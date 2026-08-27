"""
Regression tests for Stage 1: Extraction.

Verifies that the condition-centric extraction script produces the expected
CSV files with correct structure and content.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from af3inputbuilder.scripts.af3_condition_centric_extraction import extract_metrics


# ---------------------------------------------------------------------------
# Expected CSV filenames
# ---------------------------------------------------------------------------

EXPECTED_CSVS = [
    "condition_registry.csv",
    "metrics_replicates.csv",
    "metrics_conditions.csv",
    "condition_manifest.csv",
]


# ---------------------------------------------------------------------------
# Tests: output file existence
# ---------------------------------------------------------------------------


class TestExtractionProducesExpectedFiles:
    """Verify that extraction writes all four required CSVs."""

    def test_all_csvs_exist(self, synthetic_confidence_dir, tmp_output):
        """Extraction must produce all four CSV files."""
        result = extract_metrics(str(synthetic_confidence_dir), str(tmp_output), verbose=False)

        assert result["conditions"] == 1
        for csv_name in EXPECTED_CSVS:
            assert (tmp_output / csv_name).exists(), f"Missing: {csv_name}"

    def test_no_errors_on_valid_input(self, synthetic_confidence_dir, tmp_output):
        """Valid input should produce zero extraction errors."""
        result = extract_metrics(str(synthetic_confidence_dir), str(tmp_output), verbose=False)
        assert result["errors"] == [], f"Unexpected errors: {result['errors']}"


# ---------------------------------------------------------------------------
# Tests: metrics_replicates.csv structure
# ---------------------------------------------------------------------------


class TestMetricsReplicatesStructure:
    """Verify the schema of metrics_replicates.csv."""

    @pytest.fixture(autouse=True)
    def run_extraction(self, synthetic_confidence_dir, tmp_output):
        extract_metrics(str(synthetic_confidence_dir), str(tmp_output), verbose=False)
        self.df = pd.read_csv(tmp_output / "metrics_replicates.csv")

    def test_required_columns_present(self):
        """metrics_replicates.csv must have condition_id, condition_name, replicate_id."""
        for col in ["condition_id", "condition_name", "replicate_id"]:
            assert col in self.df.columns, f"Missing column: {col}"

    def test_global_metric_columns_present(self):
        """Must contain core AF3 metric columns."""
        expected = ["pLDDT_mean", "pae_mean", "contact_prob_mean"]
        for col in expected:
            assert col in self.df.columns, f"Missing metric column: {col}"

    def test_row_count_matches_replicates(self):
        """2 seeds x 2 samples = 4 rows."""
        assert len(self.df) == 4

    def test_replicate_ids_contain_seed_sample(self):
        """Each replicate_id should contain seed-N and sample-M patterns."""
        for rid in self.df["replicate_id"]:
            assert "seed-" in rid, f"replicate_id missing seed pattern: {rid}"
            assert "sample-" in rid, f"replicate_id missing sample pattern: {rid}"

    def test_metric_values_are_numeric(self):
        """All metric columns should contain valid numeric values."""
        metric_cols = [c for c in self.df.columns if c not in
                       ["condition_id", "condition_name", "replicate_id"]]
        for col in metric_cols:
            assert pd.api.types.is_numeric_dtype(self.df[col].dropna()), \
                f"Column {col} is not numeric"


# ---------------------------------------------------------------------------
# Tests: metrics_replicates.csv loads correctly
# ---------------------------------------------------------------------------


class TestMetricsReplicatesLoadable:
    """Verify metrics_replicates.csv can be loaded and is well-formed."""

    def test_loadable_as_dataframe(self, synthetic_confidence_dir, tmp_output):
        """Must load without errors."""
        extract_metrics(str(synthetic_confidence_dir), str(tmp_output), verbose=False)
        df = pd.read_csv(tmp_output / "metrics_replicates.csv")
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_no_completely_null_rows(self, synthetic_confidence_dir, tmp_output):
        """No row should be entirely NaN."""
        extract_metrics(str(synthetic_confidence_dir), str(tmp_output), verbose=False)
        df = pd.read_csv(tmp_output / "metrics_replicates.csv")
        assert df.dropna(how="all").shape[0] == df.shape[0]


# ---------------------------------------------------------------------------
# Tests: real test data (testdata/pou2/)
# ---------------------------------------------------------------------------


class TestExtractionRealData:
    """Integration tests using the real testdata/pou2/ directory.

    testdata/pou2/ contains:
    - 8 condition directories (pou_baseline, pou_dna, etc.)
    - 1 metrics_example directory (development test data)

    The extraction script treats ALL subdirectories as conditions.
    Each condition has 10 seeds x 5 samples = 50 seed-specific replicates,
    plus aggregate confidences files (pou_baseline_confidences.json,
    pou_baseline_summary_confidences.json) that also match *_confidences.json.

    Expected: 9 conditions, ~820 total replicates.
    """

    def test_pou2_extracts_conditions(self, pou_parent_dir, tmp_output):
        """testdata/pou2/ should extract multiple conditions (> 8)."""
        result = extract_metrics(str(pou_parent_dir), str(tmp_output), verbose=False)
        # 8 real conditions + metrics_example = 9 total
        assert result["conditions"] >= 8

    def test_pou2_extracts_replicates(self, pou_parent_dir, tmp_output):
        """Should extract hundreds of replicates from the real data."""
        result = extract_metrics(str(pou_parent_dir), str(tmp_output), verbose=False)
        # 8 conditions x 50 seed replicates + extras from aggregate files + metrics_example
        assert result["replicates"] >= 400

    def test_pou2_condition_names(self, pou_parent_dir, tmp_output):
        """All 8 real condition names should be discovered from _data.json files."""
        extract_metrics(str(pou_parent_dir), str(tmp_output), verbose=False)
        registry = pd.read_csv(tmp_output / "condition_registry.csv")
        names = set(registry["condition_name"])
        expected_names = {
            "pou_baseline", "pou_dna", "pou_sep102", "pou_sep102_dna",
            "pou_tpo101", "pou_tpo101_dna", "pou_tpo101_sep102", "pou_tpo101_sep102_dna",
        }
        # All expected names should be present (metrics_example may add a duplicate)
        assert expected_names.issubset(names)

    def test_pou2_replicates_csv_loadable(self, pou_parent_dir, tmp_output):
        """metrics_replicates.csv should load as a valid DataFrame."""
        extract_metrics(str(pou_parent_dir), str(tmp_output), verbose=False)
        df = pd.read_csv(tmp_output / "metrics_replicates.csv")
        assert len(df) >= 400
        assert "condition_id" in df.columns
        assert "pLDDT_mean" in df.columns

    def test_pou2_ranking_scores_populated(self, pou_parent_dir, tmp_output):
        """ranking_score column should be populated for seed/sample replicates."""
        extract_metrics(str(pou_parent_dir), str(tmp_output), verbose=False)
        df = pd.read_csv(tmp_output / "metrics_replicates.csv")
        assert "ranking_score" in df.columns
        non_null = df["ranking_score"].notna().sum()
        assert non_null >= 400, f"Expected >= 400 ranking scores, got {non_null}"


# ---------------------------------------------------------------------------
# Tests: multi-condition extraction (synthetic)
# ---------------------------------------------------------------------------


class TestMultiConditionExtraction:
    """Verify extraction handles multiple conditions."""

    def test_two_conditions(self, synthetic_multi_condition_dir, tmp_output):
        """Should discover 2 conditions."""
        result = extract_metrics(
            str(synthetic_multi_condition_dir), str(tmp_output), verbose=False
        )
        assert result["conditions"] == 2
        assert result["replicates"] == 8  # 2 conditions x 2 seeds x 2 samples

    def test_condition_registry_has_both(self, synthetic_multi_condition_dir, tmp_output):
        """condition_registry.csv should list both conditions."""
        extract_metrics(
            str(synthetic_multi_condition_dir), str(tmp_output), verbose=False
        )
        registry = pd.read_csv(tmp_output / "condition_registry.csv")
        assert len(registry) == 2
        names = set(registry["condition_name"])
        assert "condition_alpha" in names
        assert "condition_beta" in names

    def test_condition_manifest_seeds(self, synthetic_multi_condition_dir, tmp_output):
        """condition_manifest.csv should report 2 seeds per condition."""
        extract_metrics(
            str(synthetic_multi_condition_dir), str(tmp_output), verbose=False
        )
        manifest = pd.read_csv(tmp_output / "condition_manifest.csv")
        assert (manifest["seeds"] == 2).all()
