"""
Regression tests for Stage 4: Analysis (seed aggregation + descriptive stats).

Verifies that the pipeline correctly reads metrics_replicates.csv, computes
seed-level aggregation, and produces correct mean/SD values.
"""

import numpy as np
import pandas as pd
import pytest

from af3inputbuilder.scripts.af3_condition_centric_extraction import extract_metrics


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run_extraction_to_tables(input_dir, tables_dir):
    """Run extraction and return the tables directory Path."""
    extract_metrics(str(input_dir), str(tables_dir), verbose=False)
    return tables_dir


def _run_analysis_stage(tables_dir):
    """Run the pipeline's Stage 4 analysis on pre-extracted tables."""
    from af3_analysis.pipeline import _stage_run_analysis
    from af3_analysis.config import AnalysisConfig
    from pathlib import Path

    run_dir = tables_dir.parent
    config = AnalysisConfig(
        run_id="test_run",
        output_root=run_dir.parent,
        random_seed=42,
    )
    result = _stage_run_analysis(run_dir, config)
    return result


# ---------------------------------------------------------------------------
# Tests: seed_aggregated.csv
# ---------------------------------------------------------------------------


class TestSeedAggregated:
    """Verify seed_aggregated.csv structure and row count."""

    @pytest.fixture(autouse=True)
    def run_pipeline(self, synthetic_confidence_dir, tmp_output):
        tables_dir = tmp_output / "tables"
        tables_dir.mkdir()
        _run_extraction_to_tables(synthetic_confidence_dir, tables_dir)
        _run_analysis_stage(tables_dir)
        self.tables = tables_dir

    def test_seed_aggregated_exists(self):
        """seed_aggregated.csv must be produced."""
        assert (self.tables / "seed_aggregated.csv").exists()

    def test_seed_aggregated_row_count(self):
        """
        One row per condition x seed.
        synthetic_confidence_dir: 1 condition x 2 seeds = 2 rows.
        """
        df = pd.read_csv(self.tables / "seed_aggregated.csv")
        assert len(df) == 2

    def test_seed_aggregated_has_condition_columns(self):
        """Must contain condition_id, condition_name, seed."""
        df = pd.read_csv(self.tables / "seed_aggregated.csv")
        for col in ["condition_id", "condition_name", "seed"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_seed_values_are_correct(self):
        """Seed column should contain the expected seed IDs."""
        df = pd.read_csv(self.tables / "seed_aggregated.csv")
        seeds = sorted(df["seed"].unique())
        assert seeds == [1, 2]

    def test_metric_values_present(self):
        """Core metric columns should have numeric values."""
        df = pd.read_csv(self.tables / "seed_aggregated.csv")
        for col in ["pLDDT_mean", "pae_mean", "contact_prob_mean"]:
            if col in df.columns:
                assert pd.api.types.is_numeric_dtype(df[col]), f"{col} not numeric"
                assert df[col].notna().any(), f"{col} is all NaN"


# ---------------------------------------------------------------------------
# Tests: descriptive_stats.csv
# ---------------------------------------------------------------------------


class TestDescriptiveStats:
    """Verify descriptive_stats.csv is produced correctly."""

    @pytest.fixture(autouse=True)
    def run_pipeline(self, synthetic_confidence_dir, tmp_output):
        tables_dir = tmp_output / "tables"
        tables_dir.mkdir()
        _run_extraction_to_tables(synthetic_confidence_dir, tables_dir)
        _run_analysis_stage(tables_dir)
        self.tables = tables_dir

    def test_descriptive_stats_exists(self):
        """descriptive_stats.csv must be produced."""
        assert (self.tables / "descriptive_stats.csv").exists()

    def test_descriptive_stats_has_mean_and_std(self):
        """Should contain _mean and _std suffixed columns."""
        df = pd.read_csv(self.tables / "descriptive_stats.csv")
        mean_cols = [c for c in df.columns if c.endswith("_mean")]
        std_cols = [c for c in df.columns if c.endswith("_std")]
        assert len(mean_cols) > 0, "No _mean columns found"
        assert len(std_cols) > 0, "No _std columns found"

    def test_descriptive_stats_row_count(self):
        """One row per condition. synthetic has 1 condition."""
        df = pd.read_csv(self.tables / "descriptive_stats.csv")
        assert len(df) == 1


# ---------------------------------------------------------------------------
# Tests: mean/SD calculations
# ---------------------------------------------------------------------------


class TestMeanSDCalculations:
    """Verify that mean and SD calculations are numerically correct."""

    @pytest.fixture(autouse=True)
    def run_pipeline(self, synthetic_confidence_dir, tmp_output):
        tables_dir = tmp_output / "tables"
        tables_dir.mkdir()
        _run_extraction_to_tables(synthetic_confidence_dir, tables_dir)
        _run_analysis_stage(tables_dir)
        self.seed_agg = pd.read_csv(tables_dir / "seed_aggregated.csv")
        self.desc_stats = pd.read_csv(tables_dir / "descriptive_stats.csv")

    def test_condition_mean_matches_seed_means(self):
        """Condition-level mean should equal the mean of seed-level means."""
        seed_means = self.seed_agg["pLDDT_mean"]
        expected_mean = seed_means.mean()

        col = "pLDDT_mean_mean"
        if col in self.desc_stats.columns:
            actual_mean = self.desc_stats[col].iloc[0]
            assert abs(actual_mean - expected_mean) < 1e-10, \
                f"Mean mismatch: expected {expected_mean}, got {actual_mean}"

    def test_condition_std_is_seed_variance(self):
        """Condition-level SD should reflect between-seed variation."""
        seed_values = self.seed_agg["pLDDT_mean"]
        if len(seed_values) > 1:
            expected_std = seed_values.std(ddof=1)
            col = "pLDDT_mean_std"
            if col in self.desc_stats.columns:
                actual_std = self.desc_stats[col].iloc[0]
                assert abs(actual_std - expected_std) < 1e-10, \
                    f"SD mismatch: expected {expected_std}, got {actual_std}"


# ---------------------------------------------------------------------------
# Tests: multi-condition seed aggregation
# ---------------------------------------------------------------------------


class TestMultiConditionAggregation:
    """Verify seed aggregation across multiple conditions."""

    @pytest.fixture(autouse=True)
    def run_pipeline(self, synthetic_multi_condition_dir, tmp_output):
        tables_dir = tmp_output / "tables"
        tables_dir.mkdir()
        _run_extraction_to_tables(synthetic_multi_condition_dir, tables_dir)
        _run_analysis_stage(tables_dir)
        self.seed_agg = pd.read_csv(tables_dir / "seed_aggregated.csv")

    def test_row_count_two_conditions(self):
        """2 conditions x 2 seeds = 4 rows."""
        assert len(self.seed_agg) == 4

    def test_both_conditions_present(self):
        """Both condition_ids should appear in seed_aggregated."""
        condition_ids = set(self.seed_agg["condition_id"].unique())
        assert len(condition_ids) == 2

    def test_each_condition_has_two_seeds(self):
        """Each condition should have exactly 2 seed rows."""
        for cond_id in self.seed_agg["condition_id"].unique():
            cond_seeds = self.seed_agg[self.seed_agg["condition_id"] == cond_id]
            assert len(cond_seeds) == 2


# ---------------------------------------------------------------------------
# Tests: real data (testdata/pou2/)
# ---------------------------------------------------------------------------


class TestAnalysisRealData:
    """Integration tests using real testdata/pou2/ data.

    testdata/pou2/ has 9 directories (8 conditions + metrics_example).
    The extraction produces 9 conditions with varying replicate counts.
    After seed aggregation, the result should have one row per
    (condition x unique seed).
    """

    @pytest.fixture(autouse=True)
    def run_pipeline(self, pou_parent_dir, tmp_output):
        tables_dir = tmp_output / "tables"
        tables_dir.mkdir()
        _run_extraction_to_tables(pou_parent_dir, tables_dir)
        _run_analysis_stage(tables_dir)
        self.seed_agg = pd.read_csv(tables_dir / "seed_aggregated.csv")
        self.desc_stats = pd.read_csv(tables_dir / "descriptive_stats.csv")

    def test_seed_aggregated_has_rows(self):
        """Should produce at least 80 seed-level rows (8 conditions x 10 seeds)."""
        assert len(self.seed_agg) >= 80

    def test_seed_values_range(self):
        """Seeds should range from 1 to at least 10."""
        seeds = sorted(self.seed_agg["seed"].unique())
        assert 1 in seeds
        assert 10 in seeds

    def test_plddt_mean_reasonable_range(self):
        """pLDDT mean should be between 0 and 100."""
        assert (self.seed_agg["pLDDT_mean"] >= 0).all()
        assert (self.seed_agg["pLDDT_mean"] <= 100).all()

    def test_descriptive_stats_has_rows(self):
        """Should produce at least 8 condition-level summary rows."""
        assert len(self.desc_stats) >= 8

    def test_descriptive_stats_has_mean_std_pairs(self):
        """For pLDDT_mean, both _mean and _std should exist."""
        assert "pLDDT_mean_mean" in self.desc_stats.columns
        assert "pLDDT_mean_std" in self.desc_stats.columns
