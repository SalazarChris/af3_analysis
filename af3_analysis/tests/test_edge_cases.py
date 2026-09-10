"""
Regression tests for edge cases: error handling, malformed input, and seed parsing.

Particularly focuses on:
- Requirement 7: Missing/malformed input produces explicit errors
- Requirement 8: Seed parsing does not silently assign invalid replicates to seed 1
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from af3inputbuilder.scripts.af3_condition_centric_extraction import extract_metrics


# ---------------------------------------------------------------------------
# Tests: empty / missing input
# ---------------------------------------------------------------------------


class TestEmptyInput:
    """Verify behaviour with empty or missing input directories."""

    def test_empty_directory(self, tmp_output):
        """An empty directory should produce zero conditions, not crash."""
        empty_dir = tmp_output / "empty_input"
        empty_dir.mkdir()
        result = extract_metrics(str(empty_dir), str(tmp_output / "out"), verbose=False)
        assert result["conditions"] == 0
        assert result["replicates"] == 0

    def test_nonexistent_directory(self, tmp_path):
        """A nonexistent directory should not raise an unhandled exception."""
        fake_path = str(tmp_path / "does_not_exist")
        # The extraction script doesn't check existence upfront,
        # but Path.iterdir() will raise. Verify it doesn't produce
        # silent wrong output.
        with pytest.raises((FileNotFoundError, OSError)):
            extract_metrics(fake_path, str(tmp_path / "out"), verbose=False)


# ---------------------------------------------------------------------------
# Tests: malformed confidence files
# ---------------------------------------------------------------------------


class TestMalformedInput:
    """Verify that malformed JSON files produce errors, not silent corruption."""

    def test_malformed_json_in_confidences(self, tmp_path):
        """A malformed JSON file should be reported as an error."""
        cond_dir = tmp_path / "bad_json" / "test_cond"
        cond_dir.mkdir(parents=True)

        # Write _data.json
        (cond_dir / "test_cond_data.json").write_text(json.dumps({"name": "test"}))

        # Write a subdirectory with malformed confidences
        sub = cond_dir / "seed-1_sample-0"
        sub.mkdir()
        (sub / "test_cond_seed-1_sample-0_confidences.json").write_text("NOT VALID JSON {{{")

        # Also write a ranking scores file so extraction doesn't fail on that
        (cond_dir / "test_cond_ranking_scores.csv").write_text(
            "seed,sample,ranking_score\n1,0,0.5\n"
        )

        result = extract_metrics(str(tmp_path / "bad_json"), str(tmp_path / "out"), verbose=False)
        # Should report an error for the malformed file
        assert len(result["errors"]) > 0
        assert any("Failed to read" in e for e in result["errors"])

    def test_empty_confidences_file(self, tmp_path):
        """An empty confidences file should produce an error."""
        cond_dir = tmp_path / "empty_conf" / "test_cond"
        cond_dir.mkdir(parents=True)

        (cond_dir / "test_cond_data.json").write_text(json.dumps({"name": "test"}))

        sub = cond_dir / "seed-1_sample-0"
        sub.mkdir()
        (sub / "test_cond_seed-1_sample-0_confidences.json").write_text("")

        (cond_dir / "test_cond_ranking_scores.csv").write_text(
            "seed,sample,ranking_score\n1,0,0.5\n"
        )

        result = extract_metrics(str(tmp_path / "empty_conf"), str(tmp_path / "out"), verbose=False)
        assert len(result["errors"]) > 0


# ---------------------------------------------------------------------------
# Tests: missing _data.json (condition name fallback)
# ---------------------------------------------------------------------------


class TestDataJsonFallback:
    """Verify that missing _data.json uses folder name as condition name."""

    def test_fallback_to_folder_name(self, tmp_path):
        """Without _data.json, condition_name should equal folder name."""
        cond_dir = tmp_path / "no_data_json" / "my_custom_condition"
        cond_dir.mkdir(parents=True)

        # Write confidences without a _data.json
        sub = cond_dir / "seed-1_sample-0"
        sub.mkdir()
        conf = {
            "atom_plddts": [80.0],
            "contact_probs": [0.5],
            "pae": [3.0],
            "token_chain_ids": ["A"],
            "atom_chain_ids": ["A"],
        }
        (sub / "my_custom_condition_seed-1_sample-0_confidences.json").write_text(
            json.dumps(conf)
        )

        (cond_dir / "my_custom_condition_ranking_scores.csv").write_text(
            "seed,sample,ranking_score\n1,0,0.5\n"
        )

        result = extract_metrics(
            str(tmp_path / "no_data_json"), str(tmp_path / "out"), verbose=False
        )
        assert result["conditions"] == 1
        assert result["replicates"] == 1

        # Check the condition name falls back to folder name
        registry = pd.read_csv(tmp_path / "out" / "condition_registry.csv")
        assert registry["condition_name"].iloc[0] == "my_custom_condition"


# ---------------------------------------------------------------------------
# Tests: seed parsing (Requirement 8)
# ---------------------------------------------------------------------------


class TestSeedParsing:
    """
    Tests for the module-level extract_seed() function (P7 fix).

    Valid replicate IDs must produce the correct integer seed.
    Invalid/malformed IDs must raise ValueError.
    """

    def test_valid_standard_format(self):
        """Standard replicate_id formats parse to correct seeds."""
        from af3_analysis.pipeline import extract_seed

        assert extract_seed("pou_baseline_seed-1_sample-0") == 1
        assert extract_seed("pou_baseline_seed-10_sample-4") == 10
        assert extract_seed("cond_001_seed-42_sample-0") == 42
        assert extract_seed("cond_seed-3_sample-0") == 3

    def test_valid_large_seed_number(self):
        """Multi-digit seed numbers parse correctly."""
        from af3_analysis.pipeline import extract_seed

        assert extract_seed("x_seed-999_sample-0") == 999
        assert extract_seed("x_seed-12345_sample-1") == 12345

    def test_valid_seed_at_start(self):
        """seed-N can appear at the start of the replicate_id."""
        from af3_analysis.pipeline import extract_seed

        assert extract_seed("seed-5_sample-2") == 5

    def test_invalid_no_seed_pattern(self):
        """replicate_id with no seed- prefix raises ValueError."""
        from af3_analysis.pipeline import extract_seed

        with pytest.raises(ValueError, match="no 'seed-<N>' segment found"):
            extract_seed("no_seed_here")

    def test_invalid_empty_string(self):
        """Empty replicate_id raises ValueError."""
        from af3_analysis.pipeline import extract_seed

        with pytest.raises(ValueError, match="no 'seed-<N>' segment found"):
            extract_seed("")

    def test_invalid_seed_non_numeric(self):
        """seed-N where N is not a number raises ValueError."""
        from af3_analysis.pipeline import extract_seed

        with pytest.raises(ValueError, match="no 'seed-<N>' segment found"):
            extract_seed("cond_seed-abc_sample-0")

    def test_invalid_condition_only(self):
        """A bare condition name without seed raises ValueError."""
        from af3_analysis.pipeline import extract_seed

        with pytest.raises(ValueError, match="no 'seed-<N>' segment found"):
            extract_seed("condition_only")

    def test_pipeline_excludes_bad_replicate_ids(self, tmp_path):
        """
        When metrics_replicates.csv contains a mix of valid and invalid
        replicate IDs, the pipeline should exclude the invalid rows
        (with a warning) rather than silently assigning seed=1.
        """
        from af3_analysis.pipeline import _stage_run_analysis
        from af3_analysis.config import AnalysisConfig

        tables_dir = tmp_path / "tables"
        tables_dir.mkdir()

        rows = []
        # 3 valid replicates across 2 seeds
        for seed in [1, 2]:
            rows.append({
                "condition_id": "cond_001",
                "condition_name": "test",
                "replicate_id": f"test_seed-{seed}_sample-0",
                "pLDDT_mean": 70.0 + seed,
            })
        # 1 invalid replicate (no seed pattern)
        rows.append({
            "condition_id": "cond_001",
            "condition_name": "test",
            "replicate_id": "test_no_seed_pattern",
            "pLDDT_mean": 99.0,
        })
        pd.DataFrame(rows).to_csv(tables_dir / "metrics_replicates.csv", index=False)

        config = AnalysisConfig(run_id="bad_id_test", output_root=tmp_path, random_seed=42)
        result = _stage_run_analysis(tmp_path, config)

        seed_agg = pd.read_csv(tables_dir / "seed_aggregated.csv")
        # Only the 2 valid rows should appear, grouped by seed
        assert len(seed_agg) == 2
        seeds = sorted(seed_agg["seed"].unique())
        assert seeds == [1, 2]
        # The invalid row (pLDDT_mean=99) should NOT appear in any seed
        assert 99.0 not in seed_agg["pLDDT_mean"].values

    def test_pipeline_succeeds_with_all_valid_ids(self, tmp_path):
        """
        When all replicate IDs are valid, the pipeline should produce
        correct seed aggregation without warnings.
        """
        from af3_analysis.pipeline import _stage_run_analysis
        from af3_analysis.config import AnalysisConfig

        tables_dir = tmp_path / "tables"
        tables_dir.mkdir()

        rows = []
        for seed in [1, 2, 3]:
            for sample in [0, 1]:
                rows.append({
                    "condition_id": "cond_001",
                    "condition_name": "test",
                    "replicate_id": f"test_seed-{seed}_sample-{sample}",
                    "pLDDT_mean": 70.0 + seed + sample * 0.1,
                })
        pd.DataFrame(rows).to_csv(tables_dir / "metrics_replicates.csv", index=False)

        config = AnalysisConfig(run_id="all_valid_test", output_root=tmp_path, random_seed=42)
        result = _stage_run_analysis(tmp_path, config)

        seed_agg = pd.read_csv(tables_dir / "seed_aggregated.csv")
        seeds = sorted(seed_agg["seed"].unique())
        assert seeds == [1, 2, 3], f"Expected seeds [1, 2, 3], got {seeds}"
        assert len(seed_agg) == 3, f"Expected 3 rows (one per seed), got {len(seed_agg)}"


# ---------------------------------------------------------------------------
# Tests: metrics_replicates.csv with missing values
# ---------------------------------------------------------------------------


class TestMissingValues:
    """Verify pipeline handles missing metric values correctly."""

    def test_nan_metrics_produced_correctly(self, tmp_path):
        """NaN values in metrics should propagate correctly through aggregation."""
        tables_dir = tmp_path / "tables"
        tables_dir.mkdir()

        # Create metrics_replicates.csv with some NaN values
        rows = []
        for seed in [1, 2]:
            for sample in [0, 1]:
                row = {
                    "condition_id": "cond_001",
                    "condition_name": "test",
                    "replicate_id": f"test_seed-{seed}_sample-{sample}",
                    "pLDDT_mean": 80.0 if seed == 1 else float("nan"),
                    "pae_mean": 3.0,
                }
                rows.append(row)
        pd.DataFrame(rows).to_csv(tables_dir / "metrics_replicates.csv", index=False)

        from af3_analysis.pipeline import _stage_run_analysis
        from af3_analysis.config import AnalysisConfig
        config = AnalysisConfig(run_id="nan_test", output_root=tmp_path, random_seed=42)
        _stage_run_analysis(tmp_path, config)

        seed_agg = pd.read_csv(tables_dir / "seed_aggregated.csv")
        # seed-2 rows should have NaN for pLDDT_mean
        seed2_rows = seed_agg[seed_agg["seed"] == 2]
        assert seed2_rows["pLDDT_mean"].isna().all(), \
            "seed-2 should have NaN pLDDT_mean but didn't"
