"""
Tests for the wide-to-long adapter (af3_analysis.io.wide_to_long).

Verifies that the adapter correctly converts the pipeline's wide-format
seed_aggregated output into the long format expected by statistics/ and
exploratory/ modules.
"""

import numpy as np
import pandas as pd
import pytest

from af3_analysis.io.wide_to_long import (
    seed_aggregated_to_long,
    long_to_seed_aggregated,
    load_seed_aggregated_long,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_wide_df():
    """A minimal wide-format DataFrame mimicking seed_aggregated.csv."""
    return pd.DataFrame({
        "condition_id": ["c1", "c1", "c2", "c2"],
        "condition_name": ["alpha", "alpha", "beta", "beta"],
        "seed": [1, 2, 1, 2],
        "pLDDT_mean": [80.0, 82.0, 75.0, 77.0],
        "pae_mean": [3.0, 2.5, 4.0, 3.5],
        "ranking_score": [0.9, 0.85, 0.7, 0.65],
    })


@pytest.fixture
def wide_with_chain_columns():
    """Wide DataFrame with both global and chain-level metrics."""
    return pd.DataFrame({
        "condition_id": ["c1", "c2"],
        "condition_name": ["alpha", "beta"],
        "seed": [1, 1],
        "pLDDT_mean": [80.0, 75.0],
        "pae_mean": [3.0, 4.0],
        "chain_A_plddt": [82.0, 76.0],
        "chain_A_residues": [153.0, 153.0],
        "chain_B_plddt": [np.nan, 70.0],
    })


@pytest.fixture
def wide_with_nans():
    """Wide DataFrame with NaN values in some metric columns."""
    return pd.DataFrame({
        "condition_id": ["c1", "c1", "c2"],
        "condition_name": ["alpha", "alpha", "beta"],
        "seed": [1, 2, 1],
        "pLDDT_mean": [80.0, np.nan, 75.0],
        "pae_mean": [3.0, 2.5, np.nan],
    })


@pytest.fixture
def real_seed_aggregated():
    """Load the actual seed_aggregated.csv from test data."""
    from pathlib import Path
    from af3inputbuilder.scripts.af3_condition_centric_extraction import extract_metrics
    from af3_analysis.pipeline import _stage_run_analysis
    from af3_analysis.config import AnalysisConfig

    testdata = str(Path(__file__).resolve().parents[2] / "testdata" / "pou2")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        tables_dir = Path(td) / "tables"
        tables_dir.mkdir()
        extract_metrics(testdata, str(tables_dir), verbose=False)
        config = AnalysisConfig(run_id="test", output_root=Path(td), random_seed=42)
        _stage_run_analysis(Path(td), config)
        yield pd.read_csv(tables_dir / "seed_aggregated.csv")


# ---------------------------------------------------------------------------
# Tests: basic conversion
# ---------------------------------------------------------------------------


class TestBasicConversion:
    """Verify the adapter produces correct long-format output."""

    def test_output_columns(self, minimal_wide_df):
        """Long format must have exactly these columns."""
        long = seed_aggregated_to_long(minimal_wide_df)
        expected_cols = [
            "condition_id", "condition_name", "seed",
            "metric_id", "scope_type", "scope_id", "value",
        ]
        assert list(long.columns) == expected_cols

    def test_row_count(self, minimal_wide_df):
        """3 metrics x 4 rows = 12 long-format rows."""
        long = seed_aggregated_to_long(minimal_wide_df)
        assert len(long) == 12

    def test_condition_ids_preserved(self, minimal_wide_df):
        """All condition_ids from the wide format should appear."""
        long = seed_aggregated_to_long(minimal_wide_df)
        assert set(long["condition_id"].unique()) == {"c1", "c2"}

    def test_seeds_preserved(self, minimal_wide_df):
        """All seed values should appear."""
        long = seed_aggregated_to_long(minimal_wide_df)
        assert set(long["seed"].unique()) == {1, 2}

    def test_metric_ids_correct(self, minimal_wide_df):
        """metric_id should list the 3 metric column names."""
        long = seed_aggregated_to_long(minimal_wide_df)
        assert set(long["metric_id"].unique()) == {
            "pLDDT_mean", "pae_mean", "ranking_score",
        }

    def test_all_global_scope(self, minimal_wide_df):
        """All metrics in the minimal fixture are global (no chain columns)."""
        long = seed_aggregated_to_long(minimal_wide_df)
        assert (long["scope_type"] == "global").all()
        assert (long["scope_id"] == "").all()

    def test_values_match_wide(self, minimal_wide_df):
        """Each long-format value should match the corresponding wide cell."""
        long = seed_aggregated_to_long(minimal_wide_df)
        for _, row in long.iterrows():
            wide_val = minimal_wide_df.loc[
                (minimal_wide_df["condition_id"] == row["condition_id"])
                & (minimal_wide_df["seed"] == row["seed"]),
                row["metric_id"],
            ].iloc[0]
            if pd.isna(wide_val):
                assert row["value"] is None or pd.isna(row["value"])
            else:
                assert row["value"] == pytest.approx(wide_val)


# ---------------------------------------------------------------------------
# Tests: chain-level metrics
# ---------------------------------------------------------------------------


class TestChainMetrics:
    """Verify chain-level metrics are correctly classified."""

    def test_chain_metrics_detected(self, wide_with_chain_columns):
        """chain_A_plddt and chain_B_plddt should appear as chain metrics."""
        long = seed_aggregated_to_long(wide_with_chain_columns)
        chain_rows = long[long["scope_type"] == "chain"]
        metric_ids = set(chain_rows["metric_id"].unique())
        assert "plddt" in metric_ids  # chain_A_plddt -> metric_id="plddt"

    def test_chain_letters_correct(self, wide_with_chain_columns):
        """Chain A and B should be identified correctly."""
        long = seed_aggregated_to_long(wide_with_chain_columns)
        chain_rows = long[long["scope_type"] == "chain"]
        chain_ids = set(chain_rows["scope_id"].unique())
        assert "A" in chain_ids
        assert "B" in chain_ids

    def test_residues_excluded(self, wide_with_chain_columns):
        """chain_*_residues columns should NOT appear as metrics."""
        long = seed_aggregated_to_long(wide_with_chain_columns)
        # Check no metric_id contains "residues"
        assert not any("residues" in mid for mid in long["metric_id"].unique())

    def test_chain_nan_preserved(self, wide_with_chain_columns):
        """NaN values in chain metrics should be preserved as None."""
        long = seed_aggregated_to_long(wide_with_chain_columns)
        # chain_B_plddt for c1 is NaN -> should be None in long format
        c1_chain_b = long[
            (long["condition_id"] == "c1")
            & (long["scope_type"] == "chain")
            & (long["scope_id"] == "B")
            & (long["metric_id"] == "plddt")
        ]
        assert len(c1_chain_b) == 1
        val = c1_chain_b["value"].iloc[0]
        assert val is None or pd.isna(val)

    def test_global_and_chain_count(self, wide_with_chain_columns):
        """Should have 2 global metrics + 2 chain metrics = 4 distinct metric_ids."""
        long = seed_aggregated_to_long(wide_with_chain_columns)
        global_count = len(long[long["scope_type"] == "global"]["metric_id"].unique())
        chain_count = len(long[long["scope_type"] == "chain"]["metric_id"].unique())
        assert global_count == 2  # pLDDT_mean, pae_mean
        assert chain_count == 1   # plddt (from chain_A and chain_B)


# ---------------------------------------------------------------------------
# Tests: NaN handling
# ---------------------------------------------------------------------------


class TestNaNHandling:
    """Verify NaN values are handled correctly."""

    def test_nan_values_become_none(self, wide_with_nans):
        """NaN metric values should be stored as None in the long format."""
        long = seed_aggregated_to_long(wide_with_nans)
        # c1 seed-2 has NaN pLDDT_mean
        row = long[
            (long["condition_id"] == "c1")
            & (long["seed"] == 2)
            & (long["metric_id"] == "pLDDT_mean")
        ]
        assert len(row) == 1
        val = row["value"].iloc[0]
        assert val is None or pd.isna(val)


# ---------------------------------------------------------------------------
# Tests: exclude condition_name option
# ---------------------------------------------------------------------------


class TestExcludeConditionName:
    """Verify the include_condition_name parameter."""

    def test_exclude_condition_name(self, minimal_wide_df):
        """When include_condition_name=False, condition_name is not in output."""
        long = seed_aggregated_to_long(
            minimal_wide_df, include_condition_name=False
        )
        assert "condition_name" not in long.columns

    def test_include_condition_name_default(self, minimal_wide_df):
        """By default, condition_name is included."""
        long = seed_aggregated_to_long(minimal_wide_df)
        assert "condition_name" in long.columns
        assert set(long["condition_name"].unique()) == {"alpha", "beta"}


# ---------------------------------------------------------------------------
# Tests: round-trip
# ---------------------------------------------------------------------------


class TestRoundTrip:
    """Verify long -> wide -> long produces consistent results."""

    def test_round_trip_preserves_values(self, minimal_wide_df):
        """Converting long -> wide should recover the original wide format."""
        long = seed_aggregated_to_long(minimal_wide_df)
        wide_back = long_to_seed_aggregated(long)

        # Sort both for comparison
        original = minimal_wide_df.sort_values(
            ["condition_id", "seed"]
        ).reset_index(drop=True)
        recovered = wide_back.sort_values(
            ["condition_id", "seed"]
        ).reset_index(drop=True)

        # Compare metric columns
        for col in ["pLDDT_mean", "pae_mean", "ranking_score"]:
            pd.testing.assert_series_equal(
                original[col], recovered[col], check_names=False
            )

    def test_round_trip_preserves_row_count(self, wide_with_chain_columns):
        """Row count should be preserved through round-trip."""
        long = seed_aggregated_to_long(wide_with_chain_columns)
        wide_back = long_to_seed_aggregated(long)
        assert len(wide_back) == len(wide_with_chain_columns)


# ---------------------------------------------------------------------------
# Tests: input validation
# ---------------------------------------------------------------------------


class TestInputValidation:
    """Verify the adapter rejects invalid inputs."""

    def test_missing_condition_id_raises(self):
        """Missing condition_id column should raise ValueError."""
        df = pd.DataFrame({"seed": [1], "pLDDT_mean": [80.0]})
        with pytest.raises(ValueError, match="missing required columns"):
            seed_aggregated_to_long(df)

    def test_missing_seed_raises(self):
        """Missing seed column should raise ValueError."""
        df = pd.DataFrame({
            "condition_id": ["c1"],
            "condition_name": ["x"],
            "pLDDT_mean": [80.0],
        })
        with pytest.raises(ValueError, match="missing required columns"):
            seed_aggregated_to_long(df)

    def test_non_numeric_columns_ignored(self):
        """Non-numeric columns that aren't id columns should be silently ignored."""
        df = pd.DataFrame({
            "condition_id": ["c1"],
            "condition_name": ["x"],
            "seed": [1],
            "pLDDT_mean": [80.0],
            "some_text_col": ["hello"],
        })
        long = seed_aggregated_to_long(df)
        assert "some_text_col" not in long["metric_id"].values

    def test_empty_dataframe(self):
        """An empty DataFrame should produce an empty long DataFrame."""
        df = pd.DataFrame(columns=["condition_id", "condition_name", "seed", "pLDDT_mean"])
        long = seed_aggregated_to_long(df)
        assert len(long) == 0
        assert list(long.columns) == [
            "condition_id", "condition_name", "seed",
            "metric_id", "scope_type", "scope_id", "value",
        ]


# ---------------------------------------------------------------------------
# Tests: real data
# ---------------------------------------------------------------------------


class TestRealData:
    """Integration tests using the actual pipeline output."""

    def test_long_format_schema(self, real_seed_aggregated):
        """Real data should produce the expected long-format schema."""
        long = seed_aggregated_to_long(real_seed_aggregated)
        expected_cols = [
            "condition_id", "condition_name", "seed",
            "metric_id", "scope_type", "scope_id", "value",
        ]
        assert list(long.columns) == expected_cols

    def test_global_metrics_present(self, real_seed_aggregated):
        """All expected global metrics should appear."""
        long = seed_aggregated_to_long(real_seed_aggregated)
        global_metrics = set(
            long[long["scope_type"] == "global"]["metric_id"].unique()
        )
        expected_globals = {
            "pLDDT_mean", "pLDDT_max", "pLDDT_min", "pLDDT_median",
            "contact_prob_mean", "contact_prob_max", "contact_prob_min", "contact_prob_median",
            "pae_mean", "pae_max", "pae_min", "pae_median",
            "ranking_score",
        }
        assert expected_globals.issubset(global_metrics)

    def test_chain_metrics_present(self, real_seed_aggregated):
        """Chain-level metrics (plddt) should be detected."""
        long = seed_aggregated_to_long(real_seed_aggregated)
        chain_metrics = set(
            long[long["scope_type"] == "chain"]["metric_id"].unique()
        )
        assert "plddt" in chain_metrics

    def test_chain_letters(self, real_seed_aggregated):
        """Chain A should be present in chain-level data."""
        long = seed_aggregated_to_long(real_seed_aggregated)
        chain_scope_ids = set(
            long[long["scope_type"] == "chain"]["scope_id"].unique()
        )
        assert "A" in chain_scope_ids

    def test_no_residues_in_metrics(self, real_seed_aggregated):
        """chain_*_residues should not appear as metric_ids."""
        long = seed_aggregated_to_long(real_seed_aggregated)
        assert not any("residues" in mid for mid in long["metric_id"].unique())

    def test_row_count_consistent(self, real_seed_aggregated):
        """Long format should have n_rows_wide * n_metrics rows."""
        wide = real_seed_aggregated
        n_global = len([c for c in wide.columns if c not in
                        ["condition_id", "condition_name", "seed"]
                        and pd.api.types.is_numeric_dtype(wide[c])
                        and not c.startswith("chain_")])
        n_chain = len([c for c in wide.columns
                       if c.startswith("chain_") and c.endswith("_plddt")
                       and pd.api.types.is_numeric_dtype(wide[c])])
        n_rows = len(wide)
        expected = n_rows * (n_global + n_chain)
        long = seed_aggregated_to_long(real_seed_aggregated)
        assert len(long) == expected

    def test_round_trip_real_data(self, real_seed_aggregated):
        """Round-trip on real data should recover original metric values."""
        long = seed_aggregated_to_long(real_seed_aggregated)
        wide_back = long_to_seed_aggregated(long)

        original = real_seed_aggregated.sort_values(
            ["condition_id", "seed"]
        ).reset_index(drop=True)
        recovered = wide_back.sort_values(
            ["condition_id", "seed"]
        ).reset_index(drop=True)

        for col in ["pLDDT_mean", "pae_mean", "ranking_score"]:
            if col in original.columns and col in recovered.columns:
                pd.testing.assert_series_equal(
                    original[col], recovered[col], check_names=False
                )

    def test_load_function(self, real_seed_aggregated, tmp_path):
        """load_seed_aggregated_long should load and convert correctly."""
        csv_path = tmp_path / "seed_aggregated.csv"
        real_seed_aggregated.to_csv(csv_path, index=False)

        long = load_seed_aggregated_long(str(csv_path))
        assert len(long) > 0
        assert "metric_id" in long.columns
        assert "scope_type" in long.columns
