"""
Tests for the F13 similarity view family (F13A-F13D) and the
condition-level aggregation of the preserved prediction-level matrix.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from af3_analysis.visualization_v3.figures.f13_similarity_views import (
    aggregate_condition_matrix,
    within_condition_values,
    generate_f13a_condition_similarity,
    generate_f13b_within_condition_reproducibility,
    generate_f13d_prediction_matrix,
)


def _block_matrix():
    """3 conditions x 3 predictions with distinct within/between scales.

    Within-condition RMSD ~ 1.0, between-condition ~ 5.0 (A-B) / 9.0 (A-C).
    """
    n = 9
    m = np.full((n, n), np.nan)
    valid = np.zeros((n, n), dtype=bool)
    cond = ["cA"] * 3 + ["cB"] * 3 + ["cC"] * 3

    for i in range(n):
        valid[i, i] = True
        m[i, i] = 0.0
        for j in range(i + 1, n):
            same = cond[i] == cond[j]
            if same:
                d = 1.0 + 0.1 * ((i + j) % 3)
            else:
                pair = {cond[i], cond[j]}
                d = 5.0 if pair == {"cA", "cB"} else (
                    9.0 if pair == {"cA", "cC"} else 7.0)
            m[i, j] = m[j, i] = d
            valid[i, j] = valid[j, i] = True
    return m, valid, cond


class TestAggregation:
    def test_diagonal_is_within_condition_not_zero(self):
        m, valid, cond = _block_matrix()
        long_df, summary_matrix, _ = aggregate_condition_matrix(m, valid, cond)
        diag = long_df[long_df["condition_a"] == long_df["condition_b"]]
        assert (diag["comparison_type"] == "within").all()
        assert (diag["median_rmsd"] > 0).all()  # NOT auto-zero
        assert np.isclose(diag["median_rmsd"].iloc[0], 1.0, atol=0.1)

    def test_between_condition_medians(self):
        m, valid, cond = _block_matrix()
        long_df, _, _ = aggregate_condition_matrix(m, valid, cond)
        ab = long_df[(long_df.condition_a == "cA") & (long_df.condition_b == "cB")]
        assert np.isclose(ab["median_rmsd"].iloc[0], 5.0, atol=0.1)
        ac = long_df[(long_df.condition_a == "cA") & (long_df.condition_b == "cC")]
        assert np.isclose(ac["median_rmsd"].iloc[0], 9.0, atol=0.1)

    def test_symmetry_of_long_table(self):
        m, valid, cond = _block_matrix()
        long_df, _, _ = aggregate_condition_matrix(m, valid, cond)
        for _, row in long_df[long_df.comparison_type == "between"].iterrows():
            mirror = long_df[(long_df.condition_a == row.condition_b)
                             & (long_df.condition_b == row.condition_a)]
            assert len(mirror) == 1
            assert np.isclose(mirror["median_rmsd"].iloc[0],
                              row["median_rmsd"])

    def test_summary_matrix_symmetric_and_filled(self):
        m, valid, cond = _block_matrix()
        _, summary_matrix, _ = aggregate_condition_matrix(m, valid, cond)
        assert np.allclose(summary_matrix, summary_matrix.T, equal_nan=True)
        assert np.isfinite(summary_matrix).all()

    def test_self_comparison_excluded(self):
        """Self-comparisons (zero diagonal) must not leak into stats."""
        m, valid, cond = _block_matrix()
        w = within_condition_values(m, valid, cond)
        for c, vals in w.items():
            assert len(vals) == 3 * 2  # 3 preds -> 6 ordered pairs, no self
            assert (vals > 0).all()

    def test_missing_predictions_handled(self):
        """NaN/invalid entries reduce n_pairs and coverage, not correctness."""
        m, valid, cond = _block_matrix()
        # Invalidate the cA<->cB block between prediction 2 and 3.
        m[2, 3] = m[3, 2] = np.nan
        valid[2, 3] = valid[3, 2] = False
        long_df, _, _ = aggregate_condition_matrix(m, valid, cond)
        ab = long_df[(long_df.condition_a == "cA") & (long_df.condition_b == "cB")]
        row = ab.iloc[0]
        assert row["n_pairs"] == 8  # 9 ordered pairs minus the invalid one
        assert np.isclose(row["coverage"], 8 / 9)
        assert row["n_pairs_possible"] == 9

    def test_matrix_shape_mismatch_raises(self):
        m, valid, cond = _block_matrix()
        with pytest.raises(ValueError):
            aggregate_condition_matrix(m, valid, cond[:-1])

    def test_deterministic_ordering(self):
        m, valid, cond = _block_matrix()
        long1, _, _ = aggregate_condition_matrix(m, valid, cond)
        long2, _, _ = aggregate_condition_matrix(m, valid, cond)
        pd.testing.assert_frame_equal(long1, long2)


class TestFigures:
    def test_f13a_renders_and_documents_scale(self, tmp_path):
        m, valid, cond = _block_matrix()
        long_df, summary_matrix, _ = aggregate_condition_matrix(m, valid, cond)
        r = generate_f13a_condition_similarity(
            long_df, summary_matrix, sorted(set(cond)), tmp_path,
            summary_stat="median")
        assert r["status"] == "pass"
        assert r["matrix_type"] == "condition x condition"
        assert r["n_conditions"] == 3
        assert Path(r["output_path"]).is_file()

    def test_f13b_renders_with_reproducibility_terminology(self, tmp_path):
        m, valid, cond = _block_matrix()
        w = within_condition_values(m, valid, cond)
        r = generate_f13b_within_condition_reproducibility(w, tmp_path)
        assert r["status"] == "pass"
        assert any("not biological replicates" in x for x in r["warnings"])
        assert any("prediction reproducibility" in x.lower()
                   for x in r["warnings"])

    def test_f13b_skips_without_within_pairs(self, tmp_path):
        r = generate_f13b_within_condition_reproducibility(
            {"cA": np.array([])}, tmp_path)
        assert r["status"] == "skip"

    def test_f13d_renders_without_prediction_labels(self, tmp_path):
        m, valid, cond = _block_matrix()
        preds = [f"p{i}" for i in range(9)]
        seeds = [1, 2, 3] * 3
        samples = [0] * 9
        clusters = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        r = generate_f13d_prediction_matrix(
            m, valid, preds, cond, seeds, samples, tmp_path,
            cluster_labels=clusters,
            annotation_strips={
                "Seed": seeds,
                "Predicted structural cluster": (clusters + 1).tolist(),
            })
        assert r["status"] == "pass"
        assert r["matrix_type"] == "prediction x prediction"
        assert r["ordering"] == "cluster_condition_seed_sample"
        assert "existing predicted structural clusters" in r["warnings"][0]

    def test_f13d_shape_mismatch_skips(self, tmp_path):
        m, valid, cond = _block_matrix()
        r = generate_f13d_prediction_matrix(
            m, valid, ["p0"], ["cA"], [1], [0], tmp_path)
        assert r["status"] == "skip"

    def test_no_hardcoded_condition_names(self):
        """Aggregation must work for arbitrary condition identifiers."""
        m, valid, cond = _block_matrix()
        renamed = ["alpha_x", "alpha_x", "alpha_x",
                   "zz_beta", "zz_beta", "zz_beta",
                   "gamma Q9", "gamma Q9", "gamma Q9"]
        long_df, _, _ = aggregate_condition_matrix(m, valid, renamed)
        assert set(long_df["condition_a"]) == set(renamed)


class TestViewToggles:
    """F13A/F13D are hidden by default; F13C was removed entirely."""

    def _run(self, tmp_path, view_config=None):
        from af3_analysis.visualization_v3.runner import _generate_f13_views
        m, valid, cond = _block_matrix()
        return _generate_f13_views(
            dataset=None,
            matrix=m,
            valid=valid,
            predictions=[f"p{i}" for i in range(9)],
            conditions=cond,
            seeds=[1] * 9,
            samples=[0] * 9,
            cluster_labels=None,
            figures_dir=tmp_path,
            design=None,
            view_config=view_config,
        )

    def test_default_views_are_f13b_only(self, tmp_path):
        result = self._run(tmp_path)
        assert result["status"] == "pass"
        assert set(result["views"]) == {"F13B"}
        assert result["views"]["F13B"]["status"] == "pass"

    def test_hidden_views_reenable_via_config(self, tmp_path):
        result = self._run(tmp_path, view_config={"F13A": True, "F13D": True})
        assert set(result["views"]) == {"F13A", "F13B", "F13D"}

    def test_f13c_is_removed(self, tmp_path):
        """F13C must never be generated, even with everything enabled."""
        result = self._run(
            tmp_path, view_config={"F13A": True, "F13B": True, "F13D": True,
                                   "F13C": True})
        assert "F13C" not in result["views"]

    def test_all_hidden_still_preserves_artifacts(self, tmp_path):
        result = self._run(
            tmp_path,
            view_config={"F13A": False, "F13B": False, "F13D": False})
        assert result["status"] == "skip"
        assert result["views"] == {}
        # Data artifacts survive: the preserved matrix and both tables.
        assert (tmp_path.parent / "tables"
                / "pairwise_structural_distances.csv").is_file()
        assert (tmp_path.parent / "tables"
                / "prediction_metadata.csv").is_file()
        assert len(result["condition_similarity_summary"]) == 9
        assert len(result["within_condition_variability"]) == 3
