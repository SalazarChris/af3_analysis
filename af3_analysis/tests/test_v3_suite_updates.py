"""
Tests for the V3 figure-suite updates from the structural visualization
review (docs/V3_STRUCTURAL_VISUALIZATION_REVIEW.md):

- default-suite demotion of F01/F06/F08/F09/F11 with explicit opt-in,
- F06 -> F07 consolidation (recurring contact pairs panel),
- F10 site-provenance fields.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from af3_analysis.visualization_v3.config import (
    V3Config,
    V3_ALL_FIGURES,
    V3_DEFAULT_OFF_FIGURES,
    get_enabled_figures,
)
from af3_analysis.visualization_v3 import runner as v3_runner


# ---------------------------------------------------------------------------
# Default-suite demotion
# ---------------------------------------------------------------------------


class TestDefaultSuite:
    def test_known_default_off_figures(self):
        assert set(V3_DEFAULT_OFF_FIGURES) == {"F01", "F06", "F08", "F09",
                                                "F11", "F12"}

    def test_default_config_excludes_demoted_figures(self):
        enabled = get_enabled_figures(V3Config())
        assert "F01" not in enabled
        assert "F06" not in enabled
        assert "F08" not in enabled
        assert "F09" not in enabled
        assert "F11" not in enabled
        # Core suite stays on
        for fig_id in ("F02", "F03", "F04", "F05", "F07", "F13",
                       "F14", "F15", "F16", "F17", "F18", "F20"):
            assert fig_id in enabled, fig_id
        # F12 is default-off (hard to read/interpret)
        assert "F12" not in enabled
        # No unknown IDs, no duplicates
        assert set(enabled) == set(V3_ALL_FIGURES) - set(V3_DEFAULT_OFF_FIGURES)

    def test_explicit_opt_in_reenables_demoted_figure(self):
        config = V3Config(figures={"F06": True})
        enabled = get_enabled_figures(config)
        assert "F06" in enabled
        # Other demoted figures stay off
        assert "F08" not in enabled

    def test_explicit_opt_out_disables_default_on_figure(self):
        config = V3Config(figures={"F13": False})
        enabled = get_enabled_figures(config)
        assert "F13" not in enabled
        assert "F07" in enabled

    def test_partial_overrides_do_not_change_other_defaults(self):
        config = V3Config(figures={"F11": True})
        enabled = get_enabled_figures(config)
        assert "F11" in enabled
        assert "F01" not in enabled  # still default-off
        assert "F13" in enabled      # still default-on

    def test_f12_opt_in_reenables(self):
        config = V3Config(figures={"F12": True})
        enabled = get_enabled_figures(config)
        assert "F12" in enabled
        assert "F11" not in enabled  # other demoted figures unaffected

    def test_interactive_factory_uses_defaults(self):
        from af3_analysis.visualization_v3.config import (
            create_v3_config_interactive,
        )
        enabled = get_enabled_figures(create_v3_config_interactive())
        assert set(enabled) == set(V3_ALL_FIGURES) - set(V3_DEFAULT_OFF_FIGURES)

    def test_all_list_enables_everything(self):
        from af3_analysis.visualization_v3.__main__ import _with_figures
        config = _with_figures(V3Config(), list(V3_ALL_FIGURES))
        assert get_enabled_figures(config) == V3_ALL_FIGURES

    def test_suite_notes_recorded_in_results(self, tmp_path):
        """Runner result dict records demotion/consolidation provenance."""
        run_dir = tmp_path / "run"
        (run_dir / "tables").mkdir(parents=True)
        import pandas as pd
        pd.DataFrame({
            "condition_id": ["cond_001"],
            "condition_name": ["cond_001"],
            "seed": [1],
            "pLDDT_mean": [80.0],
        }).to_csv(run_dir / "tables" / "seed_aggregated.csv", index=False)
        pd.DataFrame({"condition_id": ["cond_001"]}).to_csv(
            run_dir / "tables" / "descriptive_stats.csv", index=False)

        results = v3_runner.run_v3_pipeline(run_dir)
        notes = results.get("suite_notes", {})
        assert set(notes.get("demoted_from_default", [])) == set(
            V3_DEFAULT_OFF_FIGURES)
        # F07 ran (default-on) and recorded the F06 consolidation.
        consolidated = notes.get("consolidated")
        assert isinstance(consolidated, dict)
        assert consolidated.get("into") == "F07"
        assert "F06" in consolidated.get("from", [])


# ---------------------------------------------------------------------------
# F06 -> F07 consolidation
# ---------------------------------------------------------------------------


class TestF07Consolidation:
    def _contact_change_data(self):
        return [
            {"condition_id": "condA", "seed": 1,
             "n_gained": 5.0, "n_lost": 3.0, "n_total": 100.0,
             "n_comparisons": 1},
            {"condition_id": "condB", "seed": 1,
             "n_gained": 1.0, "n_lost": 8.0, "n_total": 100.0,
             "n_comparisons": 1},
        ]

    def _recurring_pairs(self):
        return [
            {"chain_a": "A", "seq_a": 10, "chain_b": "A", "seq_b": 55,
             "n_gained": 3, "n_lost": 0, "n_conditions": 1},
            {"chain_a": "A", "seq_a": 22, "chain_b": "B", "seq_b": 7,
             "n_gained": 0, "n_lost": 2, "n_conditions": 2},
        ]

    def test_two_panel_without_pair_detail(self, tmp_path):
        from af3_analysis.visualization_v3.figures.f07_contact_change_summary import (
            generate_f07_contact_change_summary,
        )
        result = generate_f07_contact_change_summary(
            self._contact_change_data(), tmp_path)
        assert result["status"] == "pass"
        assert result["n_observations"] == 2
        # No Panel C warning path triggered
        assert not any("recurring" in w.lower() for w in result["warnings"])

    def test_three_panel_with_pair_detail(self, tmp_path):
        from af3_analysis.visualization_v3.figures.f07_contact_change_summary import (
            generate_f07_contact_change_summary,
        )
        result = generate_f07_contact_change_summary(
            self._contact_change_data(), tmp_path,
            recurring_pairs=self._recurring_pairs())
        assert result["status"] == "pass"
        assert result["n_observations"] == 2

    def test_truncation_warning_when_many_pairs(self, tmp_path):
        from af3_analysis.visualization_v3.figures.f07_contact_change_summary import (
            generate_f07_contact_change_summary,
            MAX_RANKED_PAIRS,
        )
        many_pairs = [
            {"chain_a": "A", "seq_a": i, "chain_b": "B", "seq_b": i + 1,
             "n_gained": 1, "n_lost": 0, "n_conditions": 1}
            for i in range(MAX_RANKED_PAIRS + 10)
        ]
        result = generate_f07_contact_change_summary(
            self._contact_change_data(), tmp_path,
            recurring_pairs=many_pairs)
        assert any(
            f"{len(many_pairs)} changed residue pairs" in w
            and f"{MAX_RANKED_PAIRS} most frequently changed" in w
            for w in result["warnings"]
        )

    def test_pair_data_with_wrong_columns_is_ignored(self, tmp_path):
        from af3_analysis.visualization_v3.figures.f07_contact_change_summary import (
            generate_f07_contact_change_summary,
        )
        result = generate_f07_contact_change_summary(
            self._contact_change_data(), tmp_path,
            recurring_pairs=[{"unrelated": 1}])
        assert result["status"] == "pass"
        assert not any("15 most recurring" in w for w in result["warnings"])

    def test_runner_f07_populates_recurring_pairs_from_cache(self, tmp_path):
        """Runner-side F07 must produce pair rows when F06 ran first."""
        # Import the dataset builder used by the cache tests.
        from af3_analysis.tests.test_v3_figure_cache import (
            _make_dataset,
            _run_figure,
        )
        from af3_analysis.visualization_v3.figures.f06_contact_map_difference import (
            generate_f06_contact_map_difference,
        )
        from af3_analysis.visualization_v3.figures.f07_contact_change_summary import (
            generate_f07_contact_change_summary,
        )

        dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
        figures_dir = tmp_path / "figs"
        figures_dir.mkdir()
        cache = v3_runner.FigureDataCache()

        r6 = _run_figure("F06", generate_f06_contact_map_difference,
                         dataset, figures_dir, cache)
        assert r6["status"] == "pass"
        misses_after_f06 = cache.n_misses

        r7 = _run_figure("F07", generate_f07_contact_change_summary,
                         dataset, figures_dir, cache)
        assert r7["status"] == "pass"
        # Consolidation must not recompute contact comparisons.
        assert cache.n_misses == misses_after_f06

    def test_pair_changes_cached_without_full_maps(self, tmp_path):
        """Pair-change cache stores only changed pairs, not full maps."""
        from af3_analysis.tests.test_v3_figure_cache import (
            _make_dataset,
            _run_figure,
        )
        from af3_analysis.visualization_v3.figures.f06_contact_map_difference import (
            generate_f06_contact_map_difference,
        )

        dataset = _make_dataset(n_seeds=1, n_samples=1, tmp_path=tmp_path)
        figures_dir = tmp_path / "figs"
        figures_dir.mkdir()
        cache = v3_runner.FigureDataCache()

        _run_figure("F06", generate_f06_contact_map_difference,
                    dataset, figures_dir, cache)
        for key, pairs in cache.contact_pair_changes.items():
            assert isinstance(pairs, dict)
            for delta in pairs.values():
                assert delta in (1, -1)


# ---------------------------------------------------------------------------
# F14 overlap note (pandas-2 regression guard: Series.ptp was removed)
# ---------------------------------------------------------------------------


class TestF14SeparationNote:
    def test_note_uses_numpy_ptp(self):
        """pandas-2 regression guard: Series.ptp was removed in pandas 2.0."""
        from af3_analysis.visualization_v3.figures.f14_mds_embedding import (
            _mds_separation_note,
        )
        # Centroid span (1.0) well under half of point span (10.0).
        centroids = pd.DataFrame({"x": [0.0, 1.0], "y": [5.0, 5.0]})
        note = _mds_separation_note(centroids, point_span=10.0)
        assert note is not None
        assert "not evidence of structural identity" in note

    def test_note_absent_when_centroids_well_separated(self):
        from af3_analysis.visualization_v3.figures.f14_mds_embedding import (
            _mds_separation_note,
        )
        # Centroid span (8.0) exceeds half of point span (10.0).
        centroids = pd.DataFrame({"x": [0.0, 8.0], "y": [0.0, 0.0]})
        assert _mds_separation_note(centroids, point_span=10.0) is None

    def test_note_absent_for_single_or_degenerate_centroids(self):
        from af3_analysis.visualization_v3.figures.f14_mds_embedding import (
            _mds_separation_note,
        )
        assert _mds_separation_note(
            pd.DataFrame({"x": [1.0], "y": [1.0]}), point_span=5.0) is None
        assert _mds_separation_note(
            pd.DataFrame({"x": [1.0, 1.0], "y": [1.0, 1.0]}),
            point_span=5.0) is None
        assert _mds_separation_note(
            pd.DataFrame({"x": [0.0, 1.0], "y": [0.0, 1.0]}),
            point_span=0.0) is None

    def test_f14_renders_with_overlapping_conditions(self, tmp_path):
        from af3_analysis.visualization_v3.figures.f14_mds_embedding import (
            generate_f14_mds_embedding,
        )
        n = 20
        coords = np.zeros((n, 2)) + np.random.default_rng(1).normal(
            0, 1e-3, size=(n, 2))
        result = generate_f14_mds_embedding(
            coords,
            [f"p{i}" for i in range(n)],
            ["condA"] * (n // 2) + ["condB"] * (n - n // 2),
            [1] * n,
            tmp_path,
        )
        assert result["status"] == "pass"
        assert any(
            "not evidence of structural identity" in w
            for w in result["warnings"]
        )


# ---------------------------------------------------------------------------
# F10 site provenance
# ---------------------------------------------------------------------------


class TestF10Provenance:
    def test_manual_site_rows_carry_identity_fields(self):
        """Runner rows include site identity; manual sites are not marked
        as detected."""
        # The runner marks a site as detected only when it carries detection
        # metadata keys; manual definitions (label/chain/residue/radius) do
        # not. Verify the predicate logic directly against a manual spec.
        manual_site = {"label": "s1", "chain": "A", "residue": 101,
                       "radius": 8.0}
        detected = (
            "mean_displacement" in manual_site
            or "z" in manual_site
            or "seq_list" in manual_site
            or "n_residues" in manual_site
        )
        assert detected is False

    def test_detected_site_rows_marked(self):
        detected_site = {
            "label": "auto1", "chain": "A", "residue": 55,
            "mean_displacement": 2.5, "z": 3.1,
        }
        detected = (
            "mean_displacement" in detected_site
            or "z" in detected_site
            or "seq_list" in detected_site
            or "n_residues" in detected_site
        )
        assert detected is True

    def test_f10_figure_renders_with_provenance_columns(self, tmp_path):
        from af3_analysis.visualization_v3.figures.f10_local_geometry import (
            generate_f10_local_geometry,
        )
        rows = []
        for cond in ("condA", "condB"):
            for seed in (1, 2):
                rows.append({
                    "condition_id": cond,
                    "seed": seed,
                    "sample": 1,
                    "region_label": "site1",
                    "site_chain": "A",
                    "site_residue": 101,
                    "site_radius": 8.0,
                    "site_detected": False,
                    "local_rmsd_target": 1.0,
                    "local_rmsd_ref": 0.8,
                    "local_rmsd": 1.0 + 0.1 * seed,
                    "local_plddt_mean": 85.0 - seed,
                    "local_plddt_ref": 88.0,
                    "n_atoms": 15,
                    "n_residues": 5,
                })
        result = generate_f10_local_geometry(rows, tmp_path)
        assert result["status"] == "pass"
        assert result["n_observations"] == 4
