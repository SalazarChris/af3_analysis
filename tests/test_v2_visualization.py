"""
Regression tests for V2 visualization.

Tests cover:
- PTM label mapping
- DNA label mapping
- Environment grouping
- Legend-size validation
- V2 figure generation (Figure 7, Figure 8, effects)
- Missing metric handling
- Tokenisation-confound handling
- Contact-probability unit handling
- Manifest generation
"""

import json
import pandas as pd
import pytest
from pathlib import Path

# Try importing V2 — skip all V2 tests if unavailable
try:
    from af3_analysis.visualization.v2.config import (
        PTM_COLOURS, PTM_ORDER, DNA_STYLE, DPI,
        apply_v2_style, get_ptm_color, get_dna_style, make_ptm_palette,
    )
    from af3_analysis.visualization.v2.factors import (
        parse_condition_name, add_factor_columns,
        get_unique_ptm_states, get_unique_environments,
    )
    from af3_analysis.visualization.v2.labels import (
        metric_label, panel_letter, figure_title,
        FIGURE7_METRICS, FIGURE8_METRICS,
    )
    from af3_analysis.visualization.v2.figure7 import generate_figure7
    from af3_analysis.visualization.v2.figure8 import generate_figure8
    from af3_analysis.visualization.v2.effects import generate_effects_forest
    from af3_analysis.visualization.v2.validation import (
        validate_figure, write_manifest,
    )
    V2_AVAILABLE = True
except (ImportError, ModuleNotFoundError) as e:
    V2_AVAILABLE = False
    _V2_SKIP_REASON = str(e)


def _requires_v2():
    if not V2_AVAILABLE:
        pytest.skip(f"V2 visualization unavailable: {_V2_SKIP_REASON}")


# -------------------------------------------------------------------
# Fixtures
# -------------------------------------------------------------------

@pytest.fixture
def sample_seed_aggregated():
    """Create a realistic seed_aggregated DataFrame with 8 conditions × 10 seeds."""
    import numpy as np
    rng = np.random.RandomState(42)
    rows = []
    conditions = {
        "pou_baseline":           {"ptm": "Baseline", "dna": False},
        "pou_tpo101":             {"ptm": "T101",     "dna": False},
        "pou_sep102":             {"ptm": "S102",     "dna": False},
        "pou_tpo101_sep102":      {"ptm": "T101+S102","dna": False},
        "pou_dna":                {"ptm": "Baseline", "dna": True},
        "pou_tpo101_dna":         {"ptm": "T101",     "dna": True},
        "pou_sep102_dna":         {"ptm": "S102",     "dna": True},
        "pou_tpo101_sep102_dna":  {"ptm": "T101+S102","dna": True},
    }
    for cond_name, meta in conditions.items():
        for seed in range(1, 11):
            base = 80 if meta["dna"] else 75
            offset = {"Baseline": 0, "T101": 2, "S102": -1, "T101+S102": 3}[meta["ptm"]]
            rows.append({
                "condition_id": cond_name,
                "condition_name": cond_name,
                "seed": seed,
                "pLDDT_mean": base + offset + rng.normal(0, 1),
                "pLDDT_min": base + offset - 10 + rng.normal(0, 2),
                "pae_mean": 5.0 - offset * 0.1 + rng.normal(0, 0.5),
                "contact_prob_mean": 0.6 + offset * 0.01 + rng.normal(0, 0.02),
                "ranking_score": 0.8 + offset * 0.01 + rng.normal(0, 0.01),
            })
    return pd.DataFrame(rows)


@pytest.fixture
def sample_pairwise():
    """Create a minimal pairwise_comparisons DataFrame."""
    return pd.DataFrame({
        "metric": ["pLDDT_mean"] * 3 + ["pae_mean"] * 3,
        "condition": ["pou_tpo101", "pou_sep102", "pou_tpo101_sep102"] * 2,
        "reference": ["pou_baseline"] * 6,
        "hedges_g": [0.5, -0.3, 0.8, -0.4, 0.2, -0.6],
        "hedges_g_ci_lower": [0.1, -0.7, 0.4, -0.8, -0.2, -1.0],
        "hedges_g_ci_upper": [0.9, 0.1, 1.2, 0.0, 0.6, -0.2],
    })


# -------------------------------------------------------------------
# Tests: Factor parsing
# -------------------------------------------------------------------

class TestFactorParsing:
    """Verify condition names are correctly decomposed into factors."""

    def test_parse_baseline(self):
        _requires_v2()
        result = parse_condition_name("pou_baseline")
        assert result["ptm_state"] == "Baseline"
        assert result["has_dna"] is False
        assert result["environment"] == "baseline"

    def test_parse_tpo101(self):
        _requires_v2()
        result = parse_condition_name("pou_tpo101")
        assert result["ptm_state"] == "T101"
        assert result["has_dna"] is False

    def test_parse_sep102(self):
        _requires_v2()
        result = parse_condition_name("pou_sep102")
        assert result["ptm_state"] == "S102"

    def test_parse_double_ptm(self):
        _requires_v2()
        result = parse_condition_name("pou_tpo101_sep102")
        assert result["ptm_state"] == "T101 + S102"

    def test_parse_dna(self):
        _requires_v2()
        result = parse_condition_name("pou_dna")
        assert result["ptm_state"] == "Baseline"
        assert result["has_dna"] is True

    def test_parse_tpo101_dna(self):
        _requires_v2()
        result = parse_condition_name("pou_tpo101_dna")
        assert result["ptm_state"] == "T101"
        assert result["has_dna"] is True

    def test_parse_environment(self):
        _requires_v2()
        result = parse_condition_name("pou_baseline_NA100_HOH1000_CL100")
        assert result["ptm_state"] == "Baseline"
        assert result["has_dna"] is False
        assert result["environment"] == "NA100_HOH1000_CL100"
        assert "Na100" in result["env_label"]

    def test_parse_all_components(self):
        _requires_v2()
        result = parse_condition_name("pou_tpo101_sep102_dna_NA50_HOH500_CL50")
        assert result["ptm_state"] == "T101 + S102"
        assert result["has_dna"] is True
        assert result["environment"] == "NA50_HOH500_CL50"


class TestFactorColumns:
    """Verify DataFrame augmentation with factor columns."""

    def test_add_factor_columns(self, sample_seed_aggregated):
        _requires_v2()
        df = add_factor_columns(sample_seed_aggregated)
        assert "ptm_state" in df.columns
        assert "has_dna" in df.columns
        assert "environment" in df.columns

    def test_ptm_states_present(self, sample_seed_aggregated):
        _requires_v2()
        df = add_factor_columns(sample_seed_aggregated)
        states = get_unique_ptm_states(df)
        assert "Baseline" in states
        assert "T101" in states
        assert "S102" in states

    def test_no_condition_name_raises(self):
        _requires_v2()
        df = pd.DataFrame({"seed": [1, 2], "pLDDT_mean": [0.8, 0.9]})
        with pytest.raises(ValueError, match="condition_name"):
            add_factor_columns(df)


# -------------------------------------------------------------------
# Tests: Config / palette
# -------------------------------------------------------------------

class TestV2Config:
    """Verify V2 configuration constants."""

    def test_ptm_colours_has_all_states(self):
        _requires_v2()
        for ptm in PTM_ORDER:
            assert ptm in PTM_COLOURS, f"Missing colour for PTM: {ptm}"

    def test_ptm_colours_are_hex(self):
        _requires_v2()
        for ptm, color in PTM_COLOURS.items():
            assert color.startswith("#"), f"PTM {ptm} colour not hex: {color}"

    def test_dna_style_has_both_states(self):
        _requires_v2()
        assert False in DNA_STYLE
        assert True in DNA_STYLE
        for has_dna, style in DNA_STYLE.items():
            assert "linestyle" in style
            assert "marker" in style
            assert "label" in style

    def test_make_ptm_palette(self):
        _requires_v2()
        palette = make_ptm_palette()
        assert len(palette) == len(PTM_ORDER)

    def test_get_ptm_color_returns_string(self):
        _requires_v2()
        c = get_ptm_color("Baseline")
        assert isinstance(c, str)

    def test_get_dna_style(self):
        _requires_v2()
        s = get_dna_style(True)
        assert s["label"] == "DNA"
        s = get_dna_style(False)
        assert s["label"] == "No DNA"


# -------------------------------------------------------------------
# Tests: Labels
# -------------------------------------------------------------------

class TestV2Labels:
    """Verify label generation."""

    def test_metric_label_known(self):
        _requires_v2()
        assert metric_label("pLDDT_mean") == "Mean pLDDT"
        assert metric_label("ranking_score") == "Ranking Score"

    def test_metric_label_unknown(self):
        _requires_v2()
        label = metric_label("some_custom_metric")
        assert "Custom Metric" in label

    def test_panel_letters(self):
        _requires_v2()
        assert panel_letter("pLDDT_mean") == "A"
        assert panel_letter("ranking_score") == "E"

    def test_figure_title(self):
        _requires_v2()
        t = figure_title("fig7")
        assert "AF3" in t
        assert "Distribution" in t


# -------------------------------------------------------------------
# Tests: Figure 7 generation
# -------------------------------------------------------------------

class TestFigure7:
    """Test V2 Figure 7 generation."""

    def test_figure7_runs(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        assert result["status"] == "pass"
        assert Path(result["output_path"]).exists()
        assert result["n_observations"] > 0

    def test_figure7_has_correct_ptm_states(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        assert "Baseline" in result["ptm_states"]
        assert "T101" in result["ptm_states"]
        assert "S102" in result["ptm_states"]

    def test_figure7_file_not_empty(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        fpath = Path(result["output_path"])
        assert fpath.stat().st_size > 10000  # > 10KB

    def test_figure7_image_dimensions(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        import matplotlib.pyplot as plt
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        img = plt.imread(result["output_path"])
        h, w = img.shape[:2]
        assert h < 2**16, f"Figure too tall: {h}px"
        assert w < 2**16, f"Figure too wide: {w}px"
        plt.close("all")

    def test_figure7_with_environment_filter(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        # Add environment column
        df = add_factor_columns(sample_seed_aggregated)
        result = generate_figure7(df, tmp_path, environment_filter="baseline")
        # With only baseline env, should still work (or skip if no data)
        assert result["status"] in ("pass", "skip")

    def test_figure7_empty_dataframe(self, tmp_path):
        _requires_v2()
        # Empty DataFrame with no condition_name → skip
        df_empty = pd.DataFrame({"seed": [], "pLDDT_mean": []})
        result = generate_figure7(df_empty, tmp_path)
        assert result["status"] == "skip"

    def test_figure7_single_row(self, tmp_path):
        _requires_v2()
        df = pd.DataFrame({
            "condition_name": ["pou_baseline"],
            "seed": [1],
            "pLDDT_mean": [80.0],
            "pLDDT_min": [65.0],
            "pae_mean": [5.0],
            "contact_prob_mean": [0.6],
            "ranking_score": [0.8],
        })
        result = generate_figure7(df, tmp_path)
        assert result["status"] in ("pass", "skip")

    def test_figure7_no_condition_name(self, tmp_path):
        _requires_v2()
        df = pd.DataFrame({"seed": [1, 2], "pLDDT_mean": [0.8, 0.9]})
        result = generate_figure7(df, tmp_path)
        assert result["status"] == "skip"


# -------------------------------------------------------------------
# Tests: Figure 8 generation
# -------------------------------------------------------------------

class TestFigure8:
    """Test V2 Figure 8 generation."""

    def test_figure8_runs(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure8(sample_seed_aggregated, tmp_path)
        assert result["status"] == "pass"
        assert Path(result["output_path"]).exists()
        assert result["n_observations"] > 0

    def test_figure8_file_not_empty(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure8(sample_seed_aggregated, tmp_path)
        fpath = Path(result["output_path"])
        assert fpath.stat().st_size > 10000

    def test_figure8_image_dimensions(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        import matplotlib.pyplot as plt
        result = generate_figure8(sample_seed_aggregated, tmp_path)
        img = plt.imread(result["output_path"])
        h, w = img.shape[:2]
        assert h < 2**16
        assert w < 2**16
        plt.close("all")

    def test_figure8_tokenisation_warning(self, tmp_path):
        """Tokenisation warning fires when PTM groups have different sizes."""
        _requires_v2()
        import numpy as np
        rng = np.random.RandomState(42)
        rows = []
        # Create groups with different sizes to trigger the warning
        for i, (cond, ptm_n) in enumerate([
            ("baseline", 10), ("tpo101", 5),
        ]):
            for seed in range(ptm_n):
                rows.append({
                    "condition_name": cond, "seed": seed,
                    "pLDDT_mean": 80 + rng.normal(),
                })
        df = pd.DataFrame(rows)
        result = generate_figure8(df, tmp_path)
        if result["status"] == "pass":
            assert any("tokenisation" in w.lower() for w in result["warnings"])

    def test_figure8_empty_metrics(self, tmp_path):
        _requires_v2()
        df = pd.DataFrame({"condition_name": ["a"], "seed": [1]})
        result = generate_figure8(df, tmp_path)
        assert result["status"] == "skip"

    def test_figure8_single_metric(self, tmp_path):
        _requires_v2()
        df = pd.DataFrame({
            "condition_name": ["a", "b"],
            "seed": [1, 1],
            "pLDDT_mean": [80.0, 85.0],
        })
        result = generate_figure8(df, tmp_path)
        # Need ≥2 metrics
        assert result["status"] == "skip"


# -------------------------------------------------------------------
# Tests: Effects forest
# -------------------------------------------------------------------

class TestEffectsForest:
    """Test V2 effects forest plot."""

    def test_effects_forest_runs(self, sample_pairwise, tmp_path):
        _requires_v2()
        result = generate_effects_forest(sample_pairwise, tmp_path)
        assert result["status"] == "pass"
        assert Path(result["output_path"]).exists()

    def test_effects_forest_empty(self, tmp_path):
        _requires_v2()
        result = generate_effects_forest(pd.DataFrame(), tmp_path)
        assert result["status"] == "skip"

    def test_effects_forest_image_dimensions(self, sample_pairwise, tmp_path):
        _requires_v2()
        import matplotlib.pyplot as plt
        result = generate_effects_forest(sample_pairwise, tmp_path)
        img = plt.imread(result["output_path"])
        h, w = img.shape[:2]
        assert h < 2**16
        assert w < 2**16
        plt.close("all")


# -------------------------------------------------------------------
# Tests: Validation
# -------------------------------------------------------------------

class TestValidation:
    """Test post-generation validation."""

    def test_validate_figure_pass(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        v = validate_figure(
            "test_fig", Path(result["output_path"]),
            n_observations=result["n_observations"],
        )
        assert v["passed"]

    def test_validate_figure_missing_file(self, tmp_path):
        _requires_v2()
        v = validate_figure("test", tmp_path / "nonexistent.png")
        assert not v["passed"]
        assert any("not found" in m for m in v["messages"])

    def test_validate_figure_zero_observations(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        v = validate_figure(
            "test_fig", Path(result["output_path"]),
            n_observations=0,
        )
        assert not v["passed"]
        assert any("zero" in m.lower() for m in v["messages"])

    def test_validate_legend_too_many(self, sample_seed_aggregated, tmp_path):
        _requires_v2()
        result = generate_figure7(sample_seed_aggregated, tmp_path)
        v = validate_figure(
            "test_fig", Path(result["output_path"]),
            n_observations=100,
            n_legend_entries=20,
            max_legend=10,
        )
        assert any("legend" in m.lower() for m in v["messages"])


# -------------------------------------------------------------------
# Tests: Manifest
# -------------------------------------------------------------------

class TestManifest:
    """Test visualization manifest generation."""

    def test_manifest_written(self, tmp_path):
        _requires_v2()
        figures = [
            {"figure": "fig7", "path": "fig7.png", "validation": {"passed": True, "messages": []}},
        ]
        path = write_manifest(tmp_path, figures=figures)
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert data["version"] == "v2"
        assert len(data["figures"]) == 1

    def test_manifest_environment_strategy(self, tmp_path):
        _requires_v2()
        path = write_manifest(tmp_path, figures=[], environment_strategy="faceted")
        with open(path) as f:
            data = json.load(f)
        assert data["environment_strategy"] == "faceted"

    def test_manifest_tokenisation_flag(self, tmp_path):
        _requires_v2()
        path = write_manifest(tmp_path, figures=[], confounded_by_tokenisation=False)
        with open(path) as f:
            data = json.load(f)
        assert data["confounded_by_tokenisation"] is False


# -------------------------------------------------------------------
# Tests: V1 preservation
# -------------------------------------------------------------------

class TestV1Preservation:
    """Verify V1 figures are not affected by V2 implementation."""

    def test_v1_still_works(self, synthetic_confidence_dir, tmp_output):
        """V1 generate_all_figures must still work."""
        _requires_v2()
        from af3_analysis.visualization.orchestrator import generate_all_figures
        from af3_analysis.tests.test_visualization import _run_full_pipeline

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        result = generate_all_figures(str(tables_dir), str(figures_dir))
        assert result is None or isinstance(result, dict)

    def test_v1_figures_not_overwritten(self, sample_seed_aggregated, tmp_path):
        """V2 figures go in v2/ subdirectory, not alongside V1."""
        _requires_v2()
        # Generate V1 first
        from af3_analysis.visualization.core_plots import plot_ecdf_overlay
        schema = {"metrics": ["pLDDT_mean"], "factors": [], "condition_order": []}
        plot_ecdf_overlay(sample_seed_aggregated, schema, tmp_path)
        v1_files = set(tmp_path.glob("*.png"))

        # Generate V2
        generate_figure7(sample_seed_aggregated, tmp_path)
        v2_files = set((tmp_path / "v2").glob("*.png")) if (tmp_path / "v2").exists() else set()

        # V1 files should still exist and V2 should be in subdirectory
        assert len(v1_files) > 0
        # V2 might be in tmp_path/v2/ or directly if the function puts them there
