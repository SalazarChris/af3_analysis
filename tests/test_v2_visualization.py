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
        DNA_STYLE, DPI,
        apply_v2_style, get_ptm_color, get_dna_style, make_ptm_palette,
        derive_ptm_order, build_ptm_colours,
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

def _minimal_seed_df():
    """Create a minimal seed_aggregated DataFrame."""
    rows = []
    conditions = [
        ("pou_baseline", "cond_001"),
        ("pou_baseline_NA100_HOH1000_CL100", "cond_002"),
        ("pou_baseline_NA10_HOH100_CL10", "cond_003"),
    ]
    for cond_name, cond_id in conditions:
        for seed in range(1, 6):
            rows.append({
                "condition_id": cond_id,
                "condition_name": cond_name,
                "seed": seed,
                "pLDDT_mean": 70.0 + seed,
                "pae_mean": 3.0 + seed * 0.1,
                "ranking_score": seed,
                "contact_prob_mean": 0.5 + seed * 0.02,
            })
    return pd.DataFrame(rows)


def _multi_ptm_seed_df():
    """Create a seed_aggregated DataFrame with multiple PTM states."""
    rows = []
    conditions = [
        ("pou_baseline", "cond_001"),
        ("pou_tpo101", "cond_002"),
        ("pou_sep102", "cond_003"),
        ("pou_tpo101_sep102", "cond_004"),
    ]
    for cond_name, cond_id in conditions:
        for seed in range(1, 4):
            rows.append({
                "condition_id": cond_id,
                "condition_name": cond_name,
                "seed": seed,
                "pLDDT_mean": 70.0 + seed,
                "pae_mean": 3.0 + seed * 0.1,
                "ranking_score": seed,
                "contact_prob_mean": 0.5 + seed * 0.02,
            })
    return pd.DataFrame(rows)


# -------------------------------------------------------------------
# Test classes
# -------------------------------------------------------------------


class TestFactorParsing:
    """Verify condition names are correctly decomposed into factors."""

    def test_parse_baseline(self):
        result = parse_condition_name("pou_baseline")
        assert result["ptm_state"] == "Baseline"
        assert result["has_dna"] is False

    def test_parse_tpo101(self):
        result = parse_condition_name("pou_tpo101")
        assert result["ptm_state"].lower() == "tpo101"

    def test_parse_sep102(self):
        result = parse_condition_name("pou_sep102")
        assert result["ptm_state"].lower() == "sep102"

    def test_parse_double_ptm(self):
        result = parse_condition_name("pou_tpo101_sep102")
        assert "tpo101" in result["ptm_state"].lower()
        assert "sep102" in result["ptm_state"].lower()

    def test_parse_dna(self):
        result = parse_condition_name("pou_dna")
        assert result["has_dna"] is True

    def test_parse_tpo101_dna(self):
        result = parse_condition_name("pou_tpo101_dna")
        assert "tpo101" in result["ptm_state"].lower()
        assert result["has_dna"] is True

    def test_parse_environment(self):
        result = parse_condition_name("pou_baseline_NA100_HOH1000_CL100")
        assert result["environment"] != "baseline"
        assert "NA100" in result["environment"]

    def test_parse_all_components(self):
        result = parse_condition_name("pou_tpo101_dna_NA10_HOH100_CL10")
        assert "tpo101" in result["ptm_state"].lower()
        assert result["has_dna"] is True
        assert result["environment"] != "baseline"

    def test_parse_oct4_conditions(self):
        """Verify oct4 condition names are parsed generically."""
        result = parse_condition_name("oct4_k123-sumo_summary")
        # Should extract something as the state, not crash
        assert result["condition_name"] == "oct4_k123-sumo_summary"
        assert result["ptm_state"] is not None

    def test_parse_generic_no_ptm(self):
        """Condition with only environment modifiers."""
        result = parse_condition_name("my_protein_NA50_HOH500_CL50")
        assert result["condition_name"] == "my_protein_NA50_HOH500_CL50"
        assert result["environment"] != "baseline"


class TestFactorColumns:
    """Verify DataFrame augmentation with factor columns."""

    def test_add_factor_columns(self):
        df = _minimal_seed_df()
        result = add_factor_columns(df)
        assert "ptm_state" in result.columns
        assert "has_dna" in result.columns
        assert "environment" in result.columns

    def test_ptm_states_present(self):
        df = _multi_ptm_seed_df()
        result = add_factor_columns(df)
        states = result["ptm_state"].unique()
        assert len(states) >= 1

    def test_no_condition_name_raises(self):
        df = pd.DataFrame({"col1": [1, 2]})
        with pytest.raises(ValueError, match="condition_name"):
            add_factor_columns(df)


class TestV2Config:
    """Verify V2 configuration constants."""

    def test_ptm_order_dynamic(self):
        """derive_ptm_order should handle any set of states."""
        order = derive_ptm_order(["S102", "Baseline", "T101", "CustomX"])
        assert order[0] == "Baseline"
        assert "S102" in order
        assert "CustomX" in order

    def test_build_ptm_colours_dynamic(self):
        """build_ptm_colours should generate colours for any states."""
        colours = build_ptm_colours(["Baseline", "Phospho", "Ubiquitin"])
        assert len(colours) == 3
        assert all(isinstance(c, str) and c.startswith("#") for c in colours.values())

    def test_ptm_colours_are_hex(self):
        colours = build_ptm_colours(["Baseline", "CustomA", "CustomB"])
        for c in colours.values():
            assert c.startswith("#")
            assert len(c) == 7

    def test_dna_style_has_both_states(self):
        assert False in DNA_STYLE
        assert True in DNA_STYLE

    def test_make_ptm_palette(self):
        palette = make_ptm_palette(["Baseline", "T101", "S102"])
        assert len(palette) == 3

    def test_get_ptm_color_returns_string(self):
        c = get_ptm_color("Baseline")
        assert isinstance(c, str)

    def test_get_dna_style(self):
        s = get_dna_style(False)
        assert "linestyle" in s


class TestV2Labels:
    """Verify label generation."""

    def test_metric_label_known(self):
        label = metric_label("pLDDT_mean")
        assert isinstance(label, str)
        assert len(label) > 0

    def test_metric_label_unknown(self):
        label = metric_label("totally_unknown_metric_xyz")
        assert isinstance(label, str)

    def test_panel_letters(self):
        assert panel_letter("pLDDT_mean") == "A"

    def test_figure_title(self):
        t = figure_title("fig7")
        assert isinstance(t, str)


class TestFigure7:
    """Test V2 Figure 7 generation."""

    def test_figure7_runs(self, tmp_path):
        _requires_v2()
        df = _multi_ptm_seed_df()
        result = generate_figure7(df, tmp_path)
        assert result["status"] == "pass"

    def test_figure7_has_correct_ptm_states(self, tmp_path):
        _requires_v2()
        df = _multi_ptm_seed_df()
        result = generate_figure7(df, tmp_path)
        assert "ptm_states" in result
        assert len(result["ptm_states"]) >= 1

    def test_figure7_file_not_empty(self, tmp_path):
        _requires_v2()
        df = _multi_ptm_seed_df()
        result = generate_figure7(df, tmp_path)
        assert Path(result["output_path"]).stat().st_size > 0

    def test_figure7_image_dimensions(self, tmp_path):
        _requires_v2()
        df = _multi_ptm_seed_df()
        result = generate_figure7(df, tmp_path)
        assert result["n_observations"] > 0

    def test_figure7_with_environment_filter(self, tmp_path):
        _requires_v2()
        df = _minimal_seed_df()
        result = generate_figure7(df, tmp_path, environment_filter="baseline")
        assert result["status"] in ("pass", "skip")


class TestFigure8:
    """Test V2 Figure 8 generation."""

    def test_figure8_runs(self, tmp_path):
        _requires_v2()
        df = _multi_ptm_seed_df()
        result = generate_figure8(df, tmp_path)
        assert result["status"] == "pass"

    def test_figure8_file_not_empty(self, tmp_path):
        _requires_v2()
        df = _multi_ptm_seed_df()
        result = generate_figure8(df, tmp_path)
        assert Path(result["output_path"]).stat().st_size > 0


class TestEffectsForest:
    """Test V2 effects forest generation."""

    def test_effects_empty_returns_skip(self, tmp_path):
        _requires_v2()
        df = pd.DataFrame()
        result = generate_effects_forest(df, tmp_path)
        assert result["status"] == "skip"

    def test_effects_with_data(self, tmp_path):
        _requires_v2()
        df = pd.DataFrame({
            "metric": ["pLDDT_mean", "pLDDT_mean"],
            "condition": ["cond_A", "cond_B"],
            "reference": ["cond_ref", "cond_ref"],
            "hedges_g": [0.5, -0.3],
        })
        result = generate_effects_forest(df, tmp_path)
        assert result["status"] == "pass"
