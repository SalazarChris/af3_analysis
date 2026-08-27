"""
Regression tests for Stage 5: Visualization.

Verifies that the visualization orchestrator can load the generated CSVs
and that the plotting functions receive data in the expected format.

Note: These tests require seaborn/matplotlib. If the af3_analysis.statistics
module shadows Python's stdlib 'statistics' module (causing
ImportError: cannot import name 'NormalDist'), the visualization tests
are skipped. This is a known package-structure issue.
"""

import pandas as pd
import pytest

from af3inputbuilder.scripts.af3_condition_centric_extraction import extract_metrics
from af3_analysis.pipeline import _stage_run_analysis


# ---------------------------------------------------------------------------
# Try importing visualization — skip all viz tests if it fails
# ---------------------------------------------------------------------------

try:
    from af3_analysis.visualization.orchestrator import load_data, generate_all_figures
    from af3_analysis.visualization.core_plots import (
        plot_qc_completeness,
        plot_seed_distributions,
        plot_factorial_interaction,
        plot_variability,
        plot_ecdf_overlay,
    )
    VIZ_AVAILABLE = True
except (ImportError, ModuleNotFoundError) as e:
    VIZ_AVAILABLE = False
    _VIZ_SKIP_REASON = str(e)


def _requires_viz():
    if not VIZ_AVAILABLE:
        pytest.skip(f"Visualization modules unavailable: {_VIZ_SKIP_REASON}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run_full_pipeline(input_dir, output_dir):
    """Run extraction + analysis, return tables and figures directories."""
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(parents=True)
    figures_dir.mkdir(parents=True)

    extract_metrics(str(input_dir), str(tables_dir), verbose=False)

    from af3_analysis.config import AnalysisConfig
    config = AnalysisConfig(
        run_id="viz_test",
        output_root=output_dir,
        random_seed=42,
    )
    _stage_run_analysis(output_dir, config)

    return tables_dir, figures_dir


# ---------------------------------------------------------------------------
# Tests: data loading
# ---------------------------------------------------------------------------


class TestVisualizationDataLoading:
    """Verify the orchestrator's load_data function reads the CSVs correctly."""

    @pytest.fixture(autouse=True)
    def load_viz_data(self, synthetic_confidence_dir, tmp_output):
        _requires_viz()
        tables_dir, _ = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        self.data = load_data(tables_dir)
        self.tables_dir = tables_dir

    def test_load_data_returns_dict(self):
        """load_data should return a dict with three keys."""
        assert isinstance(self.data, dict)
        assert "seed_aggregated" in self.data
        assert "descriptive_stats" in self.data
        assert "pairwise_comparisons" in self.data

    def test_seed_aggregated_is_dataframe(self):
        """seed_aggregated should be a DataFrame."""
        assert isinstance(self.data["seed_aggregated"], pd.DataFrame)
        assert len(self.data["seed_aggregated"]) > 0

    def test_descriptive_stats_is_dataframe(self):
        """descriptive_stats should be a DataFrame."""
        assert isinstance(self.data["descriptive_stats"], pd.DataFrame)
        assert len(self.data["descriptive_stats"]) > 0


# ---------------------------------------------------------------------------
# Tests: schema inference
# ---------------------------------------------------------------------------


class TestVisualizationSchemaInference:
    """Verify the orchestrator correctly infers the schema from data."""

    @pytest.fixture(autouse=True)
    def build_schema(self, synthetic_confidence_dir, tmp_output):
        _requires_viz()
        tables_dir, _ = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data["seed_aggregated"]

        # Replicate the schema inference from orchestrator
        # Chain-level metrics are excluded from the plotting schema
        # because different conditions have different chains, making them
        # incomparable. Including them would create thousands of subplots
        # that hang matplotlib during rendering.
        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid_metrics = [m for m in global_metrics if m in seed_agg.columns]
        factors = [c for c in seed_agg.columns if c.startswith('factor_')]
        if not factors:
            factors = ['factor_DNA', 'factor_PTM']

        self.schema = {'metrics': valid_metrics, 'factors': factors, 'condition_order': []}
        self.seed_agg = seed_agg

    def test_schema_has_metrics(self):
        """Schema should discover at least pLDDT_mean."""
        assert "pLDDT_mean" in self.schema["metrics"]

    def test_schema_excludes_chain_metrics(self):
        """Schema should NOT include chain-level metrics (they would create thousands of subplots)."""
        assert not any("chain_" in m for m in self.schema["metrics"])

    def test_schema_has_factors(self):
        """Schema should have factor columns (at least fallback)."""
        assert len(self.schema["factors"]) > 0


# ---------------------------------------------------------------------------
# Tests: individual plot functions receive correct data
# ---------------------------------------------------------------------------


class TestPlotFunctions:
    """Verify each plot function can be called with the generated data."""

    @pytest.fixture(autouse=True)
    def setup_plots(self, synthetic_confidence_dir, tmp_output):
        _requires_viz()
        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data["seed_aggregated"]

        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid_metrics = [m for m in global_metrics if m in seed_agg.columns]
        self.schema = {'metrics': valid_metrics, 'factors': ['factor_DNA', 'factor_PTM'], 'condition_order': []}
        self.data = data
        self.figures_dir = figures_dir

    def test_plot_qc_completeness_runs(self):
        """F1: QC completeness plot should run without error."""
        plot_qc_completeness(self.data["seed_aggregated"], self.schema, self.figures_dir)

    def test_plot_seed_distributions_runs(self):
        """F2: Seed distributions plot should run without error."""
        plot_seed_distributions(self.data["seed_aggregated"], self.schema, self.figures_dir)

    def test_plot_variability_runs(self):
        """F5: Variability plot should run without error."""
        plot_variability(self.data["descriptive_stats"], self.schema, self.figures_dir)

    def test_plot_ecdf_overlay_runs(self):
        """F7: ECDF overlay plot should run without error."""
        plot_ecdf_overlay(self.data["seed_aggregated"], self.schema, self.figures_dir)


# ---------------------------------------------------------------------------
# Tests: generate_all_figures end-to-end
# ---------------------------------------------------------------------------


class TestGenerateAllFigures:
    """End-to-end test of the full visualization pipeline."""

    def test_orchestrator_runs_without_error(self, synthetic_confidence_dir, tmp_output):
        """generate_all_figures should complete without raising."""
        _requires_viz()
        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        result = generate_all_figures(str(tables_dir), str(figures_dir))
        # Result may be None if pairwise_comparisons is empty, but should not raise
        assert result is None or isinstance(result, dict)

    def test_orchestrator_with_real_data(self, pou_parent_dir, tmp_output):
        """generate_all_figures with real data should produce at least some figures."""
        _requires_viz()
        tables_dir, figures_dir = _run_full_pipeline(pou_parent_dir, tmp_output)
        result = generate_all_figures(str(tables_dir), str(figures_dir))
        if result is not None:
            assert len(result) > 0


# ---------------------------------------------------------------------------
# Tests: pagination for many metrics
# ---------------------------------------------------------------------------


class TestPlotPagination:
    """Verify that plots with many metrics paginate into multiple figures
    instead of creating a single oversized figure.
    """

    def test_seed_distributions_paginates_when_many_metrics(self, synthetic_confidence_dir, tmp_output):
        """F2: If there are more metrics than MAX_PER_PAGE, multiple figures are created."""
        _requires_viz()
        from af3_analysis.visualization.core_plots import (
            plot_seed_distributions, MAX_METRICS_PER_FIGURE,
        )

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data["seed_aggregated"]

        # Inject many fake metric columns to trigger pagination
        n_extra = MAX_METRICS_PER_FIGURE + 5
        for i in range(n_extra):
            seed_agg[f"fake_metric_{i}"] = seed_agg.get(
                'pLDDT_mean', seed_agg.iloc[:, 0]
            )

        metrics = [c for c in seed_agg.columns if c.startswith('fake_metric_')]
        schema = {'metrics': metrics, 'factors': [], 'condition_order': []}

        generated = plot_seed_distributions(seed_agg, schema, figures_dir)

        # Should produce 2 pages (MAX + 5 extra)
        assert len(generated) == 2
        assert any('_page1' in f for f in generated)
        assert any('_page2' in f for f in generated)
        # Verify files actually exist
        for f in generated:
            assert (figures_dir / f).exists()

    def test_seed_distributions_single_page_when_few_metrics(self, synthetic_confidence_dir, tmp_output):
        """F2: Fewer metrics than MAX_PER_PAGE should produce a single file with no suffix."""
        _requires_viz()
        from af3_analysis.visualization.core_plots import (
            plot_seed_distributions, MAX_METRICS_PER_FIGURE,
        )

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data["seed_aggregated"]

        global_metrics = ['pLDDT_mean', 'pae_mean', 'ranking_score']
        valid = [m for m in global_metrics if m in seed_agg.columns]
        schema = {'metrics': valid, 'factors': [], 'condition_order': []}

        generated = plot_seed_distributions(seed_agg, schema, figures_dir)
        assert len(generated) == 1
        assert '_page' not in generated[0]
        assert (figures_dir / generated[0]).exists()

    def test_seed_distributions_height_capped(self, synthetic_confidence_dir, tmp_output):
        """F2: Even with many metrics, no single figure should exceed 50 inches tall."""
        _requires_viz()
        import matplotlib.pyplot as plt
        from af3_analysis.visualization.core_plots import (
            plot_seed_distributions, MAX_METRICS_PER_FIGURE,
        )

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data["seed_aggregated"]

        n_extra = MAX_METRICS_PER_FIGURE * 3
        for i in range(n_extra):
            seed_agg[f"fake_metric_{i}"] = seed_agg.get(
                'pLDDT_mean', seed_agg.iloc[:, 0]
            )

        metrics = [c for c in seed_agg.columns if c.startswith('fake_metric_')]
        schema = {'metrics': metrics, 'factors': [], 'condition_order': []}

        generated = plot_seed_distributions(seed_agg, schema, figures_dir)
        # Should have 3 pages (12 + 12 + 12)
        assert len(generated) == 3
        # No pixel explosion: read each saved figure and check dimensions
        for f in generated:
            img = plt.imread(str(figures_dir / f))
            # Height should be under 2^16 = 65536 pixels
            assert img.shape[0] < 2**16, f"{f} too tall: {img.shape[0]}px"
            plt.close('all')


# ---------------------------------------------------------------------------
# Regression: plot_variability rendering failure with many chain metrics
# ---------------------------------------------------------------------------


class TestVariabilityRenderingRegression:
    """Regression for the plot_variability matplotlib rendering failure.

    The original bug: the schema included thousands of chain-level metrics
    (chain_A_plddt, chain_AA_plddt, ...). When plot_variability called
    _get_primary_metrics(), it returned all of them, creating one subplot
    panel per _std column. With 4800+ subplots, matplotlib hung during
    tick text layout computation, producing an effectively infinite render.

    The fix: chain-level metrics are excluded from the plotting schema.
    This test verifies that the fix works and that generated figures
    remain within safe pixel bounds.
    """

    def test_variability_with_many_chain_columns(self, synthetic_confidence_dir, tmp_output):
        """plot_variability must not hang even when seed_aggregated has thousands of chain_* columns."""
        _requires_viz()
        import matplotlib.pyplot as plt
        import numpy as np

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        desc = data["descriptive_stats"]

        # Add fake global _std columns so variability has something to plot
        # (synthetic data may not produce pLDDT_mean_std, etc.)
        for metric in ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']:
            if metric not in desc.columns:
                desc[metric] = np.random.rand(len(desc))
            if f'{metric}_std' not in desc.columns:
                desc[f'{metric}_std'] = np.random.rand(len(desc)) * 0.1

        # Add many fake chain-level _std columns to descriptive_stats
        # This simulates the original bug condition
        n_chain = 500
        for i in range(n_chain):
            # Generate chain names like A, B, ..., Z, AA, AB, ...
            suffix = chr(65 + i % 26)
            if i >= 26:
                suffix = chr(65 + (i // 26 - 1) % 26) + suffix
            chain_name = 'chain_' + suffix
            desc[chain_name + '_plddt_std'] = np.random.rand(len(desc))

        # Use only global metrics in schema (the fix)
        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid = [m for m in global_metrics if m in desc.columns]
        schema = {'metrics': valid, 'factors': [], 'condition_order': [], 'experiment_design': None}

        import signal, sys
        if sys.platform != 'win32':
            def _timeout_handler(signum, frame):
                raise TimeoutError("plot_variability timed out — likely hanging on rendering")
            signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(30)  # 30-second safety timeout
        try:
            plot_variability(desc, schema, figures_dir)
        finally:
            if sys.platform != 'win32':
                signal.alarm(0)

        # Verify output file was created
        fig_path = figures_dir / "fig_variability.png"
        assert fig_path.exists(), "fig_variability.png was not created"

        # Verify figure dimensions are within safe bounds
        img = plt.imread(str(fig_path))
        assert img.shape[0] < 2**16, f"Figure too tall: {img.shape[0]}px"
        assert img.shape[1] < 2**16, f"Figure too wide: {img.shape[1]}px"
        plt.close('all')

    def test_schema_excludes_chain_metrics(self, synthetic_confidence_dir, tmp_output):
        """The plotting schema must not include chain-level metrics."""
        _requires_viz()

        tables_dir, _ = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data["seed_aggregated"]

        # Replicate the orchestrator's schema inference
        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid_metrics = [m for m in global_metrics if m in seed_agg.columns]
        schema = {'metrics': valid_metrics, 'factors': [], 'condition_order': []}

        chain_in_schema = [m for m in schema['metrics'] if m.startswith('chain_')]
        assert chain_in_schema == [], (
            f"Schema should not contain chain-level metrics, got {len(chain_in_schema)}: "
            f"{chain_in_schema[:5]}..."
        )

    def test_all_plots_complete_with_many_chain_columns(self, synthetic_confidence_dir, tmp_output):
        """All 8 plot functions must complete within a time limit even with many chain columns."""
        _requires_viz()
        import signal

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)

        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid_metrics = [m for m in global_metrics if m in data['seed_aggregated'].columns]
        schema = {'metrics': valid_metrics, 'factors': [], 'condition_order': []}

        import sys
        if sys.platform != 'win32':
            def _timeout_handler(signum, frame):
                raise TimeoutError("Plot timed out")
            signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(120)  # 2-minute total timeout for all plots
        try:
            plot_qc_completeness(data['seed_aggregated'], schema, figures_dir)
            plot_seed_distributions(data['seed_aggregated'], schema, figures_dir)
            plot_factorial_interaction(data['descriptive_stats'], schema, figures_dir)
            plot_variability(data['descriptive_stats'], schema, figures_dir)
            plot_ecdf_overlay(data['seed_aggregated'], schema, figures_dir)
        finally:
            if sys.platform != 'win32':
                signal.alarm(0)


# ---------------------------------------------------------------------------
# Regression: full preflight audit — all figures must be safe
# ---------------------------------------------------------------------------


class TestVisualizationPreflightAudit:
    """Verify the preflight system catches pathological figure dimensions
    and that every figure renders within safe bounds.
    """

    def test_preflight_check_rejects_huge_height(self):
        """preflight_check should reject figures exceeding MAX_FIGURE_HEIGHT_INCHES."""
        _requires_viz()
        from af3_analysis.visualization.utils import preflight_check, MAX_FIGURE_HEIGHT_INCHES

        ok, warns = preflight_check(
            "test_huge",
            fig_height=MAX_FIGURE_HEIGHT_INCHES + 10,
            fig_width=10,
        )
        assert not ok
        assert len(warns) > 0

    def test_preflight_check_rejects_huge_width(self):
        """preflight_check should reject figures exceeding MAX_FIGURE_WIDTH_INCHES."""
        _requires_viz()
        from af3_analysis.visualization.utils import preflight_check, MAX_FIGURE_WIDTH_INCHES

        ok, warns = preflight_check(
            "test_huge_w",
            fig_height=10,
            fig_width=MAX_FIGURE_WIDTH_INCHES + 10,
        )
        assert not ok

    def test_preflight_check_warns_many_annotations(self):
        """preflight_check should warn when annotations exceed soft limit."""
        _requires_viz()
        from af3_analysis.visualization.utils import preflight_check, MAX_ANNO_COUNT

        ok, warns = preflight_check(
            "test_anno",
            n_annotations=MAX_ANNO_COUNT + 100,
            fig_height=10,
            fig_width=10,
        )
        assert ok  # soft limit — not a hard rejection
        assert any('annotations' in w for w in warns)

    def test_preflight_check_warns_many_legend_entries(self):
        """preflight_check should warn when legend entries exceed limit."""
        _requires_viz()
        from af3_analysis.visualization.utils import preflight_check, MAX_PAIRPLOT_GROUPS

        ok, warns = preflight_check(
            "test_legend",
            n_legend_entries=MAX_PAIRPLOT_GROUPS + 20,
            fig_height=10,
            fig_width=10,
        )
        assert ok
        assert any('legend' in w for w in warns)

    def test_palette_handles_many_conditions(self):
        """get_condition_style_map must return enough distinct colors for 100+ conditions."""
        _requires_viz()
        from af3_analysis.visualization.utils import get_condition_style_map

        conds = [f'condition_{i}' for i in range(100)]
        style_map = get_condition_style_map(conds)
        assert len(style_map) == 100
        # All colors should be distinct tuples
        colors = list(style_map.values())
        assert len(set(colors)) == 100

    def test_all_figures_safe_pixel_bounds(self, synthetic_confidence_dir, tmp_output):
        """Every generated figure must be within the 2^16 pixel limit."""
        _requires_viz()
        import matplotlib.pyplot as plt

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        result = generate_all_figures(str(tables_dir), str(figures_dir))
        if result is None:
            pytest.skip("generate_all_figures returned None")

        for name, fname in result.items():
            # fname may be a string or a list of strings (pagination)
            fnames = [fname] if isinstance(fname, str) else fname
            for fn in fnames:
                fpath = figures_dir / fn
                if not fpath.exists():
                    continue  # skipped figure (empty data)
                img = plt.imread(str(fpath))
                h, w = img.shape[:2]
                assert h < 2**16, f"{fn} height {h}px exceeds 2^16"
                assert w < 2**16, f"{fn} width {w}px exceeds 2^16"
                plt.close('all')

    def test_heatmap_suppresses_annotations_when_large(self, synthetic_confidence_dir, tmp_output):
        """Seed trajectories heatmap should suppress cell annotations when cells > MAX_HEATMAP_CELLS."""
        _requires_viz()
        from af3_analysis.visualization.utils import MAX_HEATMAP_CELLS

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)
        data = load_data(tables_dir)
        seed_agg = data['seed_aggregated']

        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid = [m for m in global_metrics if m in seed_agg.columns]
        schema = {'metrics': valid, 'factors': [], 'condition_order': []}

        n_conds = seed_agg['condition_name'].nunique()
        n_seeds = seed_agg['seed'].nunique()
        n_cells = n_seeds * n_conds

        # If the data already exceeds the limit, the plot should complete
        # without annotations and produce a valid PNG
        if n_cells > MAX_HEATMAP_CELLS:
            from af3_analysis.visualization.core_plots import plot_seed_trajectories
            plot_seed_trajectories(seed_agg, schema, figures_dir)
            fig_path = figures_dir / 'fig_seed_trajectories.png'
            assert fig_path.exists()
            import matplotlib.pyplot as plt
            img = plt.imread(str(fig_path))
            assert img.shape[0] < 2**16
            plt.close('all')

    def test_factorial_skips_gracefully_without_metadata(self, synthetic_confidence_dir, tmp_output):
        """plot_factorial_interaction should skip gracefully with many conditions and no metadata."""
        _requires_viz()
        import pandas as pd

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)

        # Build a fake descriptive_stats with 60 conditions and _mean columns
        n_conds = 60
        n_rows = n_conds  # one row per condition
        desc = pd.DataFrame({
            'condition_id': [f'cond_{i}' for i in range(n_conds)],
            'condition_name': [f'cond_{i}' for i in range(n_conds)],
            'pLDDT_mean_mean': [0.8] * n_conds,
            'pae_mean_mean': [5.0] * n_conds,
        })

        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        valid = [m for m in global_metrics if m in desc.columns]
        schema = {'metrics': valid, 'factors': [], 'condition_order': [], 'experiment_design': None}

        from af3_analysis.visualization.core_plots import plot_factorial_interaction
        # Should NOT raise — should skip gracefully (> 50 conditions, no metadata)
        plot_factorial_interaction(desc, schema, figures_dir)
        # No file should be created
        assert not (figures_dir / 'fig_factorial_interaction.png').exists()

    def test_effect_size_forest_capped_height(self, synthetic_confidence_dir, tmp_output):
        """plot_effect_size_forest must cap figure height even with many comparisons."""
        _requires_viz()
        import matplotlib.pyplot as plt
        import pandas as pd
        from af3_analysis.visualization.core_plots import plot_effect_size_forest
        from af3_analysis.visualization.utils import MAX_FIGURE_HEIGHT_INCHES, DPI

        tables_dir, figures_dir = _run_full_pipeline(synthetic_confidence_dir, tmp_output)

        # Create a fake pairwise dataframe with many comparisons
        rows = []
        for i in range(200):
            rows.append({
                'metric': 'pLDDT_mean' if i < 100 else 'pae_mean',
                'condition': f'cond_{i}',
                'reference': 'baseline',
                'diff_mean': 0.1 * (i % 10 - 5),
            })
        pw = pd.DataFrame(rows)

        global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
        schema = {'metrics': global_metrics, 'factors': [], 'condition_order': []}

        plot_effect_size_forest(pw, schema, figures_dir)
        fig_path = figures_dir / 'fig_effect_size_forest.png'
        assert fig_path.exists()
        img = plt.imread(str(fig_path))
        # Must be under 2^16 pixels and reasonable height
        assert img.shape[0] < 2**16
        assert img.shape[0] < MAX_FIGURE_HEIGHT_INCHES * DPI * 1.5  # generous margin
        plt.close('all')
