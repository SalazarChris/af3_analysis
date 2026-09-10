import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import warnings
from typing import Dict, Any, List, Tuple
import numpy as np

from .utils import (
    get_display_label, get_condition_order, get_condition_style_map,
    format_figure_layout, adjust_legend, preflight_check,
    MAX_FIGURE_HEIGHT_INCHES, MAX_FIGURE_WIDTH_INCHES, DPI,
    MAX_ANNO_COUNT, MAX_TICK_LABELS, MAX_HEATMAP_CELLS,
    MAX_PAIRPLOT_GROUPS, MAX_Y_LABELS_FOREST, MAX_ECDF_LINES,
)

warnings.filterwarnings("ignore", category=UserWarning)

# Maximum metrics per single figure to avoid exceeding matplotlib's
# 2^16 pixel limit.  At 100 DPI and ~5 inches per panel this keeps
# the figure well under the ceiling even for tall panels.
MAX_METRICS_PER_FIGURE = 12

def get_plot_path(filename: str, figures_dir: Path) -> Path:
    return figures_dir / filename


def _paginate_metrics(metrics: List[str], max_per_page: int = MAX_METRICS_PER_FIGURE) -> List[List[str]]:
    """Split *metrics* into pages of at most *max_per_page* items.

    Returns a list of lists.  Each inner list is one page.
    """
    return [metrics[i:i + max_per_page] for i in range(0, len(metrics), max_per_page)]


def _paginate_suffix(page_idx: int, n_pages: int) -> str:
    """Return filename suffix like ``""`` or ``"_page2"``.

    Page 1 of 1 returns ``""`` so single-page output filenames
    stay unchanged.
    """
    if n_pages <= 1:
        return ""
    return f"_page{page_idx + 1}"

def _get_primary_metrics(df: pd.DataFrame, schema: Dict[str, Any] = None) -> List[str]:
    if schema and 'metrics' in schema and schema['metrics']:
        def metric_exists(m):
            return (m in df.columns) or (f"{m}_mean" in df.columns) or (f"{m}_std" in df.columns)
        return [m for m in schema['metrics'] if metric_exists(m)]
        
    # Heuristic fallback if schema is missing/empty
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    valid_metrics = [c for c in numeric_cols if c not in ['seed', 'sample', 'replicate']]
    return [c for c in valid_metrics if 'mean' in c or 'plddt' in c.lower() or 'ptm' in c.lower() or 'score' in c.lower() or 'fraction' in c.lower() or 'clash' in c.lower()]

def plot_qc_completeness(manifest_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path):
    print("--- Executing plot_qc_completeness ---")
    metrics = _get_primary_metrics(manifest_df, schema)
    if not metrics: return
    
    heatmap_data = []
    if 'condition_name' not in manifest_df.columns:
        print("WARNING: condition_name not found. Skipping QC completeness plot.")
        return
        
    conditions = get_condition_order(manifest_df['condition_name'].unique().tolist(), schema)
    for cond in conditions:
        cond_df = manifest_df[manifest_df['condition_name'] == cond]
        row = {'condition': get_display_label(cond, 'condition', schema=schema)}
        for m in metrics:
            row[get_display_label(m, 'metric')] = cond_df[m].notna().mean()
        row['n'] = len(cond_df)
        heatmap_data.append(row)
    
    if not heatmap_data: return
    
    heatmap_df = pd.DataFrame(heatmap_data).set_index('condition')
    n_counts = heatmap_df.pop('n')
    
    n_conds = len(conditions)
    fig_height = max(4, n_conds * 0.6)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    
    # Preflight: suppress annotations if matrix is very large
    n_cells = len(heatmap_df) * len(heatmap_df.columns)
    ok, warns = preflight_check(
        "plot_qc_completeness",
        n_rows=len(heatmap_df), n_cols=len(heatmap_df.columns),
        n_axes=1, fig_width=11, fig_height=fig_height,
        n_annotations=n_cells, n_tick_labels_x=len(heatmap_df.columns),
        n_tick_labels_y=n_conds,
    )
    for w in warns: print(f"  PREFLIGHT WARNING: {w}")
    use_annot = n_cells <= MAX_HEATMAP_CELLS
    
    sns.heatmap(heatmap_df, annot=use_annot, cmap='viridis', vmin=0, vmax=1, ax=ax, cbar_kws={'label': 'Fraction Non-Null'})
    
    ax.set_ylabel("")
    ax.tick_params(axis='x', rotation=30)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
    
    # Annotate sample counts clearly outside heatmap
    for i, n in enumerate(n_counts):
        ax.text(len(heatmap_df.columns) + 0.15, i + 0.5, f"n={n}", va='center', fontweight='bold')
        
    format_figure_layout(fig, "Figure 1: QC Completeness Distribution (Fraction Non-Null)")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(get_plot_path("fig_qc_completeness.png", save_path), dpi=DPI)
    plt.close()

def plot_seed_distributions(seed_aggregated_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path) -> List[str]:
    """Box + strip plot of seed-level values per condition.

    Paginates across multiple figures when the number of metrics
    exceeds MAX_METRICS_PER_FIGURE.

    Returns a list of generated filenames (relative to *save_path*).
    """
    print("--- Executing plot_seed_distributions ---")
    metrics = _get_primary_metrics(seed_aggregated_df, schema)
    if not metrics:
        return []

    pages = _paginate_metrics(metrics)
    generated: List[str] = []

    conditions = get_condition_order(seed_aggregated_df['condition_name'].unique().tolist(), schema)
    style_map = get_condition_style_map(conditions, schema)

    # Display representation (shared across pages)
    df_plot = seed_aggregated_df.copy()
    df_plot['condition_display'] = df_plot['condition_name'].map(lambda x: get_display_label(x, 'condition', schema=schema))
    display_conditions = [get_display_label(c, 'condition', schema=schema) for c in conditions]
    display_style_map = {get_display_label(c, 'condition', schema=schema): color for c, color in style_map.items()}

    for page_idx, page_metrics in enumerate(pages):
        n = len(page_metrics)
        fig, axes = plt.subplots(n, 1, figsize=(11, 4.5 * n))
        if n == 1:
            axes = [axes]

        panel_offset = page_idx * MAX_METRICS_PER_FIGURE
        for i, (ax, metric) in enumerate(zip(axes, page_metrics)):
            sns.boxplot(
                data=df_plot,
                x='condition_display',
                y=metric,
                order=display_conditions,
                ax=ax,
                color='#f5f5f5',
                showfliers=False,
                width=0.5
            )
            sns.stripplot(
                data=df_plot,
                x='condition_display',
                y=metric,
                order=display_conditions,
                hue='condition_display',
                palette=display_style_map,
                ax=ax,
                alpha=0.6,
                jitter=0.25,
                size=5,
                legend=False
            )

            ax.set_title(f"Panel {chr(65 + panel_offset + i)}: {get_display_label(metric)}")
            ax.set_xlabel("")
            ax.set_ylabel(get_display_label(metric))
            ax.tick_params(axis='x', rotation=30)
            ax.set_xticklabels(display_conditions, rotation=30, ha='right')

        suffix = _paginate_suffix(page_idx, len(pages))
        fname = f"fig_seed_distributions{suffix}.png"
        title = "Figure 2: Seed Distributions"
        if len(pages) > 1:
            title += f" (page {page_idx + 1}/{len(pages)})"
        format_figure_layout(fig, title)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(get_plot_path(fname, save_path), dpi=DPI)
        plt.close()
        generated.append(fname)

    return generated

def plot_factorial_interaction(descriptive_stats_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path):
    """Plot factorial interaction using experiment design metadata.

    Adapts to the experimental design:
    * 2+ attributes → x-axis = first attribute, hue = second attribute
    * 1 attribute → x-axis = attribute, no hue
    * No metadata → falls back to alphabetical condition display
    """
    print("--- Executing plot_factorial_interaction ---")
    if 'condition_name' not in descriptive_stats_df.columns:
        print("WARNING: condition_name not found in descriptive_stats_df. Interaction plot skipped.")
        return
    
    df = descriptive_stats_df.copy()
    n_unique_conds = df['condition_name'].nunique()
    design = schema.get('experiment_design') if schema else None
    
    # Guard: without experiment design, too many conditions makes the
    # bar/line plot unreadable and wastes rendering time.
    if design is None and n_unique_conds > 50:
        print(f"  Skipping factorial interaction: {n_unique_conds} conditions "
              f"without experiment design metadata (max 50).")
        return
    
    # --- Build factor columns from metadata (or fall back) ---
    attribute_names: List[str] = []
    if design is not None and design.n_attributes > 0:
        # Add attribute columns to the DataFrame from metadata
        attribute_names = design.attribute_names
        for attr_name in attribute_names:
            col_name = f"attr_{attr_name}"
            df[col_name] = df['condition_name'].map(
                lambda c, a=attr_name: design.conditions[c].attributes.get(a)
                if c in design.conditions else None
            )
    else:
        # No metadata: just show conditions alphabetically (no factorial interpretation)
        print("  No experiment metadata available; showing conditions without factor structure.")
        attribute_names = []
    
    # Choose x-axis and hue from attributes
    if len(attribute_names) >= 2:
        x_attr = attribute_names[0]
        hue_attr = attribute_names[1]
    elif len(attribute_names) == 1:
        x_attr = attribute_names[0]
        hue_attr = None
    else:
        x_attr = None
        hue_attr = None
    
    # Filter outcome metrics to primary metrics (exclude chain-level —
    # they are not comparable across conditions with different chains)
    primary = _get_primary_metrics(df, schema)
    primary_global = [m for m in primary if not m.startswith('chain_')]
    metrics = [m for m in primary_global if f"{m}_mean" in df.columns]
    if not metrics:
        for c in df.columns:
            if c.endswith('_mean') and not c.endswith('_std') and not c.startswith('chain_'):
                metrics.append(c[:-5])
    
    if not metrics:
        print("  No metrics with _mean columns found. Skipping.")
        return
        
    n_panels = len(metrics)
    fig_height = min(4 * n_panels, 50)
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, fig_height))
    if n_panels == 1:
        axes = [axes]
    
    for i, (ax, metric) in enumerate(zip(axes, metrics)):
        mean_col = f"{metric}_mean"
        if mean_col not in df.columns:
            continue
            
        if x_attr is not None:
            # --- Factorial-style plot using metadata attributes ---
            x_col = f"attr_{x_attr}"
            
            # Format attribute values for display
            df_plot = df.copy()
            df_plot[x_col] = df_plot[x_col].map(
                lambda v, a=x_attr: _format_attr_value(v, a, design)
            )
            
            if hue_attr is not None:
                hue_col = f"attr_{hue_attr}"
                df_plot[hue_col] = df_plot[hue_col].map(
                    lambda v, a=hue_attr: _format_attr_value(v, a, design)
                )
                
                sns.lineplot(
                    data=df_plot, x=x_col, y=mean_col, hue=hue_col,
                    marker='o', ax=ax, linewidth=2.2, markersize=8
                )
                adjust_legend(ax, title=hue_attr)
            else:
                sns.barplot(
                    data=df_plot, x=x_col, y=mean_col,
                    ax=ax, errorbar=None
                )
                
            ax.set_xlabel(x_attr)
        else:
            # --- No metadata: alphabetical condition display ---
            cond_col = 'condition_name'
            df_plot = df.copy()
            df_plot['condition_display'] = df_plot[cond_col].map(
                lambda c: get_display_label(c, 'condition', schema=schema)
            )
            cond_order = get_condition_order(df[cond_col].unique().tolist(), schema)
            display_order = [get_display_label(c, 'condition', schema=schema) for c in cond_order]
            sns.barplot(
                data=df_plot, x='condition_display', y=mean_col,
                order=display_order, ax=ax, errorbar=None
            )
            ax.tick_params(axis='x', rotation=30)
            ax.set_xticklabels(display_order, rotation=30, ha='right')
            ax.set_xlabel("Condition")
        
        ax.set_title(f"Panel {chr(65+i)}: {get_display_label(metric)}")
        ax.set_ylabel(f"Mean {get_display_label(metric)}")
    
    title = "Figure 3: Factorial Interaction Analysis"
    if design is not None:
        from experiment_metadata import inspect_design as _inspect
        insp = _inspect(design)
        if insp.is_complete_factorial:
            level_counts = ' × '.join(
                str(len(insp.attribute_levels[a])) for a in insp.attribute_names
            )
            title += f" ({level_counts} design, {insp.n_conditions} conditions)"
        else:
            title += f" ({insp.n_conditions} conditions, incomplete)"
    
    format_figure_layout(fig, title)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(get_plot_path("fig_factorial_interaction.png", save_path), dpi=DPI)
    plt.close()


def _format_attr_value(value, attr_name: str, design=None) -> str:
    """Format an attribute value for display on plot axes/legends."""
    if value is None:
        return "?"
    if isinstance(value, bool):
        return f"No {attr_name}" if not value else attr_name
    # Categorical: use the value as-is, title-cased
    return str(value).replace('_', ' ').title()


def plot_effect_size_forest(pairwise_comparisons_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path):
    print("--- Executing plot_effect_size_forest ---")
    if len(pairwise_comparisons_df) == 0:
        return
        
    metrics = pairwise_comparisons_df['metric'].unique()
    n_metrics = len(metrics)
    
    val_col = 'hedges_g' if 'hedges_g' in pairwise_comparisons_df.columns else 'diff_mean'
    
    # Cap figure height: each metric panel gets at most MAX_Y_LABELS_FOREST rows
    rows_per_metric = []
    for metric in metrics:
        subset = pairwise_comparisons_df[pairwise_comparisons_df['metric'] == metric]
        rows_per_metric.append(min(len(subset), MAX_Y_LABELS_FOREST))
    total_rows = sum(rows_per_metric)
    fig_height = max(4, min(total_rows * 0.5, MAX_FIGURE_HEIGHT_INCHES))
    
    fig, axes = plt.subplots(n_metrics, 1, figsize=(11, fig_height))
    if n_metrics == 1: axes = [axes]
    
    for i, (ax, metric) in enumerate(zip(axes, metrics)):
        subset = pairwise_comparisons_df[pairwise_comparisons_df['metric'] == metric].copy()
        
        subset = subset.sort_values(by=val_col, ascending=True).reset_index(drop=True)
        subset['contrast'] = subset['condition'].apply(lambda x: get_display_label(x, 'condition', schema=schema)) + " vs " + subset['reference'].apply(lambda x: get_display_label(x, 'condition', schema=schema))
        
        if len(subset) > MAX_Y_LABELS_FOREST:
            subset = subset.head(MAX_Y_LABELS_FOREST)
            print(f"  WARNING: {metric} has {len(pairwise_comparisons_df[pairwise_comparisons_df['metric']==metric])} comparisons; "
                  f"showing top {MAX_Y_LABELS_FOREST}")
        
        if f'{val_col}_ci_lower' in subset.columns and f'{val_col}_ci_upper' in subset.columns:
            ax.errorbar(subset[val_col], range(len(subset)), 
                        xerr=[subset[val_col] - subset[f'{val_col}_ci_lower'], subset[f'{val_col}_ci_upper'] - subset[val_col]], 
                        fmt='o', color='black', capsize=4, elinewidth=1.5)
        else:
            ax.plot(subset[val_col], range(len(subset)), 'o', color='black', markersize=7)
            
        ax.set_yticks(range(len(subset)))
        ax.set_yticklabels(subset['contrast'])
        ax.axvline(0, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
        ax.set_title(f"Panel {chr(65+i)}: {get_display_label(metric)}")
        ax.set_xlabel(f"Effect Size ({get_display_label(val_col)})")
        
    format_figure_layout(fig, f"Figure 4: Effect Size Analysis ({get_display_label(val_col)})")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(get_plot_path("fig_effect_size_forest.png", save_path), dpi=DPI)
    plt.close()

def plot_variability(descriptive_stats_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path):
    print("--- Executing plot_variability ---")
    if 'condition_name' not in descriptive_stats_df.columns:
        print("WARNING: condition_name not found. Variability plot skipped.")
        return
        
    primary = _get_primary_metrics(descriptive_stats_df, schema)
    # Exclude chain-level metrics: they are not comparable across conditions
    # (different conditions have different chains) and would create thousands
    # of subplots, causing matplotlib to hang during rendering.
    primary_global = [m for m in primary if not m.startswith('chain_')]
    std_cols = [f"{m}_std" for m in primary_global if f"{m}_std" in descriptive_stats_df.columns]
    if not std_cols:
        print("  No global std columns found. Skipping variability plot.")
        return
    
    n_panels = len(std_cols)
    fig, axes = plt.subplots(n_panels, 1, figsize=(11, min(4 * n_panels, 50)))
    if n_panels == 1: axes = [axes]
    
    conditions = get_condition_order(descriptive_stats_df['condition_name'].unique().tolist(), schema)
    style_map = get_condition_style_map(conditions, schema)
    
    df_plot = descriptive_stats_df.copy()
    df_plot['condition_display'] = df_plot['condition_name'].map(lambda x: get_display_label(x, 'condition', schema=schema))
    display_conditions = [get_display_label(c, 'condition', schema=schema) for c in conditions]
    display_style_map = {get_display_label(c, 'condition', schema=schema): color for c, color in style_map.items()}
    
    for i, (ax, col) in enumerate(zip(axes, std_cols)):
        metric = col.replace('_std', '')
        sns.barplot(
            data=df_plot, 
            x='condition_display', 
            y=col, 
            order=display_conditions, 
            hue='condition_display', 
            palette=display_style_map, 
            ax=ax, 
            legend=False
        )
        ax.set_title(f"Panel {chr(65+i)}: {get_display_label(metric)} Variability")
        ax.set_xlabel("")
        ax.set_ylabel(f"SD ({get_display_label(metric)})")
        ax.tick_params(axis='x', rotation=30)
        ax.set_xticklabels(display_conditions, rotation=30, ha='right')
        
    format_figure_layout(fig, "Figure 5: Variability (between-seed sd)")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(get_plot_path("fig_variability.png", save_path), dpi=DPI)
    plt.close()

def plot_seed_trajectories(seed_aggregated_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path):
    """Plot seed × condition heatmap for each metric.

    Replaces the previous trajectory line plot. Lines between categorical
    conditions implied false continuity; a heatmap shows the same
    seed × condition pattern without that misleading visual encoding.

    Each panel is a heatmap: rows = seeds, columns = conditions,
    color = metric value.
    """
    print("--- Executing plot_seed_trajectories (heatmap) ---")
    metrics = _get_primary_metrics(seed_aggregated_df, schema)
    if not metrics:
        return

    if 'condition_name' not in seed_aggregated_df.columns or 'seed' not in seed_aggregated_df.columns:
        print("  WARNING: condition_name or seed column missing. Skipping.")
        return

    conditions = get_condition_order(
        seed_aggregated_df['condition_name'].unique().tolist(), schema
    )
    display_conditions = [
        get_display_label(c, 'condition', schema=schema) for c in conditions
    ]

    n_conditions = len(conditions)
    n_seeds = seed_aggregated_df['seed'].nunique()
    n_metrics = len(metrics)
    
    # Cap figure width/height
    fig_width = max(10, min(n_conditions * 1.4, MAX_FIGURE_WIDTH_INCHES))
    fig_height = min(3.5 * n_metrics, MAX_FIGURE_HEIGHT_INCHES)
    
    fig, axes = plt.subplots(n_metrics, 1, figsize=(fig_width, fig_height))
    if n_metrics == 1:
        axes = [axes]

    for i, (ax, metric) in enumerate(zip(axes, metrics)):
        # Pivot: rows = seed, columns = condition, values = metric
        pivot = seed_aggregated_df.pivot_table(
            index='seed', columns='condition_name', values=metric, aggfunc='first'
        )
        # Reorder columns to match condition order
        pivot = pivot[[c for c in conditions if c in pivot.columns]]
        n_cells = pivot.shape[0] * pivot.shape[1]
        
        # Preflight check
        ok, warns = preflight_check(
            f"plot_seed_trajectories[{metric}]",
            n_rows=pivot.shape[0], n_cols=pivot.shape[1],
            n_axes=1, fig_width=fig_width, fig_height=fig_height,
            n_heatmap_cells=n_cells,
            n_tick_labels_x=pivot.shape[1], n_tick_labels_y=pivot.shape[0],
        )
        for w in warns: print(f"  PREFLIGHT WARNING: {w}")
        
        # Suppress cell annotations for large heatmaps to prevent
        # rendering explosion (each annotation is a text artist that
        # matplotlib must lay out and render individually).
        use_annot = n_cells <= MAX_HEATMAP_CELLS
        if not use_annot:
            print(f"  INFO: suppressing annotations for {metric} "
                  f"({n_cells} cells > {MAX_HEATMAP_CELLS})")
        
        # Rename columns to display labels
        pivot.columns = [
            get_display_label(c, 'condition', schema=schema)
            for c in pivot.columns
        ]

        sns.heatmap(
            pivot,
            ax=ax,
            annot=use_annot,
            fmt='.2f',
            cmap='YlOrRd',
            linewidths=0.5,
            linecolor='white',
            cbar_kws={'label': get_display_label(metric), 'shrink': 0.8},
        )
        ax.set_title(
            f"Panel {chr(65+i)}: {get_display_label(metric)} — Seed × Condition"
        )
        ax.set_xlabel("")
        ax.set_ylabel("Seed")
        ax.tick_params(axis='x', rotation=30)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')

    format_figure_layout(fig, "Figure 6: Seed × Condition Heatmap")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(get_plot_path("fig_seed_trajectories.png", save_path), dpi=DPI)
    plt.close()

def plot_ecdf_overlay(seed_aggregated_df: pd.DataFrame, schema: Dict[str, Any], save_path: Path):
    print("--- Executing plot_ecdf_overlay ---")
    metrics = _get_primary_metrics(seed_aggregated_df, schema)
    if not metrics: return
    
    n_metrics = len(metrics)
    fig_height = min(4 * n_metrics, MAX_FIGURE_HEIGHT_INCHES)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(11, fig_height))
    if n_metrics == 1: axes = [axes]
    
    conditions = get_condition_order(seed_aggregated_df['condition_name'].unique().tolist(), schema)
    
    # Limit ECDF lines to avoid unreadable overlapping curves
    if len(conditions) > MAX_ECDF_LINES:
        print(f"  WARNING: {len(conditions)} conditions exceeds ECDF limit "
              f"({MAX_ECDF_LINES}). Showing all conditions but color cycling "
              f"will occur.")
    
    style_map = get_condition_style_map(conditions, schema)
    
    df_plot = seed_aggregated_df.copy()
    df_plot['condition_display'] = df_plot['condition_name'].map(lambda x: get_display_label(x, 'condition', schema=schema))
    
    display_conditions = [get_display_label(c, 'condition', schema=schema) for c in conditions]
    display_style_map = {get_display_label(c, 'condition', schema=schema): color for c, color in style_map.items()}
    
    for i, (ax, metric) in enumerate(zip(axes, metrics)):
        sns.ecdfplot(
            data=df_plot, 
            x=metric, 
            hue='condition_display', 
            hue_order=display_conditions, 
            palette=display_style_map, 
            ax=ax,
            linewidth=1.5
        )
        ax.set_title(f"Panel {chr(65+i)}: ECDF Overlay - {get_display_label(metric)}")
        ax.set_xlabel(get_display_label(metric))
        ax.set_ylabel("Proportion")
        # Move legend outside to avoid covering data
        adjust_legend(ax, title="Condition")
        
    format_figure_layout(fig, "Figure 7: ECDF Overlay")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(get_plot_path("fig_ecdf_overlay.png", save_path), dpi=DPI)
    plt.close()


def plot_metric_relationships(
    dataframe: pd.DataFrame,
    schema: Dict[str, Any],
    save_path: Path
):
    print("--- Executing plot_metric_relationships (Schema-Driven) ---")
    metrics = schema.get('metrics', [])
    if not metrics:
        return
        
    # Limit to at most 5 key metrics for clean, readable scatter matrix (F8)
    f8_metrics = ['pLDDT_mean', 'pae_mean', 'contact_prob_mean', 'ranking_score']
    # Note: chain-level metrics are intentionally excluded here because
    # different conditions have different chains → the pairplot would be
    # meaningless and the legend would have hundreds of entries.
        
    valid_metrics = [m for m in f8_metrics if m in dataframe.columns and pd.api.types.is_numeric_dtype(dataframe[m])]
    if len(valid_metrics) < 2:
        return
        
    n_vars = len(valid_metrics)
    n_pairs = n_vars * (n_vars - 1) // 2 if True else n_vars  # corner=True
    
    try:
        conditions = get_condition_order(dataframe['condition_name'].unique().tolist(), schema)
        style_map = get_condition_style_map(conditions, schema)
        
        df_plot = dataframe.copy()
        df_plot['condition_display'] = df_plot['condition_name'].map(lambda x: get_display_label(x, 'condition', schema=schema))
        rename_map = {m: get_display_label(m, 'metric') for m in valid_metrics}
        df_plot = df_plot.rename(columns=rename_map)
        
        valid_metrics_renamed = list(rename_map.values())
        display_conditions = [get_display_label(c, 'condition', schema=schema) for c in conditions]
        display_style_map = {get_display_label(c, 'condition', schema=schema): color for c, color in style_map.items()}
        
        n_conditions = len(display_conditions)
        
        # Preflight check for pairplot
        n_axes_est = n_vars * n_vars  # pairplot creates full grid
        ok, warns = preflight_check(
            "plot_metric_relationships",
            n_rows=len(df_plot), n_cols=n_vars,
            n_axes=n_axes_est,
            fig_width=n_vars * 2.2, fig_height=n_vars * 2.2,
            n_legend_entries=n_conditions,
        )
        for w in warns: print(f"  PREFLIGHT WARNING: {w}")
        
        fig = sns.pairplot(
            df_plot, 
            vars=valid_metrics_renamed, 
            hue='condition_display', 
            hue_order=display_conditions, 
            palette=display_style_map, 
            corner=True, 
            plot_kws={'alpha': 0.5, 's': 15},
            height=2.0
        )
        
        fig.figure.suptitle("Figure 8: Metric Relationship Matrix", y=0.98, fontsize=14, fontweight='bold')
        sns.move_legend(fig, "center right", bbox_to_anchor=(0.95, 0.5), title="Condition")
        
        plt.savefig(get_plot_path("fig_scatter_matrix.png", save_path), dpi=DPI)
        plt.close()
    except Exception as e:
        print(f"ERROR during F8 plot generation: {e}")
