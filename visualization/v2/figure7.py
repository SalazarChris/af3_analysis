"""
V2 Figure 7 — Distributional Analysis of AF3 Confidence Metrics
================================================================

Five panels (A–E), each showing the distribution of one global metric
across PTM states, faceted by DNA presence when helpful.

Design
------
* **Colour** encodes PTM state (4 levels max).
* **Line-style / marker** encodes DNA presence (2 levels).
* Where multiple environments exist, individual seed-level points are
  shown with transparency; the main curves use all data.
* Ranking Score uses box+strip instead of ECDF because it is
  discrete / concentrated.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .config import (
    DPI, DOUBLE_COL_WIDTH, FONT_SIZES, MAX_FIGURE_HEIGHT_INCHES,
    PTM_COLOURS, PTM_ORDER, DNA_ORDER, DNA_STYLE, SEED_STYLE,
    apply_v2_style, get_ptm_color, get_dna_style, make_ptm_palette,
)
from .factors import add_factor_columns, get_unique_ptm_states, get_unique_environments
from .labels import (
    FIGURE7_METRICS, metric_label, panel_letter, figure_title,
)


def generate_figure7(
    seed_aggregated: pd.DataFrame,
    save_path: Path,
    *,
    design: Any = None,
    environment_filter: Optional[str] = None,
) -> Dict[str, Any]:
    """Generate the complete Figure 7 distributional panel.

    Parameters
    ----------
    seed_aggregated : DataFrame
        Must contain ``condition_name``, ``seed``, and metric columns.
    save_path : Path
        Directory to write the output PNG.
    design : ExperimentDesign, optional
        Experiment metadata for labels.
    environment_filter : str, optional
        If given, restrict to this environment value (e.g. ``"baseline"``).
        ``None`` means plot all environments pooled.

    Returns
    -------
    dict with ``status``, ``output_path``, ``metrics``, ``n_observations``,
    ``warnings``.
    """
    apply_v2_style()

    df = seed_aggregated.copy()
    if "condition_name" not in df.columns:
        return {"status": "skip", "reason": "No condition_name column"}

    # Augment with factor columns
    df = add_factor_columns(df, design=design)

    # Optional environment filter
    if environment_filter is not None:
        df = df[df["environment"] == environment_filter]
        if df.empty:
            return {"status": "skip", "reason": f"No data for environment={environment_filter}"}

    # Select available metrics
    available_metrics = [m for m in FIGURE7_METRICS if m in df.columns]
    if not available_metrics:
        return {"status": "skip", "reason": "No metrics available"}

    n_metrics = len(available_metrics)

    # Determine whether to facet by DNA or show combined
    # For the main figure, pool across DNA and encode via line style
    # Figure dimensions: thesis double-column width
    fig_height = 4.0 * n_metrics
    fig_height = min(fig_height, MAX_FIGURE_HEIGHT_INCHES)
    fig, axes = plt.subplots(n_metrics, 1, figsize=(DOUBLE_COL_WIDTH, fig_height))
    if n_metrics == 1:
        axes = [axes]

    palette = make_ptm_palette()
    ptm_states = [p for p in PTM_ORDER if p in df["ptm_state"].unique()]

    warnings_list: List[str] = []
    total_obs = 0

    for idx, (ax, metric) in enumerate(zip(axes, available_metrics)):
        if metric not in df.columns:
            ax.set_visible(False)
            continue

        letter = panel_letter(metric)
        label = metric_label(metric)

        # Check if ranking score (use box+strip instead of ECDF)
        if metric == "ranking_score":
            _plot_ranking_box(ax, df, metric, palette, ptm_states, label)
        else:
            _plot_ecdf_dual(ax, df, metric, palette, ptm_states, label)

        ax.set_title(f"  {letter}. {label}", loc="left",
                      fontsize=FONT_SIZES["panel_title"], fontweight="bold")

        total_obs += df[metric].notna().sum()

    # --- Legend ---
    # Build a combined legend: PTM colours + DNA line styles
    _add_combined_legend(fig, ptm_states, list(DNA_STYLE.keys()))

    fig.suptitle(
        figure_title("fig7"),
        fontsize=FONT_SIZES["figure_title"],
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout(rect=[0, 0, 0.88, 0.98])

    out = save_path / "figure_7_distribution_v2.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    # --- Post-generation validation ---
    _validate_figure7(df, available_metrics, ptm_states, warnings_list)

    return {
        "status": "pass",
        "output_path": str(out),
        "metrics": available_metrics,
        "n_observations": total_obs,
        "ptm_states": ptm_states,
        "warnings": warnings_list,
    }


# -----------------------------------------------------------------------
# Internal plot helpers
# -----------------------------------------------------------------------

def _plot_ecdf_dual(ax, df, metric, palette, ptm_states, label):
    """ECDF with colour=PTM, line-style=DNA."""
    for ptm in ptm_states:
        for has_dna in DNA_ORDER:
            subset = df[(df["ptm_state"] == ptm) & (df["has_dna"] == has_dna)]
            if subset.empty:
                continue
            style = get_dna_style(has_dna)
            sns.ecdfplot(
                data=subset, x=metric, ax=ax,
                color=palette.get(ptm, "grey"),
                linestyle=style["linestyle"],
                linewidth=1.8,
                label=f"{ptm} ({style['label']})",
            )
    ax.set_xlabel(label)
    ax.set_ylabel("Cumulative proportion")
    ax.set_xlim(left=None, right=None)


def _plot_ranking_box(ax, df, metric, palette, ptm_states, label):
    """Box + jittered strip for the discrete ranking score."""
    # Build a combined x label for the box plot
    df_plot = df.copy()
    df_plot["ptm_dna"] = df_plot.apply(
        lambda r: r["ptm_state"] + (" + DNA" if r["has_dna"] else ""), axis=1
    )
    order = []
    for ptm in ptm_states:
        for has_dna in DNA_ORDER:
            lbl = ptm + (" + DNA" if has_dna else "")
            if (df_plot["ptm_dna"] == lbl).any():
                order.append(lbl)

    box_palette = {}
    for ptm in ptm_states:
        for has_dna in DNA_ORDER:
            lbl = ptm + (" + DNA" if has_dna else "")
            base = palette.get(ptm, "grey")
            # Slightly darken for DNA
            box_palette[lbl] = base

    sns.boxplot(
        data=df_plot, x="ptm_dna", y=metric, order=order,
        hue="ptm_dna", hue_order=order, palette=box_palette,
        ax=ax, width=0.6, showfliers=False, linewidth=1.0,
        legend=False,
    )
    sns.stripplot(
        data=df_plot, x="ptm_dna", y=metric, order=order,
        color="black", alpha=SEED_STYLE["alpha"],
        size=SEED_STYLE["s"] * 0.4, jitter=0.2, ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel(label)
    ax.tick_params(axis="x", rotation=20, labelsize=FONT_SIZES["tick_label"])


def _add_combined_legend(fig, ptm_states, dna_values):
    """Create a two-part legend: PTM colours + DNA line styles."""
    from matplotlib.lines import Line2D
    handles = []
    labels = []

    # PTM section
    handles.append(Line2D([0], [0], color="white", linewidth=0))
    labels.append("PTM State")
    for ptm in ptm_states:
        c = PTM_COLOURS.get(ptm, "grey")
        handles.append(Line2D([0], [0], color=c, linewidth=2.5))
        labels.append(ptm)

    # DNA section
    handles.append(Line2D([0], [0], color="white", linewidth=0))
    labels.append("DNA")
    for has_dna in dna_values:
        style = DNA_STYLE[has_dna]
        handles.append(
            Line2D([0], [0], color="grey",
                   linestyle=style["linestyle"],
                   marker=style["marker"], markersize=6, linewidth=1.5)
        )
        labels.append(style["label"])

    fig.legend(
        handles, labels,
        loc="center right",
        bbox_to_anchor=(0.99, 0.5),
        frameon=True, framealpha=0.9,
        fontsize=FONT_SIZES["legend_text"],
        title_fontsize=FONT_SIZES["legend_title"],
    )


# -----------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------

def _validate_figure7(df, metrics, ptm_states, warnings_list):
    """Post-generation validation checks."""
    for metric in metrics:
        n_groups = df.groupby(["ptm_state", "has_dna"])[metric].apply(
            lambda s: s.notna().sum()
        )
        for (ptm, dna), count in n_groups.items():
            if count == 0:
                warnings_list.append(
                    f"WARNING: {metric} — {ptm} (DNA={dna}) has zero observations"
                )
    if len(ptm_states) < 2:
        warnings_list.append(
            f"WARNING: Only {len(ptm_states)} PTM state(s) present — "
            f"factor comparison is limited"
        )
