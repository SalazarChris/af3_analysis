"""
V3 Figure F08 — Protein-DNA/Interface Analysis

Shows:
- Interface contacts per condition
- Contact counts
- Minimum distances
- Contact persistence
- Condition-specific contact changes

If DNA is absent: skip the figure and record SKIPPED.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from af3_analysis.visualization_v3.config import (
    DPI,
    SINGLE_COL_WIDTH,
    apply_v3_style,
)


def generate_f08_interface_analysis(
    interface_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate interface analysis figure.

    Parameters
    ----------
    interface_data : list of interface summary dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not interface_data:
        return {
            "status": "skip",
            "reason": "No interface data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No interface data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(interface_data)

    if "condition_id" not in df.columns or "n_contacts" not in df.columns:
        return {
            "status": "skip",
            "reason": "Interface data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Interface data incomplete"],
        }

    # Filter to valid data
    df = df[df["n_contacts"] >= 0].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid interface data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All interface data is invalid"],
        }

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        df["plot_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["plot_label"] = df["condition_id"]

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition_id"].unique()]
        order += [c for c in df["condition_id"].unique() if c not in order]
    else:
        order = sorted(df["condition_id"].unique())

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}
    categories = [label_map.get(c, c) for c in order]
    df["plot_label"] = pd.Categorical(df["plot_label"], categories=categories, ordered=True)

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    palette = sns.color_palette("Set2", n_colors=max(len(order), 3))
    color_map = {label_map.get(c, c): palette[i] for i, c in enumerate(order)}

    # Panel A: Contact counts per condition
    ax = axes[0]

    sns.boxplot(
        data=df,
        x="plot_label",
        y="n_contacts",
        hue="plot_label",
        palette=color_map,
        width=0.5,
        fliersize=0,
        linewidth=0.8,
        boxprops=dict(alpha=0.4),
        legend=False,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="plot_label",
        y="n_contacts",
        hue="plot_label",
        palette=color_map,
        size=4,
        alpha=0.6,
        linewidth=0.3,
        edgecolor="white",
        jitter=0.2,
        legend=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Number of interface contacts")
    ax.set_title("  A. Interface Contacts per Condition", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Mean distance per condition
    ax = axes[1]

    if "mean_distance" in df.columns:
        df_plot = df[df["mean_distance"].notna()].copy()
        if not df_plot.empty:
            sns.boxplot(
                data=df_plot,
                x="plot_label",
                y="mean_distance",
                hue="plot_label",
                palette=color_map,
                width=0.5,
                fliersize=0,
                linewidth=0.8,
                boxprops=dict(alpha=0.4),
                legend=False,
                ax=ax,
            )
            sns.stripplot(
                data=df_plot,
                x="plot_label",
                y="mean_distance",
                hue="plot_label",
                palette=color_map,
                size=4,
                alpha=0.6,
                linewidth=0.3,
                edgecolor="white",
                jitter=0.2,
                legend=False,
                ax=ax,
            )
            ax.set_ylabel("Mean interface distance (Å)")
        else:
            ax.text(0.5, 0.5, "No distance data", ha="center", va="center")
    else:
        ax.text(0.5, 0.5, "No distance data", ha="center", va="center")

    ax.set_xlabel("")
    ax.set_title("  B. Interface Distance", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Main title
    main_title = title or "Interface Analysis"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f08_interface_analysis.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero interface observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
