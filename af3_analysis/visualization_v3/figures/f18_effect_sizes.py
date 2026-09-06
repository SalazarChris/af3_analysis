"""
V3 Figure F18 — Structural Effect Sizes

Shows structural effect sizes for condition/reference comparisons:
- RMSD effect
- local displacement effect
- contact-change effect
- interface-change effect
- domain-motion effect

Includes effect estimate, uncertainty, valid seeds, direction consistency.
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
    DOUBLE_COL_WIDTH,
    apply_v3_style,
)


def generate_f18_effect_sizes(
    effect_size_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate effect sizes figure.

    Parameters
    ----------
    effect_size_data : list of effect size dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not effect_size_data:
        return {
            "status": "skip",
            "reason": "No effect size data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No effect size data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(effect_size_data)

    required_cols = ["condition_id", "estimate", "metric_id"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Effect size data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Effect size data incomplete"],
        }

    df = df[df["estimate"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid effect size values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All effect size values are NaN"],
        }

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        df["condition_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["condition_label"] = df["condition_id"]

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition_id"].unique()]
        order += [c for c in df["condition_id"].unique() if c not in order]
    else:
        order = sorted(df["condition_id"].unique())

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}

    n_obs = len(df)

    # --- Plot ---
    # Determine figure size based on number of effects
    n_metrics = df["metric_id"].nunique()
    fig_height = min(8.0, max(4.0, n_metrics * 0.5))
    fig, ax = plt.subplots(1, 1, figsize=(DOUBLE_COL_WIDTH, fig_height))

    palette = sns.color_palette("Set2", n_colors=max(len(order), 2))

    # Create forest plot
    df_plot = df.sort_values("estimate")

    y_positions = range(len(df_plot))

    # Plot points
    colors = [palette[order.index(c) % len(palette)] if c in order else "#BBBBBB"
              for c in df_plot["condition_id"]]

    ax.scatter(
        df_plot["estimate"],
        y_positions,
        c=colors,
        s=50,
        alpha=0.7,
        edgecolors="white",
        linewidth=0.3,
        zorder=3,
    )

    # Plot error bars
    if "ci_lower" in df_plot.columns and "ci_upper" in df_plot.columns:
        for i, row in df_plot.iterrows():
            if pd.notna(row["ci_lower"]) and pd.notna(row["ci_upper"]):
                ax.plot(
                    [row["ci_lower"], row["ci_upper"]],
                    [i, i],
                    color=colors[i],
                    linewidth=1.0,
                    alpha=0.7,
                )

    # Zero line
    ax.axvline(0, color="grey", linestyle="-", linewidth=0.8, alpha=0.5)

    # Labels
    labels = [f"{row['metric_id']}\n({row['condition_label']})"
              for i, row in df_plot.iterrows()]

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Effect estimate")
    ax.set_title(
        f"  {title or 'Structural Effect Sizes (condition vs reference)'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.tick_params(axis="x", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.grid(True, alpha=0.3, axis="x")

    # Main title
    fig.suptitle("Structural Effect Sizes", fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f18_effect_sizes.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero effect size observations")
    if n_metrics > 10:
        warnings.append(f"Large number of metrics ({n_metrics}), figure may be dense")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
