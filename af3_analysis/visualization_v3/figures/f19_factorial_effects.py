"""
V3 Figure F19 — Factorial Structural Effects

If condition registry defines a factorial experiment, use those metadata fields.
Does NOT parse condition names.

Estimates structural effects corresponding to:
- main effects
- two-way interactions
- higher-order interactions where supported

Uses metadata-driven contrasts, not name parsing.
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


def generate_f19_factorial_effects(
    factorial_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate factorial effects figure.

    Parameters
    ----------
    factorial_data : list of factorial contrast dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not factorial_data:
        return {
            "status": "skip",
            "reason": "No factorial data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No factorial structural effect data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(factorial_data)

    required_cols = ["contrast_id", "estimate"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Factorial data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Factorial data incomplete"],
        }

    df = df[df["estimate"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid factorial effect values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All factorial effect values are NaN"],
        }

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    # Panel A: Effect estimates (forest plot style)
    ax = axes[0]

    df_plot = df.sort_values("estimate")

    y_positions = range(len(df_plot))

    colors = sns.color_palette("Set2", n_colors=len(df_plot))

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

    # Error bars if available
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

    ax.axvline(0, color="grey", linestyle="-", linewidth=0.8, alpha=0.5)

    labels = [row["contrast_id"] for i, row in df_plot.iterrows()]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Effect estimate")
    ax.set_title("  A. Factorial Contrast Estimates", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.grid(True, alpha=0.3, axis="x")

    # Panel B: Contrast types (bar)
    ax = axes[1]

    if "contrast_type" in df.columns:
        type_counts = df["contrast_type"].value_counts()

        bars = ax.bar(
            type_counts.index,
            type_counts.values,
            0.6,
            color=sns.color_palette("Set2", len(type_counts)),
            edgecolor="white",
            linewidth=0.8,
        )

        ax.set_xlabel("Contrast type")
        ax.set_ylabel("Number of contrasts")
        ax.set_title("  B. Contrast Types", loc="left", fontsize=11, fontweight="bold")
        ax.tick_params(axis="x", rotation=30, labelsize=9)
        ax.tick_params(axis="y", labelsize=9)

        # Add count labels
        for bar, count in zip(bars, type_counts.values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                    str(count), ha="center", va="bottom", fontsize=9)

        sns.despine(ax=ax, left=True)
        ax.yaxis.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, "No contrast types", ha="center", va="center")
        ax.set_title("  B. Contrast Types", loc="left", fontsize=11, fontweight="bold")

    # Main title
    main_title = title or "Factorial Structural Effects"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f19_factorial_effects.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero factorial effect observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
