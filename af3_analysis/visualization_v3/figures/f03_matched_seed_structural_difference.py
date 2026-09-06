"""
V3 Figure F03 — Matched-Seed Structural Difference

Shows:
- RMSD per seed
- Coverage per seed
- Valid seed count
- Mean, median, SD, IQR
- Direction consistency where applicable

Extension of F02: adds seed-level detail and reproducibility.
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


def generate_f03_matched_seed_structural_difference(
    seed_rmsd_results: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate matched-seed structural difference figure.

    Parameters
    ----------
    seed_rmsd_results : list of dicts with seed-level RMSD results
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not seed_rmsd_results:
        return {
            "status": "skip",
            "reason": "No seed-level RMSD results",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No matched-seed RMSD data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(seed_rmsd_results)

    if "seed" not in df.columns or "rmsd_mean" not in df.columns:
        return {
            "status": "skip",
            "reason": "Seed results missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Seed results incomplete"],
        }

    # Filter to valid results
    df = df[df["rmsd_mean"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid seed-level RMSD",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All seed RMSD values are NaN"],
        }

    # Build labels
    if "condition_label" in df.columns:
        label_col = "condition_label"
    elif "condition_a" in df.columns:
        label_col = "condition_a"
    else:
        return {
            "status": "skip",
            "reason": "No condition column",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No condition identifier"],
        }

    if design and hasattr(design, "get_condition_label"):
        condition_map = {
            c: design.get_condition_label(c)
            for c in df[label_col].unique()
        }
        df["plot_label"] = df[label_col].map(condition_map)
    else:
        df["plot_label"] = df[label_col]

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df[label_col].unique()]
        order += [c for c in df[label_col].unique() if c not in order]
    else:
        order = sorted(df[label_col].unique())

    if design and hasattr(design, "get_condition_label"):
        label_map = {c: design.get_condition_label(c) for c in order}
    else:
        label_map = {c: c for c in order}

    categories = [label_map.get(c, c) for c in order]
    df["plot_label"] = pd.Categorical(df["plot_label"], categories=categories, ordered=True)

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0),
                            gridspec_kw={"width_ratios": [2, 1]})

    palette = sns.color_palette("Set2", n_colors=max(len(order), 3))
    color_map = {label_map.get(c, c): palette[i] for i, c in enumerate(order)}

    # Panel A: Seed-level RMSD with individual seeds
    ax = axes[0]

    # Box plot
    sns.boxplot(
        data=df,
        x="plot_label",
        y="rmsd_mean",
        hue="plot_label",
        palette=color_map,
        width=0.5,
        fliersize=0,
        linewidth=0.8,
        boxprops=dict(alpha=0.4),
        legend=False,
        ax=ax,
    )

    # Show individual seeds
    sns.stripplot(
        data=df,
        x="plot_label",
        y="rmsd_mean",
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
    ax.set_ylabel("Mean RMSD per seed (Å)")
    ax.set_title("  A. Seed-Level RMSD", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Coverage per seed
    ax = axes[1]

    if "coverage_mean" in df.columns:
        df_plot = df[df["coverage_mean"].notna()].copy()
        sns.boxplot(
            data=df_plot,
            x="plot_label",
            y="coverage_mean",
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
            y="coverage_mean",
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

        ax.set_ylabel("Coverage")
        ax.set_title("  B. Coverage per Seed", loc="left", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1.05)
    else:
        ax.text(0.5, 0.5, "No coverage data", ha="center", va="center")
        ax.set_title("  B. Coverage per Seed", loc="left", fontsize=11, fontweight="bold")

    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Main title
    main_title = title or "Matched-Seed Structural Difference"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f03_matched_seed_structural_difference.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero seed-level observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
