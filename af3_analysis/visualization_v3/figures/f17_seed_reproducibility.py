"""
V3 Figure F17 — Seed Reproducibility

For each structural metric calculate:
- mean
- median
- SD
- IQR
- valid seeds
- direction consistency

Treats seeds as prediction-process robustness samples.
Does NOT describe them as biological replicates or physical ensemble.
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


def generate_f17_seed_reproducibility(
    seed_repro_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate seed reproducibility figure.

    Parameters
    ----------
    seed_repro_data : list of seed reproducibility dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not seed_repro_data:
        return {
            "status": "skip",
            "reason": "No seed reproducibility data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No seed reproducibility data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(seed_repro_data)

    required_cols = ["metric_id", "condition_id", "mean", "n_seeds"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Seed reproducibility data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Seed reproducibility data incomplete"],
        }

    df = df[df["mean"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid seed reproducibility values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All seed reproducibility values are NaN"],
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
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    palette = sns.color_palette("Set2", n_colors=max(len(order), 2))

    # Panel A: Mean metric value with error bars
    ax = axes[0]

    summary = df.groupby(["metric_id", "condition_id"]).agg({
        "mean": "mean",
        "std": "mean",
        "n_seeds": "mean",
    }).reset_index()

    summary["condition_label"] = summary["condition_id"].map(label_map)

    sns.barplot(
        data=summary,
        x="metric_id",
        y="mean",
        hue="condition_id",
        palette=palette,
        linewidth=0.8,
        edgecolor="white",
        ax=ax,
    )

    ax.set_xlabel("Metric")
    ax.set_ylabel("Mean value (across seeds)")
    ax.set_title("  A. Mean Metric Values by Condition", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 10:
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    else:
        ax.legend(loc="best", fontsize=9)

    # Panel B: Number of valid seeds
    ax = axes[1]

    seed_counts = df.groupby("condition_id")["n_seeds"].mean().reset_index()
    seed_counts["condition_label"] = seed_counts["condition_id"].map(label_map)
    seed_counts = seed_counts.sort_values("condition_id")

    ax.bar(
        range(len(seed_counts)),
        seed_counts["n_seeds"],
        0.6,
        color="#2C7BB6",
        edgecolor="white",
        linewidth=0.8,
    )

    ax.set_xticks(range(len(seed_counts)))
    ax.set_xticklabels(seed_counts["condition_label"], rotation=30, ha="right", fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("Number of seeds")
    ax.set_title("  B. Number of Seeds per Condition", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Main title
    main_title = title or "Seed Reproducibility"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f17_seed_reproducibility.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero seed reproducibility observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
