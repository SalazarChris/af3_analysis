"""
V3 Figure F16 — Confidence Change vs Structural Change

For each condition vs reference:
- Delta confidence (ΔpLDDT, ΔPAE)
- Delta structural (RMSD)
- Correlation between them

Purpose: understand if confidence changes accompany structural changes.
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


def generate_f16_confidence_change_vs_structural_change(
    delta_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate confidence change vs structural change figure.

    Parameters
    ----------
    delta_data : list of delta observation dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not delta_data:
        return {
            "status": "skip",
            "reason": "No delta data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No confidence change vs structural change data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(delta_data)

    required_cols = ["condition_id", "delta_rmsd", "delta_plddt"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Delta data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Delta data incomplete"],
        }

    df = df[df["delta_rmsd"].notna() & df["delta_plddt"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid delta values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All delta values are NaN"],
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

    # Panel A: Delta RMSD vs Delta pLDDT scatter
    ax = axes[0]

    sns.scatterplot(
        data=df,
        x="delta_plddt",
        y="delta_rmsd",
        hue="condition_id",
        palette=palette,
        s=50,
        alpha=0.6,
        linewidth=0.3,
        edgecolor="white",
        ax=ax,
    )

    # Add zero lines
    ax.axhline(0, color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
    ax.axvline(0, color="grey", linestyle="-", linewidth=0.5, alpha=0.5)

    ax.set_xlabel("ΔpLDDT (condition - reference)")
    ax.set_ylabel("ΔRMSD (condition - reference)")
    ax.set_title("  A. Confidence Change vs Structural Change", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.grid(True, alpha=0.3)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 10:
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    else:
        ax.legend(loc="best", fontsize=8)

    # Panel B: Summary bar chart
    ax = axes[1]

    summary = df.groupby("condition_id").agg({
        "delta_rmsd": "mean",
        "delta_plddt": "mean",
    }).reset_index()

    summary["condition_label"] = summary["condition_id"].map(label_map)

    x_pos = np.arange(len(summary))
    width = 0.35

    ax.bar(
        x_pos - width / 2,
        summary["delta_rmsd"],
        width,
        label="ΔRMSD",
        color="#2C7BB6",
        edgecolor="white",
    )
    ax.bar(
        x_pos + width / 2,
        summary["delta_plddt"],
        width,
        label="ΔpLDDT",
        color="#D7191C",
        edgecolor="white",
    )

    ax.axhline(0, color="grey", linestyle="-", linewidth=0.5, alpha=0.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(summary["condition_label"], rotation=30, ha="right", fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("Mean Δ value")
    ax.set_title("  B. Mean Changes per Condition", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="best", fontsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Main title
    main_title = title or "Confidence Change vs Structural Change"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f16_confidence_change_vs_structural_change.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero delta observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
