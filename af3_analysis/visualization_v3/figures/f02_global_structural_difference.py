"""
V3 Figure F02 — Global Structural Difference

Shows RMSD distribution across conditions vs reference.
Box + strip plot of aligned Cα RMSD.
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


def generate_f02_global_structural_difference(
    rmsd_results: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate global structural difference figure.

    Parameters
    ----------
    rmsd_results : list of dicts with RMSD values per comparison
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not rmsd_results:
        return {
            "status": "skip",
            "reason": "No RMSD results",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No structural RMSD data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(rmsd_results)

    if "condition_a" not in df.columns or "rmsd" not in df.columns:
        return {
            "status": "skip",
            "reason": "RMSD results missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["RMSD results incomplete"],
        }

    # Filter to valid RMSD values
    df = df[df["rmsd"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid RMSD values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All RMSD values are NaN"],
        }

    # Determine condition label column
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

    # Build labels
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

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}
    categories = [label_map.get(c, c) for c in order]
    df["plot_label"] = pd.Categorical(df["plot_label"], categories=categories, ordered=True)

    n_obs = len(df)

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(SINGLE_COL_WIDTH, 5.0))

    palette = sns.color_palette("Set2", n_colors=max(len(order), 3))
    color_map = {label_map.get(c, c): palette[i] for i, c in enumerate(order)}

    sns.boxplot(
        data=df,
        x="plot_label",
        y="rmsd",
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
        y="rmsd",
        hue="plot_label",
        palette=color_map,
        size=3,
        alpha=0.4,
        linewidth=0.3,
        edgecolor="white",
        jitter=0.15,
        legend=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Aligned Cα RMSD (Å)")
    ax.set_title(
        f"  {title or 'Global Structural Difference'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    out_path = save_path / "fig_v3_f02_global_structural_difference.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero RMSD observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
