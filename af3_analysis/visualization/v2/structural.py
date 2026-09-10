"""
V2 Structural Figures — Structural Geometry Visualization
==========================================================

Two publication-quality figures for structural analysis results:

* **Figure S1 (RMSD Distribution):** Box + strip plot of aligned Cα RMSD
  across conditions, coloured by condition.

* **Figure S2 (Structural Variability):** Bar chart of within-condition
  structural variability (mean ± SD RMSD) across conditions.

Design
------
* Generic — works for any AF3 project with structural_metrics.csv.
* Uses the same V2 style as confidence-metric figures.
* Condition labels come from metadata when available, else raw IDs.
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

from .config import (
    DPI,
    DOUBLE_COL_WIDTH,
    FONT_SIZES,
    SINGLE_COL_WIDTH,
    apply_v2_style,
)


# ---------------------------------------------------------------------------
# Figure S1 — RMSD Distribution
# ---------------------------------------------------------------------------

def generate_structural_rmsd_distribution(
    structural_comparisons: pd.DataFrame,
    save_path: Path,
    *,
    design: Any = None,
    reference_condition: Optional[str] = None,
    condition_order: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Generate a box + strip plot of comparison RMSD by target condition.

    Uses ``structural_comparisons.csv`` which contains meaningful
    between-prediction RMSD values (not self-RMSD).

    Parameters
    ----------
    structural_comparisons : DataFrame
        Must contain ``condition_a``, ``condition_b``, ``metric_id``,
        ``value``, ``alignment_status``.
    save_path : Path
        Directory to write the output PNG.
    design : ExperimentDesign, optional
        For human-readable condition labels.
    reference_condition : str, optional
        The reference condition to filter comparisons to.
    condition_order : list of str, optional
        Explicit condition ordering on the x-axis.

    Returns
    -------
    dict with ``status``, ``output_path``, ``n_observations``, ``warnings``.
    """
    apply_v2_style()
    warnings: List[str] = []

    if structural_comparisons is None or structural_comparisons.empty:
        return {"status": "skip", "reason": "No structural comparisons data"}

    # Filter to global Cα RMSD with valid values and comparable status
    df = structural_comparisons.copy()
    df = df[
        (df["metric_id"] == "rmsd_global_ca")
        & (df["alignment_status"] == "comparable")
        & (df["value"].notna())
    ].copy()

    if df.empty:
        return {"status": "skip", "reason": "No valid rmsd_global_ca comparison values"}

    # If a reference condition is specified, filter to comparisons against it
    if reference_condition and reference_condition in df["condition_b"].unique():
        df = df[df["condition_b"] == reference_condition].copy()
        target_col = "condition_a"
        title_suffix = f" vs {reference_condition}"
    elif reference_condition and reference_condition in df["condition_a"].unique():
        df = df[df["condition_a"] == reference_condition].copy()
        target_col = "condition_b"
        title_suffix = f" (ref: {reference_condition})"
    else:
        # Use condition_a as the target (condition_b is reference)
        target_col = "condition_a"
        ref_vals = df["condition_b"].unique()
        title_suffix = f" vs {ref_vals[0]}" if len(ref_vals) == 1 else ""

    if df.empty:
        return {"status": "skip", "reason": "No comparisons match the reference filter"}

    # Build condition labels
    all_conditions = df[target_col].unique()
    label_map = _build_condition_label_map_from_ids(all_conditions, design)
    df["condition_label"] = df[target_col].map(label_map)

    # Sort conditions
    if condition_order:
        ordered = [c for c in condition_order if c in df[target_col].unique()]
    else:
        ordered = sorted(df[target_col].unique())
    categories = [label_map.get(c, c) for c in ordered]
    df["condition_label"] = pd.Categorical(df["condition_label"], categories=categories, ordered=True)

    # Colour palette
    n_conds = len(ordered)
    palette = sns.color_palette("Set2", n_colors=max(n_conds, 3))
    colour_map = {label_map.get(c, c): palette[i] for i, c in enumerate(ordered)}

    n_obs = len(df)

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(SINGLE_COL_WIDTH, 5.0))

    sns.boxplot(
        data=df,
        x="condition_label",
        y="value",
        hue="condition_label",
        palette=colour_map,
        width=0.5,
        fliersize=0,
        linewidth=0.8,
        boxprops=dict(alpha=0.4),
        legend=False,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="condition_label",
        y="value",
        hue="condition_label",
        palette=colour_map,
        size=3,
        alpha=0.4,
        linewidth=0.3,
        edgecolor="white",
        jitter=0.15,
        legend=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Aligned C\u03b1 RMSD (\u00c5)")
    ax.set_title(
        f"Structural RMSD{title_suffix}",
        fontsize=FONT_SIZES["panel_title"],
        fontweight="bold",
    )
    ax.tick_params(axis="x", rotation=30, labelsize=FONT_SIZES["tick_label"])
    ax.tick_params(axis="y", labelsize=FONT_SIZES["tick_label"])

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    out_path = save_path / "fig_structural_rmsd_distribution.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Figure S2 — Structural Variability by Condition
# ---------------------------------------------------------------------------

def generate_structural_variability(
    structural_comparisons: pd.DataFrame,
    save_path: Path,
    *,
    design: Any = None,
    reference_condition: Optional[str] = None,
    condition_order: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Generate a bar chart of between-condition structural variability.

    For each condition (vs reference), shows:
    - mean RMSD across matched-seed comparisons (bar height)
    - SD of RMSD (error bar)

    Uses ``structural_comparisons.csv`` for meaningful comparison values.

    Parameters
    ----------
    structural_comparisons : DataFrame
        Must contain ``condition_a``, ``condition_b``, ``metric_id``,
        ``value``, ``alignment_status``.
    save_path : Path
        Directory to write the output PNG.
    design : ExperimentDesign, optional
        For human-readable condition labels.
    reference_condition : str, optional
        The reference condition to filter comparisons to.
    condition_order : list of str, optional
        Explicit condition ordering on the x-axis.

    Returns
    -------
    dict with ``status``, ``output_path``, ``n_observations``, ``warnings``.
    """
    apply_v2_style()
    warnings: List[str] = []

    if structural_comparisons is None or structural_comparisons.empty:
        return {"status": "skip", "reason": "No structural comparisons data"}

    df = structural_comparisons.copy()
    df = df[
        (df["metric_id"] == "rmsd_global_ca")
        & (df["alignment_status"] == "comparable")
        & (df["value"].notna())
    ].copy()

    if df.empty:
        return {"status": "skip", "reason": "No valid rmsd_global_ca comparison values"}

    # Filter to reference condition comparisons
    if reference_condition and reference_condition in df["condition_b"].unique():
        df = df[df["condition_b"] == reference_condition].copy()
        target_col = "condition_a"
    elif reference_condition and reference_condition in df["condition_a"].unique():
        df = df[df["condition_a"] == reference_condition].copy()
        target_col = "condition_b"
    else:
        target_col = "condition_a"

    if df.empty:
        return {"status": "skip", "reason": "No comparisons match the reference filter"}

    # Aggregate per condition
    stats = (
        df.groupby(target_col)["value"]
        .agg(["mean", "std", "count", "median"])
        .reset_index()
        .rename(columns={"mean": "mean_rmsd", "std": "sd_rmsd", "count": "n", "median": "median_rmsd", target_col: "condition_id"})
    )
    stats["sd_rmsd"] = stats["sd_rmsd"].fillna(0)

    # Build labels
    label_map = _build_condition_label_map_from_ids(stats["condition_id"].values, design)
    stats["condition_label"] = stats["condition_id"].map(label_map)

    # Sort
    if condition_order:
        ordered = [c for c in condition_order if c in stats["condition_id"].values]
    else:
        ordered = sorted(stats["condition_id"].values)
    stats = stats.set_index("condition_id").loc[ordered].reset_index()
    stats["condition_label"] = stats["condition_id"].map(label_map)

    n_obs = int(stats["n"].sum())

    # Colour palette
    n_conds = len(ordered)
    palette = sns.color_palette("Set2", n_colors=max(n_conds, 3))
    colours = [palette[i] for i in range(n_conds)]

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(SINGLE_COL_WIDTH, 5.0))

    ax.bar(
        range(n_conds),
        stats["mean_rmsd"],
        yerr=stats["sd_rmsd"],
        color=colours,
        edgecolor="white",
        linewidth=0.8,
        capsize=4,
        error_kw={"linewidth": 1.2, "capthick": 1.2},
        width=0.6,
    )

    # Annotate n on top of each bar
    for i, (_, row) in enumerate(stats.iterrows()):
        ax.text(
            i,
            row["mean_rmsd"] + row["sd_rmsd"] + 0.15,
            f"n={int(row['n'])}",
            ha="center",
            va="bottom",
            fontsize=FONT_SIZES["annotation"],
            color="0.4",
        )

    ax.set_xticks(range(n_conds))
    ax.set_xticklabels(
        stats["condition_label"],
        rotation=30,
        ha="right",
        fontsize=FONT_SIZES["tick_label"],
    )
    ax.set_ylabel("Mean RMSD (\u00c5) \u00b1 SD")
    ax.set_title(
        "Between-Condition Structural Variability",
        fontsize=FONT_SIZES["panel_title"],
        fontweight="bold",
    )
    ax.tick_params(axis="y", labelsize=FONT_SIZES["tick_label"])

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    out_path = save_path / "fig_structural_variability.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_condition_label_map_from_ids(
    condition_ids,
    design: Any = None,
) -> Dict[str, str]:
    """Map condition IDs to human-readable labels when metadata is available."""
    label_map: Dict[str, str] = {}
    if hasattr(condition_ids, "tolist"):
        ids = condition_ids.tolist()
    else:
        ids = list(condition_ids)

    if design is not None and hasattr(design, "conditions"):
        for cid in ids:
            cond = design.conditions.get(cid)
            if cond is not None and hasattr(cond, "name") and cond.name:
                label_map[cid] = cond.name
            else:
                label_map[cid] = cid
    else:
        for cid in ids:
            label_map[cid] = cid

    return label_map


def _build_condition_label_map(
    df: pd.DataFrame,
    design: Any = None,
) -> Dict[str, str]:
    """Map condition_id to human-readable labels when metadata is available."""
    label_map = {}
    ids = df["condition_id"].unique()

    if design is not None and hasattr(design, "conditions"):
        for cid in ids:
            cond = design.conditions.get(cid)
            if cond is not None and hasattr(cond, "name") and cond.name:
                label_map[cid] = cond.name
            else:
                label_map[cid] = cid
    else:
        for cid in ids:
            label_map[cid] = cid

    return label_map
