"""
V3 Figure F15 — Confidence × Geometry

Shows categories:
1. low structural difference + high confidence
2. high structural difference + high confidence
3. low structural difference + low confidence
4. high structural difference + low confidence

Purpose: distinguish confident structural differences from uncertain ones.
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


def generate_f15_confidence_geometry(
    confidence_geometry_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate confidence × geometry figure.

    Parameters
    ----------
    confidence_geometry_data : list of confidence-geometry observation dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not confidence_geometry_data:
        return {
            "status": "skip",
            "reason": "No confidence-geometry data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No confidence-geometry data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(confidence_geometry_data)

    required_cols = ["rmsd_to_reference", "plddt_mean"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Confidence-geometry data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Confidence-geometry data incomplete"],
        }

    df = df[df["rmsd_to_reference"].notna() & df["plddt_mean"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid confidence-geometry values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All confidence-geometry values are NaN"],
        }

    # Classify into categories. These cutoffs are pre-existing descriptive
    # display thresholds inherited from the original figure. They are NOT
    # validated classification boundaries; they are documented on the
    # figure itself so the categorization is transparent.
    rmsd_high = 2.0
    rmsd_low = 0.5
    plddt_high = 70.0
    plddt_low = 50.0
    threshold_note = (
        f"Descriptive display cutoffs: RMSD {rmsd_low}/{rmsd_high} Å, "
        f"pLDDT {plddt_low}/{plddt_high} (not validated boundaries)"
    )

    def classify(row):
        rmsd = row["rmsd_to_reference"]
        plddt = row["plddt_mean"]

        if rmsd <= rmsd_low:
            struct_cat = "low_diff"
        elif rmsd >= rmsd_high:
            struct_cat = "high_diff"
        else:
            struct_cat = "moderate_diff"

        if plddt >= plddt_high:
            conf_cat = "high_conf"
        elif plddt <= plddt_low:
            conf_cat = "low_conf"
        else:
            conf_cat = "moderate_conf"

        return f"{struct_cat}_{conf_cat}"

    df["category"] = df.apply(classify, axis=1)

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        df["condition_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["condition_label"] = df["condition_id"]

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    # Palette sized to the categories actually present (avoids the seaborn
    # warning when fewer categories occur).
    palette = sns.color_palette(
        "Set2", n_colors=max(df["category"].nunique(), 3))

    # Panel A: RMSD vs pLDDT scatter with categories
    ax = axes[0]

    sns.scatterplot(
        data=df,
        x="plddt_mean",
        y="rmsd_to_reference",
        hue="category",
        palette=palette,
        s=40,
        alpha=0.6,
        linewidth=0.3,
        edgecolor="white",
        ax=ax,
    )

    # Add threshold lines, labeled so the cutoffs are readable directly
    # from the panel.
    ax.axhline(rmsd_high, color="red", linestyle="--", alpha=0.5, linewidth=1)
    ax.axhline(rmsd_low, color="green", linestyle="--", alpha=0.5, linewidth=1)
    ax.axvline(plddt_high, color="blue", linestyle="--", alpha=0.5, linewidth=1)
    ax.axvline(plddt_low, color="orange", linestyle="--", alpha=0.5, linewidth=1)
    x_min, x_max = ax.get_xlim()
    ax.text(x_max, rmsd_high, f" {rmsd_high} Å", va="bottom", ha="right",
            fontsize=7, color="red", alpha=0.8)
    ax.text(x_max, rmsd_low, f" {rmsd_low} Å", va="bottom", ha="right",
            fontsize=7, color="green", alpha=0.8)
    y_min, y_max = ax.get_ylim()
    ax.text(plddt_high, y_max, f"pLDDT {plddt_high} ", va="top", ha="right",
            fontsize=7, color="blue", alpha=0.8)
    ax.text(plddt_low, y_max, f"pLDDT {plddt_low} ", va="top", ha="right",
            fontsize=7, color="orange", alpha=0.8)

    ax.set_xlabel("Mean pLDDT")
    ax.set_ylabel("RMSD to Reference (Å)")
    ax.set_title("  A. Confidence × Geometry Space", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.grid(True, alpha=0.3)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 8:
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    else:
        ax.legend(loc="best", fontsize=8)

    # Panel B: Category distribution (bar)
    ax = axes[1]

    cat_counts = df["category"].value_counts()
    cat_order = ["low_diff_high_conf", "high_diff_high_conf",
                 "low_diff_low_conf", "high_diff_low_conf",
                 "moderate_diff_high_conf", "moderate_diff_low_conf",
                 "low_diff_moderate_conf", "high_diff_moderate_conf",
                 "moderate_diff_moderate_conf"]

    existing_cats = [c for c in cat_order if c in cat_counts.index]
    values = [cat_counts.get(c, 0) for c in existing_cats]
    colors_cat = sns.color_palette("Set2", len(existing_cats))

    if values:
        ax.bar(existing_cats, values, color=colors_cat, edgecolor="white", linewidth=0.8)
        ax.set_xticks(range(len(existing_cats)))
        ax.set_xticklabels(existing_cats, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Number of predictions")

        # Add count labels
        for i, v in enumerate(values):
            ax.text(i, v + max(values) * 0.02, str(v), ha="center", va="bottom", fontsize=9)
    else:
        ax.text(0.5, 0.5, "No categories", ha="center", va="center")

    ax.set_title("  B. Category Distribution", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Main title: surface the descriptive cutoffs so the categorization is
    # transparent to the reader.
    main_title = title or f"Confidence × Geometry ({threshold_note})"
    fig.suptitle(main_title, fontsize=12, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f15_confidence_geometry.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero confidence-geometry observations")
    warnings.append(
        "Category cutoffs are descriptive display thresholds inherited from "
        "the original figure, not validated classification boundaries"
    )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
