"""
V3 Figure F10 — Local/PTM-Site Geometry

Shows local RMSD, local displacement, and local confidence for configured sites/regions.
Works for arbitrary proteins and residues (configuration-driven, not hard-coded).
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


def generate_f10_local_geometry(
    local_geometry_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate local geometry figure.

    Parameters
    ----------
    local_geometry_data : list of local geometry result dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not local_geometry_data:
        return {
            "status": "skip",
            "reason": "No local geometry data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No local geometry data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(local_geometry_data)

    required_cols = ["region_label", "local_rmsd"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Local geometry data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Local geometry data incomplete"],
        }

    df = df[df["local_rmsd"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid local RMSD values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All local RMSD values are NaN"],
        }

    # Build labels for regions
    region_labels = sorted(df["region_label"].unique())

    # Condition labels and ordering (design-driven where available; no
    # hardcoded condition identities).
    if design and hasattr(design, "get_condition_label"):
        df["condition_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["condition_label"] = df["condition_id"]

    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition_id"].unique()]
        order += [c for c in df["condition_id"].unique() if c not in order]
    else:
        order = sorted(df["condition_id"].unique())

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}
    categories = [label_map.get(c, c) for c in order]
    df["plot_label"] = pd.Categorical(df["condition_label"], categories=categories, ordered=True)

    palette = sns.color_palette("Set2", n_colors=max(len(order), 3))
    color_map = {label_map.get(c, c): palette[i] for i, c in enumerate(order)}

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    # Panel A: Local RMSD by region, grouped per condition so the figure
    # shows which condition drives each region's local difference.
    ax = axes[0]

    region_order = sorted(df["region_label"].unique())
    region_data = (
        df.groupby(["region_label", "plot_label"], observed=True)["local_rmsd"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )

    x_pos = np.arange(len(region_order))
    n_cond = max(len(order), 1)
    bar_width = 0.8 / n_cond
    labeled_conditions = set()
    for j, cond in enumerate(order):
        cond_label = label_map.get(cond, cond)
        sub = (
            region_data[region_data["plot_label"] == cond_label]
            .set_index("region_label")
        )
        offset = (j - (n_cond - 1) / 2) * bar_width
        for k, region in enumerate(region_order):
            if region not in sub.index:
                continue
            row = sub.loc[region]
            mean_val = row["mean"]
            if pd.isna(mean_val):
                continue
            yerr = row["std"] if pd.notna(row["std"]) else None
            ax.bar(
                k + offset, mean_val, bar_width * 0.9,
                yerr=yerr,
                color=color_map[cond_label],
                edgecolor="white",
                linewidth=0.8,
                capsize=2,
                label=cond_label if cond_label not in labeled_conditions else None,
            )
            labeled_conditions.add(cond_label)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(region_order, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Local RMSD (Å)")
    ax.set_title("  A. Local RMSD by Region", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=9)

    if len(labeled_conditions) > 1:
        ax.legend(fontsize=7, loc="best")

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Local confidence vs RMSD (scatter), colored by condition so
    # the continuous confidence-geometry relationship stays visible.
    ax = axes[1]

    if "local_plddt_mean" in df.columns:
        df_plot = df[df["local_plddt_mean"].notna()].copy()
        if not df_plot.empty:
            for cond in order:
                cond_label = label_map.get(cond, cond)
                sub = df_plot[df_plot["condition_id"] == cond]
                if sub.empty:
                    continue
                ax.scatter(
                    sub["local_plddt_mean"],
                    sub["local_rmsd"],
                    color=color_map[cond_label],
                    s=30,
                    alpha=0.6,
                    edgecolors="white",
                    linewidth=0.3,
                    label=cond_label,
                )
            ax.set_xlabel("Local mean pLDDT")
            ax.set_ylabel("Local RMSD (Å)")
            handles, labels = ax.get_legend_handles_labels()
            if len(handles) > 12:
                ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7)
            elif len(handles) > 1:
                ax.legend(fontsize=7, loc="best")
        else:
            ax.text(0.5, 0.5, "No confidence data", ha="center", va="center")
            ax.set_xlabel("Local mean pLDDT")
            ax.set_ylabel("Local RMSD (Å)")
    else:
        ax.text(0.5, 0.5, "No confidence data", ha="center", va="center")
        ax.set_xlabel("Local mean pLDDT")
        ax.set_ylabel("Local RMSD (Å)")

    ax.set_title("  B. Local Confidence vs RMSD", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.grid(True, alpha=0.3)

    # Main title. When the runner passes no title, make the provenance of
    # data-detected sites explicit here; runner-provided titles already
    # carry this annotation.
    detected_in_data = (
        "site_detected" in df.columns
        and bool(df["site_detected"].fillna(False).any())
    )
    if title is not None:
        main_title = title
    elif detected_in_data:
        main_title = (
            "Local/PTM-Site Geometry — predicted structural sites "
            "(detected from data; not functional annotations)"
        )
    else:
        main_title = "Local/PTM-Site Geometry"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f10_local_geometry.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero local geometry observations")
    if len(region_order) > 10:
        warnings.append(f"Large number of regions ({len(region_order)}), figure may be dense")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
