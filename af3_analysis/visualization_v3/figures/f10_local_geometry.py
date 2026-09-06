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

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    # Panel A: Local RMSD by region
    ax = axes[0]

    region_data = df.groupby("region_label")["local_rmsd"].agg(["mean", "std", "count"]).reset_index()
    region_data = region_data.sort_values("region_label")

    x_pos = np.arange(len(region_data))
    ax.bar(
        x_pos,
        region_data["mean"],
        0.6,
        yerr=region_data["std"],
        color="#2C7BB6",
        edgecolor="white",
        linewidth=0.8,
        capsize=3,
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(region_data["region_label"], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Local RMSD (Å)")
    ax.set_title("  A. Local RMSD by Region", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Local confidence vs RMSD (scatter)
    ax = axes[1]

    if "local_plddt_mean" in df.columns:
        df_plot = df[df["local_plddt_mean"].notna()].copy()
        if not df_plot.empty:
            scatter = ax.scatter(
                df_plot["local_plddt_mean"],
                df_plot["local_rmsd"],
                c=df_plot["local_rmsd"],
                cmap="RdYlBu_r",
                s=30,
                alpha=0.6,
                edgecolors="white",
                linewidth=0.3,
            )
            plt.colorbar(scatter, ax=ax, label="Local RMSD (Å)")
            ax.set_xlabel("Local mean pLDDT")
            ax.set_ylabel("Local RMSD (Å)")
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

    # Main title
    main_title = title or "Local/PTM-Site Geometry"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f10_local_geometry.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero local geometry observations")
    if len(region_data) > 10:
        warnings.append(f"Large number of regions ({len(region_data)}), figure may be dense")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
