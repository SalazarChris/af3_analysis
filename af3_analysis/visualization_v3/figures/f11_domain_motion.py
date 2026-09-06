"""
V3 Figure F11 — Domain/Region Motion

Shows:
- Region RMSD
- Centroid displacement
- Inter-region distance
- Relative orientation changes

Region definitions come from metadata/configuration.
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


def generate_f11_domain_motion(
    domain_motion_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate domain motion figure.

    Parameters
    ----------
    domain_motion_data : list of domain motion result dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not domain_motion_data:
        return {
            "status": "skip",
            "reason": "No domain motion data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No domain motion data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(domain_motion_data)

    required_cols = ["region_label", "centroid_displacement"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Domain motion data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Domain motion data incomplete"],
        }

    df = df[df["centroid_displacement"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid centroid displacement values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All centroid displacement values are NaN"],
        }

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        df["plot_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["plot_label"] = df.get("condition_id", df.get("condition_a", "unknown"))

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition_id"].unique()]
        order += [c for c in df["condition_id"].unique() if c not in order]
        if not order:
            order = sorted(df["plot_label"].unique())
    else:
        order = sorted(df["plot_label"].unique())

    # Process region labels
    region_order = sorted(df["region_label"].unique())

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    palette = sns.color_palette("Set2", n_colors=max(len(region_order), 3))

    # Panel A: Centroid displacement by region
    ax = axes[0]

    # Prepare data for grouped boxplot
    plot_data = df.pivot_table(
        index="region_label",
        columns="plot_label",
        values="centroid_displacement",
        aggfunc="mean",
    )

    if not plot_data.empty:
        # Sort by region label
        plot_data = plot_data.reindex(sorted(plot_data.index))

        sns.boxplot(
            data=plot_data,
            palette=palette[:len(plot_data.columns)],
            linewidth=0.8,
            fliersize=0,
            boxprops=dict(alpha=0.4),
            ax=ax,
        )

        ax.set_xlabel("Region")
        ax.set_ylabel("Centroid displacement (Å)")
        ax.set_title("  A. Centroid Displacement by Region", loc="left", fontsize=11, fontweight="bold")
        ax.tick_params(axis="x", rotation=30, labelsize=9)
        ax.tick_params(axis="y", labelsize=9)
    else:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_title("  A. Centroid Displacement by Region", loc="left", fontsize=11, fontweight="bold")

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Inter-region distance changes
    ax = axes[1]

    if "distance_change" in df.columns:
        df_plot = df[df["distance_change"].notna()].copy()
        if not df_plot.empty:
            sns.boxplot(
                data=df_plot,
                x="region_label",
                y="distance_change",
                hue="region_label",
                palette=palette,
                linewidth=0.8,
                fliersize=0,
                boxprops=dict(alpha=0.4),
                legend=False,
                ax=ax,
            )
            sns.stripplot(
                data=df_plot,
                x="region_label",
                y="distance_change",
                hue="region_label",
                palette=palette,
                size=3,
                alpha=0.5,
                linewidth=0.3,
                edgecolor="white",
                jitter=0.2,
                legend=False,
                ax=ax,
            )

            ax.set_xlabel("Region")
            ax.set_ylabel("Inter-region distance change (Å)")
            ax.set_title("  B. Inter-region Distance Changes", loc="left", fontsize=11, fontweight="bold")
            ax.tick_params(axis="x", rotation=30, labelsize=9)
            ax.tick_params(axis="y", labelsize=9)

            sns.despine(ax=ax, left=True)
            ax.yaxis.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, "No distance change data", ha="center", va="center")
            ax.set_title("  B. Inter-region Distance Changes", loc="left", fontsize=11, fontweight="bold")
    else:
        ax.text(0.5, 0.5, "No distance data", ha="center", va="center")
        ax.set_title("  B. Inter-region Distance Changes", loc="left", fontsize=11, fontweight="bold")

    # Main title
    main_title = title or "Domain/Region Motion"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f11_domain_motion.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero domain motion observations")
    if len(region_order) > 10:
        warnings.append(f"Large number of regions ({len(region_order)}), figure may be dense")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
