"""
V3 Figure F07 — Contact Change Summary

For each condition vs reference:
- contacts gained
- contacts lost
- percentage changed
- seed consistency
- recurring residue-pair changes

Generates summary table.
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


def generate_f07_contact_change_summary(
    contact_change_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate contact change summary figure.

    Parameters
    ----------
    contact_change_data : list of contact change summary dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not contact_change_data:
        return {
            "status": "skip",
            "reason": "No contact change data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No contact change data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(contact_change_data)

    required_cols = ["condition_id", "n_gained", "n_lost", "n_total"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Contact change data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Contact change data incomplete"],
        }

    df = df[df["n_total"] > 0].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid contact data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All contact totals are zero"],
        }

    # Calculate percentage changed
    df["pct_changed"] = (df["n_gained"] + df["n_lost"]) / df["n_total"]

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        df["plot_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["plot_label"] = df["condition_id"]

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition_id"].unique()]
        order += [c for c in df["condition_id"].unique() if c not in order]
    else:
        order = sorted(df["condition_id"].unique())

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}
    categories = [label_map.get(c, c) for c in order]
    df["plot_label"] = pd.Categorical(df["plot_label"], categories=categories, ordered=True)

    n_obs = len(df)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    palette = sns.color_palette("Set2", n_colors=max(len(order), 3))
    color_map = {label_map.get(c, c): palette[i] for i, c in enumerate(order)}

    # Panel A: Gained vs Lost contacts (grouped bar)
    ax = axes[0]
    x_pos = np.arange(len(df))
    width = 0.35

    ax.bar(
        x_pos - width / 2,
        df["n_gained"],
        width,
        label="Gained",
        color="#2C7BB6",
        edgecolor="white",
    )
    ax.bar(
        x_pos + width / 2,
        df["n_lost"],
        width,
        label="Lost",
        color="#D7191C",
        edgecolor="white",
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(df["plot_label"], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Number of residue pairs")
    ax.set_title("  A. Contact Changes per Condition", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Percentage changed (bar)
    ax = axes[1]

    sns.barplot(
        data=df,
        x="plot_label",
        y="pct_changed",
        hue="plot_label",
        palette=color_map,
        linewidth=0.8,
        edgecolor="white",
        legend=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Fraction of contacts changed")
    ax.set_title("  B. Contact Change Fraction", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Main title
    main_title = title or "Contact Change Summary"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f07_contact_change_summary.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero contact change observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
