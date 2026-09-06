"""
V3 Figure F09 — Interface Change Map

Shows which interfaces are gained, lost, or unchanged between conditions.
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


def generate_f09_interface_change_map(
    interface_changes: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate interface change map figure.

    Parameters
    ----------
    interface_changes : list of interface change dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not interface_changes:
        return {
            "status": "skip",
            "reason": "No interface change data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No interface change data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(interface_changes)

    if "interface_id" not in df.columns or "change" not in df.columns:
        return {
            "status": "skip",
            "reason": "Interface change data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Interface change data incomplete"],
        }

    # Count changes by type
    change_counts = df["change"].value_counts()

    n_gained = change_counts.get("gained", 0)
    n_lost = change_counts.get("lost", 0)
    n_unchanged = change_counts.get("unchanged", 0)
    n_increased = change_counts.get("increased", 0)
    n_decreased = change_counts.get("decreased", 0)
    n_total = len(df)

    n_obs = n_total

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 4.0))

    # Panel A: Change summary (stacked bar)
    ax = axes[0]

    categories = ["Gained", "Lost", "Increased", "Decreased", "Unchanged"]
    values = [n_gained, n_lost, n_increased, n_decreased, n_unchanged]
    colors = ["#2C7BB6", "#D7191C", "#FDAE61", "#4C72B0", "#BBBBBB"]

    # Filter to non-zero
    non_zero_idx = [i for i, v in enumerate(values) if v > 0]
    if non_zero_idx:
        cats = [categories[i] for i in non_zero_idx]
        vals = [values[i] for i in non_zero_idx]
        cols = [colors[i] for i in non_zero_idx]

        ax.bar(cats, vals, color=cols, edgecolor="white", linewidth=0.8)
        ax.set_ylabel("Number of interfaces")
        ax.set_title("  A. Interface Changes", loc="left", fontsize=11, fontweight="bold")
        ax.tick_params(axis="both", labelsize=9)

        # Add count labels
        for i, v in enumerate(vals):
            ax.text(i, v + max(vals) * 0.02, str(v), ha="center", va="bottom", fontsize=9)
    else:
        ax.text(0.5, 0.5, "No changes", ha="center", va="center")

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Change distribution (pie)
    ax = axes[1]

    changed_cats = ["Gained/Lost", "Unchanged"]
    changed_vals = [n_gained + n_lost + n_increased + n_decreased, n_unchanged]
    changed_colors = ["#FDAE61", "#BBBBBB"]

    if sum(changed_vals) > 0:
        wedges, texts, autotexts = ax.pie(
            changed_vals,
            labels=changed_cats,
            autopct="%1.1f%%",
            colors=changed_colors,
            startangle=90,
        )
        for t in autotexts:
            t.set_fontsize(9)
        for t in texts:
            t.set_fontsize(9)

    ax.set_title("  B. Changed vs Unchanged", loc="left", fontsize=11, fontweight="bold")

    # Main title
    pct_changed = (n_gained + n_lost + n_increased + n_decreased) / n_total if n_total > 0 else 0
    main_title = title or f"Interface Change Map (Δ: {pct_changed:.1%} changed)"
    fig.suptitle(main_title, fontsize=12, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f09_interface_change_map.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero interface change observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
