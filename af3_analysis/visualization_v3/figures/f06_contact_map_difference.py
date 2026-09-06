"""
V3 Figure F06 — Contact Map Difference

Shows Δcontact_map where:
+1 = gained
 0 = unchanged
-1 = lost
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


def generate_f06_contact_map_difference(
    contact_diff_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate contact map difference figure.

    Parameters
    ----------
    contact_diff_data : list of contact difference dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not contact_diff_data:
        return {
            "status": "skip",
            "reason": "No contact difference data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No contact map difference data available"],
        }

    # Build DataFrame from delta_map
    delta_entries = []
    for diff in contact_diff_data:
        delta_map = diff.get("delta_map", {})
        for pair, delta in delta_map.items():
            chain_a, seq_a, chain_b, seq_b = pair
            delta_entries.append({
                "condition": diff.get("condition", "unknown"),
                "chain_a": chain_a,
                "seq_a": seq_a,
                "chain_b": chain_b,
                "seq_b": seq_b,
                "delta": delta,  # +1, 0, -1
            })

    if not delta_entries:
        return {
            "status": "skip",
            "reason": "No delta entries",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No contact changes found"],
        }

    df = pd.DataFrame(delta_entries)

    # Summary statistics
    n_gained = sum(1 for d in delta_entries if d["delta"] == 1)
    n_lost = sum(1 for d in delta_entries if d["delta"] == -1)
    n_unchanged = sum(1 for d in delta_entries if d["delta"] == 0)
    n_total = len(delta_entries)
    pct_changed = (n_gained + n_lost) / n_total if n_total > 0 else 0.0

    n_obs = n_total

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 4.0))

    # Panel A: Contact change summary (bar chart)
    ax = axes[0]
    categories = ["Gained", "Lost", "Unchanged"]
    values = [n_gained, n_lost, n_unchanged]
    colors = ["#2C7BB6", "#D7191C", "#BBBBBB"]

    ax.bar(categories, values, color=colors, edgecolor="white", linewidth=0.8)
    ax.set_ylabel("Number of residue pairs")
    ax.set_title("  A. Contact Changes", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

    # Add percentage labels
    for i, (cat, val) in enumerate(zip(categories, values)):
        pct = val / n_total * 100 if n_total > 0 else 0
        ax.text(i, val + max(values) * 0.02, f"{pct:.1f}%",
                ha="center", va="bottom", fontsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Contact type distribution (pie)
    ax = axes[1]
    if n_gained + n_lost > 0:
        labels = ["Gained", "Lost"]
        sizes = [n_gained, n_lost]
        colors_pie = ["#2C7BB6", "#D7191C"]

        wedges, texts, autotexts = ax.pie(
            sizes,
            labels=labels,
            autopct="%1.1f%%",
            colors=colors_pie,
            startangle=90,
        )
        for t in autotexts:
            t.set_fontsize(9)
        for t in texts:
            t.set_fontsize(9)

    ax.set_title("  B. Changed Contacts", loc="left", fontsize=11, fontweight="bold")

    # Main title
    main_title = title or f"Contact Map Difference (Δ: {pct_changed:.1%} changed)"
    fig.suptitle(main_title, fontsize=12, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f06_contact_map_difference.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero contact observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
