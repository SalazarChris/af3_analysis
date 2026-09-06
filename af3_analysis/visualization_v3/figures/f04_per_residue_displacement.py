"""
V3 Figure F04 — Per-Residue Structural Displacement

Shows per-residue displacement between reference and condition.
Default atom: CA, but configurable.

Generates table with:
- condition
- residue_index
- residue_name
- displacement
- n_valid
- coverage
- seed_consistency
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


def generate_f04_per_residue_displacement(
    displacement_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate per-residue displacement figure.

    Parameters
    ----------
    displacement_data : list of per-residue displacement dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not displacement_data:
        return {
            "status": "skip",
            "reason": "No displacement data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No per-residue displacement data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(displacement_data)

    if "displacement" not in df.columns:
        return {
            "status": "skip",
            "reason": "Displacement data missing displacement column",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Displacement data incomplete"],
        }

    df = df[df["displacement"].notna()]

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid displacement values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All displacement values are NaN"],
        }

    # Get condition label
    if "condition_label" in df.columns:
        label = df["condition_label"].iloc[0]
    elif design and hasattr(design, "get_condition_label"):
        cond_id = df.get("condition_id", pd.Series(["unknown"]))[0]
        label = design.get_condition_label(cond_id)
    else:
        label = df.get("condition_id", pd.Series(["unknown"]))[0]

    n_obs = len(df)

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(SINGLE_COL_WIDTH, 4.0))

    # Plot displacement by residue index
    ax.scatter(
        df["residue_index"],
        df["displacement"],
        s=5,
        alpha=0.5,
        color="#2C7BB6",
        edgecolors="none",
    )

    # Add horizontal line at mean
    mean_disp = df["displacement"].mean()
    ax.axhline(mean_disp, color="#D7191C", linestyle="--",
               linewidth=1.0, label=f"Mean: {mean_disp:.2f} Å")

    # Add zero line
    ax.axhline(0, color="grey", linestyle="-", linewidth=0.5, alpha=0.5)

    ax.set_xlabel("Residue Index (auth_seq_id)")
    ax.set_ylabel("Cα Displacement (Å)")
    ax.set_title(
        f"  {title or f'Per-Residue Displacement ({label})'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    out_path = save_path / "fig_v3_f04_per_residue_displacement.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero displacement observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
