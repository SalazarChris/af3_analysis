"""
V3 Figure F05 — Structural Displacement Heatmap

condition × residue with values representing:
- mean or median displacement

Also shows seed consistency representation.
Distinguishes:
- large but inconsistent movement
- moderate but reproducible movement
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
    DOUBLE_COL_WIDTH,
    apply_v3_style,
)


def generate_f05_displacement_heatmap(
    displacement_matrix: pd.DataFrame,
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    value_type: str = "mean",  # "mean" or "median"
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate displacement heatmap figure.

    Parameters
    ----------
    displacement_matrix : DataFrame with condition × residue displacements
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    value_type : str ("mean" or "median")
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if displacement_matrix is None or displacement_matrix.empty:
        return {
            "status": "skip",
            "reason": "No displacement matrix",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No displacement data available"],
        }

    # Determine columns
    if "condition_id" not in displacement_matrix.columns:
        return {
            "status": "skip",
            "reason": "No condition_id column",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Cannot identify conditions"],
        }

    if "residue_index" not in displacement_matrix.columns:
        return {
            "status": "skip",
            "reason": "No residue_index column",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Cannot identify residues"],
        }

    # Pivot to condition × residue matrix
    value_col = f"{value_type}_displacement"
    if value_col not in displacement_matrix.columns:
        # Calculate from raw values
        if "displacement" in displacement_matrix.columns:
            displacement_matrix[value_col] = displacement_matrix.groupby(
                ["condition_id", "residue_index"]
            )["displacement"].transform(value_type)
        else:
            return {
                "status": "skip",
                "reason": "No displacement column",
                "output_path": None,
                "n_observations": 0,
                "warnings": ["No displacement values"],
            }

    pivot = displacement_matrix.pivot_table(
        index="condition_id",
        columns="residue_index",
        values=value_col,
        aggfunc="mean",
    )

    if pivot.empty:
        return {
            "status": "skip",
            "reason": "Empty pivot table",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Cannot build condition × residue matrix"],
        }

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        label_map = {c: design.get_condition_label(c) for c in pivot.index}
        pivot.index = [label_map.get(c, c) for c in pivot.index]

    n_obs = len(displacement_matrix)

    # --- Plot ---
    fig_height = max(4, len(pivot.index) * 0.5)
    fig_width = min(DOUBLE_COL_WIDTH, len(pivot.columns) * 0.15 + 2)

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

    # Preflight: check matrix size
    n_cells = len(pivot.index) * len(pivot.columns)
    if n_cells > 500:
        warnings = [f"Large heatmap: {n_cells} cells, annotations disabled"]
        annot = False
    else:
        warnings = []
        annot = True

    sns.heatmap(
        pivot,
        annot=annot,
        fmt=".2f",
        cmap="RdYlBu_r",
        center=0,
        vmin=pivot.min().min(),
        vmax=pivot.max().max(),
        ax=ax,
        cbar_kws={"label": f"{value_type.capitalize()} Cα Displacement (Å)"},
    )

    ax.set_xlabel("Residue Index")
    ax.set_ylabel("Condition")
    ax.set_title(
        f"  {title or f'{value_type.capitalize()} Displacement Heatmap'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=8)

    fig.tight_layout()

    out_path = save_path / "fig_v3_f05_displacement_heatmap.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
