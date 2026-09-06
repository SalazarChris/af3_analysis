"""
V3 Figure F13 — Structural Similarity Matrix

Shows pairwise structural distance matrix as a heatmap.
Preserves metadata for every prediction.
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


def generate_f13_similarity_matrix(
    distance_matrix: np.ndarray,
    predictions: List[str],
    conditions: List[str],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate structural similarity matrix figure.

    Parameters
    ----------
    distance_matrix : (N, N) array of RMSD values
    predictions : list of prediction_id
    conditions : list of condition_id
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    n = len(predictions)
    if n == 0:
        return {
            "status": "skip",
            "reason": "No structures",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No structural data available"],
        }

    if distance_matrix.shape[0] != n or distance_matrix.shape[1] != n:
        return {
            "status": "skip",
            "reason": "Distance matrix shape mismatch",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Distance matrix dimensions inconsistent"],
        }

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        labels = [f"{design.get_condition_label(c) if c in design.conditions else c}\n{p}"
                  for c, p in zip(conditions, predictions)]
    else:
        # Truncate prediction IDs for readability
        short_preds = [p[:25] + "..." if len(p) > 25 else p for p in predictions]
        labels = [f"{c}\n{p}" for c, p in zip(conditions, short_preds)]

    n_obs = n

    # --- Plot ---
    # Determine figure size based on number of structures
    fig_size = min(12.0, max(4.0, n * 0.15))
    fig, ax = plt.subplots(1, 1, figsize=(fig_size, fig_size))

    # Create annotation labels (only if not too many)
    annot = n <= 20

    sns.heatmap(
        distance_matrix,
        annot=annot,
        fmt=".2f",
        cmap="YlOrRd",
        vmin=0,
        ax=ax,
        cbar_kws={"label": "RMSD (Å)"},
    )

    ax.set_xlabel("Prediction index")
    ax.set_ylabel("Prediction index")
    ax.set_title(
        f"  {title or 'Structural Similarity Matrix (RMSD)'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=8)

    # Short tick labels
    if n > 20:
        # Show every Nth label
        tick_interval = max(1, n // 20)
        ticks = list(range(0, n, tick_interval))
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        short_labels = [f"{i}" for i in ticks]
        ax.set_xticklabels(short_labels, rotation=45, ha="right")
        ax.set_yticklabels(short_labels)
    else:
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
        ax.set_yticklabels(labels, fontsize=6)

    fig.tight_layout()

    out_path = save_path / "fig_v3_f13_similarity_matrix.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero similarity matrix observations")
    if n > 50:
        warnings.append(f"Large matrix ({n}×{n}), figure may be dense")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
