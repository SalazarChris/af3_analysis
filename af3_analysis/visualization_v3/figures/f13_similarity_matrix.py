"""
V3 Figure F13 — Structural Similarity Matrix

Shows pairwise structural distance matrix as a heatmap.
Preserves metadata for every prediction.

This is the primary direct pairwise structural-comparison figure.
Presentation (visualization only; no analysis change):
- optional row/column ordering by EXISTING predicted structural cluster
  assignments when provided (no recomputation of clustering);
- condition-boundary lines using existing condition metadata;
- sparse index ticks so the matrix stays readable at prediction scale;
- compact method note documenting the existing comparison basis.
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
    cluster_labels: Optional[np.ndarray] = None,
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
    cluster_labels : optional (N,) array of EXISTING predicted structural
        cluster assignments (e.g. from the F12 clustering step). Used only
        for display ordering; never recomputed here.
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

    # --- Display ordering (existing data only) ---
    # Prefer existing cluster assignments (from the runner's F12 step);
    # otherwise fall back to condition identity, then prediction id.
    if (
        cluster_labels is not None
        and len(cluster_labels) == n
        and len(np.unique(cluster_labels)) > 1
    ):
        order_key = sorted(
            range(n),
            key=lambda i: (int(cluster_labels[i]), str(conditions[i]),
                           str(predictions[i])),
        )
        ordering_note = "rows/columns ordered by existing predicted structural clusters"
    else:
        order_key = sorted(
            range(n),
            key=lambda i: (str(conditions[i]), str(predictions[i])),
        )
        ordering_note = "rows/columns ordered by condition identity"

    ordered_matrix = distance_matrix[np.ix_(order_key, order_key)]
    ordered_conditions = [conditions[i] for i in order_key]

    # Condition boundaries in the ordered layout (existing metadata only).
    boundaries = []
    for i in range(1, n):
        if ordered_conditions[i] != ordered_conditions[i - 1]:
            boundaries.append(i)

    n_obs = n

    # --- Plot ---
    fig_size = min(12.0, max(4.0, n * 0.15))
    fig, ax = plt.subplots(1, 1, figsize=(fig_size, fig_size))

    # Annotation labels only for small matrices (existing behavior).
    annot = n <= 20

    sns.heatmap(
        ordered_matrix,
        annot=annot,
        fmt=".2f",
        cmap="viridis",
        vmin=0,
        ax=ax,
        cbar_kws={"label": "RMSD (Å)", "fraction": 0.04, "pad": 0.02},
    )

    # Cluster boundaries (when cluster ordering is active) and condition
    # boundaries (always) as thin separators.
    if (
        cluster_labels is not None
        and len(cluster_labels) == n
        and len(np.unique(cluster_labels)) > 1
    ):
        ordered_clusters = [int(cluster_labels[i]) for i in order_key]
        for i in range(1, n):
            if ordered_clusters[i] != ordered_clusters[i - 1]:
                ax.axhline(i, color="white", linewidth=0.8, alpha=0.9)
                ax.axvline(i, color="white", linewidth=0.8, alpha=0.9)
    for i in boundaries:
        ax.axhline(i, color="black", linewidth=0.7, alpha=0.6)
        ax.axvline(i, color="black", linewidth=0.7, alpha=0.6)

    ax.set_xlabel("Prediction index")
    ax.set_ylabel("Prediction index")
    ax.set_title(
        f"  {title or 'Structural Similarity Matrix (RMSD)'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=8)

    # Compact structural-comparison method note.
    # This documents the existing comparison basis without inventing
    # methodology; any field not confirmed here should stay marked as
    # not confirmed rather than inferred.
    method_lines = [
        "Structural comparison:",
        "Metric: pairwise RMSD",
        "Design: all-vs-all over comparable predictions",
    ]
    if reference_condition:
        method_lines.append(
            f"Reference: {reference_condition} "
            f"(prediction-to-reference comparisons where applicable)"
        )
    method_lines.append("Distance units: Å")
    method_lines.append(
        "Comparison basis: existing pairwise RMSD implementation "
        "(see pipeline structural modules)"
    )
    method_lines.append(ordering_note)
    method_note = "\n".join(method_lines)

    if n <= 50:
        ax.text(
            0.5, -0.24, method_note,
            transform=ax.transAxes,
            ha="center", va="top",
            fontsize=6.0,
            style="italic",
            color="#444444",
        )

    # Sparse index ticks (unchanged numbering; display subsampling only).
    tick_interval = max(1, n // 20)
    ticks = list(range(0, n, tick_interval))
    ax.set_xticks([t + 0.5 for t in ticks])
    ax.set_yticks([t + 0.5 for t in ticks])
    tick_labels = [str(t) for t in ticks]
    ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=6)
    ax.set_yticklabels(tick_labels, fontsize=6)

    fig.tight_layout()

    out_path = save_path / "fig_v3_f13_similarity_matrix.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero similarity matrix observations")
    if n > 50:
        warnings.append(f"Large matrix ({n}×{n}), figure may be dense")
    if n > 20 and reference_condition:
        warnings.append(
            "Reference-condition-dependent comparisons use the configured "
            "reference; seed-level conclusions should use seed-matched analyses"
        )
    if cluster_labels is not None and len(cluster_labels) == n:
        warnings.append(
            "Rows/columns ordered by existing predicted structural cluster "
            "assignments (display ordering only; no reclustering performed)"
        )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
