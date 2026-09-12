"""
V3 Figure F12 — Structural Clustering

Uses pairwise structural distance matrix.
Shows predicted structural clusters using hierarchical clustering.

Do NOT claim that clusters represent biological states.
Use terminology: "predicted structural cluster"
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


def generate_f12_structural_clustering(
    distance_matrix: np.ndarray,
    conditions: List[str],
    seeds: List[int],
    predictions: List[str],
    save_path: Path,
    *,
    cluster_labels: Optional[np.ndarray] = None,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate structural clustering figure.

    Parameters
    ----------
    distance_matrix : (N, N) array of RMSD values
    conditions : list of condition_id
    seeds : list of seed
    predictions : list of prediction_id
    cluster_labels : (N,) array of cluster assignments (optional)
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

    # Generate cluster labels if not provided
    if cluster_labels is None:
        from sklearn.cluster import AgglomerativeClustering
        from scipy.spatial.distance import squareform

        n_clusters = max(2, min(5, int(np.sqrt(n))))
        try:
            condensed = squareform(distance_matrix)
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                linkage="average",
            )
            clustering.fit(condensed)
            cluster_labels = clustering.labels_
        except Exception:
            cluster_labels = np.zeros(n, dtype=int)

    n_clusters = len(np.unique(cluster_labels))

    n_obs = n

    # Reorder predictions by predicted structural cluster (then condition
    # and seed) so block structure in the distance matrix is visible.
    # This is a presentation ordering only; it does not alter the
    # clustering or the underlying distances.
    ordered_cluster_labels = np.asarray(cluster_labels).astype(int)
    row_order = sorted(
        range(n),
        key=lambda i: (
            int(ordered_cluster_labels[i]), str(conditions[i]), int(seeds[i]),
            str(predictions[i]),
        ),
    )
    ordered_matrix = distance_matrix[np.ix_(row_order, row_order)]
    ordered_predictions = [predictions[i] for i in row_order]
    ordered_conditions = [conditions[i] for i in row_order]
    ordered_cluster_labels = ordered_cluster_labels[row_order]

    # --- Plot ---
    fig_height = min(8.0, max(4.0, n * 0.15))
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL_WIDTH, fig_height))

    # One consistent color per predicted structural cluster, used in both
    # panels so assignments can be read across the figure.
    cluster_ids = sorted(np.unique(ordered_cluster_labels).tolist())
    cluster_palette_list = sns.color_palette("Set2", n_colors=max(len(cluster_ids), 2))
    cluster_palette = {
        f"Cluster {cid + 1}": cluster_palette_list[i]
        for i, cid in enumerate(cluster_ids)
    }

    # Panel A: distance matrix ordered by predicted structural cluster.
    # Perceptually uniform colormap for the continuous distances (display
    # choice only; matrix values, vmin=0, and no clipping are unchanged).
    ax = axes[0]

    sns.heatmap(
        ordered_matrix,
        annot=False,
        cmap="viridis",
        vmin=0,
        ax=ax,
        cbar_kws={"label": "RMSD (Å)", "fraction": 0.04, "pad": 0.02},
        xticklabels=False,
        yticklabels=False,
    )

    # Draw white boundaries between consecutive cluster blocks.
    for i in range(1, n):
        if ordered_cluster_labels[i] != ordered_cluster_labels[i - 1]:
            ax.axhline(i, color="white", linewidth=1.2)
            ax.axvline(i, color="white", linewidth=1.2)

    ax.set_xlabel("Prediction index (ordered by predicted cluster)")
    ax.set_ylabel("Prediction index (ordered by predicted cluster)")
    ax.set_title("  A. Structural Distance Matrix", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=8)

    # Panel B: cluster assignment as a compact ordered strip. The y-axis
    # already identifies each predicted structural cluster, so points are
    # NOT colored by cluster (no redundant 20-color legend); condition
    # identity is used for coloring instead, consistent with the other
    # figures' condition palette.
    ax = axes[1]

    viz_df = pd.DataFrame({
        "prediction_id": ordered_predictions,
        "condition": ordered_conditions,
        "cluster": ordered_cluster_labels,
        "index": range(n),
    })
    if design and hasattr(design, "get_condition_label"):
        viz_df["condition_label"] = viz_df["condition"].map(
            lambda c: design.get_condition_label(c)
            if c in design.conditions else c
        )
    else:
        viz_df["condition_label"] = viz_df["condition"].astype(str)

    cond_order = list(dict.fromkeys(viz_df["condition_label"]))
    cond_palette_list = sns.color_palette("Set2", n_colors=max(len(cond_order), 3))
    cond_palette = dict(zip(cond_order, cond_palette_list))

    for cond_label, sub in viz_df.groupby("condition_label", observed=True):
        ax.scatter(
            sub["index"], sub["cluster"],
            s=16, marker="s", alpha=0.85,
            color=cond_palette.get(cond_label, "#555555"),
            edgecolors="none", label=cond_label,
        )

    ax.set_yticks(cluster_ids)
    if len(cluster_ids) > 15:
        # Keep tick text readable on dense assignments.
        ax.set_yticklabels(
            [str(cid + 1) if cid % 2 == 0 else "" for cid in cluster_ids],
            fontsize=8,
        )
    else:
        ax.set_yticklabels([f"{cid + 1}" for cid in cluster_ids], fontsize=8)
    ax.set_ylim(min(cluster_ids) - 0.5, max(cluster_ids) + 0.5)
    ax.set_xlim(-1, n)
    ax.set_xlabel("Prediction index (ordered by predicted cluster)", fontsize=9)
    ax.set_ylabel("Predicted structural cluster", fontsize=9)
    ax.set_title("  B. Cluster Assignments", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", labelsize=8)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Legend: conditions only (cluster identity is on the y-axis). Omitted
    # entirely for a single condition.
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 1:
        ax.legend(
            loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7,
            title="Condition",
        )

    # Main title
    main_title = title or f"Structural Clustering ({n_clusters} clusters)"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f12_structural_clustering.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero clustering observations")
    if n_clusters > 10:
        warnings.append(f"Large number of clusters ({n_clusters})")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
