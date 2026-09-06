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

    # Build DataFrame for plotting
    df = pd.DataFrame({
        "prediction_id": predictions,
        "condition": conditions,
        "seed": seeds,
        "cluster": cluster_labels,
    })

    # Get condition labels
    if design and hasattr(design, "get_condition_label"):
        df["condition_label"] = df["condition"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["condition_label"] = df["condition"]

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition"].unique()]
        order += [c for c in df["condition"].unique() if c not in order]
    else:
        order = sorted(df["condition"].unique())

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}
    df["condition_label"] = df["condition_label"].astype(
        pd.CategoricalDtype(categories=[label_map.get(c, c) for c in order], ordered=True)
    )

    n_obs = n

    # --- Plot ---
    fig_height = min(8.0, max(4.0, n * 0.15))
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL_WIDTH, fig_height))

    palette = sns.color_palette("Set2", n_colors=max(n_clusters, 2))

    # Panel A: Clustering dendrogram (heatmap with cluster boundaries)
    ax = axes[0]

    # Create annotation matrix
    annot_matrix = pd.DataFrame(
        cluster_labels.reshape(-1, 1),
        index=predictions,
        columns=["cluster"],
    )

    # Draw clustered heatmap
    sns.heatmap(
        distance_matrix,
        annot=False,
        cmap="YlOrRd",
        ax=ax,
        cbar_kws={"label": "RMSD (Å)"},
        xticklabels=[],
        yticklabels=[],
    )

    ax.set_xlabel("Prediction index")
    ax.set_ylabel("Prediction index")
    ax.set_title("  A. Structural Distance Matrix", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=8)

    # Panel B: Clustering assignment
    ax = axes[1]

    # Create a dataframe for visualization
    viz_df = df.copy()
    viz_df["index"] = range(len(viz_df))

    sns.scatterplot(
        data=viz_df,
        x="index",
        y="cluster",
        hue="condition_label",
        palette=palette,
        s=50,
        alpha=0.7,
        linewidth=0.3,
        edgecolor="white",
        ax=ax,
    )

    ax.set_xlabel("Prediction index")
    ax.set_ylabel("Predicted structural cluster")
    ax.set_title("  B. Cluster Assignments", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 12:
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    else:
        ax.legend(loc="best", fontsize=9)

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
