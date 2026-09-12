"""
V3 Figure F14 — Structural MDS/PCA Embedding

Generates 2D structural-space representation using MDS.
Color/group using metadata: condition, factor, DNA presence, PTM state, seed.

Do NOT interpret axes as biological mechanisms.
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


def generate_f14_mds_embedding(
    coordinates: np.ndarray,
    predictions: List[str],
    conditions: List[str],
    seeds: List[int],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate MDS embedding figure.

    Parameters
    ----------
    coordinates : (N, 2) array from MDS embedding
    predictions : list of prediction_id
    conditions : list of condition_id
    seeds : list of seed
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

    if coordinates.shape[0] != n:
        return {
            "status": "skip",
            "reason": "Coordinates shape mismatch",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Coordinates dimensions inconsistent"],
        }

    if coordinates.shape[1] < 2:
        return {
            "status": "skip",
            "reason": "Need at least 2 dimensions for embedding",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Insufficient embedding dimensions"],
        }

    # Build DataFrame
    df = pd.DataFrame({
        "prediction_id": predictions,
        "condition": conditions,
        "seed": seeds,
        "x": coordinates[:, 0],
        "y": coordinates[:, 1],
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

    n_obs = n

    # --- Plot ---
    fig, ax = plt.subplots(1, 1, figsize=(SINGLE_COL_WIDTH, 6.0))

    palette = sns.color_palette("Set2", n_colors=max(len(order), 2))

    # Small, transparent markers reduce overplotting while keeping every
    # point visible. This is a presentation change only; the embedding
    # coordinates are not altered.
    sns.scatterplot(
        data=df,
        x="x",
        y="y",
        hue="condition",
        palette=palette,
        s=22,
        alpha=0.55,
        linewidth=0.0,
        ax=ax,
    )

    # Condition centroids derived from the existing coordinates only.
    # They summarize where each existing condition sits in the displayed
    # projection; they do not add a new embedding, clustering, or metric.
    centroid_rows = []
    for cond_id, grp in df.groupby("condition"):
        label = label_map.get(cond_id, cond_id)
        centroid_rows.append({
            "condition": cond_id,
            "condition_label": label,
            "x": float(grp["x"].mean()),
            "y": float(grp["y"].mean()),
            "n": int(grp.shape[0]),
        })
    centroid_df = pd.DataFrame(centroid_rows)
    if not centroid_df.empty:
        sns.scatterplot(
            data=centroid_df,
            x="x",
            y="y",
            hue="condition_label",
            palette=palette,
            s=70,
            alpha=0.95,
            linewidth=1.0,
            edgecolor="white",
            marker="D",
            ax=ax,
            legend=False,
        )
        for _, row in centroid_df.iterrows():
            ax.annotate(
                f"{row['condition_label']} (n={row['n']})",
                (row["x"], row["y"]),
                xytext=(6, 4),
                textcoords="offset points",
                fontsize=7,
                color="#222222",
                annotation_clip=True,
            )

    ax.set_xlabel("MDS Dimension 1")
    ax.set_ylabel("MDS Dimension 2")
    ax.set_title(
        f"  {title or 'Structural Space (MDS Embedding)'}",
        loc="left",
        fontsize=12,
        fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.grid(True, alpha=0.3)

    # Legend: keep condition legend compact; do not duplicate centroid
    # labels inside it.
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 12:
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, title="Condition")
    else:
        ax.legend(loc="best", fontsize=9, title="Condition")

    # Main title plus an always-on caveat about what projection overlap
    # does and does not mean (review readability expectation: overlap in
    # the projection is not evidence of structural identity).
    fig.suptitle(
        "Do different conditions occupy distinguishable regions of predicted structural space?",
        fontsize=10,
        fontstyle="italic",
        y=0.98,
    )
    ax.text(
        0.5, -0.12,
        "Overlap in this 2D projection is not evidence of structural identity; "
        "separation is not evidence of difference.",
        transform=ax.transAxes,
        ha="center", va="top", fontsize=7,
        style="italic", color="#444444",
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f14_mds_embedding.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero MDS observations")
    if n_obs > 1:
        point_span = _point_cloud_span(df)
        if _mds_separation_note(centroid_df, point_span) is not None:
            warnings.append(_mds_separation_note(centroid_df, point_span))

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }


def _point_cloud_span(df: pd.DataFrame) -> float:
    """Max spread of all displayed points in the projection."""
    if df.empty:
        return 0.0
    return max(
        float(np.ptp(df["x"].to_numpy())),
        float(np.ptp(df["y"].to_numpy())),
    )


def _condition_centroid_span(centroid_df: pd.DataFrame) -> Optional[float]:
    """Max spread of condition centroids, or None when not assessable.

    Spread is measured on the displayed projection only. A single condition
    (one centroid) or coincident centroids return None: nothing meaningful
    to report.
    """
    if centroid_df.shape[0] < 2:
        return None
    # Series.ptp() was removed in pandas 2.0; use np.ptp on arrays.
    ptp_x = float(np.ptp(centroid_df["x"].to_numpy()))
    ptp_y = float(np.ptp(centroid_df["y"].to_numpy()))
    span = max(ptp_x, ptp_y)
    if span <= 0:
        return None
    return span


def _mds_separation_note(
    centroid_df: pd.DataFrame,
    point_span: float,
) -> Optional[str]:
    """Conservative interpretation note for crowded MDS projections.

    Returns a note only when condition centroids occupy a limited region of
    the displayed projection relative to the overall point spread
    (centroid span at most half the point-cloud span). This ratio is a
    presentation heuristic for when the caveat caption applies; it imposes
    no threshold on the underlying structural data. Overlap in this 2D
    projection is never evidence of structural identity.
    """
    if point_span <= 0:
        return None
    span = _condition_centroid_span(centroid_df)
    if span is None:
        return None
    if span / point_span > 0.5:
        return None
    return (
        "Condition centroids occupy a limited region of the displayed 2D "
        "projection relative to the overall spread; overlap in this "
        "projection is not evidence of structural identity."
    )
