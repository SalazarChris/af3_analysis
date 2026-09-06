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

    sns.scatterplot(
        data=df,
        x="x",
        y="y",
        hue="condition",
        palette=palette,
        s=80,
        alpha=0.7,
        linewidth=0.3,
        edgecolor="white",
        ax=ax,
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

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 12:
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8, title="Condition")
    else:
        ax.legend(loc="best", fontsize=9, title="Condition")

    # Main title
    fig.suptitle(
        "Do different conditions occupy distinguishable regions of predicted structural space?",
        fontsize=10,
        fontstyle="italic",
        y=0.98,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f14_mds_embedding.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero MDS observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
