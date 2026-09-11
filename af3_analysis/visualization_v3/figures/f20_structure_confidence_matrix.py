"""
V3 Figure F20 — Structure-Confidence Relationship Matrix

Shows relationships between multiple structural metrics and confidence metrics.
Combines:
- structural geometry (RMSD, displacement)
- AF3 confidence (pLDDT, PAE, contact_prob)

Keeps metrics conceptually separate. Does NOT merge into one arbitrary score.
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
from scipy import stats

from af3_analysis.visualization_v3.config import (
    DPI,
    SINGLE_COL_WIDTH,
    apply_v3_style,
)


def generate_f20_structure_confidence_matrix(
    relationship_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate structure-confidence relationship matrix figure.

    Parameters
    ----------
    relationship_data : list of relationship observation dicts
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not relationship_data:
        return {
            "status": "skip",
            "reason": "No relationship data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No structure-confidence relationship data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(relationship_data)

    # Identify structural and confidence metrics
    structural_metrics = [col for col in df.columns
                         if col.startswith("rmsd") or col.startswith("displacement")
                         or col.startswith("contact_change") or col.startswith("interface_change")]
    confidence_metrics = [col for col in df.columns
                         if col.startswith("plddt") or col.startswith("pae")
                         or col.startswith("contact_prob") or col.startswith("ranking")]

    if not structural_metrics or not confidence_metrics:
        return {
            "status": "skip",
            "reason": "Relationship data missing metric columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Relationship data incomplete — need both structural and confidence metrics"],
        }

    # Filter to valid rows
    df = df[structural_metrics + confidence_metrics].dropna(how="any").copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid relationship values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All relationship values are NaN"],
        }

    n_obs = len(df)

    # --- Plot ---
    n_struct = len(structural_metrics)
    n_conf = len(confidence_metrics)
    fig_height = min(10.0, max(4.0, max(n_struct, n_conf) * 0.5))

    fig, axes = plt.subplots(n_struct, n_conf, figsize=(SINGLE_COL_WIDTH * 1.5, fig_height))

    # Handle single row/column cases
    if n_struct == 1 and n_conf == 1:
        axes = np.array([[axes]])
    elif n_struct == 1:
        axes = axes.reshape(1, -1)
    elif n_conf == 1:
        axes = axes.reshape(-1, 1)

    for i, struct_metric in enumerate(structural_metrics):
        for j, conf_metric in enumerate(confidence_metrics):
            ax = axes[i, j]

            x_vals = df[conf_metric]
            y_vals = df[struct_metric]

            # Scatter
            ax.scatter(
                x_vals,
                y_vals,
                s=20,
                alpha=0.5,
                color="#2C7BB6",
                edgecolors="white",
                linewidth=0.3,
            )

            # Correlation
            if len(x_vals) > 2 and np.ptp(np.asarray(x_vals)) > 0:
                corr, p_val = stats.pearsonr(x_vals, y_vals)
                ax.text(
                    0.05, 0.95,
                    f"r = {corr:.2f}\np = {p_val:.3f}",
                    transform=ax.transAxes,
                    fontsize=8,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                )

                # Add regression line
                slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
                x_line = np.linspace(x_vals.min(), x_vals.max(), 100)
                y_line = slope * x_line + intercept
                ax.plot(x_line, y_line, color="red", linewidth=1, linestyle="--", alpha=0.5)
            elif len(x_vals) > 2:
                # Zero-variance metric (e.g. constant pLDDT): correlation
                # and regression are undefined. Annotate instead of crashing.
                ax.text(
                    0.05, 0.95,
                    "r = n/a\n(constant x)",
                    transform=ax.transAxes,
                    fontsize=8,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                )

            ax.set_xlabel(conf_metric, fontsize=8)
            ax.set_ylabel(struct_metric, fontsize=8)
            ax.tick_params(axis="both", labelsize=8)

            sns.despine(ax=ax, left=True)
            ax.grid(True, alpha=0.3)

    # Main title
    main_title = title or "Structure-Confidence Relationship Matrix"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f20_structure_confidence_matrix.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero relationship observations")
    if n_struct * n_conf > 20:
        warnings.append(f"Large relationship matrix ({n_struct}×{n_conf}), figure may be dense")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
