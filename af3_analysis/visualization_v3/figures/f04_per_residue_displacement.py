"""
V3 Figure F04 — Per-Residue Structural Displacement

Shows per-residue displacement between reference and condition.
Default atom: CA, but configurable.

Presentation (visualization only; no analysis change):
- one small-multiple panel per condition when several conditions are
  present, sharing the y-axis so amplitudes are comparable;
- raw prediction-level observations at reduced marker size/opacity;
- a per-residue median line with an IQR band summarizing the SAME
  per-residue observations already plotted (descriptive display
  aggregation of existing values; no new metric);
- the existing global mean line, labeled.

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
from matplotlib.ticker import MaxNLocator

from af3_analysis.visualization_v3.config import (
    DPI,
    DOUBLE_COL_WIDTH,
    SINGLE_COL_WIDTH,
    apply_v3_style,
)


def _condition_order(df: pd.DataFrame, design: Optional[Any]) -> List[Any]:
    """Condition order from existing metadata where available."""
    conds = list(df["condition_id"].dropna().unique())
    if design and hasattr(design, "condition_names"):
        known = [c for c in design.condition_names if c in conds]
        known += [c for c in conds if c not in known]
        return known
    return sorted(conds)


def _condition_labels(
    order: List[Any],
    df: pd.DataFrame,
    design: Optional[Any],
) -> Dict[Any, str]:
    """Display labels from existing design metadata; IDs otherwise."""
    if design and hasattr(design, "get_condition_label"):
        return {
            c: (design.get_condition_label(c)
                if c in getattr(design, "conditions", []) else str(c))
            for c in order
        }
    # Fall back to any existing per-row label column (already provided by
    # the runner in some paths).
    if "condition_label" in df.columns:
        mapping = (
            df.dropna(subset=["condition_label"])
            .drop_duplicates("condition_id")
            .set_index("condition_id")["condition_label"]
            .to_dict()
        )
        return {c: str(mapping.get(c, c)) for c in order}
    return {c: str(c) for c in order}


def _per_residue_summary(group: pd.DataFrame) -> pd.DataFrame:
    """Per-residue median and IQR of the plotted displacement values.

    Pure descriptive summary of the observations already in the panel;
    no filtering, no new metric.
    """
    return (
        group.groupby("residue_index")["displacement"]
        .agg(median="median", q25=lambda s: s.quantile(0.25),
             q75=lambda s: s.quantile(0.75))
        .reset_index()
        .sort_values("residue_index")
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

    df = df[df["displacement"].notna()].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid displacement values",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All displacement values are NaN"],
        }

    if "condition_id" not in df.columns:
        df["condition_id"] = "unknown"

    n_obs = len(df)

    order = _condition_order(df, design)
    label_map = _condition_labels(order, df, design)

    # --- Plot: small multiples, one panel per condition, shared y-axis ---
    n_cond = len(order)
    fig_height = min(10.0, max(4.0, 2.2 * n_cond))
    fig_width = min(DOUBLE_COL_WIDTH, SINGLE_COL_WIDTH * 1.6)
    fig, axes = plt.subplots(
        n_cond, 1,
        figsize=(fig_width, fig_height),
        sharex=True, sharey=True, squeeze=False,
    )

    x_min = float(df["residue_index"].min())
    x_max = float(df["residue_index"].max())

    warnings = []
    for ax, cond in zip(axes[:, 0], order):
        group = df[df["condition_id"] == cond]
        if group.empty:
            ax.text(0.5, 0.5, "No observations", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9)
            ax.set_ylabel(str(label_map[cond]), fontsize=9)
            continue

        # Raw prediction-level observations: small, semi-transparent so
        # dense regions remain interpretable.
        ax.scatter(
            group["residue_index"],
            group["displacement"],
            s=4,
            alpha=0.15,
            color="#2C7BB6",
            edgecolors="none",
            rasterized=True,
        )

        # Per-residue median with IQR band: descriptive summary of the
        # SAME observations plotted above (visualization-only aggregation).
        summary = _per_residue_summary(group)
        if not summary.empty:
            ax.fill_between(
                summary["residue_index"],
                summary["q25"], summary["q75"],
                color="#2C7BB6", alpha=0.25, linewidth=0,
                label="Per-residue IQR",
            )
            ax.plot(
                summary["residue_index"], summary["median"],
                color="#2C7BB6", linewidth=1.1, label="Per-residue median",
            )

        # Existing global mean line for this condition's data (kept from
        # the original figure; display-only reference, not a threshold).
        mean_disp = group["displacement"].mean()
        ax.axhline(
            mean_disp, color="#D7191C", linestyle="--", linewidth=1.0,
            label=f"Mean: {mean_disp:.2f} Å",
        )

        ax.set_ylabel(str(label_map[cond]), fontsize=9)
        ax.yaxis.grid(True, alpha=0.3)
        sns.despine(ax=ax, left=True)

    # Shared x-axis: readable integer ticks without changing numbering.
    last_ax = axes[-1, 0]
    last_ax.set_xlabel("Residue Index (auth_seq_id)", fontsize=10)
    last_ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    last_ax.set_xlim(x_min - 1, x_max + 1)
    last_ax.tick_params(axis="x", labelsize=9)

    # Single legend (identical entries across panels) on the first panel.
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        axes[0, 0].legend(
            handles, labels, loc="upper right", fontsize=8, framealpha=0.9,
        )

    axes[0, 0].set_title(
        f"  {title or 'Per-Residue Displacement'}",
        loc="left", fontsize=12, fontweight="bold",
    )

    if n_cond > 1:
        warnings.append(
            f"Panels show {n_cond} conditions; raw prediction-level points "
            "are drawn at reduced opacity with per-residue median/IQR "
            "summaries of the same observations"
        )

    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_path = save_path / "fig_v3_f04_per_residue_displacement.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    if n_obs == 0:
        warnings.append("Zero displacement observations")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
