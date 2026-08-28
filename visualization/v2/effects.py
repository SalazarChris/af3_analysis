"""
V2 Effects — Factor Effects Forest Plot
========================================

Horizontal effect-size plot showing main effects and interactions.
Only plotted when pairwise_comparisons.csv contains non-empty data.

Encoding
--------
* Dots with confidence intervals.
* Ordered by metric group.
* Red dashed line at zero (no effect).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from .config import (
    DPI, DOUBLE_COL_WIDTH, FONT_SIZES, MAX_FIGURE_HEIGHT_INCHES,
    apply_v2_style,
)
from .labels import metric_label, figure_title


def generate_effects_forest(
    pairwise_comparisons: pd.DataFrame,
    save_path: Path,
    *,
    max_comparisons_per_metric: int = 20,
) -> Dict[str, Any]:
    """Generate the effect-size forest plot.

    Parameters
    ----------
    pairwise_comparisons : DataFrame
        Must contain ``metric``, ``condition``, ``reference``,
        ``hedges_g`` or ``diff_mean``.
    save_path : Path
        Output directory.
    max_comparisons_per_metric : int
        Cap per-metric row count.

    Returns
    -------
    dict with ``status``, ``output_path``, ``warnings``.
    """
    apply_v2_style()

    if pairwise_comparisons is None or pairwise_comparisons.empty:
        return {"status": "skip", "reason": "No pairwise comparisons available"}

    df = pairwise_comparisons.copy()

    val_col = "hedges_g" if "hedges_g" in df.columns else "diff_mean"
    if val_col not in df.columns:
        return {"status": "skip", "reason": f"No value column ({val_col}) found"}

    metrics = df["metric"].unique().tolist()
    n_metrics = len(metrics)

    # Filter to primary metrics only (exclude chain-level)
    primary = [m for m in metrics if not m.startswith("chain_")]
    if not primary:
        return {"status": "skip", "reason": "No primary metrics in pairwise data"}

    # Calculate figure height
    total_rows = 0
    for m in primary:
        subset = df[df["metric"] == m]
        total_rows += min(len(subset), max_comparisons_per_metric)
    fig_height = max(3.0, min(total_rows * 0.45, MAX_FIGURE_HEIGHT_INCHES))

    fig, axes = plt.subplots(n_metrics, 1, figsize=(DOUBLE_COL_WIDTH, fig_height))
    if n_metrics == 1:
        axes = [axes]

    warnings_list: List[str] = []

    for i, (ax, metric) in enumerate(zip(axes, primary)):
        subset = df[df["metric"] == metric].copy()
        subset = subset.sort_values(val_col).reset_index(drop=True)

        if len(subset) > max_comparisons_per_metric:
            subset = subset.head(max_comparisons_per_metric)
            warnings_list.append(
                f"WARNING: {metric} has {len(df[df['metric'] == metric])} comparisons; "
                f"showing top {max_comparisons_per_metric}"
            )

        # Build contrast label
        subset["contrast"] = (
            subset["condition"].astype(str) + "  vs  " + subset["reference"].astype(str)
        )

        # Error bars
        ci_lower_col = f"{val_col}_ci_lower"
        ci_upper_col = f"{val_col}_ci_upper"
        has_ci = ci_lower_col in subset.columns and ci_upper_col in subset.columns

        if has_ci:
            y_err = [
                subset[val_col] - subset[ci_lower_col],
                subset[ci_upper_col] - subset[val_col],
            ]
            ax.errorbar(
                subset[val_col], range(len(subset)),
                xerr=y_err, fmt="o", color="#2C3E50",
                capsize=3, elinewidth=1.2, markersize=5,
            )
        else:
            ax.plot(subset[val_col], range(len(subset)), "o",
                    color="#2C3E50", markersize=5)

        ax.set_yticks(range(len(subset)))
        ax.set_yticklabels(subset["contrast"], fontsize=FONT_SIZES["tick_label"])
        ax.axvline(0, color="#C0392B", linestyle="--", alpha=0.6, linewidth=1.0)

        label = metric_label(metric)
        ax.set_title(
            f"  {chr(65+i)}. {label}",
            loc="left", fontsize=FONT_SIZES["panel_title"], fontweight="bold",
        )
        ax.set_xlabel(
            f"Effect size ({val_col.replace('_', ' ').title()})",
            fontsize=FONT_SIZES["axis_label"],
        )
        ax.grid(axis="x", alpha=0.3)

    fig.suptitle(
        figure_title("effects"),
        fontsize=FONT_SIZES["figure_title"],
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out = save_path / "figure_effects_forest_v2.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "pass",
        "output_path": str(out),
        "metrics": primary,
        "warnings": warnings_list,
    }
