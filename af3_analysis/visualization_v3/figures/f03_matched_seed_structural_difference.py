"""
V3 Figure F03 — Matched-Seed Structural Difference

Shows:
- RMSD per seed
- Coverage per seed
- Valid seed count
- Mean, median, SD, IQR
- Direction consistency where applicable

Extension of F02: adds seed-level detail and reproducibility.
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


def generate_f03_matched_seed_structural_difference(
    seed_rmsd_results: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
    seed_reproducibility: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """
    Generate matched-seed structural difference figure.

    Parameters
    ----------
    seed_rmsd_results : list of dicts with seed-level RMSD results
    save_path : Path to output directory
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title
    seed_reproducibility : optional seed-level reproducibility rows
        (same schema as the F17 data). When provided and usable, a third
        panel summarizing seed-level metric means is added; otherwise the
        figure keeps its original two-panel layout.

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not seed_rmsd_results:
        return {
            "status": "skip",
            "reason": "No seed-level RMSD results",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No matched-seed RMSD data available"],
        }

    # Build DataFrame
    seed_df = pd.DataFrame(seed_rmsd_results)

    if "seed" not in seed_df.columns or "rmsd_mean" not in seed_df.columns:
        return {
            "status": "skip",
            "reason": "Seed results missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Seed results incomplete"],
        }

    # Filter to valid results
    seed_df = seed_df[seed_df["rmsd_mean"].notna()].copy()

    if seed_df.empty:
        return {
            "status": "skip",
            "reason": "No valid seed-level RMSD",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All seed RMSD values are NaN"],
        }

    # Build labels
    if "condition_label" in seed_df.columns:
        label_col = "condition_label"
    elif "condition_a" in seed_df.columns:
        label_col = "condition_a"
    else:
        return {
            "status": "skip",
            "reason": "No condition column",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No condition identifier"],
        }

    if design and hasattr(design, "get_condition_label"):
        condition_map = {
            c: design.get_condition_label(c)
            for c in seed_df[label_col].unique()
        }
        seed_df["plot_label"] = seed_df[label_col].map(condition_map)
    else:
        seed_df["plot_label"] = seed_df[label_col]

    # Sort conditions
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in seed_df[label_col].unique()]
        order += [c for c in seed_df[label_col].unique() if c not in order]
    else:
        order = sorted(seed_df[label_col].unique())

    if design and hasattr(design, "get_condition_label"):
        label_map = {c: design.get_condition_label(c) for c in order}
    else:
        label_map = {c: c for c in order}

    categories = [label_map.get(c, c) for c in order]
    seed_df["plot_label"] = pd.Categorical(
        seed_df["plot_label"], categories=categories, ordered=True
    )

    n_obs = len(seed_df)

    # Repro data, if any, is used only to add a compact reproducibility
    # summary alongside the matched-seed RMSD view.
    repro_df = None
    if seed_reproducibility:
        repro_df = pd.DataFrame(seed_reproducibility)
        required_repro = {"metric_id", "condition_id", "mean", "n_seeds"}
        if not required_repro.issubset(repro_df.columns):
            repro_df = None
        else:
            repro_df = repro_df[repro_df["mean"].notna()].copy()
            if "condition_label" in repro_df.columns:
                repro_df["plot_label"] = repro_df["condition_label"]
            elif "condition_id" in repro_df.columns:
                repro_df["plot_label"] = repro_df["condition_id"].map(label_map).fillna(
                    repro_df["condition_id"]
                )
            else:
                repro_df = None

    # Choose layout: three panels when reproducibility data is usable,
    # otherwise keep the existing two-panel layout.
    if repro_df is not None and repro_df.shape[0] > 0:
        fig, axes = plt.subplots(
            1, 3,
            figsize=(min(11.0, SINGLE_COL_WIDTH * 1.7), 5.0),
            gridspec_kw={"width_ratios": [2, 1.25, 1]},
        )
    else:
        fig, axes = plt.subplots(
            1, 2, figsize=(SINGLE_COL_WIDTH, 5.0),
            gridspec_kw={"width_ratios": [2, 1]},
        )

    palette = sns.color_palette("Set2", n_colors=max(len(order), 3))
    color_map = {label_map.get(c, c): palette[i] for i, c in enumerate(order)}

    # Panel A: Seed-level RMSD with individual seeds
    ax = axes[0]

    # Box plot
    sns.boxplot(
        data=seed_df,
        x="plot_label",
        y="rmsd_mean",
        hue="plot_label",
        palette=color_map,
        width=0.5,
        fliersize=0,
        linewidth=0.8,
        boxprops=dict(alpha=0.4),
        legend=False,
        ax=ax,
    )

    # Show individual seeds
    sns.stripplot(
        data=seed_df,
        x="plot_label",
        y="rmsd_mean",
        hue="plot_label",
        palette=color_map,
        size=4,
        alpha=0.6,
        linewidth=0.3,
        edgecolor="white",
        jitter=0.2,
        legend=False,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Mean RMSD per seed (Å)")
    ax.set_title("  A. Seed-Level RMSD", loc="left", fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Coverage per seed
    ax = axes[1]

    if "coverage_mean" in seed_df.columns:
        df_plot = seed_df[seed_df["coverage_mean"].notna()].copy()
        sns.boxplot(
            data=df_plot,
            x="plot_label",
            y="coverage_mean",
            hue="plot_label",
            palette=color_map,
            width=0.5,
            fliersize=0,
            linewidth=0.8,
            boxprops=dict(alpha=0.4),
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=df_plot,
            x="plot_label",
            y="coverage_mean",
            hue="plot_label",
            palette=color_map,
            size=4,
            alpha=0.6,
            linewidth=0.3,
            edgecolor="white",
            jitter=0.2,
            legend=False,
            ax=ax,
        )

        ax.set_ylabel("Coverage")
        ax.set_title("  B. Coverage per Seed", loc="left", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1.05)
    else:
        ax.text(0.5, 0.5, "No coverage data", ha="center", va="center")
        ax.set_title("  B. Coverage per Seed", loc="left", fontsize=11, fontweight="bold")

    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel C: compact seed reproducibility summary when available.
    # This keeps seed-level robustness inside the same figure as the
    # matched-seed structural difference, instead of scattering it across
    # a separate figure.
    if repro_df is not None and repro_df.shape[0] > 0:
        ax = axes[2]

        # Summarize by metric and condition, averaging means within a
        # metric/condition cell so the panel stays readable.
        repro_summary = (
            repro_df
            .groupby(["metric_id", "plot_label"], observed=True)
            .agg(
                mean=("mean", "mean"),
                n=("mean", "size"),
            )
            .reset_index()
        )

        if repro_summary.shape[0] > 0:
            sns.barplot(
                data=repro_summary,
                x="metric_id",
                y="mean",
                hue="plot_label",
                palette=color_map,
                linewidth=0.8,
                edgecolor="white",
                ax=ax,
            )

            ax.set_xlabel("Metric")
            ax.set_ylabel("Mean across seeds")
            ax.set_title(
                "  C. Seed-Level Metric Means", loc="left", fontsize=10,
                fontweight="bold",
            )
            ax.tick_params(axis="x", rotation=30, labelsize=7)
            ax.tick_params(axis="y", labelsize=8)

            sns.despine(ax=ax, left=True)
            ax.yaxis.grid(True, alpha=0.3)

            handles, labels = ax.get_legend_handles_labels()
            if len(handles) > 6:
                ax.legend(
                    loc="center left", bbox_to_anchor=(1.0, 0.5),
                    fontsize=7,
                )
        else:
            ax.text(0.5, 0.5, "No reproducibility data", ha="center", va="center")
            ax.set_title(
                "  C. Seed-Level Metric Means", loc="left", fontsize=10,
                fontweight="bold",
            )

        # Update overall title to reflect the consolidated figure.
        main_title = title or "Matched-Seed Structural Difference & Reproducibility"
    else:
        main_title = title or "Matched-Seed Structural Difference"

    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f03_matched_seed_structural_difference.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero seed-level observations")
    if repro_df is not None and repro_df.shape[0] == 0:
        warnings.append(
            "Seed reproducibility data was provided but had no usable rows; "
            "Panel C omitted"
        )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
