"""
V3 Figure F13 family — Structural Similarity Views

The full prediction-level pairwise matrix (aligned Cα RMSD between each
pair of predictions, as produced by the existing structural modules) is
preserved as a computational artifact:

- ``v3/cache/pairwise_rmsd_matrix.csv`` with a ``*_metadata.csv`` sidecar
  (prediction_id, condition_id, seed, sample) written by the runner when a
  matrix is available;

the human-readable views are derived, without recomputation, by
aggregating that same matrix:

- **F13A** Condition × condition median pairwise RMSD (diagonal =
  within-condition pairwise RMSD across different predictions, not zero).
- **F13B** Within-condition prediction reproducibility distributions.
- **F13D** Prediction-level matrix, redesigned for readability (metadata
  ordering, no prediction-index labels, annotation strips, major-group
  boundaries only).

F13A and F13D are hidden from the default suite (prediction-level and
matrix presentations are not human-interpretable at scale); F13B is the
default F13 view. F13C (within vs between) was removed: its single-number
between-condition summary collapsed heterogeneous pairwise distributions;
use ``condition_similarity_summary.csv`` for the full condition-pair
detail instead.

Terminology: predictions grouped by condition are prediction-process
repeats, not biological replicates. Figures describe predicted structural
similarity only.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from af3_analysis.visualization_v3.config import (
    DPI,
    apply_v3_style,
)

METRIC_LABEL = "aligned Cα RMSD"


def _fmt(x: float) -> str:
    """Compact annotation formatting."""
    return f"{x:.2f}" if pd.notna(x) else "n/a"


def _safe_pairs(values: np.ndarray) -> np.ndarray:
    """Finite, non-self pairwise values from a masked selection."""
    vals = values[np.isfinite(values)]
    return vals[vals >= 0]


def aggregate_condition_matrix(
    matrix: np.ndarray,
    valid: Optional[np.ndarray],
    conditions: List[str],
    summary: str = "median",
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, Dict[str, float]]]:
    """Aggregate a prediction-level distance matrix to condition level.

    Uses only valid, off-diagonal entries of the existing matrix; the
    diagonal of the result is the within-condition distribution across
    DIFFERENT predictions (never auto-zero). Returns:

    - long-format summary table (one row per ordered condition pair plus
      diagonal rows), with median/mean/SD/IQR/n_pairs/coverage;
    - the (n_conditions, n_conditions) matrix of the summary statistic;
    - per-condition within-condition distribution stats.
    """
    conditions = list(conditions)
    cond_ids = sorted(dict.fromkeys(conditions))
    idx = {c: i for i, c in enumerate(cond_ids)}
    k = len(cond_ids)
    n = len(conditions)
    if matrix.shape[0] != n or matrix.shape[1] != n:
        raise ValueError("matrix shape does not match metadata length")

    if valid is None:
        valid = np.isfinite(matrix)

    long_rows: List[Dict[str, Any]] = []
    summary_matrix = np.full((k, k), np.nan)
    within: Dict[str, Dict[str, float]] = {}

    for ci, a in enumerate(cond_ids):
        rows_a = np.array([i for i in range(n) if conditions[i] == a])

        # Diagonal: within-condition pairs across different predictions.
        sub = matrix[np.ix_(rows_a, rows_a)]
        sub_valid = valid[np.ix_(rows_a, rows_a)]
        off_diag = ~np.eye(len(rows_a), dtype=bool)
        within_vals = sub[off_diag & sub_valid]
        total_within = len(rows_a) * (len(rows_a) - 1)  # ordered pairs
        within_summary = {
            "n_predictions": int(len(rows_a)),
            "n_pairs": int(within_vals.size),
            "median": float(np.median(within_vals)) if within_vals.size else np.nan,
            "mean": float(np.mean(within_vals)) if within_vals.size else np.nan,
            "sd": float(np.std(within_vals, ddof=1)) if within_vals.size > 1 else np.nan,
            "iqr": (float(np.percentile(within_vals, 75))
                    - float(np.percentile(within_vals, 25)))
            if within_vals.size else np.nan,
        }
        within[a] = within_summary
        summary_matrix[ci, ci] = within_summary[summary]
        long_rows.append({
            "condition_a": a, "condition_b": a,
            "median_rmsd": within_summary["median"],
            "mean_rmsd": within_summary["mean"],
            "sd_rmsd": within_summary["sd"],
            "iqr_rmsd": within_summary["iqr"],
            "n_pairs": within_summary["n_pairs"],
            "n_pairs_possible": total_within,
            "coverage": (within_summary["n_pairs"] / total_within
                         if total_within else np.nan),
            "comparison_type": "within",
        })

        # Off-diagonal: between-condition pairs (symmetric; computed once
        # per unordered pair and mirrored).
        for cj in range(ci + 1, k):
            b = cond_ids[cj]
            rows_b = np.array([i for i in range(n) if conditions[i] == b])
            block = matrix[np.ix_(rows_a, rows_b)]
            block_valid = valid[np.ix_(rows_a, rows_b)]
            vals = block[block_valid]
            n_pairs = int(vals.size)
            n_possible = int(rows_a.size * rows_b.size)
            if n_pairs:
                med = float(np.median(vals))
                mean = float(np.mean(vals))
                sd = float(np.std(vals, ddof=1)) if n_pairs > 1 else np.nan
                iqr = (float(np.percentile(vals, 75))
                       - float(np.percentile(vals, 25)))
            else:
                med = mean = sd = iqr = np.nan
            row = {
                "condition_a": a, "condition_b": b,
                "median_rmsd": med, "mean_rmsd": mean,
                "sd_rmsd": sd, "iqr_rmsd": iqr,
                "n_pairs": n_pairs, "n_pairs_possible": n_possible,
                "coverage": n_pairs / n_possible if n_possible else np.nan,
                "comparison_type": "between",
            }
            long_rows.append(row)
            mirror = dict(row)
            mirror["condition_a"], mirror["condition_b"] = b, a
            long_rows.append(mirror)
            summary_matrix[ci, cj] = summary_matrix[cj, ci] = med

    long_df = pd.DataFrame(long_rows)
    return long_df, summary_matrix, within


def within_condition_values(
    matrix: np.ndarray,
    valid: Optional[np.ndarray],
    conditions: List[str],
) -> Dict[str, np.ndarray]:
    """Raw within-condition pairwise values per condition (self excluded)."""
    conditions = list(conditions)
    n = len(conditions)
    if valid is None:
        valid = np.isfinite(matrix)
    out: Dict[str, np.ndarray] = {}
    for a in sorted(dict.fromkeys(conditions)):
        rows_a = np.array([i for i in range(n) if conditions[i] == a])
        sub = matrix[np.ix_(rows_a, rows_a)]
        sub_valid = valid[np.ix_(rows_a, rows_a)]
        off_diag = ~np.eye(len(rows_a), dtype=bool)
        vals = sub[off_diag & sub_valid]
        out[a] = vals[np.isfinite(vals)]
    return out


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _dynamic_figsize(n: int, cell: float = 0.42,
                     max_size: float = 14.0) -> float:
    return float(min(max_size, max(4.5, n * cell)))


def _wrap(text: str, width: int = 14) -> str:
    words = str(text).split()
    lines, cur = [], ""
    for w in words:
        cand = f"{cur} {w}".strip()
        if len(cand) > width and cur:
            lines.append(cur)
            cur = w
        else:
            cur = cand
    if cur:
        lines.append(cur)
    return "\n".join(lines)


def _discrete_palette(n: int) -> List[Tuple[float, float, float]]:
    return sns.color_palette("tab20", n_colors=max(n, 1))


# ---------------------------------------------------------------------------
# F13A — condition-level similarity
# ---------------------------------------------------------------------------

def generate_f13a_condition_similarity(
    summary_long: pd.DataFrame,
    summary_matrix: np.ndarray,
    condition_ids: List[str],
    save_path: Path,
    *,
    condition_labels: Optional[Dict[str, str]] = None,
    summary_stat: str = "median",
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """F13A: condition × condition median pairwise RMSD heatmap."""
    apply_v3_style()

    cond_ids = list(condition_ids)
    k = len(cond_ids)
    if k == 0 or summary_matrix.shape[0] != k:
        return {
            "status": "skip",
            "reason": "No condition-level matrix available",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Condition matrix unavailable or empty"],
        }

    labels = condition_labels or {c: c for c in cond_ids}
    display = [_wrap(labels.get(c, c)) for c in cond_ids]

    fig_size = _dynamic_figsize(k, cell=0.75, max_size=13.0)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    hm = sns.heatmap(
        summary_matrix,
        annot=(k <= 12),
        fmt=".2f" if k <= 12 else "",
        cmap="viridis",
        ax=ax,
        cbar_kws={"label": f"{summary_stat.capitalize()} pairwise {METRIC_LABEL} (Å)",
                  "fraction": 0.04, "pad": 0.02},
        xticklabels=display,
        yticklabels=display,
        linewidths=0.0,
        square=True,
    )

    # Subtle separators only when few conditions (never per-cell grids).
    if k <= 20:
        for i in range(1, k):
            ax.axhline(i, color="white", linewidth=0.8, alpha=0.7)
            ax.axvline(i, color="white", linewidth=0.8, alpha=0.7)

    ax.set_xlabel("Condition")
    ax.set_ylabel("Condition")
    ax.set_title(
        f"  {title or 'Condition-Level Structural Similarity'}",
        loc="left", fontsize=13, fontweight="bold",
    )
    ax.tick_params(axis="both", labelsize=8)

    # Documented color scale (no truncation): vmin=0 by construction of
    # RMSD; vmax extends to the data maximum.
    hm.collections[0].set_clim(0, float(np.nanmax(summary_matrix)))

    fig.tight_layout()

    out_path = save_path / "fig_v3_f13a_condition_similarity.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    n_pairs_total = int(summary_long["n_pairs"].sum()) if len(summary_long) else 0
    warnings = []
    if k > 12:
        warnings.append(
            f"{k} conditions: cell annotations omitted (readability)"
        )
    low_cov = summary_long[summary_long["coverage"] < 0.5] if len(summary_long) else []
    if len(low_cov):
        warnings.append(
            f"{len(low_cov)} condition pairs with coverage < 0.5; see "
            "condition_similarity_summary.csv"
        )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": k,
        "matrix_type": "condition x condition",
        "metric": METRIC_LABEL,
        "aggregation": f"{summary_stat} over valid pairwise values",
        "n_conditions": k,
        "warnings": warnings,
        "extra": {"n_pairs_total": n_pairs_total},
    }


# ---------------------------------------------------------------------------
# F13B — within-condition reproducibility
# ---------------------------------------------------------------------------

def generate_f13b_within_condition_reproducibility(
    within_values: Dict[str, np.ndarray],
    save_path: Path,
    *,
    condition_labels: Optional[Dict[str, str]] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """F13B: within-condition pairwise RMSD distributions (box + points)."""
    apply_v3_style()

    conds = [c for c in sorted(within_values) if len(within_values[c])]
    if not conds:
        return {
            "status": "skip",
            "reason": "No within-condition pairs (conditions need >= 2 "
                      "predictions with valid comparisons)",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No within-condition structural variability data"],
        }

    labels = condition_labels or {}
    data_rows = []
    for c in conds:
        for v in within_values[c]:
            data_rows.append({"condition_id": c, "value": float(v)})
    plot_df = pd.DataFrame(data_rows)
    plot_df["plot_label"] = plot_df["condition_id"].map(
        lambda c: _wrap(labels.get(c, c), 16))
    # Deterministic display order: by median of the existing values.
    order = sorted(conds, key=lambda c: float(np.median(within_values[c])))

    fig_width = max(6.0, min(13.0, 0.9 * len(conds) + 2.5))
    fig, ax = plt.subplots(figsize=(fig_width, 5.2))

    sns.boxplot(
        data=plot_df, x="plot_label", y="value",
        order=[_wrap(labels.get(c, c), 16) for c in order],
        color="#BDBDBD", width=0.55, linewidth=0.8,
        showfliers=False, ax=ax,
    )
    sns.stripplot(
        data=plot_df, x="plot_label", y="value",
        order=[_wrap(labels.get(c, c), 16) for c in order],
        color="#2C7BB6", size=2.5, alpha=0.25, jitter=0.25, ax=ax,
    )

    ax.set_xlabel("Condition")
    ax.set_ylabel(f"Within-condition pairwise {METRIC_LABEL} (Å)")
    ax.set_title(
        f"  {title or 'Within-Condition Prediction Reproducibility'}",
        loc="left", fontsize=13, fontweight="bold",
    )
    ax.tick_params(axis="x", rotation=0, labelsize=8)
    ax.tick_params(axis="y", labelsize=9)
    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    fig.tight_layout()
    out_path = save_path / "fig_v3_f13b_within_condition_reproducibility.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    n_obs = int(sum(len(v) for v in within_values.values()))
    warnings = [
        "Within-condition values compare DIFFERENT predictions from the "
        "same condition (prediction reproducibility); predictions are not "
        "biological replicates"
    ]
    singletons = [c for c in sorted(within_values) if len(within_values[c]) == 0]
    if singletons:
        warnings.append(
            f"{len(singletons)} conditions with fewer than 2 predictions "
            "omitted from distributions"
        )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "matrix_type": "within-condition distributions",
        "metric": METRIC_LABEL,
        "aggregation": "raw within-condition pairwise values (self excluded)",
        "n_conditions": len(conds),
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# F13D — redesigned prediction-level matrix
# ---------------------------------------------------------------------------

def generate_f13d_prediction_matrix(
    distance_matrix: np.ndarray,
    valid: Optional[np.ndarray],
    predictions: List[str],
    conditions: List[str],
    seeds: List[int],
    samples: List[int],
    save_path: Path,
    *,
    cluster_labels: Optional[np.ndarray] = None,
    annotation_strips: Optional[Dict[str, List[Any]]] = None,
    ordering: str = "cluster_condition_seed_sample",
    show_dendrogram: bool = False,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """F13D: prediction-level matrix with metadata ordering and strips."""
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
    if distance_matrix.shape[0] != n:
        return {
            "status": "skip",
            "reason": "Distance matrix shape mismatch",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Matrix dimensions inconsistent with metadata"],
        }

    # --- Documented, deterministic ordering -----------------------------
    # Default: existing predicted structural cluster, then condition, then
    # seed, then sample — all from existing metadata; no re-clustering.
    if ordering == "hierarchical" and cluster_labels is None:
        ordering = "cluster_condition_seed_sample"

    if cluster_labels is not None and len(cluster_labels) == n:
        order = sorted(
            range(n),
            key=lambda i: (int(cluster_labels[i]), str(conditions[i]),
                           int(seeds[i]), int(samples[i]),
                           str(predictions[i])),
        )
    else:
        order = sorted(
            range(n),
            key=lambda i: (str(conditions[i]), int(seeds[i]),
                           int(samples[i]), str(predictions[i])),
        )

    ordered_matrix = distance_matrix[np.ix_(order, order)]
    ordered_conditions = [conditions[i] for i in order]
    ordered_valid = valid[np.ix_(order, order)] if valid is not None else None

    # Major group boundaries only (between conditions / clusters).
    boundaries = [i for i in range(1, n) if ordered_conditions[i] != ordered_conditions[i - 1]]

    # --- Annotation strips (generic; from provided metadata) -----------
    strip_specs: List[Tuple[str, List[Any], Dict[Any, Any]]] = []
    if annotation_strips:
        cond_colors = _discrete_palette(len(set(ordered_conditions)))
        cond_map = {c: cond_colors[i]
                    for i, c in enumerate(sorted(set(ordered_conditions)))}
        strip_specs.append(("Condition",
                            [cond_map[c] for c in ordered_conditions], {}))
        for name, values in annotation_strips.items():
            ordered_vals = [values[i] for i in order]
            uniq = sorted(set(ordered_vals), key=str)
            palette = _discrete_palette(max(len(uniq), 1))
            vmap = {v: palette[i % len(palette)] for i, v in enumerate(uniq)}
            strip_specs.append((name, [vmap[v] for v in ordered_vals], {}))

    n_strips = len(strip_specs)
    strip_h = 0.6
    mat_size = _dynamic_figsize(n, cell=0.055, max_size=13.0)
    fig_h = mat_size + (n_strips * strip_h if n < 200 else n_strips * strip_h * 0.8)
    fig, axes = plt.subplots(
        n_strips + 1, 1,
        figsize=(mat_size, fig_h),
        gridspec_kw={"height_ratios": [strip_h] * n_strips + [mat_size],
                     "hspace": 0.04},
        squeeze=False,
    )

    # Strips: (1, N, 3) RGB arrays drawn with pcolormesh.
    for r, (name, colors, _) in enumerate(strip_specs):
        ax = axes[r, 0]
        strip_rgb = np.array(
            [c[:3] for c in colors], dtype=float,
        ).reshape(1, -1, 3)
        ax.pcolormesh(
            np.arange(n + 1), np.arange(2), strip_rgb,
            rasterized=True,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=7)
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Heatmap panel (rasterized, no per-prediction ticks, no dense grid).
    ax = axes[n_strips, 0]
    ax.imshow(ordered_matrix, cmap="viridis", vmin=0,
              aspect="auto" if n > 120 else "equal",
              interpolation="nearest", rasterized=True)
    for i in boundaries:
        ax.axhline(i, color="white", linewidth=0.6, alpha=0.5)
        ax.axvline(i, color="white", linewidth=0.6, alpha=0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Predictions ordered by structural similarity "
                  "(cluster > condition > seed > sample)", fontsize=9)
    ax.set_ylabel("Predicted structures", fontsize=9)

    # Colorbar in documented units.
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(
            norm=plt.Normalize(vmin=0, vmax=float(np.nanmax(ordered_matrix))),
            cmap="viridis",
        ),
        ax=ax, fraction=0.04, pad=0.02,
    )
    cbar.set_label(f"Pairwise {METRIC_LABEL} (Å)", fontsize=9)

    fig.suptitle(
        title or "Prediction-Level Structural Similarity",
        fontsize=13, fontweight="bold", y=1.0,
    )

    out_path = save_path / "fig_v3_f13d_prediction_matrix.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = [
        "Ordering: " + (
            "existing predicted structural clusters, then condition, seed, "
            "sample (no re-clustering)" if cluster_labels is not None
            else "condition, seed, sample from metadata (no clustering)"
        ),
        "Each row/column is one prediction; tick labels are omitted at this "
        "scale — prediction metadata is in prediction_metadata.csv and the "
        "matrix CSV",
    ]
    if boundaries:
        warnings.append(
            f"{len(boundaries)} condition-group boundaries drawn (major "
            "groups only; no per-prediction grid)"
        )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n,
        "matrix_type": "prediction x prediction",
        "metric": METRIC_LABEL,
        "ordering": ordering,
        "clustering_method": ("existing hierarchical clustering labels"
                              if cluster_labels is not None else "none"),
        "n_boundaries": len(boundaries),
        "warnings": warnings,
    }
