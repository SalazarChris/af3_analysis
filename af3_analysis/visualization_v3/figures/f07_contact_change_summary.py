"""
V3 Figure F07 — Contact Change Summary

For each condition vs reference:
- contacts gained
- contacts lost
- percentage changed
- seed consistency
- recurring residue-pair changes

Presentation (visualization only; no analysis change):
- when the data contains multiple seed-level rows per condition, panels A/B
  show one row per condition summarizing the EXISTING per-seed values
  (mean, with the individual seed values overlaid); the seed-level rows
  themselves remain unchanged in the contact_changes output table;
- horizontal layouts with horizontal labels (no rotated tick text);
- Panel C shows the ranked residue-pair changes consolidated from the
  former F06 figure.
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
    DOUBLE_COL_WIDTH,
    SINGLE_COL_WIDTH,
    apply_v3_style,
)

# Ranked-panel display cap. The full changed-pair set is preserved in the
# runner aggregation; only the drawn rows are capped, and the figure
# reports the truncation in its warnings.
MAX_RANKED_PAIRS = 20


def _wrap_label(text: str, width: int = 18) -> str:
    """Wrap a display label on spaces (display only; IDs unchanged)."""
    words = str(text).split()
    if len(words) <= 1:
        return str(text)
    lines, current = [], ""
    for w in words:
        candidate = f"{current} {w}".strip()
        if len(candidate) > width and current:
            lines.append(current)
            current = w
        else:
            current = candidate
    if current:
        lines.append(current)
    return "\n".join(lines)


def generate_f07_contact_change_summary(
    contact_change_data: List[Dict[str, Any]],
    save_path: Path,
    *,
    recurring_pairs: Optional[List[Dict[str, Any]]] = None,
    design: Optional[Any] = None,
    reference_condition: Optional[str] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate contact change summary figure.

    Parameters
    ----------
    contact_change_data : list of contact change summary dicts
    save_path : Path to output directory
    recurring_pairs : optional list of per-pair change dicts with keys
        ``chain_a, seq_a, chain_b, seq_b, n_gained, n_lost`` (and optional
        ``n_conditions``). Aggregated changed-pair detail consolidated
        from the former F06 figure. When supplied and usable, a third
        panel shows the ranked pair changes; otherwise the figure keeps
        its two-panel layout.
    design : optional experiment design metadata
    reference_condition : str, optional
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not contact_change_data:
        return {
            "status": "skip",
            "reason": "No contact change data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No contact change data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(contact_change_data)

    required_cols = ["condition_id", "n_gained", "n_lost", "n_total"]
    if not all(col in df.columns for col in required_cols):
        return {
            "status": "skip",
            "reason": "Contact change data missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["Contact change data incomplete"],
        }

    df = df[df["n_total"] > 0].copy()

    if df.empty:
        return {
            "status": "skip",
            "reason": "No valid contact data",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["All contact totals are zero"],
        }

    # Existing per-row fraction changed (unchanged definition)
    df["pct_changed"] = (df["n_gained"] + df["n_lost"]) / df["n_total"]

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        df["plot_label"] = df["condition_id"].apply(
            lambda c: design.get_condition_label(c) if c in design.conditions else c
        )
    else:
        df["plot_label"] = df["condition_id"]

    # Condition order from existing metadata, ranked by the existing
    # per-condition mean fraction changed for display.
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in df["condition_id"].unique()]
        order += [c for c in df["condition_id"].unique() if c not in order]
    else:
        order = sorted(df["condition_id"].unique())

    label_map = {c: (design.get_condition_label(c) if design and hasattr(design, "get_condition_label") else c)
                 for c in order}

    n_obs = len(df)

    # Display-level aggregation: mean of the existing per-seed rows per
    # condition (descriptive summary; seed rows are NOT removed from the
    # data — they remain in the contact_changes output table and are
    # overlaid as points in panel B).
    cond_stats = (
        df.groupby("condition_id")
        .agg(
            mean_gained=("n_gained", "mean"),
            mean_lost=("n_lost", "mean"),
            mean_pct=("pct_changed", "mean"),
        )
        .reindex(order)
        .dropna(how="all")
        .reset_index()
    )
    cond_stats["plot_label"] = cond_stats["condition_id"].map(label_map)
    # Rank by the existing fraction-changed summary (display order only).
    cond_stats = cond_stats.sort_values("mean_pct", ascending=True)

    y_labels = [_wrap_label(l) for l in cond_stats["plot_label"]]

    # Panel C data: per-pair contact changes consolidated from F06.
    pairs_df = None
    if recurring_pairs:
        required_pair_cols = {
            "chain_a", "seq_a", "chain_b", "seq_b", "n_gained", "n_lost"
        }
        pairs_df = pd.DataFrame(recurring_pairs)
        if required_pair_cols.issubset(pairs_df.columns):
            pairs_df = pairs_df[
                (pairs_df["n_gained"] > 0) | (pairs_df["n_lost"] > 0)
            ].copy()
            # Existing ordering: by the plotted change count.
            pairs_df["total_changes"] = (
                pairs_df["n_gained"] + pairs_df["n_lost"]
            )
            pairs_df = pairs_df.sort_values("total_changes", ascending=True)
            if pairs_df.empty:
                pairs_df = None
        else:
            pairs_df = None

    # --- Plot ---
    if pairs_df is not None:
        fig, axes = plt.subplots(
            1, 3,
            figsize=(DOUBLE_COL_WIDTH, 0.55 * len(cond_stats) + 2.2),
            gridspec_kw={"width_ratios": [1.1, 1, 1.25]},
        )
    else:
        fig, axes = plt.subplots(
            1, 2,
            figsize=(SINGLE_COL_WIDTH * 1.4, 0.55 * len(cond_stats) + 2.2),
        )

    y_pos = np.arange(len(cond_stats))

    # Panel A: mean gained/lost contacts per condition (existing values,
    # averaged for display across the existing seed rows).
    ax = axes[0]
    ax.barh(
        y_pos + 0.2, cond_stats["mean_gained"], height=0.4,
        label="Gained", color="#2C7BB6", edgecolor="white",
    )
    ax.barh(
        y_pos - 0.2, cond_stats["mean_lost"], height=0.4,
        label="Lost", color="#D7191C", edgecolor="white",
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel("Contacts changed (mean across seeds)", fontsize=9)
    ax.set_title("  A. Contact Changes per Condition",
                 loc="left", fontsize=10, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    sns.despine(ax=ax, left=True)
    ax.xaxis.grid(True, alpha=0.3)

    # Panel B: fraction changed (existing per-row definition). Condition
    # means with the existing seed-level values overlaid as points.
    ax = axes[1]
    ax.barh(
        y_pos, cond_stats["mean_pct"], height=0.55,
        color="#66C2A5", edgecolor="white", alpha=0.85,
    )
    seed_y = {cid: i for i, cid in enumerate(cond_stats["condition_id"])}
    ax.scatter(
        df["pct_changed"],
        df["condition_id"].map(seed_y),
        s=12, alpha=0.55, color="#2C7BB6",
        edgecolors="white", linewidth=0.3, zorder=3,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel("Fraction of contacts changed", fontsize=9)
    ax.set_title("  B. Contact Change Fraction",
                 loc="left", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", labelsize=8)
    sns.despine(ax=ax, left=True)
    ax.xaxis.grid(True, alpha=0.3)

    # Panel C: ranked recurring residue-pair changes (consolidated F06).
    warnings = []
    if pairs_df is not None:
        ax = axes[2]
        top = pairs_df.tail(MAX_RANKED_PAIRS)
        pair_labels = [
            f"{row['chain_a']}:{row['seq_a']}–{row['chain_b']}:{row['seq_b']}"
            for _, row in top.iterrows()
        ]
        y_pos_c = np.arange(len(top))
        ax.barh(
            y_pos_c + 0.2, top["n_gained"], height=0.4,
            label="Gained", color="#2C7BB6", edgecolor="white",
        )
        ax.barh(
            y_pos_c - 0.2, top["n_lost"], height=0.4,
            label="Lost", color="#D7191C", edgecolor="white",
        )
        ax.set_yticks(y_pos_c)
        ax.set_yticklabels(pair_labels, fontsize=7)
        ax.set_xlabel("Comparisons with change", fontsize=9)
        ax.set_title("  C. Recurring Contact Changes",
                     loc="left", fontsize=10, fontweight="bold")
        ax.legend(loc="lower right", fontsize=8)
        ax.tick_params(axis="x", labelsize=8)
        sns.despine(ax=ax, left=True)
        ax.xaxis.grid(True, alpha=0.3)
        if len(pairs_df) > MAX_RANKED_PAIRS:
            warnings.append(
                f"{len(pairs_df)} changed residue pairs in total; Panel C "
                f"shows the {MAX_RANKED_PAIRS} most frequently changed "
                "(ranked display of existing values)"
            )

    # Main title
    main_title = title or "Contact Change Summary"
    fig.suptitle(main_title, fontsize=12, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f07_contact_change_summary.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    if n_obs == 0:
        warnings.append("Zero contact change observations")
    n_conds = df["condition_id"].nunique()
    if n_obs > n_conds:
        warnings.append(
            f"{n_obs} seed-level rows across {n_conds} conditions; panels "
            "show per-condition means of the same rows (seed values "
            "overlaid in Panel B; full per-seed data in contact_changes.csv)"
        )

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
