"""
V3 Figure F01 — Structural Data Completeness / QC

Shows:
- Number of structures per condition
- Parse success/failure rates
- Entity composition (protein, DNA, ligand, ion)
- Chain composition
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


def generate_f01_structural_qc(
    qc_records: List[Dict[str, Any]],
    save_path: Path,
    *,
    design: Optional[Any] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate structural QC completeness figure.

    Parameters
    ----------
    qc_records : list of dicts with structural QC info
    save_path : Path to output directory
    design : optional experiment design metadata
    title : optional custom title

    Returns
    -------
    dict with status, output_path, n_observations, warnings
    """
    apply_v3_style()

    if not qc_records:
        return {
            "status": "skip",
            "reason": "No QC records",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["No structural QC data available"],
        }

    # Build DataFrame
    df = pd.DataFrame(qc_records)

    if "condition_id" not in df.columns or "passed" not in df.columns:
        return {
            "status": "skip",
            "reason": "QC records missing required columns",
            "output_path": None,
            "n_observations": 0,
            "warnings": ["QC records incomplete"],
        }

    # Count per condition
    condition_counts = df.groupby("condition_id").agg({
        "passed": ["sum", "count"],
    })
    condition_counts.columns = ["n_passed", "n_total"]
    condition_counts["n_failed"] = condition_counts["n_total"] - condition_counts["n_passed"]
    condition_counts["pass_rate"] = condition_counts["n_passed"] / condition_counts["n_total"]

    # Sort by condition name
    if design and hasattr(design, "condition_names"):
        order = [c for c in design.condition_names if c in condition_counts.index]
        order += [c for c in condition_counts.index if c not in order]
    else:
        order = sorted(condition_counts.index)

    condition_counts = condition_counts.loc[order]

    # Build labels
    if design and hasattr(design, "get_condition_label"):
        labels = [design.get_condition_label(c) for c in condition_counts.index]
    else:
        labels = list(condition_counts.index)

    n_obs = int(condition_counts["n_total"].sum())

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(SINGLE_COL_WIDTH, 5.0))

    # Panel A: Pass/Fail bar chart
    ax = axes[0]
    x_pos = np.arange(len(condition_counts))
    width = 0.6

    ax.bar(
        x_pos,
        condition_counts["n_passed"],
        width,
        label="Passed",
        color="#2C7BB6",
        edgecolor="white",
    )
    ax.bar(
        x_pos,
        condition_counts["n_failed"],
        width,
        bottom=condition_counts["n_passed"],
        label="Failed",
        color="#D7191C",
        edgecolor="white",
    )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Number of structures")
    ax.set_title("  A. Parse Success/Failure", loc="left", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right")
    ax.tick_params(axis="y", labelsize=9)

    sns.despine(ax=ax, left=True)
    ax.yaxis.grid(True, alpha=0.3)

    # Panel B: Entity composition
    ax = axes[1]

    # Count entity types across all structures
    entity_counts = {}
    for qc in qc_records:
        if "entity_types" in qc and qc["entity_types"]:
            for etype in qc["entity_types"]:
                entity_counts[etype] = entity_counts.get(etype, 0) + 1

    if entity_counts:
        entity_items = sorted(entity_counts.items(), key=lambda x: -x[1])
        entity_names = [item[0] for item in entity_items]
        entity_values = [item[1] for item in entity_items]

        palette = sns.color_palette("Set2", len(entity_names))
        wedges, texts, autotexts = ax.pie(
            entity_values,
            labels=entity_names,
            autopct="%1.1f%%",
            colors=palette,
            startangle=90,
        )
        for t in autotexts:
            t.set_fontsize(8)
        for t in texts:
            t.set_fontsize(9)
    else:
        ax.text(0.5, 0.5, "No entity data", ha="center", va="center")
        ax.set_title("  B. Entity Composition", loc="left", fontsize=11, fontweight="bold")

    ax.set_title("  B. Entity Composition", loc="left", fontsize=11, fontweight="bold")

    # Main title
    main_title = title or "Structural Data Quality Control"
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = save_path / "fig_v3_f01_structural_qc.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

    warnings = []
    if n_obs == 0:
        warnings.append("Zero structures in dataset")

    return {
        "status": "pass",
        "output_path": str(out_path),
        "n_observations": n_obs,
        "warnings": warnings,
    }
