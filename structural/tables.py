"""
Output table schemas for structural analysis.

Produces normalised CSV/DataFrame tables from structural comparison results.
All output tables follow the same conventions as the existing confidence-metric
pipeline (tables/ directory, CSV format, provenance columns).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from af3_analysis.structural.comparison import StructuralComparison
from af3_analysis.structural.enums import ComparisonStatus, ParseStatus
from af3_analysis.structural.representation import NormalisedStructure


def build_structural_predictions_table(
    structures: List[NormalisedStructure],
    parse_errors: Optional[List[Dict]] = None,
) -> pd.DataFrame:
    """
    Build the structural_predictions.csv table.

    One row per successfully parsed CIF file.
    """
    rows = []
    for s in structures:
        rows.append({
            "prediction_id": f"{s.condition_id}_seed-{s.seed}_sample-{s.sample}",
            "condition_id": s.condition_id,
            "seed": s.seed,
            "sample": s.sample,
            "source_path": str(s.source_path),
            "source_checksum": s.source_checksum,
            "n_atoms": s.n_atoms,
            "n_chains": s.n_chains,
            "chain_ids": ";".join(s.chain_ids),
            "entity_types": ";".join(sorted(s.get_entity_types())),
            "parse_status": ParseStatus.SUCCESS.value,
            "parse_reason": "",
        })

    # Add failed parses
    if parse_errors:
        for err in parse_errors:
            rows.append({
                "prediction_id": f"{err['condition_id']}_seed-{err['seed']}_sample-{err['sample']}",
                "condition_id": err["condition_id"],
                "seed": err["seed"],
                "sample": err["sample"],
                "source_path": str(err.get("source_path", "")),
                "source_checksum": "",
                "n_atoms": 0,
                "n_chains": 0,
                "chain_ids": "",
                "entity_types": "",
                "parse_status": err.get("status", ParseStatus.PARSE_ERROR.value),
                "parse_reason": err.get("reason", ""),
            })

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "prediction_id", "condition_id", "seed", "sample",
        "source_path", "source_checksum", "n_atoms", "n_chains",
        "chain_ids", "entity_types", "parse_status", "parse_reason",
    ])


def build_structural_metrics_table(
    structures: List[NormalisedStructure],
    registry,
    config,
) -> pd.DataFrame:
    """
    Build the structural_metrics.csv table.

    One row per prediction × metric × scope.
    """
    from af3_analysis.structural.alignment import align_structures
    from af3_analysis.structural.enums import AtomSelection

    atom_sel = AtomSelection(config.atom_selection)
    rows = []

    for s in structures:
        # Self-alignment (identity alignment) for per-prediction metrics
        alignment = align_structures(
            s, s,
            atom_selection=atom_sel,
            min_common_atoms=0,
            min_sequence_identity=0.0,
        )

        metrics = registry.all_metrics()
        if config.enabled_metrics:
            metrics = registry.filter_by_ids(config.enabled_metrics)

        for metric in metrics:
            if not metric.is_applicable(s, s):
                rows.append({
                    "prediction_id": f"{s.condition_id}_seed-{s.seed}_sample-{s.sample}",
                    "condition_id": s.condition_id,
                    "seed": s.seed,
                    "metric_id": metric.metric_id,
                    "scope_type": metric.scope.value,
                    "scope_id": "",
                    "value": None,
                    "status": "not_applicable",
                    "reason": "metric not applicable to this structure",
                })
                continue

            if metric.scope.value == "chain":
                # Per-chain metrics
                for chain_id in s.chain_ids:
                    try:
                        value = metric.compute(s, s, alignment, chain_id=chain_id)
                        rows.append({
                            "prediction_id": f"{s.condition_id}_seed-{s.seed}_sample-{s.sample}",
                            "condition_id": s.condition_id,
                            "seed": s.seed,
                            "metric_id": metric.metric_id,
                            "scope_type": "chain",
                            "scope_id": chain_id,
                            "value": value,
                            "status": "present" if value is not None else "undefined",
                            "reason": "",
                        })
                    except Exception as e:
                        rows.append({
                            "prediction_id": f"{s.condition_id}_seed-{s.seed}_sample-{s.sample}",
                            "condition_id": s.condition_id,
                            "seed": s.seed,
                            "metric_id": metric.metric_id,
                            "scope_type": "chain",
                            "scope_id": chain_id,
                            "value": None,
                            "status": "computation_error",
                            "reason": str(e),
                        })
            else:
                # Global metrics
                try:
                    value = metric.compute(s, s, alignment)
                    rows.append({
                        "prediction_id": f"{s.condition_id}_seed-{s.seed}_sample-{s.sample}",
                        "condition_id": s.condition_id,
                        "seed": s.seed,
                        "metric_id": metric.metric_id,
                        "scope_type": metric.scope.value,
                        "scope_id": "",
                        "value": value,
                        "status": "present" if value is not None else "undefined",
                        "reason": "",
                    })
                except Exception as e:
                    rows.append({
                        "prediction_id": f"{s.condition_id}_seed-{s.seed}_sample-{s.sample}",
                        "condition_id": s.condition_id,
                        "seed": s.seed,
                        "metric_id": metric.metric_id,
                        "scope_type": metric.scope.value,
                        "scope_id": "",
                        "value": None,
                        "status": "computation_error",
                        "reason": str(e),
                    })

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "prediction_id", "condition_id", "seed", "metric_id",
        "scope_type", "scope_id", "value", "status", "reason",
    ])


def build_structural_comparisons_table(
    comparisons: List[StructuralComparison],
) -> pd.DataFrame:
    """
    Build the structural_comparisons.csv table.

    One row per comparison pair × metric.
    """
    rows = []
    for comp in comparisons:
        for metric_id, value in comp.metrics.items():
            status = comp.metric_status.get(metric_id, "unknown")
            rows.append({
                "comparison_id": comp.comparison_id,
                "prediction_a_id": comp.prediction_a_id,
                "prediction_b_id": comp.prediction_b_id,
                "condition_a": comp.condition_a,
                "condition_b": comp.condition_b,
                "seed_a": comp.seed_a,
                "seed_b": comp.seed_b,
                "metric_id": metric_id,
                "value": value,
                "alignment_status": comp.alignment_status.value,
                "alignment_reason": comp.alignment_reason.value,
                "n_common_atoms": comp.n_common_atoms,
                "n_common_residues": comp.n_common_residues,
                "metric_status": status,
            })

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "comparison_id", "prediction_a_id", "prediction_b_id",
        "condition_a", "condition_b", "seed_a", "seed_b",
        "metric_id", "value", "alignment_status", "alignment_reason",
        "n_common_atoms", "n_common_residues", "metric_status",
    ])


def build_structural_effects_table(
    comparisons: List[StructuralComparison],
) -> pd.DataFrame:
    """
    Build the structural_effects.csv table.

    One row per metric × comparison-condition-pair (seed-aggregated).
    """
    if not comparisons:
        return pd.DataFrame(columns=[
            "metric_id", "scope_type", "scope_id",
            "condition_a", "condition_b",
            "n", "mean", "median", "sd",
            "n_comparable", "n_incomparable",
        ])

    # Group by (condition_a, condition_b, metric_id)
    from collections import defaultdict
    import numpy as np

    groups = defaultdict(list)
    comparable_counts = defaultdict(int)
    incomparable_counts = defaultdict(int)

    for comp in comparisons:
        for metric_id, value in comp.metrics.items():
            key = (comp.condition_a, comp.condition_b, metric_id)
            if comp.alignment_status == ComparisonStatus.COMPARABLE:
                comparable_counts[key] += 1
            else:
                incomparable_counts[key] += 1
            if value is not None:
                groups[key].append(value)

    rows = []
    for (cond_a, cond_b, metric_id), values in groups.items():
        n = len(values)
        if n == 0:
            continue

        arr = np.array(values)
        rows.append({
            "metric_id": metric_id,
            "scope_type": "global",
            "scope_id": "",
            "condition_a": cond_a,
            "condition_b": cond_b,
            "n": n,
            "mean": float(np.mean(arr)),
            "median": float(np.median(arr)),
            "sd": float(np.std(arr, ddof=1)) if n > 1 else 0.0,
            "n_comparable": comparable_counts[(cond_a, cond_b, metric_id)],
            "n_incomparable": incomparable_counts[(cond_a, cond_b, metric_id)],
        })

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "metric_id", "scope_type", "scope_id",
        "condition_a", "condition_b",
        "n", "mean", "median", "sd",
        "n_comparable", "n_incomparable",
    ])


def build_structural_quality_table(
    structures: List[NormalisedStructure],
    parse_errors: Optional[List[Dict]] = None,
) -> pd.DataFrame:
    """
    Build the structural_quality.csv table.

    One row per condition × quality metric.
    """
    from collections import defaultdict
    import numpy as np

    # Group by condition
    by_condition = defaultdict(list)
    for s in structures:
        by_condition[s.condition_id].append(s)

    rows = []
    for cond_id, structs in sorted(by_condition.items()):
        atom_counts = [s.n_atoms for s in structs]
        chain_counts = [s.n_chains for s in structs]
        chains_per_pred = [len(s.chain_ids) for s in structs]

        rows.append({
            "condition_id": cond_id,
            "quality_metric": "n_predictions",
            "value": float(len(structs)),
            "n_predictions": len(structs),
            "n_failed": 0,
        })
        rows.append({
            "condition_id": cond_id,
            "quality_metric": "mean_atom_count",
            "value": float(np.mean(atom_counts)) if atom_counts else 0.0,
            "n_predictions": len(structs),
            "n_failed": 0,
        })
        rows.append({
            "condition_id": cond_id,
            "quality_metric": "n_entity_types",
            "value": float(len(set(
                e.entity_type for s in structs for e in s.entities
            ))),
            "n_predictions": len(structs),
            "n_failed": 0,
        })

        # Check consistency: do all predictions have the same chain set?
        chain_sets = [frozenset(s.chain_ids) for s in structs]
        unique_chain_sets = len(set(chain_sets))
        rows.append({
            "condition_id": cond_id,
            "quality_metric": "chain_set_consistency",
            "value": 1.0 if unique_chain_sets == 1 else 0.0,
            "n_predictions": len(structs),
            "n_failed": 0,
        })

    # Add failed parse counts
    if parse_errors:
        fail_counts = defaultdict(int)
        for err in parse_errors:
            fail_counts[err["condition_id"]] += 1
        for cond_id, count in fail_counts.items():
            rows.append({
                "condition_id": cond_id,
                "quality_metric": "n_failed_parses",
                "value": float(count),
                "n_predictions": 0,
                "n_failed": count,
            })

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        "condition_id", "quality_metric", "value", "n_predictions", "n_failed",
    ])


def write_structural_tables(
    output_dir: Path,
    structures: List[NormalisedStructure],
    comparisons: List[StructuralComparison],
    registry,
    config,
    *,
    parse_errors: Optional[List[Dict]] = None,
) -> Dict[str, Path]:
    """
    Write all structural output tables to the output directory.

    Parameters
    ----------
    output_dir : Path
        Tables directory (e.g., <run>/tables/).
    structures : list of NormalisedStructure
        All successfully parsed structures.
    comparisons : list of StructuralComparison
        All comparison results.
    registry : MetricRegistry
        Metric registry.
    config : StructuralConfig
        Analysis configuration.
    parse_errors : list of dict, optional
        Records of failed CIF parses.

    Returns
    -------
    dict mapping table name to file path.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {}

    # structural_predictions.csv
    df_pred = build_structural_predictions_table(structures, parse_errors)
    path = output_dir / "structural_predictions.csv"
    df_pred.to_csv(path, index=False)
    paths["structural_predictions"] = path

    # structural_metrics.csv
    df_met = build_structural_metrics_table(structures, registry, config)
    path = output_dir / "structural_metrics.csv"
    df_met.to_csv(path, index=False)
    paths["structural_metrics"] = path

    # structural_comparisons.csv
    df_comp = build_structural_comparisons_table(comparisons)
    path = output_dir / "structural_comparisons.csv"
    df_comp.to_csv(path, index=False)
    paths["structural_comparisons"] = path

    # structural_effects.csv
    df_eff = build_structural_effects_table(comparisons)
    path = output_dir / "structural_effects.csv"
    df_eff.to_csv(path, index=False)
    paths["structural_effects"] = path

    # structural_quality.csv
    df_qual = build_structural_quality_table(structures, parse_errors)
    path = output_dir / "structural_quality.csv"
    df_qual.to_csv(path, index=False)
    paths["structural_quality"] = path

    return paths
