"""
Cross-condition structural comparison orchestrator.

Manages which pairs of structures to compare, invokes alignment and
metrics, and produces normalised comparison results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from af3_analysis.structural.alignment import AlignmentResult, align_structures
from af3_analysis.structural.config import StructuralConfig
from af3_analysis.structural.enums import (
    AtomSelection,
    ComparisonMode,
    ComparisonReason,
    ComparisonStatus,
)
from af3_analysis.structural.metrics.registry import MetricRegistry, get_default_registry
from af3_analysis.structural.representation import NormalisedStructure


# ---------------------------------------------------------------------------
# Comparison result
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StructuralComparison:
    """One pairwise comparison between two predictions."""
    comparison_id: str
    prediction_a_id: str
    prediction_b_id: str
    condition_a: str
    condition_b: str
    seed_a: int
    seed_b: int

    alignment_status: ComparisonStatus
    alignment_reason: ComparisonReason
    n_common_atoms: int
    n_common_residues: int

    metrics: Dict[str, Optional[float]] = field(default_factory=dict)
    metric_status: Dict[str, str] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Prediction grouping
# ---------------------------------------------------------------------------

@dataclass
class PredictionGroup:
    """Group of predictions for a single condition."""
    condition_id: str
    predictions: Dict[int, List[NormalisedStructure]]  # seed -> [structures]


def group_predictions(
    structures: List[NormalisedStructure],
) -> Dict[str, PredictionGroup]:
    """Group structures by condition and seed."""
    groups: Dict[str, PredictionGroup] = {}

    for s in structures:
        if s.condition_id not in groups:
            groups[s.condition_id] = PredictionGroup(
                condition_id=s.condition_id,
                predictions={},
            )
        pg = groups[s.condition_id]
        if s.seed not in pg.predictions:
            pg.predictions[s.seed] = []
        pg.predictions[s.seed].append(s)

    return groups


# ---------------------------------------------------------------------------
# Comparison orchestrator
# ---------------------------------------------------------------------------

def _make_comparison_id(cond_a: str, seed_a: int, cond_b: str, seed_b: int) -> str:
    """Generate a deterministic comparison ID."""
    return f"{cond_a}_s{seed_a}_vs_{cond_b}_s{seed_b}"


def _make_prediction_id(condition_id: str, seed: int, sample: int) -> str:
    return f"{condition_id}_seed-{seed}_sample-{sample}"


def compare_conditions(
    structures: List[NormalisedStructure],
    config: StructuralConfig,
    *,
    registry: Optional[MetricRegistry] = None,
) -> List[StructuralComparison]:
    """
    Run structural comparisons across conditions.

    Parameters
    ----------
    structures : list of NormalisedStructure
        All parsed structures.
    config : StructuralConfig
        Configuration for the analysis.
    registry : MetricRegistry, optional
        Metric registry to use. If None, uses default.

    Returns
    -------
    list of StructuralComparison
    """
    if registry is None:
        registry = get_default_registry()

    groups = group_predictions(structures)
    atom_sel = AtomSelection(config.atom_selection)

    # Determine reference condition
    reference_condition = config.reference_condition
    if reference_condition is None and config.reference_strategy == "first_condition":
        reference_condition = sorted(groups.keys())[0] if groups else None

    # Determine which metrics to compute
    if config.enabled_metrics:
        metrics = registry.filter_by_ids(config.enabled_metrics)
    else:
        metrics = registry.all_metrics()

    # Generate comparison pairs
    pairs = _generate_pairs(groups, config, reference_condition)

    # Run comparisons
    comparisons = []
    for (cond_a, seed_a, sample_a, struct_a), (cond_b, seed_b, sample_b, struct_b) in pairs:
        # Skip self-comparisons
        if cond_a == cond_b and seed_a == seed_b and sample_a == sample_b:
            continue

        comparison = _compare_single(
            struct_a, struct_b, cond_a, seed_a, sample_a,
            cond_b, seed_b, sample_b, metrics, atom_sel, config,
        )
        comparisons.append(comparison)

    return comparisons


def _generate_pairs(
    groups: Dict[str, PredictionGroup],
    config: StructuralConfig,
    reference_condition: Optional[str],
) -> List[Tuple]:
    """
    Generate comparison pairs based on strategy.

    Returns list of ((cond_a, seed_a, sample_a, struct_a), (cond_b, seed_b, sample_b, struct_b)).
    """
    pairs = []

    if config.comparison_mode == "matched_seed" and reference_condition:
        # Compare each non-reference condition to reference at matched seeds
        ref_group = groups.get(reference_condition)
        if ref_group is None:
            return pairs

        for cond_id, group in groups.items():
            if cond_id == reference_condition:
                continue
            for seed in group.predictions:
                if seed not in ref_group.predictions:
                    continue
                for target_struct in group.predictions[seed]:
                    for ref_struct in ref_group.predictions[seed]:
                        pairs.append((
                            (cond_id, seed, target_struct.sample, target_struct),
                            (reference_condition, seed, ref_struct.sample, ref_struct),
                        ))

    elif config.comparison_mode == "pairwise_all":
        # All condition pairs, matched seeds
        cond_list = sorted(groups.keys())
        for i, cond_a in enumerate(cond_list):
            for cond_b in cond_list[i + 1:]:
                group_a = groups[cond_a]
                group_b = groups[cond_b]
                common_seeds = set(group_a.predictions.keys()) & set(group_b.predictions.keys())
                for seed in common_seeds:
                    for sa in group_a.predictions[seed]:
                        for sb in group_b.predictions[seed]:
                            pairs.append((
                                (cond_a, seed, sa.sample, sa),
                                (cond_b, seed, sb.sample, sb),
                            ))

    elif config.comparison_mode == "all_vs_all":
        # Every prediction vs every other prediction at same condition or different
        all_preds = []
        for cond_id, group in groups.items():
            for seed, structs in group.predictions.items():
                for s in structs:
                    all_preds.append((cond_id, seed, s.sample, s))
        for i, p_a in enumerate(all_preds):
            for p_b in all_preds[i + 1:]:
                pairs.append((p_a, p_b))

    return pairs


def _compare_single(
    struct_a: NormalisedStructure,
    struct_b: NormalisedStructure,
    cond_a: str,
    seed_a: int,
    sample_a: int,
    cond_b: str,
    seed_b: int,
    sample_b: int,
    metrics: list,
    atom_sel: AtomSelection,
    config: StructuralConfig,
) -> StructuralComparison:
    """Run a single pairwise comparison."""
    pred_a_id = _make_prediction_id(cond_a, seed_a, sample_a)
    pred_b_id = _make_prediction_id(cond_b, seed_b, sample_b)
    comp_id = _make_comparison_id(cond_a, seed_a, cond_b, seed_b)

    # Run alignment
    alignment = align_structures(
        struct_a, struct_b,
        atom_selection=atom_sel,
        min_common_atoms=config.min_common_atoms,
        min_sequence_identity=config.min_sequence_identity,
    )

    # Run metrics
    metric_values = {}
    metric_statuses = {}

    for metric in metrics:
        if alignment.status == ComparisonStatus.NOT_COMPARABLE:
            metric_values[metric.metric_id] = None
            metric_statuses[metric.metric_id] = "skipped_alignment_failed"
        elif not metric.is_applicable(struct_a, struct_b):
            metric_values[metric.metric_id] = None
            metric_statuses[metric.metric_id] = "not_applicable"
        else:
            try:
                value = metric.compute(
                    struct_a, struct_b, alignment,
                    contact_cutoff=config.contact_cutoff_angstrom,
                )
                metric_values[metric.metric_id] = value
                metric_statuses[metric.metric_id] = "present" if value is not None else "undefined"
            except Exception:
                metric_values[metric.metric_id] = None
                metric_statuses[metric.metric_id] = "computation_error"

    return StructuralComparison(
        comparison_id=comp_id,
        prediction_a_id=pred_a_id,
        prediction_b_id=pred_b_id,
        condition_a=cond_a,
        condition_b=cond_b,
        seed_a=seed_a,
        seed_b=seed_b,
        alignment_status=alignment.status,
        alignment_reason=alignment.reason,
        n_common_atoms=alignment.n_common_atoms,
        n_common_residues=alignment.n_common_residues,
        metrics=metric_values,
        metric_status=metric_statuses,
    )
