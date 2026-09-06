"""
V3 Reference Resolution.

Explicitly resolves the reference structure/condition from configuration.
Supports:
- condition → reference
- condition A → condition B

For matched-seed comparisons, uses:
reference(seed N) ↔ condition(seed N)
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .model import Dataset, StructureData


class ReferenceResolutionError(Exception):
    """Error resolving reference."""
    pass


def resolve_reference(
    dataset: Dataset,
    config_reference: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Resolve reference from configuration.

    Parameters
    ----------
    dataset : Dataset
        The V3 dataset.
    config_reference : dict, optional
        Reference configuration from V3Config.reference.

    Returns
    -------
    dict with:
        - reference_condition: str
        - reference_strategy: str
        - paired_conditions: list of (ref_condition, target_condition)
    """
    result = {
        "reference_condition": None,
        "reference_strategy": "none",
        "paired_conditions": [],
    }

    # An empty reference config is equivalent to no explicit reference
    # (V3Config.reference defaults to {}, not None).
    if config_reference is None or not config_reference:
        # Default: first condition alphabetically
        if dataset.conditions:
            ref_cond = sorted(dataset.conditions.keys())[0]
            result["reference_condition"] = ref_cond
            result["reference_strategy"] = "first_condition"
            result["paired_conditions"] = [
                (ref_cond, cond)
                for cond in dataset.conditions.keys()
                if cond != ref_cond
            ]
        return result

    # Explicit reference configuration
    if "condition" in config_reference:
        ref_cond = config_reference["condition"]
        if ref_cond not in dataset.conditions:
            raise ReferenceResolutionError(
                f"Reference condition '{ref_cond}' not in dataset"
            )
        result["reference_condition"] = ref_cond
        result["reference_strategy"] = "explicit_reference"
        result["paired_conditions"] = [
            (ref_cond, cond)
            for cond in dataset.conditions.keys()
            if cond != ref_cond
        ]

    elif "condition_a" in config_reference and "condition_b" in config_reference:
        # Pairwise reference
        cond_a = config_reference["condition_a"]
        cond_b = config_reference["condition_b"]
        if cond_a not in dataset.conditions:
            raise ReferenceResolutionError(
                f"Condition A '{cond_a}' not in dataset"
            )
        if cond_b not in dataset.conditions:
            raise ReferenceResolutionError(
                f"Condition B '{cond_b}' not in dataset"
            )
        result["reference_condition"] = cond_a
        result["reference_strategy"] = "pairwise_reference"
        result["paired_conditions"] = [(cond_a, cond_b)]

    elif "conditions" in config_reference:
        # Multiple conditions as reference
        ref_conds = config_reference["conditions"]
        invalid = [c for c in ref_conds if c not in dataset.conditions]
        if invalid:
            raise ReferenceResolutionError(
                f"Reference conditions not in dataset: {invalid}"
            )
        result["reference_condition"] = ref_conds[0] if ref_conds else None
        result["reference_strategy"] = "multi_reference"
        result["paired_conditions"] = [
            (ref, cond)
            for ref in ref_conds
            for cond in dataset.conditions.keys()
            if cond not in ref_conds
        ]

    return result


def get_reference_structure(
    dataset: Dataset,
    condition_id: str,
    seed: int,
    sample: Optional[int] = None,
) -> Optional[StructureData]:
    """
    Get the reference structure for a given condition/seed.

    For matched-seed comparisons, uses reference(seed N) ↔ condition(seed N).

    Parameters
    ----------
    dataset : Dataset
        The V3 dataset.
    condition_id : str
        The condition to get reference for.
    seed : int
        The seed.
    sample : int, optional
        Specific sample, or None for first available.

    Returns
    -------
    StructureData or None
    """
    ref_cond = dataset.reference_condition
    if ref_cond is None:
        return None

    # Get reference condition's structures
    ref_seeds = dataset.predictions.get(ref_cond, {})
    ref_seed_data = ref_seeds.get(seed, {})

    if not ref_seed_data:
        return None

    if sample is not None:
        return ref_seed_data.get(sample)
    else:
        # Return first available sample
        return next(iter(ref_seed_data.values())) if ref_seed_data else None


def get_matched_seed_pairs(
    dataset: Dataset,
    target_condition: str,
    reference_condition: Optional[str] = None,
) -> List[Tuple[int, StructureData, StructureData]]:
    """
    Get matched-seed pairs for comparison.

    Returns list of (seed, reference_structure, target_structure).
    For each seed, uses all sample combinations (seed_pool) or
    matched samples (matched_sample).

    Parameters
    ----------
    dataset : Dataset
        The V3 dataset.
    target_condition : str
        The target condition to compare against reference.
    reference_condition : str, optional
        The reference condition. Uses dataset.reference_condition if None.

    Returns
    -------
    list of (seed, ref_structure, target_structure)
    """
    if reference_condition is None:
        reference_condition = dataset.reference_condition

    if reference_condition is None:
        return []

    # Get seeds for both conditions
    ref_seeds = dataset.predictions.get(reference_condition, {})
    target_seeds = dataset.predictions.get(target_condition, {})

    common_seeds = sorted(set(ref_seeds.keys()) & set(target_seeds.keys()))
    pairs = []

    for seed in common_seeds:
        ref_samples = ref_seeds[seed]
        target_samples = target_seeds[seed]

        # All sample combinations (seed_pool)
        for ref_struct in ref_samples.values():
            for target_struct in target_samples.values():
                pairs.append((seed, ref_struct, target_struct))

    return pairs


def get_all_condition_pairs(
    dataset: Dataset,
    reference_condition: Optional[str] = None,
) -> List[Tuple[str, str]]:
    """
    Get all (reference, target) condition pairs for comparison.

    Parameters
    ----------
    dataset : Dataset
        The V3 dataset.
    reference_condition : str, optional
        The reference condition. Uses dataset.reference_condition if None.

    Returns
    -------
    list of (reference_condition, target_condition)
    """
    if reference_condition is None:
        reference_condition = dataset.reference_condition

    if reference_condition is None:
        return []

    pairs = []
    for cond_id in dataset.conditions.keys():
        if cond_id != reference_condition:
            pairs.append((reference_condition, cond_id))

    return pairs


def validate_reference_resolution(
    dataset: Dataset,
    resolution: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Validate reference resolution.

    Checks:
    - Reference condition exists in dataset
    - Common seeds exist for paired conditions
    - At least one comparison can be made
    """
    report = {
        "valid": True,
        "warnings": [],
        "errors": [],
    }

    ref_cond = resolution.get("reference_condition")
    if ref_cond is None:
        report["warnings"].append("No reference condition specified")
        return report

    if ref_cond not in dataset.conditions:
        report["errors"].append(
            f"Reference condition '{ref_cond}' not in dataset"
        )
        report["valid"] = False
        return report

    # Check paired conditions
    for ref, target in resolution.get("paired_conditions", []):
        if ref not in dataset.conditions:
            report["errors"].append(f"Reference '{ref}' not in dataset")
            report["valid"] = False
        if target not in dataset.conditions:
            report["errors"].append(f"Target '{target}' not in dataset")
            report["valid"] = False

        # Check common seeds
        ref_seeds = set(dataset.predictions.get(ref, {}).keys())
        target_seeds = set(dataset.predictions.get(target, {}).keys())
        common = ref_seeds & target_seeds

        if len(common) == 0:
            report["warnings"].append(
                f"No common seeds between '{ref}' and '{target}'"
            )

    return report
