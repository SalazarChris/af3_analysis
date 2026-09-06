"""
V3 Matched-Seed Structural Analysis.

For every condition/reference pair:
- seed 1 → compare
- seed 2 → compare
...
- seed N → compare

Produces:
- RMSD per seed
- coverage per seed
- valid seed count
- mean
- median
- SD
- IQR
- direction consistency where applicable
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..model import Dataset, StructureData


@dataclass
class MatchedSeedResult:
    """Result for matched-seed analysis."""

    condition_a: str
    condition_b: str
    seed: int
    n_comparisons: int
    n_valid: int
    n_incomparable: int

    # RMSD stats
    rmsd_values: List[float] = field(default_factory=list)
    rmsd_mean: Optional[float] = None
    rmsd_median: Optional[float] = None
    rmsd_std: Optional[float] = None
    rmsd_iqr: Optional[float] = None

    # Coverage stats
    coverage_values: List[float] = field(default_factory=list)
    coverage_mean: Optional[float] = None

    # Displacement stats
    displacement_mean: Optional[float] = None
    displacement_median: Optional[float] = None

    # Status
    status: str = "valid"
    reason: str = ""


def analyze_matched_seeds(
    dataset: Dataset,
    condition_a: str,
    condition_b: str,
    *,
    reference_condition: Optional[str] = None,
    alignment_atom: str = "CA",
    min_common_atoms: int = 3,
    min_coverage: float = 0.80,
    min_sequence_identity: float = 0.5,
) -> List[MatchedSeedResult]:
    """
    Run matched-seed structural analysis between two conditions.

    Parameters
    ----------
    dataset : Dataset
    condition_a, condition_b : str
        Conditions to compare.
    reference_condition : str, optional
        Reference condition. Uses dataset.reference_condition if None.
    alignment_atom : str
    min_common_atoms : int
    min_coverage : float
    min_sequence_identity : float

    Returns
    -------
    list of MatchedSeedResult, one per common seed
    """
    if reference_condition is None:
        reference_condition = dataset.reference_condition

    # Determine which is reference and which is target
    if reference_condition == condition_a:
        ref_condition = condition_a
        target_condition = condition_b
    elif reference_condition == condition_b:
        ref_condition = condition_b
        target_condition = condition_a
    else:
        # No explicit reference; use first condition
        ref_condition = condition_a
        target_condition = condition_b

    # Get common seeds
    ref_seeds = dataset.predictions.get(ref_condition, {})
    target_seeds = dataset.predictions.get(target_condition, {})

    common_seeds = sorted(set(ref_seeds.keys()) & set(target_seeds.keys()))

    if not common_seeds:
        return []

    results = []

    for seed in common_seeds:
        ref_samples = ref_seeds[seed]
        target_samples = target_seeds[seed]

        rmsd_values = []
        coverage_values = []

        n_comparisons = 0
        n_valid = 0
        n_incomparable = 0

        # Compare all sample combinations
        for ref_struct in ref_samples.values():
            for target_struct in target_samples.values():
                n_comparisons += 1

                # Import RMSD calculation
                from ..structural.rmsd import calculate_rmsd

                result = calculate_rmsd(
                    target_struct,
                    ref_struct,
                    alignment_atom=alignment_atom,
                    min_common_atoms=min_common_atoms,
                    min_sequence_identity=min_sequence_identity,
                    min_coverage=min_coverage,
                )

                if result["status"] == "comparable" and result["rmsd"] is not None:
                    rmsd_values.append(result["rmsd"])
                    coverage_values.append(result["coverage"])
                    n_valid += 1
                else:
                    n_incomparable += 1

        # Calculate stats
        if rmsd_values:
            rmsd_arr = np.array(rmsd_values)
            rmsd_mean = float(np.mean(rmsd_arr))
            rmsd_median = float(np.median(rmsd_arr))
            rmsd_std = float(np.std(rmsd_arr))
            rmsd_iqr = float(np.percentile(rmsd_arr, 75) - np.percentile(rmsd_arr, 25))
            coverage_mean = float(np.mean(coverage_values))
        else:
            rmsd_mean = rmsd_median = rmsd_std = rmsd_iqr = coverage_mean = None

        # Determine status
        if n_valid == 0:
            status = "no_valid_comparisons"
            reason = f"No valid comparisons for seed {seed}"
        elif n_valid < n_comparisons:
            status = "partial"
            reason = f"{n_valid}/{n_comparisons} valid comparisons"
        else:
            status = "valid"
            reason = ""

        results.append(MatchedSeedResult(
            condition_a=condition_a,
            condition_b=condition_b,
            seed=seed,
            n_comparisons=n_comparisons,
            n_valid=n_valid,
            n_incomparable=n_incomparable,
            rmsd_values=rmsd_values,
            rmsd_mean=rmsd_mean,
            rmsd_median=rmsd_median,
            rmsd_std=rmsd_std,
            rmsd_iqr=rmsd_iqr,
            coverage_values=coverage_values,
            coverage_mean=coverage_mean,
            status=status,
            reason=reason,
        ))

    return results


def aggregate_matched_seed_results(
    results: List[MatchedSeedResult],
) -> Dict[str, Any]:
    """
    Aggregate matched-seed results across seeds.

    Returns seed-level summary statistics.
    """
    if not results:
        return {
            "n_seeds": 0,
            "n_valid_seeds": 0,
            "status": "no_data",
        }

    # Collect seed-level RMSD means
    seed_means = []
    seed_coverages = []
    seed_valid = []

    for r in results:
        if r.rmsd_mean is not None:
            seed_means.append(r.rmsd_mean)
            seed_coverages.append(r.coverage_mean)
            seed_valid.append(r.n_valid)

    if not seed_means:
        return {
            "n_seeds": len(results),
            "n_valid_seeds": 0,
            "status": "no_valid_data",
        }

    rmsd_arr = np.array(seed_means)
    coverage_arr = np.array([c for c in seed_coverages if c is not None])

    return {
        "n_seeds": len(results),
        "n_valid_seeds": len(seed_means),
        "rmsd_mean_of_means": float(np.mean(rmsd_arr)),
        "rmsd_median_of_means": float(np.median(rmsd_arr)),
        "rmsd_std_of_means": float(np.std(rmsd_arr)),
        "rmsd_iqr_of_means": float(np.percentile(rmsd_arr, 75) - np.percentile(rmsd_arr, 25)),
        "coverage_mean": float(np.mean(coverage_arr)) if len(coverage_arr) > 0 else None,
        "total_comparisons": sum(r.n_comparisons for r in results),
        "total_valid": sum(r.n_valid for r in results),
        "status": "valid" if len(seed_means) > 0 else "no_valid_data",
    }


def matched_seed_summary_table(
    results: List[MatchedSeedResult],
) -> pd.DataFrame:
    """Create summary table from matched-seed results."""
    if not results:
        return pd.DataFrame()

    rows = []
    for r in results:
        row = {
            "condition_a": r.condition_a,
            "condition_b": r.condition_b,
            "seed": r.seed,
            "n_comparisons": r.n_comparisons,
            "n_valid": r.n_valid,
            "n_incomparable": r.n_incomparable,
            "rmsd_mean": r.rmsd_mean,
            "rmsd_median": r.rmsd_median,
            "rmsd_std": r.rmsd_std,
            "rmsd_iqr": r.rmsd_iqr,
            "coverage_mean": r.coverage_mean,
            "status": r.status,
            "reason": r.reason,
        }
        rows.append(row)

    return pd.DataFrame(rows)
