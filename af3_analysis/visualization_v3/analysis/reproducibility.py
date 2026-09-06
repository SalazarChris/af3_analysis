"""
V3 Seed Reproducibility.

For each structural metric calculate:
- mean
- median
- SD
- IQR
- valid seeds
- direction consistency

Treats seeds as prediction-process robustness samples.
Does NOT describe them as biological replicates or physical ensemble.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class SeedReproducibility:
    """Seed reproducibility statistics for one metric."""

    metric_id: str
    condition_id: str
    reference_condition: Optional[str] = None

    # Descriptive stats
    mean: Optional[float] = None
    median: Optional[float] = None
    std: Optional[float] = None
    iqr: Optional[float] = None
    q25: Optional[float] = None
    q75: Optional[float] = None
    min: Optional[float] = None
    max: Optional[float] = None

    # Sample info
    n_seeds: int = 0
    n_valid_seeds: int = 0
    n_comparisons_total: int = 0
    n_comparisons_valid: int = 0

    # Direction consistency
    direction_consistent: Optional[float] = None  # fraction of seeds with same direction
    direction: Optional[str] = None  # "positive", "negative", "mixed"

    # Coverage
    mean_coverage: Optional[float] = None

    # Status
    status: str = "valid"
    reason: str = ""


def calculate_seed_reproducibility(
    seed_values: Dict[int, List[float]],  # seed -> list of values (one per comparison)
    seed_coverages: Optional[Dict[int, List[float]]] = None,
    *,
    metric_id: str = "rmsd",
    condition_id: str = "",
    reference_condition: Optional[str] = None,
    direction_threshold: float = 0.0,
) -> SeedReproducibility:
    """
    Calculate seed reproducibility for a metric.

    Parameters
    ----------
    seed_values : dict mapping seed -> list of values
    seed_coverages : dict mapping seed -> list of coverage values
    metric_id : str
    condition_id : str
    reference_condition : str, optional
    direction_threshold : float
        Minimum value to consider as non-zero direction.

    Returns
    -------
    SeedReproducibility
    """
    # Collect all values with seed info
    all_values = []
    all_seeds = []
    all_coverages = []

    for seed, values in seed_values.items():
        if not values:
            continue

        seed_mean = np.mean(values)
        all_values.append(seed_mean)
        all_seeds.append(seed)

        if seed_coverages and seed in seed_coverages:
            coverages = seed_coverages[seed]
            if coverages:
                all_coverages.append(np.mean(coverages))

    n_seeds = len(seed_values)
    n_valid_seeds = len(all_values)
    n_comparisons_total = sum(len(v) for v in seed_values.values())
    n_comparisons_valid = n_comparisons_total  # Simplified

    # Calculate statistics
    if all_values:
        arr = np.array(all_values)
        mean_val = float(np.mean(arr))
        median_val = float(np.median(arr))
        std_val = float(np.std(arr)) if len(arr) > 1 else None
        q25 = float(np.percentile(arr, 25))
        q75 = float(np.percentile(arr, 75))
        iqr_val = q75 - q25
        min_val = float(np.min(arr))
        max_val = float(np.max(arr))

        # Direction consistency
        if direction_threshold is not None:
            positive = np.sum(arr > direction_threshold)
            negative = np.sum(arr < -direction_threshold)
            total_directional = positive + negative

            if total_directional > 0:
                frac_consistent = max(positive, negative) / total_directional
                direction = "positive" if positive > negative else "negative"
            else:
                frac_consistent = None
                direction = "mixed"
        else:
            frac_consistent = None
            direction = "mixed" if mean_val > 0 else "positive"

        mean_coverage = float(np.mean(all_coverages)) if all_coverages else None

        return SeedReproducibility(
            metric_id=metric_id,
            condition_id=condition_id,
            reference_condition=reference_condition,
            mean=mean_val,
            median=median_val,
            std=std_val,
            iqr=iqr_val,
            q25=q25,
            q75=q75,
            min=min_val,
            max=max_val,
            n_seeds=n_seeds,
            n_valid_seeds=n_valid_seeds,
            n_comparisons_total=n_comparisons_total,
            n_comparisons_valid=n_comparisons_valid,
            direction_consistent=frac_consistent,
            direction=direction,
            mean_coverage=mean_coverage,
            status="valid",
            reason="",
        )
    else:
        return SeedReproducibility(
            metric_id=metric_id,
            condition_id=condition_id,
            reference_condition=reference_condition,
            n_seeds=n_seeds,
            n_valid_seeds=0,
            status="no_valid_seeds",
            reason="No valid seed-level values",
        )


def seed_reproducibility_table(
    reproducibilities: List[SeedReproducibility],
) -> pd.DataFrame:
    """Create DataFrame from seed reproducibility results."""
    rows = []
    for r in reproducibilities:
        row = {
            "metric_id": r.metric_id,
            "condition_id": r.condition_id,
            "reference_condition": r.reference_condition,
            "mean": r.mean,
            "median": r.median,
            "std": r.std,
            "iqr": r.iqr,
            "q25": r.q25,
            "q75": r.q75,
            "min": r.min,
            "max": r.max,
            "n_seeds": r.n_seeds,
            "n_valid_seeds": r.n_valid_seeds,
            "n_comparisons_total": r.n_comparisons_total,
            "n_comparisons_valid": r.n_comparisons_valid,
            "direction_consistent": r.direction_consistent,
            "direction": r.direction,
            "mean_coverage": r.mean_coverage,
            "status": r.status,
        }
        rows.append(row)

    return pd.DataFrame(rows)


def summarize_across_conditions(
    reproducibilities: List[SeedReproducibility],
    *,
    group_by: str = "metric_id",
) -> Dict[str, Any]:
    """
    Summarize seed reproducibility across conditions.

    Parameters
    ----------
    reproducibilities : list of SeedReproducibility
    group_by : str ("metric_id" or "condition_id")

    Returns
    -------
    dict with summary statistics
    """
    df = pd.DataFrame([
        {
            "group": getattr(r, group_by),
            "mean": r.mean,
            "median": r.median,
            "std": r.std,
            "n_seeds": r.n_seeds,
            "n_valid_seeds": r.n_valid_seeds,
            "direction_consistent": r.direction_consistent,
            "direction": r.direction,
        }
        for r in reproducibilities
        if r.mean is not None
    ])

    if df.empty:
        return {}

    grouped = df.groupby("group").agg({
        "mean": ["mean", "std"],
        "median": ["mean", "std"],
        "n_seeds": "sum",
        "n_valid_seeds": "sum",
        "direction_consistent": "mean",
    })

    return {
        "group_by": group_by,
        "n_groups": len(df["group"].unique()),
        "summary": grouped.to_dict(),
    }
