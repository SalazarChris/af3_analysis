"""
V3 Structural Effect Sizes.

For every condition/reference comparison calculate:
- RMSD effect
- local displacement effect
- contact-change effect
- interface-change effect
- domain-motion effect

Reuses existing statistical framework where possible.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class StructuralEffectSize:
    """Effect size for one structural metric."""

    metric_id: str
    condition_a: str
    condition_b: str
    seed: int

    # Point estimate
    estimate: Optional[float] = None

    # Uncertainty
    std_error: Optional[float] = None
    confidence_interval_lower: Optional[float] = None
    confidence_interval_upper: Optional[float] = None
    confidence_level: float = 0.95

    # Direction
    direction: Optional[str] = None  # "increased", "decreased", "unchanged"
    direction_consistent: Optional[bool] = None  # across seeds?

    # Coverage
    n_valid: int = 0
    n_comparisons: int = 0
    coverage_mean: Optional[float] = None

    # Effect size measures
    cohens_d: Optional[float] = None
    hedges_g: Optional[float] = None

    # Status
    status: str = "valid"
    reason: str = ""


def calculate_structural_effect_sizes(
    rmsd_values: List[float],
    coverage_values: List[float],
    *,
    condition_a: str,
    condition_b: str,
    reference_condition: Optional[str] = None,
    confidence_level: float = 0.95,
) -> List[StructuralEffectSize]:
    """
    Calculate structural effect sizes.

    Parameters
    ----------
    rmsd_values : list of RMSD values (one per comparison)
    coverage_values : list of coverage values
    condition_a, condition_b : str
    reference_condition : str, optional
    confidence_level : float

    Returns
    -------
    list of StructuralEffectSize
    """
    if not rmsd_values:
        return []

    # Calculate basic statistics
    rmsd_arr = np.array(rmsd_values)
    coverage_arr = np.array(coverage_values) if coverage_values else np.array([])

    n = len(rmsd_arr)
    mean_rmsd = float(np.mean(rmsd_arr))
    std_rmsd = float(np.std(rmsd_arr)) if n > 1 else None
    median_rmsd = float(np.median(rmsd_arr))

    # Calculate effect size (Cohen's d relative to zero)
    # For RMSD, the "effect" is the deviation from zero (no difference)
    if std_rmsd is not None and std_rmsd > 0:
        cohens_d = mean_rmsd / std_rmsd
        # Hedges' g correction for small samples
        if n > 1:
            correction = 1 - 3 / (4 * n - 9)
            hedges_g = cohens_d * correction
        else:
            hedges_g = cohens_d
    else:
        cohens_d = None
        hedges_g = None

    # Standard error
    std_error = std_rmsd / np.sqrt(n) if std_rmsd is not None and n > 1 else None

    # Confidence interval
    ci_lower = None
    ci_upper = None
    if std_error is not None:
        from scipy import stats
        t_critical = stats.t.ppf((1 + confidence_level) / 2, df=n - 1)
        ci_lower = mean_rmsd - t_critical * std_error
        ci_upper = mean_rmsd + t_critical * std_error

    # Effect size object
    effect = StructuralEffectSize(
        metric_id="rmsd_global_ca",
        condition_a=condition_a,
        condition_b=condition_b,
        seed=0,  # Aggregated across seeds
        estimate=mean_rmsd,
        std_error=std_error,
        confidence_interval_lower=ci_lower,
        confidence_interval_upper=ci_upper,
        confidence_level=confidence_level,
        direction="increased" if mean_rmsd > 0 else ("decreased" if mean_rmsd < 0 else "unchanged"),
        direction_consistent=None,  # Would need seed-level direction
        n_valid=n,
        n_comparisons=n,
        coverage_mean=float(np.mean(coverage_arr)) if len(coverage_arr) > 0 else None,
        cohens_d=cohens_d,
        hedges_g=hedges_g,
        status="valid" if n > 0 else "no_data",
        reason="" if n > 0 else "no_rmsd_values",
    )

    return [effect]


def calculate_local_displacement_effect(
    displacement_values: List[float],
    *,
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """
    Calculate effect size for local displacement.

    Parameters
    ----------
    displacement_values : list of mean displacement values
    threshold : float
        Minimum displacement to consider "real"

    Returns
    -------
    dict with effect statistics
    """
    if not displacement_values:
        return {"status": "no_data"}

    arr = np.array(displacement_values)
    n = len(arr)
    mean_disp = float(np.mean(arr))
    std_disp = float(np.std(arr)) if n > 1 else None

    # Fraction of residues with displacement > threshold
    n_above = np.sum(arr > threshold)
    frac_above = n_above / n if n > 0 else 0.0

    # Cohen's d
    cohens_d = mean_disp / std_disp if std_disp and std_disp > 0 else None

    return {
        "mean_displacement": mean_disp,
        "std_displacement": std_disp,
        "median_displacement": float(np.median(arr)),
        "n_residues": n,
        "n_above_threshold": int(n_above),
        "fraction_above_threshold": frac_above,
        "cohens_d": cohens_d,
        "status": "valid" if n > 0 else "no_data",
    }


def calculate_contact_change_effect(
    n_gained: int,
    n_lost: int,
    n_total: int,
    *,
    pct_threshold: float = 0.05,
) -> Dict[str, Any]:
    """
    Calculate effect size for contact changes.

    Parameters
    ----------
    n_gained : int
    n_lost : int
    n_total : int
    pct_threshold : float

    Returns
    -------
    dict with effect statistics
    """
    if n_total == 0:
        return {"status": "no_data"}

    pct_gained = n_gained / n_total
    pct_lost = n_lost / n_total
    pct_changed = (n_gained + n_lost) / n_total

    # Net change
    net_change = n_gained - n_lost

    return {
        "n_gained": n_gained,
        "n_lost": n_lost,
        "n_total": n_total,
        "pct_gained": pct_gained,
        "pct_lost": pct_lost,
        "pct_changed": pct_changed,
        "net_change": net_change,
        "significant": pct_changed > pct_threshold,
        "status": "valid",
    }
