"""
V3 Factorial Structural Analysis.

If condition registry defines a factorial experiment, use those metadata fields.
Does NOT parse condition names.
Estimates structural effects corresponding to:
- main effects
- two-way interactions
- higher-order interactions where supported

Reuses existing statistical methodology.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class FactorialContrast:
    """One factorial contrast result."""

    contrast_id: str  # e.g., "main_DNA", "interaction_DNA_pTPO101"
    contrast_type: str  # "main", "interaction"
    factors: List[str]  # involved factors
    estimate: Optional[float] = None
    std_error: Optional[float] = None
    ci_lower: Optional[float] = None
    ci_upper: Optional[float] = None
    n_observations: int = 0
    n_conditions: int = 0

    # Status
    status: str = "valid"
    reason: str = ""


@dataclass
class FactorialStructuralAnalysis:
    """Complete factorial structural analysis result."""

    contrasts: List[FactorialContrast] = field(default_factory=list)
    factors: List[str] = field(default_factory=list)
    factor_levels: Dict[str, List[Any]] = field(default_factory=dict)

    # Metadata
    n_conditions: int = 0
    n_observations: int = 0
    design_type: str = ""  # "complete_factorial", "incomplete", "unknown"

    # Warnings
    warnings: List[str] = field(default_factory=list)


def factorial_structural_analysis(
    rmsd_values: Dict[str, List[float]],  # condition_id -> list of RMSD values
    metadata: Dict[str, Any],  # experiment metadata with factors
    *,
    metric_name: str = "rmsd",
) -> FactorialStructuralAnalysis:
    """
    Perform factorial analysis on structural metric.

    Parameters
    ----------
    rmsd_values : dict mapping condition_id -> list of RMSD values
    metadata : dict with experiment design
    metric_name : str

    Returns
    -------
    FactorialStructuralAnalysis
    """
    warnings = []
    contrasts = []

    # Extract factors from metadata
    factors = list(metadata.get("attributes", {}).keys())
    factor_levels = {}

    for attr_name, attr_def in metadata.get("attributes", {}).items():
        attr_type = attr_def.get("type", "binary")
        if attr_type == "binary":
            factor_levels[attr_name] = [False, True]
        else:
            # Categorical: get from conditions
            levels = set()
            for cond_attrs in metadata.get("conditions", {}).values():
                if attr_name in cond_attrs.get("attributes", {}):
                    levels.add(cond_attrs["attributes"][attr_name])
            factor_levels[attr_name] = sorted(levels, key=str)

    # Get conditions with their factor levels
    condition_factors = {}
    for cond_name, cond_data in metadata.get("conditions", {}).items():
        cond_factors = {}
        for attr_name in factors:
            if attr_name in cond_data.get("attributes", {}):
                cond_factors[attr_name] = cond_data["attributes"][attr_name]
        condition_factors[cond_name] = cond_factors

    # Calculate main effects for each factor
    for factor in factors:
        levels = factor_levels.get(factor, [])
        if not levels:
            continue

        # Get RMSD values for each level
        level_values = {}
        for level in levels:
            level_values[level] = []

        for cond_id, values in rmsd_values.items():
            cond_factors = condition_factors.get(cond_id, {})
            level = cond_factors.get(factor)
            if level is not None and values:
                level_values[level].extend(values)

        # Calculate main effect as difference between levels
        if len(levels) == 2:
            level_a = levels[0]
            level_b = levels[1]
            values_a = level_values.get(level_a, [])
            values_b = level_values.get(level_b, [])

            if values_a and values_b:
                mean_a = np.mean(values_a)
                mean_b = np.mean(values_b)
                effect = mean_b - mean_a

                n_a = len(values_a)
                n_b = len(values_b)

                # Standard error
                var_a = np.var(values_a, ddof=1) if n_a > 1 else 0
                var_b = np.var(values_b, ddof=1) if n_b > 1 else 0
                se = np.sqrt(var_a / n_a + var_b / n_b) if n_a > 0 and n_b > 0 else None

                contrast = FactorialContrast(
                    contrast_id=f"main_{factor}",
                    contrast_type="main",
                    factors=[factor],
                    estimate=float(effect),
                    std_error=float(se) if se else None,
                    n_observations=n_a + n_b,
                    n_conditions=len(set(
                        cond_id for cond_id, vals in rmsd_values.items()
                        if condition_factors.get(cond_id, {}).get(factor) in [level_a, level_b]
                    )),
                )

                if se:
                    from scipy import stats
                    df = n_a + n_b - 2
                    t_crit = stats.t.ppf(0.975, df)
                    contrast.ci_lower = float(contrast.estimate - t_crit * se)
                    contrast.ci_upper = float(contrast.estimate + t_crit * se)

                contrasts.append(contrast)

    # Calculate two-way interactions
    if len(factors) >= 2:
        from itertools import combinations
        for factor_a, factor_b in combinations(factors, 2):
            # Simple interaction: difference of differences
            # This is a simplified approach
            levels_a = factor_levels.get(factor_a, [])
            levels_b = factor_levels.get(factor_b, [])

            if len(levels_a) == 2 and len(levels_b) == 2:
                # Get 4 combinations
                combos = {}
                for cond_id, values in rmsd_values.items():
                    cf = condition_factors.get(cond_id, {})
                    key = (cf.get(factor_a), cf.get(factor_b))
                    if key not in combos:
                        combos[key] = []
                    combos[key].extend(values)

                # Check all 4 combos exist
                all_keys = [(la, lb) for la in levels_a for lb in levels_b]
                if all(k in combos and combos[k] for k in all_keys):
                    # Calculate interaction as (A1B1 - A2B1) - (A1B2 - A2B2)
                    a1b1 = np.mean(combos[(levels_a[0], levels_b[0])])
                    a2b1 = np.mean(combos[(levels_a[1], levels_b[0])])
                    a1b2 = np.mean(combos[(levels_a[0], levels_b[1])])
                    a2b2 = np.mean(combos[(levels_a[1], levels_b[1])])

                    interaction = (a1b1 - a2b1) - (a1b2 - a2b2)

                    contrast = FactorialContrast(
                        contrast_id=f"interaction_{factor_a}_{factor_b}",
                        contrast_type="interaction",
                        factors=[factor_a, factor_b],
                        estimate=float(interaction),
                        n_observations=sum(len(combos[k]) for k in all_keys),
                        n_conditions=len([k for k in all_keys if combos[k]]),
                    )
                    contrasts.append(contrast)

    return FactorialStructuralAnalysis(
        contrasts=contrasts,
        factors=factors,
        factor_levels=factor_levels,
        n_conditions=len(rmsd_values),
        n_observations=sum(len(v) for v in rmsd_values.values()),
        design_type="complete_factorial" if len(contrasts) > 0 else "unknown",
        warnings=warnings,
    )


def factorial_contrasts_table(
    analysis: FactorialStructuralAnalysis,
) -> pd.DataFrame:
    """Create DataFrame from factorial analysis."""
    rows = []
    for contrast in analysis.contrasts:
        row = {
            "contrast_id": contrast.contrast_id,
            "contrast_type": contrast.contrast_type,
            "factors": "+".join(contrast.factors),
            "estimate": contrast.estimate,
            "std_error": contrast.std_error,
            "ci_lower": contrast.ci_lower,
            "ci_upper": contrast.ci_upper,
            "n_observations": contrast.n_observations,
            "n_conditions": contrast.n_conditions,
            "status": contrast.status,
        }
        rows.append(row)

    return pd.DataFrame(rows)
