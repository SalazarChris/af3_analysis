"""
Variance decomposition for AF3 Confidence Analysis Pipeline.

Computes variance components: within-seed, between-seed, between-condition.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np


@dataclass
class VarianceComponents:
    """Variance decomposition components."""
    metric_id: str
    scope: str
    between_condition_var: float
    between_seed_var: float
    within_seed_var: float
    total_var: float
    fraction_between_condition: float
    fraction_between_seed: float
    fraction_within_seed: float


@dataclass
class VarianceSummary:
    """Summary of variance decomposition."""
    metric_id: str
    scope: str
    n_conditions: int
    n_seeds: int
    between_seed_sd: float
    within_seed_sd: float
    total_sd: float


class VarianceDecomposer:
    """
    Decompose variance into components.
    
    - Within-seed sampling variation
    - Between-seed variation within condition
    - Factor-model explained vs residual variation
    """
    
    def __init__(self):
        self._components: List[Dict[str, Any]] = []
    
    def decompose_variance(
        self,
        analysis_seed_df: pd.DataFrame,
        metric_id: str,
        scope_type: str = "global",
        scope_id: str = "",
    ) -> VarianceComponents:
        """
        Decompose variance for a metric and scope.
        
        Args:
            analysis_seed_df: Analysis seed table
            metric_id: Metric to analyze
            scope_type: Scope type
            scope_id: Scope ID
        
        Returns:
            VarianceComponents
        """
        # Handle empty DataFrame
        if analysis_seed_df.empty:
            return VarianceComponents(
                metric_id=metric_id,
                scope=f"{scope_type}_{scope_id}",
                between_condition_var=np.nan,
                between_seed_var=np.nan,
                within_seed_var=np.nan,
                total_var=np.nan,
                fraction_between_condition=np.nan,
                fraction_between_seed=np.nan,
                fraction_within_seed=np.nan,
            )
        
        # Filter to the specific metric and scope
        subset = analysis_seed_df[
            (analysis_seed_df["metric_id"] == metric_id) &
            (analysis_seed_df["scope_type"] == scope_type) &
            (analysis_seed_df["scope_id"] == scope_id)
        ]
        
        if subset.empty:
            return VarianceComponents(
                metric_id=metric_id,
                scope=f"{scope_type}_{scope_id}",
                between_condition_var=np.nan,
                between_seed_var=np.nan,
                within_seed_var=np.nan,
                total_var=np.nan,
                fraction_between_condition=np.nan,
                fraction_between_seed=np.nan,
                fraction_within_seed=np.nan,
            )
        
        # Get seed means per condition
        seed_means = subset.groupby("condition_id")["value"].mean()
        
        # Between-condition variance
        between_condition_var = float(seed_means.var(ddof=1)) if len(seed_means) > 1 else 0
        
        # Between-seed variance (from within-condition seed variance)
        # For now, use overall between-seed variance
        values = subset["value"].dropna()
        between_seed_var = float(values.var(ddof=1)) if len(values) > 1 else 0
        
        # Within-seed variance (from sample_sd in analysis_seed)
        sample_sds = subset["sample_sd"].dropna()
        within_seed_var = float((sample_sds ** 2).mean()) if len(sample_sds) > 0 else 0
        
        # Total variance
        total_var = between_condition_var + between_seed_var + within_seed_var
        
        # Fractions
        if total_var > 0:
            fraction_between_condition = between_condition_var / total_var
            fraction_between_seed = between_seed_var / total_var
            fraction_within_seed = within_seed_var / total_var
        else:
            fraction_between_condition = 0
            fraction_between_seed = 0
            fraction_within_seed = 0
        
        return VarianceComponents(
            metric_id=metric_id,
            scope=f"{scope_type}_{scope_id}",
            between_condition_var=between_condition_var,
            between_seed_var=between_seed_var,
            within_seed_var=within_seed_var,
            total_var=total_var,
            fraction_between_condition=fraction_between_condition,
            fraction_between_seed=fraction_between_seed,
            fraction_within_seed=fraction_within_seed,
        )
    
    def compute_variance_summary(
        self,
        analysis_seed_df: pd.DataFrame,
        metric_id: str,
        scope_type: str = "global",
        scope_id: str = "",
    ) -> VarianceSummary:
        """
        Compute variance summary for a metric.
        
        Args:
            analysis_seed_df: Analysis seed table
            metric_id: Metric to analyze
            scope_type: Scope type
            scope_id: Scope ID
        
        Returns:
            VarianceSummary
        """
        subset = analysis_seed_df[
            (analysis_seed_df["metric_id"] == metric_id) &
            (analysis_seed_df["scope_type"] == scope_type) &
            (analysis_seed_df["scope_id"] == scope_id)
        ]
        
        values = subset["value"].dropna()
        sample_sds = subset["sample_sd"].dropna()
        
        between_seed_sd = float(values.std(ddof=1)) if len(values) > 1 else 0
        within_seed_sd = float(sample_sds.mean()) if len(sample_sds) > 0 else 0
        total_sd = float(np.sqrt(between_seed_sd**2 + within_seed_sd**2))
        
        return VarianceSummary(
            metric_id=metric_id,
            scope=f"{scope_type}_{scope_id}",
            n_conditions=subset["condition_id"].nunique(),
            n_seeds=len(values),
            between_seed_sd=between_seed_sd,
            within_seed_sd=within_seed_sd,
            total_sd=total_sd,
        )


def decompose_variance(
    analysis_seed_df: pd.DataFrame,
    metric_id: str,
    scope_type: str = "global",
    scope_id: str = "",
) -> VarianceComponents:
    """Convenience function for variance decomposition."""
    decomposer = VarianceDecomposer()
    return decomposer.decompose_variance(analysis_seed_df, metric_id, scope_type, scope_id)


def compute_variance_summary(
    analysis_seed_df: pd.DataFrame,
    metric_id: str,
    scope_type: str = "global",
    scope_id: str = "",
) -> VarianceSummary:
    """Convenience function for variance summary."""
    decomposer = VarianceDecomposer()
    return decomposer.compute_variance_summary(analysis_seed_df, metric_id, scope_type, scope_id)


__all__ = [
    "VarianceComponents",
    "VarianceSummary",
    "VarianceDecomposer",
    "decompose_variance",
    "compute_variance_summary",
]
