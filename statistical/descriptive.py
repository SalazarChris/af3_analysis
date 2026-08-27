"""
Descriptive statistics for AF3 Confidence Analysis Pipeline.

Provides condition and variability summaries from seed-level data.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np


@dataclass
class ConditionSummary:
    """Summary for a single condition."""
    condition_id: str
    n_seeds: int
    n_samples: int
    seed_means: Dict[str, float]
    seed_sds: Dict[str, float]
    sample_means: Dict[str, float]
    sample_sds: Dict[str, float]


@dataclass
class VariabilitySummary:
    """Variability decomposition summary."""
    metric_id: str
    scope: str
    between_condition_var: float
    between_seed_var: float
    within_seed_var: float
    total_var: float
    cv_between_seed: float  # Coefficient of variation
    cv_within_seed: float


class DescriptiveStats:
    """
    Compute descriptive statistics from analysis_seed data.
    """
    
    def __init__(self):
        self._condition_summaries: List[Dict[str, Any]] = []
        self._variability_summaries: List[Dict[str, Any]] = []
    
    def compute_condition_summaries(
        self,
        analysis_seed_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Compute condition-level descriptive summaries.
        
        Args:
            analysis_seed_df: Analysis seed table
        
        Returns:
            DataFrame with condition summaries
        """
        summaries = []
        
        for (cond_id, metric_id, scope_type, scope_id), group in analysis_seed_df.groupby(
            ["condition_id", "metric_id", "scope_type", "scope_id"]
        ):
            values = group["value"].dropna()
            
            summaries.append({
                "condition_id": cond_id,
                "metric_id": metric_id,
                "scope_type": scope_type,
                "scope_id": scope_id,
                "n_seeds": len(values),
                "mean": values.mean() if len(values) > 0 else None,
                "sd": values.std(ddof=1) if len(values) > 1 else None,
                "median": values.median() if len(values) > 0 else None,
                "min": values.min() if len(values) > 0 else None,
                "max": values.max() if len(values) > 0 else None,
                "seeds": group["seed"].unique().tolist(),
            })
        
        return pd.DataFrame(summaries)
    
    def compute_variability_summaries(
        self,
        analysis_seed_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Compute variability decomposition summaries.
        
        Decomposes variance into:
        - Between-condition variance
        - Between-seed variance
        - Within-seed variance
        
        Args:
            analysis_seed_df: Analysis seed table
        
        Returns:
            DataFrame with variability summaries
        """
        summaries = []
        
        for (metric_id, scope_type, scope_id), group in analysis_seed_df.groupby(
            ["metric_id", "scope_type", "scope_id"]
        ):
            # Between-seed variance
            seed_means = group.groupby("condition_id")["value"].mean()
            between_seed_var = seed_means.var(ddof=1) if len(seed_means) > 1 else 0
            
            # Between-condition variance (from seed means)
            between_cond_var = seed_means.var(ddof=1) if len(seed_means) > 1 else 0
            
            # Within-seed variance (from sample_sd in analysis_seed)
            sample_sds = group["sample_sd"].dropna()
            within_seed_var = (sample_sds ** 2).mean() if len(sample_sds) > 0 else 0
            
            # Total variance
            total_var = between_seed_var + within_seed_var
            
            # Coefficient of variation
            mean_val = seed_means.mean() if len(seed_means) > 0 else 0
            cv_between_seed = (np.sqrt(between_seed_var) / abs(mean_val)) if mean_val != 0 else 0
            cv_within_seed = (np.sqrt(within_seed_var) / abs(mean_val)) if mean_val != 0 else 0
            
            summaries.append({
                "metric_id": metric_id,
                "scope_type": scope_type,
                "scope_id": scope_id,
                "n_conditions": group["condition_id"].nunique(),
                "n_seeds": len(group),
                "between_seed_var": between_seed_var,
                "between_cond_var": between_cond_var,
                "within_seed_var": within_seed_var,
                "total_var": total_var,
                "cv_between_seed": cv_between_seed,
                "cv_within_seed": cv_within_seed,
            })
        
        return pd.DataFrame(summaries)
    
    def compute_bootstrap_ci(
        self,
        values: np.ndarray,
        n_bootstrap: int = 2000,
        confidence_level: float = 0.95,
        random_state: Optional[int] = None,
    ) -> Dict[str, float]:
        """
        Compute bootstrap confidence intervals.
        
        Args:
            values: Sample values
            n_bootstrap: Number of bootstrap iterations
            confidence_level: Confidence level (0-1)
            random_state: Random seed for reproducibility
        
        Returns:
            Dict with lower, upper, and width
        """
        if len(values) == 0:
            return {"lower": None, "upper": None, "width": None}
        
        rng = np.random.RandomState(random_state)
        
        # Bootstrap distribution of mean
        bootstrap_means = []
        for _ in range(n_bootstrap):
            sample = rng.choice(values, size=len(values), replace=True)
            bootstrap_means.append(sample.mean())
        
        bootstrap_means = np.array(bootstrap_means)
        
        # Compute confidence interval
        alpha = 1 - confidence_level
        lower = float(np.percentile(bootstrap_means, alpha / 2 * 100))
        upper = float(np.percentile(bootstrap_means, (1 - alpha / 2) * 100))
        
        return {
            "lower": lower,
            "upper": upper,
            "width": upper - lower,
            "mean": float(np.mean(bootstrap_means)),
            "ci_lower": lower,
            "ci_upper": upper,
        }


def compute_condition_summaries(
    analysis_seed_df: pd.DataFrame,
) -> pd.DataFrame:
    """Convenience function for condition summaries."""
    stats = DescriptiveStats()
    return stats.compute_condition_summaries(analysis_seed_df)


def compute_variability_summaries(
    analysis_seed_df: pd.DataFrame,
) -> pd.DataFrame:
    """Convenience function for variability summaries."""
    stats = DescriptiveStats()
    return stats.compute_variability_summaries(analysis_seed_df)


__all__ = [
    "ConditionSummary",
    "VariabilitySummary",
    "DescriptiveStats",
    "compute_condition_summaries",
    "compute_variability_summaries",
]
