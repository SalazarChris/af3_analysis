"""Variance decomposition and reference spread for AF3 Confidence Analysis Pipeline."""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any
from dataclasses import dataclass


def seed_means(long_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute mean over samples within seed for each condition/metric/scope.
    
    Args:
        long_df: Long-format measurements with columns:
            condition_id, seed, sample, metric_id, scope_type, scope_id, value
    
    Returns:
        DataFrame with seed-level means, n_samples, and sample SD
    """
    grouped = long_df.groupby(["condition_id", "seed", "metric_id", "scope_type", "scope_id"])
    
    results = []
    for key, group in grouped:
        condition_id, seed, metric_id, scope_type, scope_id = key
        values = group["value"].dropna()
        
        results.append({
            "condition_id": condition_id,
            "seed": seed,
            "metric_id": metric_id,
            "scope_type": scope_type,
            "scope_id": scope_id,
            "value": values.mean() if len(values) > 0 else None,
            "n_samples": len(values),
            "sample_sd": values.std(ddof=1) if len(values) > 1 else None,
        })
    
    return pd.DataFrame(results)


def variance_components(long_df: pd.DataFrame, by: tuple = ("condition_id", "metric_id", "scope_type", "scope_id")) -> pd.DataFrame:
    """
    Decompose variance into within-seed and between-seed components.
    
    Args:
        long_df: Long-format measurements
        by: Grouping columns for aggregation
    
    Returns:
        DataFrame with variance components per group
    """
    # First aggregate to seed level
    seed_df = seed_means(long_df)
    
    results = []
    for key, group in seed_df.groupby(list(by)):
        if isinstance(key, tuple) and len(by) == 1:
            key = (key,)
        
        values = group["value"].dropna()
        
        if len(values) < 2:
            results.append({
                **dict(zip(by, key if isinstance(key, tuple) else (key,))),
                "n_seeds": len(values),
                "within_seed_var": 0.0,
                "between_seed_var": 0.0,
                "between_seed_share": 0.0,
            })
            continue
        
        # Within-seed variance (from sample_sd column)
        within_vars = (group["sample_sd"].dropna() ** 2).dropna()
        within_seed_var = float(within_vars.mean()) if len(within_vars) > 0 else 0.0
        
        # Between-seed variance
        between_seed_var = float(values.var(ddof=1))
        
        # Total variance
        total_var = between_seed_var + within_seed_var
        
        results.append({
            **dict(zip(by, key if isinstance(key, tuple) else (key,))),
            "n_seeds": len(values),
            "within_seed_var": within_seed_var,
            "between_seed_var": between_seed_var,
            "between_seed_share": between_seed_var / total_var if total_var > 0 else 0.0,
        })
    
    return pd.DataFrame(results)


def reference_spread(seed_means_df: pd.DataFrame, reference_id: str) -> pd.DataFrame:
    """
    Compute reference condition's seed-mean SD and IQR per metric/scope.
    
    Args:
        seed_means_df: Seed-level values from seed_means()
        reference_id: Reference condition ID
    
    Returns:
        DataFrame with reference condition's SD and IQR per metric/scope
    """
    ref_df = seed_means_df[seed_means_df["condition_id"] == reference_id]
    
    results = []
    for key, group in ref_df.groupby(["metric_id", "scope_type", "scope_id"]):
        metric_id, scope_type, scope_id = key
        values = group["value"].dropna()
        
        if len(values) == 0:
            results.append({
                "metric_id": metric_id,
                "scope_type": scope_type,
                "scope_id": scope_id,
                "seed_mean_sd": None,
                "seed_mean_iqr": None,
                "n_seeds": 0,
            })
            continue
        
        sd = values.std(ddof=1)
        q75 = values.quantile(0.75)
        q25 = values.quantile(0.25)
        iqr = q75 - q25
        
        results.append({
            "metric_id": metric_id,
            "scope_type": scope_type,
            "scope_id": scope_id,
            "seed_mean_sd": sd,
            "seed_mean_iqr": iqr,
            "n_seeds": len(values),
        })
    
    return pd.DataFrame(results)


def min_detectable_difference(
    spread_df: pd.DataFrame,
    n_seeds: int,
    alpha: float = 0.05,
    power: float = 0.8,
) -> pd.DataFrame:
    """
    Estimate minimum detectable difference per metric/scope.
    
    Uses approximation: MDD ≈ 2.9 * SD / sqrt(n) for rank-test case.
    
    Args:
        spread_df: Reference spread from reference_spread()
        n_seeds: Number of seeds per condition
        alpha: Significance level
        power: Statistical power
    
    Returns:
        DataFrame with MDD per metric/scope
    """
    results = []
    for _, row in spread_df.iterrows():
        sd = row["seed_mean_sd"]
        if sd is None or np.isnan(sd):
            mdd = None
        else:
            # Approximate MDD for two-sided rank test
            mdd = 2.9 * sd / np.sqrt(n_seeds)
        
        results.append({
            "metric_id": row["metric_id"],
            "scope_type": row["scope_type"],
            "scope_id": row["scope_id"],
            "n_seeds": n_seeds,
            "min_detectable_difference": mdd,
        })
    
    return pd.DataFrame(results)


def degenerate_metrics(seed_means_df: pd.DataFrame) -> pd.DataFrame:
    """
    Find metrics with zero or near-zero spread across the design.
    
    Args:
        seed_means_df: Seed-level values from seed_means()
    
    Returns:
        DataFrame with degenerate metrics and their spread
    """
    results = []
    for key, group in seed_means_df.groupby(["metric_id", "scope_type", "scope_id"]):
        metric_id, scope_type, scope_id = key
        values = group["value"].dropna()
        
        if len(values) < 2:
            spread = 0.0
        else:
            spread = values.std(ddof=1)
        
        if spread < 1e-10:  # Near-zero threshold
            results.append({
                "metric_id": metric_id,
                "scope_type": scope_type,
                "scope_id": scope_id,
                "spread": spread,
                "n_values": len(values),
            })
    
    return pd.DataFrame(results) if results else pd.DataFrame(columns=["metric_id", "scope_type", "scope_id", "spread", "n_values"])


__all__ = [
    "seed_means",
    "variance_components",
    "reference_spread",
    "min_detectable_difference",
    "degenerate_metrics",
]
