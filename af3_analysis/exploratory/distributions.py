"""Distribution summaries for AF3 Confidence Analysis Pipeline."""

import pandas as pd
from typing import Dict, Any, List


def summarise_distributions(long_df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarise distributions at both sample and seed-mean levels.
    
    Uses median and IQR to avoid assumptions about distribution shape.
    
    Args:
        long_df: Long-format measurements
    
    Returns:
        DataFrame with distribution summaries per condition/metric
    """
    results = []
    
    for (cond_id, metric_id, scope_type, scope_id), group in long_df.groupby(
        ["condition_id", "metric_id", "scope_type", "scope_id"]
    ):
        values = group["value"].dropna()
        
        if len(values) == 0:
            continue
        
        # Sample-level statistics
        sample_median = values.median()
        sample_iqr = values.quantile(0.75) - values.quantile(0.25)
        sample_sd = values.std(ddof=1) if len(values) > 1 else None
        
        # Seed-level aggregation first
        seed_means = group.groupby("seed")["value"].mean()
        seed_median = seed_means.median()
        seed_iqr = seed_means.quantile(0.75) - seed_means.quantile(0.25)
        seed_count = len(seed_means)
        
        results.append({
            "condition_id": cond_id,
            "metric_id": metric_id,
            "scope_type": scope_type,
            "scope_id": scope_id,
            "n_samples": len(values),
            "n_seeds": seed_count,
            "sample_median": sample_median,
            "sample_iqr": sample_iqr,
            "sample_sd": sample_sd,
            "seed_median": seed_median,
            "seed_iqr": seed_iqr,
        })
    
    return pd.DataFrame(results) if results else pd.DataFrame(columns=[
        "condition_id", "metric_id", "scope_type", "scope_id",
        "n_samples", "n_seeds", "sample_median", "sample_iqr", "sample_sd",
        "seed_median", "seed_iqr",
    ])


__all__ = ["summarise_distributions"]
