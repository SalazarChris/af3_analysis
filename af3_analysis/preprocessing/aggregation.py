"""
Aggregation for AF3 Confidence Analysis Pipeline.

Builds analysis_sample and analysis_seed tables with seed-level aggregation.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np

from af3_analysis.schemas.enums import Definedness


@dataclass
class AggregationResult:
    """Result of aggregation."""
    analysis_sample: pd.DataFrame
    analysis_seed: pd.DataFrame
    condition_summaries: pd.DataFrame
    variability_summaries: pd.DataFrame


class Aggregator:
    """
    Aggregate sample-level data to seed-level.
    
    Samples are nested within seeds. Default inferential analyses use
    one seed-level observation per condition × metric × scope.
    """
    
    def __init__(self):
        self._analysis_samples: List[Dict[str, Any]] = []
        self._analysis_seeds: List[Dict[str, Any]] = []
    
    def aggregate(
        self,
        measurements_df: pd.DataFrame,
        replicates_df: pd.DataFrame,
    ) -> AggregationResult:
        """
        Aggregate sample-level measurements to seed-level.
        
        Args:
            measurements_df: Measurements table (long format)
            replicates_df: Replicates table
        
        Returns:
            AggregationResult with analysis_sample, analysis_seed, and summaries
        """
        # Build analysis_sample (one row per sample)
        self._build_analysis_samples(measurements_df, replicates_df)
        
        # Build analysis_seed (one row per seed)
        self._build_analysis_seeds()
        
        # Create condition summaries
        condition_summaries = self._build_condition_summaries()
        
        # Create variability summaries
        variability_summaries = self._build_variability_summaries()
        
        return AggregationResult(
            analysis_sample=pd.DataFrame(self._analysis_samples),
            analysis_seed=pd.DataFrame(self._analysis_seeds),
            condition_summaries=condition_summaries,
            variability_summaries=variability_summaries,
        )
    
    def _build_analysis_samples(
        self,
        measurements_df: pd.DataFrame,
        replicates_df: pd.DataFrame,
    ) -> None:
        """Build analysis_sample table."""
        # Merge measurements with replicates to get condition_id and seed
        merged = measurements_df.merge(
            replicates_df[["prediction_id", "condition_id", "seed", "sample"]],
            on="prediction_id",
            how="inner",
        )
        
        for _, row in merged.iterrows():
            self._analysis_samples.append({
                "condition_id": row["condition_id"],
                "seed": row["seed"],
                "sample": row["sample"],
                "metric_id": row["metric_id"],
                "scope_type": row.get("scope_type", "global"),
                "scope_id": row.get("scope_id", ""),
                "value": row["value"],
                "definedness": row["definedness"],
            })
    
    def _build_analysis_seeds(self) -> None:
        """Build analysis_seed table from analysis_samples."""
        if not self._analysis_samples:
            return
        
        samples_df = pd.DataFrame(self._analysis_samples)
        
        # Group by condition_id, seed, metric_id, scope_type, scope_id
        grouped = samples_df.groupby(
            ["condition_id", "seed", "metric_id", "scope_type", "scope_id"]
        )
        
        for key, group in grouped:
            cond_id, seed, metric_id, scope_type, scope_id = key
            values = group["value"].dropna()
            
            if len(values) == 0:
                # No valid values - use first record's definedness
                definedness = group.iloc[0]["definedness"]
                self._analysis_seeds.append({
                    "condition_id": cond_id,
                    "seed": seed,
                    "metric_id": metric_id,
                    "scope_type": scope_type,
                    "scope_id": scope_id,
                    "value": None,
                    "definedness": definedness,
                    "n_samples": 0,
                    "sample_sd": None,
                    "sample_min": None,
                    "sample_max": None,
                })
            elif len(values) == 1:
                # Single sample - use as-is
                self._analysis_seeds.append({
                    "condition_id": cond_id,
                    "seed": seed,
                    "metric_id": metric_id,
                    "scope_type": scope_type,
                    "scope_id": scope_id,
                    "value": values.iloc[0],
                    "definedness": group.iloc[0]["definedness"],
                    "n_samples": 1,
                    "sample_sd": None,
                    "sample_min": values.iloc[0],
                    "sample_max": values.iloc[0],
                })
            else:
                # Multiple samples - compute summary stats
                self._analysis_seeds.append({
                    "condition_id": cond_id,
                    "seed": seed,
                    "metric_id": metric_id,
                    "scope_type": scope_type,
                    "scope_id": scope_id,
                    "value": values.mean(),
                    "definedness": group.iloc[0]["definedness"],
                    "n_samples": len(values),
                    "sample_sd": values.std(ddof=1),
                    "sample_min": values.min(),
                    "sample_max": values.max(),
                })
    
    def _build_condition_summaries(self) -> pd.DataFrame:
        """Build condition-level descriptive summaries."""
        if not self._analysis_seeds:
            return pd.DataFrame()
        
        seeds_df = pd.DataFrame(self._analysis_seeds)
        
        # Group by condition_id, metric_id, scope_type, scope_id
        grouped = seeds_df.groupby(["condition_id", "metric_id", "scope_type", "scope_id"])
        
        summaries = []
        for (cond_id, metric_id, scope_type, scope_id), group in grouped:
            values = group["value"].dropna()
            
            if len(values) == 0:
                summaries.append({
                    "condition_id": cond_id,
                    "metric_id": metric_id,
                    "scope_type": scope_type,
                    "scope_id": scope_id,
                    "n_seeds": 0,
                    "mean": None,
                    "sd": None,
                    "median": None,
                    "min": None,
                    "max": None,
                    "seeds": [],
                })
            else:
                summaries.append({
                    "condition_id": cond_id,
                    "metric_id": metric_id,
                    "scope_type": scope_type,
                    "scope_id": scope_id,
                    "n_seeds": len(values),
                    "mean": values.mean(),
                    "sd": values.std(ddof=1),
                    "median": values.median(),
                    "min": values.min(),
                    "max": values.max(),
                    "seeds": group["seed"].unique().tolist(),
                })
        
        return pd.DataFrame(summaries)
    
    def _build_variability_summaries(self) -> pd.DataFrame:
        """Build variability decomposition summaries."""
        if not self._analysis_seeds:
            return pd.DataFrame()
        
        seeds_df = pd.DataFrame(self._analysis_seeds)
        
        summaries = []
        for key, group in seeds_df.groupby(["metric_id", "scope_type", "scope_id"]):
            metric_id, scope_type, scope_id = key
            # Between-seed variance
            between_seed = group.groupby("condition_id")["value"].mean()
            between_seed_var = between_seed.var(ddof=1) if len(between_seed) > 1 else 0
            
            # Within-seed variance (sample_sd from analysis_seed)
            within_seed_var = group["sample_sd"].mean() ** 2 if group["sample_sd"].notna().any() else 0
            
            summaries.append({
                "metric_id": metric_id,
                "scope_type": scope_type,
                "scope_id": scope_id,
                "n_conditions": group["condition_id"].nunique(),
                "n_seeds": len(group),
                "between_seed_var": between_seed_var,
                "within_seed_var": within_seed_var,
                "total_var": between_seed_var + within_seed_var,
            })
        
        return pd.DataFrame(summaries)


def aggregate_to_seed_level(
    measurements_df: pd.DataFrame,
    replicates_df: pd.DataFrame,
) -> AggregationResult:
    """
    Convenience function to aggregate to seed-level.
    
    Args:
        measurements_df: Measurements table
        replicates_df: Replicates table
    
    Returns:
        AggregationResult with analysis_sample and analysis_seed
    """
    aggregator = Aggregator()
    return aggregator.aggregate(measurements_df, replicates_df)


__all__ = [
    "AggregationResult",
    "Aggregator",
    "aggregate_to_seed_level",
]
