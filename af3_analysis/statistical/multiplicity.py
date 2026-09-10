"""
Multiplicity correction for AF3 Confidence Analysis Pipeline.

Implements Holm correction for multiple testing control.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np


@dataclass
class HolmResult:
    """Result of Holm correction for a single test."""
    metric_id: str
    scope: str
    raw_p_value: float
    holm_adjusted_p_value: float
    rejected: bool
    family_id: str
    rank: int


class HolmCorrector:
    """
    Apply Holm-Bonferroni correction for multiplicity control.
    
    Control family-wise error rate while being more powerful than Bonferroni.
    """
    
    def __init__(self, alpha: float = 0.05):
        self.alpha = alpha
    
    def correct_holm(
        self,
        p_values: List[float],
        metric_ids: List[str],
        scopes: List[str],
        family_ids: List[str],
    ) -> List[HolmResult]:
        """
        Apply Holm correction to a set of p-values.
        
        Args:
            p_values: Raw p-values
            metric_ids: Corresponding metric IDs
            scopes: Corresponding scope strings
            family_ids: Family IDs for grouping
        
        Returns:
            List of HolmResult with adjusted p-values
        """
        if len(p_values) == 0:
            return []
        
        # Combine into DataFrame for processing
        df = pd.DataFrame({
            "metric_id": metric_ids,
            "scope": scopes,
            "family_id": family_ids,
            "raw_p_value": p_values,
        })
        
        # Process each family separately
        results = []
        for family_id, family_df in df.groupby("family_id"):
            # Sort by p-value
            sorted_df = family_df.sort_values("raw_p_value")
            
            m = len(sorted_df)
            for i, (idx, row) in enumerate(sorted_df.iterrows()):
                rank = i + 1
                
                # Holm adjustment: p * (m - rank + 1)
                adjusted_p = row["raw_p_value"] * (m - rank + 1)
                
                # Ensure monotonically increasing
                if i > 0:
                    prev_result = results[-1]
                    adjusted_p = max(adjusted_p, prev_result.holm_adjusted_p_value)
                
                # Cap at 1
                adjusted_p = min(adjusted_p, 1.0)
                
                # Decision
                rejected = adjusted_p <= self.alpha
                
                results.append(HolmResult(
                    metric_id=row["metric_id"],
                    scope=row["scope"],
                    raw_p_value=float(row["raw_p_value"]),
                    holm_adjusted_p_value=float(adjusted_p),
                    rejected=rejected,
                    family_id=family_id,
                    rank=rank,
                ))
        
        return results
    
    def compute_family_wise_error(
        self,
        results: List[HolmResult],
    ) -> float:
        """
        Compute family-wise error rate from Holm results.
        
        Args:
            results: List of HolmResult
        
        Returns:
            FWER (0 = no error, 1 = all tests rejected incorrectly)
        """
        if not results:
            return 0.0
        
        # If any test in a family is rejected, FWER = 1
        # (This is a simplified estimate)
        families = set(r.family_id for r in results)
        fwer = 0
        
        for family_id in families:
            family_results = [r for r in results if r.family_id == family_id]
            any_rejected = any(r.rejected for r in family_results)
            fwer += 1 if any_rejected else 0
        
        return fwer / len(families) if families else 0


def holm_correction(
    p_values: List[float],
    metric_ids: List[str],
    scopes: List[str],
    family_ids: List[str],
    alpha: float = 0.05,
) -> List[HolmResult]:
    """Convenience function for Holm correction."""
    corrector = HolmCorrector(alpha=alpha)
    return corrector.correct_holm(p_values, metric_ids, scopes, family_ids)


__all__ = [
    "HolmResult",
    "HolmCorrector",
    "holm_correction",
]
