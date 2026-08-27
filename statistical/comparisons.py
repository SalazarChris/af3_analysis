"""
Comparisons for AF3 Confidence Analysis Pipeline.

Implements paired sign-flip tests and unpaired label-permutation tests.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np

from af3_analysis.statistical.resampling import Resampler


@dataclass
class TwoConditionResult:
    """Result of two-condition comparison."""
    condition_a: str
    condition_b: str
    native_diff: float  # Native-unit difference (a - b)
    native_ci_lower: float
    native_ci_upper: float
    native_ci_width: float
    effect_size: Optional[float]  # Standardized effect (Hedges' g or d_z)
    effect_ci_lower: Optional[float]
    effect_ci_upper: Optional[float]
    p_value: float
    n_seeds_a: int
    n_seeds_b: int
    paired: bool
    test_type: str


@dataclass
class MultiConditionResult:
    """Result of multi-condition comparison."""
    metric_id: str
    scope: str
    omnibus_p_value: float
    pairwise_results: List[Dict[str, Any]]
    holm_adjusted_p_values: List[float]


class TwoConditionComparisons:
    """
    Compare two conditions using seed-level data.
    """
    
    def __init__(self, seed: Optional[int] = None):
        self.resampler = Resampler(seed=seed)
    
    def compare_two_conditions(
        self,
        cond_a: str,
        cond_b: str,
        values_a: np.ndarray,
        values_b: np.ndarray,
        paired: bool = False,
        confidence_level: float = 0.95,
        n_bootstrap: int = 2000,
        n_permutation: int = 2000,
    ) -> TwoConditionResult:
        """
        Compare two conditions.
        
        Args:
            cond_a: First condition name
            cond_b: Second condition name
            values_a: Seed-level values for condition A
            values_b: Seed-level values for condition B
            paired: Whether conditions are paired by design
            confidence_level: Confidence level for CI
            n_bootstrap: Number of bootstrap iterations
            n_permutation: Number of permutation iterations
        
        Returns:
            TwoConditionResult
        """
        # Compute native difference
        native_diff = float(np.mean(values_a) - np.mean(values_b))
        
        # Bootstrap CI for the difference
        bootstrap_res = self.resampler.bootstrap_mean(
            values_a - values_b if paired else values_a,
            confidence_level=confidence_level,
            n_bootstrap=n_bootstrap,
        )
        
        # Permutation test
        if paired:
            # Paired sign-flip test
            perm_res = self._paired_sign_flip_test(values_a, values_b, n_permutation)
        else:
            # Unpaired label permutation test
            perm_res = self.resampler.permutation_test(
                values_a, values_b, n_permutations=n_permutation
            )
        
        # Compute effect size
        effect_size = self._compute_effect_size(values_a, values_b, paired)
        
        return TwoConditionResult(
            condition_a=cond_a,
            condition_b=cond_b,
            native_diff=native_diff,
            native_ci_lower=bootstrap_res.ci_lower,
            native_ci_upper=bootstrap_res.ci_upper,
            native_ci_width=bootstrap_res.ci_width,
            effect_size=effect_size,
            effect_ci_lower=None,  # Can be computed via bootstrap if needed
            effect_ci_upper=None,
            p_value=perm_res.p_value,
            n_seeds_a=len(values_a),
            n_seeds_b=len(values_b),
            paired=paired,
            test_type="paired_sign_flip" if paired else "unpaired_permutation",
        )
    
    def _paired_sign_flip_test(
        self,
        values_a: np.ndarray,
        values_b: np.ndarray,
        n_permutations: int,
    ) -> Any:
        """Paired sign-flip permutation test."""
        # Compute paired differences
        diffs = values_a - values_b
        
        # Permutation distribution under null (no difference)
        perm_stats = []
        for _ in range(n_permutations):
            # Randomly flip signs
            signs = np.random.choice([-1, 1], size=len(diffs))
            perm_stat = np.mean(diffs * signs)
            perm_stats.append(perm_stat)
        
        perm_stats = np.array(perm_stats)
        
        # Two-sided p-value
        observed = np.mean(diffs)
        p_value = float(np.mean(np.abs(perm_stats) >= np.abs(observed)))
        
        return type("PermutationResult", (), {
            "observed_stat": observed,
            "p_value": p_value,
            "n_permutations": n_permutations,
        })()
    
    def _compute_effect_size(
        self,
        values_a: np.ndarray,
        values_b: np.ndarray,
        paired: bool,
    ) -> Optional[float]:
        """Compute standardized effect size."""
        if paired:
            # Paired standardized mean difference (d_z)
            diffs = values_a - values_b
            if len(diffs) < 2:
                return None
            return float(np.mean(diffs) / np.std(diffs, ddof=1))
        else:
            # Hedges' g (small-sample corrected)
            n1, n2 = len(values_a), len(values_b)
            if n1 < 2 or n2 < 2:
                return None
            
            pooled_sd = np.sqrt(
                ((n1 - 1) * np.var(values_a, ddof=1) + (n2 - 1) * np.var(values_b, ddof=1)) 
                / (n1 + n2 - 2)
            )
            
            if pooled_sd == 0:
                return None
            
            mean_diff = np.mean(values_a) - np.mean(values_b)
            
            # Hedges' g correction
            g = mean_diff / pooled_sd
            j = 1 - (3 / (4 * (n1 + n2) - 9))
            return float(g * j)


def compare_two_conditions(
    cond_a: str,
    cond_b: str,
    values_a: np.ndarray,
    values_b: np.ndarray,
    paired: bool = False,
    confidence_level: float = 0.95,
    n_bootstrap: int = 2000,
    n_permutation: int = 2000,
    seed: Optional[int] = None,
) -> TwoConditionResult:
    """Convenience function for two-condition comparison."""
    comparer = TwoConditionComparisons(seed=seed)
    return comparer.compare_two_conditions(
        cond_a, cond_b, values_a, values_b, paired, confidence_level, n_bootstrap, n_permutation
    )


__all__ = [
    "TwoConditionResult",
    "MultiConditionResult",
    "TwoConditionComparisons",
    "compare_two_conditions",
]
