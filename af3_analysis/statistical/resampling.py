"""
Resampling for AF3 Confidence Analysis Pipeline.

Implements whole-seed bootstrap intervals using configured coverage/iterations.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd
import numpy as np


@dataclass
class BootstrapResult:
    """Result of bootstrap resampling."""
    estimate: float
    ci_lower: float
    ci_upper: float
    ci_width: float
    ci_coverage: float
    n_resamples: int
    se: float


@dataclass
class PermutationResult:
    """Result of permutation test."""
    observed_stat: float
    p_value: float
    n_permutations: int
    alternative: str
    stat_distribution: np.ndarray


class Resampler:
    """
    Perform whole-seed bootstrap resampling.
    
    Resamples whole seed-level observations, not individual samples.
    """
    
    def __init__(self, n_bootstrap: int = 2000, n_permutation: int = 2000, seed: Optional[int] = None):
        self.n_bootstrap = n_bootstrap
        self.n_permutation = n_permutation
        self.rng = np.random.RandomState(seed)
    
    def bootstrap_mean(
        self,
        values: np.ndarray,
        confidence_level: float = 0.95,
        n_bootstrap: Optional[int] = None,
    ) -> BootstrapResult:
        """
        Bootstrap confidence interval for mean.
        
        Args:
            values: Sample values
            confidence_level: Confidence level (0-1)
            n_bootstrap: Number of bootstrap iterations (uses instance default if None)
        
        Returns:
            BootstrapResult with estimate and CI
        """
        if len(values) == 0:
            return BootstrapResult(
                estimate=np.nan,
                ci_lower=np.nan,
                ci_upper=np.nan,
                ci_width=np.nan,
                ci_coverage=confidence_level,
                n_resamples=0,
                se=np.nan,
            )
        
        n_bootstrap = n_bootstrap or self.n_bootstrap
        
        # Bootstrap distribution of mean
        bootstrap_means = []
        for _ in range(n_bootstrap):
            sample = self.rng.choice(values, size=len(values), replace=True)
            bootstrap_means.append(sample.mean())
        
        bootstrap_means = np.array(bootstrap_means)
        
        # Compute CI
        alpha = 1 - confidence_level
        lower = float(np.percentile(bootstrap_means, alpha / 2 * 100))
        upper = float(np.percentile(bootstrap_means, (1 - alpha / 2) * 100))
        
        return BootstrapResult(
            estimate=float(np.mean(values)),
            ci_lower=lower,
            ci_upper=upper,
            ci_width=upper - lower,
            ci_coverage=confidence_level,
            n_resamples=n_bootstrap,
            se=float(np.std(bootstrap_means, ddof=1)),
        )
    
    def permutation_test(
        self,
        values_a: np.ndarray,
        values_b: np.ndarray,
        stat_func: str = "mean_diff",
        alternative: str = "two-sided",
        n_permutations: Optional[int] = None,
    ) -> PermutationResult:
        """
        Permutation test for difference in means.
        
        Args:
            values_a: First group values
            values_b: Second group values
            stat_func: Statistic to compute ("mean_diff", "mean_ratio")
            alternative: Alternative hypothesis ("two-sided", "less", "greater")
            n_permutations: Number of permutations (uses instance default if None)
        
        Returns:
            PermutationResult with p-value
        """
        n_permutations = n_permutations or self.n_permutation
        
        # Compute observed statistic
        if stat_func == "mean_diff":
            observed_stat = float(np.mean(values_a) - np.mean(values_b))
        else:
            raise ValueError(f"Unknown stat_func: {stat_func}")
        
        # Combine values for permutation
        combined = np.concatenate([values_a, values_b])
        n_a = len(values_a)
        
        # Permutation distribution
        perm_stats = []
        for _ in range(n_permutations):
            self.rng.shuffle(combined)
            perm_a = combined[:n_a]
            perm_b = combined[n_a:]
            
            if stat_func == "mean_diff":
                perm_stat = perm_a.mean() - perm_b.mean()
            perm_stats.append(perm_stat)
        
        perm_stats = np.array(perm_stats)
        
        # Compute p-value
        if alternative == "two-sided":
            p_value = np.mean(np.abs(perm_stats) >= np.abs(observed_stat))
        elif alternative == "greater":
            p_value = np.mean(perm_stats >= observed_stat)
        elif alternative == "less":
            p_value = np.mean(perm_stats <= observed_stat)
        else:
            raise ValueError(f"Unknown alternative: {alternative}")
        
        return PermutationResult(
            observed_stat=observed_stat,
            p_value=float(p_value),
            n_permutations=n_permutations,
            alternative=alternative,
            stat_distribution=perm_stats,
        )


def bootstrap_confidence_interval(
    values: np.ndarray,
    confidence_level: float = 0.95,
    n_bootstrap: int = 2000,
    seed: Optional[int] = None,
) -> BootstrapResult:
    """Convenience function for bootstrap CI."""
    resampler = Resampler(n_bootstrap=n_bootstrap, seed=seed)
    return resampler.bootstrap_mean(values, confidence_level)


def permutation_test(
    values_a: np.ndarray,
    values_b: np.ndarray,
    stat_func: str = "mean_diff",
    alternative: str = "two-sided",
    n_permutations: int = 2000,
    seed: Optional[int] = None,
) -> PermutationResult:
    """Convenience function for permutation test."""
    resampler = Resampler(n_permutation=n_permutations, seed=seed)
    return resampler.permutation_test(values_a, values_b, stat_func, alternative, n_permutations)


__all__ = [
    "BootstrapResult",
    "PermutationResult",
    "Resampler",
    "bootstrap_confidence_interval",
    "permutation_test",
]
