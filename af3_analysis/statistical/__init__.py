"""
Statistics module for AF3 Confidence Analysis Pipeline.

Provides descriptive, resampling, comparisons, factorial models, variance, and multiplicity analysis.
"""

# Eligibility and gating
from af3_analysis.statistical.eligibility import (
    EligibilityGates,
    check_stage_eligibility,
    get_analysis_eligibility,
)

# Descriptive statistics
from af3_analysis.statistical.descriptive import (
    DescriptiveStats,
    compute_condition_summaries,
    compute_variability_summaries,
)

# Resampling
from af3_analysis.statistical.resampling import (
    Resampler,
    BootstrapResult,
    PermutationResult,
    bootstrap_confidence_interval,
    permutation_test,
)

# Comparisons
from af3_analysis.statistical.comparisons import (
    TwoConditionComparisons,
    TwoConditionResult,
    compare_two_conditions,
)

# Variance decomposition
from af3_analysis.statistical.variance import (
    VarianceDecomposer,
    VarianceComponents,
    VarianceSummary,
    decompose_variance,
    compute_variance_summary,
)

# Multiplicity correction
from af3_analysis.statistical.multiplicity import (
    HolmCorrector,
    HolmResult,
    holm_correction,
)

__all__ = [
    # Eligibility
    "EligibilityGates",
    "check_stage_eligibility",
    "get_analysis_eligibility",
    # Descriptive
    "DescriptiveStats",
    "compute_condition_summaries",
    "compute_variability_summaries",
    # Resampling
    "Resampler",
    "BootstrapResult",
    "PermutationResult",
    "bootstrap_confidence_interval",
    "permutation_test",
    # Comparisons
    "TwoConditionComparisons",
    "TwoConditionResult",
    "compare_two_conditions",
    # Variance
    "VarianceDecomposer",
    "VarianceComponents",
    "VarianceSummary",
    "decompose_variance",
    "compute_variance_summary",
    # Multiplicity
    "HolmCorrector",
    "HolmResult",
    "holm_correction",
]