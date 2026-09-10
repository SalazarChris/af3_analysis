"""Exploratory analysis for AF3 Confidence Analysis Pipeline."""

from af3_analysis.exploratory.inventory import DesignInventory, is_fully_crossed, inventory_markdown
from af3_analysis.exploratory.variance import seed_means, variance_components, reference_spread, min_detectable_difference, degenerate_metrics
from af3_analysis.exploratory.distributions import summarise_distributions
from af3_analysis.exploratory.factors import factor_marginals, monotonicity_flag

__all__ = [
    "DesignInventory",
    "is_fully_crossed",
    "inventory_markdown",
    "seed_means",
    "variance_components",
    "reference_spread",
    "min_detectable_difference",
    "degenerate_metrics",
    "summarise_distributions",
    "factor_marginals",
    "monotonicity_flag",
]
