"""
Structural metrics sub-package.

Provides pluggable metric implementations for structural comparison.
"""

from af3_analysis.structural.metrics.base import StructuralMetric
from af3_analysis.structural.metrics.registry import MetricRegistry, get_default_registry

__all__ = ["StructuralMetric", "MetricRegistry", "get_default_registry"]
