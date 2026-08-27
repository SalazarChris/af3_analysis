"""
Registry module for AF3 Confidence Analysis Pipeline.

Provides metric registry loading and resolution functionality.
"""

from af3_analysis.registry.metric_registry import (
    MetricRegistry,
    load_registry,
    resolve_metric_availability,
)

from af3_analysis.registry.resolution import (
    MetricResolution,
    ResolutionResult,
    resolve_metrics,
)

__all__ = [
    "MetricRegistry",
    "load_registry",
    "resolve_metric_availability",
    "MetricResolution",
    "ResolutionResult",
    "resolve_metrics",
]
