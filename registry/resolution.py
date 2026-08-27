"""
Metric resolution for AF3 Confidence Analysis Pipeline.

Creates resolved availability tables from the registry and available data.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd

from af3_analysis.registry.metric_registry import MetricRegistry, resolve_metric_availability
from af3_analysis.schemas.enums import (
    MetricResolutionStatus,
    MetricPortfolioCategory,
    MetricScope,
)


@dataclass(frozen=True)
class MetricResolution:
    """
    Resolution record for a single metric.
    
    Tracks the status and reason for each metric's availability.
    """
    metric_id: str
    resolution_status: MetricResolutionStatus
    resolution_reason: str
    source: str
    scope: MetricScope
    portfolio_category: MetricPortfolioCategory


@dataclass
class ResolutionResult:
    """
    Result of metric resolution.
    
    Contains resolved metrics, exclusions, and metadata.
    """
    resolutions: List[MetricResolution]
    available_metrics: List[str]
    excluded_metrics: List[str]
    undefined_metrics: List[str]
    unavailable_metrics: List[str]
    
    @property
    def available_count(self) -> int:
        return len(self.available_metrics)
    
    @property
    def excluded_count(self) -> int:
        return len(self.excluded_metrics)
    
    @property
    def undefined_count(self) -> int:
        return len(self.undefined_metrics)
    
    @property
    def unavailable_count(self) -> int:
        return len(self.unavailable_metrics)


def resolve_metrics(
    registry: MetricRegistry,
    available_columns: List[str],
    condition_composition: Optional[Dict[str, Any]] = None,
    exclude_constants: bool = True,
) -> ResolutionResult:
    """
    Resolve metric availability and create resolution records.
    
    Args:
        registry: The metric registry to resolve against
        available_columns: List of column names available in the data
        condition_composition: Condition-specific composition info
        exclude_constants: Whether to exclude constant metrics
    
    Returns:
        ResolutionResult with available, excluded, undefined, and unavailable metrics
    """
    resolutions = []
    available_metrics = []
    excluded_metrics = []
    undefined_metrics = []
    unavailable_metrics = []
    
    for metric_id in registry.metric_ids:
        entry = registry.get_entry(metric_id)
        status = resolve_metric_availability(
            registry, available_columns, condition_composition
        ).get(metric_id, MetricResolutionStatus.UNAVAILABLE)
        
        # Determine resolution reason
        reason = _determine_resolution_reason(status, metric_id, entry, exclude_constants)
        
        resolution = MetricResolution(
            metric_id=metric_id,
            resolution_status=status,
            resolution_reason=reason,
            source=entry.source,
            scope=entry.scope,
            portfolio_category=entry.portfolio_category,
        )
        resolutions.append(resolution)
        
        # Categorize metric
        if status == MetricResolutionStatus.RESOLVED_AVAILABLE:
            available_metrics.append(metric_id)
        elif status == MetricResolutionStatus.UNDEFINED_BY_COMPOSITION:
            undefined_metrics.append(metric_id)
            excluded_metrics.append(metric_id)
        elif status == MetricResolutionStatus.UNAVAILABLE:
            unavailable_metrics.append(metric_id)
            excluded_metrics.append(metric_id)
        elif status == MetricResolutionStatus.RESOLVED_EXCLUDED:
            excluded_metrics.append(metric_id)
    
    return ResolutionResult(
        resolutions=resolutions,
        available_metrics=available_metrics,
        excluded_metrics=excluded_metrics,
        undefined_metrics=undefined_metrics,
        unavailable_metrics=unavailable_metrics,
    )


def _determine_resolution_reason(
    status: MetricResolutionStatus,
    metric_id: str,
    entry: Any,
    exclude_constants: bool,
) -> str:
    """Determine the reason for a metric's resolution status."""
    if status == MetricResolutionStatus.RESOLVED_AVAILABLE:
        return "Available in data and passes all eligibility checks"
    elif status == MetricResolutionStatus.UNDEFINED_BY_COMPOSITION:
        return f"Metric '{metric_id}' is undefined for this composition (e.g., interface metric in monomer)"
    elif status == MetricResolutionStatus.UNAVAILABLE:
        return f"Metric '{metric_id}' not found in available columns"
    elif status == MetricResolutionStatus.RESOLVED_EXCLUDED:
        if exclude_constants:
            return f"Metric '{metric_id}' excluded: constant or near-zero variance"
        return f"Metric '{metric_id}' excluded for other reasons"
    else:
        return f"Unknown resolution status: {status}"


def create_resolved_availability_table(
    resolution_result: ResolutionResult,
    condition_id: Optional[str] = None,
) -> pd.DataFrame:
    """
    Create a DataFrame of resolved metric availability.
    
    Includes all metrics with their resolution status and reason.
    """
    data = []
    for resolution in resolution_result.resolutions:
        row = {
            "metric_id": resolution.metric_id,
            "resolution_status": resolution.resolution_status.value,
            "resolution_reason": resolution.resolution_reason,
            "source": resolution.source,
            "scope": resolution.scope.value,
            "portfolio_category": resolution.portfolio_category.value,
        }
        if condition_id:
            row["condition_id"] = condition_id
        data.append(row)
    
    return pd.DataFrame(data)


__all__ = [
    "MetricResolution",
    "ResolutionResult",
    "resolve_metrics",
    "_determine_resolution_reason",
    "create_resolved_availability_table",
]
