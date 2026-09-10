"""
Metric registry for AF3 Confidence Analysis Pipeline.

Loads immutable metric registry source data and creates resolved availability tables.
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional
import pandas as pd

from af3_analysis.schemas.records import RegistryEntry
from af3_analysis.schemas.enums import (
    MetricResolutionStatus,
    MetricPortfolioCategory,
    MetricDirection,
    MetricScope,
)


class MetricRegistry:
    """
    Immutable metric registry.
    
    Contains definitions for all AF3 confidence metrics with metadata
    about their source, scope, direction, and portfolio category.
    """
    
    def __init__(self, entries: List[RegistryEntry]):
        self._entries = {e.metric_id: e for e in entries}
        self._by_category: Dict[MetricPortfolioCategory, List[str]] = {}
        self._by_scope: Dict[MetricScope, List[str]] = {}
        
        for entry in entries:
            # Index by portfolio category
            if entry.portfolio_category not in self._by_category:
                self._by_category[entry.portfolio_category] = []
            self._by_category[entry.portfolio_category].append(entry.metric_id)
            
            # Index by scope
            if entry.scope not in self._by_scope:
                self._by_scope[entry.scope] = []
            self._by_scope[entry.scope].append(entry.metric_id)
    
    @property
    def entries(self) -> List[RegistryEntry]:
        """Return all registry entries."""
        return list(self._entries.values())
    
    @property
    def metric_ids(self) -> List[str]:
        """Return all metric IDs."""
        return list(self._entries.keys())
    
    def get_entry(self, metric_id: str) -> Optional[RegistryEntry]:
        """Get a registry entry by metric ID."""
        return self._entries.get(metric_id)
    
    def get_by_category(self, category: MetricPortfolioCategory) -> List[str]:
        """Get metric IDs in a portfolio category."""
        return self._by_category.get(category, [])
    
    def get_by_scope(self, scope: MetricScope) -> List[str]:
        """Get metric IDs for a scope."""
        return self._by_scope.get(scope, [])
    
    def get_all_categories(self) -> List[MetricPortfolioCategory]:
        """Get all portfolio categories in the registry."""
        return list(self._by_category.keys())
    
    def get_all_scopes(self) -> List[MetricScope]:
        """Get all scopes in the registry."""
        return list(self._by_scope.keys())
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert registry to DataFrame for inspection."""
        data = []
        for entry in self._entries.values():
            data.append({
                "metric_id": entry.metric_id,
                "display_name": entry.display_name,
                "source": entry.source,
                "scope": entry.scope.value,
                "units": entry.units,
                "direction": entry.direction.value,
                "portfolio_category": entry.portfolio_category.value,
                "precision": entry.precision,
                "expected_range": entry.expected_range,
                "caveats": entry.caveats,
                "resolution_status": entry.resolution_status,
                "resolution_reason": entry.resolution_reason,
            })
        return pd.DataFrame(data)
    
    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> "MetricRegistry":
        """Create a registry from a DataFrame."""
        entries = []
        for _, row in df.iterrows():
            entry = RegistryEntry(
                metric_id=row["metric_id"],
                display_name=row["display_name"],
                source=row["source"],
                scope=MetricScope(row["scope"]),
                units=row["units"],
                direction=MetricDirection(row["direction"]),
                portfolio_category=MetricPortfolioCategory(row["portfolio_category"]),
                precision=row.get("precision"),
                expected_range=row.get("expected_range"),
                caveats=row.get("caveats"),
                resolution_status=row.get("resolution_status", "unresolved"),
                resolution_reason=row.get("resolution_reason"),
            )
            entries.append(entry)
        return cls(entries)
    
    @classmethod
    def load(cls, source: Path) -> "MetricRegistry":
        """Load registry from a JSON file."""
        with open(source, "r") as f:
            data = json.load(f)
        
        entries = []
        for item in data["metrics"]:
            entry = RegistryEntry(
                metric_id=item["metric_id"],
                display_name=item["display_name"],
                source=item["source"],
                scope=MetricScope(item["scope"]),
                units=item["units"],
                direction=MetricDirection(item["direction"]),
                portfolio_category=MetricPortfolioCategory(item["portfolio_category"]),
                precision=item.get("precision"),
                expected_range=item.get("expected_range"),
                caveats=item.get("caveats"),
                resolution_status=item.get("resolution_status", "unresolved"),
                resolution_reason=item.get("resolution_reason"),
            )
            entries.append(entry)
        
        return cls(entries)


def load_registry(source: Path) -> MetricRegistry:
    """
    Load a metric registry from a JSON file.
    
    Expected format:
    {
      "metrics": [
        {
          "metric_id": "ranking_score",
          "display_name": "Ranking Score",
          ...
        }
      ]
    }
    """
    return MetricRegistry.load(source)


def resolve_metric_availability(
    registry: MetricRegistry,
    available_columns: List[str],
    condition_composition: Optional[Dict[str, Any]] = None,
) -> Dict[str, MetricResolutionStatus]:
    """
    Resolve metric availability based on available columns and condition composition.
    
    Returns a dict mapping metric_id -> resolution status.
    """
    resolutions = {}
    
    for metric_id in registry.metric_ids:
        entry = registry.get_entry(metric_id)
        
        # Check if metric is in available columns
        if metric_id not in available_columns:
            resolutions[metric_id] = MetricResolutionStatus.UNAVAILABLE
            continue
        
        # Check composition constraints
        if condition_composition:
            scope = entry.scope
            if scope == MetricScope.CHAIN_PAIR and not condition_composition.get("is_complex", False):
                resolutions[metric_id] = MetricResolutionStatus.UNDEFINED_BY_COMPOSITION
                continue
        
        # Default to available
        resolutions[metric_id] = MetricResolutionStatus.RESOLVED_AVAILABLE
    
    return resolutions


__all__ = [
    "MetricRegistry",
    "load_registry",
    "resolve_metric_availability",
]
