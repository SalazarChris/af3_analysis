"""
Metric registry for structural analysis.

Maintains a registry of available structural metrics and provides
lookup and filtering by scope, entity requirements, and metric ID.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from af3_analysis.structural.metrics.base import StructuralMetric


class MetricRegistry:
    """
    Registry of structural metrics.

    Metrics are registered as instances; the registry provides
    lookup and filtering.
    """

    def __init__(self):
        self._metrics: Dict[str, StructuralMetric] = {}

    def register(self, metric: StructuralMetric) -> None:
        """Register a metric instance."""
        self._metrics[metric.metric_id] = metric

    def get(self, metric_id: str) -> Optional[StructuralMetric]:
        """Get a metric by ID."""
        return self._metrics.get(metric_id)

    def all_metrics(self) -> List[StructuralMetric]:
        """Return all registered metrics."""
        return list(self._metrics.values())

    def metric_ids(self) -> List[str]:
        """Return all metric IDs."""
        return list(self._metrics.keys())

    def applicable_metrics(
        self,
        structure_a,
        structure_b,
    ) -> List[StructuralMetric]:
        """Return metrics applicable to the given structural pair."""
        return [m for m in self._metrics.values() if m.is_applicable(structure_a, structure_b)]

    def filter_by_ids(self, metric_ids: List[str]) -> List[StructuralMetric]:
        """Return only metrics whose IDs are in the given list."""
        return [m for m in self._metrics.values() if m.metric_id in metric_ids]


def get_default_registry() -> MetricRegistry:
    """Create and return the default metric registry with all built-in metrics."""
    from af3_analysis.structural.metrics.rmsd import (
        GlobalCaRMSD,
        GlobalBackboneRMSD,
        ChainCaRMSD,
        MaxAtomDisplacement,
        MeanAtomDisplacement,
    )
    from af3_analysis.structural.metrics.centroid import (
        CentroidDisplacement,
        ChainCentroidDisplacement,
    )
    from af3_analysis.structural.metrics.contacts import (
        ContactMapDifference,
        InterfaceContactCount,
        InterfaceMinDistance,
    )
    from af3_analysis.structural.metrics.distance_matrix import (
        PairwiseDistanceChange,
    )
    from af3_analysis.structural.metrics.interface import (
        InterfaceRMSD,
    )

    registry = MetricRegistry()

    for metric_cls in [
        GlobalCaRMSD,
        GlobalBackboneRMSD,
        ChainCaRMSD,
        MaxAtomDisplacement,
        MeanAtomDisplacement,
        CentroidDisplacement,
        ChainCentroidDisplacement,
        ContactMapDifference,
        InterfaceContactCount,
        InterfaceMinDistance,
        PairwiseDistanceChange,
        InterfaceRMSD,
    ]:
        registry.register(metric_cls())

    return registry
