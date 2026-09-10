"""
Centroid displacement metrics for structural comparison.

Computes geometric centre displacement at global, per-chain, and
per-entity levels.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from af3_analysis.structural.alignment import AlignmentResult, apply_alignment
from af3_analysis.structural.enums import MetricDirection, MetricScope
from af3_analysis.structural.metrics.base import StructuralMetric
from af3_analysis.structural.representation import NormalisedStructure


class CentroidDisplacement(StructuralMetric):
    """C-alpha centroid displacement between two structures."""

    @property
    def metric_id(self) -> str:
        return "centroid_displacement"

    @property
    def display_name(self) -> str:
        return "Centroid displacement"

    @property
    def units(self) -> str:
        return "Å"

    @property
    def scope(self) -> MetricScope:
        return MetricScope.GLOBAL

    @property
    def direction(self) -> MetricDirection:
        return MetricDirection.LOWER_BETTER

    def requires(self) -> List[str]:
        return []

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        **kwargs,
    ) -> Optional[float]:
        ca_a = structure_a.get_ca_atoms()
        ca_b = structure_b.get_ca_atoms()

        if not ca_a or not ca_b:
            return None

        centroid_a = np.mean([a.coords for a in ca_a], axis=0)
        centroid_b = np.mean([b.coords for b in ca_b], axis=0)

        # Apply alignment to centroid_a if available
        if alignment.rotation is not None and alignment.translation is not None:
            centroid_a = alignment.rotation @ centroid_a + alignment.translation

        return float(np.sqrt(np.sum((centroid_a - centroid_b) ** 2)))


class ChainCentroidDisplacement(StructuralMetric):
    """Per-chain centroid displacement."""

    @property
    def metric_id(self) -> str:
        return "centroid_displacement_chain"

    @property
    def display_name(self) -> str:
        return "Per-chain centroid displacement"

    @property
    def units(self) -> str:
        return "Å"

    @property
    def scope(self) -> MetricScope:
        return MetricScope.CHAIN

    @property
    def direction(self) -> MetricDirection:
        return MetricDirection.LOWER_BETTER

    def requires(self) -> List[str]:
        return []

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        *,
        chain_id: Optional[str] = None,
        **kwargs,
    ) -> Optional[float]:
        if chain_id is None:
            return None

        ca_a = structure_a.get_ca_atoms([chain_id])
        ca_b = structure_b.get_ca_atoms([chain_id])

        if not ca_a or not ca_b:
            return None

        centroid_a = np.mean([a.coords for a in ca_a], axis=0)
        centroid_b = np.mean([b.coords for b in ca_b], axis=0)

        if alignment.rotation is not None and alignment.translation is not None:
            centroid_a = alignment.rotation @ centroid_a + alignment.translation

        return float(np.sqrt(np.sum((centroid_a - centroid_b) ** 2)))
