"""
Pairwise distance metrics for structural comparison.

Computes changes in pairwise residue-residue distances.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from af3_analysis.structural.alignment import AlignmentResult, apply_alignment
from af3_analysis.structural.enums import MetricDirection, MetricScope
from af3_analysis.structural.metrics.base import StructuralMetric
from af3_analysis.structural.representation import NormalisedStructure


class PairwiseDistanceChange(StructuralMetric):
    """Mean change in pairwise C-alpha distances after alignment."""

    @property
    def metric_id(self) -> str:
        return "pairwise_distance_change"

    @property
    def display_name(self) -> str:
        return "Pairwise distance change"

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
        """Mean absolute change in pairwise C-alpha distances."""
        if not alignment.paired_atom_ids:
            return None

        atom_map_a = {a.atom_id: a for a in structure_a.atoms}
        atom_map_b = {b.atom_id: b for b in structure_b.atoms}

        coords_a = []
        coords_b = []
        for aid, bid in alignment.paired_atom_ids:
            if aid in atom_map_a and bid in atom_map_b:
                if atom_map_a[aid].atom_name == "CA" and atom_map_b[bid].atom_name == "CA":
                    coords_a.append(atom_map_a[aid].coords)
                    coords_b.append(atom_map_b[bid].coords)

        if len(coords_a) < 3:
            return None

        coords_a = np.array(coords_a)
        coords_b = np.array(coords_b)

        if alignment.rotation is not None:
            coords_a = apply_alignment(coords_a, alignment)

        # Compute pairwise distance matrices
        n = len(coords_a)
        if n > 200:
            # For large structures, subsample to keep computation tractable
            indices = np.linspace(0, n - 1, 200, dtype=int)
            coords_a = coords_a[indices]
            coords_b = coords_b[indices]
            n = 200

        dist_a = np.sqrt(np.sum((coords_a[:, None, :] - coords_a[None, :, :]) ** 2, axis=2))
        dist_b = np.sqrt(np.sum((coords_b[:, None, :] - coords_b[None, :, :]) ** 2, axis=2))

        # Upper triangle only
        triu_indices = np.triu_indices(n, k=1)
        diff = np.abs(dist_a[triu_indices] - dist_b[triu_indices])

        return float(np.mean(diff))
