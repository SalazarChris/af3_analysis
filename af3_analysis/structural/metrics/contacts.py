"""
Contact-map metrics for structural comparison.

Computes contact-map differences using a configurable distance cutoff.
"""

from __future__ import annotations

from typing import List, Optional, Set, Tuple

import numpy as np

from af3_analysis.structural.alignment import AlignmentResult
from af3_analysis.structural.enums import MetricDirection, MetricScope
from af3_analysis.structural.metrics.base import StructuralMetric
from af3_analysis.structural.representation import NormalisedStructure


def _compute_contact_set(
    structure: NormalisedStructure,
    cutoff: float,
    chain_ids: Optional[List[str]] = None,
) -> Set[Tuple[int, int]]:
    """
    Compute set of residue-residue contacts within cutoff.

    Returns set of (res_i, res_j) pairs where i < j and
    minimum heavy-atom distance < cutoff.
    """
    # Build residue centres
    from collections import defaultdict

    residue_coords = defaultdict(list)
    for atom in structure.atoms:
        if chain_ids and atom.chain_id not in chain_ids:
            continue
        if atom.atom_type == "H":
            continue
        key = (atom.chain_id, atom.auth_seq_id or 0)
        residue_coords[key].append(np.array(atom.coords))

    centres = {}
    for key, coords in residue_coords.items():
        centres[key] = np.mean(coords, axis=0)

    keys = sorted(centres.keys())
    n = len(keys)
    if n < 2:
        return set()

    coords = np.array([centres[k] for k in keys])
    from scipy.spatial.distance import pdist
    dists = pdist(coords)
    mask = dists < cutoff

    rows, cols = np.triu_indices(n, k=1)
    return set(zip(rows[mask].tolist(), cols[mask].tolist()))


class ContactMapDifference(StructuralMetric):
    """Contact-map Jaccard distance between two structures."""

    @property
    def metric_id(self) -> str:
        return "contact_map_diff"

    @property
    def display_name(self) -> str:
        return "Contact map difference"

    @property
    def units(self) -> str:
        return "dimensionless"

    @property
    def scope(self) -> MetricScope:
        return MetricScope.GLOBAL

    @property
    def direction(self) -> MetricDirection:
        return MetricDirection.LOWER_BETTER

    def requires(self) -> List[str]:
        return ["protein"]

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        *,
        contact_cutoff: float = 8.0,
        **kwargs,
    ) -> Optional[float]:
        """Compute Jaccard distance between contact maps."""
        contacts_a = _compute_contact_set(structure_a, contact_cutoff)
        contacts_b = _compute_contact_set(structure_b, contact_cutoff)

        if not contacts_a and not contacts_b:
            return None

        if not contacts_a or not contacts_b:
            return 1.0

        intersection = contacts_a & contacts_b
        union = contacts_a | contacts_b

        if not union:
            return None

        jaccard_similarity = len(intersection) / len(union)
        return 1.0 - jaccard_similarity


class InterfaceContactCount(StructuralMetric):
    """Change in inter-chain contact count."""

    @property
    def metric_id(self) -> str:
        return "interface_contacts"

    @property
    def display_name(self) -> str:
        return "Interface contact count"

    @property
    def units(self) -> str:
        return "count"

    @property
    def scope(self) -> MetricScope:
        return MetricScope.INTERFACE

    @property
    def direction(self) -> MetricDirection:
        return MetricDirection.LOWER_BETTER

    def requires(self) -> List[str]:
        return ["protein", "protein"]

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        *,
        contact_cutoff: float = 8.0,
        **kwargs,
    ) -> Optional[float]:
        """Compute difference in inter-chain contact count."""
        from collections import defaultdict

        def _inter_chain_contacts(struct, cutoff):
            residue_coords = defaultdict(list)
            for atom in struct.atoms:
                if atom.atom_type == "H":
                    continue
                key = (atom.chain_id, atom.auth_seq_id or 0)
                residue_coords[key].append(np.array(atom.coords))

            centres = {k: np.mean(v, axis=0) for k, v in residue_coords.items()}
            keys = sorted(centres.keys())
            n = len(keys)
            if n < 2:
                return 0

            coords = np.array([centres[k] for k in keys])
            chain_arr = np.array([k[0] for k in keys])

            from scipy.spatial.distance import pdist
            dists = pdist(coords)
            rows, cols = np.triu_indices(n, k=1)
            inter_chain = chain_arr[rows] != chain_arr[cols]
            return int(np.sum((dists < cutoff) & inter_chain))

        count_a = _inter_chain_contacts(structure_a, contact_cutoff)
        count_b = _inter_chain_contacts(structure_b, contact_cutoff)

        return float(count_a - count_b)


class InterfaceMinDistance(StructuralMetric):
    """Minimum inter-chain heavy-atom distance."""

    @property
    def metric_id(self) -> str:
        return "interface_distance"

    @property
    def display_name(self) -> str:
        return "Minimum interface distance"

    @property
    def units(self) -> str:
        return "Å"

    @property
    def scope(self) -> MetricScope:
        return MetricScope.INTERFACE

    @property
    def direction(self) -> MetricDirection:
        return MetricDirection.LOWER_BETTER

    def requires(self) -> List[str]:
        return []  # works with any multi-chain structure

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        **kwargs,
    ) -> Optional[float]:
        """Minimum inter-chain distance in structure_b (reference)."""
        from scipy.spatial.distance import cdist

        chains = structure_b.chain_ids
        if len(chains) < 2:
            return None

        min_dist = float("inf")
        for i, c1 in enumerate(chains):
            coords_c1 = np.array([
                a.coords for a in structure_b.get_chain(c1) if a.atom_type != "H"
            ])
            if len(coords_c1) == 0:
                continue
            for c2 in chains[i + 1:]:
                coords_c2 = np.array([
                    a.coords for a in structure_b.get_chain(c2) if a.atom_type != "H"
                ])
                if len(coords_c2) == 0:
                    continue
                # Vectorized pairwise distance computation
                dist_matrix = cdist(coords_c1, coords_c2)
                pair_min = float(dist_matrix.min())
                if pair_min < min_dist:
                    min_dist = pair_min

        if min_dist == float("inf"):
            return None
        return min_dist
