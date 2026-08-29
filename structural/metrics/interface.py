"""
Interface metrics for multi-chain structural comparison.

Computes interface-specific metrics for multi-chain complexes.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from af3_analysis.structural.alignment import AlignmentResult
from af3_analysis.structural.enums import MetricDirection, MetricScope
from af3_analysis.structural.metrics.base import StructuralMetric
from af3_analysis.structural.representation import NormalisedStructure


class InterfaceRMSD(StructuralMetric):
    """RMSD of interface residues between two structures."""

    @property
    def metric_id(self) -> str:
        return "interface_rmsd"

    @property
    def display_name(self) -> str:
        return "Interface RMSD"

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
        *,
        contact_cutoff: float = 8.0,
        **kwargs,
    ) -> Optional[float]:
        """Compute RMSD of interface C-alpha atoms."""
        from collections import defaultdict

        def _get_interface_residues(struct, cutoff):
            residue_coords = defaultdict(list)
            for atom in struct.atoms:
                if atom.atom_name != "CA":
                    continue
                key = (atom.chain_id, atom.auth_seq_id or 0)
                residue_coords[key] = np.array(atom.coords)

            chains = list({k[0] for k in residue_coords.keys()})
            if len(chains) < 2:
                return set()

            interface_residues = set()
            chain_groups = defaultdict(list)
            for key, coords in residue_coords.items():
                chain_groups[key[0]].append((key, coords))

            for i, c1 in enumerate(chains):
                for c2 in chains[i + 1:]:
                    for k1, coords1 in chain_groups.get(c1, []):
                        for k2, coords2 in chain_groups.get(c2, []):
                            dist = np.sqrt(np.sum((coords1 - coords2) ** 2))
                            if dist < cutoff:
                                interface_residues.add(k1)
                                interface_residues.add(k2)
            return interface_residues

        iface_a = _get_interface_residues(structure_a, contact_cutoff)
        iface_b = _get_interface_residues(structure_b, contact_cutoff)

        # Use intersection of interface residues
        common = iface_a & iface_b
        if len(common) < 3:
            return None

        # Get CA coords for common interface residues
        def _get_coords(struct, residues):
            coords = []
            for atom in struct.atoms:
                if atom.atom_name == "CA":
                    key = (atom.chain_id, atom.auth_seq_id or 0)
                    if key in residues:
                        coords.append(atom.coords)
            return np.array(coords)

        coords_a = _get_coords(structure_a, common)
        coords_b = _get_coords(structure_b, common)

        if len(coords_a) < 3:
            return None

        # Apply alignment
        if alignment.rotation is not None:
            from af3_analysis.structural.alignment import apply_alignment
            coords_a = apply_alignment(coords_a, alignment)

        diff = coords_a - coords_b
        return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
