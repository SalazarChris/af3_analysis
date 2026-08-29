"""
RMSD metrics for structural comparison.

Provides global, per-chain, and region-specific RMSD calculations.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from af3_analysis.structural.alignment import AlignmentResult, apply_alignment
from af3_analysis.structural.enums import MetricDirection, MetricScope
from af3_analysis.structural.metrics.base import StructuralMetric
from af3_analysis.structural.representation import NormalisedStructure


class GlobalCaRMSD(StructuralMetric):
    """Global C-alpha RMSD after optimal superposition."""

    @property
    def metric_id(self) -> str:
        return "rmsd_global_ca"

    @property
    def display_name(self) -> str:
        return "Global Cα RMSD"

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
        return ["protein"]

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        **kwargs,
    ) -> Optional[float]:
        if alignment.rmsd_post_alignment is not None:
            return alignment.rmsd_post_alignment
        return None


class GlobalBackboneRMSD(StructuralMetric):
    """Global backbone (N, CA, C, O) RMSD after optimal superposition."""

    @property
    def metric_id(self) -> str:
        return "rmsd_backbone"

    @property
    def display_name(self) -> str:
        return "Global backbone RMSD"

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
        return ["protein"]

    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        **kwargs,
    ) -> Optional[float]:
        if alignment.rmsd_post_alignment is not None:
            return alignment.rmsd_post_alignment
        return None


class ChainCaRMSD(StructuralMetric):
    """Per-chain C-alpha RMSD."""

    @property
    def metric_id(self) -> str:
        return "rmsd_chain"

    @property
    def display_name(self) -> str:
        return "Per-chain Cα RMSD"

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
        """Compute RMSD for a specific chain. Requires chain_id in kwargs."""
        if chain_id is None:
            return None

        from af3_analysis.structural.alignment import (
            _build_residue_index,
            _match_residues,
            _match_atoms,
            _kabsch,
        )

        atoms_a = structure_a.get_ca_atoms([chain_id])
        atoms_b = structure_b.get_ca_atoms([chain_id])

        if not atoms_a or not atoms_b:
            return None

        index_a = _build_residue_index(atoms_a)
        index_b = _build_residue_index(atoms_b)

        common = _match_residues(index_a, index_b)
        if not common:
            return None

        atom_pairs = _match_atoms(index_a, index_b, common)
        if len(atom_pairs) < 3:
            return None

        atom_map_a = {a.atom_id: a for a in atoms_a}
        atom_map_b = {b.atom_id: b for b in atoms_b}

        coords_a = np.array([atom_map_a[aid].coords for aid, _ in atom_pairs])
        coords_b = np.array([atom_map_b[bid].coords for _, bid in atom_pairs])

        try:
            _, _, rmsd = _kabsch(coords_a, coords_b)
            return rmsd
        except np.linalg.LinAlgError:
            return None


class MaxAtomDisplacement(StructuralMetric):
    """Maximum per-atom displacement after alignment."""

    @property
    def metric_id(self) -> str:
        return "max_atom_displacement"

    @property
    def display_name(self) -> str:
        return "Maximum atom displacement"

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
        if not alignment.paired_atom_ids:
            return None

        atom_map_a = {a.atom_id: a for a in structure_a.atoms}
        atom_map_b = {b.atom_id: b for b in structure_b.atoms}

        coords_a = []
        coords_b = []
        for aid, bid in alignment.paired_atom_ids:
            if aid in atom_map_a and bid in atom_map_b:
                coords_a.append(atom_map_a[aid].coords)
                coords_b.append(atom_map_b[bid].coords)

        if not coords_a:
            return None

        coords_a = np.array(coords_a)
        coords_b = np.array(coords_b)

        # Apply alignment
        if alignment.rotation is not None:
            coords_a = apply_alignment(coords_a, alignment)

        displacements = np.sqrt(np.sum((coords_a - coords_b) ** 2, axis=1))
        return float(np.max(displacements))


class MeanAtomDisplacement(StructuralMetric):
    """Mean per-atom displacement after alignment."""

    @property
    def metric_id(self) -> str:
        return "mean_atom_displacement"

    @property
    def display_name(self) -> str:
        return "Mean atom displacement"

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
        if not alignment.paired_atom_ids:
            return None

        atom_map_a = {a.atom_id: a for a in structure_a.atoms}
        atom_map_b = {b.atom_id: b for b in structure_b.atoms}

        coords_a = []
        coords_b = []
        for aid, bid in alignment.paired_atom_ids:
            if aid in atom_map_a and bid in atom_map_b:
                coords_a.append(atom_map_a[aid].coords)
                coords_b.append(atom_map_b[bid].coords)

        if not coords_a:
            return None

        coords_a = np.array(coords_a)
        coords_b = np.array(coords_b)

        if alignment.rotation is not None:
            coords_a = apply_alignment(coords_a, alignment)

        displacements = np.sqrt(np.sum((coords_a - coords_b) ** 2, axis=1))
        return float(np.mean(displacements))
