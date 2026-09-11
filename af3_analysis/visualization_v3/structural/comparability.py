"""
V3 Structural Comparability.

Determines the common structural space before calculating RMSD or displacement.
For every comparison calculates:
- Common chains
- Common residues
- Common atoms
- Common sequence
- Coverage

If coverage is below threshold, marks comparison as LOW_COVERAGE.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..model import StructureData


class ComparabilityError(Exception):
    """Comparability error."""
    pass


@dataclass
class CommonStructuralSpace:
    """The common structural space between two structures."""

    # Chains
    common_chain_ids: List[str]
    chain_to_entity: Dict[str, int]  # chain_id -> entity_id

    # Residues (by chain, by auth_seq_id)
    common_residues: Dict[str, List[int]]  # chain_id -> [auth_seq_id, ...]

    # Reference residues (for coverage)
    reference_residues: Dict[str, List[int]]  # chain_id -> [auth_seq_id, ...]

    # Common atoms (for alignment)
    common_ca_atoms: List[Tuple[int, int, int]]  # (chain_idx, auth_seq_id_a, auth_seq_id_b)
    common_backbone_atoms: List[Tuple[int, int, int]]

    # Sequence identity
    sequence_identity: float

    # Coverage
    coverage: float  # common / reference

    # Status
    status: str  # comparable, partially_comparable, not_comparable
    reason: str = ""

    def __post_init__(self):
        if not hasattr(self, 'common_chain_ids'):
            object.__setattr__(self, 'common_chain_ids', [])
        if not hasattr(self, 'chain_to_entity'):
            object.__setattr__(self, 'chain_to_entity', {})
        if not hasattr(self, 'common_residues'):
            object.__setattr__(self, 'common_residues', {})
        if not hasattr(self, 'reference_residues'):
            object.__setattr__(self, 'reference_residues', {})
        if not hasattr(self, 'common_ca_atoms'):
            object.__setattr__(self, 'common_ca_atoms', [])
        if not hasattr(self, 'common_backbone_atoms'):
            object.__setattr__(self, 'common_backbone_atoms', [])
        if not hasattr(self, 'sequence_identity'):
            object.__setattr__(self, 'sequence_identity', 0.0)
        if not hasattr(self, 'coverage'):
            object.__setattr__(self, 'coverage', 0.0)


def find_common_structural_space(
    structure_a: StructureData,
    structure_b: StructureData,
    *,
    alignment_atom: str = "CA",
    min_sequence_identity: float = 0.5,
    prefer_protein_chains: bool = True,
) -> CommonStructuralSpace:
    """
    Find the common structural space between two structures.

    Parameters
    ----------
    structure_a, structure_b : StructureData
        The two structures to compare.
    alignment_atom : str
        Which atom to use for alignment ('CA', 'backbone', 'all_heavy').
    min_sequence_identity : float
        Minimum fraction of common residues.
    prefer_protein_chains : bool
        If True, prefer protein chains for comparison.

    Returns
    -------
    CommonStructuralSpace
    """
    # Determine chains to compare
    chains_a = set(structure_a.chain_ids)
    chains_b = set(structure_b.chain_ids)

    if prefer_protein_chains:
        protein_a = set(structure_a.get_protein_chains())
        protein_b = set(structure_b.get_protein_chains())
        if protein_a and protein_b:
            common_chains = sorted(protein_a & protein_b)
        else:
            common_chains = sorted(chains_a & chains_b)
    else:
        common_chains = sorted(chains_a & chains_b)

    if not common_chains:
        return CommonStructuralSpace(
            common_chain_ids=[],
            chain_to_entity={},
            common_residues={},
            reference_residues={},
            common_ca_atoms=[],
            common_backbone_atoms=[],
            sequence_identity=0.0,
            coverage=0.0,
            status="not_comparable",
            reason="no_common_chains",
        )

    # Build common residue mapping per chain
    common_residues = {}
    reference_residues = {}
    chain_to_entity = {}
    total_common = 0
    total_ref = 0

    for chain_id in common_chains:
        chain_a = structure_a.get_chain(chain_id)
        chain_b = structure_b.get_chain(chain_id)

        if chain_a is None or chain_b is None:
            continue

        # Get entity_id for this chain
        if chain_a:
            chain_to_entity[chain_id] = chain_a.entity_id if hasattr(chain_a, 'entity_id') else 0

        # Get residues by auth_seq_id
        residues_a = {
            r.auth_seq_id: r
            for r in chain_a.residues
            if r.auth_seq_id is not None
        }
        residues_b = {
            r.auth_seq_id: r
            for r in chain_b.residues
            if r.auth_seq_id is not None
        }

        # Common residues
        common_ids = sorted(set(residues_a.keys()) & set(residues_b.keys()))
        common_residues[chain_id] = common_ids
        reference_residues[chain_id] = sorted(residues_a.keys())

        total_common += len(common_ids)
        total_ref += len(residues_a)

        # Track common CA atoms
        for auth_seq_id in common_ids:
            res_a = residues_a[auth_seq_id]
            res_b = residues_b[auth_seq_id]

            if alignment_atom == "CA":
                if res_a.ca_coords is not None and res_b.ca_coords is not None:
                    pass  # Will be counted in common_ca_atoms

    # Calculate sequence identity
    sequence_identity = (
        total_common / total_ref if total_ref > 0 else 0.0
    )

    # Check minimum sequence identity
    if sequence_identity < min_sequence_identity:
        return CommonStructuralSpace(
            common_chain_ids=common_chains,
            chain_to_entity=chain_to_entity,
            common_residues=common_residues,
            reference_residues=reference_residues,
            common_ca_atoms=[],
            common_backbone_atoms=[],
            sequence_identity=sequence_identity,
            coverage=total_common / total_ref if total_ref > 0 else 0.0,
            status="not_comparable",
            reason=f"sequence_identity {sequence_identity:.2%} < {min_sequence_identity:.2%}",
        )

    # Build common CA atom pairs
    common_ca_atoms = []
    common_backbone_atoms = []

    for chain_id in common_chains:
        chain_a = structure_a.get_chain(chain_id)
        chain_b = structure_b.get_chain(chain_id)

        if chain_a is None or chain_b is None:
            continue

        for auth_seq_id in common_residues.get(chain_id, []):
            # O(1) lookup via the chain's prebuilt residue index.
            res_a = chain_a.get_residue(auth_seq_id)
            res_b = chain_b.get_residue(auth_seq_id)

            if res_a and res_b:
                # CA atoms
                if res_a.ca_coords is not None and res_b.ca_coords is not None:
                    common_ca_atoms.append((chain_id, auth_seq_id, auth_seq_id))
                # Backbone atoms
                backbone_a = all([
                    res_a.n_coords is not None,
                    res_a.ca_coords is not None,
                    res_a.c_coords is not None,
                    res_a.o_coords is not None,
                ])
                backbone_b = all([
                    res_b.n_coords is not None,
                    res_b.ca_coords is not None,
                    res_b.c_coords is not None,
                    res_b.o_coords is not None,
                ])
                if backbone_a and backbone_b:
                    common_backbone_atoms.append((chain_id, auth_seq_id, auth_seq_id))

    # Calculate coverage
    coverage = total_common / total_ref if total_ref > 0 else 0.0

    # Determine status
    if len(common_ca_atoms) < 3:
        status = "not_comparable"
        reason = "insufficient_common_atoms"
    elif coverage < 0.80:
        status = "partially_comparable"
        reason = f"low_coverage {coverage:.2%}"
    else:
        status = "comparable"
        reason = "valid"

    return CommonStructuralSpace(
        common_chain_ids=common_chains,
        chain_to_entity=chain_to_entity,
        common_residues=common_residues,
        reference_residues=reference_residues,
        common_ca_atoms=common_ca_atoms,
        common_backbone_atoms=common_backbone_atoms,
        sequence_identity=sequence_identity,
        coverage=coverage,
        status=status,
        reason=reason,
    )


def has_sufficient_coverage(
    space: CommonStructuralSpace,
    threshold: float = 0.80,
) -> bool:
    """Check if coverage meets threshold."""
    return space.coverage >= threshold and space.status == "comparable"


# ---------------------------------------------------------------------------
# Prepared-structure fast path
#
# One-vs-many comparisons (the pairwise matrix, reference displacements)
# repeat the same chain/residue bookkeeping for every pair. The helpers
# below precompute each structure's chain selections and per-chain residue
# index once so the per-pair work is set intersection only. Results are
# identical to the non-prepared path (validated by tests).
# ---------------------------------------------------------------------------

@dataclass
class PreparedStructure:
    """Precomputed per-chain residue/CA inventories for repeated comparisons."""

    structure: StructureData
    # chain_id -> {auth_seq_id: ca coords or None}; protein chains only
    # (empty when the structure has no protein chains)
    protein_ca_by_res: Dict[str, Dict[int, Optional[Tuple[float, float, float]]]]
    # chain_id -> {auth_seq_id: ca coords or None}; all chains
    all_ca_by_res: Dict[str, Dict[int, Optional[Tuple[float, float, float]]]]


def prepare_structure(
    structure: StructureData,
) -> PreparedStructure:
    """Precompute residue inventories for a structure.

    Both protein-only and all-chain inventories are kept so the pair-time
    chain selection can replicate ``find_common_structural_space`` exactly:
    protein chains of both structures when both have them, otherwise all
    chains of both.
    """
    all_ca_by_res: Dict[str, Dict[int, Optional[Tuple[float, float, float]]]] = {}
    for chain_id in structure.chain_ids:
        chain = structure.get_chain(chain_id)
        if chain is None:
            continue
        all_ca_by_res[chain_id] = {
            r.auth_seq_id: r.ca_coords
            for r in chain.residues
            if r.auth_seq_id is not None
        }

    protein_ids = set(structure.get_protein_chains())
    protein_ca_by_res = {
        cid: m for cid, m in all_ca_by_res.items() if cid in protein_ids
    }
    return PreparedStructure(
        structure=structure,
        protein_ca_by_res=protein_ca_by_res,
        all_ca_by_res=all_ca_by_res,
    )


def _pair_ca_maps(
    prepared_a: PreparedStructure,
    prepared_b: PreparedStructure,
) -> Tuple[
    Dict[str, Dict[int, Optional[Tuple[float, float, float]]]],
    Dict[str, Dict[int, Optional[Tuple[float, float, float]]]],
]:
    """Chain selection mirroring find_common_structural_space."""
    if prepared_a.protein_ca_by_res and prepared_b.protein_ca_by_res:
        return prepared_a.protein_ca_by_res, prepared_b.protein_ca_by_res
    return prepared_a.all_ca_by_res, prepared_b.all_ca_by_res


def find_common_space_prepared(
    prepared_a: PreparedStructure,
    prepared_b: PreparedStructure,
    *,
    min_sequence_identity: float = 0.5,
) -> Tuple[str, float, Dict[str, List[int]], int]:
    """Common-space determination between two prepared structures.

    Returns ``(status, coverage, common_residues_by_chain, n_common_residues)``
    where status is ``comparable`` / ``partially_comparable`` /
    ``not_comparable`` using the same thresholds and reasons as
    ``find_common_structural_space`` (sequence identity gate, then the 0.80
    coverage boundary and the >=3 atom minimum).
    """
    ca_a_all, ca_b_all = _pair_ca_maps(prepared_a, prepared_b)
    common_chains = sorted(set(ca_a_all) & set(ca_b_all))
    if not common_chains:
        return "not_comparable", 0.0, {}, 0

    common_residues: Dict[str, List[int]] = {}
    total_common = 0
    total_ref = 0
    n_ca_pairs = 0
    for chain_id in common_chains:
        ca_a = ca_a_all[chain_id]
        ca_b = ca_b_all[chain_id]
        common_ids = sorted(set(ca_a) & set(ca_b))
        common_residues[chain_id] = common_ids
        total_common += len(common_ids)
        total_ref += len(ca_a)
        n_ca_pairs += sum(
            1 for rid in common_ids
            if ca_a[rid] is not None and ca_b[rid] is not None
        )

    if total_ref == 0:
        return "not_comparable", 0.0, common_residues, 0

    sequence_identity = total_common / total_ref
    coverage = sequence_identity

    if sequence_identity < min_sequence_identity:
        return "not_comparable", coverage, common_residues, n_ca_pairs
    if n_ca_pairs < 3:
        return "not_comparable", coverage, common_residues, n_ca_pairs
    if coverage < 0.80:
        return "partially_comparable", coverage, common_residues, n_ca_pairs
    return "comparable", coverage, common_residues, n_ca_pairs


def kabsch_rmsd_prepared(
    prepared_a: PreparedStructure,
    prepared_b: PreparedStructure,
    common_residues: Dict[str, List[int]],
) -> Optional[float]:
    """CA Kabsch RMSD over the given common residues, or None if too few.

    Coordinates are gathered with O(1) dict lookups from the prepared
    inventories; alignment uses the same ``_kabsch_rmsd`` implementation as
    the standard path.
    """
    coords_a: List[Tuple[float, float, float]] = []
    coords_b: List[Tuple[float, float, float]] = []
    ca_a_all, ca_b_all = _pair_ca_maps(prepared_a, prepared_b)
    for chain_id, ids in common_residues.items():
        ca_a = ca_a_all[chain_id]
        ca_b = ca_b_all[chain_id]
        for rid in ids:
            ca = ca_a[rid]
            cb = ca_b[rid]
            if ca is not None and cb is not None:
                coords_a.append(ca)
                coords_b.append(cb)

    if len(coords_a) < 3:
        return None

    from .rmsd import _kabsch_rmsd

    return _kabsch_rmsd(
        np.array(coords_a, dtype=np.float64),
        np.array(coords_b, dtype=np.float64),
    )


def get_common_ca_coords(
    space: CommonStructuralSpace,
    structure_a: StructureData,
    structure_b: StructureData,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get paired CA coordinates for common residues.

    Returns two (N, 3) arrays of coordinates.
    """
    coords_a = []
    coords_b = []

    for chain_id, auth_seq_a, auth_seq_b in space.common_ca_atoms:
        chain_a = structure_a.get_chain(chain_id)
        chain_b = structure_b.get_chain(chain_id)

        if chain_a is None or chain_b is None:
            continue

        # O(1) lookup via the chain's prebuilt residue index.
        res_a = chain_a.get_residue(auth_seq_a)
        res_b = chain_b.get_residue(auth_seq_b)

        if res_a and res_b:
            if res_a.ca_coords is not None and res_b.ca_coords is not None:
                coords_a.append(res_a.ca_coords)
                coords_b.append(res_b.ca_coords)

    if not coords_a:
        return np.array([]).reshape(0, 3), np.array([]).reshape(0, 3)

    return (
        np.array(coords_a, dtype=np.float64),
        np.array(coords_b, dtype=np.float64),
    )


def get_common_backbone_coords(
    space: CommonStructuralSpace,
    structure_a: StructureData,
    structure_b: StructureData,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get paired backbone coordinates for common residues.

    Returns two (N, 3) arrays of coordinates.
    """
    coords_a = []
    coords_b = []

    for chain_id, auth_seq_a, auth_seq_b in space.common_backbone_atoms:
        chain_a = structure_a.get_chain(chain_id)
        chain_b = structure_b.get_chain(chain_id)

        if chain_a is None or chain_b is None:
            continue

        # O(1) lookup via the chain's prebuilt residue index.
        res_a = chain_a.get_residue(auth_seq_a)
        res_b = chain_b.get_residue(auth_seq_b)

        if res_a and res_b:
            # Order: N, CA, C, O
            for atom_name in ["n_coords", "ca_coords", "c_coords", "o_coords"]:
                if getattr(res_a, atom_name) is not None and getattr(res_b, atom_name) is not None:
                    coords_a.append(getattr(res_a, atom_name))
                    coords_b.append(getattr(res_b, atom_name))

    if not coords_a:
        return np.array([]).reshape(0, 3), np.array([]).reshape(0, 3)

    return (
        np.array(coords_a, dtype=np.float64),
        np.array(coords_b, dtype=np.float64),
    )
