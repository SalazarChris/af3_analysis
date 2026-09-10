"""
Structural alignment engine.

Provides pairwise alignment of NormalisedStructure objects using the
Kabsch algorithm for optimal rigid-body superposition. Handles missing
residues, unequal chain lengths, and reports explicit comparison status.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from af3_analysis.structural.enums import (
    AlignmentMethod,
    AtomSelection,
    ComparisonReason,
    ComparisonStatus,
)
from af3_analysis.structural.representation import AtomRecord, NormalisedStructure


# ---------------------------------------------------------------------------
# Alignment result
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AlignmentResult:
    """Result of aligning two NormalisedStructure objects."""
    status: ComparisonStatus
    reason: ComparisonReason

    # Atom pairing: list of (atom_id_A, atom_id_B)
    paired_atom_ids: Tuple[Tuple[int, int], ...] = ()
    n_common_atoms: int = 0

    # Common residue info: (chain_id, auth_seq_A, auth_seq_B)
    common_residues: Tuple[Tuple[str, int, int], ...] = ()
    n_common_residues: int = 0

    # Alignment transform
    rotation: Optional[np.ndarray] = None  # 3x3
    translation: Optional[np.ndarray] = None  # (3,)
    rmsd_pre_alignment: Optional[float] = None
    rmsd_post_alignment: Optional[float] = None

    # Missing atoms for diagnostics
    n_missing_in_b: int = 0
    n_missing_in_a: int = 0


# ---------------------------------------------------------------------------
# Atom selection
# ---------------------------------------------------------------------------

def _select_atoms(
    structure: NormalisedStructure,
    selection: AtomSelection,
    chain_ids: Optional[List[str]] = None,
) -> List[AtomRecord]:
    """Select atoms from a structure based on selection mode."""
    if selection == AtomSelection.CA:
        return structure.get_ca_atoms(chain_ids)
    elif selection == AtomSelection.BACKBONE:
        return structure.get_backbone_atoms(chain_ids)
    elif selection == AtomSelection.ALL_HEAVY:
        return structure.get_all_heavy_atoms(chain_ids)
    else:  # ALL
        if chain_ids is None:
            return list(structure.atoms)
        result = []
        for cid in chain_ids:
            result.extend(structure.get_chain(cid))
        return result


# ---------------------------------------------------------------------------
# Residue-level matching
# ---------------------------------------------------------------------------

def _build_residue_index(
    atoms: List[AtomRecord],
) -> Dict[Tuple[str, int], List[AtomRecord]]:
    """Build index: (chain_id, auth_seq_id) -> [AtomRecord, ...]."""
    index: Dict[Tuple[str, int], List[AtomRecord]] = {}
    for atom in atoms:
        key = (atom.chain_id, atom.auth_seq_id)
        if key not in index:
            index[key] = []
        index[key].append(atom)
    return index


def _match_residues(
    index_a: Dict[Tuple[str, int], List[AtomRecord]],
    index_b: Dict[Tuple[str, int], List[AtomRecord]],
) -> List[Tuple[str, int, int]]:
    """
    Find common residues between two structures.

    Returns list of (chain_id, seq_a, seq_b) where seq_a == seq_b
    (matched by chain_id and auth_seq_id).
    """
    common = []
    for key_a in index_a:
        if key_a in index_b:
            chain_id, seq_id = key_a
            common.append((chain_id, seq_id, seq_id))
    return sorted(common)


def _match_atoms(
    index_a: Dict[Tuple[str, int], List[AtomRecord]],
    index_b: Dict[Tuple[str, int], List[AtomRecord]],
    common_residues: List[Tuple[str, int, int]],
) -> List[Tuple[int, int]]:
    """
    Match atoms between two structures at matched residues.

    For each common residue, match atoms by atom_name.
    Returns list of (atom_id_a, atom_id_b).
    """
    pairs = []
    for chain_id, seq_a, seq_b in common_residues:
        atoms_a = {a.atom_name: a for a in index_a.get((chain_id, seq_a), [])}
        atoms_b = {a.atom_name: a for a in index_b.get((chain_id, seq_b), [])}
        for atom_name in atoms_a:
            if atom_name in atoms_b:
                pairs.append((atoms_a[atom_name].atom_id, atoms_b[atom_name].atom_id))
    return pairs


# ---------------------------------------------------------------------------
# Kabsch algorithm
# ---------------------------------------------------------------------------

def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Compute optimal rigid-body superposition of P onto Q using Kabsch algorithm.

    Parameters
    ----------
    P : (N, 3) array - coordinates to be superposed
    Q : (N, 3) array - reference coordinates

    Returns
    -------
    rotation : (3, 3) array
    translation : (3,) array
    rmsd : float - RMSD after superposition
    """
    assert P.shape == Q.shape
    assert P.shape[1] == 3

    # Centre both point sets
    centroid_p = P.mean(axis=0)
    centroid_q = Q.mean(axis=0)
    P_centered = P - centroid_p
    Q_centered = Q - centroid_q

    # Covariance matrix
    H = P_centered.T @ Q_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Ensure right-handed coordinate system
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.eye(3)
    sign_matrix[2, 2] = np.sign(d)

    # Optimal rotation
    rotation = Vt.T @ sign_matrix @ U.T

    # Optimal translation
    translation = centroid_q - rotation @ centroid_p

    # Compute RMSD after alignment (on centred coordinates)
    P_aligned = (rotation @ P_centered.T).T
    rmsd = np.sqrt(np.mean(np.sum((P_aligned - Q_centered) ** 2, axis=1)))

    return rotation, translation, rmsd


# ---------------------------------------------------------------------------
# Pre-alignment RMSD (no superposition)
# ---------------------------------------------------------------------------

def _compute_rmsd_no_alignment(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Compute RMSD without superposition (just centroid-aligned)."""
    centroid_a = coords_a.mean(axis=0)
    centroid_b = coords_b.mean(axis=0)
    diff = (coords_a - centroid_a) - (coords_b - centroid_b)
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


# ---------------------------------------------------------------------------
# Public alignment API
# ---------------------------------------------------------------------------

def align_structures(
    structure_a: NormalisedStructure,
    structure_b: NormalisedStructure,
    *,
    atom_selection: AtomSelection = AtomSelection.CA,
    chain_ids_a: Optional[List[str]] = None,
    chain_ids_b: Optional[List[str]] = None,
    min_common_atoms: int = 3,
    min_sequence_identity: float = 0.5,
) -> AlignmentResult:
    """
    Align two NormalisedStructure objects and return alignment result.

    Parameters
    ----------
    structure_a : NormalisedStructure
        First structure (target to be superposed).
    structure_b : NormalisedStructure
        Second structure (reference).
    atom_selection : AtomSelection
        Which atoms to use for alignment.
    chain_ids_a, chain_ids_b : list of str, optional
        Restrict alignment to these chains. If None, use all.
    min_common_atoms : int
        Minimum number of common atoms for a valid alignment.
    min_sequence_identity : float
        Minimum fraction of common residues relative to smaller structure.

    Returns
    -------
    AlignmentResult
        Alignment result with status, pairing, and transform.
    """
    # Select atoms
    atoms_a = _select_atoms(structure_a, atom_selection, chain_ids_a)
    atoms_b = _select_atoms(structure_b, atom_selection, chain_ids_b)

    if not atoms_a or not atoms_b:
        return AlignmentResult(
            status=ComparisonStatus.NOT_COMPARABLE,
            reason=ComparisonReason.INSUFFICIENT_ATOMS,
        )

    # Build residue indices
    index_a = _build_residue_index(atoms_a)
    index_b = _build_residue_index(atoms_b)

    # Find common residues
    common_residues = _match_residues(index_a, index_b)

    if not common_residues:
        return AlignmentResult(
            status=ComparisonStatus.NOT_COMPARABLE,
            reason=ComparisonReason.NO_COMMON_RESIDUES,
            n_missing_in_b=len(index_a),
            n_missing_in_a=len(index_b),
        )

    # Check sequence identity
    n_common = len(common_residues)
    n_smaller = min(len(index_a), len(index_b))
    identity = n_common / n_smaller if n_smaller > 0 else 0.0

    if identity < min_sequence_identity:
        return AlignmentResult(
            status=ComparisonStatus.NOT_COMPARABLE,
            reason=ComparisonReason.SEQUENCE_MISMATCH,
            n_common_residues=n_common,
        )

    # Match atoms
    atom_pairs = _match_atoms(index_a, index_b, common_residues)

    if len(atom_pairs) < min_common_atoms:
        return AlignmentResult(
            status=ComparisonStatus.NOT_COMPARABLE,
            reason=ComparisonReason.INSUFFICIENT_ATOMS,
            n_common_atoms=len(atom_pairs),
            n_common_residues=n_common,
        )

    # Extract paired coordinates
    atom_map_a = {a.atom_id: a for a in atoms_a}
    atom_map_b = {b.atom_id: b for b in atoms_b}

    coords_a = np.array([atom_map_a[aid].coords for aid, _ in atom_pairs])
    coords_b = np.array([atom_map_b[bid].coords for _, bid in atom_pairs])

    # Pre-alignment RMSD
    rmsd_pre = _compute_rmsd_no_alignment(coords_a, coords_b)

    # Kabsch alignment
    try:
        rotation, translation, rmsd_post = _kabsch(coords_a, coords_b)
    except np.linalg.LinAlgError:
        return AlignmentResult(
            status=ComparisonStatus.NOT_COMPARABLE,
            reason=ComparisonReason.INSUFFICIENT_ATOMS,
            n_common_atoms=len(atom_pairs),
            n_common_residues=n_common,
        )

    # Determine overall status
    n_missing_a = len(index_b) - n_common
    n_missing_b = len(index_a) - n_common

    # Check if all chains are represented in both structures
    chains_a = {key[0] for key in index_a}
    chains_b = {key[0] for key in index_b}
    all_chains_covered = chains_a == chains_b

    if n_common == min(len(index_a), len(index_b)) and all_chains_covered:
        status = ComparisonStatus.COMPARABLE
        reason = ComparisonReason.VALID
    else:
        status = ComparisonStatus.PARTIALLY_COMPARABLE
        reason = ComparisonReason.VALID

    return AlignmentResult(
        status=status,
        reason=reason,
        paired_atom_ids=tuple(atom_pairs),
        n_common_atoms=len(atom_pairs),
        common_residues=tuple(common_residues),
        n_common_residues=n_common,
        rotation=rotation,
        translation=translation,
        rmsd_pre_alignment=rmsd_pre,
        rmsd_post_alignment=rmsd_post,
        n_missing_in_b=n_missing_b,
        n_missing_in_a=n_missing_a,
    )


def apply_alignment(
    coords: np.ndarray,
    alignment: AlignmentResult,
) -> np.ndarray:
    """
    Apply an alignment transform to coordinates.

    Parameters
    ----------
    coords : (N, 3) array
    alignment : AlignmentResult with rotation and translation

    Returns
    -------
    (N, 3) array of transformed coordinates
    """
    if alignment.rotation is None or alignment.translation is None:
        return coords.copy()
    return (alignment.rotation @ coords.T).T + alignment.translation
