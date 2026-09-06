"""
V3 RMSD Calculation.

Implements global structural RMSD calculation using Kabsch alignment.
Supports Cα RMSD and configurable alignment atom.

Returns both RMSD and coverage.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .comparability import (
    find_common_structural_space,
    CommonStructuralSpace,
    get_common_ca_coords,
    get_common_backbone_coords,
)


def calculate_rmsd(
    structure_a: Any,
    structure_b: Any,
    *,
    alignment_atom: str = "CA",
    min_common_atoms: int = 3,
    min_sequence_identity: float = 0.5,
    min_coverage: float = 0.80,
) -> Dict[str, Any]:
    """
    Calculate RMSD between two structures.

    Parameters
    ----------
    structure_a, structure_b : StructureData
        The two structures to compare.
    alignment_atom : str
        Which atom to use for alignment ('CA', 'backbone').
    min_common_atoms : int
        Minimum number of common atoms for valid alignment.
    min_sequence_identity : float
        Minimum fraction of common residues.
    min_coverage : float
        Minimum coverage fraction.

    Returns
    -------
    dict with:
        - rmsd: float or None
        - coverage: float
        - status: str
        - reason: str
        - n_common_atoms: int
        - n_common_residues: int
    """
    # Find common structural space
    space = find_common_structural_space(
        structure_a,
        structure_b,
        alignment_atom=alignment_atom,
        min_sequence_identity=min_sequence_identity,
    )

    result = {
        "rmsd": None,
        "coverage": space.coverage,
        "status": space.status,
        "reason": space.reason,
        "n_common_atoms": len(space.common_ca_atoms),
        "n_common_residues": sum(
            len(residues) for residues in space.common_residues.values()
        ),
    }

    # Check if we have enough atoms
    if len(space.common_ca_atoms) < min_common_atoms:
        result["reason"] = f"insufficient_common_atoms: {len(space.common_ca_atoms)} < {min_common_atoms}"
        return result

    # Get coordinates
    if alignment_atom == "CA":
        coords_a, coords_b = get_common_ca_coords(space, structure_a, structure_b)
    elif alignment_atom == "backbone":
        coords_a, coords_b = get_common_backbone_coords(space, structure_a, structure_b)
    else:
        coords_a, coords_b = get_common_ca_coords(space, structure_a, structure_b)

    if len(coords_a) < min_common_atoms:
        result["reason"] = "insufficient_atom_pairs"
        return result

    # Calculate RMSD with Kabsch alignment
    try:
        rmsd_value = _kabsch_rmsd(coords_a, coords_b)
        result["rmsd"] = float(rmsd_value)
        result["status"] = "comparable"
        result["reason"] = "valid"
    except Exception as e:
        result["reason"] = f"alignment_error: {str(e)[:100]}"

    return result


def _kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """
    Calculate RMSD after Kabsch optimal superposition.

    Parameters
    ----------
    coords_a, coords_b : (N, 3) arrays

    Returns
    -------
    float : RMSD in Angstroms
    """
    assert coords_a.shape == coords_b.shape
    assert coords_a.shape[1] == 3

    # Center both point sets
    centroid_a = coords_a.mean(axis=0)
    centroid_b = coords_b.mean(axis=0)
    centered_a = coords_a - centroid_a
    centered_b = coords_b - centroid_b

    # Covariance matrix
    H = centered_a.T @ centered_b

    # SVD
    try:
        U, S, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        # Fallback: no rotation
        aligned_a = centered_a
    else:
        # Ensure right-handed coordinate system
        d = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.eye(3)
        sign_matrix[2, 2] = np.sign(d)

        # Optimal rotation
        rotation = Vt.T @ sign_matrix @ U.T
        aligned_a = (rotation @ centered_a.T).T

    # Calculate RMSD
    diff = aligned_a - centered_b
    rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

    return float(rmsd)


def calculate_rmsd_matrix(
    structures: List[Any],
    *,
    alignment_atom: str = "CA",
    reference_structure: Optional[Any] = None,
    min_common_atoms: int = 3,
    min_sequence_identity: float = 0.5,
) -> np.ndarray:
    """
    Calculate pairwise RMSD matrix for all structures.

    Parameters
    ----------
    structures : list of StructureData
    alignment_atom : str
        Which atom to use.
    reference_structure : StructureData, optional
        If provided, calculate RMSD to this reference only.
    min_common_atoms : int
    min_sequence_identity : float

    Returns
    -------
    (N, N) array of RMSD values, or (N,) if reference provided.
    N = len(structures)
    """
    n = len(structures)

    if reference_structure is not None:
        # Calculate RMSD to reference for each structure
        rmsd_values = []
        for struct in structures:
            result = calculate_rmsd(
                struct,
                reference_structure,
                alignment_atom=alignment_atom,
                min_common_atoms=min_common_atoms,
                min_sequence_identity=min_sequence_identity,
            )
            rmsd_values.append(result["rmsd"] if result["rmsd"] is not None else np.nan)
        return np.array(rmsd_values)

    # Full pairwise matrix
    matrix = np.full((n, n), np.nan)

    for i in range(n):
        for j in range(i + 1, n):
            result = calculate_rmsd(
                structures[i],
                structures[j],
                alignment_atom=alignment_atom,
                min_common_atoms=min_common_atoms,
                min_sequence_identity=min_sequence_identity,
            )
            if result["rmsd"] is not None:
                matrix[i, j] = result["rmsd"]
                matrix[j, i] = result["rmsd"]

    return matrix


def calculate_rmsd_stats(
    rmsd_values: List[float],
) -> Dict[str, float]:
    """
    Calculate descriptive statistics for RMSD values.

    Parameters
    ----------
    rmsd_values : list of float

    Returns
    -------
    dict with mean, median, std, min, max, q25, q75
    """
    valid = [v for v in rmsd_values if v is not None and not np.isnan(v)]

    if not valid:
        return {
            "mean": np.nan,
            "median": np.nan,
            "std": np.nan,
            "min": np.nan,
            "max": np.nan,
            "q25": np.nan,
            "q75": np.nan,
            "n": 0,
        }

    arr = np.array(valid)
    return {
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "q25": float(np.percentile(arr, 25)),
        "q75": float(np.percentile(arr, 75)),
        "n": len(valid),
    }
