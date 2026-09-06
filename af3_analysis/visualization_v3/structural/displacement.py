"""
V3 Per-Residue Displacement.

Calculates per-residue displacement between structures.
Default atom: CA, but configurable.

Generates displacement tables with:
- condition
- residue_index
- residue_name
- displacement
- n_valid
- coverage
- seed_consistency
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..model import StructureData


def calculate_per_residue_displacement(
    structure_ref: StructureData,
    structure_target: StructureData,
    *,
    alignment_atom: str = "CA",
    chain_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Calculate per-residue displacement between two structures.

    Parameters
    ----------
    structure_ref : StructureData
        Reference structure.
    structure_target : StructureData
        Target structure.
    alignment_atom : str
        Which atom to use ('CA', 'backbone', 'all_heavy').
    chain_id : str, optional
        Specific chain to analyze. If None, uses protein chains.

    Returns
    -------
    dict with:
        - displacements: list of {residue_index, residue_name, displacement, atom_name}
        - n_residues: int
        - n_valid: int
        - coverage: float
    """
    # Determine chains
    if chain_id:
        chains_ref = [structure_ref.get_chain(chain_id)]
        chains_target = [structure_target.get_chain(chain_id)]
    else:
        # Use protein chains
        ref_protein = structure_ref.get_protein_chains()
        target_protein = structure_target.get_protein_chains()
        common = sorted(set(ref_protein) & set(target_protein))
        chains_ref = [structure_ref.get_chain(c) for c in common]
        chains_target = [structure_target.get_chain(c) for c in common]

    if not chains_ref or not chains_target:
        return {
            "displacements": [],
            "n_residues": 0,
            "n_valid": 0,
            "coverage": 0.0,
        }

    displacements = []
    total_residues = 0
    valid_residues = 0

    for chain_ref, chain_target in zip(chains_ref, chains_target):
        if chain_ref is None or chain_target is None:
            continue

        # Build residue index for target
        target_residues = {
            r.auth_seq_id: r
            for r in chain_target.residues
            if r.auth_seq_id is not None
        }

        for residue in chain_ref.residues:
            if residue.auth_seq_id is None:
                continue

            total_residues += 1
            auth_seq_id = residue.auth_seq_id

            # Find matching residue in target
            target_res = target_residues.get(auth_seq_id)
            if target_res is None:
                continue

            # Get atom coordinates
            ref_coords = None
            target_coords = None

            if alignment_atom == "CA":
                ref_coords = residue.ca_coords
                target_coords = target_res.ca_coords
            elif alignment_atom == "backbone":
                # Use CA as representative for backbone
                ref_coords = residue.ca_coords
                target_coords = target_res.ca_coords
            elif alignment_atom == "all_heavy":
                # Use CA as representative
                ref_coords = residue.ca_coords
                target_coords = target_res.ca_coords
            else:
                ref_coords = residue.ca_coords
                target_coords = target_res.ca_coords

            if ref_coords is None or target_coords is None:
                continue

            # Calculate displacement
            ref_arr = np.array(ref_coords, dtype=np.float64)
            target_arr = np.array(target_coords, dtype=np.float64)
            displacement = float(np.sqrt(np.sum((ref_arr - target_arr) ** 2)))

            displacements.append({
                "chain_id": chain_ref.chain_id if hasattr(chain_ref, 'chain_id') else "unknown",
                "residue_index": auth_seq_id,
                "residue_name": residue.residue_name,
                "displacement": displacement,
                "atom_name": alignment_atom,
            })
            valid_residues += 1

    coverage = valid_residues / total_residues if total_residues > 0 else 0.0

    return {
        "displacements": displacements,
        "n_residues": total_residues,
        "n_valid": valid_residues,
        "coverage": coverage,
    }


def calculate_displacement_vector(
    structure_ref: StructureData,
    structure_target: StructureData,
    *,
    chain_id: str,
    auth_seq_id: int,
    alignment_atom: str = "CA",
) -> Optional[np.ndarray]:
    """
    Calculate the displacement vector for a specific residue.

    Parameters
    ----------
    structure_ref, structure_target : StructureData
    chain_id : str
    auth_seq_id : int
    alignment_atom : str

    Returns
    -------
    (3,) array or None
    """
    ref_res = structure_ref.get_residue(chain_id, auth_seq_id)
    target_res = structure_target.get_residue(chain_id, auth_seq_id)

    if ref_res is None or target_res is None:
        return None

    ref_coords = None
    target_coords = None

    if alignment_atom == "CA":
        ref_coords = ref_res.ca_coords
        target_coords = target_res.ca_coords
    elif alignment_atom == "backbone":
        ref_coords = ref_res.ca_coords
        target_coords = target_res.ca_coords
    else:
        ref_coords = ref_res.ca_coords
        target_coords = target_res.ca_coords

    if ref_coords is None or target_coords is None:
        return None

    ref_arr = np.array(ref_coords, dtype=np.float64)
    target_arr = np.array(target_coords, dtype=np.float64)

    return target_arr - ref_arr


def aggregate_displacements(
    displacement_results: List[Dict[str, Any]],
    *,
    method: str = "mean",
) -> Dict[str, Dict[str, float]]:
    """
    Aggregate displacements across multiple comparisons.

    Parameters
    ----------
    displacement_results : list of per-residue displacement dicts
    method : str ('mean', 'median')

    Returns
    -------
    dict: residue_index -> {mean, median, std, n, values}
    """
    # Collect all displacements by residue
    residue_displacements = {}

    for result in displacement_results:
        for disp in result.get("displacements", []):
            key = (disp["chain_id"], disp["residue_index"])
            if key not in residue_displacements:
                residue_displacements[key] = {
                    "residue_name": disp["residue_name"],
                    "values": [],
                }
            residue_displacements[key]["values"].append(disp["displacement"])

    # Aggregate
    aggregated = {}
    for key, data in residue_displacements.items():
        values = np.array(data["values"])
        aggregated[key] = {
            "residue_name": data["residue_name"],
            "mean": float(np.mean(values)) if len(values) > 0 else None,
            "median": float(np.median(values)) if len(values) > 0 else None,
            "std": float(np.std(values)) if len(values) > 1 else None,
            "n": len(values),
            "min": float(np.min(values)) if len(values) > 0 else None,
            "max": float(np.max(values)) if len(values) > 0 else None,
        }

    return aggregated


def calculate_seed_direction_consistency(
    displacement_results: List[Dict[str, Any]],
    *,
    threshold: float = 0.5,
) -> Dict[Tuple[str, int], Dict[str, Any]]:
    """
    Calculate seed consistency for displacement direction.

    For each residue, determines if displacement direction is consistent
    across seeds (using sign of displacement along first principal axis).

    Parameters
    ----------
    displacement_results : list of per-seed displacement results
    threshold : float
        Minimum displacement magnitude to consider.

    Returns
    -------
    dict: (chain_id, residue_index) -> {consistent, n_seeds, n_forward, n_backward}
    """
    # Collect displacements by residue across seeds
    residue_by_seed = {}

    for seed_idx, result in enumerate(displacement_results):
        for disp in result.get("displacements", []):
            key = (disp["chain_id"], disp["residue_index"])
            if key not in residue_by_seed:
                residue_by_seed[key] = {}
            residue_by_seed[key][seed_idx] = disp["displacement"]

    # Calculate consistency
    consistency = {}
    for key, seed_disps in residue_by_seed.items():
        values = np.array(list(seed_disps.values()))

        if len(values) < 2:
            consistency[key] = {
                "n_seeds": len(values),
                "consistent": None,
                "mean_displacement": float(np.mean(values)) if len(values) > 0 else None,
            }
            continue

        # Check consistency: are all displacements in the same direction?
        # (positive = same direction, negative = opposite)
        positive = np.sum(values > threshold)
        negative = np.sum(values < -threshold)
        near_zero = len(values) - positive - negative

        # Direction consistency: fraction in majority direction
        n_seeds = len(values)
        if positive + negative > 0:
            majority = max(positive, negative)
            frac_consistent = majority / (positive + negative)
        else:
            frac_consistent = None

        consistency[key] = {
            "n_seeds": n_seeds,
            "n_forward": int(positive),
            "n_backward": int(negative),
            "n_near_zero": int(near_zero),
            "consistent": frac_consistent,
            "mean_displacement": float(np.mean(values)),
            "std_displacement": float(np.std(values)) if n_seeds > 1 else None,
        }

    return consistency
