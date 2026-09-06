"""
V3 Local/Region Geometry.

Supports configurable local regions for:
- Local RMSD
- Local displacement
- Neighboring distances
- Contact changes
- Local confidence
- Interface proximity

Region definitions come from metadata/configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..model import StructureData


@dataclass
class RegionResult:
    """Results for a structural region."""

    region_id: str
    region_label: str
    condition_id: str
    prediction_id: str
    seed: int
    sample: int

    # Geometry
    n_atoms: int = 0
    n_residues: int = 0
    local_rmsd: Optional[float] = None  # RMSD relative to reference
    centroid: Optional[np.ndarray] = None  # (3,) centroid coordinates
    radius_gyration: Optional[float] = None

    # Confidence
    local_plddt_mean: Optional[float] = None

    # Interface proximity
    interface_distance_min: Optional[float] = None
    interface_distance_mean: Optional[float] = None
    n_interface_contacts: int = 0

    # Status
    status: str = "valid"
    reason: str = ""


def calculate_local_geometry(
    structure: StructureData,
    region_definition: Dict[str, Any],
    *,
    reference_structure: Optional[StructureData] = None,
    alignment_atom: str = "CA",
) -> RegionResult:
    """
    Calculate local geometry for a region.

    Parameters
    ----------
    structure : StructureData
    region_definition : dict with region specification
        - label: str
        - chain_ids: list of str
        - seq_range: (start, end)
        - seq_list: list of int
        - radius: float (for spherical region around residue)
        - center_residue: (chain_id, auth_seq_id)
    reference_structure : StructureData, optional
        Reference for calculating local RMSD.
    alignment_atom : str

    Returns
    -------
    RegionResult
    """
    label = region_definition.get("label", "unnamed")
    chain_ids = region_definition.get("chain_ids")
    seq_range = region_definition.get("seq_range")
    seq_list = region_definition.get("seq_list")
    radius = region_definition.get("radius")
    center_residue = region_definition.get("center_residue")  # (chain_id, auth_seq_id)

    # Collect atoms in region
    region_atoms = []

    # Determine which chains to use
    target_chains = chain_ids
    if target_chains is None:
        target_chains = structure.chain_ids

    # Determine which residues to include
    include_residues = None  # (chain_id, auth_seq_id) pairs

    if seq_list is not None:
        include_residues = set()
        for cid in target_chains:
            for seq_id in seq_list:
                include_residues.add((cid, seq_id))
    elif seq_range is not None:
        start, end = seq_range
        include_residues = set()
        for cid in target_chains:
            for seq_id in range(start, end + 1):
                include_residues.add((cid, seq_id))

    # If radius specified, find residues within radius of center
    if radius is not None and center_residue is not None:
        center_chain, center_seq = center_residue
        center_res = structure.get_residue(center_chain, center_seq)

        if center_res is not None:
            center_coords = None
            if alignment_atom == "CA":
                center_coords = center_res.ca_coords
            elif alignment_atom == "backbone":
                center_coords = center_res.ca_coords

            if center_coords is not None:
                center_arr = np.array(center_coords, dtype=np.float64)

                # Find all residues within radius
                neighbor_residues = set()
                for cid in target_chains:
                    chain = structure.get_chain(cid)
                    if chain:
                        for res in chain.residues:
                            if res.auth_seq_id is None:
                                continue
                            res_coords = None
                            if alignment_atom == "CA":
                                res_coords = res.ca_coords
                            if res_coords is not None:
                                res_arr = np.array(res_coords, dtype=np.float64)
                                dist = float(np.sqrt(np.sum((res_arr - center_arr) ** 2)))
                                if dist <= radius:
                                    neighbor_residues.add((cid, res.auth_seq_id))

                include_residues = neighbor_residues

    # Collect atoms
    for cid in target_chains:
        chain = structure.get_chain(cid)
        if chain is None:
            continue

        for residue in chain.residues:
            if residue.auth_seq_id is None:
                continue

            # Check if this residue is in the region
            in_region = False
            if include_residues is None:
                # Include all residues in selected chains
                in_region = True
            elif (cid, residue.auth_seq_id) in include_residues:
                in_region = True

            if not in_region:
                continue

            # Get atom coordinates
            coords = None
            if alignment_atom == "CA":
                coords = residue.ca_coords
            elif alignment_atom == "backbone":
                coords = residue.ca_coords
            # Add other atoms as needed

            if coords is not None:
                region_atoms.append({
                    "chain_id": cid,
                    "auth_seq_id": residue.auth_seq_id,
                    "residue_name": residue.residue_name,
                    "coords": np.array(coords),
                    "plddt": None,  # Would come from AF3 data
                })

    # Calculate region properties
    n_atoms = len(region_atoms)
    n_residues = len(set(
        (a["chain_id"], a["auth_seq_id"]) for a in region_atoms
    ))

    # Centroid
    if n_atoms > 0:
        coords_arr = np.array([a["coords"] for a in region_atoms])
        centroid = coords_arr.mean(axis=0)

        # Radius of gyration
        diff = coords_arr - centroid
        radius_gyration = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    else:
        centroid = None
        radius_gyration = None

    # Local RMSD (if reference provided)
    local_rmsd = None
    if reference_structure is not None:
        # Get corresponding atoms in reference
        ref_atoms = []
        for atom_info in region_atoms:
            cid = atom_info["chain_id"]
            seq_id = atom_info["auth_seq_id"]

            ref_chain = reference_structure.get_chain(cid)
            if ref_chain is None:
                continue

            ref_res = ref_chain.get_residue(seq_id)
            if ref_res is None:
                continue

            ref_coords = None
            if alignment_atom == "CA":
                ref_coords = ref_res.ca_coords
            elif alignment_atom == "backbone":
                ref_coords = ref_res.ca_coords

            if ref_coords is not None:
                ref_atoms.append(np.array(ref_coords))

        if len(ref_atoms) > 0 and n_atoms > 0:
            coords_arr = np.array([a["coords"] for a in region_atoms])
            if len(ref_atoms) == len(coords_arr):
                diff = coords_arr - np.array(ref_atoms)
                local_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    # Local pLDDT (would come from AF3 data)
    local_plddt_mean = None
    # if structure.plddt_per_residue is available:
    #     ...

    # Interface proximity
    interface_distance_min = None
    interface_distance_mean = None
    n_interface_contacts = 0

    # Determine if this region is near an interface
    # This would require identifying interfaces first

    status = "valid" if n_atoms > 0 else "empty"
    reason = "" if n_atoms > 0 else "no_atoms_in_region"

    return RegionResult(
        region_id=f"{structure.prediction_id}:{label}",
        region_label=label,
        condition_id=structure.condition_id,
        prediction_id=structure.prediction_id,
        seed=structure.seed,
        sample=structure.sample,
        n_atoms=n_atoms,
        n_residues=n_residues,
        local_rmsd=local_rmsd,
        centroid=centroid,
        radius_gyration=radius_gyration,
        local_plddt_mean=local_plddt_mean,
        interface_distance_min=interface_distance_min,
        interface_distance_mean=interface_distance_mean,
        n_interface_contacts=n_interface_contacts,
        status=status,
        reason=reason,
    )


def calculate_region_comparison(
    ref_region: RegionResult,
    target_region: RegionResult,
) -> Dict[str, Any]:
    """
    Compare two region results.

    Parameters
    ----------
    ref_region : RegionResult
    target_region : RegionResult

    Returns
    -------
    dict with comparison
    """
    if ref_region.status != "valid" or target_region.status != "valid":
        return {
            "status": "invalid",
            "reason": "one or both regions invalid",
        }

    n_atoms = min(ref_region.n_atoms, target_region.n_atoms)
    n_residues = min(ref_region.n_residues, target_region.n_residues)

    # RMSD difference
    rmsd_diff = None
    if ref_region.local_rmsd is not None and target_region.local_rmsd is not None:
        rmsd_diff = target_region.local_rmsd - ref_region.local_rmsd

    # Centroid displacement
    centroid_displacement = None
    if ref_region.centroid is not None and target_region.centroid is not None:
        centroid_displacement = float(np.sqrt(np.sum(
            (target_region.centroid - ref_region.centroid) ** 2
        )))

    # Radius of gyration change
    rg_change = None
    if ref_region.radius_gyration is not None and target_region.radius_gyration is not None:
        rg_change = target_region.radius_gyration - ref_region.radius_gyration

    return {
        "region_id": ref_region.region_id,
        "region_label": ref_region.region_label,
        "n_atoms": n_atoms,
        "n_residues": n_residues,
        "local_rmsd_ref": ref_region.local_rmsd,
        "local_rmsd_target": target_region.local_rmsd,
        "local_rmsd_diff": rmsd_diff,
        "centroid_displacement": centroid_displacement,
        "rg_change": rg_change,
        "local_plddt_ref": ref_region.local_plddt_mean,
        "local_plddt_target": target_region.local_plddt_mean,
        "status": "valid",
    }
