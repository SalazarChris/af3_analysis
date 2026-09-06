"""
V3 Contact Maps.

Creates residue-level contact maps.
For each structure:
- contact(i,j) = distance(i,j) < configured threshold

Calculates:
- reference contact map
- condition contact map
- gained contacts
- lost contacts
- unchanged contacts

Creates delta_contact_map where:
+1 = gained
 0 = unchanged
-1 = lost
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from ..model import StructureData


@dataclass
class ContactMap:
    """Residue-level contact map."""
    structure_id: str
    condition_id: str
    seed: int
    sample: int
    contacts: Dict[Tuple[str, int, str, int], float]  # ((chain_a, seq_a), (chain_b, seq_b)) -> distance
    threshold: float
    n_contacts: int

    def is_contact(self, chain_a: str, seq_a: int, chain_b: str, seq_b: int) -> bool:
        """Check if two residues are in contact."""
        key = (chain_a, seq_a, chain_b, seq_b)
        if key in self.contacts:
            return self.contacts[key] < self.threshold
        # Try reverse order
        rev_key = (chain_b, seq_b, chain_a, seq_a)
        if rev_key in self.contacts:
            return self.contacts[rev_key] < self.threshold
        return False


def calculate_contact_map(
    structure: StructureData,
    *,
    threshold: float = 8.0,
    chain_id: Optional[str] = None,
    atom_name: str = "CA",
) -> ContactMap:
    """
    Calculate residue-level contact map for a structure.

    Parameters
    ----------
    structure : StructureData
    threshold : float
        Distance threshold for contact (Angstrom).
    chain_id : str, optional
        Specific chain. If None, uses all chains.
    atom_name : str
        Atom to use for distance calculation.

    Returns
    -------
    ContactMap
    """
    # Get atoms
    if chain_id:
        chains = [structure.get_chain(chain_id)]
    else:
        chains = [structure.get_chain(cid) for cid in structure.chain_ids]

    # Get CA atoms per residue
    residue_atoms = {}  # (chain_id, auth_seq_id) -> coordinates or None

    for chain in chains:
        if chain is None:
            continue

        for residue in chain.residues:
            if residue.auth_seq_id is None:
                continue

            coords = None
            if atom_name == "CA":
                coords = residue.ca_coords
            elif atom_name == "N":
                coords = residue.n_coords
            # Add other atoms as needed

            if coords is not None:
                residue_atoms[(chain.chain_id, residue.auth_seq_id)] = np.array(coords)

    # Calculate contacts
    contacts = {}
    keys = list(residue_atoms.keys())
    n_residues = len(keys)

    for i in range(n_residues):
        for j in range(i + 1, n_residues):
            key_i = keys[i]
            key_j = keys[j]

            # Skip same residue
            if key_i[0] == key_j[0] and key_i[1] == key_j[1]:
                continue

            coords_i = residue_atoms[key_i]
            coords_j = residue_atoms[key_j]

            if coords_i is None or coords_j is None:
                continue

            dist = float(np.sqrt(np.sum((coords_i - coords_j) ** 2)))
            contacts[(key_i[0], key_i[1], key_j[0], key_j[1])] = dist

    n_contacts = sum(1 for d in contacts.values() if d < threshold)

    return ContactMap(
        structure_id=structure.prediction_id,
        condition_id=structure.condition_id,
        seed=structure.seed,
        sample=structure.sample,
        contacts=contacts,
        threshold=threshold,
        n_contacts=n_contacts,
    )


def calculate_contact_difference(
    ref_map: ContactMap,
    target_map: ContactMap,
    *,
    chain_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Calculate contact map difference between two structures.

    Returns:
    - delta_map: dict of (chain_a, seq_a, chain_b, seq_b) -> +1/-1/0
    - gained: list of gained contacts
    - lost: list of lost contacts
    - unchanged: list of unchanged contacts
    - n_gained: int
    - n_lost: int
    - n_unchanged: int
    - pct_changed: float
    """
    # Get all residue pairs from both maps
    all_pairs = set()

    for key in ref_map.contacts.keys():
        all_pairs.add(key)
    for key in target_map.contacts.keys():
        all_pairs.add(key)

    delta_map = {}
    gained = []
    lost = []
    unchanged = []

    for pair in all_pairs:
        chain_a, seq_a, chain_b, seq_b = pair

        # Apply chain filter
        if chain_id and chain_a != chain_id and chain_b != chain_id:
            continue

        ref_dist = ref_map.contacts.get(pair)
        target_dist = target_map.contacts.get(pair)

        # Handle missing distances
        if ref_dist is None:
            ref_contact = False
        else:
            ref_contact = ref_dist < ref_map.threshold

        if target_dist is None:
            target_contact = False
        else:
            target_contact = target_dist < target_map.threshold

        if ref_contact and target_contact:
            delta_map[pair] = 0
            unchanged.append(pair)
        elif ref_contact and not target_contact:
            delta_map[pair] = -1
            lost.append(pair)
        elif not ref_contact and target_contact:
            delta_map[pair] = 1
            gained.append(pair)
        else:
            # Neither in contact
            delta_map[pair] = 0
            unchanged.append(pair)

    n_gained = len(gained)
    n_lost = len(lost)
    n_unchanged = len(unchanged)
    n_total = len(all_pairs)

    pct_changed = (n_gained + n_lost) / n_total if n_total > 0 else 0.0

    return {
        "delta_map": delta_map,
        "gained": gained,
        "lost": lost,
        "unchanged": unchanged,
        "n_gained": n_gained,
        "n_lost": n_lost,
        "n_unchanged": n_unchanged,
        "n_total": n_total,
        "pct_changed": pct_changed,
        "ref_n_contacts": ref_map.n_contacts,
        "target_n_contacts": target_map.n_contacts,
    }


def aggregate_contact_changes(
    contact_differences: List[Dict[str, Any]],
    *,
    min_occurrences: int = 1,
) -> Dict[str, Any]:
    """
    Aggregate contact changes across multiple comparisons.

    Parameters
    ----------
    contact_differences : list of contact difference dicts
    min_occurrences : int
        Minimum number of times a change must occur to be included.

    Returns
    -------
    dict with:
        - recurring_gained: residue pairs that are frequently gained
        - recurring_lost: residue pairs that are frequently lost
        - consistent_gained: pairs gained in all comparisons
        - consistent_lost: pairs lost in all comparisons
    """
    gained_counts = {}
    lost_counts = {}

    n_comparisons = len(contact_differences)

    for diff in contact_differences:
        for pair in diff.get("gained", []):
            gained_counts[pair] = gained_counts.get(pair, 0) + 1
        for pair in diff.get("lost", []):
            lost_counts[pair] = lost_counts.get(pair, 0) + 1

    # Filter by minimum occurrences
    recurring_gained = [
        (pair, count)
        for pair, count in gained_counts.items()
        if count >= min_occurrences
    ]
    recurring_lost = [
        (pair, count)
        for pair, count in lost_counts.items()
        if count >= min_occurrences
    ]

    # Consistent changes (occur in all comparisons)
    consistent_gained = [pair for pair, count in gained_counts.items() if count == n_comparisons]
    consistent_lost = [pair for pair, count in lost_counts.items() if count == n_comparisons]

    return {
        "recurring_gained": sorted(recurring_gained, key=lambda x: -x[1]),
        "recurring_lost": sorted(recurring_lost, key=lambda x: -x[1]),
        "consistent_gained": consistent_gained,
        "consistent_lost": consistent_lost,
        "n_gained_total": sum(1 for diff in contact_differences for _ in diff.get("gained", [])),
        "n_lost_total": sum(1 for diff in contact_differences for _ in diff.get("lost", [])),
    }
