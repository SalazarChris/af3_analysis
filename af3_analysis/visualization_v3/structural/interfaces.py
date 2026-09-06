"""
V3 Interface Analysis.

Entity-aware interface analysis supporting:
- protein-DNA
- protein-RNA
- protein-ligand
- protein-protein
- protein-ion

For each interface calculates where possible:
- contacts
- distances
- contact frequency
- residue identity
- partner identity
- seed consistency
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from ..model import StructureData


@dataclass
class InterfaceContacts:
    """Interface contacts between two entities."""

    # Identity
    interface_id: str  # e.g., "protein-DNA:A:B"
    entity_a_type: str
    entity_b_type: str
    chain_a: str
    chain_b: str

    # Contacts
    contacts: List[Dict[str, Any]] = field(default_factory=list)
    # Each contact: {residue_a, residue_b, distance, atom_a, atom_b}

    # Summary
    n_contacts: int = 0
    min_distance: Optional[float] = None
    mean_distance: Optional[float] = None
    contact_frequency: float = 0.0  # fraction of possible contacts

    # Residues involved
    residues_a: Set[Tuple[str, int]] = field(default_factory=set)  # (chain_id, auth_seq_id)
    residues_b: Set[Tuple[str, int]] = field(default_factory=set)

    # Status
    status: str = "valid"  # valid, no_contacts, unsupported

    def __post_init__(self):
        if not hasattr(self, 'contacts'):
            object.__setattr__(self, 'contacts', [])
        if not hasattr(self, 'residues_a'):
            object.__setattr__(self, 'residues_a', set())
        if not hasattr(self, 'residues_b'):
            object.__setattr__(self, 'residues_b', set())
        self.n_contacts = len(self.contacts)

        if self.contacts:
            distances = [c["distance"] for c in self.contacts]
            self.min_distance = float(np.min(distances))
            self.mean_distance = float(np.mean(distances))
            self.contact_frequency = self.n_contacts / (
                len(self.residues_a) * len(self.residues_b)
            ) if self.residues_a and self.residues_b else 0.0

            for c in self.contacts:
                self.residues_a.add((c["chain_a"], c["residue_a"]))
                self.residues_b.add((c["chain_b"], c["residue_b"]))


def find_interface_contacts(
    structure: StructureData,
    *,
    chain_a: str,
    chain_b: str,
    threshold: float = 8.0,
    atom_a: str = "CA",
    atom_b: str = "CA",
) -> Optional[InterfaceContacts]:
    """
    Find interface contacts between two chains.

    Parameters
    ----------
    structure : StructureData
    chain_a, chain_b : str
        Chain IDs to compare.
    threshold : float
        Distance threshold for contact.
    atom_a, atom_b : str
        Atoms to use for distance.

    Returns
    -------
    InterfaceContacts or None if chains not found.
    """
    chain_obj_a = structure.get_chain(chain_a)
    chain_obj_b = structure.get_chain(chain_b)

    if chain_obj_a is None or chain_obj_b is None:
        return None

    # Determine entity types
    entity_a_type = chain_obj_a.entity_type if hasattr(chain_obj_a, 'entity_type') else "unknown"
    entity_b_type = chain_obj_b.entity_type if hasattr(chain_obj_b, 'entity_type') else "unknown"

    # Get atoms for each chain
    atoms_a = []
    atoms_b = []

    for residue in chain_obj_a.residues:
        if residue.auth_seq_id is None:
            continue

        coords = None
        if atom_a == "CA":
            coords = residue.ca_coords
        elif atom_a == "N":
            coords = residue.n_coords
        # Add other atoms as needed

        if coords is not None:
            atoms_a.append({
                "chain": chain_a,
                "residue": residue.auth_seq_id,
                "residue_name": residue.residue_name,
                "coords": np.array(coords),
            })

    for residue in chain_obj_b.residues:
        if residue.auth_seq_id is None:
            continue

        coords = None
        if atom_b == "CA":
            coords = residue.ca_coords
        elif atom_b == "N":
            coords = residue.n_coords

        if coords is not None:
            atoms_b.append({
                "chain": chain_b,
                "residue": residue.auth_seq_id,
                "residue_name": residue.residue_name,
                "coords": np.array(coords),
            })

    # Calculate contacts
    contacts = []
    for atom_a in atoms_a:
        for atom_b in atoms_b:
            dist = float(np.sqrt(np.sum((atom_a["coords"] - atom_b["coords"]) ** 2)))
            if dist < threshold:
                contacts.append({
                    "chain_a": atom_a["chain"],
                    "chain_b": atom_b["chain"],
                    "residue_a": atom_a["residue"],
                    "residue_b": atom_b["residue"],
                    "residue_name_a": atom_a["residue_name"],
                    "residue_name_b": atom_b["residue_name"],
                    "distance": dist,
                })

    # Determine interface type label
    if "polypeptide" in str(entity_a_type).lower() and "polydeoxyribonucleotide" in str(entity_b_type).lower():
        interface_type = "protein-DNA"
    elif "polypeptide" in str(entity_a_type).lower() and "polyribonucleotide" in str(entity_b_type).lower():
        interface_type = "protein-RNA"
    elif "polypeptide" in str(entity_a_type).lower() and "polypeptide" in str(entity_b_type).lower():
        interface_type = "protein-protein"
    else:
        interface_type = f"{entity_a_type}-{entity_b_type}"

    interface_id = f"{interface_type}:{chain_a}:{chain_b}"

    return InterfaceContacts(
        interface_id=interface_id,
        entity_a_type=entity_a_type,
        entity_b_type=entity_b_type,
        chain_a=chain_a,
        chain_b=chain_b,
        contacts=contacts,
        status="valid" if contacts else "no_contacts",
    )


def find_all_interfaces(
    structure: StructureData,
    *,
    threshold: float = 8.0,
) -> List[InterfaceContacts]:
    """
    Find all interfaces in a structure.

    Parameters
    ----------
    structure : StructureData
    threshold : float

    Returns
    -------
    list of InterfaceContacts
    """
    interfaces = []
    chains = structure.chain_ids

    # Group chains by entity type
    chain_types = {}
    for chain_id in chains:
        chain = structure.get_chain(chain_id)
        if chain:
            # Determine polymer type
            if hasattr(chain, 'polymer_type') and chain.polymer_type:
                if "polypeptide" in chain.polymer_type:
                    chain_types[chain_id] = "protein"
                elif "polydeoxyribonucleotide" in chain.polymer_type:
                    chain_types[chain_id] = "dna"
                elif "polyribonucleotide" in chain.polymer_type:
                    chain_types[chain_id] = "rna"
                else:
                    chain_types[chain_id] = "other"
            elif hasattr(chain, 'entity_type') and chain.entity_type:
                if chain.entity_type == "non-polymer":
                    chain_types[chain_id] = "ligand"
                elif chain.entity_type == "waters" or chain.entity_type == "branched":
                    chain_types[chain_id] = "ion"
                else:
                    chain_types[chain_id] = "unknown"
            else:
                chain_types[chain_id] = "unknown"

    # Find interfaces between chains of different types
    for i, chain_a in enumerate(chains):
        for chain_b in chains[i + 1:]:
            type_a = chain_types.get(chain_a, "unknown")
            type_b = chain_types.get(chain_b, "unknown")

            # Only look for specific interface types
            if type_a == type_b:
                continue

            if type_a in ("protein",) and type_b in ("dna", "rna", "ligand", "protein", "ion"):
                pass  # Will check
            elif type_b in ("protein",) and type_a in ("dna", "rna", "ligand", "ion"):
                pass  # Will check
            else:
                continue

            interface = find_interface_contacts(
                structure,
                chain_a=chain_a,
                chain_b=chain_b,
                threshold=threshold,
            )
            if interface:
                interfaces.append(interface)

    return interfaces


def calculate_interface_comparison(
    ref_interfaces: List[InterfaceContacts],
    target_interfaces: List[InterfaceContacts],
    *,
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """
    Compare interfaces between two structures.

    Parameters
    ----------
    ref_interfaces : list of InterfaceContacts from reference
    target_interfaces : list of InterfaceContacts from target
    threshold : float
        Contact frequency change threshold.

    Returns
    -------
    dict with interface changes
    """
    # Build interface lookup by ID
    ref_by_id = {iface.interface_id: iface for iface in ref_interfaces}
    target_by_id = {iface.interface_id: iface for iface in target_interfaces}

    all_ids = sorted(set(ref_by_id.keys()) | set(target_by_id.keys()))

    changes = []
    for iface_id in all_ids:
        ref = ref_by_id.get(iface_id)
        target = target_by_id.get(iface_id)

        if ref is None and target is None:
            continue

        if ref is None:
            changes.append({
                "interface_id": iface_id,
                "change": "gained",
                "ref_n_contacts": 0,
                "target_n_contacts": target.n_contacts if target else 0,
            })
        elif target is None:
            changes.append({
                "interface_id": iface_id,
                "change": "lost",
                "ref_n_contacts": ref.n_contacts,
                "target_n_contacts": 0,
            })
        else:
            # Compare contact counts
            n_ref = ref.n_contacts
            n_target = target.n_contacts

            if n_ref == 0 and n_target == 0:
                continue

            if n_target > n_ref * (1 + threshold):
                changes.append({
                    "interface_id": iface_id,
                    "change": "increased",
                    "ref_n_contacts": n_ref,
                    "target_n_contacts": n_target,
                    "fold_change": n_target / n_ref if n_ref > 0 else float("inf"),
                })
            elif n_target < n_ref * (1 - threshold) and n_ref > 0:
                changes.append({
                    "interface_id": iface_id,
                    "change": "decreased",
                    "ref_n_contacts": n_ref,
                    "target_n_contacts": n_target,
                    "fold_change": n_target / n_ref if n_ref > 0 else float("inf"),
                })
            else:
                changes.append({
                    "interface_id": iface_id,
                    "change": "unchanged",
                    "ref_n_contacts": n_ref,
                    "target_n_contacts": n_target,
                })

    return {
        "changes": changes,
        "n_gained": sum(1 for c in changes if c["change"] == "gained"),
        "n_lost": sum(1 for c in changes if c["change"] == "lost"),
        "n_changed": sum(1 for c in changes if c["change"] in ("increased", "decreased")),
    }
