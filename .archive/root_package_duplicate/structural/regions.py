"""
Region definition and selection for structural analysis.

Provides a generic mechanism for users to define structural regions
and built-in region generators for entity-type-based selections.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from af3_analysis.structural.config import RegionDefinition
from af3_analysis.structural.representation import AtomRecord, NormalisedStructure


@dataclass(frozen=True)
class RegionSpec:
    """
    Resolved region specification for atom selection.

    Created from a RegionDefinition or built-in generator.
    """
    label: Optional[str]
    chain_ids: Optional[List[str]]
    entity_types: Optional[List[str]]
    seq_range: Optional[tuple]
    seq_list: Optional[List[int]]
    atom_names: Optional[List[str]]

    def select(self, structure: NormalisedStructure) -> List[AtomRecord]:
        """Return atoms matching this region in the given structure."""
        atoms = list(structure.atoms)

        if self.chain_ids is not None:
            atoms = [a for a in atoms if a.chain_id in self.chain_ids]

        if self.entity_types is not None:
            # Map entity types to chain IDs
            valid_chains = set()
            for chain in structure.chains:
                if chain.entity_type in self.entity_types:
                    valid_chains.add(chain.chain_id)
                elif chain.polymer_type:
                    if "polypeptide" in chain.polymer_type and "protein" in self.entity_types:
                        valid_chains.add(chain.chain_id)
                    elif "polydeoxyribonucleotide" in chain.polymer_type and "dna" in self.entity_types:
                        valid_chains.add(chain.chain_id)
                    elif "polyribonucleotide" in chain.polymer_type and "rna" in self.entity_types:
                        valid_chains.add(chain.chain_id)
            atoms = [a for a in atoms if a.chain_id in valid_chains]

        if self.seq_range is not None:
            start, end = self.seq_range
            atoms = [
                a for a in atoms
                if a.auth_seq_id is not None and start <= a.auth_seq_id <= end
            ]

        if self.seq_list is not None:
            seq_set = set(self.seq_list)
            atoms = [
                a for a in atoms
                if a.auth_seq_id is not None and a.auth_seq_id in seq_set
            ]

        if self.atom_names is not None:
            name_set = set(self.atom_names)
            atoms = [a for a in atoms if a.atom_name in name_set]

        return atoms


def from_config_region(region_def: RegionDefinition) -> RegionSpec:
    """Convert a RegionDefinition config to a RegionSpec."""
    return RegionSpec(
        label=region_def.label,
        chain_ids=region_def.chain_ids,
        entity_types=region_def.entity_types,
        seq_range=region_def.seq_range,
        seq_list=region_def.seq_list,
        atom_names=region_def.atom_names,
    )


def protein_chains_region(structure: NormalisedStructure) -> RegionSpec:
    """Select all protein-type chains."""
    return RegionSpec(
        label="protein_chains",
        chain_ids=structure.get_protein_chains(),
        entity_types=None,
        seq_range=None,
        seq_list=None,
        atom_names=None,
    )


def nucleic_acid_chains_region(structure: NormalisedStructure) -> RegionSpec:
    """Select all nucleic-acid-type chains."""
    return RegionSpec(
        label="nucleic_acid_chains",
        chain_ids=structure.get_nucleic_acid_chains(),
        entity_types=None,
        seq_range=None,
        seq_list=None,
        atom_names=None,
    )


def interface_region(
    structure_a: NormalisedStructure,
    structure_b: NormalisedStructure,
    chain_a: str,
    chain_b: str,
    cutoff_angstrom: float = 8.0,
) -> RegionSpec:
    """
    Detect interface residues within cutoff between two chains.

    Uses structure_b (reference) for detection.
    """
    from collections import defaultdict
    import numpy as np

    # Get CA atoms for both chains in structure_b
    ca_atoms = {}
    for atom in structure_b.atoms:
        if atom.atom_name == "CA" and atom.chain_id in (chain_a, chain_b):
            key = (atom.chain_id, atom.auth_seq_id or 0)
            ca_atoms[key] = atom

    # Find interface residues
    chain_a_keys = [k for k in ca_atoms if k[0] == chain_a]
    chain_b_keys = [k for k in ca_atoms if k[0] == chain_b]

    interface_keys = set()
    for ka in chain_a_keys:
        for kb in chain_b_keys:
            dist = np.sqrt(np.sum((ca_atoms[ka].xyz - ca_atoms[kb].xyz) ** 2))
            if dist < cutoff_angstrom:
                interface_keys.add(ka)
                interface_keys.add(kb)

    seq_list = [k[1] for k in interface_keys]
    chain_ids = list({k[0] for k in interface_keys})

    return RegionSpec(
        label=f"interface_{chain_a}_{chain_b}",
        chain_ids=chain_ids if chain_ids else None,
        entity_types=None,
        seq_range=None,
        seq_list=seq_list if seq_list else None,
        atom_names=None,
    )
