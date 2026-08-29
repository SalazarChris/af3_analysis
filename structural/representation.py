"""
Normalized structural representation for AF3 predicted structures.

Provides a project-agnostic internal model for parsed mmCIF coordinates.
One NormalisedStructure per CIF file; all downstream code operates on this model.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Core atom-level record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AtomRecord:
    """Single atom in a normalised structure."""
    atom_id: int
    atom_name: str
    atom_type: str  # element symbol
    comp_id: str  # residue / monomer name
    chain_id: str  # label_asym_id
    entity_id: int
    seq_id: Optional[int]  # label_seq_id (None if missing)
    auth_seq_id: Optional[int]
    auth_asym_id: str
    coords: Tuple[float, float, float]
    b_factor: Optional[float]
    occupancy: Optional[float]
    model_num: int

    @property
    def xyz(self) -> np.ndarray:
        return np.array(self.coords, dtype=np.float64)


# ---------------------------------------------------------------------------
# Residue-level summary
# ---------------------------------------------------------------------------

_PROTEIN_BACKBONE = frozenset({"N", "CA", "C", "O"})
_PROTEIN_CA = frozenset({"CA"})


@dataclass(frozen=True)
class ResidueRecord:
    """Residue-level summary aggregated from atoms."""
    comp_id: str
    chain_id: str
    entity_id: int
    seq_id: Optional[int]
    auth_seq_id: Optional[int]
    atom_names: Tuple[str, ...]
    n_atoms: int
    is_complete: bool  # all backbone atoms present (proteins only)


# ---------------------------------------------------------------------------
# Entity and chain metadata
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EntityInfo:
    """Entity metadata extracted from mmCIF _entity + _entity_poly."""
    entity_id: int
    chain_id: str  # label_asym_id
    auth_asym_id: str
    entity_type: str  # polymer, non-polymer, water, branched
    polymer_type: Optional[str]  # polypeptide(L), polydeoxyribonucleotide, etc.
    description: Optional[str]
    n_residues: int
    n_atoms: int


@dataclass(frozen=True)
class ChainInfo:
    """Per-chain metadata."""
    chain_id: str
    auth_asym_id: str
    entity_id: int
    entity_type: str
    polymer_type: Optional[str]
    n_residues: int
    n_atoms: int
    first_seq_id: Optional[int]
    last_seq_id: Optional[int]


# ---------------------------------------------------------------------------
# NormalisedStructure
# ---------------------------------------------------------------------------

@dataclass
class NormalisedStructure:
    """
    Complete normalised representation of one AF3 prediction.

    One instance per (condition, seed, sample) CIF file.
    All downstream code (alignment, metrics, comparison) operates on this.
    """
    # Identity
    source_path: Path
    source_checksum: str
    condition_id: str
    seed: int
    sample: int

    # Contents
    atoms: List[AtomRecord] = field(default_factory=list)
    residues: List[ResidueRecord] = field(default_factory=list)
    entities: List[EntityInfo] = field(default_factory=list)
    chains: List[ChainInfo] = field(default_factory=list)

    # Derived indices (built on construction)
    _by_chain: Dict[str, List[AtomRecord]] = field(default_factory=dict)
    _by_entity: Dict[int, List[AtomRecord]] = field(default_factory=dict)
    _by_seq: Dict[Tuple[str, int], List[AtomRecord]] = field(default_factory=dict)
    _by_residue: Dict[Tuple[str, int], List[AtomRecord]] = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Index construction
    # ------------------------------------------------------------------

    def build_indices(self) -> None:
        """Build derived indices from atoms. Call after atoms are populated."""
        self._by_chain = defaultdict(list)
        self._by_entity = defaultdict(list)
        self._by_seq = defaultdict(list)
        self._by_residue = defaultdict(list)

        for atom in self.atoms:
            self._by_chain[atom.chain_id].append(atom)
            self._by_entity[atom.entity_id].append(atom)
            if atom.auth_seq_id is not None:
                key = (atom.chain_id, atom.auth_seq_id)
                self._by_seq[key].append(atom)
                self._by_residue[key].append(atom)

        # Convert defaultdicts to regular dicts for hashability
        self._by_chain = dict(self._by_chain)
        self._by_entity = dict(self._by_entity)
        self._by_seq = dict(self._by_seq)
        self._by_residue = dict(self._by_residue)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def chain_ids(self) -> List[str]:
        return [c.chain_id for c in self.chains]

    @property
    def entity_ids(self) -> List[int]:
        return [e.entity_id for e in self.entities]

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_chains(self) -> int:
        return len(self.chains)

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get_chain(self, chain_id: str) -> List[AtomRecord]:
        """Return all atoms in the given chain."""
        return list(self._by_chain.get(chain_id, []))

    def get_entity(self, entity_id: int) -> List[AtomRecord]:
        """Return all atoms in the given entity."""
        return list(self._by_entity.get(entity_id, []))

    def get_ca_atoms(self, chain_ids: Optional[List[str]] = None) -> List[AtomRecord]:
        """Return C-alpha atoms, optionally restricted to specific chains."""
        result = []
        chains = chain_ids if chain_ids is not None else list(self._by_chain.keys())
        for cid in chains:
            for atom in self._by_chain.get(cid, []):
                if atom.atom_name == "CA":
                    result.append(atom)
        return result

    def get_backbone_atoms(self, chain_ids: Optional[List[str]] = None) -> List[AtomRecord]:
        """Return backbone atoms (N, CA, C, O), optionally restricted to chains."""
        result = []
        chains = chain_ids if chain_ids is not None else list(self._by_chain.keys())
        for cid in chains:
            for atom in self._by_chain.get(cid, []):
                if atom.atom_name in _PROTEIN_BACKBONE:
                    result.append(atom)
        return result

    def get_all_heavy_atoms(self, chain_ids: Optional[List[str]] = None) -> List[AtomRecord]:
        """Return all non-hydrogen atoms, optionally restricted to chains."""
        result = []
        chains = chain_ids if chain_ids is not None else list(self._by_chain.keys())
        for cid in chains:
            for atom in self._by_chain.get(cid, []):
                if atom.atom_type != "H":
                    result.append(atom)
        return result

    def get_residue(self, chain_id: str, auth_seq_id: int) -> List[AtomRecord]:
        """Return all atoms in a specific residue."""
        return list(self._by_residue.get((chain_id, auth_seq_id), []))

    def get_entity_types(self) -> Set[str]:
        """Return set of entity types present."""
        return {e.entity_type for e in self.entities}

    def get_polymer_types(self) -> Set[str]:
        """Return set of polymer types present."""
        return {e.polymer_type for e in self.entities if e.polymer_type is not None}

    def get_protein_chains(self) -> List[str]:
        """Return chain IDs whose entity is polypeptide."""
        return [
            c.chain_id for c in self.chains
            if c.polymer_type is not None and "polypeptide" in c.polymer_type
        ]

    def get_nucleic_acid_chains(self) -> List[str]:
        """Return chain IDs whose entity is nucleic acid."""
        return [
            c.chain_id for c in self.chains
            if c.polymer_type is not None
            and ("polydeoxyribonucleotide" in c.polymer_type or "polyribonucleotide" in c.polymer_type)
        ]
