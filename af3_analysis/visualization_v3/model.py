"""
V3 Normalized Data Model.

Defines the normalized hierarchy for V3 analysis:
Dataset → Condition → Seed → Prediction → StructureData → Entity → Chain → Residue → Atom

This model is constructed by the adapter from existing pipeline outputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Core hierarchy
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AtomPosition:
    """A single atom position with metadata."""
    atom_name: str
    chain_id: str
    entity_id: int
    auth_seq_id: Optional[int]
    coords: Tuple[float, float, float]
    residue_name: str

    @property
    def x(self) -> float:
        return self.coords[0]

    @property
    def y(self) -> float:
        return self.coords[1]

    @property
    def z(self) -> float:
        return self.coords[2]

    def distance_to(self, other: "AtomPosition") -> float:
        """Calculate Euclidean distance to another atom."""
        return np.sqrt(
            (self.x - other.x) ** 2 +
            (self.y - other.y) ** 2 +
            (self.z - other.z) ** 2
        )


@dataclass(frozen=True)
class ResidueGeometry:
    """Geometric data for a single residue."""
    residue_name: str
    chain_id: str
    auth_seq_id: int
    entity_id: int
    # Backbone atoms
    n_coords: Optional[Tuple[float, float, float]] = None
    ca_coords: Optional[Tuple[float, float, float]] = None
    c_coords: Optional[Tuple[float, float, float]] = None
    o_coords: Optional[Tuple[float, float, float]] = None
    # All atoms (for sidechain analysis)
    all_atoms: List[AtomPosition] = field(default_factory=list)
    # Confidence metrics (if available)
    plddt_mean: Optional[float] = None

    @property
    def has_backbone(self) -> bool:
        """Check if backbone atoms (N, CA, C, O) are all present."""
        return all([
            self.n_coords is not None,
            self.ca_coords is not None,
            self.c_coords is not None,
            self.o_coords is not None,
        ])

    @property
    def has_ca(self) -> bool:
        """Check if CA atom is present."""
        return self.ca_coords is not None


@dataclass(frozen=True)
class ChainGeometry:
    """Geometric data for a single chain."""
    chain_id: str
    entity_id: int
    entity_type: str  # polymer, non-polymer, water, etc.
    polymer_type: Optional[str]  # polypeptide(L), polydeoxyribonucleotide, etc.
    residues: List[ResidueGeometry] = field(default_factory=list)
    # Sorted by auth_seq_id
    _residue_index: Dict[int, ResidueGeometry] = field(default_factory=dict)

    def __post_init__(self):
        # Build residue index
        object.__setattr__(self, "_residue_index", {
            r.auth_seq_id: r for r in self.residues if r.auth_seq_id is not None
        })

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    def get_residue(self, auth_seq_id: int) -> Optional[ResidueGeometry]:
        """Get residue by auth_seq_id."""
        return self._residue_index.get(auth_seq_id)

    def get_ca_atoms(self) -> List[AtomPosition]:
        """Get all CA atoms in this chain."""
        return [
            AtomPosition(
                atom_name="CA",
                chain_id=self.chain_id,
                entity_id=r.entity_id,
                auth_seq_id=r.auth_seq_id,
                coords=r.ca_coords,
                residue_name=r.residue_name,
            )
            for r in self.residues
            if r.ca_coords is not None
        ]

    def get_backbone_atoms(self) -> List[AtomPosition]:
        """Get all backbone atoms (N, CA, C, O)."""
        atoms = []
        for r in self.residues:
            if r.n_coords:
                atoms.append(AtomPosition(
                    atom_name="N", chain_id=self.chain_id,
                    entity_id=r.entity_id, auth_seq_id=r.auth_seq_id,
                    coords=r.n_coords, residue_name=r.residue_name,
                ))
            if r.ca_coords:
                atoms.append(AtomPosition(
                    atom_name="CA", chain_id=self.chain_id,
                    entity_id=r.entity_id, auth_seq_id=r.auth_seq_id,
                    coords=r.ca_coords, residue_name=r.residue_name,
                ))
            if r.c_coords:
                atoms.append(AtomPosition(
                    atom_name="C", chain_id=self.chain_id,
                    entity_id=r.entity_id, auth_seq_id=r.auth_seq_id,
                    coords=r.c_coords, residue_name=r.residue_name,
                ))
            if r.o_coords:
                atoms.append(AtomPosition(
                    atom_name="O", chain_id=self.chain_id,
                    entity_id=r.entity_id, auth_seq_id=r.auth_seq_id,
                    coords=r.o_coords, residue_name=r.residue_name,
                ))
        return atoms


@dataclass(frozen=True)
class EntityGeometry:
    """Geometric data for a single entity."""
    entity_id: int
    entity_type: str
    polymer_type: Optional[str]
    description: Optional[str]
    chains: List[ChainGeometry] = field(default_factory=list)
    _chain_index: Dict[str, ChainGeometry] = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "_chain_index", {
            c.chain_id: c for c in self.chains
        })

    def get_chain(self, chain_id: str) -> Optional[ChainGeometry]:
        """Get chain by chain_id."""
        return self._chain_index.get(chain_id)


@dataclass(frozen=True)
class StructureData:
    """
    Normalized structural data for one prediction.

    This is the core data structure for V3 geometric analysis.
    It contains the full hierarchy: Entity → Chain → Residue → Atom.
    """
    prediction_id: str
    condition_id: str
    seed: int
    sample: int
    source_path: Path

    # Confidence metrics (from AF3 output)
    plddt_mean: Optional[float] = None
    plddt_min: Optional[float] = None
    plddt_max: Optional[float] = None
    plddt_median: Optional[float] = None
    pae_mean: Optional[float] = None
    contact_prob_mean: Optional[float] = None

    # Structural data
    entities: List[EntityGeometry] = field(default_factory=list)
    _entity_index: Dict[int, EntityGeometry] = field(default_factory=dict)
    _chain_index: Dict[str, ChainGeometry] = field(default_factory=dict)

    # QC information
    parse_status: str = "success"
    parse_reason: str = ""

    def __post_init__(self):
        # Build indices
        object.__setattr__(self, "_entity_index", {
            e.entity_id: e for e in self.entities
        })
        object.__setattr__(self, "_chain_index", {
            c.chain_id: c for e in self.entities for c in e.chains
        })

    @property
    def chain_ids(self) -> List[str]:
        """Get all chain IDs."""
        return list(self._chain_index.keys())

    @property
    def entity_ids(self) -> List[int]:
        """Get all entity IDs."""
        return list(self._entity_index.keys())

    @property
    def n_chains(self) -> int:
        return len(self._chain_index)

    @property
    def n_entities(self) -> int:
        return len(self._entity_index)

    def get_chain(self, chain_id: str) -> Optional[ChainGeometry]:
        """Get chain by chain_id."""
        return self._chain_index.get(chain_id)

    def get_entity(self, entity_id: int) -> Optional[EntityGeometry]:
        """Get entity by entity_id."""
        return self._entity_index.get(entity_id)

    def get_protein_chains(self) -> List[str]:
        """Get chain IDs for protein/polypeptide chains."""
        return [
            c.chain_id for c in self._chain_index.values()
            if c.polymer_type and "polypeptide" in c.polymer_type
        ]

    def get_nucleic_acid_chains(self) -> List[str]:
        """Get chain IDs for DNA/RNA chains."""
        return [
            c.chain_id for c in self._chain_index.values()
            if c.polymer_type and (
                "polydeoxyribonucleotide" in c.polymer_type or
                "polyribonucleotide" in c.polymer_type
            )
        ]

    @property
    def has_dna(self) -> bool:
        """Whether the structure contains a DNA/RNA chain."""
        return bool(self.get_nucleic_acid_chains())

    def get_ca_atoms(self, chain_ids: Optional[List[str]] = None) -> List[AtomPosition]:
        """Get all CA atoms, optionally filtered by chain."""
        chains = chain_ids if chain_ids is not None else list(self._chain_index.keys())
        atoms = []
        for cid in chains:
            chain = self._chain_index.get(cid)
            if chain:
                atoms.extend(chain.get_ca_atoms())
        return atoms

    def get_backbone_atoms(self, chain_ids: Optional[List[str]] = None) -> List[AtomPosition]:
        """Get all backbone atoms, optionally filtered by chain."""
        chains = chain_ids if chain_ids is not None else list(self._chain_index.keys())
        atoms = []
        for cid in chains:
            chain = self._chain_index.get(cid)
            if chain:
                atoms.extend(chain.get_backbone_atoms())
        return atoms

    def get_residue(self, chain_id: str, auth_seq_id: int) -> Optional[ResidueGeometry]:
        """Get a specific residue."""
        chain = self._chain_index.get(chain_id)
        if chain:
            return chain.get_residue(auth_seq_id)
        return None

    def get_common_residues(self, other: "StructureData", chain_id: Optional[str] = None) -> List[Tuple[int, int]]:
        """
        Find common residues between this structure and another.

        Returns list of (auth_seq_id_this, auth_seq_id_other) pairs.
        Assumes same chain_id if specified.
        """
        common = []
        if chain_id:
            chain_a = self._chain_index.get(chain_id)
            chain_b = other._chain_index.get(chain_id)
            if chain_a and chain_b:
                ids_a = {r.auth_seq_id for r in chain_a.residues if r.auth_seq_id is not None}
                ids_b = {r.auth_seq_id for r in chain_b.residues if r.auth_seq_id is not None}
                common_ids = ids_a & ids_b
                common = [(cid, cid) for cid in sorted(common_ids)]
        else:
            # Compare across all chains - more complex, skip for now
            pass
        return common

    def get_common_ca_atoms(self, other: "StructureData", chain_id: Optional[str] = None) -> List[Tuple[AtomPosition, AtomPosition]]:
        """Get paired CA atoms for common residues."""
        pairs = []
        if chain_id:
            chain_a = self._chain_index.get(chain_id)
            chain_b = other._chain_index.get(chain_id)
            if chain_a and chain_b:
                for res_a in chain_a.residues:
                    if res_a.ca_coords is None or res_a.auth_seq_id is None:
                        continue
                    res_b = chain_b.get_residue(res_a.auth_seq_id)
                    if res_b and res_b.ca_coords is not None:
                        atoms_a = AtomPosition(
                            atom_name="CA", chain_id=chain_id,
                            entity_id=res_a.entity_id, auth_seq_id=res_a.auth_seq_id,
                            coords=res_a.ca_coords, residue_name=res_a.residue_name,
                        )
                        atoms_b = AtomPosition(
                            atom_name="CA", chain_id=chain_id,
                            entity_id=res_b.entity_id, auth_seq_id=res_b.auth_seq_id,
                            coords=res_b.ca_coords, residue_name=res_b.residue_name,
                        )
                        pairs.append((atoms_a, atoms_b))
        return pairs


@dataclass(frozen=True)
class MetricObservation:
    """A single metric observation at a specific scope."""
    metric_id: str
    value: Optional[float]
    scope_type: str  # global, chain, interface, region
    scope_id: str = ""  # chain_id, interface_id, region_label
    condition_id: str = ""
    seed: int = 0
    sample: int = 0
    prediction_id: str = ""


@dataclass(frozen=True)
class ConditionSummary:
    """Summary statistics for one condition."""
    condition_id: str
    condition_name: str
    n_seeds: int
    n_predictions: int
    metrics: Dict[str, Dict[str, float]] = field(default_factory=dict)
    # metric_id -> {mean, std, median, min, max, n}


@dataclass(frozen=True)
class SeedSummary:
    """Summary statistics for one seed within a condition."""
    condition_id: str
    seed: int
    n_samples: int
    metrics: Dict[str, float] = field(default_factory=dict)
    # metric_id -> mean value across samples


@dataclass(frozen=True)
class Dataset:
    """
    Complete V3 dataset.

    Contains all structures, organized by condition → seed → prediction.
    """
    name: str
    conditions: Dict[str, ConditionSummary] = field(default_factory=dict)
    seeds: Dict[str, Dict[int, SeedSummary]] = field(default_factory=dict)
    # condition_id -> seed -> SeedSummary
    predictions: Dict[str, Dict[int, Dict[int, StructureData]]] = field(default_factory=dict)
    # condition_id -> seed -> sample -> StructureData

    # Reference resolution
    reference_condition: Optional[str] = None

    # Metadata
    experiment_metadata: Optional[Dict[str, Any]] = None

    def get_structure(self, condition_id: str, seed: int, sample: int) -> Optional[StructureData]:
        """Get a specific structure."""
        cond_seeds = self.predictions.get(condition_id, {})
        seed_preds = cond_seeds.get(seed, {})
        return seed_preds.get(sample)

    def get_structures_for_condition(self, condition_id: str) -> List[StructureData]:
        """Get all structures for a condition."""
        result = []
        cond_seeds = self.predictions.get(condition_id, {})
        for seed_preds in cond_seeds.values():
            for structure in seed_preds.values():
                result.append(structure)
        return result

    def get_seeds_for_condition(self, condition_id: str) -> List[int]:
        """Get all seeds for a condition."""
        cond_seeds = self.predictions.get(condition_id, {})
        return sorted(cond_seeds.keys())

    def get_common_seeds(self, condition_a: str, condition_b: str) -> List[int]:
        """Get seeds common to both conditions."""
        seeds_a = set(self.get_seeds_for_condition(condition_a))
        seeds_b = set(self.get_seeds_for_condition(condition_b))
        return sorted(seeds_a & seeds_b)


# ---------------------------------------------------------------------------
# QC records
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StructuralQC:
    """QC record for one structure."""
    prediction_id: str
    condition_id: str
    seed: int
    sample: int
    status: str  # success, parse_error, empty, etc.
    reason: str = ""
    n_atoms: int = 0
    n_chains: int = 0
    n_residues: int = 0
    entity_types: List[str] = field(default_factory=list)
    chain_ids: List[str] = field(default_factory=list)
    has_protein: bool = False
    has_dna: bool = False
    has_rna: bool = False
    has_ligand: bool = False
    has_ion: bool = False


@dataclass(frozen=True)
class ComparisonQC:
    """QC record for one structural comparison."""
    comparison_id: str
    condition_a: str
    condition_b: str
    seed_a: int
    seed_b: int
    status: str  # comparable, partially_comparable, not_comparable
    reason: str = ""
    n_common_atoms: int = 0
    n_common_residues: int = 0
    coverage: float = 0.0  # fraction of reference residues covered
    low_coverage: bool = False
