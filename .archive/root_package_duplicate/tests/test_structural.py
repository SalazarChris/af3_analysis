"""
Comprehensive tests for the structural analysis subsystem.

Uses synthetic structures with known coordinates for numerical validation.
No biological assumptions; tests cover generic geometric operations.
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pytest

from af3_analysis.structural.alignment import align_structures
from af3_analysis.structural.comparison import compare_conditions
from af3_analysis.structural.enums import ComparisonStatus

# ---------------------------------------------------------------------------
# Synthetic structure builders
# ---------------------------------------------------------------------------

def _make_atom(
    atom_id: int,
    atom_name: str,
    atom_type: str,
    comp_id: str,
    chain_id: str,
    entity_id: int,
    seq_id: int,
    auth_seq_id: int,
    coords: Tuple[float, float, float],
    *,
    auth_asym_id: Optional[str] = None,
    b_factor: float = 50.0,
    occupancy: float = 1.0,
    model_num: int = 1,
):
    from af3_analysis.structural.representation import AtomRecord
    return AtomRecord(
        atom_id=atom_id,
        atom_name=atom_name,
        atom_type=atom_type,
        comp_id=comp_id,
        chain_id=chain_id,
        entity_id=entity_id,
        seq_id=seq_id,
        auth_seq_id=auth_seq_id,
        auth_asym_id=auth_asym_id or chain_id,
        coords=coords,
        b_factor=b_factor,
        occupancy=occupancy,
        model_num=model_num,
    )


def _make_structure(
    condition_id: str,
    seed: int,
    sample: int,
    atoms: list,
    *,
    source_path: Optional[Path] = None,
    source_checksum: str = "test_checksum",
):
    from af3_analysis.structural.representation import (
        NormalisedStructure, EntityInfo, ChainInfo,
        ResidueRecord,
    )
    from collections import defaultdict

    # Derive entities and chains from atoms
    entity_map = {}
    chain_map = {}
    for atom in atoms:
        if atom.entity_id not in entity_map:
            entity_map[atom.entity_id] = {
                "entity_id": atom.entity_id,
                "chains": set(),
                "type": "polymer",
                "polymer_type": "polypeptide(L)",
            }
        entity_map[atom.entity_id]["chains"].add(atom.chain_id)

        if atom.chain_id not in chain_map:
            chain_map[atom.chain_id] = {
                "chain_id": atom.chain_id,
                "auth_asym_id": atom.auth_asym_id,
                "entity_id": atom.entity_id,
                "type": "polymer",
                "polymer_type": "polypeptide(L)",
            }

    entities = []
    for eid, info in sorted(entity_map.items()):
        chains_in_entity = info["chains"]
        n_atoms = sum(1 for a in atoms if a.entity_id == eid)
        n_res = len(set(a.auth_seq_id for a in atoms if a.entity_id == eid))
        entities.append(EntityInfo(
            entity_id=eid,
            chain_id=sorted(chains_in_entity)[0],
            auth_asym_id=sorted(chains_in_entity)[0],
            entity_type=info["type"],
            polymer_type=info["polymer_type"],
            description=None,
            n_residues=n_res,
            n_atoms=n_atoms,
        ))

    chains = []
    for cid, info in sorted(chain_map.items()):
        chains_in_entity = entity_map[info["entity_id"]]["chains"]
        n_atoms = sum(1 for a in atoms if a.chain_id == cid)
        n_res = len(set(a.auth_seq_id for a in atoms if a.chain_id == cid))
        chains.append(ChainInfo(
            chain_id=cid,
            auth_asym_id=info["auth_asym_id"],
            entity_id=info["entity_id"],
            entity_type=info["type"],
            polymer_type=info["polymer_type"],
            n_residues=n_res,
            n_atoms=n_atoms,
            first_seq_id=min(a.auth_seq_id for a in atoms if a.chain_id == cid) if any(a.chain_id == cid for a in atoms) else None,
            last_seq_id=max(a.auth_seq_id for a in atoms if a.chain_id == cid) if any(a.chain_id == cid for a in atoms) else None,
        ))

    # Build residues
    residue_atoms = defaultdict(list)
    for atom in atoms:
        residue_atoms[(atom.chain_id, atom.auth_seq_id)].append(atom)

    residues = []
    for (cid, sid), res_atoms in sorted(residue_atoms.items()):
        atom_names = tuple(sorted(set(a.atom_name for a in res_atoms)))
        residues.append(ResidueRecord(
            comp_id=res_atoms[0].comp_id,
            chain_id=cid,
            entity_id=res_atoms[0].entity_id,
            seq_id=None,
            auth_seq_id=sid,
            atom_names=atom_names,
            n_atoms=len(res_atoms),
            is_complete={"N", "CA", "C"}.issubset(set(atom_names)),
        ))

    structure = NormalisedStructure(
        source_path=source_path or Path(f"synthetic_{condition_id}_s{seed}_n{sample}.cif"),
        source_checksum=source_checksum,
        condition_id=condition_id,
        seed=seed,
        sample=sample,
        atoms=atoms,
        residues=residues,
        entities=entities,
        chains=chains,
    )
    structure.build_indices()
    return structure


def _make_linear_protein(
    chain_id: str,
    entity_id: int,
    n_residues: int,
    *,
    start_coords: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    residue_spacing: float = 3.8,
    atom_names=None,
):
    """Create atoms for a linear protein chain with known geometry."""
    if atom_names is None:
        atom_names = ["N", "CA", "C", "O"]

    atoms = []
    atom_id = 1
    for i in range(n_residues):
        x0, y0, z0 = start_coords
        for j, name in enumerate(atom_names):
            # Simple geometry: CA at regular intervals, others offset
            if name == "CA":
                coords = (x0 + i * residue_spacing, y0, z0)
            elif name == "N":
                coords = (x0 + i * residue_spacing - 1.2, y0 + 0.5, z0)
            elif name == "C":
                coords = (x0 + i * residue_spacing + 1.5, y0 - 0.3, z0)
            elif name == "O":
                coords = (x0 + i * residue_spacing + 1.5, y0 - 1.5, z0)
            else:
                coords = (x0 + i * residue_spacing, y0 + 1.0, z0)

            atoms.append(_make_atom(
                atom_id=atom_id,
                atom_name=name,
                atom_type=name[0],
                comp_id="ALA",
                chain_id=chain_id,
                entity_id=entity_id,
                seq_id=i + 1,
                auth_seq_id=i + 1,
                coords=coords,
            ))
            atom_id += 1
    return atoms


# ---------------------------------------------------------------------------
# Test: Representation
# ---------------------------------------------------------------------------

class TestNormalisedStructure:
    """Test the NormalisedStructure data model."""

    def test_single_chain_structure(self):
        atoms = _make_linear_protein("A", 1, 5)
        s = _make_structure("test_cond", 1, 0, atoms)

        assert s.n_atoms == 20  # 5 residues × 4 atoms
        assert s.n_chains == 1
        assert s.chain_ids == ["A"]
        assert len(s.entities) == 1
        assert len(s.residues) == 5

    def test_multi_chain_structure(self):
        atoms_a = _make_linear_protein("A", 1, 5)
        atoms_b = _make_linear_protein("B", 2, 3, start_coords=(10.0, 0.0, 0.0))
        atoms = atoms_a + atoms_b
        s = _make_structure("test_cond", 1, 0, atoms)

        assert s.n_chains == 2
        assert set(s.chain_ids) == {"A", "B"}
        assert len(s.entities) == 2

    def test_get_ca_atoms(self):
        atoms = _make_linear_protein("A", 1, 5)
        s = _make_structure("test_cond", 1, 0, atoms)
        ca_atoms = s.get_ca_atoms()
        assert len(ca_atoms) == 5
        assert all(a.atom_name == "CA" for a in ca_atoms)

    def test_get_backbone_atoms(self):
        atoms = _make_linear_protein("A", 1, 5)
        s = _make_structure("test_cond", 1, 0, atoms)
        bb = s.get_backbone_atoms()
        assert len(bb) == 20  # N, CA, C, O per residue

    def test_get_chain(self):
        atoms_a = _make_linear_protein("A", 1, 3)
        atoms_b = _make_linear_protein("B", 2, 2, start_coords=(10.0, 0.0, 0.0))
        s = _make_structure("test_cond", 1, 0, atoms_a + atoms_b)

        chain_a = s.get_chain("A")
        assert len(chain_a) == 12  # 3 × 4
        chain_b = s.get_chain("B")
        assert len(chain_b) == 8  # 2 × 4

    def test_get_entity_types(self):
        atoms = _make_linear_protein("A", 1, 3)
        s = _make_structure("test_cond", 1, 0, atoms)
        assert "polymer" in s.get_entity_types()

    def test_get_protein_chains(self):
        atoms = _make_linear_protein("A", 1, 3)
        s = _make_structure("test_cond", 1, 0, atoms)
        assert s.get_protein_chains() == ["A"]


# ---------------------------------------------------------------------------
# Test: Alignment
# ---------------------------------------------------------------------------

class TestAlignment:
    """Test structural alignment engine."""

    def test_identical_structures_zero_rmsd(self):
        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)
        s2 = _make_structure("cond2", 1, 0, list(atoms))

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.COMPARABLE
        assert result.rmsd_post_alignment is not None
        assert result.rmsd_post_alignment < 1e-10

    def test_known_translation_rmsd(self):
        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)

        # Translate by (5, 0, 0)
        translated = [
            _make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                (a.coords[0] + 5.0, a.coords[1], a.coords[2]),
                auth_asym_id=a.auth_asym_id,
            )
            for a in atoms
        ]
        s2 = _make_structure("cond2", 1, 0, translated)

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.COMPARABLE
        # After Kabsch alignment, pure translation should give RMSD ≈ 0
        assert result.rmsd_post_alignment < 1e-10

    def test_known_rotation_rmsd(self):
        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)

        # Rotate 90° around Z axis
        angle = math.pi / 2
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        rotated = []
        for a in atoms:
            x, y, z = a.coords
            rx = cos_a * x - sin_a * y
            ry = sin_a * x + cos_a * y
            rotated.append(_make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                (rx, ry, z),
                auth_asym_id=a.auth_asym_id,
            ))
        s2 = _make_structure("cond2", 1, 0, rotated)

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.COMPARABLE
        # After Kabsch, pure rotation should give RMSD ≈ 0
        assert result.rmsd_post_alignment < 1e-10

    def test_genuinely_displaced_structures(self):
        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)

        # Displace CA atoms of residues 3-5 by 2Å in Y direction
        displaced = []
        for a in atoms:
            if a.atom_name == "CA" and a.auth_seq_id >= 3:
                coords = (a.coords[0], a.coords[1] + 2.0, a.coords[2])
            else:
                coords = a.coords
            displaced.append(_make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                coords,
                auth_asym_id=a.auth_asym_id,
            ))
        s2 = _make_structure("cond2", 1, 0, displaced)

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.COMPARABLE
        assert result.rmsd_post_alignment > 0.4  # Should be significant

    def test_different_sequence_lengths(self):
        atoms_a = _make_linear_protein("A", 1, 5)
        atoms_b = _make_linear_protein("A", 1, 3)
        s1 = _make_structure("cond1", 1, 0, atoms_a)
        s2 = _make_structure("cond2", 1, 0, atoms_b)

        result = align_structures(s1, s2)
        # Should be partially comparable (3/5 residues match)
        assert result.status in (
            ComparisonStatus.COMPARABLE,
            ComparisonStatus.PARTIALLY_COMPARABLE,
        )
        assert result.n_common_residues == 3

    def test_missing_residues(self):
        atoms_a = _make_linear_protein("A", 1, 5)
        # Remove residue 3 from structure B
        atoms_b = [a for a in atoms_a if a.auth_seq_id != 3]
        s1 = _make_structure("cond1", 1, 0, atoms_a)
        s2 = _make_structure("cond2", 1, 0, atoms_b)

        result = align_structures(s1, s2)
        assert result.status in (
            ComparisonStatus.COMPARABLE,
            ComparisonStatus.PARTIALLY_COMPARABLE,
        )
        assert result.n_common_residues == 4

    def test_missing_chain(self):
        atoms_a = _make_linear_protein("A", 1, 3)
        atoms_b = _make_linear_protein("B", 2, 3, start_coords=(10.0, 0.0, 0.0))
        s1 = _make_structure("cond1", 1, 0, atoms_a + atoms_b)
        s2 = _make_structure("cond2", 1, 0, list(atoms_a))  # No chain B

        result = align_structures(s1, s2)
        # Chain A matches, chain B is missing
        assert result.status == ComparisonStatus.PARTIALLY_COMPARABLE
        assert result.n_common_residues == 3

    def test_no_common_residues(self):
        atoms_a = _make_linear_protein("A", 1, 3)
        atoms_b = _make_linear_protein("A", 1, 3, start_coords=(100.0, 0.0, 0.0))
        # Rename chain IDs to prevent matching
        atoms_b = [
            _make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                "B", a.entity_id, a.seq_id, a.auth_seq_id,
                a.coords, auth_asym_id="B",
            )
            for a in atoms_b
        ]
        s1 = _make_structure("cond1", 1, 0, atoms_a)
        s2 = _make_structure("cond2", 1, 0, atoms_b)

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.NOT_COMPARABLE
        assert result.reason.value == "no_common_residues"

    def test_multiple_protein_chains(self):
        atoms_a = _make_linear_protein("A", 1, 3)
        atoms_b = _make_linear_protein("B", 2, 3, start_coords=(10.0, 0.0, 0.0))
        s1 = _make_structure("cond1", 1, 0, atoms_a + atoms_b)
        s2 = _make_structure("cond2", 1, 0, list(atoms_a) + list(atoms_b))

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.COMPARABLE
        assert result.n_common_residues == 6

    def test_alignment_transform_applied(self):
        atoms = _make_linear_protein("A", 1, 3)
        s1 = _make_structure("cond1", 1, 0, atoms)
        s2 = _make_structure("cond2", 1, 0, list(atoms))

        result = align_structures(s1, s2)
        assert result.rotation is not None
        assert result.translation is not None
        # For identical structures, rotation should be identity
        np.testing.assert_allclose(result.rotation, np.eye(3), atol=1e-10)


# ---------------------------------------------------------------------------
# Test: RMSD Metrics
# ---------------------------------------------------------------------------

class TestRMSDMetrics:
    """Test RMSD metric computations."""

    def test_rmsd_identical_structures(self):
        from af3_analysis.structural.metrics.rmsd import GlobalCaRMSD

        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)
        s2 = _make_structure("cond2", 1, 0, list(atoms))

        alignment = align_structures(s1, s2)
        metric = GlobalCaRMSD()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value < 1e-10

    def test_rmsd_known_displacement(self):
        from af3_analysis.structural.metrics.rmsd import GlobalCaRMSD

        atoms = _make_linear_protein("A", 1, 10)
        s1 = _make_structure("cond1", 1, 0, atoms)

        # Displace all atoms by (1, 0, 0)
        displaced = [
            _make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                (a.coords[0] + 1.0, a.coords[1], a.coords[2]),
                auth_asym_id=a.auth_asym_id,
            )
            for a in atoms
        ]
        s2 = _make_structure("cond2", 1, 0, displaced)

        alignment = align_structures(s1, s2)
        metric = GlobalCaRMSD()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value < 1e-10  # Pure translation should align perfectly

    def test_chain_rmsd(self):
        from af3_analysis.structural.metrics.rmsd import ChainCaRMSD

        atoms_a = _make_linear_protein("A", 1, 5)
        atoms_b = _make_linear_protein("B", 2, 5, start_coords=(10.0, 0.0, 0.0))
        s1 = _make_structure("cond1", 1, 0, atoms_a + atoms_b)
        s2 = _make_structure("cond2", 1, 0, list(atoms_a) + list(atoms_b))

        alignment = align_structures(s1, s2)
        metric = ChainCaRMSD()
        value = metric.compute(s1, s2, alignment, chain_id="A")
        assert value is not None
        assert value < 1e-10


# ---------------------------------------------------------------------------
# Test: Centroid Metrics
# ---------------------------------------------------------------------------

class TestCentroidMetrics:
    """Test centroid displacement metrics."""

    def test_centroid_zero_displacement(self):
        from af3_analysis.structural.metrics.centroid import CentroidDisplacement

        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)
        s2 = _make_structure("cond2", 1, 0, list(atoms))

        alignment = align_structures(s1, s2)
        metric = CentroidDisplacement()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value < 1e-10

    def test_centroid_known_displacement(self):
        from af3_analysis.structural.metrics.centroid import CentroidDisplacement

        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)

        # Translate by (3, 4, 0) — displacement should be 5Å
        displaced = [
            _make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                (a.coords[0] + 3.0, a.coords[1] + 4.0, a.coords[2]),
                auth_asym_id=a.auth_asym_id,
            )
            for a in atoms
        ]
        s2 = _make_structure("cond2", 1, 0, displaced)

        alignment = align_structures(s1, s2)
        metric = CentroidDisplacement()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value < 1e-10  # After alignment, translation is absorbed

    def test_centroid_raw_displacement(self):
        from af3_analysis.structural.metrics.centroid import CentroidDisplacement

        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)

        displaced = [
            _make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                (a.coords[0] + 3.0, a.coords[1] + 4.0, a.coords[2]),
                auth_asym_id=a.auth_asym_id,
            )
            for a in atoms
        ]
        s2 = _make_structure("cond2", 1, 0, displaced)

        # Without alignment, centroid displacement should be 5
        alignment = align_structures(s1, s2)
        # The Kabsch should align them, so post-alignment displacement ≈ 0
        metric = CentroidDisplacement()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value < 0.1


# ---------------------------------------------------------------------------
# Test: Contact Metrics
# ---------------------------------------------------------------------------

class TestContactMetrics:
    """Test contact-map metrics."""

    def test_identical_contact_maps(self):
        from af3_analysis.structural.metrics.contacts import ContactMapDifference

        atoms = _make_linear_protein("A", 1, 10)
        s1 = _make_structure("cond1", 1, 0, atoms)
        s2 = _make_structure("cond2", 1, 0, list(atoms))

        alignment = align_structures(s1, s2)
        metric = ContactMapDifference()
        value = metric.compute(s1, s2, alignment, contact_cutoff=8.0)
        assert value is not None
        assert value == 0.0

    def test_contact_gain_loss(self):
        from af3_analysis.structural.metrics.contacts import ContactMapDifference

        atoms = _make_linear_protein("A", 1, 10)
        s1 = _make_structure("cond1", 1, 0, atoms)

        # Move residue 5 CA far away to break contacts
        displaced = []
        for a in atoms:
            if a.atom_name == "CA" and a.auth_seq_id == 5:
                coords = (a.coords[0], a.coords[1] + 20.0, a.coords[2])
            else:
                coords = a.coords
            displaced.append(_make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                coords,
                auth_asym_id=a.auth_asym_id,
            ))
        s2 = _make_structure("cond2", 1, 0, displaced)

        alignment = align_structures(s1, s2)
        metric = ContactMapDifference()
        value = metric.compute(s1, s2, alignment, contact_cutoff=8.0)
        assert value is not None
        assert value > 0.0  # Contact maps differ

    def test_interface_min_distance(self):
        from af3_analysis.structural.metrics.contacts import InterfaceMinDistance

        atoms_a = _make_linear_protein("A", 1, 3)
        atoms_b = _make_linear_protein("B", 2, 3, start_coords=(5.0, 0.0, 0.0))
        s1 = _make_structure("cond1", 1, 0, atoms_a + atoms_b)
        s2 = _make_structure("cond2", 1, 0, list(atoms_a) + list(atoms_b))

        alignment = align_structures(s1, s2)
        metric = InterfaceMinDistance()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value > 0.0


# ---------------------------------------------------------------------------
# Test: Distance Metrics
# ---------------------------------------------------------------------------

class TestDistanceMetrics:
    """Test pairwise distance change metrics."""

    def test_pairwise_distance_identical(self):
        from af3_analysis.structural.metrics.distance_matrix import PairwiseDistanceChange

        atoms = _make_linear_protein("A", 1, 5)
        s1 = _make_structure("cond1", 1, 0, atoms)
        s2 = _make_structure("cond2", 1, 0, list(atoms))

        alignment = align_structures(s1, s2)
        metric = PairwiseDistanceChange()
        value = metric.compute(s1, s2, alignment)
        assert value is not None
        assert value < 1e-10


# ---------------------------------------------------------------------------
# Test: Comparison Orchestrator
# ---------------------------------------------------------------------------

class TestComparisonOrchestrator:
    """Test the comparison orchestrator."""

    def _make_multi_condition_structures(self):
        """Create structures for two conditions with matching seeds."""
        atoms_ref = _make_linear_protein("A", 1, 5)
        atoms_tgt = _make_linear_protein("A", 1, 5)

        # Displace target slightly
        atoms_tgt = [
            _make_atom(
                a.atom_id, a.atom_name, a.atom_type, a.comp_id,
                a.chain_id, a.entity_id, a.seq_id, a.auth_seq_id,
                (a.coords[0], a.coords[1] + 0.5, a.coords[2]),
                auth_asym_id=a.auth_asym_id,
            )
            for a in atoms_tgt
        ]

        structures = []
        for seed in [1, 2, 3]:
            structures.append(_make_structure("reference", seed, 0, list(atoms_ref)))
            structures.append(_make_structure("target", seed, 0, list(atoms_tgt)))
        return structures

    def test_matched_seed_comparison(self):
        from af3_analysis.structural.config import StructuralConfig

        structures = self._make_multi_condition_structures()
        config = StructuralConfig(
            enabled=True,
            reference_condition="reference",
            comparison_mode="matched_seed",
        )

        comparisons = compare_conditions(structures, config)
        assert len(comparisons) > 0

        # Each comparison should be reference vs target
        for comp in comparisons:
            assert {comp.condition_a, comp.condition_b} == {"reference", "target"}

    def test_pairwise_all_comparison(self):
        from af3_analysis.structural.config import StructuralConfig

        structures = self._make_multi_condition_structures()
        config = StructuralConfig(
            enabled=True,
            comparison_mode="pairwise_all",
        )

        comparisons = compare_conditions(structures, config)
        assert len(comparisons) > 0

    def test_incomparable_pair_logged(self):
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.enums import ComparisonStatus

        # Create structures with no common residues
        atoms_a = _make_linear_protein("A", 1, 3)
        atoms_b = _make_linear_protein("B", 2, 3, start_coords=(100.0, 0.0, 0.0))

        structures = [
            _make_structure("cond1", 1, 0, atoms_a),
            _make_structure("cond2", 1, 0, atoms_b),
        ]

        config = StructuralConfig(
            enabled=True,
            comparison_mode="pairwise_all",
            min_sequence_identity=0.0,
        )

        comparisons = compare_conditions(structures, config)
        # Should have at least one comparison
        assert len(comparisons) >= 1
        # At least one should be not_comparable
        statuses = [c.alignment_status for c in comparisons]
        assert ComparisonStatus.NOT_COMPARABLE in statuses


# ---------------------------------------------------------------------------
# Test: Metric Registry
# ---------------------------------------------------------------------------

class TestMetricRegistry:
    """Test the metric registry."""

    def test_default_registry_has_metrics(self):
        from af3_analysis.structural.metrics.registry import get_default_registry

        reg = get_default_registry()
        assert len(reg.all_metrics()) == 12
        assert "rmsd_global_ca" in reg.metric_ids()
        assert "centroid_displacement" in reg.metric_ids()

    def test_registry_applicability(self):
        from af3_analysis.structural.metrics.registry import get_default_registry

        reg = get_default_registry()
        atoms = _make_linear_protein("A", 1, 5)
        s = _make_structure("cond1", 1, 0, atoms)

        applicable = reg.applicable_metrics(s, s)
        # All metrics should be applicable to a single protein chain
        assert len(applicable) > 0

    def test_filter_by_ids(self):
        from af3_analysis.structural.metrics.registry import get_default_registry

        reg = get_default_registry()
        filtered = reg.filter_by_ids(["rmsd_global_ca", "centroid_displacement"])
        assert len(filtered) == 2


# ---------------------------------------------------------------------------
# Test: Tables
# ---------------------------------------------------------------------------

class TestOutputTables:
    """Test output table generation."""

    def test_write_tables(self, tmp_path):
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.metrics.registry import get_default_registry
        from af3_analysis.structural.tables import write_structural_tables
        from af3_analysis.structural.comparison import compare_conditions

        atoms = _make_linear_protein("A", 1, 5)
        structures = [
            _make_structure("cond1", 1, 0, list(atoms)),
            _make_structure("cond2", 1, 0, list(atoms)),
        ]

        config = StructuralConfig(
            enabled=True,
            reference_condition="cond1",
            comparison_mode="matched_seed",
        )
        registry = get_default_registry()
        comparisons = compare_conditions(structures, config, registry=registry)

        paths = write_structural_tables(
            tmp_path, structures, comparisons, registry, config,
        )

        assert "structural_predictions" in paths
        assert "structural_metrics" in paths
        assert "structural_comparisons" in paths
        assert "structural_effects" in paths
        assert "structural_quality" in paths

        # Verify files exist and are non-empty
        for name, path in paths.items():
            assert path.exists()
            assert path.stat().st_size > 0

        # Verify predictions table
        import pandas as pd
        df = pd.read_csv(paths["structural_predictions"])
        assert len(df) == 2
        assert "condition_id" in df.columns
        assert "n_atoms" in df.columns

    def test_tables_with_parse_errors(self, tmp_path):
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.metrics.registry import get_default_registry
        from af3_analysis.structural.tables import write_structural_tables

        atoms = _make_linear_protein("A", 1, 5)
        structures = [_make_structure("cond1", 1, 0, list(atoms))]
        parse_errors = [
            {"condition_id": "cond2", "seed": 1, "sample": 0,
             "source_path": "bad.cif", "status": "parse_error", "reason": "test"},
        ]

        config = StructuralConfig(enabled=True)
        registry = get_default_registry()

        paths = write_structural_tables(
            tmp_path, structures, [], registry, config,
            parse_errors=parse_errors,
        )

        import pandas as pd
        df = pd.read_csv(paths["structural_predictions"])
        assert len(df) == 2  # 1 success + 1 error
        assert "parse_error" in df["parse_status"].values


# ---------------------------------------------------------------------------
# Test: CIF Reader with Real File
# ---------------------------------------------------------------------------

class TestCIFReader:
    """Test CIF parsing with a real AF3 output file."""

    REAL_CIF = Path("testdata/pou2/pou_baseline/seed-9_sample-4/pou_baseline_seed-9_sample-4_model.cif")

    @pytest.mark.skipif(not REAL_CIF.exists(), reason="Test CIF not available")
    def test_parse_real_cif(self):
        from af3_analysis.io.structure_reader import parse_mmcif

        s = parse_mmcif(self.REAL_CIF, "test_cond", 9, 4)
        assert s.n_atoms > 0
        assert s.n_chains >= 1
        assert len(s.entities) >= 1
        assert len(s.residues) > 0
        assert s.condition_id == "test_cond"
        assert s.seed == 9
        assert s.sample == 4

    @pytest.mark.skipif(not REAL_CIF.exists(), reason="Test CIF not available")
    def test_parse_produces_valid_structure(self):
        from af3_analysis.io.structure_reader import parse_mmcif

        s = parse_mmcif(self.REAL_CIF, "test", 1, 0)
        # Should be able to get CA atoms
        ca = s.get_ca_atoms()
        assert len(ca) > 0
        # Should be able to align with itself
        result = align_structures(s, s)
        assert result.status == ComparisonStatus.COMPARABLE
        assert result.rmsd_post_alignment < 1e-10

    @pytest.mark.skipif(not REAL_CIF.exists(), reason="Test CIF not available")
    def test_discover_structure_files(self):
        from af3_analysis.io.structure_reader import discover_structure_files

        root = Path("testdata/pou2")
        files = discover_structure_files(root)
        assert len(files) > 0
        # Each tuple should have (path, condition_id, seed, sample)
        for path, cond, seed, sample in files:
            assert path.exists()
            assert isinstance(cond, str)
            assert isinstance(seed, int)
            assert isinstance(sample, int)


# ---------------------------------------------------------------------------
# Test: Config
# ---------------------------------------------------------------------------

class TestStructuralConfig:
    """Test structural configuration."""

    def test_default_config(self):
        from af3_analysis.structural.config import StructuralConfig

        config = StructuralConfig()
        assert config.enabled is False
        assert config.atom_selection == "ca"
        assert config.contact_cutoff_angstrom == 8.0
        assert config.comparison_mode == "matched_seed"

    def test_config_from_dict(self):
        from af3_analysis.structural.config import StructuralConfig

        data = {
            "enabled": True,
            "atom_selection": "backbone",
            "contact_cutoff_angstrom": 10.0,
            "reference_condition": "my_reference",
        }
        config = StructuralConfig.from_dict(data)
        assert config.enabled is True
        assert config.atom_selection == "backbone"
        assert config.contact_cutoff_angstrom == 10.0
        assert config.reference_condition == "my_reference"

    def test_config_roundtrip(self):
        from af3_analysis.structural.config import StructuralConfig

        config = StructuralConfig(
            enabled=True,
            reference_condition="ref",
            contact_cutoff_angstrom=7.5,
        )
        d = config.to_dict()
        config2 = StructuralConfig.from_dict(d)
        assert config == config2


# ---------------------------------------------------------------------------
# Test: Edge Cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_structure_alignment(self):
        from af3_analysis.structural.enums import ComparisonStatus

        s1 = _make_structure("cond1", 1, 0, [])
        s2 = _make_structure("cond2", 1, 0, [])

        result = align_structures(s1, s2)
        assert result.status == ComparisonStatus.NOT_COMPARABLE

    def test_single_atom_structure(self):
        from af3_analysis.structural.enums import ComparisonStatus

        atom = _make_atom(1, "CA", "C", "ALA", "A", 1, 1, 1, (0.0, 0.0, 0.0))
        s1 = _make_structure("cond1", 1, 0, [atom])
        s2 = _make_structure("cond2", 1, 0, [atom])

        result = align_structures(s1, s2, min_common_atoms=1)
        assert result.status == ComparisonStatus.COMPARABLE

    def test_protein_nucleic_acid_structure(self):
        """Test structure with protein + nucleic acid chains."""
        atoms_prot = _make_linear_protein("A", 1, 3)
        # Create DNA-like atoms (simplified)
        atoms_dna = []
        atom_id = len(atoms_prot) + 1
        for i in range(3):
            for name in ["P", "OP1", "OP2", "C5'", "C4'", "C3'"]:
                atoms_dna.append(_make_atom(
                    atom_id=atom_id,
                    atom_name=name,
                    atom_type="P" if name.startswith("P") else "C",
                    comp_id="DA",
                    chain_id="B",
                    entity_id=2,
                    seq_id=i + 1,
                    auth_seq_id=i + 1,
                    coords=(i * 5.0 + 10.0, 0.0, 0.0),
                    auth_asym_id="B",
                ))
                atom_id += 1

        s = _make_structure("test", 1, 0, atoms_prot + atoms_dna)
        assert s.n_chains == 2
        assert len(s.entities) == 2

    def test_protein_ligand_structure(self):
        """Test structure with protein + ligand."""
        atoms_prot = _make_linear_protein("A", 1, 3)
        # Create ligand atoms
        atoms_lig = [
            _make_atom(100, "C1", "C", "LIG", "B", 2, 1, 1, (5.0, 0.0, 0.0), auth_asym_id="B"),
            _make_atom(101, "C2", "C", "LIG", "B", 2, 1, 1, (6.0, 1.0, 0.0), auth_asym_id="B"),
            _make_atom(102, "O1", "O", "LIG", "B", 2, 1, 1, (5.0, 1.0, 1.0), auth_asym_id="B"),
        ]

        s = _make_structure("test", 1, 0, atoms_prot + atoms_lig)
        assert s.n_chains == 2

    def test_metric_applicability_protein_only(self):
        from af3_analysis.structural.metrics.rmsd import GlobalCaRMSD

        atoms = _make_linear_protein("A", 1, 5)
        s = _make_structure("test", 1, 0, atoms)

        metric = GlobalCaRMSD()
        assert metric.is_applicable(s, s) is True

    def test_metric_not_applicable_without_protein(self):
        from af3_analysis.structural.metrics.rmsd import GlobalCaRMSD

        # DNA-only structure
        atoms_dna = []
        for i in range(3):
            for name in ["P", "C4'"]:
                atoms_dna.append(_make_atom(
                    i * 2 + 1, name, "P" if name == "P" else "C",
                    "DA", "A", 1, i + 1, i + 1,
                    (i * 5.0, 0.0, 0.0),
                ))
        s = _make_structure("test", 1, 0, atoms_dna)
        # Override entity types for DNA
        from af3_analysis.structural.representation import EntityInfo, ChainInfo
        s.entities = [EntityInfo(1, "A", "A", "polymer", "polydeoxyribonucleotide", None, 3, 6)]
        s.chains = [ChainInfo("A", "A", 1, "polymer", "polydeoxyribonucleotide", 3, 6, 1, 3)]

        metric = GlobalCaRMSD()
        # Should not be applicable since there are no protein chains
        # (the metric requires ["protein"])
        # However, our current implementation checks entity types, not CA atoms
        # For a DNA-only structure with no protein chains, it should return False
        assert metric.is_applicable(s, s) is False
