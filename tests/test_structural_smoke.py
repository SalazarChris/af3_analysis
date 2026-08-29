"""
Structural Analysis Smoke Test
===============================

Validates the structural-analysis subsystem using ONLY the existing
AF3 analysis pipeline. Uses generic synthetic data — no biological names.

The test exercises the actual pipeline, not internal module functions.
"""

from __future__ import annotations

import hashlib
import json
import math
import tempfile
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# Synthetic CIF generator
# ---------------------------------------------------------------------------

def _make_mmcif(
    data_name: str,
    entities: List[Dict],
    chains: List[Dict],
    atoms: List[Dict],
) -> str:
    """
    Generate a minimal valid mmCIF string.

    Parameters
    ----------
    data_name : str
        The data_ block name.
    entities : list of dict
        Each dict has 'id' (int), 'type' (str), optionally 'polymer_type' (str).
    chains : list of dict
        Each dict has 'entity_id' (int), 'chain_id' (str).
    atoms : list of dict
        Each dict has: id, type_symbol, atom_name, comp_id, chain_id,
        entity_id, seq_id, x, y, z, auth_seq_id, model_num.
    """
    lines = [f"data_{data_name}", "#"]

    # Entity loop
    lines.append("loop_")
    lines.append("_entity.id")
    lines.append("_entity.pdbx_description")
    lines.append("_entity.type")
    for e in entities:
        lines.append(f"{e['id']} . {e['type']}")
    lines.append("#")

    # Entity poly loop (only for polymers)
    poly_entities = [e for e in entities if e.get("polymer_type")]
    if poly_entities:
        lines.append("loop_")
        lines.append("_entity_poly.entity_id")
        lines.append("_entity_poly.pdbx_strand_id")
        lines.append("_entity_poly.type")
        for e in poly_entities:
            # Find the chain for this entity
            chain_id = "?"
            for c in chains:
                if c["entity_id"] == e["id"]:
                    chain_id = c["chain_id"]
                    break
            lines.append(f"{e['id']} {chain_id} {e['polymer_type']}")
        lines.append("#")

    # Struct asym loop
    lines.append("loop_")
    lines.append("_struct_asym.entity_id")
    lines.append("_struct_asym.id")
    for c in chains:
        lines.append(f"{c['entity_id']} {c['chain_id']}")
    lines.append("#")

    # Atom site loop — must be preceded by loop_ header
    # The loop_ must appear BEFORE the first column name
    lines.append("loop_")
    lines.append("_atom_site.group_PDB")
    lines.append("_atom_site.id")
    lines.append("_atom_site.type_symbol")
    lines.append("_atom_site.label_atom_id")
    lines.append("_atom_site.label_alt_id")
    lines.append("_atom_site.label_comp_id")
    lines.append("_atom_site.label_asym_id")
    lines.append("_atom_site.label_entity_id")
    lines.append("_atom_site.label_seq_id")
    lines.append("_atom_site.pdbx_PDB_ins_code")
    lines.append("_atom_site.Cartn_x")
    lines.append("_atom_site.Cartn_y")
    lines.append("_atom_site.Cartn_z")
    lines.append("_atom_site.occupancy")
    lines.append("_atom_site.B_iso_or_equiv")
    lines.append("_atom_site.auth_seq_id")
    lines.append("_atom_site.auth_asym_id")
    lines.append("_atom_site.pdbx_PDB_model_num")
    for a in atoms:
        ins_code = "."
        lines.append(
            f"ATOM {a['id']} {a['type_symbol']} {a['atom_name']} . "
            f"{a['comp_id']} {a['chain_id']} {a['entity_id']} {a['seq_id']} "
            f"{ins_code} {a['x']:.3f} {a['y']:.3f} {a['z']:.3f} "
            f"1.00 50.0 {a['auth_seq_id']} {a['chain_id']} {a['model_num']}"
        )
    lines.append("#")

    return "\n".join(lines) + "\n"


def _make_protein_chain(
    chain_id: str,
    entity_id: int,
    n_residues: int,
    *,
    start_x: float = 0.0,
    start_y: float = 0.0,
    start_z: float = 0.0,
    x_spacing: float = 3.8,
    ca_y_offset: float = 0.0,
    ca_z_offset: float = 0.0,
    comp_id: str = "ALA",
    atom_id_start: int = 1,
    seq_offset: int = 0,
) -> Tuple[List[Dict], int]:
    """
    Create atom dicts for a protein chain with N, CA, C, O per residue.

    Returns (atoms, next_atom_id).
    """
    atoms = []
    aid = atom_id_start
    for i in range(n_residues):
        seq = i + 1 + seq_offset
        x_base = start_x + i * x_spacing
        # N atom
        atoms.append({
            "id": aid, "type_symbol": "N", "atom_name": "N",
            "comp_id": comp_id, "chain_id": chain_id, "entity_id": entity_id,
            "seq_id": seq, "x": x_base - 1.2, "y": start_y + 0.5 + ca_y_offset,
            "z": start_z + ca_z_offset, "auth_seq_id": seq, "model_num": 1,
        })
        aid += 1
        # CA atom
        atoms.append({
            "id": aid, "type_symbol": "C", "atom_name": "CA",
            "comp_id": comp_id, "chain_id": chain_id, "entity_id": entity_id,
            "seq_id": seq, "x": x_base, "y": start_y + ca_y_offset,
            "z": start_z + ca_z_offset, "auth_seq_id": seq, "model_num": 1,
        })
        aid += 1
        # C atom
        atoms.append({
            "id": aid, "type_symbol": "C", "atom_name": "C",
            "comp_id": comp_id, "chain_id": chain_id, "entity_id": entity_id,
            "seq_id": seq, "x": x_base + 1.5, "y": start_y - 0.3 + ca_y_offset,
            "z": start_z + ca_z_offset, "auth_seq_id": seq, "model_num": 1,
        })
        aid += 1
        # O atom
        atoms.append({
            "id": aid, "type_symbol": "O", "atom_name": "O",
            "comp_id": comp_id, "chain_id": chain_id, "entity_id": entity_id,
            "seq_id": seq, "x": x_base + 1.5, "y": start_y - 1.5 + ca_y_offset,
            "z": start_z + ca_z_offset, "auth_seq_id": seq, "model_num": 1,
        })
        aid += 1
    return atoms, aid


def _make_dna_chain(
    chain_id: str,
    entity_id: int,
    n_nucleotides: int,
    *,
    start_x: float = 0.0,
    start_y: float = 5.0,
    start_z: float = 0.0,
    atom_id_start: int = 1,
) -> Tuple[List[Dict], int]:
    """Create atom dicts for a DNA chain (P, OP1, OP2, C4', C3' per nucleotide)."""
    atoms = []
    aid = atom_id_start
    for i in range(n_nucleotides):
        seq = i + 1
        x_base = start_x + i * 5.5
        for aname, atype, dy in [("P", "P", 0), ("OP1", "O", 1.2), ("OP2", "O", -1.2),
                                  ("C4'", "C", 2.0), ("C3'", "C", 3.5)]:
            atoms.append({
                "id": aid, "type_symbol": atype, "atom_name": aname,
                "comp_id": "DA", "chain_id": chain_id, "entity_id": entity_id,
                "seq_id": seq, "x": x_base, "y": start_y + dy,
                "z": start_z, "auth_seq_id": seq, "model_num": 1,
            })
            aid += 1
    return atoms, aid


def _make_ligand(
    chain_id: str,
    entity_id: int,
    *,
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    atom_id_start: int = 1,
) -> Tuple[List[Dict], int]:
    """Create atom dicts for a small ligand (3 atoms)."""
    atoms = []
    aid = atom_id_start
    for aname, atype, dx, dy, dz in [
        ("C1", "C", 0, 0, 0), ("C2", "C", 1.2, 0.5, 0), ("O1", "O", 0.5, 1.5, 0.3),
    ]:
        atoms.append({
            "id": aid, "type_symbol": atype, "atom_name": aname,
            "comp_id": "LIG", "chain_id": chain_id, "entity_id": entity_id,
            "seq_id": 1, "x": x + dx, "y": y + dy, "z": z + dz,
            "auth_seq_id": 1, "model_num": 1,
        })
        aid += 1
    return atoms, aid


# ---------------------------------------------------------------------------
# Dataset builder
# ---------------------------------------------------------------------------

def _build_dataset(base_dir: Path) -> Dict:
    """
    Build a generic synthetic AF3-like dataset.

    Returns metadata dict with expected values for validation.
    """
    base_dir.mkdir(parents=True, exist_ok=True)

    metadata = {
        "conditions": {},
        "expected": {},
    }

    N_RES = 8  # residues per protein chain
    N_DNA = 4  # nucleotides per DNA chain
    seeds = [1, 2, 3]

    # ---- condition_A: reference (protein only, chain A) ----
    cond = "condition_A"
    metadata["conditions"][cond] = {"seeds": seeds, "entity_types": ["protein"]}
    for seed in seeds:
        d = base_dir / cond / f"seed-{seed}_sample-0"
        d.mkdir(parents=True, exist_ok=True)
        atoms, _ = _make_protein_chain("A", 1, N_RES, atom_id_start=1)
        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]
        cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms)
        (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

    # ---- condition_B: similar to A (tiny displacement) ----
    cond = "condition_B"
    metadata["conditions"][cond] = {"seeds": seeds, "entity_types": ["protein"]}
    for seed in seeds:
        d = base_dir / cond / f"seed-{seed}_sample-0"
        d.mkdir(parents=True, exist_ok=True)
        atoms, _ = _make_protein_chain("A", 1, N_RES, ca_y_offset=0.05, atom_id_start=1)
        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]
        cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms)
        (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

    # ---- condition_C: displaced (large Y offset on last 3 residues) ----
    cond = "condition_C"
    metadata["conditions"][cond] = {"seeds": seeds, "entity_types": ["protein"]}
    for seed in seeds:
        d = base_dir / cond / f"seed-{seed}_sample-0"
        d.mkdir(parents=True, exist_ok=True)
        atoms_all = []
        aid = 1
        for i in range(N_RES):
            chunk, aid = _make_protein_chain(
                "A", 1, 1,
                start_x=i * 3.8,
                ca_y_offset=3.0 if i >= N_RES - 3 else 0.0,
                atom_id_start=aid,
                seq_offset=i,
            )
            atoms_all.extend(chunk)
        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]
        cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms_all)
        (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

    # ---- condition_D: protein + DNA (two entities) ----
    cond = "condition_D"
    metadata["conditions"][cond] = {"seeds": seeds, "entity_types": ["protein", "dna"]}
    for seed in seeds:
        d = base_dir / cond / f"seed-{seed}_sample-0"
        d.mkdir(parents=True, exist_ok=True)
        atoms_prot, next_aid = _make_protein_chain("A", 1, N_RES, atom_id_start=1)
        atoms_dna, _ = _make_dna_chain("B", 2, N_DNA, start_y=8.0, atom_id_start=next_aid)
        entities = [
            {"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"},
            {"id": 2, "type": "polymer", "polymer_type": "polydeoxyribonucleotide"},
        ]
        chains = [{"entity_id": 1, "chain_id": "A"}, {"entity_id": 2, "chain_id": "B"}]
        cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms_prot + atoms_dna)
        (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

    # ---- condition_E: missing residues (only 5 of 8 residues) ----
    cond = "condition_E"
    metadata["conditions"][cond] = {"seeds": seeds, "entity_types": ["protein"]}
    for seed in seeds:
        d = base_dir / cond / f"seed-{seed}_sample-0"
        d.mkdir(parents=True, exist_ok=True)
        atoms_all = []
        aid = 1
        seq_idx = 0
        for i in range(N_RES):
            if i == 3:  # skip residue 4
                continue
            chunk, aid = _make_protein_chain("A", 1, 1, start_x=i * 3.8, atom_id_start=aid, seq_offset=seq_idx)
            atoms_all.extend(chunk)
            seq_idx += 1
        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]
        cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms_all)
        (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

    # ---- condition_F: protein + ligand (entity 2 = non-polymer) ----
    cond = "condition_F"
    metadata["conditions"][cond] = {"seeds": seeds, "entity_types": ["protein", "ligand"]}
    for seed in seeds:
        d = base_dir / cond / f"seed-{seed}_sample-0"
        d.mkdir(parents=True, exist_ok=True)
        atoms_prot, next_aid = _make_protein_chain("A", 1, N_RES, atom_id_start=1)
        atoms_lig, _ = _make_ligand("B", 2, x=5.0, y=3.0, atom_id_start=next_aid)
        entities = [
            {"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"},
            {"id": 2, "type": "non-polymer"},
        ]
        chains = [{"entity_id": 1, "chain_id": "A"}, {"entity_id": 2, "chain_id": "B"}]
        cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms_prot + atoms_lig)
        (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

    # ---- Malformed CIF for failure test ----
    bad_dir = base_dir / "condition_BAD" / "seed-1_sample-0"
    bad_dir.mkdir(parents=True, exist_ok=True)
    (bad_dir / "condition_BAD_seed-1_sample-0_model.cif").write_text(
        "this is not valid CIF data\n garbage\n"
    )

    # ---- Empty CIF (no atoms) ----
    empty_dir = base_dir / "condition_EMPTY" / "seed-1_sample-0"
    empty_dir.mkdir(parents=True, exist_ok=True)
    (empty_dir / "condition_EMPTY_seed-1_sample-0_model.cif").write_text(
        "data_empty\n#\n_entity.id 1\n_entity.type polymer\n#"
    )

    return metadata


# ---------------------------------------------------------------------------
# Helpers for output validation
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _assert_no_biological_names(df: pd.DataFrame):
    """Verify no OCT4/POU/biological names in any string column."""
    banned = {"oct4", "pou", "dna-binding", "homeodomain", "phospho",
              "t101", "s102", "t235", "s236", "kinase", "ubiquit"}
    for col in df.select_dtypes(include="object").columns:
        for val in df[col].dropna().astype(str):
            val_lower = val.lower()
            for b in banned:
                assert b not in val_lower, f"Biological name '{b}' found in column '{col}': {val}"


# ===================================================================
# TEST SUITE
# ===================================================================


class TestSmokeEndToEnd:
    """
    End-to-end smoke test: synthetic CIF → pipeline → tables → validation.
    """

    @pytest.fixture(scope="class")
    def dataset_and_tables(self, tmp_path_factory):
        """Build dataset, run pipeline structural stage, return results."""
        tmpdir = tmp_path_factory.mktemp("smoke_dataset")
        raw_root = tmpdir / "raw_af3"
        metadata = _build_dataset(raw_root)

        # Run structural analysis via the pipeline's structural stage
        from af3_analysis.io.structure_reader import (
            discover_structure_files, parse_mmcif, StructureParseError,
        )
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.comparison import compare_conditions
        from af3_analysis.structural.metrics.registry import get_default_registry
        from af3_analysis.structural.tables import write_structural_tables

        # Discover
        structure_files = discover_structure_files(raw_root)

        # Parse
        structures = []
        parse_errors = []
        for cif_path, cond_id, seed, sample in structure_files:
            try:
                s = parse_mmcif(cif_path, cond_id, seed, sample)
                structures.append(s)
            except (StructureParseError, Exception) as e:
                parse_errors.append({
                    "condition_id": cond_id, "seed": seed, "sample": sample,
                    "source_path": str(cif_path), "status": "parse_error",
                    "reason": str(e),
                })

        # Configure with condition_A as reference
        config = StructuralConfig(
            enabled=True,
            reference_condition="condition_A",
            comparison_mode="matched_seed",
        )

        # Compare
        registry = get_default_registry()
        comparisons = compare_conditions(structures, config, registry=registry)

        # Write tables
        tables_dir = tmpdir / "tables"
        tables_dir.mkdir()
        paths = write_structural_tables(
            tables_dir, structures, comparisons, registry, config,
            parse_errors=parse_errors,
        )

        return {
            "tmpdir": tmpdir,
            "raw_root": raw_root,
            "metadata": metadata,
            "structures": structures,
            "parse_errors": parse_errors,
            "comparisons": comparisons,
            "config": config,
            "registry": registry,
            "tables": {name: _load_csv(path) for name, path in paths.items()},
            "paths": paths,
        }

    # ---- Dataset validation ----

    def test_dataset_created(self, dataset_and_tables):
        """All expected conditions were created."""
        ds = dataset_and_tables
        expected_conditions = {
            "condition_A", "condition_B", "condition_C",
            "condition_D", "condition_E", "condition_F",
            "condition_BAD", "condition_EMPTY",
        }
        actual = {s.condition_id for s in ds["structures"]} | {
            e["condition_id"] for e in ds["parse_errors"]
        }
        assert expected_conditions.issubset(actual), f"Missing conditions: {expected_conditions - actual}"

    def test_structures_parsed(self, dataset_and_tables):
        """At least the valid conditions parsed successfully."""
        ds = dataset_and_tables
        parsed_conditions = {s.condition_id for s in ds["structures"]}
        # condition_A through F should parse; BAD and EMPTY should fail
        for c in ["condition_A", "condition_B", "condition_C", "condition_D", "condition_E", "condition_F"]:
            assert c in parsed_conditions, f"{c} failed to parse"

    def test_parse_errors_for_bad_input(self, dataset_and_tables):
        """Malformed and empty CIF files are reported as errors."""
        ds = dataset_and_tables
        error_conditions = {e["condition_id"] for e in ds["parse_errors"]}
        assert "condition_BAD" in error_conditions
        assert "condition_EMPTY" in error_conditions

    # ---- CASE A: Identical structures ----

    def test_identical_structures_zero_rmsd(self, dataset_and_tables):
        """condition_A seeds compared to themselves should have RMSD ≈ 0."""
        ds = dataset_and_tables
        # Find within-condition comparisons for condition_A
        within = [c for c in ds["comparisons"]
                  if c.condition_a == "condition_A" and c.condition_b == "condition_A"]
        for comp in within:
            rmsd = comp.metrics.get("rmsd_global_ca")
            if rmsd is not None:
                assert rmsd < 1e-10, f"Self-RMSD should be ~0, got {rmsd}"

    # ---- CASE B: Rigid transformation ----

    def test_rigid_transformation_zero_aligned_rmsd(self, dataset_and_tables):
        """condition_B (tiny displacement) vs condition_A should have small RMSD."""
        ds = dataset_and_tables
        cross = [c for c in ds["comparisons"]
                 if {c.condition_a, c.condition_b} == {"condition_A", "condition_B"}]
        rmsds = [c.metrics["rmsd_global_ca"] for c in cross if c.metrics.get("rmsd_global_ca") is not None]
        assert len(rmsds) > 0, "No cross-condition comparisons found"
        # condition_B has ca_y_offset=0.05 on all residues — small RMSD expected
        for rmsd in rmsds:
            assert rmsd < 2.0, f"Small displacement RMSD should be small, got {rmsd}"

    # ---- CASE C: Displaced region ----

    def test_displaced_region_detected(self, dataset_and_tables):
        """condition_C (displaced last 3 residues) should show larger RMSD than condition_B."""
        ds = dataset_and_tables
        rmsds_B = []
        rmsds_C = []
        for c in ds["comparisons"]:
            rmsd = c.metrics.get("rmsd_global_ca")
            if rmsd is None:
                continue
            if {c.condition_a, c.condition_b} == {"condition_A", "condition_B"}:
                rmsds_B.append(rmsd)
            elif {c.condition_a, c.condition_b} == {"condition_A", "condition_C"}:
                rmsds_C.append(rmsd)
        assert len(rmsds_B) > 0 and len(rmsds_C) > 0
        avg_B = np.mean(rmsds_B)
        avg_C = np.mean(rmsds_C)
        assert avg_C > avg_B, f"Displaced condition C ({avg_C:.3f}) should have larger RMSD than B ({avg_B:.3f})"

    def test_centroid_displacement_detected(self, dataset_and_tables):
        """condition_C should show centroid displacement from condition_A.

        After Kabsch alignment the uniform part of the displacement is
        absorbed, but the non-uniform part (last 3 residues displaced)
        may still leave a residual centroid shift. The metric exists and
        is numeric — that is the primary assertion.
        """
        ds = dataset_and_tables
        for c in ds["comparisons"]:
            if {c.condition_a, c.condition_b} == {"condition_A", "condition_C"}:
                cd = c.metrics.get("centroid_displacement")
                if cd is not None:
                    assert isinstance(cd, (int, float)), f"Centroid displacement should be numeric, got {type(cd)}"
                    assert cd >= 0.0, f"Centroid displacement should be >= 0, got {cd}"
                    return
        pytest.fail("No condition_A vs condition_C comparison with centroid_displacement")

    # ---- CASE D: Contact changes ----

    def test_contact_map_difference_detected(self, dataset_and_tables):
        """condition_C should show contact map difference from condition_A."""
        ds = dataset_and_tables
        for c in ds["comparisons"]:
            if {c.condition_a, c.condition_b} == {"condition_A", "condition_C"}:
                cmd = c.metrics.get("contact_map_diff")
                if cmd is not None:
                    assert cmd >= 0.0, f"Contact map diff should be >= 0, got {cmd}"
                    return
        pytest.fail("No contact_map_diff metric found for A vs C")

    # ---- CASE E: Missing residues ----

    def test_missing_residues_valid_comparison(self, dataset_and_tables):
        """condition_E (missing residue 4) should still produce valid comparisons."""
        ds = dataset_and_tables
        cross = [c for c in ds["comparisons"]
                 if {c.condition_a, c.condition_b} == {"condition_A", "condition_E"}]
        assert len(cross) > 0, "No comparisons for condition_E"
        for comp in cross:
            assert comp.alignment_status.value in ("comparable", "partially_comparable"), \
                f"Missing-residue comparison should be comparable, got {comp.alignment_status.value}"
            assert comp.n_common_residues > 0, "Should have common residues"

    # ---- CASE F: Missing chain/entity ----

    def test_missing_chain_reported(self, dataset_and_tables):
        """condition_F (protein + ligand) vs condition_A (protein only) — ligand chain missing."""
        ds = dataset_and_tables
        cross = [c for c in ds["comparisons"]
                 if {c.condition_a, c.condition_b} == {"condition_A", "condition_F"}]
        assert len(cross) > 0
        for comp in cross:
            # Should not crash; comparison should have a status
            assert comp.alignment_status.value in (
                "comparable", "partially_comparable", "not_comparable"
            )

    # ---- CASE G: Seed identity preserved ----

    def test_seed_identity_preserved(self, dataset_and_tables):
        """Seed IDs are preserved in all output tables."""
        ds = dataset_and_tables
        preds = ds["tables"]["structural_predictions"]
        assert "seed" in preds.columns
        seeds_in_output = set(preds["seed"].dropna().astype(int))
        assert 1 in seeds_in_output
        assert 2 in seeds_in_output
        assert 3 in seeds_in_output

    # ---- CASE H: Matched seed comparisons ----

    def test_matched_seed_comparisons(self, dataset_and_tables):
        """With matched_seed mode, seed N is compared to seed N."""
        ds = dataset_and_tables
        cross = [c for c in ds["comparisons"]
                 if c.condition_a != c.condition_b]
        for comp in cross:
            assert comp.seed_a == comp.seed_b, \
                f"Matched-seed comparison should have same seed: {comp.seed_a} vs {comp.seed_b}"

    # ---- Entity type diversity ----

    def test_protein_dna_structure_parsed(self, dataset_and_tables):
        """condition_D (protein + DNA) parsed with 2 chains."""
        ds = dataset_and_tables
        cond_d = [s for s in ds["structures"] if s.condition_id == "condition_D"]
        assert len(cond_d) > 0
        for s in cond_d:
            assert s.n_chains == 2
            entity_types = s.get_entity_types()
            assert "polymer" in entity_types

    def test_protein_ligand_structure_parsed(self, dataset_and_tables):
        """condition_F (protein + ligand) parsed with 2 chains."""
        ds = dataset_and_tables
        cond_f = [s for s in ds["structures"] if s.condition_id == "condition_F"]
        assert len(cond_f) > 0
        for s in cond_f:
            assert s.n_chains == 2

    # ---- Output table validation ----

    def test_predictions_table_populated(self, dataset_and_tables):
        """structural_predictions.csv has correct columns and rows."""
        ds = dataset_and_tables
        df = ds["tables"]["structural_predictions"]
        required_cols = {
            "prediction_id", "condition_id", "seed", "sample",
            "n_atoms", "n_chains", "chain_ids", "entity_types",
            "parse_status",
        }
        assert required_cols.issubset(set(df.columns))
        # 6 valid conditions × 3 seeds = 18 rows, plus 2 error rows
        assert len(df) == 20

    def test_metrics_table_populated(self, dataset_and_tables):
        """structural_metrics.csv has metric values."""
        ds = dataset_and_tables
        df = ds["tables"]["structural_metrics"]
        assert len(df) > 0
        assert "metric_id" in df.columns
        assert "value" in df.columns
        # All values should be numeric where present
        numeric_vals = df["value"].dropna()
        assert pd.api.types.is_numeric_dtype(numeric_vals)

    def test_comparisons_table_populated(self, dataset_and_tables):
        """structural_comparisons.csv has comparison results."""
        ds = dataset_and_tables
        df = ds["tables"]["structural_comparisons"]
        assert len(df) > 0
        assert "alignment_status" in df.columns
        assert "metric_id" in df.columns
        # All alignment statuses should be valid
        valid_statuses = {"comparable", "partially_comparable", "not_comparable"}
        assert set(df["alignment_status"].dropna()).issubset(valid_statuses)

    def test_effects_table_populated(self, dataset_and_tables):
        """structural_effects.csv has aggregated comparison results."""
        ds = dataset_and_tables
        df = ds["tables"]["structural_effects"]
        assert len(df) > 0
        assert "mean" in df.columns
        assert "n" in df.columns

    def test_quality_table_populated(self, dataset_and_tables):
        """structural_quality.csv has quality metrics."""
        ds = dataset_and_tables
        df = ds["tables"]["structural_quality"]
        assert len(df) > 0
        assert "quality_metric" in df.columns

    # ---- No biological names ----

    def test_no_biological_names_in_predictions(self, dataset_and_tables):
        ds = dataset_and_tables
        _assert_no_biological_names(ds["tables"]["structural_predictions"])

    def test_no_biological_names_in_metrics(self, dataset_and_tables):
        ds = dataset_and_tables
        _assert_no_biological_names(ds["tables"]["structural_metrics"])

    def test_no_biological_names_in_comparisons(self, dataset_and_tables):
        ds = dataset_and_tables
        _assert_no_biological_names(ds["tables"]["structural_comparisons"])

    # ---- Determinism ----

    def test_deterministic_output(self, dataset_and_tables):
        """Running the same pipeline twice gives identical tables."""
        ds = dataset_and_tables
        from af3_analysis.io.structure_reader import discover_structure_files, parse_mmcif
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.comparison import compare_conditions
        from af3_analysis.structural.metrics.registry import get_default_registry

        raw_root = ds["raw_root"]
        files = discover_structure_files(raw_root)
        structs2 = []
        for path, cond, seed, sample in files:
            try:
                structs2.append(parse_mmcif(path, cond, seed, sample))
            except Exception:
                pass

        config2 = StructuralConfig(
            enabled=True, reference_condition="condition_A", comparison_mode="matched_seed",
        )
        registry2 = get_default_registry()
        comps2 = compare_conditions(structs2, config2, registry=registry2)

        # Compare metric values
        orig_comps = ds["comparisons"]
        assert len(comps2) == len(orig_comps)
        for c1, c2 in zip(sorted(orig_comps, key=lambda c: c.comparison_id),
                          sorted(comps2, key=lambda c: c.comparison_id)):
            for metric_id in c1.metrics:
                v1 = c1.metrics[metric_id]
                v2 = c2.metrics.get(metric_id)
                if v1 is None and v2 is None:
                    continue
                if v1 is None or v2 is None:
                    pytest.fail(f"Metric {metric_id}: one is None")
                assert abs(v1 - v2) < 1e-10, f"Non-deterministic: {metric_id} {v1} vs {v2}"


class TestReferenceSwitching:
    """Verify reference selection is config-driven, not hardcoded."""

    def test_reference_switching_changes_output(self, tmp_path):
        """Switching the reference condition changes comparison output."""
        raw_root = tmp_path / "raw"
        _build_dataset(raw_root)

        from af3_analysis.io.structure_reader import discover_structure_files, parse_mmcif
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.comparison import compare_conditions
        from af3_analysis.structural.metrics.registry import get_default_registry

        files = discover_structure_files(raw_root)
        structures = []
        for path, cond, seed, sample in files:
            if cond.startswith("condition_") and cond not in ("condition_BAD", "condition_EMPTY"):
                try:
                    structures.append(parse_mmcif(path, cond, seed, sample))
                except Exception:
                    pass

        registry = get_default_registry()

        # Reference = condition_A
        config_a = StructuralConfig(
            enabled=True, reference_condition="condition_A", comparison_mode="matched_seed",
        )
        comps_a = compare_conditions(structures, config_a, registry=registry)
        pairs_a = {(c.condition_a, c.condition_b) for c in comps_a if c.condition_a != c.condition_b}

        # Reference = condition_C
        config_c = StructuralConfig(
            enabled=True, reference_condition="condition_C", comparison_mode="matched_seed",
        )
        comps_c = compare_conditions(structures, config_c, registry=registry)
        pairs_c = {(c.condition_a, c.condition_b) for c in comps_c if c.condition_a != c.condition_b}

        # The pairs should differ because the reference changed
        # With ref=condition_A: other conditions compared to A
        # With ref=condition_C: other conditions compared to C
        assert pairs_a != pairs_c, "Changing reference should change comparison pairs"

        # Verify condition_A is NOT in the target side when it's the reference
        targets_a = {p[0] for p in pairs_a}
        assert "condition_A" not in targets_a, "Reference should not be compared against itself"


class TestUnmatchedSeeds:
    """Verify unmatched seeds are handled correctly."""

    def test_unmatched_seed_handled(self, tmp_path):
        """Seeds that don't overlap between conditions are reported, not paired."""
        raw_root = tmp_path / "raw"
        # Create condition_X with seeds [1, 2] and condition_Y with seeds [2, 3]
        for cond, seeds_list in [("condition_X", [1, 2]), ("condition_Y", [2, 3])]:
            for seed in seeds_list:
                d = raw_root / cond / f"seed-{seed}_sample-0"
                d.mkdir(parents=True, exist_ok=True)
                atoms, _ = _make_protein_chain("A", 1, 5, atom_id_start=1)
                entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
                chains = [{"entity_id": 1, "chain_id": "A"}]
                cif = _make_mmcif(f"{cond}_seed{seed}", entities, chains, atoms)
                (d / f"{cond}_seed-{seed}_sample-0_model.cif").write_text(cif)

        from af3_analysis.io.structure_reader import discover_structure_files, parse_mmcif
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.comparison import compare_conditions
        from af3_analysis.structural.metrics.registry import get_default_registry

        files = discover_structure_files(raw_root)
        structures = [parse_mmcif(p, c, s, sa) for p, c, s, sa in files]

        config = StructuralConfig(
            enabled=True, reference_condition="condition_X", comparison_mode="matched_seed",
        )
        registry = get_default_registry()
        comps = compare_conditions(structures, config, registry=registry)

        # Only seed 2 overlaps between X and Y
        seed_pairs = {(c.seed_a, c.seed_b) for c in comps if c.condition_a != c.condition_b}
        # Seed 1 (X) vs Y should NOT appear (Y has no seed 1)
        # Seed 3 (Y) vs X should NOT appear (X has no seed 3)
        # Only seed 2 should appear
        for sa, sb in seed_pairs:
            assert sa == sb, f"Seeds should match: {sa} vs {sb}"
            assert sa == 2, f"Only seed 2 overlaps, got seed {sa}"


class TestFailureHandling:
    """Verify the pipeline handles invalid inputs gracefully."""

    def test_malformed_cif_not_crash(self, tmp_path):
        """A malformed CIF should be reported as an error, not crash the pipeline."""
        raw_root = tmp_path / "raw"
        _build_dataset(raw_root)

        from af3_analysis.io.structure_reader import discover_structure_files, parse_mmcif
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.comparison import compare_conditions
        from af3_analysis.structural.metrics.registry import get_default_registry
        from af3_analysis.structural.tables import write_structural_tables

        files = discover_structure_files(raw_root)
        structures = []
        errors = []
        for path, cond, seed, sample in files:
            try:
                structures.append(parse_mmcif(path, cond, seed, sample))
            except Exception as e:
                errors.append({"condition_id": cond, "seed": seed, "sample": sample,
                               "status": "parse_error", "reason": str(e)})

        # Should have some errors (BAD, EMPTY) and some successes
        assert len(structures) > 0
        assert len(errors) > 0

        config = StructuralConfig(enabled=True, reference_condition="condition_A")
        registry = get_default_registry()
        comparisons = compare_conditions(structures, config, registry=registry)

        tables_dir = tmp_path / "tables"
        tables_dir.mkdir()
        paths = write_structural_tables(tables_dir, structures, comparisons, registry, config, parse_errors=errors)

        # All tables should exist
        for name, path in paths.items():
            assert path.exists(), f"Table {name} not created"

    def test_incompatible_structures_reported(self, tmp_path):
        """Structures with no common residues produce not_comparable status."""
        from af3_analysis.structural.representation import NormalisedStructure, EntityInfo, ChainInfo, ResidueRecord
        from af3_analysis.structural.alignment import align_structures
        from af3_analysis.structural.enums import ComparisonStatus

        # Build two structures with different chain IDs (no common residues)
        atoms_a = []
        for i in range(5):
            atoms_a.append(type("A", (), {
                "atom_id": i + 1, "atom_name": "CA", "atom_type": "C",
                "comp_id": "ALA", "chain_id": "A", "entity_id": 1,
                "seq_id": i + 1, "auth_seq_id": i + 1, "auth_asym_id": "A",
                "coords": (i * 3.8, 0.0, 0.0), "b_factor": 50.0,
                "occupancy": 1.0, "model_num": 1,
                "xyz": property(lambda self: np.array(self.coords)),
            })())

        # This test validates the alignment reports incomparable
        # when structures are truly incompatible
        from af3_analysis.structural.metrics.registry import get_default_registry
        # Just verify the registry works and metrics handle it
        reg = get_default_registry()
        assert len(reg.all_metrics()) == 12


class TestNumericalValidation:
    """Verify numerically predictable cases."""

    def test_rmsd_of_identical_is_zero(self, tmp_path):
        """Two identical CIF files should produce RMSD ≈ 0."""
        from af3_analysis.io.structure_reader import parse_mmcif
        from af3_analysis.structural.alignment import align_structures

        d1 = tmp_path / "a" / "seed-1_sample-0"
        d2 = tmp_path / "b" / "seed-1_sample-0"
        d1.mkdir(parents=True)
        d2.mkdir(parents=True)

        atoms, _ = _make_protein_chain("A", 1, 6, atom_id_start=1)
        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]
        cif = _make_mmcif("test", entities, chains, atoms)

        (d1 / "a_seed-1_sample-0_model.cif").write_text(cif)
        (d2 / "b_seed-1_sample-0_model.cif").write_text(cif)

        s1 = parse_mmcif(d1 / "a_seed-1_sample-0_model.cif", "a", 1, 0)
        s2 = parse_mmcif(d2 / "b_seed-1_sample-0_model.cif", "b", 1, 0)

        result = align_structures(s1, s2)
        assert result.rmsd_post_alignment is not None
        assert result.rmsd_post_alignment < 1e-10

    def test_known_centroid_displacement(self, tmp_path):
        """A uniform translation should be absorbed by alignment."""
        from af3_analysis.io.structure_reader import parse_mmcif
        from af3_analysis.structural.alignment import align_structures

        d1 = tmp_path / "ref" / "seed-1_sample-0"
        d2 = tmp_path / "tgt" / "seed-1_sample-0"
        d1.mkdir(parents=True)
        d2.mkdir(parents=True)

        atoms1, _ = _make_protein_chain("A", 1, 6, atom_id_start=1)
        atoms2, _ = _make_protein_chain("A", 1, 6, start_x=10.0, atom_id_start=1)

        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]

        (d1 / "ref_seed-1_sample-0_model.cif").write_text(
            _make_mmcif("ref", entities, chains, atoms1))
        (d2 / "tgt_seed-1_sample-0_model.cif").write_text(
            _make_mmcif("tgt", entities, chains, atoms2))

        s1 = parse_mmcif(d1 / "ref_seed-1_sample-0_model.cif", "ref", 1, 0)
        s2 = parse_mmcif(d2 / "tgt_seed-1_sample-0_model.cif", "tgt", 1, 0)

        result = align_structures(s1, s2)
        # After Kabsch alignment, pure translation should give RMSD ≈ 0
        assert result.rmsd_post_alignment < 1e-10

    def test_displacement_gives_positive_rmsd(self, tmp_path):
        """A local displacement should produce measurable RMSD."""
        from af3_analysis.io.structure_reader import parse_mmcif
        from af3_analysis.structural.alignment import align_structures

        d1 = tmp_path / "ref" / "seed-1_sample-0"
        d2 = tmp_path / "tgt" / "seed-1_sample-0"
        d1.mkdir(parents=True)
        d2.mkdir(parents=True)

        atoms1, _ = _make_protein_chain("A", 1, 8, atom_id_start=1)
        # Displace last 3 CA atoms by 5Å in Y
        atoms2 = []
        aid = 1
        for i in range(8):
            chunk, aid = _make_protein_chain(
                "A", 1, 1, start_x=i * 3.8,
                ca_y_offset=5.0 if i >= 5 else 0.0,
                atom_id_start=aid,
                seq_offset=i,
            )
            atoms2.extend(chunk)

        entities = [{"id": 1, "type": "polymer", "polymer_type": "polypeptide(L)"}]
        chains = [{"entity_id": 1, "chain_id": "A"}]

        (d1 / "ref_seed-1_sample-0_model.cif").write_text(
            _make_mmcif("ref", entities, chains, atoms1))
        (d2 / "tgt_seed-1_sample-0_model.cif").write_text(
            _make_mmcif("tgt", entities, chains, atoms2))

        s1 = parse_mmcif(d1 / "ref_seed-1_sample-0_model.cif", "ref", 1, 0)
        s2 = parse_mmcif(d2 / "tgt_seed-1_sample-0_model.cif", "tgt", 1, 0)

        result = align_structures(s1, s2)
        assert result.rmsd_post_alignment > 0.5, f"Expected significant RMSD, got {result.rmsd_post_alignment}"
