"""
Unit tests for the smoke-test runner.

These test the runner logic (condition selection, JSON generation,
metadata recording, failure handling) with mocked AF3 execution.
"""

import json
import csv
import pytest
from pathlib import Path

from af3_builder.condition_manifest import (
    load_master_manifest,
    load_protein_registry,
    load_construct_registry,
    load_modification_registry,
    load_nucleic_acid_registry,
    load_ligand_registry,
    load_ion_registry,
    load_af3_compatibility_registry,
    load_covalent_bond_registry,
)
from af3_builder.condition_manifest.builder import (
    build_job,
    resolve_condition,
    _validate_spec_for_build,
)
from af3_builder.validation.validator import AF3Validator, ValidationError


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).resolve().parents[2]  # FINALVERSIONTHESIS/
POU_DIR = _REPO_ROOT / "testdata" / "pou2" / "registries"
PX_DIR = _REPO_ROOT / "testdata" / "protein_x" / "registries"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def pou_registries():
    return {
        "protein_registry": load_protein_registry(POU_DIR / "protein_registry.csv"),
        "construct_registry": load_construct_registry(POU_DIR / "construct_registry.csv"),
        "modification_registry": load_modification_registry(POU_DIR / "modification_registry.csv"),
        "nucleic_acid_registry": load_nucleic_acid_registry(POU_DIR / "nucleic_acid_registry.csv"),
        "ligand_registry": load_ligand_registry(POU_DIR / "ligand_registry.csv"),
        "ion_registry": load_ion_registry(POU_DIR / "ion_registry.csv"),
        "af3_compatibility_registry": load_af3_compatibility_registry(POU_DIR / "af3_compatibility_registry.csv"),
        "covalent_bond_registry": load_covalent_bond_registry(POU_DIR / "covalent_bond_registry.csv"),
    }


@pytest.fixture
def pou_manifest():
    return load_master_manifest(
        POU_DIR / "master_condition_manifest.csv",
        modifications_path=POU_DIR / "condition_modifications.csv",
        entities_path=POU_DIR / "condition_entities.csv",
        factors_path=POU_DIR / "condition_factors.csv",
    )


# ---------------------------------------------------------------------------
# Test: Full pipeline for each smoke-test category
# ---------------------------------------------------------------------------

class TestSmokeTestPipeline:
    """Verify each smoke-test category produces valid JSON."""

    @pytest.mark.parametrize("condition_id,category", [
        ("pou_baseline", "A"),
        ("pou_tpo101", "B"),
        ("pou_tpo101_sep102", "C"),
        ("pou_dna", "D"),
        ("pou_tpo101_dna", "E"),
        ("pou_tpo101_sep102_dna", "E_extended"),
    ])
    def test_pou_condition_generates_valid_json(
        self, condition_id, category, pou_manifest, pou_registries
    ):
        """Each POU condition should produce valid AF3 JSON."""
        spec = resolve_condition(pou_manifest, condition_id, **pou_registries)
        _validate_spec_for_build(spec)
        jb = build_job(pou_manifest, condition_id, seeds=[1], **pou_registries)
        d = jb.to_dict()

        # Structural validation
        assert d["dialect"] == "alphafold3"
        assert d["modelSeeds"] == [1]
        assert len(d["sequences"]) >= 1

        # AF3Validator
        AF3Validator.validate_job(d, require_files=False)

    def test_pou_baseline_has_one_protein(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_baseline", **pou_registries)
        assert len(spec.proteins) == 1
        assert spec.proteins[0].entity_id == "POU_DOMAIN"

    def test_pou_baseline_has_ion(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_baseline", **pou_registries)
        assert len(spec.ions) >= 1
        assert spec.ions[0].ccd_code == "MG"

    def test_pou_tpo101_has_one_modification(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_tpo101", **pou_registries)
        mods = spec.proteins[0].modifications
        assert len(mods) == 1
        assert mods[0]["ccd_code"] == "TPO"

    def test_pou_tpo101_sep102_has_two_modifications(
        self, pou_manifest, pou_registries
    ):
        spec = resolve_condition(pou_manifest, "pou_tpo101_sep102", **pou_registries)
        mods = spec.proteins[0].modifications
        assert len(mods) == 2

    def test_pou_dna_has_dna_entity(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_dna", **pou_registries)
        assert len(spec.dna) == 1
        assert len(spec.dna[0].sequence) > 0

    def test_pou_tpo101_dna_has_both(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_tpo101_dna", **pou_registries)
        assert len(spec.proteins[0].modifications) == 1
        assert len(spec.dna) == 1

    def test_json_modification_format(self, pou_manifest, pou_registries):
        """Modifications should use AF3 camelCase format."""
        jb = build_job(pou_manifest, "pou_tpo101", seeds=[1], **pou_registries)
        d = jb.to_dict()
        prot = [s for s in d["sequences"] if "protein" in s][0]["protein"]
        mod = prot["modifications"][0]
        assert "ccdCode" in mod
        assert "position" in mod
        assert "modification_id" not in mod
        assert "ccd_code" not in mod
        assert "af3_status" not in mod


# ---------------------------------------------------------------------------
# Module-level Protein X fixtures (shared across test classes)
# ---------------------------------------------------------------------------

@pytest.fixture
def px_registries():
    return {
        "protein_registry": load_protein_registry(PX_DIR / "protein_registry.csv"),
        "construct_registry": load_construct_registry(PX_DIR / "construct_registry.csv"),
        "modification_registry": load_modification_registry(PX_DIR / "modification_registry.csv"),
        "nucleic_acid_registry": load_nucleic_acid_registry(PX_DIR / "nucleic_acid_registry.csv"),
        "ligand_registry": load_ligand_registry(PX_DIR / "ligand_registry.csv"),
        "ion_registry": load_ion_registry(PX_DIR / "ion_registry.csv"),
        "af3_compatibility_registry": load_af3_compatibility_registry(PX_DIR / "af3_compatibility_registry.csv"),
        "covalent_bond_registry": load_covalent_bond_registry(PX_DIR / "covalent_bond_registry.csv"),
    }


@pytest.fixture
def px_manifest():
    return load_master_manifest(
        PX_DIR / "master_condition_manifest.csv",
        modifications_path=PX_DIR / "condition_modifications.csv",
        entities_path=PX_DIR / "condition_entities.csv",
        factors_path=PX_DIR / "condition_factors.csv",
    )


# ---------------------------------------------------------------------------
# Test: Protein X (genericity)
# ---------------------------------------------------------------------------

class TestProteinXSmokeTest:
    """Verify Protein X works through the same pipeline."""

    @pytest.mark.parametrize("condition_id,category", [
        ("kinase_baseline", "A_px"),
        ("kinase_ligand_A", "F"),
        ("kinase_phospho", "B_px"),
        ("kinase_phospho_ligand", "H"),
        ("kinase_multi", "H_extended"),
    ])
    def test_px_condition_generates_valid_json(
        self, condition_id, category, px_manifest, px_registries
    ):
        spec = resolve_condition(px_manifest, condition_id, **px_registries)
        _validate_spec_for_build(spec)
        jb = build_job(px_manifest, condition_id, seeds=[1], **px_registries)
        d = jb.to_dict()
        assert d["dialect"] == "alphafold3"
        AF3Validator.validate_job(d, require_files=False)

    def test_px_ligand_has_atp(self, px_manifest, px_registries):
        spec = resolve_condition(px_manifest, "kinase_ligand_A", **px_registries)
        lig = [l for l in spec.ligands if l.entity_id == "example_ligand_A"]
        assert len(lig) == 1
        assert lig[0].ccd_code == "ATP"

    def test_px_ligand_json_format(self, px_manifest, px_registries):
        jb = build_job(px_manifest, "kinase_ligand_A", seeds=[1], **px_registries)
        d = jb.to_dict()
        ligs = [s for s in d["sequences"] if "ligand" in s
                if "ccdCodes" in s.get("ligand", {})]
        assert len(ligs) >= 1
        assert ligs[0]["ligand"]["ccdCodes"] == ["ATP"]


# ---------------------------------------------------------------------------
# Test: Error handling
# ---------------------------------------------------------------------------

class TestSmokeTestErrorHandling:
    """Verify the runner handles errors gracefully."""

    def test_unknown_condition_fails(self, pou_manifest, pou_registries):
        with pytest.raises(ValueError, match="not found"):
            resolve_condition(pou_manifest, "nonexistent", **pou_registries)

    def test_unsupported_modification_fails(self):
        from af3_builder.condition_manifest.manifest import (
            MasterManifest, ConditionRecord, ConditionModificationRecord,
            ConditionEntityRecord,
        )
        from af3_builder.condition_manifest.registries import (
            ModificationRecord, AF3CompatibilityRecord,
        )

        manifest = MasterManifest()
        manifest.conditions["c1"] = ConditionRecord(
            condition_id="c1", condition_name="Bad"
        )
        manifest.modifications["m1"] = ConditionModificationRecord(
            condition_id="c1", modification_id="bad_mod",
            sequence_position="10", construct_id="POU_DOMAIN",
        )
        manifest.entities["e1"] = ConditionEntityRecord(
            condition_id="c1", entity_type="protein",
            entity_id="POU_DOMAIN", stoichiometry="1",
        )

        mod_reg = {"bad_mod": ModificationRecord(modification_id="bad_mod")}
        af3_reg = {
            "rep": AF3CompatibilityRecord(
                representation_id="rep",
                modification_id="bad_mod",
                af3_status="unsupported",
            )
        }
        construct_reg = load_construct_registry(POU_DIR / "construct_registry.csv")

        with pytest.raises(ValueError, match="UNSUPPORTED"):
            build_job(
                manifest, "c1", seeds=[1],
                construct_registry=construct_reg,
                modification_registry=mod_reg,
                af3_compatibility_registry=af3_reg,
            )


# ---------------------------------------------------------------------------
# Test: JSON component completeness (no silent loss)
# ---------------------------------------------------------------------------

class TestNoSilentComponentLoss:
    """Verify no required biological component disappears during serialization."""

    def test_dna_condition_retains_dna(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_dna", **pou_registries)
        jb = build_job(pou_manifest, "pou_dna", seeds=[1], **pou_registries)
        d = jb.to_dict()

        types = [list(s.keys())[0] for s in d["sequences"]]
        assert "dna" in types

    def test_modification_condition_retains_mods(self, pou_manifest, pou_registries):
        spec = resolve_condition(pou_manifest, "pou_tpo101", **pou_registries)
        jb = build_job(pou_manifest, "pou_tpo101", seeds=[1], **pou_registries)
        d = jb.to_dict()

        prot = [s for s in d["sequences"] if "protein" in s][0]["protein"]
        assert "modifications" in prot
        assert len(prot["modifications"]) == 1

    def test_complex_condition_retains_all(self, pou_manifest, pou_registries):
        spec = resolve_condition(
            pou_manifest, "pou_tpo101_sep102_dna", **pou_registries
        )
        jb = build_job(
            pou_manifest, "pou_tpo101_sep102_dna", seeds=[1], **pou_registries
        )
        d = jb.to_dict()

        types = [list(s.keys())[0] for s in d["sequences"]]
        assert "protein" in types
        assert "dna" in types

        prot = [s for s in d["sequences"] if "protein" in s][0]["protein"]
        assert len(prot["modifications"]) == 2

    def test_ligand_condition_retains_ligand(self, px_manifest, px_registries):
        spec = resolve_condition(px_manifest, "kinase_ligand_A", **px_registries)
        jb = build_job(px_manifest, "kinase_ligand_A", seeds=[1], **px_registries)
        d = jb.to_dict()

        types = [list(s.keys())[0] for s in d["sequences"]]
        assert "ligand" in types


# ---------------------------------------------------------------------------
# Test: Determinism
# ---------------------------------------------------------------------------

class TestDeterminism:
    """Same inputs should produce same outputs."""

    def test_deterministic_output(self, pou_manifest, pou_registries):
        jb1 = build_job(pou_manifest, "pou_baseline", seeds=[42], **pou_registries)
        jb2 = build_job(pou_manifest, "pou_baseline", seeds=[42], **pou_registries)
        assert jb1.to_dict() == jb2.to_dict()

    def test_different_seeds_different_json(self, pou_manifest, pou_registries):
        jb1 = build_job(pou_manifest, "pou_baseline", seeds=[1], **pou_registries)
        jb2 = build_job(pou_manifest, "pou_baseline", seeds=[2], **pou_registries)
        assert jb1.to_dict()["modelSeeds"] != jb2.to_dict()["modelSeeds"]
