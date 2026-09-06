"""
V3 adapter unit tests.

Covers the read-only translation layer between existing pipeline outputs
and the V3 data model, specifically:

- _replicate_to_stem: replicate id -> CIF condition stem reduction
- _build_condition_stem_map: registry/table condition_name -> condition_id
- _convert_to_structure_data: NormalisedStructure -> StructureData
  (including prediction_id construction, which previously crashed with
  AttributeError on real CIFs)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from af3_analysis.visualization_v3.adapter import (
    _replicate_to_stem,
    _build_condition_stem_map,
    _convert_to_structure_data,
)
from af3_analysis.io.structure_reader import parse_mmcif

TESTDATA = Path("testdata/pou2")
BASELINE_CIF = (
    TESTDATA / "pou_baseline" / "seed-10_sample-0"
    / "pou_baseline_seed-10_sample-0_model.cif"
)


# ---------------------------------------------------------------------------
# _replicate_to_stem
# ---------------------------------------------------------------------------

class TestReplicateToStem:
    @pytest.mark.parametrize(
        ("replicate_id", "expected"),
        [
            ("pou_baseline", "pou_baseline"),
            ("pou_baseline_summary", "pou_baseline"),
            ("pou_baseline_seed-10_sample-0", "pou_baseline"),
            ("pou_baseline_seed-10_sample-0_summary", "pou_baseline"),
            ("oct4__k123-sumo_seed-3_sample-4", "oct4__k123-sumo"),
            ("oct4__k123-sumo__dna_summary", "oct4__k123-sumo__dna"),
            ("pou_baseline_NA100_HOH1000_CL100_seed-2_sample-1",
             "pou_baseline_NA100_HOH1000_CL100"),
        ],
    )
    def test_reductions(self, replicate_id, expected):
        assert _replicate_to_stem(replicate_id) == expected

    def test_empty(self):
        assert _replicate_to_stem("_summary") is None


# ---------------------------------------------------------------------------
# _build_condition_stem_map
# ---------------------------------------------------------------------------

class TestBuildConditionStemMap:
    def test_maps_registry_stems_to_condition_ids(self, tmp_path):
        tables = tmp_path / "tables"
        tables.mkdir()
        (tables / "condition_registry.csv").write_text(
            'condition_id,condition_name,n_replicates,replicate_ids\n'
            'cond_001,pou_baseline,3,"[\'\'pou_baseline\'\', '
            '\'\'pou_baseline_summary\'\', '
            '\'\'pou_baseline_seed-10_sample-0\'\']"\n',
            encoding="utf-8",
        )
        seed_aggregated = pd.DataFrame({
            "condition_id": ["cond_001"],
            "condition_name": ["pou_baseline"],
            "seed": [10],
            "pLDDT_mean": [80.0],
        })
        stem_map = _build_condition_stem_map(tmp_path, seed_aggregated)
        assert stem_map["pou_baseline"] == "cond_001"

    def test_condition_name_fallback_without_registry(self, tmp_path):
        seed_aggregated = pd.DataFrame({
            "condition_id": ["cond_002", "cond_002"],
            "condition_name": ["pou_sep102", "pou_sep102"],
            "seed": [1, 2],
            "pLDDT_mean": [80.0, 81.0],
        })
        stem_map = _build_condition_stem_map(tmp_path, seed_aggregated)
        assert stem_map["pou_sep102"] == "cond_002"

    def test_missing_registry_is_silent(self, tmp_path):
        seed_aggregated = pd.DataFrame({
            "condition_id": ["cond_001"],
            "condition_name": ["pou_baseline"],
        })
        stem_map = _build_condition_stem_map(tmp_path, seed_aggregated)
        assert stem_map["pou_baseline"] == "cond_001"


# ---------------------------------------------------------------------------
# _convert_to_structure_data
# ---------------------------------------------------------------------------

class TestConvertToStructureData:
    @pytest.mark.skipif(
        not BASELINE_CIF.is_file(),
        reason="testdata/pou2 CIF not available",
    )
    def test_converts_real_cif_without_crash(self):
        structure = parse_mmcif(
            BASELINE_CIF,
            condition_id="cond_001",
            seed=10,
            sample=0,
        )
        data = _convert_to_structure_data(structure, None)

        assert data.prediction_id == "cond_001_seed-10_sample-0"
        assert data.condition_id == "cond_001"
        assert data.seed == 10
        assert data.sample == 0
        assert data.parse_status == "success"
        assert data.entities, "expected at least one entity"
        assert any(
            e.polymer_type in ("polypeptide(L)", "polydeoxyribonucleotide")
            for e in data.entities
        )

    def test_prediction_id_follows_project_convention(self, tmp_path):
        """Regression: NormalisedStructure has no prediction_id attribute."""
        # Build via the real path when possible; otherwise assert the
        # conversion builds prediction_id from condition/seed/sample.
        from af3_analysis.visualization_v3 import model as v3_model

        class _FakeEntity:
            entity_id = 1
            entity_type = "polypeptide(L)"
            polymer_type = "polypeptide"
            description = "toy"

        class _FakeChain:
            chain_id = "A"
            entity_id = 1
            entity_type = "polypeptide(L)"
            polymer_type = "polypeptide"
            residues = []

        class _FakeStructure:
            source_path = tmp_path / "toy.cif"
            condition_id = "cond_001"
            seed = 3
            sample = 4
            entities = [_FakeEntity()]
            chains = [_FakeChain()]
            residues = []

            def get_chain(self, chain_id):
                return []

            def get_protein_chains(self):
                return [_FakeChain()]

            def get_nucleic_acid_chains(self):
                return []

        data = _convert_to_structure_data(_FakeStructure(), None)
        assert data.prediction_id == "cond_001_seed-3_sample-4"
        assert data.condition_id == "cond_001"