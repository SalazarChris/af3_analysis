"""
V3 visualization pipeline smoke tests.

Uses a small synthetic dataset built from two toy structures, so tests run
without real AF3 data. Structural correctness is exercised elsewhere;
these tests verify wiring, unit-of-analysis pairing, and table/manifest
outputs of the runner.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from af3_analysis.visualization_v3.config import V3Config
from af3_analysis.visualization_v3.model import (
    AtomPosition,
    ChainGeometry,
    EntityGeometry,
    ResidueGeometry,
    StructureData,
    Dataset,
    ConditionSummary,
    SeedSummary,
)
from af3_analysis.visualization_v3 import runner as v3_runner
from af3_analysis.visualization_v3.reference import resolve_reference
from af3_analysis.visualization_v3.figures import (
    generate_f02_global_structural_difference,
    generate_f03_matched_seed_structural_difference,
    generate_f10_local_geometry,
    generate_f11_domain_motion,
    generate_f12_structural_clustering,
    generate_f13_similarity_matrix,
    generate_f14_mds_embedding,
    generate_f16_confidence_change_vs_structural_change,
    generate_f19_factorial_effects,
)


# ---------------------------------------------------------------------------
# Synthetic structure helpers
# ---------------------------------------------------------------------------

def _make_residue(chain_id: str, entity_id: int, seq_id: int, xyz) -> ResidueGeometry:
    return ResidueGeometry(
        residue_name="ALA",
        chain_id=chain_id,
        auth_seq_id=seq_id,
        entity_id=entity_id,
        ca_coords=tuple(xyz),
    )


def _make_structure(
    condition_id: str,
    seed: int,
    sample: int,
    *,
    offset: float = 0.0,
    source_path: str = "synthetic.cif",
) -> StructureData:
    """Protein chain A (10 residues) + DNA chain B (4 residues), shifted by
    `offset` along x. Chain B sits 6 A from chain A so a protein-DNA interface
    exists."""
    residues = [
        _make_residue("A", 1, i + 1, (float(i) + offset, 0.0, 0.0))
        for i in range(10)
    ]
    chain = ChainGeometry(
        chain_id="A",
        entity_id=1,
        entity_type="polymer",
        polymer_type="polypeptide(L)",
        residues=residues,
    )
    dna_residues = [
        ResidueGeometry(
            residue_name="DA",
            chain_id="B",
            auth_seq_id=j,
            entity_id=2,
            ca_coords=(float(j - 1) + offset, 6.0, 0.0),
        )
        for j in range(1, 5)
    ]
    dna_chain = ChainGeometry(
        chain_id="B",
        entity_id=2,
        entity_type="polymer",
        polymer_type="polydeoxyribonucleotide",
        residues=dna_residues,
    )
    entity = EntityGeometry(
        entity_id=1,
        entity_type="polymer",
        polymer_type="polypeptide(L)",
        description="synthetic protein",
        chains=[chain],
    )
    dna_entity = EntityGeometry(
        entity_id=2,
        entity_type="polymer",
        polymer_type="polydeoxyribonucleotide",
        description="synthetic DNA",
        chains=[dna_chain],
    )
    return StructureData(
        prediction_id=f"{condition_id}_seed-{seed}_sample-{sample}",
        condition_id=condition_id,
        seed=seed,
        sample=sample,
        source_path=Path(source_path),
        entities=[entity, dna_entity],
        parse_status="success",
    )


def _make_dataset(
    n_seeds: int = 2,
    n_samples: int = 1,
    tmp_path: Path = None,
    conditions: tuple = ("ref_cond", "target_cond"),
):
    """Build a Dataset with synthetic structures and matching seed summaries."""
    from af3_analysis.visualization_v3.model import Dataset, SeedSummary, ConditionSummary

    ref_condition, target_condition = conditions
    predictions = {}
    seeds = {}

    for cond in conditions:
        predictions[cond] = {}
        seeds[cond] = {}
        for seed in range(1, n_seeds + 1):
            predictions[cond][seed] = {}
            for sample in range(n_samples):
                # target structures are shifted so RMSD > 0 deterministically
                offset = 0.0 if cond == ref_condition else 1.5
                predictions[cond][seed][sample] = _make_structure(
                    cond, seed, sample, offset=offset,
                    source_path=str(tmp_path) if tmp_path else "synthetic.cif",
                )
            seeds[cond][seed] = SeedSummary(
                condition_id=cond,
                seed=seed,
                n_samples=n_samples,
                metrics={
                    "pLDDT_mean": 80.0 if cond == ref_condition else 75.0,
                },
            )

    condition_summaries = {
        cond: ConditionSummary(
            condition_id=cond,
            condition_name=cond,
            n_seeds=n_seeds,
            n_predictions=n_seeds * n_samples,
        )
        for cond in conditions
    }

    return Dataset(
        name="synthetic",
        conditions=condition_summaries,
        predictions=predictions,
        seeds=seeds,
        reference_condition=ref_condition,
    )


# ---------------------------------------------------------------------------
# Unit tests: helpers
# ---------------------------------------------------------------------------

class TestV3RunnerHelpers:
    def test_iter_predictions_sorted_and_complete(self, tmp_path):
        dataset = _make_dataset(n_seeds=2, n_samples=2, tmp_path=tmp_path)
        ids = [
            (c, s, m)
            for c, s, m, _ in v3_runner._iter_predictions(dataset)
        ]
        assert len(ids) == 8
        # sorted by condition then seed then sample
        assert ids == sorted(ids)
        assert ids[0] == ("ref_cond", 1, 0)
        assert ids[-1] == ("target_cond", 2, 1)

    def test_get_reference_structure_matched_sample(self, tmp_path):
        dataset = _make_dataset(tmp_path=tmp_path)
        ref = v3_runner._get_reference_structure(dataset, 2, 0, "ref_cond")
        assert ref is not None
        assert ref.condition_id == "ref_cond"
        assert ref.seed == 2 and ref.sample == 0
        # unmatched seed returns None (no invention)
        assert v3_runner._get_reference_structure(dataset, 99, 0, "ref_cond") is None

    def test_safe_float_converts_nonfinite(self):
        assert v3_runner._safe_float(float("nan")) is None
        assert v3_runner._safe_float(float("inf")) is None
        assert v3_runner._safe_float(1.5) == 1.5

    def test_pairwise_cache_key_deterministic(self, tmp_path):
        structures_a = [
            _make_structure("c1", 1, 0, source_path=str(tmp_path / "a.cif")),
            _make_structure("c2", 1, 0, source_path=str(tmp_path / "b.cif")),
        ]
        structures_b = list(reversed(structures_a))
        key_a = v3_runner._pairwise_cache_key(structures_a, V3Config())
        key_b = v3_runner._pairwise_cache_key(structures_b, V3Config())
        assert key_a == key_b  # order-independent

        structures_c = [
            _make_structure("c1", 1, 0, source_path=str(tmp_path / "a.cif")),
            _make_structure("c3", 1, 0, source_path=str(tmp_path / "c.cif")),
        ]
        key_c = v3_runner._pairwise_cache_key(structures_c, V3Config())
        assert key_c != key_a


# ---------------------------------------------------------------------------
# Figure data preparation
# ---------------------------------------------------------------------------

class TestFigureDataPreparation:
    def test_f02_matched_sample_pairing(self, tmp_path):
        dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        data = v3_runner._generate_figure_with_data(
            generate_f02_global_structural_difference, "F02", dataset, [], None,
            ref_resolution, V3Config(), tmp_path,
        )
        assert data["status"] == "pass"
        # 1 target condition x 2 seeds x 1 sample
        assert data["n_observations"] == 2
        assert Path(data["output_path"]).is_file()

    def test_f03_uses_matched_seed_analysis(self, tmp_path):
        dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        data = v3_runner._generate_figure_with_data(
            generate_f03_matched_seed_structural_difference, "F03", dataset,
            [], None, ref_resolution, V3Config(), tmp_path,
        )
        assert data["status"] == "pass"
        # one row per matched seed
        assert data["n_observations"] == 2

    def test_f16_preserves_missingness(self, tmp_path):
        """Rows without reference pLDDT must be excluded, not imputed."""
        dataset = _make_dataset(n_seeds=1, n_samples=1, tmp_path=tmp_path)
        # remove reference seed summaries -> no reference pLDDT
        dataset.seeds["ref_cond"] = {}
        ref_resolution = {"reference_condition": "ref_cond"}
        data = v3_runner._generate_figure_with_data(
            generate_f16_confidence_change_vs_structural_change, "F16",
            dataset, [], None, ref_resolution, V3Config(), tmp_path,
        )
        # all rows lacked delta_plddt -> no imputation -> empty input
        assert data["status"] == "skip"

    def test_f10_skips_without_configured_sites(self, tmp_path):
        dataset = _make_dataset(tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        data = v3_runner._generate_figure_with_data(
            generate_f10_local_geometry, "F10", dataset, [], None,
            ref_resolution, V3Config(), tmp_path,
        )
        assert data["status"] == "skip"
        assert "sites" in data["reason"].lower()

    def test_f11_skips_without_configured_regions(self, tmp_path):
        dataset = _make_dataset(tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        data = v3_runner._generate_figure_with_data(
            generate_f11_domain_motion, "F11", dataset, [], None,
            ref_resolution, V3Config(), tmp_path,
        )
        assert data["status"] == "skip"
        assert "regions" in data["reason"].lower()

    def test_f19_skips_without_metadata(self, tmp_path):
        dataset = _make_dataset(tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        data = v3_runner._generate_figure_with_data(
            generate_f19_factorial_effects, "F19", dataset, [], None,
            ref_resolution, V3Config(), tmp_path,
        )
        assert data["status"] == "skip"
        assert "metadata" in data["reason"].lower()

    def test_f12_f13_f14_skip_without_pairwise_matrix(self, tmp_path):
        dataset = _make_dataset(tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        generators = {
            "F12": generate_f12_structural_clustering,
            "F13": generate_f13_similarity_matrix,
            "F14": generate_f14_mds_embedding,
        }
        for fig_id, gen in generators.items():
            data = v3_runner._generate_figure_with_data(
                gen, fig_id, dataset, [], None, ref_resolution,
                V3Config(), tmp_path,
            )
            assert data["status"] == "skip", fig_id


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

class TestTableGeneration:
    def test_tables_written_with_expected_grain(self, tmp_path):
        dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        tables_dir = tmp_path / "v3" / "tables"
        tables_dir.mkdir(parents=True)

        results: dict = {"tables": {}, "warnings": [], "errors": []}
        v3_runner._generate_v3_tables(
            dataset, [], None, ref_resolution, V3Config(),
            tables_dir, results,
        )

        written = results["tables"]
        for name in [
            "structural_summary",
            "per_residue_displacement",
            "contact_changes",
            "interface_contacts",
            "confidence_geometry",
        ]:
            assert name in written, f"missing table: {name}"
            assert Path(written[name]).is_file()

        # structural_summary grain: one row per prediction
        df = pd.read_csv(written["structural_summary"])
        assert len(df) == 4  # 2 conditions x 2 seeds
        assert set(df["condition_id"]) == {"ref_cond", "target_cond"}

        # per-residue displacement grain: one row per residue per prediction
        df = pd.read_csv(written["per_residue_displacement"])
        assert len(df) == 20  # 10 residues x 2 target predictions
        assert (df["displacement"] > 0).all()

        # contact_changes grain: one row per target prediction
        df = pd.read_csv(written["contact_changes"])
        assert len(df) == 2
        assert (df["reference_condition"] == "ref_cond").all()

        # missingness preserved: rmsd columns may be absent only if calc
        # failed; here they must be present and finite
        df = pd.read_csv(written["confidence_geometry"])
        assert df["rmsd_to_reference"].notna().all()

    def test_no_tables_when_no_reference(self, tmp_path):
        dataset = _make_dataset(tmp_path=tmp_path)
        ref_resolution = {"reference_condition": None}
        tables_dir = tmp_path / "v3" / "tables"
        tables_dir.mkdir(parents=True)

        results: dict = {"tables": {}, "warnings": [], "errors": []}
        v3_runner._generate_v3_tables(
            dataset, [], None, ref_resolution, V3Config(),
            tables_dir, results,
        )
        # structural_summary is reference-independent and must still exist
        assert "structural_summary" in results["tables"]
        assert "per_residue_displacement" not in results["tables"]


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

class TestManifest:
    def test_manifest_roundtrip(self, tmp_path):
        results = {
            "status": "complete",
            "summary": {"n_figures_success": 1},
            "figures": {
                "F01": {
                    "status": "pass",
                    "output_path": "x.png",
                    "n_observations": 3,
                    "warnings": [],
                },
            },
            "structural_qc": {"n_structures": np.int64(4)},
            "warnings": [],
            "errors": [],
            "skipped": [],
        }
        metadata_dir = tmp_path / "metadata"
        metadata_dir.mkdir()
        v3_runner._write_v3_manifest(metadata_dir, results)

        manifest = json.loads((metadata_dir / "figure_manifest.json").read_text())
        assert manifest["status"] == "complete"
        assert manifest["figures"][0]["figure_id"] == "F01"
        assert manifest["structural_qc"]["n_structures"] == 4  # np.int64 serialized


# ---------------------------------------------------------------------------
# End-to-end smoke test through the runner (no raw CIFs; geometry-light)
# ---------------------------------------------------------------------------

class TestRunnerSmoke:
    def test_run_pipeline_completes_without_structures(self, tmp_path):
        """With no structures, the pipeline must complete and skip figures."""
        run_dir = tmp_path / "run"
        (run_dir / "tables").mkdir(parents=True)
        # Minimal tables so the adapter can load them
        pd.DataFrame({
            "condition_id": ["cond_001"],
            "condition_name": ["cond_001"],
            "seed": [1],
            "pLDDT_mean": [80.0],
        }).to_csv(run_dir / "tables" / "seed_aggregated.csv", index=False)
        pd.DataFrame({"condition_id": ["cond_001"]}).to_csv(
            run_dir / "tables" / "descriptive_stats.csv", index=False)

        results = v3_runner.run_v3_pipeline(run_dir)

        assert results["status"] in ("complete", "completed_with_errors")
        assert results["dataset_summary"]["n_conditions"] == 1
        # every figure must be pass/skip/failed — never unknown
        for fig_id, fig in results["figures"].items():
            assert fig["status"] in ("pass", "skip", "failed"), fig_id
        # manifest + report written
        assert (run_dir / "v3" / "metadata" / "figure_manifest.json").is_file()
        assert (run_dir / "v3" / "report" / "V3_VISUALIZATION_REPORT.md").is_file()


# ---------------------------------------------------------------------------
# Reference resolution
# ---------------------------------------------------------------------------

class TestReferenceResolution:
    def _dataset(self):
        conditions = {
            f"cond_{i:03d}": ConditionSummary(
                condition_id=f"cond_{i:03d}",
                condition_name=f"cond_{i:03d}",
                n_seeds=1,
                n_predictions=1,
                metrics={},
            )
            for i in range(1, 4)
        }
        return Dataset(
            name="toy",
            conditions=conditions,
            seeds={},
            predictions={},
            reference_condition=None,
            experiment_metadata=None,
        )

    def test_empty_config_falls_back_to_first_condition(self):
        """V3Config.reference defaults to {} — must behave like no reference."""
        dataset = self._dataset()
        result = resolve_reference(dataset, {})
        assert result["reference_condition"] == "cond_001"
        assert result["reference_strategy"] == "first_condition"
        assert len(result["paired_conditions"]) == 2

    def test_none_config_falls_back_to_first_condition(self):
        dataset = self._dataset()
        result = resolve_reference(dataset, None)
        assert result["reference_condition"] == "cond_001"

    def test_explicit_reference_used(self):
        dataset = self._dataset()
        result = resolve_reference(dataset, {"condition": "cond_002"})
        assert result["reference_condition"] == "cond_002"
        assert result["reference_strategy"] == "explicit_reference"
