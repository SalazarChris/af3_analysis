"""
Tests for the V3 cross-figure structural calculation cache.

Verifies that sharing structural calculations across figures via
``FigureDataCache`` produces results identical to computing each figure
independently, and that repeated calculations are actually reused
(hit/miss accounting). Uses a small synthetic dataset so no real AF3
data is required.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from af3_analysis.visualization_v3 import runner as v3_runner
from af3_analysis.visualization_v3.config import V3Config
from af3_analysis.visualization_v3.model import (
    ChainGeometry,
    ConditionSummary,
    Dataset,
    EntityGeometry,
    ResidueGeometry,
    SeedSummary,
    StructureData,
)
from af3_analysis.visualization_v3.figures import (
    generate_f02_global_structural_difference,
    generate_f03_matched_seed_structural_difference,
    generate_f04_per_residue_displacement,
    generate_f05_displacement_heatmap,
    generate_f06_contact_map_difference,
    generate_f07_contact_change_summary,
    generate_f08_interface_analysis,
    generate_f09_interface_change_map,
    generate_f15_confidence_geometry,
    generate_f20_structure_confidence_matrix,
)
from af3_analysis.visualization_v3.structural.rmsd import calculate_rmsd


# ---------------------------------------------------------------------------
# Synthetic dataset (same pattern as test_v3_runner.py)
# ---------------------------------------------------------------------------

def _make_structure(
    condition_id: str,
    seed: int,
    sample: int,
    *,
    offset: float = 0.0,
    wobble: float = 0.0,
    dna_y: float = 6.0,
    plddt: Optional[float] = None,
    pae: Optional[float] = None,
    contact_prob: Optional[float] = None,
    source_path: str = "synthetic.cif",
) -> StructureData:
    """Protein chain A (10 residues) + DNA chain B (4 residues), shifted by
    `offset` along x. `wobble` adds a per-residue sine distortion (RMSD is
    translation-invariant, so shape distortion is needed for nonzero,
    seed-varying RMSD). Chain B sits `dna_y` A from chain A so a protein-DNA
    interface exists; varying dna_y between conditions changes the interface
    so F08/F09 have interface changes to report."""
    residues = [
        ResidueGeometry(
            residue_name="ALA",
            chain_id="A",
            auth_seq_id=i + 1,
            entity_id=1,
            ca_coords=(float(i) + offset, wobble * float(np.sin(i)), 0.0),
        )
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
            ca_coords=(float(j - 1) + offset, dna_y, 0.0),
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
        plddt_mean=plddt,
        pae_mean=pae,
        contact_prob_mean=contact_prob,
        parse_status="success",
    )


def _make_dataset(n_seeds: int = 2, n_samples: int = 1, tmp_path: Path = None):
    ref_condition, target_condition = "ref_cond", "target_cond"
    predictions = {}
    seeds = {}
    for cond in (ref_condition, target_condition):
        predictions[cond] = {}
        seeds[cond] = {}
        for seed in range(1, n_seeds + 1):
            predictions[cond][seed] = {}
            for sample in range(n_samples):
                # Vary the rigid shift per seed so RMSD values differ across
                # predictions (F15's regression needs non-degenerate x).
                offset = (0.0 if cond == ref_condition
                          else 1.0 + 0.1 * seed + 0.01 * sample)
                # Seed-dependent shape distortion: RMSD is translation
                # invariant, so without this every comparison would give
                # RMSD ~ 0 and F15's regression would be degenerate.
                wobble = 0.0 if cond == ref_condition else 0.25 * seed
                # Target's DNA is displaced far from the protein: the
                # protein geometry (RMSD/displacement) changes by the rigid
                # shift, and the protein-DNA interface loses its contacts
                # so F08/F09 have interface changes beyond the fold-change
                # band to report.
                dna_y = 6.0 if cond == ref_condition else 20.0
                is_ref = cond == ref_condition
                predictions[cond][seed][sample] = _make_structure(
                    cond, seed, sample, offset=offset, wobble=wobble,
                    dna_y=dna_y,
                    # Seed-varying confidence so F20's regression panels
                    # have non-degenerate x.
                    plddt=(85.0 if is_ref else 80.0) - 0.25 * seed,
                    pae=(4.0 if is_ref else 6.0) + 0.1 * seed,
                    contact_prob=0.5 if is_ref else 0.45,
                    source_path=str(tmp_path) if tmp_path else "synthetic.cif",
                )
            seeds[cond][seed] = SeedSummary(
                condition_id=cond,
                seed=seed,
                n_samples=n_samples,
                metrics={"pLDDT_mean": 80.0 if cond == ref_condition else 75.0},
            )
    condition_summaries = {
        cond: ConditionSummary(
            condition_id=cond,
            condition_name=cond,
            n_seeds=n_seeds,
            n_predictions=n_seeds * n_samples,
        )
        for cond in (ref_condition, target_condition)
    }
    return Dataset(
        name="synthetic",
        conditions=condition_summaries,
        predictions=predictions,
        seeds=seeds,
        reference_condition=ref_condition,
    )


REF_RESOLUTION = {"reference_condition": "ref_cond"}


def _run_figure(fig_id, generator, dataset, figures_dir, cache=None):
    return v3_runner._generate_figure_with_data(
        generator, fig_id, dataset, [], None, REF_RESOLUTION,
        V3Config(), figures_dir, cache,
    )


# ---------------------------------------------------------------------------
# Helper-level equivalence: cached vs direct computation
# ---------------------------------------------------------------------------


def test_cached_rmsd_matches_direct_calculation(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    ref = dataset.predictions["ref_cond"][1][0]
    tgt = dataset.predictions["target_cond"][1][0]
    cfg = V3Config()

    cache = v3_runner.FigureDataCache()
    first = v3_runner._cached_rmsd(cache, tgt, ref, cfg)
    second = v3_runner._cached_rmsd(cache, tgt, ref, cfg)

    direct = calculate_rmsd(
        tgt, ref,
        alignment_atom=cfg.structure.alignment_atom,
        min_common_atoms=cfg.structure.min_common_atoms,
        min_sequence_identity=cfg.structure.min_sequence_identity,
        min_coverage=cfg.structure.minimum_coverage,
    )
    assert first == (direct["rmsd"], direct["coverage"], direct["status"])
    assert second == first
    assert cache.n_misses == 1
    assert cache.n_hits == 1


def test_cached_rmsd_without_cache_computes_directly(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    ref = dataset.predictions["ref_cond"][1][0]
    tgt = dataset.predictions["target_cond"][1][0]
    cfg = V3Config()

    value = v3_runner._cached_rmsd(None, tgt, ref, cfg)
    direct = calculate_rmsd(
        tgt, ref, alignment_atom=cfg.structure.alignment_atom,
        min_common_atoms=cfg.structure.min_common_atoms,
        min_sequence_identity=cfg.structure.min_sequence_identity,
        min_coverage=cfg.structure.minimum_coverage,
    )
    assert value == (direct["rmsd"], direct["coverage"], direct["status"])


def test_cached_displacement_matches_direct_calculation(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    ref = dataset.predictions["ref_cond"][1][0]
    tgt = dataset.predictions["target_cond"][1][0]
    cfg = V3Config()

    cache = v3_runner.FigureDataCache()
    first = v3_runner._cached_displacement(cache, ref, tgt, cfg)
    second = v3_runner._cached_displacement(cache, ref, tgt, cfg)
    direct = v3_runner.calculate_per_residue_displacement(
        ref, tgt, alignment_atom=cfg.structure.alignment_atom)

    assert first["n_valid"] == direct["n_valid"]
    assert first["coverage"] == pytest.approx(direct["coverage"])
    assert [d["displacement"] for d in first["displacements"]] == pytest.approx(
        [d["displacement"] for d in direct["displacements"]])
    assert second is first  # shared object, not a copy
    assert cache.n_hits == 1


def test_cached_contact_summary_matches_direct_calculation(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    ref = dataset.predictions["ref_cond"][1][0]
    tgt = dataset.predictions["target_cond"][1][0]
    cfg = V3Config()

    cache = v3_runner.FigureDataCache()
    first = v3_runner._cached_contact_diff_summary(cache, ref, tgt, cfg)
    second = v3_runner._cached_contact_diff_summary(cache, ref, tgt, cfg)

    ref_map = v3_runner.calculate_contact_map(
        ref, threshold=cfg.structure.contact_distance)
    tgt_map = v3_runner.calculate_contact_map(
        tgt, threshold=cfg.structure.contact_distance)
    direct = v3_runner.calculate_contact_difference(ref_map, tgt_map)

    for key in v3_runner._contact_diff_summary_keys():
        assert first[key] == direct[key], key
    assert second == first
    assert cache.n_hits == 1


def test_full_contact_diff_populates_summary_cache(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    ref = dataset.predictions["ref_cond"][1][0]
    tgt = dataset.predictions["target_cond"][1][0]
    cfg = V3Config()

    cache = v3_runner.FigureDataCache()
    diff = v3_runner._full_contact_diff(cache, ref, tgt, cfg)
    assert (ref.prediction_id, tgt.prediction_id) in cache.contact_diff_summary

    summary = v3_runner._cached_contact_diff_summary(cache, ref, tgt, cfg)
    for key in v3_runner._contact_diff_summary_keys():
        assert summary[key] == diff[key], key
    # The map-bearing dicts returned to F06 remain fresh per call.
    diff2 = v3_runner._full_contact_diff(cache, ref, tgt, cfg)
    assert diff2 is not diff
    assert diff2["n_gained"] == diff["n_gained"]


def test_cached_interfaces_keyed_per_distance(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    struct = dataset.predictions["target_cond"][1][0]
    cfg = V3Config()
    cfg_6a = replace(cfg, structure=replace(cfg.structure, contact_distance=6.0))

    cache = v3_runner.FigureDataCache()
    a1 = v3_runner._cached_interfaces(cache, struct, cfg)
    a2 = v3_runner._cached_interfaces(cache, struct, cfg)
    b1 = v3_runner._cached_interfaces(cache, struct, cfg_6a)

    assert a1 is a2          # same key -> reused object
    assert a1 is not b1      # different threshold -> separate entry
    assert len(cache.interfaces) == 2
    # Interface results equivalent to direct computation
    direct = v3_runner.find_all_interfaces(
        struct, threshold=cfg.structure.contact_distance)
    assert [i.interface_id for i in a1] == [i.interface_id for i in direct]
    assert [i.n_contacts for i in a1] == [i.n_contacts for i in direct]


def test_matched_seed_results_cached_and_share_rmsd_entries(tmp_path):
    dataset = _make_dataset(n_seeds=2, n_samples=2, tmp_path=tmp_path)
    cfg = V3Config()
    cache = v3_runner.FigureDataCache()

    m1 = v3_runner._get_matched_seed_results(dataset, REF_RESOLUTION, cfg, cache)
    n_rmsd_entries = len(cache.rmsd)
    assert n_rmsd_entries > 0
    m2 = v3_runner._get_matched_seed_results(dataset, REF_RESOLUTION, cfg, cache)

    assert m1["target_cond"] is m2["target_cond"]
    assert cache.n_hits >= 1
    assert len(cache.rmsd) == n_rmsd_entries  # no recomputation

    # The rmsd entries are keyed in calculate_rmsd argument order so F02's
    # direct path can consume them.
    tgt = dataset.predictions["target_cond"][1][0]
    ref = dataset.predictions["ref_cond"][1][0]
    key = (tgt.prediction_id, ref.prediction_id)
    assert key in cache.rmsd


def test_local_geometry_cache_keys_on_definition_and_reference(tmp_path):
    dataset = _make_dataset(tmp_path=tmp_path)
    ref = dataset.predictions["ref_cond"][1][0]
    tgt = dataset.predictions["target_cond"][1][0]

    cache = v3_runner.FigureDataCache()
    site_a = {"label": "s1", "chain": "A", "residue": 3, "radius": 5.0}
    site_b = {"label": "s1", "chain": "A", "residue": 4, "radius": 5.0}

    r1 = v3_runner._cached_local_geometry(cache, tgt, site_a, ref, "CA")
    r1_again = v3_runner._cached_local_geometry(cache, tgt, site_a, ref, "CA")
    r2 = v3_runner._cached_local_geometry(cache, tgt, site_b, ref, "CA")
    r3 = v3_runner._cached_local_geometry(cache, ref, site_a, ref, "CA")

    assert r1 is r1_again
    assert r2 is not r1       # different definition
    assert r3 is not r1       # different structure (and reference role)
    assert len(cache.local_geometry) == 3

    direct = v3_runner.calculate_local_geometry(
        tgt, site_a, reference_structure=ref, alignment_atom="CA")
    assert r1.local_rmsd == pytest.approx(direct.local_rmsd)
    assert r1.n_atoms == direct.n_atoms


# ---------------------------------------------------------------------------
# Figure-level reuse: second consumer adds hits, not misses
# ---------------------------------------------------------------------------


def test_f05_reuses_f04_displacements(tmp_path):
    dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
    figures_dir = tmp_path / "figs"
    figures_dir.mkdir()
    cache = v3_runner.FigureDataCache()

    r4 = _run_figure("F04", generate_f04_per_residue_displacement,
                     dataset, figures_dir, cache)
    assert r4["status"] == "pass"
    misses_after_f04 = cache.n_misses
    hits_after_f04 = cache.n_hits

    r5 = _run_figure("F05", generate_f05_displacement_heatmap,
                     dataset, figures_dir, cache)
    assert r5["status"] == "pass"
    assert cache.n_misses == misses_after_f04  # no new calculations
    assert cache.n_hits > hits_after_f04       # reused every comparison


def test_f07_reuses_f06_contact_comparisons(tmp_path):
    dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
    figures_dir = tmp_path / "figs"
    figures_dir.mkdir()
    cache = v3_runner.FigureDataCache()

    r6 = _run_figure("F06", generate_f06_contact_map_difference,
                     dataset, figures_dir, cache)
    assert r6["status"] == "pass"
    misses_after_f06 = cache.n_misses

    r7 = _run_figure("F07", generate_f07_contact_change_summary,
                     dataset, figures_dir, cache)
    assert r7["status"] == "pass"
    assert cache.n_misses == misses_after_f06
    assert cache.n_hits >= 2


def test_f09_reuses_f08_interface_calculations(tmp_path):
    dataset = _make_dataset(n_seeds=2, n_samples=1, tmp_path=tmp_path)
    figures_dir = tmp_path / "figs"
    figures_dir.mkdir()
    cache = v3_runner.FigureDataCache()

    r8 = _run_figure("F08", generate_f08_interface_analysis,
                     dataset, figures_dir, cache)
    assert r8["status"] == "pass"
    misses_after_f08 = cache.n_misses

    r9 = _run_figure("F09", generate_f09_interface_change_map,
                     dataset, figures_dir, cache)
    assert r9["status"] == "pass"
    assert cache.n_misses == misses_after_f08
    assert cache.n_hits >= 1  # per-structure interface lists reused


def test_f02_reuses_matched_seed_rmsds(tmp_path):
    dataset = _make_dataset(n_seeds=2, n_samples=2, tmp_path=tmp_path)
    figures_dir = tmp_path / "figs"
    figures_dir.mkdir()
    cfg = V3Config()
    cache = v3_runner.FigureDataCache()

    v3_runner._get_matched_seed_results(dataset, REF_RESOLUTION, cfg, cache)
    misses_after_matched = cache.n_misses

    r02 = _run_figure("F02", generate_f02_global_structural_difference,
                      dataset, figures_dir, cache)
    assert r02["status"] == "pass"
    assert cache.n_misses == misses_after_matched
    assert cache.n_hits >= 2


# ---------------------------------------------------------------------------
# End-to-end equivalence: full figure sweep + tables with shared cache
# must produce byte-identical tables to a fully independent run
# ---------------------------------------------------------------------------

_FIGURE_SWEEP = [
    ("F02", generate_f02_global_structural_difference),
    ("F03", generate_f03_matched_seed_structural_difference),
    ("F04", generate_f04_per_residue_displacement),
    ("F05", generate_f05_displacement_heatmap),
    ("F06", generate_f06_contact_map_difference),
    ("F07", generate_f07_contact_change_summary),
    ("F08", generate_f08_interface_analysis),
    ("F09", generate_f09_interface_change_map),
    ("F15", generate_f15_confidence_geometry),
    ("F20", generate_f20_structure_confidence_matrix),
]


def test_full_run_tables_identical_with_and_without_cache(tmp_path):
    dataset = _make_dataset(n_seeds=2, n_samples=2, tmp_path=tmp_path)
    cfg = V3Config()

    dir_cached = tmp_path / "cached"
    dir_independent = tmp_path / "independent"
    figs_cached = tmp_path / "figs_cached"
    figs_independent = tmp_path / "figs_independent"
    for d in (dir_cached, dir_independent, figs_cached, figs_independent):
        d.mkdir()

    # Shared-cache run: figure sweep populates the cache, tables reuse it.
    cache = v3_runner.FigureDataCache()
    for fig_id, gen in _FIGURE_SWEEP:
        result = _run_figure(fig_id, gen, dataset, figs_cached, cache)
        assert result["status"] == "pass", fig_id
    results_cached: dict = {"tables": {}, "warnings": [], "errors": []}
    v3_runner._generate_v3_tables(
        dataset, [], None, REF_RESOLUTION, cfg, dir_cached,
        results_cached, cache,
    )

    # Independent run: no cache anywhere.
    for fig_id, gen in _FIGURE_SWEEP:
        _run_figure(fig_id, gen, dataset, figs_independent, None)
    results_independent: dict = {"tables": {}, "warnings": [], "errors": []}
    v3_runner._generate_v3_tables(
        dataset, [], None, REF_RESOLUTION, cfg, dir_independent,
        results_independent, None,
    )

    assert set(results_cached["tables"]) == set(results_independent["tables"])
    for name, path_a in results_cached["tables"].items():
        path_b = results_independent["tables"][name]
        assert Path(path_a).read_bytes() == Path(path_b).read_bytes(), name


def test_cache_summary_shape():
    cache = v3_runner.FigureDataCache()
    summary = cache.as_summary()
    assert summary["hits"] == 0
    assert summary["misses"] == 0
    for key in (
        "rmsd_entries", "displacement_entries",
        "contact_diff_summary_entries", "interface_entries",
        "matched_seed_entries", "local_geometry_entries",
    ):
        assert summary[key] == 0
