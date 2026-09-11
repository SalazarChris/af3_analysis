"""Tests for V3 automatic site/region detection (structural/detection.py)
and its integration into the runner (F10/F11 effective definitions)."""

from dataclasses import replace
from pathlib import Path

import pytest

from af3_analysis.visualization_v3 import runner as v3_runner
from af3_analysis.visualization_v3.config import V3Config
from af3_analysis.visualization_v3.model import (
    ChainGeometry,
    Dataset,
    EntityGeometry,
    ResidueGeometry,
    StructureData,
)
from af3_analysis.visualization_v3.structural.detection import (
    detect_contact_regions,
    detect_displacement_sites,
)
from af3_analysis.visualization_v3.figures import (
    generate_f10_local_geometry,
    generate_f11_domain_motion,
)


# ---------------------------------------------------------------------------
# Synthetic structures
# ---------------------------------------------------------------------------

def _linear_chain_residues(n: int, bump_residues=(), bump: float = 0.0):
    """n CA residues spaced 1 A along x; optional subset shifted by `bump`."""
    residues = []
    for i in range(n):
        seq_id = i + 1
        x = float(i) + (bump if seq_id in bump_residues else 0.0)
        residues.append(ResidueGeometry(
            residue_name="ALA",
            chain_id="A",
            auth_seq_id=seq_id,
            entity_id=1,
            ca_coords=(x, 0.0, 0.0),
        ))
    return residues


def _protein_structure(residues, condition_id, seed, sample, tmp_path):
    chain = ChainGeometry(
        chain_id="A",
        entity_id=1,
        entity_type="polymer",
        polymer_type="polypeptide(L)",
        residues=list(residues),
    )
    entity = EntityGeometry(
        entity_id=1,
        entity_type="polymer",
        polymer_type="polypeptide(L)",
        description="synthetic protein",
        chains=[chain],
    )
    return StructureData(
        prediction_id=f"{condition_id}_seed-{seed}_sample-{sample}",
        condition_id=condition_id,
        seed=seed,
        sample=sample,
        source_path=Path(tmp_path / "synthetic.cif"),
        entities=[entity],
        parse_status="success",
    )


def _bump_dataset(tmp_path, n_seeds=2, bump=2.0, bump_residues=(5, 6)):
    """Two conditions; target has a localized displacement at bump_residues."""
    preds = {}
    for cond, shift in (("ref_cond", 0.0), ("target_cond", bump)):
        preds[cond] = {}
        for seed in range(1, n_seeds + 1):
            preds[cond][seed] = {
                1: _protein_structure(
                    _linear_chain_residues(10, bump_residues, shift),
                    cond, seed, 1, tmp_path,
                )
            }
    return Dataset(name="test", predictions=preds)


# ---------------------------------------------------------------------------
# Site detection
# ---------------------------------------------------------------------------

class TestDetectDisplacementSites:
    def test_finds_localized_displacement(self, tmp_path):
        dataset = _bump_dataset(tmp_path, bump=2.0)
        sites = detect_displacement_sites(dataset, "ref_cond")
        assert len(sites) == 1
        site = sites[0]
        assert site["chain"] == "A"
        assert site["residue"] in (5, 6)
        assert site["mean_displacement"] == pytest.approx(2.0)
        assert site["radius"] == pytest.approx(8.0)
        assert site["label"] == "site_1"

    def test_uniform_displacement_yields_no_sites(self, tmp_path):
        # Rigid translation: every residue displaced equally -> no
        # concentration -> no sites.
        dataset = _bump_dataset(tmp_path, bump=0.0)
        # Give target a uniform +2 shift on all residues
        for seed in dataset.predictions["target_cond"]:
            dataset.predictions["target_cond"][seed][1] = _protein_structure(
                [ResidueGeometry(
                    residue_name="ALA", chain_id="A", auth_seq_id=i + 1,
                    entity_id=1, ca_coords=(float(i) + 2.0, 0.0, 0.0),
                ) for i in range(10)],
                "target_cond", seed, 1, tmp_path,
            )
        sites = detect_displacement_sites(dataset, "ref_cond")
        assert sites == []

    def test_empty_without_reference(self, tmp_path):
        dataset = _bump_dataset(tmp_path)
        assert detect_displacement_sites(dataset, None) == []

    def test_separates_distinct_peaks_and_orders_by_score(self, tmp_path):
        # Two separated bumps with different magnitudes: residues 2-3 (1.5 A)
        # and residues 7-8 (3.0 A). Requires > 10 residues.
        preds = {}
        for cond in ("ref_cond", "target_cond"):
            preds[cond] = {}
            for seed in (1, 2):
                if cond == "ref_cond":
                    residues = _linear_chain_residues(12)
                else:
                    residues = _linear_chain_residues(
                        12, bump_residues=(2, 3), bump=1.5)
                    residues = [
                        r if r.auth_seq_id not in (7, 8)
                        else ResidueGeometry(
                            residue_name=r.residue_name,
                            chain_id=r.chain_id,
                            auth_seq_id=r.auth_seq_id,
                            entity_id=r.entity_id,
                            ca_coords=(r.ca_coords[0] + 3.0, 0.0, 0.0),
                        )
                        for r in residues
                    ]
                preds[cond][seed] = {
                    1: _protein_structure(residues, cond, seed, 1, tmp_path)
                }
        dataset = Dataset(name="test", predictions=preds)
        sites = detect_displacement_sites(dataset, "ref_cond")
        assert len(sites) == 2
        # Higher-displacement peak labeled first
        assert sites[0]["residue"] in (7, 8)
        assert sites[0]["mean_displacement"] == pytest.approx(3.0)
        assert sites[1]["residue"] in (2, 3)
        assert sites[1]["mean_displacement"] == pytest.approx(1.5)

    def test_degenerate_all_tied_distribution(self, tmp_path):
        # One bump in a sea of exact zeros: MAD == 0 degenerate case; the
        # differing residues must still be flagged.
        dataset = _bump_dataset(tmp_path, bump=1.0, bump_residues=(3,))
        sites = detect_displacement_sites(dataset, "ref_cond")
        assert len(sites) == 1
        assert sites[0]["residue"] == 3
        assert sites[0]["mean_displacement"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Region detection
# ---------------------------------------------------------------------------

def _two_cluster_residues():
    """Two 10-residue clusters 100 A apart; each cluster internally
    contact-linked (non-adjacent residues within 8 A)."""
    residues = []
    for cluster_offset in (0.0, 100.0):
        for k in range(10):
            seq_id = len(residues) + 1
            # Interleaved pattern so non-adjacent residues stay in contact
            x = cluster_offset + (k % 2) * 2.0
            y = (k // 2) * 2.5
            residues.append(ResidueGeometry(
                residue_name="ALA", chain_id="A", auth_seq_id=seq_id,
                entity_id=1, ca_coords=(x, y, 0.0),
            ))
    return residues


class TestDetectContactRegions:
    def test_finds_separated_clusters(self, tmp_path):
        structure = _protein_structure(
            _two_cluster_residues(), "ref_cond", 1, 1, tmp_path)
        regions = detect_contact_regions(structure)
        assert len(regions) == 2
        spans = sorted((r["start"], r["end"]) for r in regions)
        assert spans == [(1, 10), (11, 20)]
        for region in regions:
            assert region["chain"] == "A"
            assert region["n_residues"] == 10
            assert len(region["seq_list"]) == 10
            assert region["label"] in ("region_1", "region_2")

    def test_dense_chain_yields_single_region(self, tmp_path):
        # 1 A spacing: many non-adjacent residues fall within the 8 A
        # threshold, so the whole chain forms one contact-linked region.
        structure = _protein_structure(
            _linear_chain_residues(6), "ref_cond", 1, 1, tmp_path)
        regions = detect_contact_regions(structure)
        assert len(regions) == 1
        assert (regions[0]["start"], regions[0]["end"]) == (1, 6)

    def test_empty_structure_is_safe(self, tmp_path):
        structure = _protein_structure([], "ref_cond", 1, 1, tmp_path)
        assert detect_contact_regions(structure) == []


# ---------------------------------------------------------------------------
# Runner integration
# ---------------------------------------------------------------------------

class TestRunnerIntegration:
    def test_f10_auto_detected_sites_run(self, tmp_path):
        dataset = _bump_dataset(tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        config = replace(V3Config(), auto_sites=True)
        data = v3_runner._generate_figure_with_data(
            generate_f10_local_geometry, "F10", dataset, [], None,
            ref_resolution, config, tmp_path,
        )
        assert data["status"] == "pass"
        assert data["n_observations"] > 0
        # Provenance must be explicit in the result warnings (which reach
        # the manifest)
        warnings_text = " ".join(data.get("warnings") or [])
        assert "predicted structural" in warnings_text
        assert "not functional" in warnings_text

    def test_f11_auto_detected_regions_run(self, tmp_path):
        dataset = _bump_dataset(tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        config = replace(V3Config(), auto_regions=True)
        data = v3_runner._generate_figure_with_data(
            generate_f11_domain_motion, "F11", dataset, [], None,
            ref_resolution, config, tmp_path,
        )
        assert data["status"] == "pass"
        assert data["n_observations"] > 0
        warnings_text = " ".join(data.get("warnings") or [])
        assert "predicted structural" in warnings_text
        assert "not functional" in warnings_text

    def test_f10_f11_skip_without_manual_or_auto(self, tmp_path):
        dataset = _bump_dataset(tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        for fig_id, generator in (("F10", generate_f10_local_geometry),
                                  ("F11", generate_f11_domain_motion)):
            data = v3_runner._generate_figure_with_data(
                generator, fig_id, dataset, [], None,
                ref_resolution, V3Config(), tmp_path,
            )
            assert data["status"] == "skip"

    def test_manual_definitions_take_priority(self, tmp_path):
        dataset = _bump_dataset(tmp_path)
        ref_resolution = {"reference_condition": "ref_cond"}
        manual = [{"label": "manual_site", "chain": "A", "residue": 1,
                   "radius": 5.0}]
        sites = v3_runner.effective_sites(
            dataset, ref_resolution, replace(V3Config(), sites=manual,
                                             auto_sites=True))
        assert sites == manual

    def test_effective_sites_empty_when_auto_disabled(self, tmp_path):
        dataset = _bump_dataset(tmp_path)
        assert v3_runner.effective_sites(
            dataset, {"reference_condition": "ref_cond"}, V3Config()) == []
        assert v3_runner.effective_regions(
            dataset, [], {"reference_condition": "ref_cond"}, V3Config()) == []
