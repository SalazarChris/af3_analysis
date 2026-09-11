"""
Tests for the V3 pairwise RMSD matrix fast path.

Verifies that the prepared-structure path used by
``calculate_pairwise_rmsd_matrix`` gives results identical to the general
per-pair ``calculate_rmsd`` path, including chain-selection fallbacks,
identity/coverage gates, and matrix metadata save/load round trips.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from af3_analysis.visualization_v3.model import (
    ChainGeometry,
    EntityGeometry,
    ResidueGeometry,
    StructureData,
)
from af3_analysis.visualization_v3.structural.comparability import (
    find_common_space_prepared,
    kabsch_rmsd_prepared,
    prepare_structure,
)
from af3_analysis.visualization_v3.structural.pairwise import (
    calculate_pairwise_rmsd_matrix,
    load_pairwise_matrix,
    save_pairwise_matrix,
)
from af3_analysis.visualization_v3.structural.rmsd import calculate_rmsd


def _make_structure(
    prediction_id: str,
    *,
    n_protein: int = 12,
    n_dna: int = 5,
    offset: float = 0.0,
) -> StructureData:
    """Protein chain A + optional DNA chain B, shifted by `offset` along x."""

    def _coords(i: int, y: float) -> tuple:
        return (float(i) + offset, y * np.sin(i), 0.3 * np.cos(i))

    residues = [
        ResidueGeometry(
            residue_name="ALA",
            chain_id="A",
            auth_seq_id=i + 1,
            entity_id=1,
            ca_coords=_coords(i, 0.5),
        )
        for i in range(n_protein)
    ]
    chains = [
        ChainGeometry(
            chain_id="A",
            entity_id=1,
            entity_type="polymer",
            polymer_type="polypeptide(L)",
            residues=residues,
        )
    ]
    entities = [
        EntityGeometry(
            entity_id=1,
            entity_type="polymer",
            polymer_type="polypeptide(L)",
            description=None,
            chains=chains,
        )
    ]
    if n_dna > 0:
        dna_residues = [
            ResidueGeometry(
                residue_name="DA",
                chain_id="B",
                auth_seq_id=i + 1,
                entity_id=2,
                ca_coords=_coords(i, 10.0),
            )
            for i in range(n_dna)
        ]
        chains.append(
            ChainGeometry(
                chain_id="B",
                entity_id=2,
                entity_type="polymer",
                polymer_type="polydeoxyribonucleotide",
                residues=dna_residues,
            )
        )
        entities.append(
            EntityGeometry(
                entity_id=2,
                entity_type="polymer",
                polymer_type="polydeoxyribonucleotide",
                description=None,
                chains=[chains[1]],
            )
        )

    return StructureData(
        prediction_id=prediction_id,
        condition_id="cond_x",
        seed=1,
        sample=1,
        source_path="synthetic.cif",
        entities=entities,
    )


def _fixture_structures():
    return [
        _make_structure("p1", offset=0.0),
        _make_structure("p2", offset=0.4),
        _make_structure("p3", offset=2.0, n_protein=10),  # partial coverage
        _make_structure("p4", n_dna=0),                   # protein-only
        _make_structure("p5", n_protein=4),               # below identity gate vs p1
    ]


# ---------------------------------------------------------------------------
# Prepared fast path vs general path equivalence
# ---------------------------------------------------------------------------


def test_prepared_path_matches_general_path_per_pair():
    structures = _fixture_structures()
    prepared = [prepare_structure(s) for s in structures]

    for i in range(len(structures)):
        for j in range(i + 1, len(structures)):
            slow = calculate_rmsd(
                structures[i],
                structures[j],
                alignment_atom="CA",
                min_common_atoms=3,
                min_sequence_identity=0.5,
                min_coverage=0.80,
            )
            status, coverage, common_residues, n_ca = find_common_space_prepared(
                prepared[i], prepared[j], min_sequence_identity=0.5
            )
            fast = (
                kabsch_rmsd_prepared(prepared[i], prepared[j], common_residues)
                if status != "not_comparable" and n_ca >= 3
                else None
            )
            if slow["rmsd"] is None:
                assert fast is None, f"pair ({i}, {j})"
            else:
                assert fast is not None, f"pair ({i}, {j})"
                assert slow["rmsd"] == pytest.approx(fast, abs=1e-9), (
                    f"pair ({i}, {j})"
                )


def test_prepared_path_identity_gate_blocks_low_identity_pairs():
    # p5 (4 residues) vs p1 (12): identity 4/12 < 0.5 -> not comparable,
    # even though 4 common CA atoms exist.
    s1 = _make_structure("p1")
    s5 = _make_structure("p5", n_protein=4)
    prepared = [prepare_structure(s1), prepare_structure(s5)]
    status, coverage, _common, n_ca = find_common_space_prepared(
        prepared[0], prepared[1], min_sequence_identity=0.5
    )
    assert status == "not_comparable"
    assert n_ca == 4
    # Even with residues in common, the gate blocks the comparison; calling
    # the RMSD step with the returned (non-empty) residues must never be
    # reached by the matrix loop. Sanity: empty space -> no RMSD.
    assert kabsch_rmsd_prepared(prepared[0], prepared[1], {}) is None


def test_matrix_is_symmetric_with_zero_diagonal():
    structures = _fixture_structures()
    matrix = calculate_pairwise_rmsd_matrix(structures, alignment_atom="CA")
    assert matrix.n_structures == len(structures)
    assert np.allclose(matrix.matrix, matrix.matrix.T, equal_nan=True)
    assert np.all(np.diag(matrix.matrix) == 0.0)
    assert np.all(np.diag(matrix.valid))


def test_matrix_excludes_pairs_below_identity_gate():
    structures = _fixture_structures()
    matrix = calculate_pairwise_rmsd_matrix(structures, alignment_atom="CA")
    # p5 (index 4) has 4 protein residues: identity 4/12 < 0.5 vs every
    # full-length structure, so all of its cross pairs are invalid. Only the
    # diagonal (self-comparison) is valid for p5.
    assert matrix.valid[4].sum() == 1
    assert matrix.valid[4, 4]
    assert not matrix.valid[:4, 4].any()
    assert not matrix.valid[4, :4].any()


def test_progress_logging_emitted(caplog):
    structures = _fixture_structures()
    with caplog.at_level(logging.INFO, logger="af3_analysis.visualization_v3.structural.pairwise"):
        calculate_pairwise_rmsd_matrix(structures, alignment_atom="CA")
    messages = [r.getMessage() for r in caplog.records]
    assert any("Computing" in m and "pairwise RMSD" in m for m in messages)
    assert any("Done" in m and "pairs in" in m for m in messages)


# ---------------------------------------------------------------------------
# Save / load round trip (cache filename consistency)
# ---------------------------------------------------------------------------


def test_save_writes_keyed_metadata_and_round_trips(tmp_path):
    structures = _fixture_structures()
    matrix = calculate_pairwise_rmsd_matrix(structures, alignment_atom="CA")

    matrix_path = save_pairwise_matrix(
        matrix, tmp_path, filename="cachekey.csv"
    )
    metadata_path = tmp_path / "cachekey_metadata.csv"
    assert matrix_path == tmp_path / "cachekey.csv"
    assert metadata_path.exists()

    loaded = load_pairwise_matrix(matrix_path, metadata_path)
    assert loaded.predictions == matrix.predictions
    assert loaded.conditions == matrix.conditions
    assert loaded.seeds == matrix.seeds
    assert loaded.samples == matrix.samples
    assert np.allclose(loaded.matrix, matrix.matrix, equal_nan=True)
    assert np.array_equal(loaded.valid, matrix.valid)


def test_round_tripped_matrix_feeds_loader_cache_path(tmp_path):
    """Runner cache contract: {key}.csv + {key}_metadata.csv round trip."""
    structures = _fixture_structures()
    matrix = calculate_pairwise_rmsd_matrix(structures, alignment_atom="CA")
    key = "abc123"
    save_pairwise_matrix(matrix, tmp_path, filename=f"{key}.csv")
    loaded = load_pairwise_matrix(
        tmp_path / f"{key}.csv", tmp_path / f"{key}_metadata.csv"
    )
    assert loaded.n_structures == matrix.n_structures
    valid_distances = loaded.matrix[loaded.valid]
    valid_distances = valid_distances[valid_distances > 0]
    if len(valid_distances):
        assert loaded.mean_distance == pytest.approx(matrix.mean_distance)
