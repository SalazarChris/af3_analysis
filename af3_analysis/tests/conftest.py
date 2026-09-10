"""
Shared fixtures for AF3 analysis pipeline regression tests.

Uses real test data from testdata/pou2/ and temporary directories for outputs.
"""

import shutil
import tempfile
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Path constants
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTDATA = REPO_ROOT / "testdata" / "pou2"


# ---------------------------------------------------------------------------
# Fixtures: input paths
# ---------------------------------------------------------------------------


@pytest.fixture
def repo_root():
    """Repository root directory."""
    return REPO_ROOT


@pytest.fixture
def pou_parent_dir():
    """
    Path to testdata/pou2/ — the correct input directory for the extraction
    script. Contains 8 condition subdirectories, each with seed subdirs.
    """
    assert TESTDATA.is_dir(), f"Test data not found: {TESTDATA}"
    return TESTDATA


# ---------------------------------------------------------------------------
# Fixtures: temporary output directories
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_output(tmp_path):
    """Temporary output directory, cleaned up after the test."""
    return tmp_path


@pytest.fixture
def extraction_output(tmp_path):
    """Temporary directory for extraction CSV outputs."""
    out = tmp_path / "tables"
    out.mkdir()
    return out


# ---------------------------------------------------------------------------
# Fixtures: synthetic data for fast unit tests
# ---------------------------------------------------------------------------


@pytest.fixture
def synthetic_confidence_dir(tmp_path):
    """
    Create a minimal synthetic input with the STRUCTURE the extraction
    script expects: parent > condition_folder > seed_subdirs.

    Returns the parent directory (pass this as input_dir to extract_metrics).
    """
    import json
    import numpy as np

    parent = tmp_path / "synthetic_input"
    cond_dir = parent / "test_condition"
    cond_dir.mkdir(parents=True)

    # Write _data.json at condition level
    data_json = {"name": "test_condition", "sequences": []}
    (cond_dir / "test_condition_data.json").write_text(json.dumps(data_json))

    # Write ranking_scores.csv at condition level
    ranking_lines = ["seed,sample,ranking_score"]
    for seed in [1, 2]:
        for sample in [0, 1]:
            ranking_lines.append(f"{seed},{sample},{0.5 + seed * 0.01 + sample * 0.001:.6f}")
    (cond_dir / "test_condition_ranking_scores.csv").write_text("\n".join(ranking_lines))

    # Write confidences files in seed subdirectories
    rng = np.random.RandomState(42)
    for seed in [1, 2]:
        for sample in [0, 1]:
            sub = cond_dir / f"seed-{seed}_sample-{sample}"
            sub.mkdir()

            conf = {
                "atom_plddts": rng.uniform(60, 95, size=100).tolist(),
                "contact_probs": rng.uniform(0, 1, size=20).tolist(),
                "pae": rng.uniform(1, 10, size=20).tolist(),
                "token_chain_ids": ["A"] * 20,
                "atom_chain_ids": ["A"] * 100,
            }
            fname = f"test_condition_seed-{seed}_sample-{sample}_confidences.json"
            (sub / fname).write_text(json.dumps(conf))

    return parent


@pytest.fixture
def synthetic_multi_condition_dir(tmp_path):
    """
    Create a synthetic input with 2 conditions, 2 seeds x 2 samples each.
    Returns the parent directory.
    """
    import json
    import numpy as np

    parent = tmp_path / "multi_cond_input"

    conditions = {
        "cond_A": {"name": "condition_alpha"},
        "cond_B": {"name": "condition_beta"},
    }

    for cond_folder, meta in conditions.items():
        cond_dir = parent / cond_folder
        cond_dir.mkdir(parents=True)

        # _data.json
        (cond_dir / f"{cond_folder}_data.json").write_text(
            json.dumps({"name": meta["name"], "sequences": []})
        )

        # ranking_scores.csv
        lines = ["seed,sample,ranking_score"]
        for seed in [1, 2]:
            for sample in [0, 1]:
                lines.append(f"{seed},{sample},{0.5 + seed * 0.01:.6f}")
        (cond_dir / f"{cond_folder}_ranking_scores.csv").write_text("\n".join(lines))

        # confidences
        rng = np.random.RandomState(hash(cond_folder) % 2**31)
        for seed in [1, 2]:
            for sample in [0, 1]:
                sub = cond_dir / f"seed-{seed}_sample-{sample}"
                sub.mkdir()
                conf = {
                    "atom_plddts": rng.uniform(60, 95, size=50).tolist(),
                    "contact_probs": rng.uniform(0, 1, size=10).tolist(),
                    "pae": rng.uniform(1, 10, size=10).tolist(),
                    "token_chain_ids": ["A"] * 10,
                    "atom_chain_ids": ["A"] * 50,
                }
                fname = f"{cond_folder}_seed-{seed}_sample-{sample}_confidences.json"
                (sub / fname).write_text(json.dumps(conf))

    return parent
