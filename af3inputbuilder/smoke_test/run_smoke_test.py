#!/usr/bin/env python3
"""
AF3 Smoke Test Runner
=====================

Minimal end-to-end smoke test of the condition_manifest → JobBuilder →
AF3 JSON pipeline.

Generates JSON files using EXISTING pipeline, validates them, and records
results.  Does NOT submit to AF3 automatically — that requires manual
submission or a separate AF3 client.

Usage::

    python smoke_test/run_smoke_test.py

Or programmatically::

    from smoke_test.run_smoke_test import run_smoke_tests
    results = run_smoke_tests(seed=1)
"""

from __future__ import annotations

import csv
import json
import os
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Bootstrap path
_HERE = Path(__file__).resolve().parent          # smoke_test/
_AF3BUILDER = _HERE.parent                       # af3inputbuilder/
_REPO_ROOT = _AF3BUILDER.parent                  # FINALVERSIONTHESIS/
if str(_AF3BUILDER) not in sys.path:
    sys.path.insert(0, str(_AF3BUILDER))

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
from af3_builder.utils.io import save_json
from af3_builder.validation.validator import AF3Validator, ValidationError


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

POU_DIR = _REPO_ROOT / "testdata" / "pou2" / "registries"
PX_DIR = _REPO_ROOT / "testdata" / "protein_x" / "registries"
OUTPUT_DIR = _HERE / "json"
RESULTS_DIR = _HERE / "results"


# ---------------------------------------------------------------------------
# Smoke-test matrix definition
# ---------------------------------------------------------------------------

@dataclass
class SmokeTestEntry:
    """A single row in the smoke-test matrix."""
    dataset: str
    condition_id: str
    category: str
    description: str
    seed: int = 1
    # Populated during execution
    generated_json: str = ""
    entity_count: int = 0
    protein_count: int = 0
    dna_count: int = 0
    ligand_count: int = 0
    ion_count: int = 0
    mod_count: int = 0
    bond_count: int = 0
    prevalidation_status: str = "NOT_GENERATED"
    execution_status: str = "NOT_ATTEMPTED"
    af3_status: str = "NOT_ATTEMPTED"
    error_type: str = ""
    error_message: str = ""
    timestamp: str = ""


# The 8 mandatory categories, mapped to EXISTING conditions
SMOKE_TEST_MATRIX: List[SmokeTestEntry] = [
    # POU/OCT4 dataset
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_baseline",
        category="A",
        description="Unmodified protein (baseline)",
    ),
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_tpo101",
        category="B",
        description="Protein + 1 verified PTM (pTPO101)",
    ),
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_tpo101_sep102",
        category="C",
        description="Protein + 2 verified PTMs (pTPO101 + pSEP102)",
    ),
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_dna",
        category="D",
        description="Protein + DNA",
    ),
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_tpo101_dna",
        category="E",
        description="Protein + PTM + DNA",
    ),
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_tpo101_sep102_dna",
        category="E_extended",
        description="Protein + 2 PTMs + DNA (bonus)",
    ),
    SmokeTestEntry(
        dataset="pou2",
        condition_id="pou_baseline",
        category="G",
        description="Protein + ion (MG, from baseline)",
    ),
    # Protein X dataset
    SmokeTestEntry(
        dataset="protein_x",
        condition_id="kinase_baseline",
        category="A_px",
        description="Protein-X baseline (genericity check)",
    ),
    SmokeTestEntry(
        dataset="protein_x",
        condition_id="kinase_ligand_A",
        category="F",
        description="Protein + ligand (ATP)",
    ),
    SmokeTestEntry(
        dataset="protein_x",
        condition_id="kinase_phospho",
        category="B_px",
        description="Protein-X + PTM (phospho_K42)",
    ),
    SmokeTestEntry(
        dataset="protein_x",
        condition_id="kinase_phospho_ligand",
        category="H",
        description="Protein + PTM + ligand",
    ),
    SmokeTestEntry(
        dataset="protein_x",
        condition_id="kinase_multi",
        category="H_extended",
        description="Protein + 2 PTMs + ligand (bonus)",
    ),
]


# ---------------------------------------------------------------------------
# Registry loading (load once)
# ---------------------------------------------------------------------------

def _load_registries(dataset: str) -> Dict[str, Any]:
    """Load all registries for a dataset."""
    if dataset == "pou2":
        d = POU_DIR
    elif dataset == "protein_x":
        d = PX_DIR
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    return {
        "protein_registry": load_protein_registry(d / "protein_registry.csv"),
        "construct_registry": load_construct_registry(d / "construct_registry.csv"),
        "modification_registry": load_modification_registry(d / "modification_registry.csv"),
        "nucleic_acid_registry": load_nucleic_acid_registry(d / "nucleic_acid_registry.csv"),
        "ligand_registry": load_ligand_registry(d / "ligand_registry.csv"),
        "ion_registry": load_ion_registry(d / "ion_registry.csv"),
        "af3_compatibility_registry": load_af3_compatibility_registry(d / "af3_compatibility_registry.csv"),
        "covalent_bond_registry": load_covalent_bond_registry(d / "covalent_bond_registry.csv"),
    }


def _load_manifest(dataset: str):
    """Load the master manifest for a dataset."""
    if dataset == "pou2":
        d = POU_DIR
    elif dataset == "protein_x":
        d = PX_DIR
    else:
        raise ValueError(f"Unknown dataset: {dataset}")

    return load_master_manifest(
        d / "master_condition_manifest.csv",
        modifications_path=d / "condition_modifications.csv",
        entities_path=d / "condition_entities.csv",
        factors_path=d / "condition_factors.csv",
    )


# ---------------------------------------------------------------------------
# Core smoke-test logic
# ---------------------------------------------------------------------------

def _run_single_test(
    entry: SmokeTestEntry,
    manifest,
    registries: Dict[str, Any],
) -> SmokeTestEntry:
    """Run a single smoke-test entry through the full pipeline."""
    entry.timestamp = datetime.now(timezone.utc).isoformat()
    cid = entry.condition_id
    seed = entry.seed

    # Step 1: Resolve condition
    try:
        spec = resolve_condition(manifest, cid, **registries)
    except Exception as e:
        entry.prevalidation_status = "RESOLUTION_FAILED"
        entry.error_type = "PIPELINE_BUG"
        entry.error_message = str(e)
        return entry

    # Step 2: Pre-submission validation
    try:
        _validate_spec_for_build(spec, allow_uncertain=False)
    except ValueError as e:
        entry.prevalidation_status = "PRE_SUBMISSION_VALIDATION_FAILED"
        entry.error_type = "REPRESENTATION_LIMITATION"
        entry.error_message = str(e)
        return entry

    entry.prevalidation_status = "PASSED"

    # Step 3: Build job
    try:
        jb = build_job(manifest, cid, seeds=[seed], **registries)
    except Exception as e:
        entry.prevalidation_status = "BUILD_FAILED"
        entry.error_type = "PIPELINE_BUG"
        entry.error_message = str(e)
        return entry

    # Step 4: Serialize
    try:
        job_dict = jb.to_dict()
    except Exception as e:
        entry.error_type = "PIPELINE_BUG"
        entry.error_message = f"to_dict() failed: {e}"
        return entry

    # Step 5: JSON inspection — verify no silent component loss
    sequences = job_dict.get("sequences", [])
    entry.entity_count = len(sequences)

    for s in sequences:
        key = list(s.keys())[0]
        if key == "protein":
            entry.protein_count += 1
            prot = s["protein"]
            mods = prot.get("modifications", [])
            entry.mod_count += len(mods)
        elif key == "dna":
            entry.dna_count += 1
        elif key == "ligand":
            entry.ligand_count += 1
        elif key == "rna":
            pass  # count as entity but not separately tracked here

    entry.bond_count = len(job_dict.get("bondedAtomPairs", []))

    # Step 6: AF3Validator
    try:
        AF3Validator.validate_job(job_dict, require_files=False)
        validation_ok = True
    except ValidationError as e:
        validation_ok = False
        entry.error_type = "JSON_FORMAT_ERROR"
        entry.error_message = "; ".join(e.messages)

    if not validation_ok:
        entry.prevalidation_status = "AF3_VALIDATOR_FAILED"
        return entry

    # Step 7: Write JSON
    json_filename = f"{cid}_seed{seed}.json"
    json_path = OUTPUT_DIR / json_filename
    try:
        save_json(str(json_path), job_dict)
        entry.generated_json = str(json_path.relative_to(_HERE))
    except Exception as e:
        entry.error_type = "IO_ERROR"
        entry.error_message = str(e)
        return entry

    # Step 8: AF3 execution status
    # No automated submission mechanism exists.
    entry.execution_status = "GENERATED_NOT_SUBMITTED"
    entry.af3_status = "AF3_EXECUTION_NOT_ATTEMPTED"

    return entry


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_smoke_tests(
    seed: int = 1,
    matrix: Optional[List[SmokeTestEntry]] = None,
) -> List[SmokeTestEntry]:
    """Run the full smoke-test matrix.

    Parameters
    ----------
    seed : int
        Seed value for all jobs.
    matrix : list, optional
        Override the default smoke-test matrix.

    Returns
    -------
    list of SmokeTestEntry
        Results for each tested condition.
    """
    if matrix is None:
        matrix = SMOKE_TEST_MATRIX

    # Ensure output dirs exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Cache loaded registries/manifests per dataset
    cache: Dict[str, Tuple[Any, Dict[str, Any]]] = {}

    results: List[SmokeTestEntry] = []

    print(f"\n{'='*60}")
    print(f"  AF3 Smoke Test")
    print(f"  Seed: {seed}")
    print(f"  Conditions: {len(matrix)}")
    print(f"  Output: {OUTPUT_DIR}")
    print(f"{'='*60}\n")

    for entry in matrix:
        entry.seed = seed
        ds = entry.dataset

        # Load data (cached per dataset)
        if ds not in cache:
            print(f"  Loading {ds} registries...")
            cache[ds] = (_load_manifest(ds), _load_registries(ds))

        manifest, registries = cache[ds]

        print(f"  [{entry.category}] {entry.condition_id} ({ds})...", end=" ")
        result = _run_single_test(entry, manifest, registries)
        results.append(result)

        status = result.prevalidation_status
        if status == "PASSED":
            print(f"OK  (json={result.generated_json})")
        else:
            print(f"FAIL  ({status}: {result.error_message[:60]})")

    # Write results CSV
    csv_path = RESULTS_DIR / "smoke_test_results.csv"
    _write_results_csv(results, csv_path)

    # Print summary
    print(f"\n{'='*60}")
    print(f"  Summary")
    print(f"{'='*60}")

    n_total = len(results)
    n_generated = sum(1 for r in results if r.generated_json)
    n_failed = sum(1 for r in results if not r.generated_json)

    print(f"  Total:    {n_total}")
    print(f"  Generated: {n_generated}")
    print(f"  Failed:   {n_failed}")
    print(f"  Results:  {csv_path}")

    if n_failed:
        print(f"\n  Failures:")
        for r in results:
            if not r.generated_json:
                print(f"    [{r.category}] {r.condition_id}: "
                      f"{r.error_type} — {r.error_message[:80]}")

    print(f"\n  AF3 execution: NOT_ATTEMPTED (no automated submission)")
    print(f"  To submit manually, use the JSON files in: {OUTPUT_DIR}")
    print(f"{'='*60}\n")

    return results


def _write_results_csv(
    results: List[SmokeTestEntry],
    path: Path,
) -> None:
    """Write results to CSV."""
    fieldnames = [
        "dataset", "condition_id", "category", "description", "seed",
        "generated_json", "entity_count", "protein_count", "dna_count",
        "ligand_count", "ion_count", "mod_count", "bond_count",
        "prevalidation_status", "execution_status", "af3_status",
        "error_type", "error_message", "timestamp",
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            row = {k: getattr(r, k, "") for k in fieldnames}
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_smoke_tests(seed=1)
