#!/usr/bin/env python3
"""
AF3 Smoke Test Runner
=====================

Minimal end-to-end smoke test of the condition_manifest → JobBuilder →
AF3 JSON pipeline.

Conditions are DISCOVERED from the master condition manifests under the
registries root — nothing about any specific experiment is hard-coded here.
Generates JSON files using the EXISTING pipeline, validates them, and
records results.  Does NOT submit to AF3 automatically — that requires
manual submission or a separate AF3 client.

Usage::

    python smoke_test/run_smoke_test.py [--registries-root DIR] [--seed N]
                                        [--conditions ID1,ID2,...]

Or programmatically::

    from smoke_test.run_smoke_test import run_smoke_tests
    results = run_smoke_tests(seed=1)
"""

from __future__ import annotations

import argparse
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
_REPO_ROOT = _AF3BUILDER.parent                  # repository root
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

# Registry manifests live in per-dataset subdirectories of this root.
# The root itself is data; this script contains no dataset names.
DEFAULT_REGISTRIES_ROOT = _REPO_ROOT / "testdata" / "registries"
OUTPUT_DIR = _HERE / "json"
RESULTS_DIR = _HERE / "results"


# ---------------------------------------------------------------------------
# Smoke-test entry
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


def _registry_dirs(registries_root: Path) -> List[Path]:
    """Enumerate registry directories under the root.

    Supports both layouts without any dataset names:
    - flat: <root>/master_condition_manifest.csv
    - nested: <root>/<dataset>/master_condition_manifest.csv and
      <root>/<dataset>/registries/master_condition_manifest.csv
    """
    if not registries_root.is_dir():
        raise FileNotFoundError(
            f"Registries root not found: {registries_root}"
        )
    found: Dict[Path, None] = {}
    if (registries_root / "master_condition_manifest.csv").is_file():
        found[registries_root.resolve()] = None
    # Bounded scan: one and two levels below the root.
    for pattern in ("*/master_condition_manifest.csv",
                    "*/*/master_condition_manifest.csv"):
        for manifest_path in registries_root.glob(pattern):
            found[manifest_path.parent.resolve()] = None
    dirs = sorted(found.keys())
    if not dirs:
        raise FileNotFoundError(
            f"No directories containing "
            f"'master_condition_manifest.csv' found under {registries_root}"
        )
    return dirs


def _dataset_label(registry_dir: Path, registries_root: Path) -> str:
    """Stable display label for a registry directory (data-derived)."""
    if registry_dir.name != "registries" and registry_dir != registries_root.resolve():
        return registry_dir.name
    return registry_dir.parent.name


def _discover_entries(
    registries_root: Path,
    condition_filter: Optional[List[str]] = None,
) -> Tuple[List[SmokeTestEntry], Dict[str, Path]]:
    """Discover smoke-test entries from manifest CSVs.

    Every condition in every discovered manifest gets one entry. No
    condition names, categories, or descriptions are hard-coded; category
    labels are generated positionally and are stable only within one run.

    Returns (entries, registry_dirs_by_label).
    """
    entries: List[SmokeTestEntry] = []
    dirs_by_label: Dict[str, Path] = {}
    for registry_dir in _registry_dirs(registries_root):
        label = _dataset_label(registry_dir, registries_root)
        dirs_by_label[label] = registry_dir
        manifest = load_master_manifest(
            registry_dir / "master_condition_manifest.csv",
            modifications_path=registry_dir / "condition_modifications.csv",
            entities_path=registry_dir / "condition_entities.csv",
            factors_path=registry_dir / "condition_factors.csv",
        )
        for k, cid in enumerate(manifest.condition_ids, start=1):
            if condition_filter and cid not in condition_filter:
                continue
            entries.append(SmokeTestEntry(
                dataset=label,
                condition_id=cid,
                category=f"C{k:02d}",
                description="",
            ))
    if condition_filter:
        discovered_ids = {e.condition_id for e in entries}
        missing = sorted(set(condition_filter) - discovered_ids)
        if missing:
            raise ValueError(
                f"Requested condition(s) not found in any manifest: {missing}"
            )
    return entries, dirs_by_label


# ---------------------------------------------------------------------------
# Registry loading (per dataset directory)
# ---------------------------------------------------------------------------

def _load_registries(registry_dir: Path) -> Dict[str, Any]:
    """Load all registries from one dataset registry directory."""
    return {
        "protein_registry": load_protein_registry(registry_dir / "protein_registry.csv"),
        "construct_registry": load_construct_registry(registry_dir / "construct_registry.csv"),
        "modification_registry": load_modification_registry(registry_dir / "modification_registry.csv"),
        "nucleic_acid_registry": load_nucleic_acid_registry(registry_dir / "nucleic_acid_registry.csv"),
        "ligand_registry": load_ligand_registry(registry_dir / "ligand_registry.csv"),
        "ion_registry": load_ion_registry(registry_dir / "ion_registry.csv"),
        "af3_compatibility_registry": load_af3_compatibility_registry(registry_dir / "af3_compatibility_registry.csv"),
        "covalent_bond_registry": load_covalent_bond_registry(registry_dir / "covalent_bond_registry.csv"),
    }


def _load_manifest(registry_dir: Path):
    """Load the master manifest for one dataset registry directory."""
    return load_master_manifest(
        registry_dir / "master_condition_manifest.csv",
        modifications_path=registry_dir / "condition_modifications.csv",
        entities_path=registry_dir / "condition_entities.csv",
        factors_path=registry_dir / "condition_factors.csv",
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
    registries_root: Path = None,
    condition_filter: Optional[List[str]] = None,
) -> List[SmokeTestEntry]:
    """Run the smoke-test matrix.

    Parameters
    ----------
    seed : int
        Seed value for all jobs.
    matrix : list, optional
        Explicit entries. When omitted, entries are discovered from the
        registries root.
    registries_root : Path, optional
        Root containing per-dataset registry directories. Defaults to
        DEFAULT_REGISTRIES_ROOT.
    condition_filter : list of str, optional
        Restrict to these condition ids (discovery mode only).

    Returns
    -------
    list of SmokeTestEntry
        Results for each tested condition.
    """
    if matrix is None:
        matrix, dirs_by_label = _discover_entries(
            registries_root or DEFAULT_REGISTRIES_ROOT,
            condition_filter,
        )
    else:
        dirs_by_label = {}

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
            registry_dir = dirs_by_label.get(ds)
            if registry_dir is None:
                registry_dir = (
                    (registries_root or DEFAULT_REGISTRIES_ROOT) / ds
                )
            print(f"  Loading {ds} registries...")
            cache[ds] = (_load_manifest(registry_dir), _load_registries(registry_dir))

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

def _main() -> int:
    parser = argparse.ArgumentParser(
        description="AF3 condition-manifest smoke test (conditions are "
                    "discovered from manifest CSVs)"
    )
    parser.add_argument(
        "--registries-root",
        type=str,
        default=str(DEFAULT_REGISTRIES_ROOT),
        help="Root directory containing per-dataset registry folders "
             "(each with master_condition_manifest.csv)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Seed for generated jobs (default: 1)",
    )
    parser.add_argument(
        "--conditions",
        type=str,
        default=None,
        help="Comma-separated condition ids to restrict the run "
             "(default: all discovered conditions)",
    )
    args = parser.parse_args()

    condition_filter = None
    if args.conditions:
        condition_filter = [
            c.strip() for c in args.conditions.split(",") if c.strip()
        ]

    results = run_smoke_tests(
        seed=args.seed,
        registries_root=Path(args.registries_root),
        condition_filter=condition_filter,
    )
    return 0 if all(r.generated_json for r in results) else 1


if __name__ == "__main__":
    sys.exit(_main())
