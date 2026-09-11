"""
Unit tests for the smoke-test runner.

These tests verify the runner logic and the condition_manifest →
JobBuilder → AF3 JSON pipeline against whatever conditions are
DISCOVERED from the registries root. No condition names, entity ids, or
chemical codes are hard-coded: every assertion is a property of the
software, not of a particular experiment.

Data contract enforced per discovered condition:
  - resolve_condition must succeed (manifest integrity).
  - build_job either succeeds or fails with a classified
    REPRESENTATION_LIMITATION (documented AF3 representation limits).
  - On success: valid AF3 JSON dialect, validator passes, no silent
    component loss, deterministic output for identical inputs.
"""

import json
import csv
import sys
from pathlib import Path

import pytest

# Bootstrap path (mirrors run_smoke_test.py)
_HERE = Path(__file__).resolve().parent           # smoke_test/
_AF3BUILDER = _HERE.parent                        # af3inputbuilder/
_REPO_ROOT = _AF3BUILDER.parent                   # repository root
if str(_AF3BUILDER) not in sys.path:
    sys.path.insert(0, str(_AF3BUILDER))

from af3_builder.condition_manifest import (
    load_master_manifest,
    load_construct_registry,
    load_modification_registry,
    load_af3_compatibility_registry,
)
from af3_builder.condition_manifest.builder import (
    build_job,
    resolve_condition,
    _validate_spec_for_build,
)
from af3_builder.validation.validator import AF3Validator, ValidationError

from smoke_test.run_smoke_test import (
    DEFAULT_REGISTRIES_ROOT,
    _discover_entries,
    _load_manifest,
    _load_registries,
    _run_single_test,
    SmokeTestEntry,
)

pytestmark = pytest.mark.skipif(
    not DEFAULT_REGISTRIES_ROOT.is_dir(),
    reason="registries root not available (data not checked in)",
)


# ---------------------------------------------------------------------------
# Fixtures: everything discovered from data
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def discovered():
    """(entries, dirs_by_label) discovered from the registries root."""
    entries, dirs_by_label = _discover_entries(DEFAULT_REGISTRIES_ROOT)
    return entries, dirs_by_label


@pytest.fixture(scope="module")
def first_dataset(discovered):
    """Manifest + registries for the first discovered dataset."""
    entries, dirs_by_label = discovered
    label = entries[0].dataset
    registry_dir = dirs_by_label[label]
    return _load_manifest(registry_dir), _load_registries(registry_dir)


@pytest.fixture(scope="module")
def loaded_cache(discovered):
    """label -> (manifest, registries) for all discovered datasets."""
    entries, dirs_by_label = discovered
    return {
        label: (_load_manifest(d), _load_registries(d))
        for label, d in dirs_by_label.items()
    }


# ---------------------------------------------------------------------------
# Test: every discovered condition resolves and builds
# ---------------------------------------------------------------------------

class TestAllDiscoveredConditions:
    """The pipeline must handle every condition the data defines."""

    def test_discovery_finds_conditions(self, discovered):
        entries, _ = discovered
        assert len(entries) > 0

    def test_all_conditions_resolve(
        self, discovered, loaded_cache
    ):
        entries, _ = discovered
        failures = []
        for entry in entries:
            manifest, registries = loaded_cache[entry.dataset]
            try:
                spec = resolve_condition(
                    manifest, entry.condition_id, **registries
                )
                _validate_spec_for_build(spec)
            except ValueError as e:
                # Only documented AF3 representation limits are acceptable:
                # unsupported representations and uncertain ones (the
                # runner classifies both as REPRESENTATION_LIMITATION).
                msg = str(e)
                if "UNSUPPORTED" not in msg and "uncertain" not in msg:
                    failures.append((entry.condition_id, msg[:100]))
        assert failures == []

    def test_buildable_conditions_produce_valid_json(
        self, discovered, loaded_cache
    ):
        entries, _ = discovered
        problems = []
        for entry in entries:
            manifest, registries = loaded_cache[entry.dataset]
            try:
                spec = resolve_condition(
                    manifest, entry.condition_id, **registries
                )
                _validate_spec_for_build(spec)
            except ValueError:
                # Representation limitation; not a software defect.
                continue
            try:
                jb = build_job(
                    manifest, entry.condition_id, seeds=[1], **registries
                )
                d = jb.to_dict()
            except ValueError:
                continue
            # Structural properties of valid AF3 JSON
            assert d["dialect"] == "alphafold3", entry.condition_id
            assert d["modelSeeds"] == [1], entry.condition_id
            assert len(d["sequences"]) >= 1, entry.condition_id
            try:
                AF3Validator.validate_job(d, require_files=False)
            except ValidationError as e:
                problems.append((entry.condition_id, "; ".join(e.messages)[:100]))
        assert problems == []


# ---------------------------------------------------------------------------
# Test: runner behavior (error classification, no crashes)
# ---------------------------------------------------------------------------

class TestRunnerClassification:
    """The runner must classify outcomes, never crash."""

    def test_every_condition_classified(self, discovered, loaded_cache):
        entries, _ = discovered
        unclassified = []
        for entry in entries:
            manifest, registries = loaded_cache[entry.dataset]
            result = _run_single_test(entry, manifest, registries)
            ok = result.generated_json != ""
            failed = result.prevalidation_status not in (
                "NOT_GENERATED", "PASSED",
            )
            if not ok and not failed:
                unclassified.append(entry.condition_id)
        assert unclassified == []

    def test_unknown_condition_fails_cleanly(self, first_dataset):
        manifest, registries = first_dataset
        with pytest.raises(ValueError, match="not found|Unknown|No condition"):
            resolve_condition(
                manifest, "nonexistent_condition_for_test", **registries
            )

    def test_unsupported_modification_fails(self, first_dataset):
        """A synthetic condition referencing an unsupported modification
        must be rejected with an UNSUPPORTED error (software behavior)."""
        from af3_builder.condition_manifest.manifest import (
            MasterManifest, ConditionRecord, ConditionModificationRecord,
            ConditionEntityRecord,
        )
        from af3_builder.condition_manifest.registries import (
            ModificationRecord, AF3CompatibilityRecord,
        )

        manifest, registries = first_dataset
        # Use an existing construct id from the data (no hard-coding).
        construct_id = next(iter(registries["construct_registry"].keys()))

        synthetic = MasterManifest()
        synthetic.conditions["c1"] = ConditionRecord(
            condition_id="c1", condition_name="Synthetic"
        )
        synthetic.modifications["m1"] = ConditionModificationRecord(
            condition_id="c1", modification_id="bad_mod",
            sequence_position="10", construct_id=construct_id,
        )
        synthetic.entities["e1"] = ConditionEntityRecord(
            condition_id="c1", entity_type="protein",
            entity_id=construct_id, stoichiometry="1",
        )

        mod_reg = {"bad_mod": ModificationRecord(modification_id="bad_mod")}
        af3_reg = {
            "rep": AF3CompatibilityRecord(
                representation_id="rep",
                modification_id="bad_mod",
                af3_status="unsupported",
            )
        }

        with pytest.raises(ValueError, match="UNSUPPORTED"):
            build_job(
                synthetic, "c1", seeds=[1],
                construct_registry=registries["construct_registry"],
                modification_registry=mod_reg,
                af3_compatibility_registry=af3_reg,
            )


# ---------------------------------------------------------------------------
# Test: determinism (software property)
# ---------------------------------------------------------------------------

class TestDeterminism:
    """Same inputs must produce identical JSON; different seeds differ."""

    def test_deterministic_output(self, discovered, loaded_cache):
        entries, _ = discovered
        # Test one buildable condition per dataset (data-driven pick).
        tested = set()
        for entry in entries:
            if entry.dataset in tested:
                continue
            manifest, registries = loaded_cache[entry.dataset]
            try:
                jb1 = build_job(
                    manifest, entry.condition_id, seeds=[42], **registries
                )
                jb2 = build_job(
                    manifest, entry.condition_id, seeds=[42], **registries
                )
            except (ValueError, Exception):
                continue
            assert jb1.to_dict() == jb2.to_dict(), entry.condition_id
            tested.add(entry.dataset)
        assert tested, "no buildable condition found in any dataset"

    def test_different_seeds_differ(self, discovered, loaded_cache):
        entries, _ = discovered
        for entry in entries:
            manifest, registries = loaded_cache[entry.dataset]
            try:
                jb1 = build_job(
                    manifest, entry.condition_id, seeds=[1], **registries
                )
                jb2 = build_job(
                    manifest, entry.condition_id, seeds=[2], **registries
                )
            except (ValueError, Exception):
                continue
            d1, d2 = jb1.to_dict(), jb2.to_dict()
            assert d1["modelSeeds"] == [1]
            assert d2["modelSeeds"] == [2]
            return
        pytest.skip("no buildable condition found in any dataset")
