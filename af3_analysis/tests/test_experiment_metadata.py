"""
Tests for experiment design metadata system.

Covers:
1. Current POU dataset loads correctly
2. All 8 POU conditions map correctly
3. Arbitrary attribute names work
4. Multiple components work
5. Different number of attributes works
6. Incomplete combinations detected correctly
7. Missing condition metadata produces explicit error
8. Duplicate definitions produce explicit error
9. Validation against observed data works
"""

import json
import os
import tempfile

import pytest

# Module under test
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from experiment_metadata import (
    AttributeDefinition,
    ConditionDefinition,
    ExperimentDesign,
    DesignInspection,
    load_experiment_design,
    validate_metadata_against_data,
    inspect_design,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

PUO_METADATA_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "testdata", "pou2", "experiment_metadata.json"
)


@pytest.fixture
def pou_design():
    """Load the real POU experiment metadata."""
    return load_experiment_design(PUO_METADATA_PATH)


def _write_metadata(tmp_path, data):
    """Helper to write metadata JSON and return path."""
    p = tmp_path / "experiment_metadata.json"
    p.write_text(json.dumps(data), encoding="utf-8")
    return str(p)


# ---------------------------------------------------------------------------
# 1. Current POU dataset loads correctly
# ---------------------------------------------------------------------------

class TestPOUDataset:
    def test_loads_without_error(self, pou_design):
        assert isinstance(pou_design, ExperimentDesign)

    def test_experiment_id(self, pou_design):
        assert pou_design.experiment_id == "pou_2024"

    def test_eight_conditions(self, pou_design):
        assert pou_design.n_conditions == 8

    def test_three_attributes(self, pou_design):
        assert pou_design.n_attributes == 3

    def test_attribute_names(self, pou_design):
        assert set(pou_design.attribute_names) == {"DNA", "pTPO101", "pSEP102"}

    def test_all_binary(self, pou_design):
        for name in pou_design.attribute_names:
            assert pou_design.attributes[name].attr_type == "binary"


# ---------------------------------------------------------------------------
# 2. All 8 POU conditions map correctly
# ---------------------------------------------------------------------------

class TestPOUConditionMapping:
    def test_baseline_has_no_attributes(self, pou_design):
        cond = pou_design.conditions["pou_baseline"]
        assert cond.attributes["DNA"] is False
        assert cond.attributes["pTPO101"] is False
        assert cond.attributes["pSEP102"] is False

    def test_dna_only(self, pou_design):
        cond = pou_design.conditions["pou_dna"]
        assert cond.attributes["DNA"] is True
        assert cond.attributes["pTPO101"] is False
        assert cond.attributes["pSEP102"] is False

    def test_tpo101_only(self, pou_design):
        cond = pou_design.conditions["pou_tpo101"]
        assert cond.attributes["DNA"] is False
        assert cond.attributes["pTPO101"] is True
        assert cond.attributes["pSEP102"] is False

    def test_sep102_only(self, pou_design):
        cond = pou_design.conditions["pou_sep102"]
        assert cond.attributes["DNA"] is False
        assert cond.attributes["pTPO101"] is False
        assert cond.attributes["pSEP102"] is True

    def test_all_present(self, pou_design):
        cond = pou_design.conditions["pou_tpo101_sep102_dna"]
        assert cond.attributes["DNA"] is True
        assert cond.attributes["pTPO101"] is True
        assert cond.attributes["pSEP102"] is True

    def test_labels_present(self, pou_design):
        for name in pou_design.condition_names:
            cond = pou_design.conditions[name]
            assert cond.label, f"Condition {name} has empty label"

    def test_display_label_lookup(self, pou_design):
        label = pou_design.get_condition_label("pou_baseline")
        assert label == "Baseline (POU)"

    def test_unknown_condition_gets_generic_label(self, pou_design):
        label = pou_design.get_condition_label("unknown_condition")
        assert label == "Unknown Condition"


# ---------------------------------------------------------------------------
# 3. Arbitrary attribute names work
# ---------------------------------------------------------------------------

class TestArbitraryAttributes:
    def test_custom_attribute_names(self, tmp_path):
        data = {
            "experiment_id": "custom_test",
            "description": "Test with arbitrary names",
            "attributes": {
                "ligand_A": {"type": "binary", "description": "Ligand A"},
                "temperature": {"type": "categorical", "description": "Temperature"},
                "ion_Mg": {"type": "binary", "description": "Magnesium ion"},
            },
            "conditions": {
                "cond_control": {
                    "label": "Control",
                    "attributes": {"ligand_A": False, "temperature": "25C", "ion_Mg": False},
                },
                "cond_ligand": {
                    "label": "With Ligand A",
                    "attributes": {"ligand_A": True, "temperature": "25C", "ion_Mg": False},
                },
            },
        }
        path = _write_metadata(tmp_path, data)
        design = load_experiment_design(path)

        assert design.n_attributes == 3
        assert set(design.attribute_names) == {"ion_Mg", "ligand_A", "temperature"}

        # Categorical attribute
        assert design.attributes["temperature"].attr_type == "categorical"
        assert design.get_attribute_values("temperature") == {"25C"}

        # Binary attributes
        assert design.get_attribute_values("ligand_A") == {True, False}


# ---------------------------------------------------------------------------
# 4. Multiple components work
# ---------------------------------------------------------------------------

class TestMultipleComponents:
    def test_five_attributes(self, tmp_path):
        data = {
            "experiment_id": "multi_test",
            "description": "Five attributes",
            "attributes": {
                "A": {"type": "binary"},
                "B": {"type": "binary"},
                "C": {"type": "binary"},
                "D": {"type": "binary"},
                "E": {"type": "binary"},
            },
            "conditions": {
                "c1": {"label": "C1", "attributes": {"A": False, "B": False, "C": False, "D": False, "E": False}},
                "c2": {"label": "C2", "attributes": {"A": True, "B": True, "C": True, "D": True, "E": True}},
            },
        }
        path = _write_metadata(tmp_path, data)
        design = load_experiment_design(path)
        assert design.n_attributes == 5
        assert design.n_conditions == 2


# ---------------------------------------------------------------------------
# 5. Different number of attributes works
# ---------------------------------------------------------------------------

class TestDifferentAttributeCounts:
    def test_single_attribute(self, tmp_path):
        data = {
            "experiment_id": "single",
            "description": "One attribute",
            "attributes": {"treatment": {"type": "binary"}},
            "conditions": {
                "control": {"label": "Control", "attributes": {"treatment": False}},
                "treated": {"label": "Treated", "attributes": {"treatment": True}},
            },
        }
        path = _write_metadata(tmp_path, data)
        design = load_experiment_design(path)
        assert design.n_attributes == 1
        assert design.n_conditions == 2

    def test_six_attributes(self, tmp_path):
        attrs = {f"attr_{i}": {"type": "binary"} for i in range(6)}
        conds = {
            "c_low": {"label": "Low", "attributes": {k: False for k in attrs}},
            "c_high": {"label": "High", "attributes": {k: True for k in attrs}},
        }
        data = {
            "experiment_id": "six_attrs",
            "description": "Six attributes",
            "attributes": attrs,
            "conditions": conds,
        }
        path = _write_metadata(tmp_path, data)
        design = load_experiment_design(path)
        assert design.n_attributes == 6


# ---------------------------------------------------------------------------
# 6. Incomplete combinations detected correctly
# ---------------------------------------------------------------------------

class TestIncompleteCombinations:
    def test_incomplete_factorial_detected(self, tmp_path):
        data = {
            "experiment_id": "incomplete",
            "description": "Missing combinations",
            "attributes": {
                "A": {"type": "binary"},
                "B": {"type": "binary"},
            },
            "conditions": {
                "c00": {"label": "C00", "attributes": {"A": False, "B": False}},
                "c11": {"label": "C11", "attributes": {"A": True, "B": True}},
                # Missing: (A=False, B=True) and (A=True, B=False)
            },
        }
        path = _write_metadata(tmp_path, data)
        design = load_experiment_design(path)
        insp = inspect_design(design)

        assert insp.n_conditions == 2
        assert insp.expected_n_conditions == 4
        assert not insp.is_complete_factorial
        assert len(insp.missing_combinations) == 2

    def test_complete_factorial_detected(self, pou_design):
        insp = inspect_design(pou_design)
        assert insp.is_complete_factorial
        assert insp.expected_n_conditions == 8
        assert len(insp.missing_combinations) == 0


# ---------------------------------------------------------------------------
# 7. Missing condition metadata produces explicit error
# ---------------------------------------------------------------------------

class TestMissingMetadata:
    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_experiment_design("/nonexistent/path.json")

    def test_empty_conditions_raises(self, tmp_path):
        data = {
            "experiment_id": "empty",
            "attributes": {"A": {"type": "binary"}},
            "conditions": {},
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="no 'conditions' defined"):
            load_experiment_design(path)

    def test_missing_experiment_id_raises(self, tmp_path):
        data = {
            "description": "No ID",
            "attributes": {"A": {"type": "binary"}},
            "conditions": {"c1": {"label": "C1", "attributes": {"A": True}}},
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="missing required field 'experiment_id'"):
            load_experiment_design(path)

    def test_missing_attributes_in_condition_raises(self, tmp_path):
        data = {
            "experiment_id": "test",
            "attributes": {"A": {"type": "binary"}, "B": {"type": "binary"}},
            "conditions": {
                "c1": {"label": "C1", "attributes": {"A": True}},  # missing B
            },
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="missing attributes"):
            load_experiment_design(path)

    def test_undefined_attribute_in_condition_raises(self, tmp_path):
        data = {
            "experiment_id": "test",
            "attributes": {"A": {"type": "binary"}},
            "conditions": {
                "c1": {"label": "C1", "attributes": {"A": True, "Z": False}},
            },
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="undefined attribute"):
            load_experiment_design(path)


# ---------------------------------------------------------------------------
# 8. Duplicate definitions produce explicit error
# ---------------------------------------------------------------------------

class TestDuplicateDefinitions:
    def test_duplicate_labels_raise(self, tmp_path):
        data = {
            "experiment_id": "test",
            "attributes": {"A": {"type": "binary"}},
            "conditions": {
                "c1": {"label": "Same Label", "attributes": {"A": True}},
                "c2": {"label": "Same Label", "attributes": {"A": False}},
            },
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="Duplicate label"):
            load_experiment_design(path)

    def test_invalid_attribute_type_raises(self, tmp_path):
        data = {
            "experiment_id": "test",
            "attributes": {"A": {"type": "invalid"}},
            "conditions": {
                "c1": {"label": "C1", "attributes": {"A": True}},
            },
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="invalid type"):
            load_experiment_design(path)

    def test_wrong_type_for_binary_attribute(self, tmp_path):
        data = {
            "experiment_id": "test",
            "attributes": {"A": {"type": "binary"}},
            "conditions": {
                "c1": {"label": "C1", "attributes": {"A": "yes"}},  # should be bool
            },
        }
        path = _write_metadata(tmp_path, data)
        with pytest.raises(ValueError, match="expected boolean"):
            load_experiment_design(path)


# ---------------------------------------------------------------------------
# 9. Validation against observed data
# ---------------------------------------------------------------------------

class TestValidationAgainstData:
    def test_matching_data_passes(self, pou_design):
        observed = [f"pou_{s}" for s in [
            "baseline", "dna", "sep102", "sep102_dna",
            "tpo101", "tpo101_dna", "tpo101_sep102", "tpo101_sep102_dna",
        ]]
        # Add missing prefix to match actual names
        observed = [
            "pou_baseline", "pou_dna", "pou_sep102", "pou_sep102_dna",
            "pou_tpo101", "pou_tpo101_dna", "pou_tpo101_sep102", "pou_tpo101_sep102_dna",
        ]
        messages = validate_metadata_against_data(pou_design, observed)
        # No errors expected
        assert not any("ERROR" in m for m in messages)

    def test_extra_condition_in_data_warns(self, pou_design):
        observed = [
            "pou_baseline", "pou_dna", "pou_unknown_extra",
        ]
        messages = validate_metadata_against_data(pou_design, observed)
        assert any("pou_unknown_extra" in m for m in messages)

    def test_missing_condition_in_data_warns(self, pou_design):
        observed = ["pou_baseline"]  # only 1 of 8
        messages = validate_metadata_against_data(pou_design, observed)
        assert any("not in data" in m for m in messages)


# ---------------------------------------------------------------------------
# Design inspection summary
# ---------------------------------------------------------------------------

class TestDesignInspection:
    def test_summary_contains_experiment_id(self, pou_design):
        insp = inspect_design(pou_design)
        summary = insp.summary()
        assert "pou_2024" in summary

    def test_summary_contains_condition_count(self, pou_design):
        insp = inspect_design(pou_design)
        summary = insp.summary()
        assert "8" in summary  # 8 conditions

    def test_summary_reports_complete_factorial(self, pou_design):
        insp = inspect_design(pou_design)
        summary = insp.summary()
        assert "Complete factorial" in summary


# ---------------------------------------------------------------------------
# Condition attributes matrix
# ---------------------------------------------------------------------------

class TestConditionMatrix:
    def test_matrix_shape(self, pou_design):
        matrix = pou_design.condition_attributes_matrix()
        assert len(matrix) == 8
        for cond_name in matrix:
            assert len(matrix[cond_name]) == 3

    def test_matrix_values(self, pou_design):
        matrix = pou_design.condition_attributes_matrix()
        assert matrix["pou_baseline"]["DNA"] is False
        assert matrix["pou_dna"]["DNA"] is True


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_get_conditions_with_attribute(self, pou_design):
        dna_conditions = pou_design.get_conditions_with_attribute("DNA", True)
        assert len(dna_conditions) == 4
        assert "pou_dna" in dna_conditions
        assert "pou_baseline" not in dna_conditions

    def test_get_attribute_values(self, pou_design):
        values = pou_design.get_attribute_values("DNA")
        assert values == {True, False}

    def test_unknown_attribute_raises(self, pou_design):
        with pytest.raises(KeyError, match="not defined"):
            pou_design.get_attribute_values("nonexistent")

    def test_categorical_attributes(self, tmp_path):
        data = {
            "experiment_id": "cat_test",
            "description": "Categorical attributes",
            "attributes": {
                "buffer": {"type": "categorical", "description": "Buffer type"},
            },
            "conditions": {
                "pbs": {"label": "PBS", "attributes": {"buffer": "PBS"}},
                "tris": {"label": "Tris", "attributes": {"buffer": "Tris"}},
                "hepes": {"label": "HEPES", "attributes": {"buffer": "HEPES"}},
            },
        }
        path = _write_metadata(tmp_path, data)
        design = load_experiment_design(path)

        assert design.get_attribute_values("buffer") == {"PBS", "Tris", "HEPES"}
        assert design.get_attribute_type("buffer") == "categorical"
