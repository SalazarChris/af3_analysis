"""Tests for V3 reference resolution fallbacks.

The af3.py menu supplies the raw AF3 output folder name (condition stem,
e.g. 'pou_baseline') which is not a canonical condition id. resolve_reference
must resolve it transparently with explicit provenance, and fail loudly
when unresolvable or ambiguous."""

import pytest

from af3_analysis.visualization_v3.model import (
    ConditionSummary,
    Dataset,
)
from af3_analysis.visualization_v3.reference import (
    ReferenceResolutionError,
    resolve_reference,
)


def _dataset(
    condition_names: dict,
    stem_candidates: dict = None,
    n_predictions: dict = None,
) -> Dataset:
    conditions = {
        cid: ConditionSummary(
            condition_id=cid,
            condition_name=name,
            n_seeds=1,
            n_predictions=(n_predictions or {}).get(cid, 10),
        )
        for cid, name in condition_names.items()
    }
    return Dataset(
        name="test",
        conditions=conditions,
        condition_stem_candidates=stem_candidates or {},
    )


class TestResolveReference:
    def test_exact_id_still_works(self):
        ds = _dataset({"cond_001": "pou_baseline", "cond_002": "pou_dna"})
        result = resolve_reference(ds, {"condition": "cond_001"})
        assert result["reference_condition"] == "cond_001"
        assert result["reference_strategy"] == "explicit_reference"

    def test_resolves_condition_name(self):
        ds = _dataset({"cond_002": "pou_baseline", "cond_003": "pou_dna"})
        result = resolve_reference(ds, {"condition": "pou_baseline"})
        assert result["reference_condition"] == "cond_002"
        assert result["reference_strategy"] == "explicit_reference_by_name"
        assert result["resolved_from_name"] == "pou_baseline"

    def test_resolves_unique_stem(self):
        # Stem whose condition_name differs from the stem: must resolve via
        # the stem-candidate map.
        ds = _dataset(
            {"cond_005": "OCT4 + K133-Ub"},
            stem_candidates={"pou_sep102": ["cond_005"]},
        )
        result = resolve_reference(ds, {"condition": "pou_sep102"})
        assert result["reference_condition"] == "cond_005"
        assert result["reference_strategy"] == "explicit_reference_by_stem"
        assert result["resolved_from_stem"] == "pou_sep102"

    def test_ambiguous_stem_disambiguated_by_prediction_count(self):
        # Legacy fragment (4 predictions) vs full condition (102): the full
        # condition must win, with the decision recorded explicitly.
        ds = _dataset(
            {"cond_001": "pou_baseline", "cond_002": "pou_baseline"},
            stem_candidates={"pou_baseline": ["cond_001", "cond_002"]},
            n_predictions={"cond_001": 4, "cond_002": 102},
        )
        result = resolve_reference(ds, {"condition": "pou_baseline"})
        assert result["reference_condition"] == "cond_002"
        assert result["reference_strategy"] == (
            "explicit_reference_by_stem_disambiguated"
        )
        assert result["disambiguated_alternatives"] == ["cond_001"]
        assert "cond_001" in result["disambiguation_note"]

    def test_truly_ambiguous_stem_raises(self):
        ds = _dataset(
            {"cond_001": "pou_baseline", "cond_002": "pou_baseline"},
            stem_candidates={"pou_baseline": ["cond_001", "cond_002"]},
            n_predictions={"cond_001": 102, "cond_002": 102},
        )
        with pytest.raises(ReferenceResolutionError, match="ambiguous"):
            resolve_reference(ds, {"condition": "pou_baseline"})

    def test_unknown_reference_raises_with_candidates_listed(self):
        ds = _dataset({"cond_001": "pou_baseline", "cond_002": "pou_dna"})
        with pytest.raises(ReferenceResolutionError) as excinfo:
            resolve_reference(ds, {"condition": "nonexistent"})
        msg = str(excinfo.value)
        assert "nonexistent" in msg
        assert "cond_001" in msg and "pou_baseline" in msg

    def test_none_reference_defaults_to_first_condition(self):
        ds = _dataset({"cond_002": "pou_baseline", "cond_001": "pou_dna"})
        result = resolve_reference(ds, {})
        assert result["reference_condition"] == "cond_001"
        assert result["reference_strategy"] == "first_condition"
