"""Tests for the V3 pipeline stage wiring (pipeline.py)."""

import inspect
from pathlib import Path

from af3_analysis import pipeline as af3_pipeline


def test_run_pipeline_accepts_v3_enabled_flag():
    """run_pipeline must expose v3_enabled (default True) so the menu can
    toggle the V3 stage."""
    sig = inspect.signature(af3_pipeline.run_pipeline)
    assert "v3_enabled" in sig.parameters
    assert sig.parameters["v3_enabled"].default is True


def test_v3_stage_function_exists():
    assert hasattr(af3_pipeline, "_stage_v3_visualization")


def test_v3_stage_fails_cleanly_without_run_dir(tmp_path):
    """An empty run directory (no tables/) must produce a clean 'fail'
    StageResult, not an exception."""
    result = af3_pipeline._stage_v3_visualization(
        tmp_path,
        raw_af3_root=None,
        metadata_path=None,
        reference_condition=None,
    )
    assert result.name == "v3_visualization"
    assert result.status == "fail"
    assert result.message  # non-empty diagnostic


def test_v3_stage_sets_reference_condition(tmp_path):
    """With a reference condition and a missing run dir the stage still
    fails cleanly (validation happens inside the runner), proving the
    config path executes without import/scope errors."""
    result = af3_pipeline._stage_v3_visualization(
        tmp_path,
        raw_af3_root=None,
        metadata_path=None,
        reference_condition="cond_001",
    )
    assert result.name == "v3_visualization"
    assert result.status == "fail"
