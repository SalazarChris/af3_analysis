# -*- coding: utf-8 -*-
"""AF3 Analysis Pipeline Orchestrator
===================================================
Provides a single `run_pipeline()` function that executes the complete
M0–M9 analysis workflow. The implementation stitches together the existing
condition‑centric extraction script, basic validation, manifest creation, the
placeholder statistical analysis stage, figure generation, and a minimal
report placeholder.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

import re

# Local imports from the AF3 analysis package
from .config import AnalysisConfig, create_run_directory


# ---------------------------------------------------------------------------
# Seed extraction (P7 fix)
# ---------------------------------------------------------------------------

_SEED_PATTERN = re.compile(r"seed-(\d+)")


def extract_seed(rep_id: str) -> int:
    """Extract the integer seed from a replicate ID string.

    Expected format: ``<anything>_seed-<N>_sample-<M>``

    Returns
    -------
    int
        The parsed seed value.

    Raises
    ------
    ValueError
        If *rep_id* does not contain a valid ``seed-<N>`` segment.
    """
    match = _SEED_PATTERN.search(str(rep_id))
    if match is None:
        raise ValueError(
            f"Cannot extract seed from replicate_id '{rep_id}': "
            f"no 'seed-<N>' segment found."
        )
    return int(match.group(1))

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class StageResult:
    name: str
    status: str  # "pass" or "fail"
    duration_s: float = 0.0
    message: str = ""
    records: int = 0

@dataclass
class PipelineResult:
    success: bool = False
    stages: List[StageResult] = field(default_factory=list)
    n_conditions: int = 0
    n_replicates: int = 0
    n_manifest_rows: int = 0
    output_dir: str = ""
    run_id: str = ""
    errors: List[str] = field(default_factory=list)
    elapsed_s: float = 0.0

# ---------------------------------------------------------------------------
# Stage 1 – Extraction (uses the existing script)
# ---------------------------------------------------------------------------

def _stage_extract_metrics(
    config: AnalysisConfig, *, raw_json: bool = False, summary_json: bool = False
) -> StageResult:
    """Run the condition‑centric extraction script.

    The script lives under ``af3inputbuilder/scripts/af3_condition_centric_extraction.py``
    and provides an ``extract_metrics`` function.
    """
    t0 = time.time()
    try:
        from importlib.util import spec_from_file_location, module_from_spec
        # Walk up from pipeline.py to find the sibling af3inputbuilder/ directory.
        # The exact depth varies depending on install layout (src/ vs flat).
        _pipeline_dir = Path(__file__).resolve().parent
        _af3inputbuilder_script = (
            "af3inputbuilder", "scripts", "af3_condition_centric_extraction.py"
        )
        script_path = None
        for ancestor in [_pipeline_dir] + list(_pipeline_dir.parents):
            candidate = ancestor.joinpath(*_af3inputbuilder_script)
            if candidate.is_file():
                script_path = candidate
                break
        if script_path is None:
            raise FileNotFoundError(
                "Cannot locate af3_condition_centric_extraction.py "
                f"(searched up from {_pipeline_dir})"
            )
        spec = spec_from_file_location("extraction_module", script_path)
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        extract_metrics = module.extract_metrics
    except Exception as e:
        return StageResult(
            name="extract_metrics",
            status="fail",
            duration_s=time.time() - t0,
            message=f"Failed to load extraction module: {e}",
        )

    try:
        result = extract_metrics(
            input_dir=str(config.raw_af3_root),
            output_dir=str(config.output_paths.tables_dir),
            verbose=False,
        )
    except Exception as e:
        return StageResult(
            name="extract_metrics",
            status="fail",
            duration_s=time.time() - t0,
            message=str(e),
        )

    # Optional JSON outputs – user requested to skip, so we ignore flags.
    return StageResult(
        name="extract_metrics",
        status="pass",
        duration_s=time.time() - t0,
        message="Extraction completed",
        records=result.get("replicates", 0),
    )

# ---------------------------------------------------------------------------
# Stage 2 – Validation (placeholder)
# ---------------------------------------------------------------------------

def _stage_validate(config: AnalysisConfig, previous_result: Dict[str, Any]) -> StageResult:
    return StageResult(name="validate", status="pass", message="Validation placeholder")

# ---------------------------------------------------------------------------
# Stage 3 – Manifest creation (writes a simple JSON manifest)
# ---------------------------------------------------------------------------

def _stage_build_manifest(config: AnalysisConfig, previous_result: Dict[str, Any]) -> StageResult:
    manifest = config.to_manifest()
    manifest_path = config.output_paths.manifest_dir / "manifest.json"
    try:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, sort_keys=True)
    except Exception as e:
        return StageResult(name="build_manifest", status="fail", message=str(e))
    return StageResult(name="build_manifest", status="pass", message=str(manifest_path))

# ---------------------------------------------------------------------------
# Placeholder stages that already exist in the original file
# ---------------------------------------------------------------------------

def _stage_run_analysis(output_dir: Path, config: "AnalysisConfig") -> StageResult:
    import pandas as pd
    import numpy as np
    import time
    t0 = time.time()
    
    tables_dir = output_dir / "tables"
    replicates_csv = tables_dir / "metrics_replicates.csv"
    
    if not replicates_csv.exists():
        return StageResult(name="run_analysis", status="fail", message="metrics_replicates.csv not found")
        
    try:
        df = pd.read_csv(replicates_csv)
        
        # Extract seed from replicate_id (e.g. cond1_seed-1_sample-0 -> 1)
        # P7 fix: raise on unparseable IDs instead of silently falling back to seed=1
        bad_ids = []
        seeds = []
        for rep_id in df['replicate_id']:
            try:
                seeds.append(extract_seed(rep_id))
            except ValueError as exc:
                bad_ids.append((rep_id, str(exc)))
                seeds.append(None)
        df['seed'] = seeds
        
        if bad_ids:
            # Drop rows with unparseable seed IDs so aggregation stays correct
            n_bad = len(bad_ids)
            first_msg = bad_ids[0][1]
            df = df.dropna(subset=['seed']).copy()
            df['seed'] = df['seed'].astype(int)
            print(f"  WARNING: {n_bad} replicate_id(s) had no valid seed pattern and were excluded. "
                  f"First: {first_msg}")
        else:
            df['seed'] = df['seed'].astype(int)
        
        # Numeric columns to aggregate
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if 'seed' in numeric_cols:
            numeric_cols.remove('seed')
            
        # Group by condition and seed, taking the mean of all numeric columns
        groupby_cols = ['condition_id', 'condition_name', 'seed']
        available_groupby = [c for c in groupby_cols if c in df.columns]
        
        seed_aggregated = df.groupby(available_groupby)[numeric_cols].mean().reset_index()
        seed_aggregated.to_csv(tables_dir / "seed_aggregated.csv", index=False)
        
        # Compute condition-level descriptive stats (group by condition_id)
        if 'condition_id' in seed_aggregated.columns:
            desc_cols = [c for c in numeric_cols if c in seed_aggregated.columns]
            descriptive_stats = seed_aggregated.groupby('condition_id')[desc_cols].agg(['mean', 'std']).reset_index()
            # Flatten columns for simple CSV
            descriptive_stats.columns = ['_'.join(col).strip('_') for col in descriptive_stats.columns.values]
            descriptive_stats.to_csv(tables_dir / "descriptive_stats.csv", index=False)
        else:
            pd.DataFrame().to_csv(tables_dir / "descriptive_stats.csv", index=False)
            
        # Dummy pairwise comparisons for now since real implementation requires specifying conditions
        pd.DataFrame({'reference': [], 'condition': [], 'metric': [], 'diff_mean': []}).to_csv(tables_dir / "pairwise_comparisons.csv", index=False)

        return StageResult(
            name="run_analysis", 
            status="pass", 
            duration_s=time.time() - t0,
            message="Seed-level aggregation and stats completed",
            records=len(seed_aggregated)
        )
    except Exception as e:
        return StageResult(name="run_analysis", status="fail", duration_s=time.time() - t0, message=str(e))

def _stage_generate_figures(output_dir: Path) -> StageResult:
    """Generate figures using the visualisation orchestrator.

    The orchestrator is imported lazily to avoid import‑time side effects.
    """
    t0 = time.time()
    try:
        from af3_analysis.visualization.orchestrator import generate_all_figures
        tables_dir = output_dir / "tables"
        figures_dir = output_dir / "figures"
        results = generate_all_figures(str(tables_dir), str(figures_dir))
        status = "pass" if results else "fail"
        return StageResult(
            name="generate_figures",
            status=status,
            duration_s=time.time() - t0,
            message="Figure generation completed",
            records=len(results) if results else 0,
        )
    except Exception as e:
        return StageResult(
            name="generate_figures",
            status="fail",
            duration_s=time.time() - t0,
            message=str(e),
        )

def _stage_generate_reports(output_dir: Path, result: PipelineResult) -> StageResult:
    # Minimal placeholder implementation.
    return StageResult(
        name="generate_reports",
        status="pass",
        message="Report generation placeholder",
    )

# ---------------------------------------------------------------------------
# Public entry‑point
# ---------------------------------------------------------------------------

def run_pipeline(
    config: AnalysisConfig, *, save_raw_json: bool = False, save_summary_json: bool = False
) -> PipelineResult:
    """Execute the full AF3 analysis pipeline.

    Returns a :class:`PipelineResult` containing stage information and overall
    success status.
    """
    run_dir = create_run_directory(config)
    pipeline = PipelineResult(output_dir=str(run_dir), run_id=config.run_id)

    # Stage 1 – extraction
    s1 = _stage_extract_metrics(config, raw_json=save_raw_json, summary_json=save_summary_json)
    pipeline.stages.append(s1)
    if s1.status != "pass":
        pipeline.errors.append(s1.message)
        return pipeline

    # Stage 2 – validation
    s2 = _stage_validate(config, {})
    pipeline.stages.append(s2)

    # Stage 3 – manifest
    s3 = _stage_build_manifest(config, {})
    pipeline.stages.append(s3)
    if s3.status != "pass":
        pipeline.errors.append(s3.message)
        return pipeline

    # Stage 4 – analysis (placeholder)
    s4 = _stage_run_analysis(Path(run_dir), config)
    pipeline.stages.append(s4)

    # Stage 5 – figures
    if getattr(config, "generate_figures", True):
        s5 = _stage_generate_figures(Path(run_dir))
        pipeline.stages.append(s5)
        if s5.status == "pass":
            pipeline.n_manifest_rows += s5.records

    # Stage 6 – reports
    s6 = _stage_generate_reports(Path(run_dir), pipeline)
    pipeline.stages.append(s6)

    # Determine overall success
    pipeline.success = all(st.status == "pass" for st in pipeline.stages)
    return pipeline

__all__ = ["StageResult", "PipelineResult", "run_pipeline"]