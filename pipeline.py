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
from typing import Any, Dict, List, Optional

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

def _stage_generate_figures(
    output_dir: Path,
    *,
    visualization_version: str = "v1",
    metadata_path: Optional[str] = None,
    environment_filter: Optional[str] = None,
    reference_condition: Optional[str] = None,
) -> StageResult:
    """Generate figures using the visualisation orchestrator.

    Parameters
    ----------
    output_dir : Path
        The run directory containing ``tables/`` and ``figures/``.
    visualization_version : str
        ``"v1"`` for the existing 8 figures, ``"v2"`` for the
        publication-quality factor-aware figures, or ``"both"``
        to generate both.
    metadata_path : str, optional
        Path to ``experiment_metadata.json`` for V2.
    environment_filter : str, optional
        Restrict V2 to a single environment.
    reference_condition : str, optional
        Reference condition for structural figures.
    """
    t0 = time.time()
    try:
        from af3_analysis.visualization.orchestrator import (
            generate_all_figures,
            generate_all_figures_v2,
        )
        tables_dir = output_dir / "tables"
        figures_dir = output_dir / "figures"
        results = {}

        # V1 figures
        if visualization_version in ("v1", "both"):
            v1_results = generate_all_figures(
                str(tables_dir), str(figures_dir),
                metadata_path=metadata_path,
            )
            if v1_results:
                results.update(v1_results)

        # V2 figures
        if visualization_version in ("v2", "both"):
            v2_results = generate_all_figures_v2(
                str(tables_dir), str(figures_dir),
                metadata_path=metadata_path,
                environment_filter=environment_filter,
                reference_condition=reference_condition,
            )
            if v2_results:
                results.update(v2_results)

        status = "pass" if results else "fail"
        return StageResult(
            name=f"generate_figures_{visualization_version}",
            status=status,
            duration_s=time.time() - t0,
            message=f"Figure generation completed ({visualization_version})",
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
    config: AnalysisConfig,
    *,
    save_raw_json: bool = False,
    save_summary_json: bool = False,
    visualization_version: str = "both",
    environment_filter: Optional[str] = None,
) -> PipelineResult:
    """Execute the full AF3 analysis pipeline.

    Parameters
    ----------
    config : AnalysisConfig
        Pipeline configuration.
    save_raw_json : bool
        If *True*, save raw extraction JSON.
    save_summary_json : bool
        If *True*, save summary extraction JSON.
    visualization_version : str
        ``"v1"`` for existing figures, ``"v2"`` for publication-quality
        figures, ``"both"`` for both.  Default is ``"both"``.
    environment_filter : str, optional
        Restrict V2 to a single environment.

    Returns
    -------
    PipelineResult
    """
    t_start = time.time()
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
    if s2.status != "pass":
        pipeline.errors.append(s2.message)

    # Stage 3 – manifest
    s3 = _stage_build_manifest(config, {})
    pipeline.stages.append(s3)
    if s3.status != "pass":
        pipeline.errors.append(s3.message)
        return pipeline

    # Stage 4 – analysis (placeholder)
    s4 = _stage_run_analysis(Path(run_dir), config)
    pipeline.stages.append(s4)
    if s4.status != "pass":
        pipeline.errors.append(s4.message)

    # Stage 4b – structural analysis (optional)
    if getattr(config, "coordinate_analysis_enabled", False):
        s4b = _stage_structural_analysis(Path(run_dir), config)
        pipeline.stages.append(s4b)
        if s4b.status != "pass":
            pipeline.errors.append(s4b.message)

    # Stage 5 – figures
    if getattr(config, "generate_figures", True):
        # Resolve metadata path for V2
        meta_path = None
        raw_root = getattr(config, "raw_af3_root", None)
        if raw_root is not None:
            candidate = Path(raw_root) / "experiment_metadata.json"
            if candidate.is_file():
                meta_path = str(candidate)

        s5 = _stage_generate_figures(
            Path(run_dir),
            visualization_version=visualization_version,
            metadata_path=meta_path,
            environment_filter=environment_filter,
            reference_condition=config.reference_condition,
        )
        pipeline.stages.append(s5)
        if s5.status == "pass":
            pipeline.n_manifest_rows += s5.records
        if s5.status != "pass":
            pipeline.errors.append(s5.message)

    # Stage 6 – reports
    s6 = _stage_generate_reports(Path(run_dir), pipeline)
    pipeline.stages.append(s6)
    if s6.status != "pass":
        pipeline.errors.append(s6.message)

    # Determine overall success and elapsed time
    pipeline.success = all(st.status == "pass" for st in pipeline.stages)
    pipeline.elapsed_s = time.time() - t_start
    return pipeline

def _stage_structural_analysis(
    run_dir: Path,
    config: "AnalysisConfig",
) -> StageResult:
    """Run structural/geometric analysis on AF3 predicted structures.

    This stage is optional and only runs when
    ``config.coordinate_analysis_enabled`` is True.

    Uses metadata-driven discovery, validation, and QC reporting.
    """
    import json
    import time
    t0 = time.time()

    try:
        from af3_analysis.io.structure_reader import parse_mmcif, StructureParseError
        from af3_analysis.structural.config import StructuralConfig
        from af3_analysis.structural.discovery import (
            discover_structures,
            populate_composition,
            get_all_records,
        )
        from af3_analysis.structural.comparison import compare_conditions
        from af3_analysis.structural.metrics.registry import get_default_registry
        from af3_analysis.structural.tables import write_structural_tables
    except ImportError as e:
        return StageResult(
            name="structural_analysis",
            status="fail",
            duration_s=time.time() - t0,
            message=f"Structural analysis dependencies missing: {e}",
        )

    raw_root = getattr(config, "raw_af3_root", None)
    if raw_root is None:
        return StageResult(
            name="structural_analysis",
            status="fail",
            duration_s=time.time() - t0,
            message="No raw_af3_root configured for structural analysis",
        )

    raw_root = Path(raw_root)
    if not raw_root.exists():
        return StageResult(
            name="structural_analysis",
            status="fail",
            duration_s=time.time() - t0,
            message=f"Raw AF3 root does not exist: {raw_root}",
        )

    print("[Structural] Stage 4b: Structural analysis")
    print(f"[Structural] Root: {raw_root}")

    # ------------------------------------------------------------------
    # Phase 1: Discover structures
    # ------------------------------------------------------------------
    print("[Structural] Phase 1: Discovering structure files...")

    # Load expected conditions from experiment metadata if available
    expected_conditions = None
    metadata_path = raw_root / "experiment_metadata.json"
    if metadata_path.is_file():
        try:
            with open(metadata_path) as f:
                metadata = json.load(f)
            expected_conditions = list(metadata.get("conditions", {}).keys())
            print(f"[Structural] Loaded metadata: {len(expected_conditions)} expected conditions")
        except Exception as e:
            print(f"[Structural] Warning: could not load metadata: {e}")

    report = discover_structures(raw_root, expected_conditions=expected_conditions)
    print(f"[Structural] Discovered {report.n_discovered} structure files across {report.n_conditions} conditions")
    print(f"[Structural] Orphans: {report.n_orphan}")

    if report.n_discovered == 0:
        return StageResult(
            name="structural_analysis",
            status="pass",
            duration_s=time.time() - t0,
            message="No structure files found; structural analysis skipped",
            records=0,
        )

    # ------------------------------------------------------------------
    # Phase 2: Parse structures
    # ------------------------------------------------------------------
    print("[Structural] Phase 2: Parsing CIF files...")

    structures = []
    parse_errors = []
    all_records = get_all_records(report)
    n_total = len(all_records)

    for i, rec in enumerate(all_records):
        try:
            struct = parse_mmcif(rec.path, rec.condition_id, rec.seed, rec.sample)
            structures.append(struct)
        except StructureParseError as e:
            parse_errors.append({
                "condition_id": rec.condition_id,
                "seed": rec.seed,
                "sample": rec.sample,
                "source_path": str(rec.path),
                "status": "parse_error",
                "reason": str(e),
            })
        except Exception as e:
            parse_errors.append({
                "condition_id": rec.condition_id,
                "seed": rec.seed,
                "sample": rec.sample,
                "source_path": str(rec.path),
                "status": "parse_error",
                "reason": f"Unexpected error: {e}",
            })
        # Progress logging every 50 files
        if (i + 1) % 50 == 0 or (i + 1) == n_total:
            print(f"[Structural]   Parsed {i + 1}/{n_total}")

    report.n_parse_success = len(structures)
    report.n_parse_failure = len(parse_errors)
    report.parse_errors = parse_errors

    print(f"[Structural] Parsed: {len(structures)} success, {len(parse_errors)} failures")

    if not structures:
        return StageResult(
            name="structural_analysis",
            status="fail",
            duration_s=time.time() - t0,
            message=f"All structure files failed to parse ({len(parse_errors)} errors)",
        )

    # ------------------------------------------------------------------
    # Phase 3: Populate composition information
    # ------------------------------------------------------------------
    print("[Structural] Phase 3: Analyzing entity composition...")
    populate_composition(report, structures)

    # ------------------------------------------------------------------
    # Phase 4: Write QC report
    # ------------------------------------------------------------------
    print("[Structural] Phase 4: Writing QC report...")
    qc_path = run_dir / "tables" / "structural_qc_report.txt"
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    with open(qc_path, "w", encoding="utf-8") as f:
        f.write(report.summary())
    print(f"[Structural] QC report: {qc_path}")

    # ------------------------------------------------------------------
    # Phase 5: Run comparisons
    # ------------------------------------------------------------------
    print("[Structural] Phase 5: Running structural comparisons...")

    # Build structural config from analysis config
    ref_condition = getattr(config, "reference_condition", None)
    struct_config = StructuralConfig(
        enabled=True,
        reference_condition=ref_condition,
        reference_strategy="explicit_reference" if ref_condition else "first_condition",
    )

    registry = get_default_registry()
    comparisons = compare_conditions(structures, struct_config, registry=registry)
    print(f"[Structural] Generated {len(comparisons)} structural comparisons")

    # ------------------------------------------------------------------
    # Phase 6: Write output tables
    # ------------------------------------------------------------------
    print("[Structural] Phase 6: Writing output tables...")
    tables_dir = run_dir / "tables"
    try:
        output_paths = write_structural_tables(
            tables_dir, structures, comparisons, registry, struct_config,
            parse_errors=parse_errors,
        )
        for name, path in output_paths.items():
            print(f"[Structural]   {name}: {path}")
    except Exception as e:
        return StageResult(
            name="structural_analysis",
            status="fail",
            duration_s=time.time() - t0,
            message=f"Failed to write structural tables: {e}",
        )

    elapsed = time.time() - t0
    print(f"[Structural] Structural analysis completed in {elapsed:.1f}s")
    print(f"[Structural]   Structures: {len(structures)}")
    print(f"[Structural]   Comparisons: {len(comparisons)}")
    print(f"[Structural]   Conditions: {report.n_conditions}")
    print(f"[Structural]   Parse errors: {report.n_parse_failure}")

    return StageResult(
        name="structural_analysis",
        status="pass",
        duration_s=elapsed,
        message=(
            f"Structural analysis: {len(structures)} structures, "
            f"{len(comparisons)} comparisons, {report.n_parse_failure} parse errors"
        ),
        records=len(structures) + len(comparisons),
    )


__all__ = ["StageResult", "PipelineResult", "run_pipeline"]