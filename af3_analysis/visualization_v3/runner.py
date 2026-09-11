"""
V3 Runner Module.

Orchestrates the V3 visualization pipeline with:
- [V3] stage logging
- Failure-safe figure generation
- Pathological figure protection
- Caching

Unit-of-analysis notes (see ARCHITECTURE.md):
- prediction level: one (condition, seed, sample) CIF — pairwise matrix,
  contact maps, displacement vectors
- seed level: matched-seed comparisons, effect sizes, reproducibility
- condition level: summaries, factorial contrasts

Reference pairing for per-prediction figures (F02, F15, F20) uses
"matched_sample" mode by default: a target prediction is compared to the
reference prediction with identical (seed, sample) indices. This is
deterministic but arbitrary among samples; seed-level conclusions must use
the matched-seed analyses (F03, F17, F18), which aggregate over samples.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
import traceback
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .config import V3Config, get_enabled_figures
from .adapter import adapt_from_pipeline
from .validation import (
    validate_condition_mapping,
    validate_seed_mapping,
    validate_structural_inputs,
    summarize_structural_qc,
)
from .reference import resolve_reference, validate_reference_resolution
from .analysis.matched_seed import (
    analyze_matched_seeds,
    aggregate_matched_seed_results,
    matched_seed_summary_table,
)
from .structural.rmsd import calculate_rmsd
from .structural.displacement import calculate_per_residue_displacement
from .structural.contacts import calculate_contact_map, calculate_contact_difference
from .structural.interfaces import find_all_interfaces, calculate_interface_comparison
from .structural.pairwise import (
    calculate_pairwise_rmsd_matrix,
    save_pairwise_matrix,
    load_pairwise_matrix,
)
from .structural.clustering import hierarchical_clustering, mds_embedding
from .structural.regions import calculate_local_geometry
from .analysis.confidence import integrate_confidence_geometry
from .analysis.effects import calculate_structural_effect_sizes
from .analysis.factorial import factorial_structural_analysis, factorial_contrasts_table
from .analysis.reproducibility import calculate_seed_reproducibility

from .figures import (
    generate_f01_structural_qc,
    generate_f02_global_structural_difference,
    generate_f03_matched_seed_structural_difference,
    generate_f04_per_residue_displacement,
    generate_f05_displacement_heatmap,
    generate_f06_contact_map_difference,
    generate_f07_contact_change_summary,
    generate_f08_interface_analysis,
    generate_f09_interface_change_map,
    generate_f10_local_geometry,
    generate_f11_domain_motion,
    generate_f12_structural_clustering,
    generate_f13_similarity_matrix,
    generate_f14_mds_embedding,
    generate_f15_confidence_geometry,
    generate_f16_confidence_change_vs_structural_change,
    generate_f17_seed_reproducibility,
    generate_f18_effect_sizes,
    generate_f19_factorial_effects,
    generate_f20_structure_confidence_matrix,
)

logger = logging.getLogger(__name__)


class V3PipelineError(Exception):
    """V3 pipeline error."""
    pass


def run_v3_pipeline(
    run_dir: Path,
    raw_af3_root: Optional[Path] = None,
    experiment_metadata_path: Optional[Path] = None,
    v3_config: Optional[V3Config] = None,
    *,
    overwrite: bool = False,
) -> Dict[str, Any]:
    """
    Run the complete V3 visualization pipeline.

    Parameters
    ----------
    run_dir : Path
        Existing pipeline run directory.
    raw_af3_root : Path, optional
        Raw AF3 output root with CIF files.
    experiment_metadata_path : Path, optional
        Path to experiment_metadata.json.
    v3_config : V3Config, optional
        V3 configuration. Uses defaults if None.
    overwrite : bool
        Whether to overwrite existing V3 output.

    Returns
    -------
    dict with pipeline results
    """
    logger.info("[V3] Starting V3 visualization pipeline")

    # Initialize config
    if v3_config is None:
        v3_config = V3Config()

    if not v3_config.enabled:
        logger.info("[V3] V3 pipeline disabled")
        return {"status": "disabled", "message": "V3 pipeline disabled in config"}

    # Create V3 output directory
    v3_dir = run_dir / "v3"
    v3_figures_dir = v3_dir / "figures"
    v3_tables_dir = v3_dir / "tables"
    v3_metadata_dir = v3_dir / "metadata"
    v3_logs_dir = v3_dir / "logs"
    v3_report_dir = v3_dir / "report"
    v3_cache_dir = v3_dir / "cache"

    for d in [v3_figures_dir, v3_tables_dir, v3_metadata_dir, v3_logs_dir, v3_report_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Initialize result tracking
    results: Dict[str, Any] = {
        "status": "running",
        "start_time": time.time(),
        "config": _config_to_dict(v3_config),
        "figures": {},
        "tables": {},
        "warnings": [],
        "errors": [],
        "skipped": [],
    }

    # Phase 1: Adapt dataset from pipeline
    logger.info("[V3] Phase 1: Adapting dataset from pipeline")
    try:
        dataset = adapt_from_pipeline(
            run_dir,
            raw_af3_root=raw_af3_root,
            experiment_metadata_path=experiment_metadata_path,
            v3_config=v3_config,
        )
        results["dataset_summary"] = {
            "n_conditions": len(dataset.conditions),
            "n_predictions": sum(
                sum(len(sp.values()) for sp in cond_preds.values())
                for cond_preds in dataset.predictions.values()
            ),
        }
    except Exception as e:
        logger.error("[V3] Failed to adapt dataset: %s", e)
        results["status"] = "failed"
        results["errors"].append(f"Dataset adaptation failed: {e}")
        return results

    # Phase 2: Validate dataset
    logger.info("[V3] Phase 2: Validating dataset")
    try:
        validation = validate_condition_mapping(dataset)
        if not validation.passed:
            results["warnings"].extend(validation.messages)

        validation = validate_seed_mapping(dataset)
        if not validation.passed:
            results["warnings"].extend(validation.messages)

        # Structural validation (if structures available)
        if raw_af3_root and raw_af3_root.exists() and dataset.predictions:
            validation = validate_structural_inputs(dataset)
            if not validation.passed:
                results["warnings"].extend(validation.messages)

            qc_summary = summarize_structural_qc(
                validation.details.get("qc_records", [])
            )
            results["structural_qc"] = qc_summary

    except Exception as e:
        logger.error("[V3] Validation error: %s", e)
        results["warnings"].append(f"Validation issue: {e}")

    # Phase 3: Resolve reference
    logger.info("[V3] Phase 3: Resolving reference")
    try:
        ref_resolution = resolve_reference(dataset, v3_config.reference)
        validation = validate_reference_resolution(dataset, ref_resolution)
        if not validation["valid"]:
            results["warnings"].extend(validation["warnings"])
            results["errors"].extend(validation["errors"])

        results["reference_resolution"] = ref_resolution
    except Exception as e:
        logger.error("[V3] Reference resolution error: %s", e)
        results["errors"].append(f"Reference resolution failed: {e}")
        return results

    # Phase 4: Run structural calculations
    logger.info("[V3] Phase 4: Running structural calculations")

    # Get all structures
    all_structures: List[Any] = []
    for condition_id, seeds_dict in dataset.predictions.items():
        for seed, samples_dict in seeds_dict.items():
            for sample, structure in samples_dict.items():
                all_structures.append(structure)

    if not all_structures:
        logger.warning("[V3] No structures available for analysis")
        results["warnings"].append("No structures available")

    # Calculate pairwise RMSD matrix (cached on disk)
    pairwise_matrix = None
    if len(all_structures) > 1:
        cache_key = _pairwise_cache_key(all_structures, v3_config)
        cache_path = v3_cache_dir / f"{cache_key}.csv"
        meta_path = v3_cache_dir / f"{cache_key}_metadata.csv"
        v3_cache_dir.mkdir(parents=True, exist_ok=True)

        if cache_path.exists() and meta_path.exists():
            logger.info("[V3] Loading cached pairwise RMSD matrix (%s)", cache_key)
            try:
                pairwise_matrix = load_pairwise_matrix(cache_path, meta_path)
            except Exception as e:
                logger.warning("[V3] Cache load failed, recalculating: %s", e)
                pairwise_matrix = None

        if pairwise_matrix is None:
            logger.info("[V3] Calculating pairwise RMSD matrix (%d structures)",
                        len(all_structures))
            try:
                pairwise_matrix = calculate_pairwise_rmsd_matrix(
                    all_structures,
                    alignment_atom=v3_config.structure.alignment_atom,
                    min_common_atoms=v3_config.structure.min_common_atoms,
                    min_sequence_identity=v3_config.structure.min_sequence_identity,
                    min_coverage=v3_config.structure.minimum_coverage,
                )
                save_pairwise_matrix(pairwise_matrix, v3_cache_dir,
                                     filename=f"{cache_key}.csv")
            except Exception as e:
                logger.error("[V3] Pairwise matrix calculation failed: %s", e)
                results["errors"].append(f"Pairwise matrix failed: {e}")

        if pairwise_matrix is not None and pairwise_matrix.n_structures > 0:
            results["pairwise_matrix_stats"] = {
                "n_structures": pairwise_matrix.n_structures,
                "mean_distance": _safe_float(pairwise_matrix.mean_distance),
                "median_distance": _safe_float(pairwise_matrix.median_distance),
                "min_distance": _safe_float(pairwise_matrix.min_distance),
                "max_distance": _safe_float(pairwise_matrix.max_distance),
            }

    # Phase 5: Generate figures
    logger.info("[V3] Phase 5: Generating figures")

    enabled_figures = get_enabled_figures(v3_config)
    figure_generators = {
        "F01": generate_f01_structural_qc,
        "F02": generate_f02_global_structural_difference,
        "F03": generate_f03_matched_seed_structural_difference,
        "F04": generate_f04_per_residue_displacement,
        "F05": generate_f05_displacement_heatmap,
        "F06": generate_f06_contact_map_difference,
        "F07": generate_f07_contact_change_summary,
        "F08": generate_f08_interface_analysis,
        "F09": generate_f09_interface_change_map,
        "F10": generate_f10_local_geometry,
        "F11": generate_f11_domain_motion,
        "F12": generate_f12_structural_clustering,
        "F13": generate_f13_similarity_matrix,
        "F14": generate_f14_mds_embedding,
        "F15": generate_f15_confidence_geometry,
        "F16": generate_f16_confidence_change_vs_structural_change,
        "F17": generate_f17_seed_reproducibility,
        "F18": generate_f18_effect_sizes,
        "F19": generate_f19_factorial_effects,
        "F20": generate_f20_structure_confidence_matrix,
    }

    for fig_id in enabled_figures:
        if fig_id not in figure_generators:
            logger.warning("[V3] Unknown figure ID: %s", fig_id)
            results["skipped"].append({"figure_id": fig_id, "reason": "Unknown figure"})
            continue

        generator = figure_generators[fig_id]
        logger.info("[V3] Starting %s", fig_id)

        try:
            fig_result = _generate_figure_with_data(
                generator,
                fig_id,
                dataset,
                all_structures,
                pairwise_matrix,
                ref_resolution,
                v3_config,
                v3_figures_dir,
            )

            results["figures"][fig_id] = fig_result

            if fig_result.get("status") == "pass":
                logger.info("[V3] Completed %s (n=%d)", fig_id, fig_result.get("n_observations", 0))
            elif fig_result.get("status") == "skip":
                logger.info("[V3] Skipped %s: %s", fig_id, fig_result.get("reason", "Unknown"))
                results["skipped"].append({
                    "figure_id": fig_id,
                    "reason": fig_result.get("reason", "Unknown"),
                })
            else:
                logger.warning("[V3] Failed %s: %s", fig_id, fig_result.get("error", "Unknown"))
                results["errors"].append(f"{fig_id}: {fig_result.get('error', 'Unknown')}")

        except Exception as e:
            error_msg = f"{fig_id} failed: {str(e)[:200]}"
            logger.error("[V3] %s", error_msg, exc_info=True)
            results["figures"][fig_id] = {
                "status": "failed",
                "error": str(e),
                "traceback": traceback.format_exc()[:500],
                "output_path": None,
                "n_observations": 0,
                "warnings": [],
            }
            results["errors"].append(error_msg)
            # Continue with other figures
            continue

    # Phase 6: Generate tables
    logger.info("[V3] Phase 6: Generating tables")
    try:
        _generate_v3_tables(
            dataset,
            all_structures,
            pairwise_matrix,
            ref_resolution,
            v3_config,
            v3_tables_dir,
            results,
        )
    except Exception as e:
        logger.error("[V3] Table generation failed: %s", e)
        results["errors"].append(f"Table generation failed: {e}")

    # Finalize results summary and status
    elapsed = time.time() - results["start_time"]
    results["elapsed_s"] = elapsed

    n_success = sum(1 for f in results["figures"].values() if f.get("status") == "pass")
    n_skipped = sum(1 for f in results["figures"].values() if f.get("status") == "skip")
    n_failed = sum(1 for f in results["figures"].values() if f.get("status") == "failed")

    results["summary"] = {
        "n_figures_requested": len(enabled_figures),
        "n_figures_success": n_success,
        "n_figures_skipped": n_skipped,
        "n_figures_failed": n_failed,
        "elapsed_s": elapsed,
    }

    if n_failed == 0:
        results["status"] = "complete"
    else:
        results["status"] = "completed_with_errors"

    # Phase 7: Write manifest
    logger.info("[V3] Phase 7: Writing manifest")
    try:
        _write_v3_manifest(v3_metadata_dir, results)
    except Exception as e:
        logger.error("[V3] Manifest writing failed: %s", e)
        results["errors"].append(f"Manifest writing failed: {e}")

    # Phase 8: Generate report
    logger.info("[V3] Phase 8: Generating report")
    try:
        _generate_v3_report(v3_report_dir, results)
    except Exception as e:
        logger.error("[V3] Report generation failed: %s", e)
        results["errors"].append(f"Report generation failed: {e}")

    logger.info("[V3] Pipeline complete: %d success, %d skipped, %d failed (%.1fs)",
                n_success, n_skipped, n_failed, elapsed)

    return results


# ---------------------------------------------------------------------------
# Figure data preparation
# ---------------------------------------------------------------------------

def _iter_predictions(dataset: Any):
    """Iterate all (condition_id, seed, sample, structure) tuples in stable order."""
    for condition_id in sorted(dataset.predictions.keys()):
        for seed in sorted(dataset.predictions[condition_id].keys()):
            for sample in sorted(dataset.predictions[condition_id][seed].keys()):
                yield (
                    condition_id,
                    seed,
                    sample,
                    dataset.predictions[condition_id][seed][sample],
                )


def _get_reference_structure(
    dataset: Any,
    seed: int,
    sample: int,
    ref_condition: str,
) -> Optional[Any]:
    """Get reference structure with identical (seed, sample) indices (matched_sample)."""
    return dataset.predictions.get(ref_condition, {}).get(seed, {}).get(sample)


def _get_matched_seed_results(
    dataset: Any,
    ref_resolution: Dict[str, Any],
    v3_config: V3Config,
) -> Dict[str, List[Any]]:
    """
    Run matched-seed analysis for every target condition vs the reference.

    Returns condition_id -> list of MatchedSeedResult (one per common seed).
    """
    ref_condition = ref_resolution.get("reference_condition")
    matched: Dict[str, List[Any]] = {}
    if ref_condition is None:
        return matched

    for condition_id in sorted(dataset.predictions.keys()):
        if condition_id == ref_condition:
            continue
        try:
            matched[condition_id] = analyze_matched_seeds(
                dataset,
                ref_condition,
                condition_id,
                reference_condition=ref_condition,
                alignment_atom=v3_config.structure.alignment_atom,
                min_common_atoms=v3_config.structure.min_common_atoms,
                min_coverage=v3_config.structure.minimum_coverage,
                min_sequence_identity=v3_config.structure.min_sequence_identity,
            )
        except Exception as e:
            logger.warning("[V3] Matched-seed analysis failed for %s: %s",
                           condition_id, e)
            matched[condition_id] = []
    return matched


def _generate_figure_with_data(
    generator: Any,
    fig_id: str,
    dataset: Any,
    all_structures: List[Any],
    pairwise_matrix: Any,
    ref_resolution: Dict[str, Any],
    v3_config: V3Config,
    figures_dir: Path,
) -> Dict[str, Any]:
    """Generate a figure with the appropriate data."""

    # Build common parameters
    common_params = {
        "save_path": figures_dir,
        "design": dataset.experiment_metadata,
        "reference_condition": ref_resolution.get("reference_condition"),
    }

    ref_condition = ref_resolution.get("reference_condition")

    # ------------------------------------------------------------------
    if fig_id == "F01":
        # Structural QC — reuse validation-level QC when available.
        # F01 has no reference_condition parameter, so pass only the params
        # it accepts.
        qc_records = []
        for condition_id, seed, sample, structure in _iter_predictions(dataset):
            qc_records.append({
                "prediction_id": structure.prediction_id,
                "condition_id": condition_id,
                "seed": seed,
                "sample": sample,
                "passed": structure.parse_status == "success",
                "reason": structure.parse_reason,
                "entity_types": [e.entity_type for e in structure.entities],
            })
        f01_params = {k: v for k, v in common_params.items() if k != "reference_condition"}
        return generator(qc_records=qc_records, **f01_params)

    # ------------------------------------------------------------------
    if fig_id == "F02":
        # Global structural difference (prediction-level, matched_sample pairing)
        rmsd_results = []
        if ref_condition and dataset.predictions:
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
                if ref_struct is None:
                    continue
                result = calculate_rmsd(
                    target_struct,
                    ref_struct,
                    alignment_atom=v3_config.structure.alignment_atom,
                    min_common_atoms=v3_config.structure.min_common_atoms,
                    min_sequence_identity=v3_config.structure.min_sequence_identity,
                    min_coverage=v3_config.structure.minimum_coverage,
                )
                if result["rmsd"] is not None:
                    rmsd_results.append({
                        "condition_a": condition_id,
                        "condition_b": ref_condition,
                        "seed": seed,
                        "sample": sample,
                        "rmsd": result["rmsd"],
                        "coverage": result["coverage"],
                    })
        return generator(rmsd_results=rmsd_results, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F03":
        # Matched-seed structural difference (seed-level unit of analysis)
        seed_results = []
        matched = _get_matched_seed_results(dataset, ref_resolution, v3_config)
        for condition_id in sorted(matched.keys()):
            for msr in matched[condition_id]:
                if msr.rmsd_mean is None:
                    continue
                seed_results.append({
                    "condition_a": condition_id,
                    "condition_b": ref_condition,
                    "seed": msr.seed,
                    "n_seeds": 1,
                    "rmsd_mean": msr.rmsd_mean,
                    "rmsd_median": msr.rmsd_median,
                    "rmsd_std": msr.rmsd_std,
                    "coverage_mean": msr.coverage_mean,
                })
        return generator(seed_rmsd_results=seed_results, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F04":
        # Per-residue displacement (prediction level, matched_sample pairing)
        displacement_data = []
        if ref_condition and dataset.predictions:
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
                if ref_struct is None:
                    continue
                disp = calculate_per_residue_displacement(
                    ref_struct,
                    target_struct,
                    alignment_atom=v3_config.structure.alignment_atom,
                )
                for row in disp["displacements"]:
                    displacement_data.append({
                        "condition_id": condition_id,
                        "condition_label": condition_id,
                        "seed": seed,
                        "sample": sample,
                        "chain_id": row["chain_id"],
                        "residue_index": row["residue_index"],
                        "residue_name": row["residue_name"],
                        "displacement": row["displacement"],
                        "coverage": disp["coverage"],
                        "n_valid": disp["n_valid"],
                    })
        return generator(displacement_data=displacement_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F05":
        # Displacement heatmap: long-format condition × residue frame
        if ref_condition is None:
            displacement_matrix = pd.DataFrame()
        else:
            rows = []
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
                if ref_struct is None:
                    continue
                disp = calculate_per_residue_displacement(
                    ref_struct,
                    target_struct,
                    alignment_atom=v3_config.structure.alignment_atom,
                )
                for row in disp["displacements"]:
                    rows.append({
                        "condition_id": condition_id,
                        "seed": seed,
                        "sample": sample,
                        "chain_id": row["chain_id"],
                        "residue_index": row["residue_index"],
                        "displacement": row["displacement"],
                    })
            displacement_matrix = pd.DataFrame(rows)
        return generator(displacement_matrix=displacement_matrix, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F06":
        # Contact map difference (prediction level, matched_sample pairing)
        contact_diff_data = []
        if ref_condition and dataset.predictions:
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
                if ref_struct is None:
                    continue
                try:
                    ref_map = calculate_contact_map(
                        ref_struct, threshold=v3_config.structure.contact_distance)
                    target_map = calculate_contact_map(
                        target_struct, threshold=v3_config.structure.contact_distance)
                    diff = calculate_contact_difference(ref_map, target_map)
                    diff["condition"] = condition_id
                    diff["seed"] = seed
                    diff["sample"] = sample
                    contact_diff_data.append(diff)
                except Exception as e:
                    logger.warning("[V3] F06 contact comparison failed for %s: %s",
                                   condition_id, e)
        return generator(contact_diff_data=contact_diff_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F07":
        # Contact change summary (seed-level aggregation of per-seed diffs)
        contact_change_data = []
        if ref_condition and dataset.predictions:
            for condition_id in sorted(dataset.predictions.keys()):
                if condition_id == ref_condition:
                    continue
                seeds_diffs = []
                for seed in sorted(dataset.predictions[condition_id].keys()):
                    seed_diffs = []
                    for sample, target_struct in sorted(
                            dataset.predictions[condition_id][seed].items()):
                        ref_struct = _get_reference_structure(
                            dataset, seed, sample, ref_condition)
                        if ref_struct is None:
                            continue
                        try:
                            ref_map = calculate_contact_map(
                                ref_struct, threshold=v3_config.structure.contact_distance)
                            target_map = calculate_contact_map(
                                target_struct, threshold=v3_config.structure.contact_distance)
                            diff = calculate_contact_difference(ref_map, target_map)
                            diff["condition"] = condition_id
                            seed_diffs.append(diff)
                        except Exception as e:
                            logger.warning("[V3] F07 contact diff failed (%s seed %s): %s",
                                           condition_id, seed, e)
                    if seed_diffs:
                        seeds_diffs.append(seed_diffs)

                for seed_idx, seed_diffs in enumerate(seeds_diffs):
                    contact_change_data.append({
                        "condition_id": condition_id,
                        "seed": seed_idx,
                        "n_gained": float(np.mean([d["n_gained"] for d in seed_diffs])),
                        "n_lost": float(np.mean([d["n_lost"] for d in seed_diffs])),
                        "n_total": float(np.mean([d["n_total"] for d in seed_diffs])),
                        "n_comparisons": len(seed_diffs),
                    })
        return generator(contact_change_data=contact_change_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F08":
        # Interface analysis (prediction level, matched_sample pairing)
        interface_data = []
        if ref_condition and dataset.predictions:
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
                if ref_struct is None:
                    continue
                try:
                    ref_ifaces = find_all_interfaces(
                        ref_struct, threshold=v3_config.structure.contact_distance)
                    tgt_ifaces = find_all_interfaces(
                        target_struct, threshold=v3_config.structure.contact_distance)
                    comparison = calculate_interface_comparison(ref_ifaces, tgt_ifaces)
                    for change in comparison["changes"]:
                        interface_data.append({
                            "condition_id": condition_id,
                            "seed": seed,
                            "sample": sample,
                            "interface_id": change["interface_id"],
                            "change": change["change"],
                            "ref_n_contacts": change.get("ref_n_contacts", 0),
                            "target_n_contacts": change.get("target_n_contacts", 0),
                            "n_contacts": change.get("target_n_contacts",
                                                     change.get("ref_n_contacts", 0)),
                        })
                except Exception as e:
                    logger.warning("[V3] F08 interface analysis failed for %s: %s",
                                   condition_id, e)
        return generator(interface_data=interface_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F09":
        # Interface change map (prediction level, matched_sample pairing)
        interface_changes = []
        if ref_condition and dataset.predictions:
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
                if ref_struct is None:
                    continue
                try:
                    ref_ifaces = find_all_interfaces(
                        ref_struct, threshold=v3_config.structure.contact_distance)
                    tgt_ifaces = find_all_interfaces(
                        target_struct, threshold=v3_config.structure.contact_distance)
                    comparison = calculate_interface_comparison(ref_ifaces, tgt_ifaces)
                    for change in comparison["changes"]:
                        if change["change"] == "unchanged":
                            continue
                        interface_changes.append({
                            "condition_id": condition_id,
                            "seed": seed,
                            "sample": sample,
                            "interface_id": change["interface_id"],
                            "change": change["change"],
                            "ref_n_contacts": change.get("ref_n_contacts", 0),
                            "target_n_contacts": change.get("target_n_contacts", 0),
                        })
                except Exception as e:
                    logger.warning("[V3] F09 interface change failed for %s: %s",
                                   condition_id, e)
        return generator(interface_changes=interface_changes, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F10":
        # Local geometry: config-driven regions/sites only (never invented)
        if not v3_config.sites:
            return {
                "status": "skip",
                "reason": "No sites configured (sites are config-driven)",
                "output_path": None,
                "n_observations": 0,
                "warnings": ["F10 requires v3_config.sites definitions"],
            }
        local_geometry_data = []
        if ref_condition and dataset.predictions:
            for site in v3_config.sites:
                site_label = site.get("label", "unnamed")
                for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                    if condition_id == ref_condition:
                        continue
                    ref_struct = _get_reference_structure(
                        dataset, seed, sample, ref_condition)
                    if ref_struct is None:
                        continue
                    try:
                        target_region = calculate_local_geometry(
                            target_struct, site, reference_structure=ref_struct)
                        ref_region = calculate_local_geometry(
                            ref_struct, site, reference_structure=ref_struct)
                        if target_region.status != "valid" or ref_region.status != "valid":
                            continue
                        local_geometry_data.append({
                            "condition_id": condition_id,
                            "seed": seed,
                            "sample": sample,
                            "region_label": site_label,
                            "local_rmsd_target": target_region.local_rmsd,
                            "local_rmsd_ref": ref_region.local_rmsd,
                            "local_rmsd": target_region.local_rmsd,
                            "local_plddt_target": target_region.local_plddt_mean,
                            "local_plddt_ref": ref_region.local_plddt_mean,
                            "n_atoms": target_region.n_atoms,
                            "n_residues": target_region.n_residues,
                        })
                    except Exception as e:
                        logger.warning("[V3] F10 local geometry failed (%s/%s): %s",
                                       site_label, condition_id, e)
        return generator(local_geometry_data=local_geometry_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F11":
        # Domain motion: config-driven regions only (never invented)
        if not v3_config.regions:
            return {
                "status": "skip",
                "reason": "No regions configured (regions are config-driven)",
                "output_path": None,
                "n_observations": 0,
                "warnings": ["F11 requires v3_config.regions definitions"],
            }
        domain_motion_data = []
        if ref_condition and dataset.predictions:
            for region in v3_config.regions:
                region_label = region.get("label", "unnamed")
                for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                    if condition_id == ref_condition:
                        continue
                    ref_struct = _get_reference_structure(
                        dataset, seed, sample, ref_condition)
                    if ref_struct is None:
                        continue
                    try:
                        target_region = calculate_local_geometry(
                            target_struct, region, reference_structure=ref_struct)
                        ref_region = calculate_local_geometry(
                            ref_struct, region, reference_structure=ref_struct)
                        if target_region.status != "valid" or ref_region.status != "valid":
                            continue
                        centroid_displacement = None
                        if (target_region.centroid is not None
                                and ref_region.centroid is not None):
                            centroid_displacement = float(np.linalg.norm(
                                np.asarray(target_region.centroid)
                                - np.asarray(ref_region.centroid)))
                        rg_change = None
                        if (target_region.radius_gyration is not None
                                and ref_region.radius_gyration is not None
                                and ref_region.radius_gyration != 0):
                            rg_change = (
                                target_region.radius_gyration
                                - ref_region.radius_gyration
                            )
                        domain_motion_data.append({
                            "condition_id": condition_id,
                            "seed": seed,
                            "sample": sample,
                            "region_label": region_label,
                            "centroid_displacement": centroid_displacement,
                            "local_rmsd_target": target_region.local_rmsd,
                            "local_rmsd_ref": ref_region.local_rmsd,
                            "rg_change": rg_change,
                        })
                    except Exception as e:
                        logger.warning("[V3] F11 domain motion failed (%s/%s): %s",
                                       region_label, condition_id, e)
        return generator(domain_motion_data=domain_motion_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id in ("F12", "F13", "F14"):
        # Clustering / similarity / MDS from the pairwise matrix
        if pairwise_matrix is None or pairwise_matrix.n_structures == 0:
            return {
                "status": "skip",
                "reason": "No pairwise matrix available",
                "output_path": None,
                "n_observations": 0,
                "warnings": ["Pairwise structural matrix unavailable"],
            }

        matrix = np.asarray(pairwise_matrix.matrix, dtype=float)
        # Fill invalid entries with the finite maximum so clustering can run;
        # validity is retained separately for reporting.
        if pairwise_matrix.valid is not None and pairwise_matrix.valid.any():
            finite_vals = matrix[pairwise_matrix.valid & ~np.eye(matrix.shape[0], dtype=bool)]
            finite_vals = finite_vals[np.isfinite(finite_vals)]
        else:
            finite_vals = matrix[np.isfinite(matrix)]
        fill_value = float(np.max(finite_vals)) if finite_vals.size else 0.0
        matrix_filled = np.where(np.isfinite(matrix), matrix, fill_value)
        np.fill_diagonal(matrix_filled, 0.0)

        predictions = list(pairwise_matrix.predictions)
        conditions = list(pairwise_matrix.conditions)
        seeds = list(pairwise_matrix.seeds)

        if fig_id == "F12":
            clustering = hierarchical_clustering(
                matrix_filled,
                n_clusters=v3_config.clustering.n_clusters,
                linkage=v3_config.clustering.linkage,
            )
            return generator(
                distance_matrix=matrix_filled,
                conditions=conditions,
                seeds=seeds,
                predictions=predictions,
                cluster_labels=clustering.get("labels"),
                **common_params,
            )

        if fig_id == "F13":
            return generator(
                distance_matrix=matrix_filled,
                predictions=predictions,
                conditions=conditions,
                **common_params,
            )

        # F14
        embedding = mds_embedding(
            matrix_filled,
            n_components=2,
            dissimilarity="precomputed",
        )
        coordinates = np.asarray(embedding.get("coordinates"))
        if coordinates.size == 0:
            return {
                "status": "skip",
                "reason": "MDS produced no coordinates",
                "output_path": None,
                "n_observations": 0,
                "warnings": ["MDS embedding failed"],
            }
        return generator(
            coordinates=coordinates,
            predictions=predictions,
            conditions=conditions,
            seeds=seeds,
            **common_params,
        )

    # ------------------------------------------------------------------
    if fig_id == "F15":
        # Confidence × geometry (prediction level, matched_sample pairing)
        confidence_geometry_data = []
        if ref_condition and dataset.predictions:
            structures = []
            rmsd_values = []
            displacement_values = []
            plddt_values = []
            pae_values = []
            contact_prob_values = []
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(
                    dataset, seed, sample, ref_condition)
                rmsd = None
                mean_disp = None
                if ref_struct is not None:
                    result = calculate_rmsd(
                        target_struct,
                        ref_struct,
                        alignment_atom=v3_config.structure.alignment_atom,
                        min_common_atoms=v3_config.structure.min_common_atoms,
                        min_sequence_identity=v3_config.structure.min_sequence_identity,
                        min_coverage=v3_config.structure.minimum_coverage,
                    )
                    rmsd = result["rmsd"]
                    if rmsd is not None:
                        disp = calculate_per_residue_displacement(
                            ref_struct,
                            target_struct,
                            alignment_atom=v3_config.structure.alignment_atom,
                        )
                        if disp["displacements"]:
                            mean_disp = float(np.mean(
                                [d["displacement"] for d in disp["displacements"]]))
                structures.append(target_struct)
                rmsd_values.append(rmsd)
                displacement_values.append(mean_disp)
                plddt_values.append(target_struct.plddt_mean)
                pae_values.append(target_struct.pae_mean)
                contact_prob_values.append(target_struct.contact_prob_mean)

            analysis = integrate_confidence_geometry(
                structures,
                rmsd_values=rmsd_values,
                displacement_values=displacement_values,
                plddt_values=plddt_values,
                pae_values=pae_values,
                contact_prob_values=contact_prob_values,
            )
            for obs in analysis.observations:
                confidence_geometry_data.append({
                    "prediction_id": obs.prediction_id,
                    "condition_id": obs.condition_id,
                    "seed": obs.seed,
                    "sample": obs.sample,
                    "rmsd_to_reference": obs.rmsd_to_reference,
                    "mean_displacement": obs.mean_displacement,
                    "plddt_mean": obs.plddt_mean,
                    "pae_mean": obs.pae_mean,
                    "contact_prob_mean": obs.contact_prob_mean,
                    "structural_category": obs.structural_category,
                    "confidence_category": obs.confidence_category,
                })
        return generator(confidence_geometry_data=confidence_geometry_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F16":
        # Confidence change vs structural change (seed-level deltas vs reference)
        delta_data = []
        matched = _get_matched_seed_results(dataset, ref_resolution, v3_config)
        ref_seed_metrics = dataset.seeds.get(ref_condition, {}) if ref_condition else {}
        for condition_id in sorted(matched.keys()):
            for msr in matched[condition_id]:
                if msr.rmsd_mean is None:
                    continue
                ref_metrics = ref_seed_metrics.get(msr.seed)
                if ref_metrics is None:
                    continue
                ref_plddt = ref_metrics.metrics.get("pLDDT_mean")
                if ref_plddt is None:
                    # Missing reference confidence: record structural delta only.
                    # pLDDT delta is left absent (no imputation).
                    delta_data.append({
                        "condition_id": condition_id,
                        "seed": msr.seed,
                        "delta_rmsd": msr.rmsd_mean,
                        "delta_plddt": None,
                    })
                else:
                    tgt_plddt = dataset.seeds.get(condition_id, {}).get(msr.seed)
                    tgt_plddt_val = tgt_plddt.metrics.get("pLDDT_mean") if tgt_plddt else None
                    delta_data.append({
                        "condition_id": condition_id,
                        "seed": msr.seed,
                        "delta_rmsd": msr.rmsd_mean,
                        "delta_plddt": (
                            (tgt_plddt_val - ref_plddt)
                            if tgt_plddt_val is not None else None
                        ),
                    })
        # F16 requires both deltas present per row; drop incomplete rows here
        # and record the exclusion count rather than imputing.
        complete = [d for d in delta_data if d["delta_plddt"] is not None]
        n_dropped = len(delta_data) - len(complete)
        if n_dropped:
            logger.warning("[V3] F16: %d rows excluded (missing pLDDT delta; "
                           "missingness preserved, no imputation)", n_dropped)
        return generator(delta_data=complete, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F17":
        # Seed reproducibility (seed-level unit of analysis)
        seed_repro_data = []
        matched = _get_matched_seed_results(dataset, ref_resolution, v3_config)
        for condition_id in sorted(matched.keys()):
            msr_list = matched[condition_id]
            if not msr_list:
                continue
            seed_values: Dict[int, List[float]] = {}
            seed_coverages: Dict[int, List[float]] = {}
            for msr in msr_list:
                vals = [v for v in msr.rmsd_values if v is not None and np.isfinite(v)]
                if vals:
                    seed_values[msr.seed] = vals
                if msr.coverage_values:
                    seed_coverages[msr.seed] = list(msr.coverage_values)
            repro = calculate_seed_reproducibility(
                seed_values,
                seed_coverages if seed_coverages else None,
                metric_id="rmsd_global_ca",
                condition_id=condition_id,
                reference_condition=ref_condition,
            )
            seed_repro_data.append({
                "metric_id": repro.metric_id,
                "condition_id": repro.condition_id,
                "reference_condition": repro.reference_condition,
                "mean": repro.mean,
                "median": repro.median,
                "std": repro.std,
                "iqr": repro.iqr,
                "n_seeds": repro.n_seeds,
                "n_valid_seeds": repro.n_valid_seeds,
                "n_comparisons_total": repro.n_comparisons_total,
                "n_comparisons_valid": repro.n_comparisons_valid,
                "direction_consistent": repro.direction_consistent,
                "direction": repro.direction,
                "mean_coverage": repro.mean_coverage,
                "status": repro.status,
            })
        return generator(seed_repro_data=seed_repro_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F18":
        # Effect sizes (seed-level via matched-seed results)
        effect_size_data = []
        matched = _get_matched_seed_results(dataset, ref_resolution, v3_config)
        for condition_id in sorted(matched.keys()):
            msr_list = matched[condition_id]
            if not msr_list:
                continue
            rmsd_values = []
            coverage_values = []
            for msr in msr_list:
                for i, v in enumerate(msr.rmsd_values):
                    if v is None or not np.isfinite(v):
                        continue
                    rmsd_values.append(v)
                    cov = (msr.coverage_values[i]
                           if i < len(msr.coverage_values) else None)
                    if cov is not None and np.isfinite(cov):
                        coverage_values.append(cov)
            effects = calculate_structural_effect_sizes(
                rmsd_values,
                coverage_values,
                condition_a=condition_id,
                condition_b=ref_condition,
                reference_condition=ref_condition,
            )
            for eff in effects:
                effect_size_data.append({
                    "metric_id": eff.metric_id,
                    "condition_id": eff.condition_a,
                    "reference_condition": eff.condition_b,
                    "estimate": eff.estimate,
                    "std_error": eff.std_error,
                    "ci_lower": eff.confidence_interval_lower,
                    "ci_upper": eff.confidence_interval_upper,
                    "n_valid": eff.n_valid,
                    "n_comparisons": eff.n_comparisons,
                    "hedges_g": eff.hedges_g,
                    "coverage_mean": eff.coverage_mean,
                    "status": eff.status,
                })
        return generator(effect_size_data=effect_size_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F19":
        # Factorial contrasts (metadata-driven; requires experiment metadata)
        metadata = dataset.experiment_metadata
        if not metadata or not metadata.get("attributes") or not metadata.get("conditions"):
            return {
                "status": "skip",
                "reason": "No experiment metadata (factor contrasts are "
                          "metadata-driven and were not defined)",
                "output_path": None,
                "n_observations": 0,
                "warnings": ["F19 requires experiment metadata with attributes"],
            }
        rmsd_values: Dict[str, List[float]] = {}
        matched = _get_matched_seed_results(dataset, ref_resolution, v3_config)
        meta_conds = metadata.get("conditions", {})
        for condition_id in sorted(matched.keys()):
            values = []
            for msr in matched[condition_id]:
                values.extend(
                    v for v in msr.rmsd_values
                    if v is not None and np.isfinite(v)
                )
            if values:
                cond_obj = dataset.conditions.get(condition_id)
                cond_name = cond_obj.condition_name if cond_obj else condition_id
                key = condition_id
                if condition_id not in meta_conds and cond_name in meta_conds:
                    key = cond_name
                rmsd_values[key] = values

        if ref_condition:
            ref_obj = dataset.conditions.get(ref_condition)
            ref_name = ref_obj.condition_name if ref_obj else ref_condition
            ref_key = ref_condition
            if ref_condition not in meta_conds and ref_name in meta_conds:
                ref_key = ref_name
            if ref_key not in rmsd_values and ref_key in meta_conds:
                n_ref_seeds = len(dataset.seeds.get(ref_condition, {})) or 1
                rmsd_values[ref_key] = [0.0] * n_ref_seeds
        if not rmsd_values:
            return {
                "status": "skip",
                "reason": "No RMSD values for factor contrasts",
                "output_path": None,
                "n_observations": 0,
                "warnings": [],
            }
        try:
            analysis = factorial_structural_analysis(
                rmsd_values, metadata, metric_name="rmsd_global_ca")
        except Exception as e:
            return {
                "status": "failed",
                "error": f"Factorial analysis failed: {e}",
                "output_path": None,
                "n_observations": 0,
                "warnings": [],
            }
        factorial_data = []
        for contrast in analysis.contrasts:
            factorial_data.append({
                "contrast_id": contrast.contrast_id,
                "contrast_type": contrast.contrast_type,
                "factors": "+".join(contrast.factors),
                "estimate": contrast.estimate,
                "std_error": contrast.std_error,
                "ci_lower": contrast.ci_lower,
                "ci_upper": contrast.ci_upper,
                "n_observations": contrast.n_observations,
                "n_conditions": contrast.n_conditions,
            })
        return generator(factorial_data=factorial_data, **common_params)

    # ------------------------------------------------------------------
    if fig_id == "F20":
        # Structure-confidence relationship matrix (prediction level)
        relationship_data = []
        if ref_condition and dataset.predictions:
            for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
                if condition_id == ref_condition:
                    continue
                ref_struct = _get_reference_structure(
                    dataset, seed, sample, ref_condition)
                rmsd = None
                mean_disp = None
                if ref_struct is not None:
                    result = calculate_rmsd(
                        target_struct,
                        ref_struct,
                        alignment_atom=v3_config.structure.alignment_atom,
                        min_common_atoms=v3_config.structure.min_common_atoms,
                        min_sequence_identity=v3_config.structure.min_sequence_identity,
                        min_coverage=v3_config.structure.minimum_coverage,
                    )
                    rmsd = result["rmsd"]
                    if rmsd is not None:
                        disp = calculate_per_residue_displacement(
                            ref_struct,
                            target_struct,
                            alignment_atom=v3_config.structure.alignment_atom,
                        )
                        if disp["displacements"]:
                            mean_disp = float(np.mean(
                                [d["displacement"] for d in disp["displacements"]]))
                row: Dict[str, Any] = {
                    "prediction_id": target_struct.prediction_id,
                    "condition_id": condition_id,
                    "seed": seed,
                    "sample": sample,
                }
                if rmsd is not None:
                    row["rmsd_to_reference"] = rmsd
                if mean_disp is not None:
                    row["displacement_mean"] = mean_disp
                if target_struct.plddt_mean is not None:
                    row["plddt_mean"] = target_struct.plddt_mean
                if target_struct.pae_mean is not None:
                    row["pae_mean"] = target_struct.pae_mean
                if target_struct.contact_prob_mean is not None:
                    row["contact_prob_mean"] = target_struct.contact_prob_mean
                relationship_data.append(row)
        return generator(relationship_data=relationship_data, **common_params)

    # ------------------------------------------------------------------
    return {
        "status": "skip",
        "reason": f"Figure {fig_id} not yet implemented",
        "output_path": None,
        "n_observations": 0,
        "warnings": [],
    }


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def _generate_v3_tables(
    dataset: Any,
    all_structures: List[Any],
    pairwise_matrix: Any,
    ref_resolution: Dict[str, Any],
    v3_config: V3Config,
    tables_dir: Path,
    results: Dict[str, Any],
) -> None:
    """Generate V3 output tables."""
    ref_condition = ref_resolution.get("reference_condition")
    written: Dict[str, str] = {}

    # --- structural_summary.csv (prediction-level confidence metrics) ---
    rows = []
    for condition_id, seed, sample, structure in _iter_predictions(dataset):
        rows.append({
            "prediction_id": structure.prediction_id,
            "condition_id": condition_id,
            "condition_name": structure.condition_id,
            "seed": seed,
            "sample": sample,
            "plddt_mean": structure.plddt_mean,
            "plddt_min": structure.plddt_min,
            "plddt_max": structure.plddt_max,
            "plddt_median": structure.plddt_median,
            "pae_mean": structure.pae_mean,
            "contact_prob_mean": structure.contact_prob_mean,
            "n_chains": structure.n_chains,
            "n_entities": structure.n_entities,
            "parse_status": structure.parse_status,
        })
    if rows:
        path = tables_dir / "structural_summary.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        written["structural_summary"] = str(path)

    # --- per_residue_displacement.csv (prediction level, matched_sample) ---
    rows = []
    if ref_condition and dataset.predictions:
        for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
            if condition_id == ref_condition:
                continue
            ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
            if ref_struct is None:
                continue
            try:
                disp = calculate_per_residue_displacement(
                    ref_struct,
                    target_struct,
                    alignment_atom=v3_config.structure.alignment_atom,
                )
            except Exception as e:
                logger.warning("[V3] Table displacement failed for %s: %s",
                               condition_id, e)
                continue
            for row in disp["displacements"]:
                rows.append({
                    "condition_id": condition_id,
                    "seed": seed,
                    "sample": sample,
                    "chain_id": row["chain_id"],
                    "residue_index": row["residue_index"],
                    "residue_name": row["residue_name"],
                    "displacement": row["displacement"],
                    "atom_name": row.get("atom_name",
                                         v3_config.structure.alignment_atom),
                    "n_valid": disp["n_valid"],
                    "coverage": disp["coverage"],
                })
    if rows:
        path = tables_dir / "per_residue_displacement.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        written["per_residue_displacement"] = str(path)

    # --- contact_changes.csv (prediction level, matched_sample) ---
    rows = []
    if ref_condition and dataset.predictions:
        for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
            if condition_id == ref_condition:
                continue
            ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
            if ref_struct is None:
                continue
            try:
                ref_map = calculate_contact_map(
                    ref_struct, threshold=v3_config.structure.contact_distance)
                target_map = calculate_contact_map(
                    target_struct, threshold=v3_config.structure.contact_distance)
                diff = calculate_contact_difference(ref_map, target_map)
            except Exception as e:
                logger.warning("[V3] Table contact diff failed for %s: %s",
                               condition_id, e)
                continue
            rows.append({
                "condition_id": condition_id,
                "seed": seed,
                "sample": sample,
                "reference_condition": ref_condition,
                "n_gained": diff["n_gained"],
                "n_lost": diff["n_lost"],
                "n_unchanged": diff["n_unchanged"],
                "n_total": diff["n_total"],
                "pct_changed": diff["pct_changed"],
                "ref_n_contacts": diff["ref_n_contacts"],
                "target_n_contacts": diff["target_n_contacts"],
                "contact_threshold": v3_config.structure.contact_distance,
            })
    if rows:
        path = tables_dir / "contact_changes.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        written["contact_changes"] = str(path)

    # --- interface_contacts.csv (prediction level, matched_sample) ---
    rows = []
    if ref_condition and dataset.predictions:
        for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
            if condition_id == ref_condition:
                continue
            ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
            if ref_struct is None:
                continue
            try:
                ref_ifaces = find_all_interfaces(
                    ref_struct, threshold=v3_config.structure.contact_distance)
                tgt_ifaces = find_all_interfaces(
                    target_struct, threshold=v3_config.structure.contact_distance)
                comparison = calculate_interface_comparison(ref_ifaces, tgt_ifaces)
            except Exception as e:
                logger.warning("[V3] Table interface diff failed for %s: %s",
                               condition_id, e)
                continue
            for change in comparison["changes"]:
                rows.append({
                    "condition_id": condition_id,
                    "seed": seed,
                    "sample": sample,
                    "reference_condition": ref_condition,
                    "interface_id": change["interface_id"],
                    "change": change["change"],
                    "ref_n_contacts": change.get("ref_n_contacts", 0),
                    "target_n_contacts": change.get("target_n_contacts", 0),
                })
    if rows:
        path = tables_dir / "interface_contacts.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        written["interface_contacts"] = str(path)

    # --- structural_clusters.csv (prediction level; saved from pairwise matrix) ---
    if pairwise_matrix is not None and pairwise_matrix.n_structures > 0:
        try:
            matrix = np.asarray(pairwise_matrix.matrix, dtype=float)
            if pairwise_matrix.valid is not None and pairwise_matrix.valid.any():
                finite_vals = matrix[pairwise_matrix.valid
                                     & ~np.eye(matrix.shape[0], dtype=bool)]
                finite_vals = finite_vals[np.isfinite(finite_vals)]
            else:
                finite_vals = matrix[np.isfinite(matrix)]
            fill_value = float(np.max(finite_vals)) if finite_vals.size else 0.0
            matrix_filled = np.where(np.isfinite(matrix), matrix, fill_value)
            np.fill_diagonal(matrix_filled, 0.0)
            clustering = hierarchical_clustering(
                matrix_filled,
                n_clusters=v3_config.clustering.n_clusters,
                linkage=v3_config.clustering.linkage,
            )
            cluster_rows = []
            for i, pred in enumerate(pairwise_matrix.predictions):
                cluster_rows.append({
                    "prediction_id": pred,
                    "condition_id": pairwise_matrix.conditions[i],
                    "seed": pairwise_matrix.seeds[i],
                    "sample": pairwise_matrix.samples[i],
                    "cluster": int(clustering["labels"][i]),
                    "n_clusters": int(clustering["n_clusters"]),
                    "clustering_method": clustering.get("method", "hierarchical"),
                    "linkage": clustering.get("linkage", v3_config.clustering.linkage),
                })
            path = tables_dir / "structural_clusters.csv"
            pd.DataFrame(cluster_rows).to_csv(path, index=False)
            written["structural_clusters"] = str(path)
        except Exception as e:
            logger.warning("[V3] structural_clusters table failed: %s", e)

    # --- confidence_geometry.csv (prediction level, matched_sample) ---
    rows = []
    if ref_condition and dataset.predictions:
        for condition_id, seed, sample, target_struct in _iter_predictions(dataset):
            if condition_id == ref_condition:
                continue
            ref_struct = _get_reference_structure(dataset, seed, sample, ref_condition)
            rmsd = None
            mean_disp = None
            if ref_struct is not None:
                result = calculate_rmsd(
                    target_struct,
                    ref_struct,
                    alignment_atom=v3_config.structure.alignment_atom,
                    min_common_atoms=v3_config.structure.min_common_atoms,
                    min_sequence_identity=v3_config.structure.min_sequence_identity,
                    min_coverage=v3_config.structure.minimum_coverage,
                )
                rmsd = result["rmsd"]
                if rmsd is not None:
                    disp = calculate_per_residue_displacement(
                        ref_struct,
                        target_struct,
                        alignment_atom=v3_config.structure.alignment_atom,
                    )
                    if disp["displacements"]:
                        mean_disp = float(np.mean(
                            [d["displacement"] for d in disp["displacements"]]))
            rows.append({
                "prediction_id": target_struct.prediction_id,
                "condition_id": condition_id,
                "seed": seed,
                "sample": sample,
                "reference_condition": ref_condition,
                "rmsd_to_reference": rmsd,
                "mean_displacement": mean_disp,
                "plddt_mean": target_struct.plddt_mean,
                "pae_mean": target_struct.pae_mean,
                "contact_prob_mean": target_struct.contact_prob_mean,
            })
    if rows:
        path = tables_dir / "confidence_geometry.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        written["confidence_geometry"] = str(path)

    results["tables"] = written
    logger.info("[V3] Wrote %d tables to %s", len(written), tables_dir)


# ---------------------------------------------------------------------------
# Manifest / report
# ---------------------------------------------------------------------------

def _config_to_dict(v3_config: V3Config) -> Dict[str, Any]:
    """Serialize V3Config for the manifest without requiring to_dict()."""
    try:
        return v3_config.to_dict()
    except Exception:
        return {
            "enabled": v3_config.enabled,
            "reference": v3_config.reference,
            "figures": v3_config.figures,
            "structure": vars(v3_config.structure),
            "clustering": vars(v3_config.clustering),
            "sites": v3_config.sites,
            "regions": v3_config.regions,
        }


def _safe_float(value: Any) -> Optional[float]:
    """Convert NaN/inf to None for JSON serialization."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _pairwise_cache_key(structures: List[Any], v3_config: V3Config) -> str:
    """Deterministic cache key from prediction identities + structural config."""
    h = hashlib.sha256()
    h.update(json.dumps(vars(v3_config.structure), sort_keys=True, default=str)
             .encode("utf-8"))
    for s in sorted(structures, key=lambda x: x.prediction_id):
        h.update(s.prediction_id.encode("utf-8"))
        try:
            h.update(str(Path(s.source_path).stat().st_mtime_ns).encode("utf-8"))
        except OSError:
            pass
    return h.hexdigest()[:16]


def _write_v3_manifest(
    metadata_dir: Path,
    results: Dict[str, Any],
) -> None:
    """Write V3 figure manifest."""
    manifest = {
        "version": "3.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "run_id": (results.get("config", {}) or {}).get("run_id", ""),
        "status": results.get("status", "unknown"),
        "figures": [],
        "summary": results.get("summary", {}),
        "dataset_summary": results.get("dataset_summary", {}),
        "structural_qc": _jsonable(results.get("structural_qc", {})),
        "pairwise_matrix_stats": results.get("pairwise_matrix_stats", {}),
        "reference_resolution": results.get("reference_resolution", {}),
        "tables": results.get("tables", {}),
        "warnings": results.get("warnings", []),
        "errors": results.get("errors", []),
        "skipped": results.get("skipped", []),
    }

    for fig_id, fig_result in results.get("figures", {}).items():
        entry = {
            "figure_id": fig_id,
            "status": fig_result.get("status", "unknown"),
            "output_path": fig_result.get("output_path"),
            "n_observations": fig_result.get("n_observations", 0),
            "warnings": fig_result.get("warnings", []),
            "skipped_reason": (fig_result.get("reason")
                               if fig_result.get("status") == "skip" else None),
            "error": (fig_result.get("error")
                      if fig_result.get("status") == "failed" else None),
        }
        manifest["figures"].append(entry)

    manifest_path = metadata_dir / "figure_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, default=str)


def _jsonable(obj: Any) -> Any:
    """Best-effort conversion of numpy types for JSON serialization."""
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return _safe_float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _generate_v3_report(
    report_dir: Path,
    results: Dict[str, Any],
) -> None:
    """Generate V3 visualization report."""
    report_path = report_dir / "V3_VISUALIZATION_REPORT.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# V3 Visualization Report\n\n")
        f.write(f"**Generated:** {datetime.now(timezone.utc).isoformat()}\n\n")

        f.write("## Pipeline Status\n\n")
        f.write(f"- Status: `{results.get('status', 'unknown')}`\n")
        f.write(f"- Elapsed: {results.get('elapsed_s', 0):.1f}s\n\n")

        dataset_summary = results.get("dataset_summary", {})
        if dataset_summary:
            f.write("## Dataset\n\n")
            f.write(f"- Conditions: {dataset_summary.get('n_conditions', 0)}\n")
            f.write(f"- Predictions: {dataset_summary.get('n_predictions', 0)}\n\n")

        ref = results.get("reference_resolution", {})
        if ref:
            f.write("## Reference\n\n")
            f.write(f"- Reference condition: "
                    f"`{ref.get('reference_condition', 'n/a')}`\n")
            f.write(f"- Paired conditions: "
                    f"{len(ref.get('paired_conditions', []))}\n\n")

        summary = results.get("summary", {})
        f.write("## Figure Summary\n\n")
        f.write(f"- Requested: {summary.get('n_figures_requested', 0)}\n")
        f.write(f"- Success: {summary.get('n_figures_success', 0)}\n")
        f.write(f"- Skipped: {summary.get('n_figures_skipped', 0)}\n")
        f.write(f"- Failed: {summary.get('n_figures_failed', 0)}\n\n")

        f.write("## Figures\n\n")
        f.write("| Figure ID | Status | Observations | Output |\n")
        f.write("|-----------|--------|-------------|--------|\n")

        for fig_id, fig_result in results.get("figures", {}).items():
            status = fig_result.get("status", "unknown")
            n_obs = fig_result.get("n_observations", 0)
            output = fig_result.get("output_path", "")

            f.write(f"| {fig_id} | {status} | {n_obs} | {output or '—'} |\n")

        f.write("\n## Tables\n\n")
        tables = results.get("tables", {})
        if tables:
            for name, path in tables.items():
                f.write(f"- `{name}`: {path}\n")
        else:
            f.write("- None written\n")

        f.write("\n## Warnings\n\n")
        for warning in results.get("warnings", []):
            f.write(f"- {warning}\n")

        if results.get("warnings"):
            f.write("\n")

        f.write("## Errors\n\n")
        for error in results.get("errors", []):
            f.write(f"- {error}\n")

        if results.get("errors"):
            f.write("\n")

        f.write("## Skipped Figures\n\n")
        for skipped in results.get("skipped", []):
            f.write(f"- {skipped.get('figure_id', '?')}: "
                    f"{skipped.get('reason', '?')}\n")

        if results.get("skipped"):
            f.write("\n")

    logger.info("[V3] Report written: %s", report_path)
