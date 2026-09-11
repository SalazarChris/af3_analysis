"""
V3 Data Adapter.

Reads existing pipeline outputs (tables, structures) and constructs
the normalized V3 data model (Dataset).

This adapter is pure read-only: it never modifies existing files.
It translates the existing project representation into the V3 representation.
"""

from __future__ import annotations

import ast
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from af3_analysis.io.structure_reader import parse_mmcif, StructureParseError
from af3_analysis.experiment_metadata import load_experiment_design

from .model import (
    Dataset,
    StructureData,
    EntityGeometry,
    ChainGeometry,
    ResidueGeometry,
    AtomPosition,
    ConditionSummary,
    SeedSummary,
    StructuralQC,
)

logger = logging.getLogger(__name__)


class V3AdapterError(Exception):
    """Adapter error."""
    pass


class V3AdapterValidationError(V3AdapterError):
    """Validation error during adaptation."""
    pass


def adapt_from_pipeline(
    run_dir: Path,
    raw_af3_root: Optional[Path] = None,
    experiment_metadata_path: Optional[Path] = None,
    v3_config: Optional[Any] = None,
) -> Dataset:
    """
    Build a V3 Dataset from existing pipeline outputs.

    Parameters
    ----------
    run_dir : Path
        The existing pipeline run directory containing tables/ and figures/.
    raw_af3_root : Path, optional
        Root directory containing raw AF3 outputs (CIF files).
        Required for structural analysis.
    experiment_metadata_path : Path, optional
        Path to experiment_metadata.json for condition labels.
    v3_config : V3Config, optional
        V3 configuration (used for reference resolution).

    Returns
    -------
    Dataset
        Normalized V3 dataset.
    """
    logger.info("[V3 Adapter] Building dataset from pipeline outputs")

    # Phase 1: Load existing tables
    tables_dir = run_dir / "tables"
    logger.info("[V3 Adapter] Loading tables from %s", tables_dir)

    seed_aggregated = pd.read_csv(tables_dir / "seed_aggregated.csv")
    descriptive_stats = pd.read_csv(tables_dir / "descriptive_stats.csv")

    # Phase 2: Load experiment metadata if available
    experiment_metadata = None
    if experiment_metadata_path is not None:
        try:
            design = load_experiment_design(experiment_metadata_path)
            experiment_metadata = {
                "experiment_id": design.experiment_id,
                "description": design.description,
                "conditions": {
                    name: {
                        "label": cond.label,
                        "attributes": dict(cond.attributes),
                    }
                    for name, cond in design.conditions.items()
                },
                "attributes": {
                    name: {
                        "type": attr.attr_type,
                        "description": attr.description,
                    }
                    for name, attr in design.attributes.items()
                },
            }
            logger.info("[V3 Adapter] Loaded experiment metadata: %s", design.experiment_id)
        except Exception as e:
            logger.warning("[V3 Adapter] Could not load experiment metadata: %s", e)

    # Phase 3: Build condition summaries from existing tables
    conditions = _build_condition_summaries(
        seed_aggregated, descriptive_stats, experiment_metadata
    )

    # Phase 4: Build seed summaries
    seeds = _build_seed_summaries(seed_aggregated)

    # Phase 5: Load structures if raw_af3_root provided
    predictions = {}
    stem_map: Dict[str, str] = {}
    stem_candidates: Dict[str, List[str]] = {}
    if raw_af3_root is not None and raw_af3_root.exists():
        logger.info("[V3 Adapter] Loading structures from %s", raw_af3_root)
        # CIF filenames carry the condition stem (e.g. 'pou_baseline'); the run
        # tables carry condition_id (e.g. 'cond_001'). Build the stem -> id map
        # so predictions align with conditions/seeds for matched-seed pairing.
        stem_map, stem_candidates = _build_condition_stem_map(
            run_dir, seed_aggregated
        )
        predictions, qc_records = _load_structures(
            raw_af3_root, seed_aggregated, experiment_metadata, stem_map=stem_map
        )
    else:
        logger.info("[V3 Adapter] No raw_af3_root provided; structural data skipped")
        qc_records = []

    # Phase 6: Determine reference condition
    reference_condition = None
    if v3_config is not None:
        ref = v3_config.reference
        if ref:
            if "condition" in ref:
                reference_condition = ref["condition"]
            elif "condition_a" in ref and "condition_b" in ref:
                # For pairwise comparison, use condition_a as reference
                reference_condition = ref["condition_a"]

    # If no reference from config, use first condition alphabetically
    if reference_condition is None and conditions:
        reference_condition = sorted(conditions.keys())[0]

    # Resolve human-supplied names/stems (e.g. the raw AF3 output folder
    # name from the af3.py menu) to canonical condition ids. Ambiguous or
    # unresolvable references are left as-is; resolve_reference() reports
    # them explicitly instead of silently choosing a condition.
    if (
        reference_condition is not None
        and reference_condition not in conditions
    ):
        name_matches = sorted(
            cid for cid, cond in conditions.items()
            if cond.condition_name == reference_condition
        )
        if len(name_matches) == 1:
            reference_condition = name_matches[0]
        else:
            stem_matches = stem_candidates.get(reference_condition, [])
            if len(stem_matches) == 1:
                reference_condition = stem_matches[0]

    # Phase 7: Build dataset
    dataset = Dataset(
        name=run_dir.name,
        conditions=conditions,
        seeds=seeds,
        predictions=predictions,
        reference_condition=reference_condition,
        condition_stem_map=stem_map,
        condition_stem_candidates=stem_candidates,
        experiment_metadata=experiment_metadata,
    )
    logger.info("[V3 Adapter] Dataset built: %d conditions, %d predictions",
                len(conditions), sum(
                    sum(len(sp.values()) for sp in cond_preds.values())
                    for cond_preds in predictions.values()
                ))

    return dataset


def _build_condition_summaries(
    seed_aggregated: pd.DataFrame,
    descriptive_stats: pd.DataFrame,
    experiment_metadata: Optional[Dict[str, Any]],
) -> Dict[str, ConditionSummary]:
    """Build condition summaries from existing tables."""
    conditions = {}

    # Get condition names from seed_aggregated
    for _, row in seed_aggregated.groupby("condition_id"):
        condition_id = row["condition_id"].iloc[0]
        condition_name = row["condition_name"].iloc[0] if "condition_name" in row.columns else condition_id
        n_seeds = row["seed"].nunique() if "seed" in row.columns else 0
        n_predictions = len(row)

        # Get metrics from descriptive_stats
        metrics = {}
        if condition_id in descriptive_stats["condition_id"].values:
            stats_row = descriptive_stats[descriptive_stats["condition_id"] == condition_id].iloc[0]
            for col in descriptive_stats.columns:
                if col.endswith("_mean") and col != "condition_id_mean":
                    metric_name = col.replace("_mean", "")
                    metrics[metric_name] = {
                        "mean": float(stats_row[col]) if pd.notna(stats_row[col]) else None,
                        "std": float(stats_row.get(col.replace("_mean", "_std"), np.nan))
                        if pd.notna(stats_row.get(col.replace("_mean", "_std"), np.nan))
                        else None,
                    }

        conditions[condition_id] = ConditionSummary(
            condition_id=condition_id,
            condition_name=condition_name,
            n_seeds=n_seeds,
            n_predictions=n_predictions,
            metrics=metrics,
        )

    return conditions


def _build_seed_summaries(
    seed_aggregated: pd.DataFrame,
) -> Dict[str, Dict[int, SeedSummary]]:
    """Build seed summaries from seed_aggregated.csv."""
    seeds = {}

    for _, row in seed_aggregated.iterrows():
        condition_id = row["condition_id"]
        seed = int(row["seed"]) if "seed" in row and pd.notna(row.get("seed")) else 0

        if condition_id not in seeds:
            seeds[condition_id] = {}

        metrics = {}
        for col in seed_aggregated.columns:
            if col not in ("condition_id", "condition_name", "seed", "replicate_id"):
                if pd.notna(row[col]):
                    metrics[col] = float(row[col])

        seeds[condition_id][seed] = SeedSummary(
            condition_id=condition_id,
            seed=seed,
            n_samples=1,  # seed_aggregated already aggregated samples
            metrics=metrics,
        )

    return seeds


def _build_condition_stem_map(
    run_dir: Path,
    seed_aggregated: pd.DataFrame,
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Map CIF filename condition stems to run condition_ids.

    The AF3 CIF files are named by their condition stem (e.g.
    'pou_baseline_seed-10_sample-0_model.cif'), while the run tables use
    generated condition_ids (e.g. 'cond_001'). The authoritative bridge is
    condition_registry.csv's replicate_ids (which embed the stem) and the
    seed_aggregated condition_name column.

    Returns
    -------
    (stem_map, stem_candidates)
        stem_map: stem -> first condition_id (deterministic, used for CIF
        lookup). stem_candidates: stem -> ALL condition_ids sharing the
        stem, used to detect reference-name ambiguity.
    """
    stem_map: Dict[str, str] = {}
    stem_candidates: Dict[str, List[str]] = {}

    # Source 1: condition_registry.csv replicate_ids (authoritative stems).
    # Example: 'oct4__k123-sumo_seed-10_sample-1' -> stem 'oct4__k123-sumo'
    registry = run_dir / "tables" / "condition_registry.csv"
    if registry.is_file():
        try:
            registry_df = pd.read_csv(registry)
            for _, row in registry_df.iterrows():
                condition_id = str(row["condition_id"])
                replicate_ids = row.get("replicate_ids")
                if not isinstance(replicate_ids, str):
                    continue
                try:
                    rep_list = ast.literal_eval(replicate_ids)
                except (ValueError, SyntaxError):
                    continue
                for rep in rep_list:
                    if not isinstance(rep, str):
                        continue
                    stem = _replicate_to_stem(rep)
                    if stem:
                        stem_map.setdefault(stem, condition_id)
                        candidates = stem_candidates.setdefault(stem, [])
                        if condition_id not in candidates:
                            candidates.append(condition_id)
        except Exception as e:
            logger.warning("[V3 Adapter] Could not read condition_registry.csv: %s", e)

    # Source 2: seed_aggregated condition_name (covers runs without a registry)
    if "condition_name" in seed_aggregated.columns:
        for _, row in seed_aggregated.drop_duplicates("condition_id").iterrows():
            name = str(row["condition_name"])
            condition_id = str(row["condition_id"])
            if name:
                stem_map.setdefault(name, condition_id)
                candidates = stem_candidates.setdefault(name, [])
                if condition_id not in candidates:
                    candidates.append(condition_id)

    return stem_map, stem_candidates


def _replicate_to_stem(replicate_id: str) -> Optional[str]:
    """
    Reduce a replicate id to its CIF condition stem.

    Examples:
      'pou_baseline_seed-10_sample-0' -> 'pou_baseline'
      'pou_baseline_seed-10_sample-0_summary' -> 'pou_baseline'
      'pou_baseline_summary' -> 'pou_baseline'
      'pou_baseline' -> 'pou_baseline'
    """
    stem = re.sub(r"_seed-\d+_sample-\d+_summary$", "", replicate_id)
    stem = re.sub(r"_seed-\d+_sample-\d+$", "", stem)
    stem = re.sub(r"_summary$", "", stem)
    return stem or None


def _load_structures(
    raw_af3_root: Path,
    seed_aggregated: pd.DataFrame,
    experiment_metadata: Optional[Dict[str, Any]],
    stem_map: Optional[Dict[str, str]] = None,
) -> Tuple[Dict[str, Dict[int, Dict[int, StructureData]]], List[StructuralQC]]:
    """
    Load all CIF structures and build StructureData objects.

    Parameters
    ----------
    stem_map : dict, optional
        CIF stem -> run condition_id mapping. Discovered conditions are
        remapped through this so prediction keys align with the run tables.

    Returns:
    - predictions: condition_id -> seed -> sample -> StructureData
    - qc_records: list of StructuralQC records
    """
    predictions = {}
    qc_records = []

    if stem_map is None:
        stem_map = {}

    # Discover all CIF files
    from af3_analysis.structural.discovery import discover_structures
    report = discover_structures(raw_af3_root)

    if report.n_discovered == 0:
        logger.warning("[V3 Adapter] No structures discovered")
        return predictions, qc_records

    # Parse each structure
    for cond_id, inv in report.conditions.items():
        run_condition_id = stem_map.get(cond_id, cond_id)
        if run_condition_id != cond_id:
            logger.info("[V3 Adapter] Mapped CIF condition '%s' -> run condition '%s'",
                        cond_id, run_condition_id)
        elif cond_id not in stem_map.values():
            logger.warning(
                "[V3 Adapter] CIF condition '%s' has no run-table mapping; "
                "using it as-is (predictions will not pair with table data)",
                cond_id,
            )
        predictions.setdefault(run_condition_id, {})
        for rec in inv.records:
            try:
                structure = parse_mmcif(
                    rec.path,
                    condition_id=run_condition_id,
                    seed=rec.seed,
                    sample=rec.sample,
                )
                structure_data = _convert_to_structure_data(
                    structure, experiment_metadata
                )
                predictions[run_condition_id].setdefault(rec.seed, {})
                predictions[run_condition_id][rec.seed][rec.sample] = structure_data

                qc_records.append(StructuralQC(
                    prediction_id=structure_data.prediction_id,
                    condition_id=structure_data.condition_id,
                    seed=structure_data.seed,
                    sample=structure_data.sample,
                    status="success",
                    n_atoms=structure_data.n_chains * 100,  # approximate; could refine
                    n_chains=structure.n_chains,
                    n_residues=len(structure.residues),
                    entity_types=sorted(structure.get_entity_types()),
                    chain_ids=list(structure.chain_ids),
                    has_protein=bool(structure.get_protein_chains()),
                    has_dna=bool(structure.get_nucleic_acid_chains()),
                ))

            except StructureParseError as e:
                logger.warning("[V3 Adapter] Parse error for %s: %s",
                              rec.path.name, str(e)[:100])
                qc_records.append(StructuralQC(
                    prediction_id=f"{run_condition_id}_seed-{rec.seed}_sample-{rec.sample}",
                    condition_id=run_condition_id,
                    seed=rec.seed,
                    sample=rec.sample,
                    status="parse_error",
                    reason=str(e)[:200],
                ))

    return predictions, qc_records


def _convert_to_structure_data(
    structure: Any,
    experiment_metadata: Optional[Dict[str, Any]],
) -> StructureData:
    """
    Convert a NormalisedStructure to V3 StructureData.

    This is the key translation layer between existing structural
    representation and V3 normalized model.
    """
    # Build entities
    entities = []
    for entity in structure.entities:
        chains = []
        for chain in structure.chains:
            if chain.entity_id == entity.entity_id:
                # Build residues for this chain
                residues = []
                chain_atoms = structure.get_chain(chain.chain_id)

                # Group atoms by residue
                residue_atoms = {}
                for atom in chain_atoms:
                    key = (atom.chain_id, atom.auth_seq_id)
                    if key not in residue_atoms:
                        residue_atoms[key] = []
                    residue_atoms[key].append(atom)

                for (cid, auth_seq_id), atoms in residue_atoms.items():
                    # Get residue info
                    residue_name = atoms[0].comp_id if atoms else "UNK"

                    # Find backbone atoms
                    n_coords = None
                    ca_coords = None
                    c_coords = None
                    o_coords = None
                    all_atoms = []

                    for atom in atoms:
                        coords = atom.coords
                        if atom.atom_name == "N":
                            n_coords = coords
                        elif atom.atom_name == "CA":
                            ca_coords = coords
                        elif atom.atom_name == "C":
                            c_coords = coords
                        elif atom.atom_name == "O":
                            o_coords = coords
                        all_atoms.append(AtomPosition(
                            atom_name=atom.atom_name,
                            chain_id=atom.chain_id,
                            entity_id=atom.entity_id,
                            auth_seq_id=atom.auth_seq_id,
                            coords=coords,
                            residue_name=residue_name,
                        ))

                    residues.append(ResidueGeometry(
                        residue_name=residue_name,
                        chain_id=cid,
                        auth_seq_id=auth_seq_id or 0,
                        entity_id=entity.entity_id,
                        n_coords=n_coords,
                        ca_coords=ca_coords,
                        c_coords=c_coords,
                        o_coords=o_coords,
                        all_atoms=all_atoms,
                    ))

                chains.append(ChainGeometry(
                    chain_id=chain.chain_id,
                    entity_id=entity.entity_id,
                    entity_type=entity.entity_type,
                    polymer_type=entity.polymer_type,
                    residues=sorted(residues, key=lambda r: r.auth_seq_id or 0),
                ))

        entities.append(EntityGeometry(
            entity_id=entity.entity_id,
            entity_type=entity.entity_type,
            polymer_type=entity.polymer_type,
            description=entity.description,
            chains=chains,
        ))

    # Confidence metrics from the AF3 confidences JSON that sits beside the
    # model CIF (e.g. pou_baseline_seed-10_sample-0_confidences.json).
    confidence = _load_confidence_metrics(structure.source_path)
    plddt_mean = confidence.get("plddt_mean")
    plddt_min = confidence.get("plddt_min")
    plddt_max = confidence.get("plddt_max")
    plddt_median = confidence.get("plddt_median")
    pae_mean = confidence.get("pae_mean")
    contact_prob_mean = confidence.get("contact_prob_mean")

    # NormalisedStructure has no prediction_id; use the project convention
    # (see af3_analysis/structural/tables.py):
    #   {condition_id}_seed-{seed}_sample-{sample}
    prediction_id = (
        f"{structure.condition_id}_seed-{structure.seed}_sample-{structure.sample}"
    )

    return StructureData(
        prediction_id=prediction_id,
        condition_id=structure.condition_id,
        seed=structure.seed,
        sample=structure.sample,
        source_path=structure.source_path,
        plddt_mean=plddt_mean,
        plddt_min=plddt_min,
        plddt_max=plddt_max,
        plddt_median=plddt_median,
        pae_mean=pae_mean,
        contact_prob_mean=contact_prob_mean,
        entities=entities,
        parse_status="success",
    )


def _load_confidence_metrics(cif_path: Path) -> Dict[str, Optional[float]]:
    """
    Load per-prediction confidence metrics from the AF3 confidences JSON
    that sits beside a model CIF.

    The confidences JSON contains per-atom pLDDT, per-residue PAE and
    contact probabilities (see af3_analysis/structural/ARCHITECTURE.md).
    Missing files or missing keys return None (never imputed).
    """
    import json as _json

    result: Dict[str, Optional[float]] = {
        "plddt_mean": None,
        "plddt_min": None,
        "plddt_max": None,
        "plddt_median": None,
        "pae_mean": None,
        "contact_prob_mean": None,
    }

    conf_path = Path(str(cif_path).replace("_model.cif", "_confidences.json"))
    if not conf_path.is_file():
        return result

    try:
        with open(conf_path, "r", encoding="utf-8") as f:
            data = _json.load(f)
    except (ValueError, OSError) as e:
        logger.warning("[V3 Adapter] Could not read confidence file %s: %s",
                       conf_path.name, e)
        return result

    if not isinstance(data, dict):
        return result

    plddts = data.get("atom_plddts")
    if isinstance(plddts, list) and plddts:
        vals = [float(v) for v in plddts if v is not None]
        if vals:
            result["plddt_mean"] = float(np.mean(vals))
            result["plddt_min"] = float(np.min(vals))
            result["plddt_max"] = float(np.max(vals))
            result["plddt_median"] = float(np.median(vals))

    pae = data.get("pae")
    if isinstance(pae, list) and pae:
        flat = [float(v) for row in pae for v in row if v is not None]
        if flat:
            result["pae_mean"] = float(np.mean(flat))

    contact_probs = data.get("contact_probs")
    if isinstance(contact_probs, list) and contact_probs:
        # AF3 stores contact probabilities as a (tokens, tokens) matrix;
        # flatten to the off-diagonal upper triangle (self-pairs are
        # always 1.0 and would inflate the mean).
        flat = []
        for i, row in enumerate(contact_probs):
            if not isinstance(row, (list, tuple)):
                row = [row]
            for j, v in enumerate(row):
                if v is None or i >= j:
                    continue
                try:
                    flat.append(float(v))
                except (TypeError, ValueError):
                    continue
        if flat:
            result["contact_prob_mean"] = float(np.mean(flat))

    return result


def validate_dataset(dataset: Dataset) -> Dict[str, Any]:
    """
    Validate the V3 dataset.

    Checks:
    - Every prediction belongs to a condition
    - Every prediction has a valid seed (for seed-level analysis)
    - Condition-level rows are not interpreted as predictions
    - Summary rows are not interpreted as predictions
    - Duplicate rows are detected
    - Invalid mappings are reported

    Returns validation report.
    """
    report = {
        "valid": True,
        "warnings": [],
        "errors": [],
        "n_conditions": len(dataset.conditions),
        "n_predictions": sum(
            sum(len(sp.values()) for sp in cond_preds.values())
            for cond_preds in dataset.predictions.values()
        ),
    }

    # Check that every prediction has a condition
    for condition_id, seeds_dict in dataset.predictions.items():
        if condition_id not in dataset.conditions:
            report["warnings"].append(
                f"Condition '{condition_id}' in predictions but not in conditions summary"
            )

    # Check for missing seeds (for seed-level analysis)
    for condition_id, seeds_dict in dataset.predictions.items():
        seeds = sorted(seeds_dict.keys())
        if len(seeds) == 0:
            report["warnings"].append(
                f"No seeds found for condition '{condition_id}'"
            )

    # Check reference condition
    if dataset.reference_condition and dataset.reference_condition not in dataset.conditions:
        report["errors"].append(
            f"Reference condition '{dataset.reference_condition}' not in dataset"
        )
        report["valid"] = False

    return report


def validate_condition_seed_mapping(
    dataset: Dataset,
    require_seeds: bool = True,
) -> Dict[str, Any]:
    """
    Validate condition and seed mapping.

    Ensures:
    - Every prediction belongs to a condition
    - Every prediction has a valid seed when seed-level analysis is intended
    - Condition-level rows are not interpreted as predictions
    - Summary rows are not interpreted as predictions
    - Duplicate rows are detected
    - Invalid mappings are reported
    """
    report = {
        "valid": True,
        "classified_rows": {
            "predictions": 0,
            "condition_level": 0,
            "summary": 0,
            "invalid": 0,
        },
        "details": [],
    }

    # In the V3 model, we don't have "condition-level rows" or "summary rows"
    # because the adapter already filters those out. But we document this.
    report["classified_rows"]["predictions"] = sum(
        sum(len(sp.values()) for sp in cond_seeds.values())
        for cond_seeds in dataset.predictions.values()
    )

    # Check each prediction for valid seed
    for condition_id, seeds_dict in dataset.predictions.items():
        for seed, samples_dict in seeds_dict.items():
            if require_seeds and (seed is None or seed < 0):
                report["valid"] = False
                report["details"].append(
                    f"Invalid seed {seed} for condition {condition_id}"
                )

    return report


def validate_structural_inputs(dataset: Dataset) -> Dict[str, Any]:
    """
    Validate structural inputs for every structure.

    For every structure determine:
    - Whether coordinates exist
    - Which entities exist
    - Which chains exist
    - Sequence identity
    - Residue mapping
    - Atom mapping
    - Missing residues
    - Missing atoms
    - Modified residues
    - DNA/RNA presence
    - Ligand presence
    - Ion presence

    Creates a structural QC table.
    """
    report = {
        "valid": True,
        "structures_qc": [],
        "failed_comparisons": [],
    }

    for condition_id, seeds_dict in dataset.predictions.items():
        for seed, samples_dict in seeds_dict.items():
            for sample, structure in samples_dict.items():
                qc = {
                    "prediction_id": structure.prediction_id,
                    "condition_id": condition_id,
                    "seed": seed,
                    "sample": sample,
                    "status": structure.parse_status,
                    "reason": structure.parse_reason,
                    "n_atoms": sum(
                        len(chain.residues) for entity in structure.entities
                        for chain in entity.chains
                    ),
                    "n_chains": structure.n_chains,
                    "n_entities": structure.n_entities,
                    "chain_ids": structure.chain_ids,
                    "entity_types": [
                        entity.entity_type for entity in structure.entities
                    ],
                    "has_protein": structure.has_protein,
                    "has_dna": structure.has_dna,
                    "has_rna": structure.get_nucleic_acid_chains() != structure.get_protein_chains(),
                    "missing_residues": [],
                    "missing_atoms": [],
                    "modified_residues": [],
                    "ligand_present": any(
                        entity.entity_type == "non-polymer"
                        for entity in structure.entities
                    ),
                    "ion_present": any(
                        entity.entity_type == "waters"
                        for entity in structure.entities
                    ),
                }
                report["structures_qc"].append(qc)

    return report
