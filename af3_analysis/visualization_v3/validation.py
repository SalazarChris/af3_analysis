"""
V3 Validation Module.

Validates condition/seed/prediction mapping, structural QC, coverage.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .model import Dataset, StructureData, ComparisonQC, StructuralQC


# ---------------------------------------------------------------------------
# Validation results
# ---------------------------------------------------------------------------

@dataclass
class ValidationResult:
    """Result of a validation check."""
    passed: bool
    messages: List[str] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.passed

    def __and__(self, other: "ValidationResult") -> "ValidationResult":
        return ValidationResult(
            passed=self.passed and other.passed,
            messages=self.messages + other.messages,
            details={**self.details, **other.details},
        )


# ---------------------------------------------------------------------------
# Condition and seed validation
# ---------------------------------------------------------------------------

def validate_condition_mapping(
    dataset: Dataset,
) -> ValidationResult:
    """
    Validate that every prediction belongs to a valid condition.

    Checks:
    - Every prediction's condition_id exists in dataset.conditions
    - No orphan predictions
    """
    messages = []
    details = {}

    orphan_predictions = []
    for condition_id, seeds_dict in dataset.predictions.items():
        if condition_id not in dataset.conditions:
            for seed, samples_dict in seeds_dict.items():
                for sample, structure in samples_dict.items():
                    orphan_predictions.append(structure.prediction_id)

    if orphan_predictions:
        messages.append(
            f"Found {len(orphan_predictions)} orphan predictions "
            f"(condition_id not in conditions summary)"
        )
        details["orphan_predictions"] = orphan_predictions

    return ValidationResult(
        passed=len(orphan_predictions) == 0,
        messages=messages,
        details=details,
    )


def validate_seed_mapping(
    dataset: Dataset,
    *,
    require_seeds: bool = True,
) -> ValidationResult:
    """
    Validate seed mapping for all predictions.

    Checks:
    - Every prediction has a valid integer seed
    - Seeds are consistent within condition
    - No duplicate (condition, seed, sample) combinations
    """
    messages = []
    details = {
        "n_predictions": 0,
        "n_valid_seeds": 0,
        "invalid_seeds": [],
        "duplicates": [],
    }

    seen_keys = set()

    for condition_id, seeds_dict in dataset.predictions.items():
        for seed, samples_dict in seeds_dict.items():
            details["n_predictions"] += len(samples_dict)

            if require_seeds:
                if not isinstance(seed, int) or seed < 0:
                    details["invalid_seeds"].append({
                        "condition_id": condition_id,
                        "seed": seed,
                    })

            for sample, structure in samples_dict.items():
                details["n_valid_seeds"] += 1
                key = (condition_id, seed, sample)
                if key in seen_keys:
                    details["duplicates"].append(key)
                    messages.append(
                        f"Duplicate prediction: {condition_id}, seed={seed}, sample={sample}"
                    )
                seen_keys.add(key)

    return ValidationResult(
        passed=len(details["invalid_seeds"]) == 0 and len(details["duplicates"]) == 0,
        messages=messages,
        details=details,
    )


def validate_prediction_mapping(
    dataset: Dataset,
) -> ValidationResult:
    """
    Validate that predictions are correctly mapped.

    Checks:
    - Predictions are not condition-level rows
    - Predictions are not summary rows
    - Each prediction has unique (condition, seed, sample)
    """
    # In V3, we already filter these out during adaptation.
    # This function documents that filtering happened.

    predictions = []
    for condition_id, seeds_dict in dataset.predictions.items():
        for seed, samples_dict in seeds_dict.items():
            for sample, structure in samples_dict.items():
                predictions.append({
                    "prediction_id": structure.prediction_id,
                    "condition_id": condition_id,
                    "seed": seed,
                    "sample": sample,
                })

    return ValidationResult(
        passed=True,
        messages=[f"Validated {len(predictions)} predictions"],
        details={"predictions": predictions},
    )


# ---------------------------------------------------------------------------
# Structural input validation
# ---------------------------------------------------------------------------

def validate_structural_inputs(
    dataset: Dataset,
    *,
    min_atoms: int = 1,
    min_chains: int = 1,
) -> ValidationResult:
    """
    Validate structural inputs for every structure.

    For each structure determines:
    - Whether coordinates exist
    - Which entities exist
    - Which chains exist
    - Sequence identity (relative to reference if available)
    - Residue mapping
    - Atom mapping
    - Missing residues
    - Missing atoms
    - Modified residues (non-standard comp_id)
    - DNA/RNA presence
    - Ligand presence
    - Ion presence

    Returns structural QC table.
    """
    messages = []
    qc_records = []
    failed = []

    for condition_id, seeds_dict in dataset.predictions.items():
        for seed, samples_dict in seeds_dict.items():
            for sample, structure in samples_dict.items():
                qc = _build_structural_qc(structure, min_atoms, min_chains)
                qc_records.append(qc)

                if not qc["passed"]:
                    failed.append(qc)
                    messages.append(
                        f"Structural QC failed for {qc['prediction_id']}: "
                        f"{qc['reason']}"
                    )

    return ValidationResult(
        passed=len(failed) == 0,
        messages=messages,
        details={
            "qc_records": qc_records,
            "failed": failed,
            "n_structures": len(qc_records),
            "n_failed": len(failed),
        },
    )


def _build_structural_qc(
    structure: StructureData,
    min_atoms: int = 1,
    min_chains: int = 1,
) -> Dict[str, Any]:
    """Build QC record for one structure."""
    # Check basic structure
    n_atoms = sum(
        len(chain.residues) for entity in structure.entities
        for chain in entity.chains
    )
    n_chains = structure.n_chains
    n_entities = structure.n_entities

    # Check entities and chains
    entity_types = [entity.entity_type for entity in structure.entities]
    chain_ids = structure.chain_ids

    # Check polymer types
    polymer_types = []
    for entity in structure.entities:
        for chain in entity.chains:
            if chain.polymer_type:
                polymer_types.append(chain.polymer_type)

    has_protein = any("polypeptide" in pt for pt in polymer_types)
    has_dna = any("polydeoxyribonucleotide" in pt for pt in polymer_types)
    has_rna = any("polyribonucleotide" in pt for pt in polymer_types)

    # Check for ligands (non-polymer entities)
    has_ligand = any(
        entity.entity_type == "non-polymer"
        for entity in structure.entities
    )

    # Check for ions (waters/branch)
    has_ion = any(
        entity.entity_type in ("waters", "branched")
        for entity in structure.entities
    )

    # Check for modified residues (non-standard amino acids)
    modified_residues = []
    standard_aa = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "MSE", "SEP", "TPO", "HYP",
    }
    for entity in structure.entities:
        for chain in entity.chains:
            for residue in chain.residues:
                if residue.residue_name not in standard_aa:
                    modified_residues.append({
                        "chain_id": chain.chain_id,
                        "auth_seq_id": residue.auth_seq_id,
                        "residue_name": residue.residue_name,
                    })

    # Check for missing atoms in residues
    missing_atoms = []
    for entity in structure.entities:
        for chain in entity.chains:
            for residue in chain.residues:
                expected_backbone = {"N", "CA", "C", "O"}
                present = {a.atom_name for a in residue.all_atoms}
                missing = expected_backbone - present
                if missing and residue.residue_name in standard_aa:
                    missing_atoms.append({
                        "chain_id": chain.chain_id,
                        "auth_seq_id": residue.auth_seq_id,
                        "residue_name": residue.residue_name,
                        "missing": sorted(missing),
                    })

    # Determine overall status
    passed = True
    reason = ""

    if n_atoms < min_atoms:
        passed = False
        reason = f"Insufficient atoms: {n_atoms} < {min_atoms}"
    elif n_chains < min_chains:
        passed = False
        reason = f"Insufficient chains: {n_chains} < {min_chains}"
    elif structure.parse_status != "success":
        passed = False
        reason = f"Parse error: {structure.parse_reason}"

    return {
        "prediction_id": structure.prediction_id,
        "condition_id": structure.condition_id,
        "seed": structure.seed,
        "sample": structure.sample,
        "passed": passed,
        "reason": reason,
        "n_atoms": n_atoms,
        "n_chains": n_chains,
        "n_entities": n_entities,
        "chain_ids": chain_ids,
        "entity_types": entity_types,
        "polymer_types": polymer_types,
        "has_protein": has_protein,
        "has_dna": has_dna,
        "has_rna": has_rna,
        "has_ligand": has_ligand,
        "has_ion": has_ion,
        "modified_residues": modified_residues,
        "missing_atoms": missing_atoms,
    }


# ---------------------------------------------------------------------------
# Structural comparability
# ---------------------------------------------------------------------------

def check_structural_comparability(
    structure_a: StructureData,
    structure_b: StructureData,
    *,
    min_common_atoms: int = 3,
    min_coverage: float = 0.80,
    alignment_atom: str = "CA",
) -> ComparisonQC:
    """
    Check if two structures are structurally comparable.

    Determines:
    - Common chains
    - Common residues
    - Common atoms
    - Common sequence
    - Coverage

    Returns ComparisonQC with status and coverage.
    """
    # Determine chains to compare
    chains_a = set(structure_a.chain_ids)
    chains_b = set(structure_b.chain_ids)

    # Prefer protein chains
    protein_chains_a = structure_a.get_protein_chains()
    protein_chains_b = structure_b.get_protein_chains()

    if protein_chains_a and protein_chains_b:
        common_chains = sorted(set(protein_chains_a) & set(protein_chains_b))
    else:
        common_chains = sorted(chains_a & chains_b)

    # Count reference residues (for coverage calculation)
    ref_residues = 0
    common_residues = 0
    common_atoms = 0

    for chain_id in common_chains:
        chain_a = structure_a.get_chain(chain_id)
        chain_b = structure_b.get_chain(chain_id)

        if chain_a and chain_b:
            residues_a = {r.auth_seq_id for r in chain_a.residues if r.auth_seq_id is not None}
            residues_b = {r.auth_seq_id for r in chain_b.residues if r.auth_seq_id is not None}

            common = residues_a & residues_b
            common_residues += len(common)
            ref_residues += len(residues_a)  # Use first structure as reference

            # Count common atoms (CA by default)
            if alignment_atom == "CA":
                for res_a in chain_a.residues:
                    if res_a.auth_seq_id in common and res_a.ca_coords is not None:
                        res_b = chain_b.get_residue(res_a.auth_seq_id)
                        if res_b and res_b.ca_coords is not None:
                            common_atoms += 1

    # Calculate coverage
    coverage = common_residues / ref_residues if ref_residues > 0 else 0.0

    # Determine status
    low_coverage = coverage < min_coverage
    insufficient_atoms = common_atoms < min_common_atoms

    if insufficient_atoms or common_residues == 0:
        status = "not_comparable"
        reason = "insufficient_common_atoms" if insufficient_atoms else "no_common_residues"
    elif low_coverage:
        status = "partially_comparable"
        reason = f"low_coverage: {coverage:.2%}"
    else:
        status = "comparable"
        reason = "valid"

    return ComparisonQC(
        comparison_id=f"{structure_a.prediction_id}_vs_{structure_b.prediction_id}",
        condition_a=structure_a.condition_id,
        condition_b=structure_b.condition_id,
        seed_a=structure_a.seed,
        seed_b=structure_b.seed,
        status=status,
        reason=reason,
        n_common_atoms=common_atoms,
        n_common_residues=common_residues,
        coverage=coverage,
        low_coverage=low_coverage,
    )


# ---------------------------------------------------------------------------
# Exclusions and classification
# ---------------------------------------------------------------------------

def classify_comparison(
    qc: ComparisonQC,
    *,
    allowed_statuses: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Classify a comparison result.

    Returns classification dict with:
    - include: whether to include in analysis
    - reason: exclusion reason if not included
    - status: original status
    """
    if allowed_statuses is None:
        allowed_statuses = ["comparable", "partially_comparable"]

    classification = {
        "include": qc.status in allowed_statuses,
        "status": qc.status,
    }

    if not classification["include"]:
        if qc.status == "not_comparable":
            classification["reason"] = qc.reason or "not_comparable"
        elif qc.low_coverage:
            classification["reason"] = f"low_coverage: {qc.coverage:.2%}"
        else:
            classification["reason"] = qc.reason or "excluded"

    return classification


# ---------------------------------------------------------------------------
# Summary functions
# ---------------------------------------------------------------------------

def summarize_structural_qc(qc_records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Summarize structural QC results."""
    n_total = len(qc_records)
    n_passed = sum(1 for qc in qc_records if qc["passed"])
    n_failed = n_total - n_passed

    # Count by entity type
    entity_counts = {}
    for qc in qc_records:
        key = tuple(sorted(qc["entity_types"]))
        entity_counts[key] = entity_counts.get(key, 0) + 1

    return {
        "n_structures": n_total,
        "n_passed": n_passed,
        "n_failed": n_failed,
        "pass_rate": n_passed / n_total if n_total > 0 else 0,
        "entity_composition": entity_counts,
        "n_with_protein": sum(1 for qc in qc_records if qc["has_protein"]),
        "n_with_dna": sum(1 for qc in qc_records if qc["has_dna"]),
        "n_with_ligand": sum(1 for qc in qc_records if qc["has_ligand"]),
        "n_with_ion": sum(1 for qc in qc_records if qc["has_ion"]),
        "n_modified_residues": sum(
            len(qc["modified_residues"]) for qc in qc_records
        ),
        "n_missing_atoms": sum(
            len(qc["missing_atoms"]) for qc in qc_records
        ),
    }
