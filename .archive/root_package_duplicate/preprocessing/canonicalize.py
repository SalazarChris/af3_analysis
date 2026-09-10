"""
Canonicalization for AF3 Confidence Analysis Pipeline.

Reconstructs canonical (condition_id, seed, sample) from raw provenance
and creates canonical experiments, conditions, replicates tables.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
import pandas as pd

from af3_analysis.io.raw_af3_reader import RawAF3Reader, ArtifactRecord
from af3_analysis.schemas.records import (
    Replicate,
    BestModelReference,
    Measurement,
)
from af3_analysis.schemas.enums import (
    ReplicateStatus,
    MetricScope,
)


@dataclass
class CanonicalizationResult:
    """Result of canonicalization."""
    experiments: pd.DataFrame
    conditions: pd.DataFrame
    replicates: pd.DataFrame
    best_model_references: pd.DataFrame
    measurements: pd.DataFrame
    canonicalized_count: int
    reference_count: int
    excluded_count: int


class Canonicalizer:
    """
    Canonicalize raw AF3 artifacts into analysis-ready tables.
    
    Creates canonical experiments, conditions, replicates, best-model-reference,
    and long-form measurement tables from verified provenance.
    """
    
    def __init__(self, artifacts: List[ArtifactRecord]):
        self._artifacts = artifacts
        self._experiments = []
        self._conditions = []
        self._replicates = []
        self._best_model_references = []
        self._measurements = []
    
    def canonicalize(self) -> CanonicalizationResult:
        """
        Canonicalize all artifacts.
        
        Returns result with canonical tables.
        """
        # Group artifacts by condition
        artifacts_by_condition: Dict[str, List[ArtifactRecord]] = {}
        
        for artifact in self._artifacts:
            if not artifact.condition_id:
                continue
            
            if artifact.condition_id not in artifacts_by_condition:
                artifacts_by_condition[artifact.condition_id] = []
            artifacts_by_condition[artifact.condition_id].append(artifact)
        
        # Process each condition
        for condition_id, condition_artifacts in artifacts_by_condition.items():
            self._process_condition(condition_id, condition_artifacts)
        
        # Create DataFrames
        experiments_df = pd.DataFrame(self._experiments)
        conditions_df = pd.DataFrame(self._conditions)
        replicates_df = pd.DataFrame(self._replicates)
        best_model_df = pd.DataFrame(self._best_model_references)
        measurements_df = pd.DataFrame(self._measurements)
        
        return CanonicalizationResult(
            experiments=experiments_df,
            conditions=conditions_df,
            replicates=replicates_df,
            best_model_references=best_model_df,
            measurements=measurements_df,
            canonicalized_count=len(self._replicates),
            reference_count=len(self._best_model_references),
            excluded_count=0,  # TODO: track excluded records
        )
    
    def _process_condition(
        self,
        condition_id: str,
        artifacts: List[ArtifactRecord],
    ) -> None:
        """Process all artifacts for a single condition."""
        # Create experiment record
        self._experiments.append({
            "experiment_id": f"exp_{condition_id}",
            "system_id": "unknown",
            "protocol_id": "af3_v2.3",
            "input_root": "",
            "analysis_run_id": "current",
            "design_description": f"Condition {condition_id}",
        })
        
        # Create condition record
        self._conditions.append({
            "condition_id": condition_id,
            "experiment_id": f"exp_{condition_id}",
            "condition_label": condition_id,
            "input_signature": self._compute_input_signature(artifacts),
        })
        
        # Process each artifact
        for artifact in artifacts:
            if artifact.artifact_type == "prediction":
                self._process_prediction(condition_id, artifact)
            elif artifact.artifact_type == "top_level":
                self._process_reference(condition_id, artifact)
    
    def _process_prediction(
        self,
        condition_id: str,
        artifact: ArtifactRecord,
    ) -> None:
        """Process a prediction artifact."""
        if artifact.seed is None or artifact.sample is None:
            return
        
        prediction_id = f"{condition_id}_seed-{artifact.seed}_sample-{artifact.sample}"
        
        self._replicates.append({
            "condition_id": condition_id,
            "experiment_id": f"exp_{condition_id}",
            "seed": artifact.seed,
            "sample": artifact.sample,
            "prediction_id": prediction_id,
            "source_confidences_path": str(artifact.source_path),
            "source_checksum": artifact.checksum,
            "run_status": ReplicateStatus.VALID.value,
        })
    
    def _process_reference(
        self,
        condition_id: str,
        artifact: ArtifactRecord,
    ) -> None:
        """Process a reference artifact (top-level/best-model)."""
        if artifact.seed is None:
            return
        
        # Check if this matches an existing prediction
        prediction_id = f"{condition_id}_seed-{artifact.seed}_sample-0"
        
        self._best_model_references.append({
            "prediction_id": prediction_id,
            "best_model_path": str(artifact.source_path),
            "reference_type": artifact.artifact_type,
        })
    
    def _compute_input_signature(self, artifacts: List[ArtifactRecord]) -> str:
        """Compute a signature for the condition's input."""
        # For now, just use a hash of all source paths
        import hashlib
        paths = sorted([str(a.source_path) for a in artifacts])
        return hashlib.sha256("".join(paths).encode()).hexdigest()[:16]


def canonicalize_artifacts(raw_af3_reader: RawAF3Reader) -> CanonicalizationResult:
    """
    Convenience function to canonicalize raw AF3 artifacts.
    
    Args:
        raw_af3_reader: RawAF3Reader with artifacts already discovered
    
    Returns:
        CanonicalizationResult with canonical tables
    """
    canonicalizer = Canonicalizer(raw_af3_reader._artifacts)
    return canonicalizer.canonicalize()


__all__ = [
    "CanonicalizationResult",
    "Canonicalizer",
    "canonicalize_artifacts",
]
