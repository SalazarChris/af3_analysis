"""
Raw AF3 artifact reader for AF3 Confidence Analysis Pipeline.

Discovers and checksums raw confidence, summary, model, ranking, and mapping artifacts.
"""

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass

from af3_analysis.io.provenance import Provenance, Checksum, SourceFile


class RawAF3Error(Exception):
    """Error reading raw AF3 artifacts."""
    pass


@dataclass(frozen=True)
class ArtifactClassification:
    """Classification of an AF3 artifact."""
    artifact_type: str  # "prediction", "summary", "top_level", "best_model", "malformed", "unclassified"
    condition_id: Optional[str] = None
    seed: Optional[int] = None
    sample: Optional[int] = None
    model_rank: Optional[int] = None
    reference_type: Optional[str] = None


@dataclass(frozen=True)
class ArtifactRecord:
    """Record for an AF3 artifact."""
    source_path: Path
    checksum: str
    file_size: int
    artifact_type: str
    condition_id: Optional[str] = None
    seed: Optional[int] = None
    sample: Optional[int] = None
    model_rank: Optional[int] = None
    af3_version: Optional[str] = None


class RawAF3Reader:
    """
    Read and classify raw AF3 artifacts.
    
    Discovers artifacts, computes checksums, classifies them,
    and extracts canonical (condition_id, seed, sample) identity.
    """
    
    # File patterns for different artifact types
    # ORDER MATTERS: More specific patterns must come first
    # Use [^_]+ for condition to match up to first underscore after condition name
    SUMMARY_PATTERNS = [
        re.compile(r"^(?P<cond>.+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_summary_confidences\.json$"),
    ]
    
    PREDICTION_PATTERNS = [
        re.compile(r"^(?P<cond>.+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_confidences\.json$"),
    ]
    
    TOP_LEVEL_PATTERNS = [
        re.compile(r"^(?P<cond>.+)_seed-(?P<seed>\d+)_confidences\.json$"),
        re.compile(r"^(?P<cond>.+)_seed-(?P<seed>\d+)_best_model_confidences\.json$"),
    ]
    
    MODEL_PATTERNS = [
        re.compile(r"^(?P<cond>.+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model\.cif$"),
    ]
    
    RANKING_PATTERNS = [
        re.compile(r"^(?P<cond>.+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_ranking_scores\.csv$"),
    ]
    
    def __init__(self, raw_af3_root: Path):
        self.raw_af3_root = raw_af3_root
        self._provenance = Provenance()
        self._artifacts: List[ArtifactRecord] = []
        self._classification_counts: Dict[str, int] = {}
    
    def discover_artifacts(self) -> List[ArtifactRecord]:
        """
        Discover all AF3 artifacts in the raw root.
        
        Returns list of ArtifactRecord instances.
        """
        artifacts = []
        
        # Walk the directory tree
        for path in self.raw_af3_root.rglob("*"):
            if not path.is_file():
                continue
            
            # Skip hidden files and directories
            if any(part.startswith(".") for part in path.parts):
                continue
            
            artifact = self._process_file(path)
            if artifact:
                artifacts.append(artifact)
        
        self._artifacts = artifacts
        return artifacts
    
    def _process_file(self, path: Path) -> Optional[ArtifactRecord]:
        """Process a single file and create an ArtifactRecord."""
        try:
            checksum = Checksum.from_file(path).value
            file_size = path.stat().st_size
            
            # Classify the artifact
            classification = self._classify_file(path)
            
            return ArtifactRecord(
                source_path=path,
                checksum=checksum,
                file_size=file_size,
                artifact_type=classification.artifact_type,
                condition_id=classification.condition_id,
                seed=classification.seed,
                sample=classification.sample,
                model_rank=classification.model_rank,
            )
        except Exception as e:
            # Log error but don't fail completely
            print(f"Warning: Could not process {path}: {e}")
            return None
    
    def _classify_file(self, path: Path) -> ArtifactClassification:
        """Classify an AF3 file based on its name and content."""
        filename = path.name.lower()
        
        # Try prediction patterns first
        for pattern in self.PREDICTION_PATTERNS:
            match = pattern.match(filename)
            if match:
                return ArtifactClassification(
                    artifact_type="prediction",
                    condition_id=match.group("cond"),
                    seed=int(match.group("seed")),
                    sample=int(match.group("sample")),
                )
        
        # Try summary patterns
        for pattern in self.SUMMARY_PATTERNS:
            match = pattern.match(filename)
            if match:
                return ArtifactClassification(
                    artifact_type="summary",
                    condition_id=match.group("cond"),
                    seed=int(match.group("seed")),
                    sample=int(match.group("sample")),
                )
        
        # Try top-level patterns
        for pattern in self.TOP_LEVEL_PATTERNS:
            match = pattern.match(filename)
            if match:
                return ArtifactClassification(
                    artifact_type="top_level",
                    condition_id=match.group("cond"),
                    seed=int(match.group("seed")),
                    reference_type=self._determine_reference_type(path),
                )
        
        # Try model patterns
        for pattern in self.MODEL_PATTERNS:
            match = pattern.match(filename)
            if match:
                return ArtifactClassification(
                    artifact_type="best_model",
                    condition_id=match.group("cond"),
                    seed=int(match.group("seed")),
                    sample=int(match.group("sample")),
                )
        
        # Try ranking patterns
        for pattern in self.RANKING_PATTERNS:
            match = pattern.match(filename)
            if match:
                return ArtifactClassification(
                    artifact_type="ranking",
                    condition_id=match.group("cond"),
                    seed=int(match.group("seed")),
                    sample=int(match.group("sample")),
                )
        
        # Try to read JSON and extract metadata
        if filename.endswith("_confidences.json") or filename.endswith("_data.json"):
            try:
                with open(path, "r") as f:
                    data = json.load(f)
                
                # Try to extract seed/sample from data
                seed = data.get("seed")
                sample = data.get("sample")
                condition = data.get("condition")
                
                if seed is not None and sample is not None:
                    return ArtifactClassification(
                        artifact_type="prediction",
                        condition_id=str(condition) if condition else None,
                        seed=int(seed),
                        sample=int(sample),
                    )
            except (json.JSONDecodeError, IOError):
                pass
        
        # Unclassified
        return ArtifactClassification(
            artifact_type="unclassified",
            condition_id=None,
            seed=None,
            sample=None,
        )
    
    def _determine_reference_type(self, path: Path) -> Optional[str]:
        """Determine the type of reference artifact."""
        filename = path.name.lower()
        
        if "best_model" in filename:
            return "best_ranked"
        elif "top_level" in filename:
            return "top_level"
        
        return None
    
    def get_artifacts_by_type(self, artifact_type: str) -> List[ArtifactRecord]:
        """Get artifacts of a specific type."""
        return [a for a in self._artifacts if a.artifact_type == artifact_type]
    
    def get_artifacts_by_condition(self, condition_id: str) -> List[ArtifactRecord]:
        """Get artifacts for a specific condition."""
        return [a for a in self._artifacts if a.condition_id == condition_id]
    
    def get_predictions(self) -> List[ArtifactRecord]:
        """Get all prediction artifacts."""
        return self.get_artifacts_by_type("prediction")
    
    def get_summaries(self) -> List[ArtifactRecord]:
        """Get all summary artifacts."""
        return self.get_artifacts_by_type("summary")
    
    def get_top_levels(self) -> List[ArtifactRecord]:
        """Get all top-level reference artifacts."""
        return self.get_artifacts_by_type("top_level")
    
    def to_inventory(self) -> Dict[str, Any]:
        """Create an inventory of discovered artifacts."""
        inventory = {
            "raw_af3_root": str(self.raw_af3_root),
            "total_artifacts": len(self._artifacts),
            "by_type": {},
        }
        
        for artifact in self._artifacts:
            artifact_type = artifact.artifact_type
            if artifact_type not in inventory["by_type"]:
                inventory["by_type"][artifact_type] = 0
            inventory["by_type"][artifact_type] += 1
        
        return inventory


def discover_raw_af3_artifacts(raw_af3_root: Path) -> RawAF3Reader:
    """
    Convenience function to discover all AF3 artifacts.
    
    Args:
        raw_af3_root: Path to raw AF3 output directory
    
    Returns:
        RawAF3Reader instance with artifacts discovered
    """
    reader = RawAF3Reader(raw_af3_root)
    reader.discover_artifacts()
    return reader


__all__ = [
    "RawAF3Error",
    "ArtifactClassification",
    "ArtifactRecord",
    "RawAF3Reader",
    "discover_raw_af3_artifacts",
]
