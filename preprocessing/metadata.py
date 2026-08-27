"""
Metadata loading for AF3 Confidence Analysis Pipeline.

Loads and validates condition metadata from legacy sources.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
import pandas as pd


class MetadataLoadError(Exception):
    """Error loading condition metadata."""
    pass


def load_condition_metadata(
    legacy_csv_root: Path,
) -> Dict[str, Dict[str, Any]]:
    """
    Load condition metadata from legacy CSVs.
    
    Only loads reviewed condition metadata from condition_manifest.csv.
    Does NOT infer factors from condition names.
    
    Args:
        legacy_csv_root: Path to directory containing legacy CSVs
    
    Returns:
        Dict mapping condition_id to metadata dict
    """
    manifest_path = legacy_csv_root / "condition_manifest.csv"
    registry_path = legacy_csv_root / "condition_registry.csv"
    
    if not manifest_path.exists():
        raise MetadataLoadError(
            f"condition_manifest.csv not found at {manifest_path}"
        )
    
    # Load manifest
    manifest_df = pd.read_csv(manifest_path)
    
    # Load registry for additional metadata
    condition_metadata = {}
    if registry_path.exists():
        registry_df = pd.read_csv(registry_path)
        
        # Merge or add registry data
        for _, row in manifest_df.iterrows():
            cond_id = row["condition_id"]
            metadata = row.to_dict()
            
            # Try to find matching registry entry
            reg_row = registry_df[
                registry_df["condition_id"] == cond_id
            ]
            if len(reg_row) == 1:
                reg_dict = reg_row.iloc[0].to_dict()
                metadata["replicate_list"] = reg_dict.get("replicate_list", [])
                metadata["reported_replicate_count"] = reg_dict.get(
                    "reported_replicate_count", 0
                )
            
            condition_metadata[cond_id] = metadata
    else:
        # Use manifest only
        for _, row in manifest_df.iterrows():
            condition_metadata[row["condition_id"]] = row.to_dict()
    
    return condition_metadata


def validate_condition_metadata(
    metadata: Dict[str, Dict[str, Any]],
) -> List[str]:
    """
    Validate condition metadata and return any issues.
    
    Args:
        metadata: Dict mapping condition_id to metadata
    
    Returns:
        List of issue descriptions (empty if valid)
    """
    issues = []
    
    if not metadata:
        issues.append("No conditions found in metadata")
        return issues
    
    # Check for duplicate condition names
    condition_names = [
        m.get("condition_name", m.get("condition_id", "unknown"))
        for m in metadata.values()
    ]
    if len(condition_names) != len(set(condition_names)):
        issues.append("Duplicate condition names found")
    
    # Check for missing condition IDs
    for cond_id, meta in metadata.items():
        if not cond_id:
            issues.append("Condition with empty ID found")
        if "condition_name" not in meta:
            issues.append(f"Condition {cond_id} missing condition_name")
    
    return issues


__all__ = [
    "MetadataLoadError",
    "load_condition_metadata",
    "validate_condition_metadata",
]
