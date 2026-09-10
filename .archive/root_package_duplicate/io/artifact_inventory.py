"""
Artifact inventory for AF3 Confidence Analysis Pipeline.

Creates immutable input inventory from Phase 1 legacy CSVs.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
import pandas as pd


@dataclass(frozen=True)
class CSVInventory:
    """Inventory entry for a single CSV file."""
    name: str
    source_path: Path
    checksum: str
    row_count: int
    dynamic_chain_columns: List[str]
    source_row_ids: List[int]


@dataclass(frozen=True)
class LegacyInputInventory:
    """Immutable inventory of all Phase 1 legacy inputs."""
    legacy_csv_root: Path
    csv_files: Dict[str, CSVInventory]
    total_replicates: int
    total_conditions: int
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert inventory to dictionary for JSON serialization."""
        return {
            "legacy_csv_root": str(self.legacy_csv_root),
            "csv_files": {
                name: {
                    "name": inv.name,
                    "source_path": str(inv.source_path),
                    "checksum": inv.checksum,
                    "row_count": inv.row_count,
                    "dynamic_chain_columns": inv.dynamic_chain_columns,
                    "source_row_ids": inv.source_row_ids,
                }
                for name, inv in self.csv_files.items()
            },
            "total_replicates": self.total_replicates,
            "total_conditions": self.total_conditions,
        }


class ArtifactInventory:
    """
    Create immutable input inventory from Phase 1 legacy CSVs.
    
    Does NOT perform canonicalization - this is just an inventory.
    """
    
    def __init__(self, legacy_csv_root: Path, dataframes: Dict[str, pd.DataFrame]):
        self.legacy_csv_root = legacy_csv_root
        self._dataframes = dataframes
        self._inventory = self._build_inventory()
    
    def _build_inventory(self) -> LegacyInputInventory:
        """Build the immutable inventory."""
        csv_files = {}
        
        for csv_name, df in self._dataframes.items():
            # Compute checksum
            csv_path = self.legacy_csv_root / csv_name
            checksum = self._compute_checksum(csv_path)
            
            # Extract dynamic chain columns
            dynamic_chain_cols = self._extract_dynamic_chain_columns(df)
            
            # Create inventory entry
            csv_files[csv_name] = CSVInventory(
                name=csv_name,
                source_path=csv_path,
                checksum=checksum,
                row_count=len(df),
                dynamic_chain_columns=dynamic_chain_cols,
                source_row_ids=df["_source_row_id"].tolist(),
            )
        
        # Calculate totals
        total_replicates = self._count_replicates()
        total_conditions = self._count_conditions()
        
        return LegacyInputInventory(
            legacy_csv_root=self.legacy_csv_root,
            csv_files=csv_files,
            total_replicates=total_replicates,
            total_conditions=total_conditions,
        )
    
    def _compute_checksum(self, path: Path) -> str:
        """Compute SHA256 checksum of a file."""
        import hashlib
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    
    def _extract_dynamic_chain_columns(self, df: pd.DataFrame) -> List[str]:
        """Extract chain-related columns from the DataFrame."""
        chain_cols = []
        for col in df.columns:
            if col.startswith("chain_") and col.endswith(("residues", "plddt")):
                chain_cols.append(col)
        return chain_cols
    
    def _count_replicates(self) -> int:
        """Count total replicates from metrics_replicates.csv."""
        if "metrics_replicates.csv" in self._dataframes:
            return len(self._dataframes["metrics_replicates.csv"])
        return 0
    
    def _count_conditions(self) -> int:
        """Count unique conditions from condition_manifest.csv."""
        if "condition_manifest.csv" in self._dataframes:
            return len(self._dataframes["condition_manifest.csv"])
        return 0
    
    @property
    def inventory(self) -> LegacyInputInventory:
        """Get the immutable inventory."""
        return self._inventory
    
    def get_dataframe(self, csv_name: str) -> Optional[pd.DataFrame]:
        """Get the DataFrame for a CSV."""
        return self._dataframes.get(csv_name)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return self._inventory.to_dict()


def create_artifact_inventory(
    legacy_csv_root: Path,
    dataframes: Dict[str, pd.DataFrame],
) -> ArtifactInventory:
    """
    Convenience function to create artifact inventory.
    
    Args:
        legacy_csv_root: Path to legacy CSV directory
        dataframes: Dict mapping CSV names to DataFrames
    
    Returns:
        ArtifactInventory instance
    """
    return ArtifactInventory(legacy_csv_root, dataframes)


__all__ = [
    "CSVInventory",
    "LegacyInputInventory",
    "ArtifactInventory",
    "create_artifact_inventory",
]
