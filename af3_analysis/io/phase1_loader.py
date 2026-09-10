"""
Phase 1 CSV loader for AF3 Confidence Analysis Pipeline.

Loads the four option-6 CSVs as an auditable legacy intake artifact.
"""

import hashlib
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass

from af3_analysis.io.provenance import Provenance, Checksum


class LegacyCSVError(Exception):
    """Error loading legacy CSV files."""
    pass


@dataclass(frozen=True)
class LegacyInput:
    """Input from a single legacy CSV file."""
    source_path: Path
    checksum: str
    source_row_id: int
    dynamic_chain_columns: List[str]
    raw_data: Dict[str, Any]


class Phase1Loader:
    """
    Load the current four-CSV option-6 output as an auditable legacy intake artifact.
    
    Preserves source row IDs, source paths, checksums, dynamic chain columns,
    and legacy condition identifiers.
    """
    
    EXPECTED_CSV_FILES = [
        "condition_manifest.csv",
        "condition_registry.csv",
        "metrics_replicates.csv",
        "metrics_conditions.csv",
    ]
    
    def __init__(self, legacy_csv_root: Path):
        self.legacy_csv_root = legacy_csv_root
        self._provenance = Provenance()
        self._inputs: Dict[str, LegacyInput] = {}
        self._raw_dataframes: Dict[str, pd.DataFrame] = {}
    
    def load(self) -> Dict[str, pd.DataFrame]:
        """
        Load all four expected CSV files.
        
        Returns dict mapping CSV name to DataFrame.
        """
        dataframes = {}
        
        for csv_name in self.EXPECTED_CSV_FILES:
            csv_path = self.legacy_csv_root / csv_name
            if not csv_path.exists():
                raise LegacyCSVError(f"Missing expected CSV file: {csv_name}")
            
            df = self._load_csv(csv_path)
            dataframes[csv_name] = df
            
            # Store provenance
            checksum = self._compute_checksum(csv_path)
            self._inputs[csv_name] = LegacyInput(
                source_path=csv_path,
                checksum=checksum,
                source_row_id=0,  # Will be set per row
                dynamic_chain_columns=self._extract_dynamic_chain_columns(df),
                raw_data={},
            )
        
        return dataframes
    
    def _load_csv(self, path: Path) -> pd.DataFrame:
        """Load a CSV file with source tracking."""
        df = pd.read_csv(path)
        
        # Add source path column for provenance
        df["_source_path"] = str(path)
        df["_source_checksum"] = self._compute_checksum(path)
        
        # Add source row ID column
        df["_source_row_id"] = df.index
        
        return df
    
    def _compute_checksum(self, path: Path) -> str:
        """Compute SHA256 checksum of a file."""
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    
    def _extract_dynamic_chain_columns(self, df: pd.DataFrame) -> List[str]:
        """Extract chain-related columns from the DataFrame."""
        chain_cols = []
        for col in df.columns:
            if col.startswith("chain_") and col.endswith(("residues", "plddt")):
                chain_cols.append(col)
        return chain_cols
    
    def get_provenance(self, csv_name: str) -> Optional[LegacyInput]:
        """Get provenance information for a CSV."""
        return self._inputs.get(csv_name)
    
    def get_raw_dataframe(self, csv_name: str) -> Optional[pd.DataFrame]:
        """Get the raw DataFrame for a CSV."""
        return self._raw_dataframes.get(csv_name)
    
    def to_inventory(self) -> Dict[str, Any]:
        """Create an immutable input inventory."""
        inventory = {
            "legacy_csv_root": str(self.legacy_csv_root),
            "csv_files": {},
        }
        
        for csv_name, legacy_input in self._inputs.items():
            inventory["csv_files"][csv_name] = {
                "source_path": str(legacy_input.source_path),
                "checksum": legacy_input.checksum,
                "dynamic_chain_columns": legacy_input.dynamic_chain_columns,
                "row_count": len(self._raw_dataframes.get(csv_name, pd.DataFrame())),
            }
        
        return inventory


def load_phase1_legacy(legacy_csv_root: Path) -> Dict[str, pd.DataFrame]:
    """
    Convenience function to load all Phase 1 legacy CSVs.
    
    Args:
        legacy_csv_root: Path to directory containing the four expected CSV files
    
    Returns:
        Dict mapping CSV names to DataFrames
    """
    loader = Phase1Loader(legacy_csv_root)
    return loader.load()


__all__ = [
    "LegacyCSVError",
    "LegacyInput",
    "Phase1Loader",
    "load_phase1_legacy",
]
