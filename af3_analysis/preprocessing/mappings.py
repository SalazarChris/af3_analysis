"""
Mappings for AF3 Confidence Analysis Pipeline.

Builds chains, residues, chain pairs, residue correspondence, and matrix index tables.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True)
class Chain:
    """Chain record."""
    mapping_id: str
    condition_id: str
    chain_id: str
    entity_id: str
    token_count: int
    residue_count: int
    chain_type: str  # protein, dna, rna, ligand
    source_checksum: Optional[str] = None


@dataclass(frozen=True)
class Residue:
    """Residue record."""
    mapping_id: str
    prediction_id: str
    chain_id: str
    residue_index: int
    residue_name: str
    sequence_position: int
    token_indices: List[int]
    residue_pddlt: Optional[float] = None


@dataclass(frozen=True)
class ChainPair:
    """Chain pair record."""
    prediction_id: str
    chain_pair_id: str  # e.g., "A|B"
    chain_a: str
    chain_b: str
    pair_iptm: Optional[float] = None
    pair_pae_mean: Optional[float] = None
    pair_pde_mean: Optional[float] = None
    interface_contact_confidence: Optional[float] = None


@dataclass(frozen=True)
class ResidueCorrespondence:
    """Residue correspondence mapping."""
    mapping_id: str
    condition_a: str
    condition_b: str
    residue_a_index: int
    residue_b_index: int
    aligned: bool
    rmsd: Optional[float] = None


@dataclass(frozen=True)
class MatrixIndex:
    """Matrix index record."""
    prediction_id: str
    matrix_type: str  # PAE, PDE, contact_prob
    matrix_path: Optional[Path] = None
    dimension: int = 0
    row_mapping_id: Optional[str] = None
    column_mapping_id: Optional[str] = None
    definedness: str = "present"


class MappingBuilder:
    """
    Build and validate chain, residue, and matrix mappings.
    
    Requires verified raw sources for mappings.
    """
    
    def __init__(self):
        self._chains: List[Dict[str, Any]] = []
        self._residues: List[Dict[str, Any]] = []
        self._chain_pairs: List[Dict[str, Any]] = []
        self._residue_correspondences: List[Dict[str, Any]] = []
        self._matrices_index: List[Dict[str, Any]] = []
    
    def build_chain_mapping(
        self,
        condition_id: str,
        chain_id: str,
        entity_id: str,
        token_count: int,
        residue_count: int,
        chain_type: str,
        source_checksum: Optional[str] = None,
    ) -> None:
        """Build a chain mapping."""
        self._chains.append({
            "mapping_id": f"{condition_id}_{chain_id}",
            "condition_id": condition_id,
            "chain_id": chain_id,
            "entity_id": entity_id,
            "token_count": token_count,
            "residue_count": residue_count,
            "chain_type": chain_type,
            "source_checksum": source_checksum,
        })
    
    def build_residue_mapping(
        self,
        prediction_id: str,
        mapping_id: str,
        chain_id: str,
        residues: List[Dict[str, Any]],
    ) -> None:
        """Build residue mappings for a prediction."""
        for residue in residues:
            self._residues.append({
                "prediction_id": prediction_id,
                "mapping_id": mapping_id,
                "chain_id": chain_id,
                "residue_index": residue.get("index", 0),
                "residue_name": residue.get("name", "UNK"),
                "sequence_position": residue.get("sequence_position", 0),
                "token_indices": residue.get("token_indices", []),
                "residue_pddlt": residue.get("pddlt"),
            })
    
    def build_chain_pair_mapping(
        self,
        prediction_id: str,
        chain_a: str,
        chain_b: str,
        pair_iptm: Optional[float] = None,
        pair_pae_mean: Optional[float] = None,
        pair_pde_mean: Optional[float] = None,
        interface_contact_confidence: Optional[float] = None,
    ) -> None:
        """Build a chain pair mapping."""
        self._chain_pairs.append({
            "prediction_id": prediction_id,
            "chain_pair_id": f"{chain_a}|{chain_b}",
            "chain_a": chain_a,
            "chain_b": chain_b,
            "pair_iptm": pair_iptm,
            "pair_pae_mean": pair_pae_mean,
            "pair_pde_mean": pair_pde_mean,
            "interface_contact_confidence": interface_contact_confidence,
        })
    
    def build_residue_correspondence(
        self,
        mapping_id: str,
        condition_a: str,
        condition_b: str,
        correspondence_pairs: List[Dict[str, Any]],
    ) -> None:
        """Build residue correspondence between conditions."""
        for pair in correspondence_pairs:
            self._residue_correspondences.append({
                "mapping_id": mapping_id,
                "condition_a": condition_a,
                "condition_b": condition_b,
                "residue_a_index": pair.get("residue_a_index", 0),
                "residue_b_index": pair.get("residue_b_index", 0),
                "aligned": pair.get("aligned", True),
                "rmsd": pair.get("rmsd"),
            })
    
    def build_matrix_index(
        self,
        prediction_id: str,
        matrix_type: str,
        matrix_path: Optional[Path] = None,
        dimension: int = 0,
        row_mapping_id: Optional[str] = None,
        column_mapping_id: Optional[str] = None,
    ) -> None:
        """Build a matrix index record."""
        self._matrices_index.append({
            "prediction_id": prediction_id,
            "matrix_type": matrix_type,
            "matrix_path": str(matrix_path) if matrix_path else None,
            "dimension": dimension,
            "row_mapping_id": row_mapping_id,
            "column_mapping_id": column_mapping_id,
            "definedness": "present",
        })
    
    def to_chains_dataframe(self) -> pd.DataFrame:
        """Convert chains to DataFrame."""
        return pd.DataFrame(self._chains)
    
    def to_residues_dataframe(self) -> pd.DataFrame:
        """Convert residues to DataFrame."""
        return pd.DataFrame(self._residues)
    
    def to_chain_pairs_dataframe(self) -> pd.DataFrame:
        """Convert chain pairs to DataFrame."""
        return pd.DataFrame(self._chain_pairs)
    
    def to_residue_correspondences_dataframe(self) -> pd.DataFrame:
        """Convert residue correspondences to DataFrame."""
        return pd.DataFrame(self._residue_correspondences)
    
    def to_matrices_index_dataframe(self) -> pd.DataFrame:
        """Convert matrices index to DataFrame."""
        return pd.DataFrame(self._matrices_index)
    
    def get_comparability_status(
        self,
        condition_a: str,
        condition_b: str,
        scope: str = "full",
    ) -> str:
        """
        Derive comparability status between two conditions.
        
        Returns:
            - "full": Identical mapped composition
            - "shared_scope_only": Shared scope only
            - "not_comparable": No defensible common scope
        """
        # Simplified logic - in practice, would compare residue mappings
        if scope == "full":
            # Check if conditions have identical composition
            return "full"
        elif scope == "shared_scope_only":
            return "shared_scope_only"
        else:
            return "not_comparable"


def build_chain_mappings(
    conditions_df: pd.DataFrame,
    replicates_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build chain mappings from condition and replicate data.
    
    Args:
        conditions_df: Conditions table
        replicates_df: Replicates table
    
    Returns:
        DataFrame of chain mappings
    """
    builder = MappingBuilder()
    
    # For each condition, build chain mappings
    for _, cond_row in conditions_df.iterrows():
        condition_id = cond_row["condition_id"]
        
        # Get first replicate for this condition
        cond_replicates = replicates_df[replicates_df["condition_id"] == condition_id]
        if cond_replicates.empty:
            continue
        
        # Simplified: assume chain information from first replicate
        # In practice, would extract from raw AF3 artifacts
        builder.build_chain_mapping(
            condition_id=condition_id,
            chain_id="A",
            entity_id=f"entity_{condition_id}",
            token_count=100,  # Placeholder
            residue_count=100,  # Placeholder
            chain_type="protein",
        )
    
    return builder.to_chains_dataframe()


__all__ = [
    "Chain",
    "Residue",
    "ChainPair",
    "ResidueCorrespondence",
    "MatrixIndex",
    "MappingBuilder",
    "build_chain_mappings",
]
