"""
Configuration for the structural analysis subsystem.

All fields have project-agnostic defaults. No biological assumptions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class RegionDefinition:
    """
    User-supplied structural region definition.

    All fields optional; intersection of specified filters selects atoms.
    """
    label: Optional[str] = None
    chain_ids: Optional[List[str]] = None
    entity_ids: Optional[List[int]] = None
    entity_types: Optional[List[str]] = None
    seq_range: Optional[tuple] = None  # (start, end) inclusive auth_seq_id
    seq_list: Optional[List[int]] = None  # explicit auth_seq_ids
    atom_names: Optional[List[str]] = None


@dataclass(frozen=True)
class StructuralConfig:
    """
    Configuration for structural analysis.

    All defaults are project-agnostic.
    """
    # Enable/disable
    enabled: bool = False

    # Reference
    reference_strategy: str = "explicit_reference"  # explicit_reference, first_condition, pairwise_all, pooled_reference
    reference_condition: Optional[str] = None
    reference_path: Optional[Path] = None

    # Atom selection
    atom_selection: str = "ca"  # ca, backbone, all_heavy, all
    min_common_atoms: int = 3

    # Alignment
    alignment_method: str = "kabsch"
    min_sequence_identity: float = 0.5

    # Contact detection
    contact_cutoff_angstrom: float = 8.0

    # Region definitions (optional)
    regions: List[RegionDefinition] = field(default_factory=list)

    # Metrics to compute (empty = all applicable)
    enabled_metrics: List[str] = field(default_factory=list)

    # Missing-data policy
    drop_incomparable: bool = False

    # Output settings
    output_csv: bool = True
    output_parquet: bool = True

    # Seed settings
    comparison_mode: str = "matched_seed"
    min_seeds: int = 2

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StructuralConfig":
        """Create from a dictionary (e.g., parsed JSON config)."""
        regions = []
        for r in data.get("regions", []):
            regions.append(RegionDefinition(**r))

        ref_path = data.get("reference_path")
        if ref_path is not None:
            ref_path = Path(ref_path)

        return cls(
            enabled=bool(data.get("enabled", False)),
            reference_strategy=str(data.get("reference_strategy", "explicit_reference")),
            reference_condition=data.get("reference_condition"),
            reference_path=ref_path,
            atom_selection=str(data.get("atom_selection", "ca")),
            min_common_atoms=int(data.get("min_common_atoms", 10)),
            alignment_method=str(data.get("alignment_method", "kabsch")),
            min_sequence_identity=float(data.get("min_sequence_identity", 0.5)),
            contact_cutoff_angstrom=float(data.get("contact_cutoff_angstrom", 8.0)),
            regions=regions,
            enabled_metrics=list(data.get("enabled_metrics", [])),
            drop_incomparable=bool(data.get("drop_incomparable", False)),
            output_csv=bool(data.get("output_csv", True)),
            output_parquet=bool(data.get("output_parquet", True)),
            comparison_mode=str(data.get("comparison_mode", "matched_seed")),
            min_seeds=int(data.get("min_seeds", 2)),
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {
            "enabled": self.enabled,
            "reference_strategy": self.reference_strategy,
            "reference_condition": self.reference_condition,
            "reference_path": str(self.reference_path) if self.reference_path else None,
            "atom_selection": self.atom_selection,
            "min_common_atoms": self.min_common_atoms,
            "alignment_method": self.alignment_method,
            "min_sequence_identity": self.min_sequence_identity,
            "contact_cutoff_angstrom": self.contact_cutoff_angstrom,
            "regions": [],
            "enabled_metrics": self.enabled_metrics,
            "drop_incomparable": self.drop_incomparable,
            "output_csv": self.output_csv,
            "output_parquet": self.output_parquet,
            "comparison_mode": self.comparison_mode,
            "min_seeds": self.min_seeds,
        }
        for r in self.regions:
            d = {}
            if r.label is not None:
                d["label"] = r.label
            if r.chain_ids is not None:
                d["chain_ids"] = r.chain_ids
            if r.entity_ids is not None:
                d["entity_ids"] = r.entity_ids
            if r.entity_types is not None:
                d["entity_types"] = r.entity_types
            if r.seq_range is not None:
                d["seq_range"] = list(r.seq_range)
            if r.seq_list is not None:
                d["seq_list"] = r.seq_list
            if r.atom_names is not None:
                d["atom_names"] = r.atom_names
            result["regions"].append(d)
        return result
