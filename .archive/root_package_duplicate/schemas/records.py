"""
Record definitions for AF3 Confidence Analysis Pipeline.

Validated row contracts for experiments, conditions, replicates, measurements,
registry entries, QC findings, exclusions, and mappings.
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
from af3_analysis.schemas.enums import (
    ConditionStatus,
    ReplicateStatus,
    ComparabilityStatus,
    MetricScope,
    MetricPortfolioCategory,
    MetricDirection,
)


@dataclass(frozen=True)
class Experiment:
    """
    One row per scientifically coherent experiment/system series.
    
    Primary key: experiment_id
    """
    experiment_id: str
    system_id: str
    protocol_id: str
    input_root: Path
    analysis_run_id: str
    design_description: str
    metadata_schema_version: str = "1.0"


@dataclass(frozen=True)
class Condition:
    """
    One row per unique input condition within an experiment.
    
    Primary key: (condition_id, experiment_id)
    """
    condition_id: str
    experiment_id: str
    condition_label: str
    input_signature: str
    condition_name: Optional[str] = None
    condition_folder: Optional[str] = None
    expected_composition: Optional[str] = None
    comparability_group: Optional[str] = None
    status: ConditionStatus = ConditionStatus.PLANNED
    factors_json: Optional[str] = None  # JSON string for factor metadata
    source_row_id: Optional[str] = None  # Original CSV row ID for provenance


@dataclass(frozen=True)
class Replicate:
    """
    One row per genuine prediction, keyed by (condition_id, seed, sample).
    
    Primary key: (condition_id, seed, sample)
    """
    condition_id: str
    seed: int
    sample: int
    prediction_id: str
    source_confidences_path: Optional[Path] = None
    source_summary_path: Optional[Path] = None
    source_model_path: Optional[Path] = None
    source_ranking_path: Optional[Path] = None
    run_status: ReplicateStatus = ReplicateStatus.VALID
    model_rank: Optional[int] = None
    qc_flags: Optional[str] = None  # Semicolon-delimited QC finding IDs
    source_checksum: Optional[str] = None
    af3_version: Optional[str] = None


@dataclass(frozen=True)
class BestModelReference:
    """
    Reference to best model for a prediction.
    
    Primary key: prediction_id
    """
    prediction_id: str
    best_model_path: Path
    reference_type: str  # e.g., "top_level", "best_ranked", "manual_override"


@dataclass(frozen=True)
class Measurement:
    """
    Long-form metric measurement table.
    
    Primary key: (prediction_id, metric_id, scope_type, scope_id)
    """
    prediction_id: str
    metric_id: str
    scope_type: MetricScope
    scope_id: str  # Empty for global, chain ID, ordered pair like "A|B", or canonical residue key
    value: Optional[float] = None
    definedness: str = "present"  # From Definedness enum
    source_precision: Optional[int] = None
    source_name: Optional[str] = None  # Raw JSON/CSV field name
    source_checksum: Optional[str] = None


@dataclass(frozen=True)
class RegistryEntry:
    """
    Metric registry entry defining a metric's properties.
    
    Primary key: metric_id
    """
    metric_id: str
    display_name: str
    source: str  # e.g., "confidence_json", "summary_json", "ranking_csv"
    scope: MetricScope
    units: str
    direction: MetricDirection
    portfolio_category: MetricPortfolioCategory
    precision: Optional[int] = None
    expected_range: Optional[str] = None  # e.g., "0-100", "0-1"
    caveats: Optional[str] = None  # Important notes about the metric
    resolution_status: str = "unresolved"
    resolution_reason: Optional[str] = None


@dataclass(frozen=True)
class QCFinding:
    """
    QC finding with rule ID and evidence.
    
    Primary key: finding_id
    """
    finding_id: str
    rule_id: str
    severity: str  # From QCSeverity enum
    affected_keys: str  # Comma-separated primary keys
    evidence: str
    action: str  # e.g., "exclude", "flag", "ignore"
    timestamp: datetime = field(default_factory=datetime.utcnow)


@dataclass(frozen=True)
class Exclusion:
    """
    Record of an excluded prediction with rule and evidence.
    
    Primary key: exclusion_id
    """
    exclusion_id: str
    primary_key: str  # The primary key of the excluded record
    rule_id: str
    source_path: str
    evidence: str
    resolution: str  # e.g., "duplicate", "malformed", "unidentifiable_identity"
    timestamp: datetime = field(default_factory=datetime.utcnow)


@dataclass(frozen=True)
class Mapping:
    """
    Mapping record for chain/residue correspondence.
    
    Primary key: mapping_id
    """
    mapping_id: str
    condition_id: str
    chain_id: str
    entity_id: str
    token_count: int
    residue_count: int
    chain_type: str  # e.g., "protein", "dna", "rna", "ligand"
    source_checksum: Optional[str] = None
    resolved: bool = True
    resolution_reason: Optional[str] = None


@dataclass(frozen=True)
class ChainPair:
    """
    Chain pair record for interface analysis.
    
    Primary key: (prediction_id, chain_pair_id)
    """
    prediction_id: str
    chain_pair_id: str  # e.g., "A|B"
    chain_a: str
    chain_b: str
    pair_iptm: Optional[float] = None
    pair_pae_mean: Optional[float] = None
    pair_pde_mean: Optional[float] = None
    interface_contact_confidence: Optional[float] = None
    definedness: str = "present"


__all__ = [
    "Experiment",
    "Condition",
    "Replicate",
    "BestModelReference",
    "Measurement",
    "RegistryEntry",
    "QCFinding",
    "Exclusion",
    "Mapping",
    "ChainPair",
]
