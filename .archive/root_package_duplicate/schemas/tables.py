"""
Table schemas and validation for AF3 Confidence Analysis Pipeline.

Implements table contracts with primary keys, foreign keys, types, and validators.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Set
from dataclasses import dataclass, field
import pandas as pd

from af3_analysis.schemas.records import (
    Experiment,
    Condition,
    Replicate,
    BestModelReference,
    Measurement,
    RegistryEntry,
    QCFinding,
    Exclusion,
    Mapping,
    ChainPair,
)
from af3_analysis.schemas.enums import (
    Definedness,
    ConditionStatus,
    ReplicateStatus,
    ComparabilityStatus,
    MetricScope,
    MetricPortfolioCategory,
    MetricDirection,
    AnalysisLevel,
)


class ValidationError(Exception):
    """Validation error for table records."""
    
    def __init__(self, message: str, record: Optional[Dict[str, Any]] = None):
        self.message = message
        self.record = record or {}
        super().__init__(self.message)


class TableSchema:
    """Base class for table schemas with validation."""
    
    @property
    def primary_keys(self) -> List[str]:
        """Return primary key columns."""
        raise NotImplementedError
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        """Return foreign key relationships."""
        raise NotImplementedError
    
    @property
    def required_columns(self) -> Set[str]:
        """Return required column names."""
        raise NotImplementedError
    
    @property
    def allowed_columns(self) -> Set[str]:
        """Return all allowed column names."""
        raise NotImplementedError
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        """Validate a single record."""
        raise NotImplementedError
    
    def validate_table(self, df: pd.DataFrame) -> None:
        """Validate a DataFrame table."""
        raise NotImplementedError


@dataclass
class ExperimentTable(TableSchema):
    """Table of experiments."""
    records: List[Experiment] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["experiment_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {"experiment_id": "conditions.experiment_id"}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"experiment_id", "system_id", "protocol_id", "input_root", 
                "analysis_run_id", "design_description"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"metadata_schema_version"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "experiment_id" not in record:
            raise ValidationError("Missing experiment_id", record)
        if "system_id" not in record:
            raise ValidationError("Missing system_id", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())
        # Check for duplicate experiment_ids
        if df.duplicated(subset=["experiment_id"]).any():
            raise ValidationError("Duplicate experiment_id found", {})


@dataclass
class ConditionTable(TableSchema):
    """Table of conditions."""
    records: List[Condition] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["condition_id", "experiment_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {"experiment_id": "experiments.experiment_id"}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"condition_id", "experiment_id", "condition_label", "input_signature"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {
            "condition_name", "condition_folder", "expected_composition",
            "comparability_group", "status", "factors_json", "source_row_id"
        }
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "condition_id" not in record:
            raise ValidationError("Missing condition_id", record)
        if "experiment_id" not in record:
            raise ValidationError("Missing experiment_id", record)
        if "input_signature" not in record:
            raise ValidationError("Missing input_signature", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())
        # Check referential integrity
        # Check for duplicate (condition_id, experiment_id)


@dataclass
class ReplicateTable(TableSchema):
    """Table of replicates."""
    records: List[Replicate] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["condition_id", "seed", "sample"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {
            "condition_id": "conditions.condition_id",
            "experiment_id": "conditions.experiment_id",
        }
    
    @property
    def required_columns(self) -> Set[str]:
        return {"condition_id", "seed", "sample", "prediction_id"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {
            "source_confidences_path", "source_summary_path", "source_model_path",
            "source_ranking_path", "run_status", "model_rank", "qc_flags",
            "source_checksum", "af3_version"
        }
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "condition_id" not in record:
            raise ValidationError("Missing condition_id", record)
        if "seed" not in record:
            raise ValidationError("Missing seed", record)
        if "sample" not in record:
            raise ValidationError("Missing sample", record)
        if "prediction_id" not in record:
            raise ValidationError("Missing prediction_id", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())
        # Check for duplicate (condition_id, seed, sample)
        # Check for orphaned condition_ids


@dataclass
class MeasurementTable(TableSchema):
    """Table of measurements."""
    records: List[Measurement] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["prediction_id", "metric_id", "scope_type", "scope_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {
            "prediction_id": "replicates.prediction_id",
            "metric_id": "registry.metric_id",
        }
    
    @property
    def required_columns(self) -> Set[str]:
        return {"prediction_id", "metric_id", "scope_type", "definedness"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"value", "source_precision", "source_name", "source_checksum"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "prediction_id" not in record:
            raise ValidationError("Missing prediction_id", record)
        if "metric_id" not in record:
            raise ValidationError("Missing metric_id", record)
        if "definedness" not in record:
            raise ValidationError("Missing definedness", record)
        definedness = record.get("definedness")
        value = record.get("value")
        
        # Check definedness/value consistency
        if definedness in ["present"] and value is None:
            raise ValidationError("Value must be present when definedness is 'present'", record)
        if definedness not in ["present"] and value is not None:
            raise ValidationError("Value must be None when definedness is not 'present'", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())
        # Check referential integrity


@dataclass
class RegistryTable(TableSchema):
    """Table of registry entries."""
    records: List[RegistryEntry] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["metric_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"metric_id", "display_name", "source", "scope", "units", "direction", "portfolio_category"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"precision", "expected_range", "caveats", "resolution_status", "resolution_reason"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "metric_id" not in record:
            raise ValidationError("Missing metric_id", record)
        if "display_name" not in record:
            raise ValidationError("Missing display_name", record)
        if "scope" not in record:
            raise ValidationError("Missing scope", record)
        if "units" not in record:
            raise ValidationError("Missing units", record)
        if "direction" not in record:
            raise ValidationError("Missing direction", record)
        if "portfolio_category" not in record:
            raise ValidationError("Missing portfolio_category", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())
        # Check for duplicate metric_ids


@dataclass
class QCFindingsTable(TableSchema):
    """Table of QC findings."""
    records: List[QCFinding] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["finding_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"finding_id", "rule_id", "severity", "affected_keys", "action"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"evidence", "timestamp"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "finding_id" not in record:
            raise ValidationError("Missing finding_id", record)
        if "rule_id" not in record:
            raise ValidationError("Missing rule_id", record)
        if "severity" not in record:
            raise ValidationError("Missing severity", record)
        if "affected_keys" not in record:
            raise ValidationError("Missing affected_keys", record)
        if "action" not in record:
            raise ValidationError("Missing action", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())


@dataclass
class ExclusionTable(TableSchema):
    """Table of exclusions."""
    records: List[Exclusion] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["exclusion_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"exclusion_id", "primary_key", "rule_id", "source_path", "evidence", "resolution"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"timestamp"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "exclusion_id" not in record:
            raise ValidationError("Missing exclusion_id", record)
        if "primary_key" not in record:
            raise ValidationError("Missing primary_key", record)
        if "rule_id" not in record:
            raise ValidationError("Missing rule_id", record)
        if "source_path" not in record:
            raise ValidationError("Missing source_path", record)
        if "evidence" not in record:
            raise ValidationError("Missing evidence", record)
        if "resolution" not in record:
            raise ValidationError("Missing resolution", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())


@dataclass
class MappingTable(TableSchema):
    """Table of mappings."""
    records: List[Mapping] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["mapping_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"mapping_id", "condition_id", "chain_id", "entity_id", "token_count", "residue_count", "chain_type"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"source_checksum", "resolved", "resolution_reason"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "mapping_id" not in record:
            raise ValidationError("Missing mapping_id", record)
        if "condition_id" not in record:
            raise ValidationError("Missing condition_id", record)
        if "chain_id" not in record:
            raise ValidationError("Missing chain_id", record)
        if "entity_id" not in record:
            raise ValidationError("Missing entity_id", record)
        if "token_count" not in record:
            raise ValidationError("Missing token_count", record)
        if "residue_count" not in record:
            raise ValidationError("Missing residue_count", record)
        if "chain_type" not in record:
            raise ValidationError("Missing chain_type", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())


@dataclass
class ChainPairTable(TableSchema):
    """Table of chain pairs."""
    records: List[ChainPair] = field(default_factory=list)
    
    @property
    def primary_keys(self) -> List[str]:
        return ["prediction_id", "chain_pair_id"]
    
    @property
    def foreign_keys(self) -> Dict[str, str]:
        return {"prediction_id": "replicates.prediction_id"}
    
    @property
    def required_columns(self) -> Set[str]:
        return {"prediction_id", "chain_pair_id", "chain_a", "chain_b"}
    
    @property
    def allowed_columns(self) -> Set[str]:
        return self.required_columns | {"pair_iptm", "pair_pae_mean", "pair_pde_mean", 
                                        "interface_contact_confidence", "definedness"}
    
    def validate_record(self, record: Dict[str, Any]) -> None:
        if "prediction_id" not in record:
            raise ValidationError("Missing prediction_id", record)
        if "chain_pair_id" not in record:
            raise ValidationError("Missing chain_pair_id", record)
        if "chain_a" not in record:
            raise ValidationError("Missing chain_a", record)
        if "chain_b" not in record:
            raise ValidationError("Missing chain_b", record)
    
    def validate_table(self, df: pd.DataFrame) -> None:
        if df.empty:
            return
        for _, row in df.iterrows():
            self.validate_record(row.to_dict())


# Predefined metric registry for AF3 confidence metrics
AF3_METRIC_REGISTRY = {
    "ranking_score": RegistryEntry(
        metric_id="ranking_score",
        display_name="Ranking Score",
        source="ranking_csv",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.GLOBAL_CONFIDENCE,
        precision=None,
        expected_range="0-10",
        caveats="AF3 ranking score, higher = more confident prediction",
        resolution_status="resolved_available",
        resolution_reason="Extracted from ranking_scores.csv files",
    ),
    "ptm": RegistryEntry(
        metric_id="ptm",
        display_name="pTM Score",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.GLOBAL_CONFIDENCE,
        precision=4,
        expected_range="0-1",
        caveats="Predicted TM-score, higher = more confident global fold",
        resolution_status="resolved_available",
        resolution_reason="Available in confidence JSON files",
    ),
    "iptm": RegistryEntry(
        metric_id="iptm",
        display_name="ipTM Score",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.INTERFACE_CONFIDENCE,
        precision=4,
        expected_range="0-1",
        caveats="Predicted interface TM-score, higher = more confident interface",
        resolution_status="resolved_available",
        resolution_reason="Available in confidence JSON files",
    ),
    "plddt_mean": RegistryEntry(
        metric_id="plddt_mean",
        display_name="pLDDT Mean",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.GLOBAL_CONFIDENCE,
        precision=2,
        expected_range="0-100",
        caveats="Mean pLDDT across all atoms, higher = more confident",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
    "plddt_max": RegistryEntry(
        metric_id="plddt_max",
        display_name="pLDDT Max",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.GLOBAL_CONFIDENCE,
        precision=2,
        expected_range="0-100",
        caveats="Maximum pLDDT value",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
    "plddt_min": RegistryEntry(
        metric_id="plddt_min",
        display_name="pLDDT Min",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.LOWER_BETTER,
        portfolio_category=MetricPortfolioCategory.STRUCTURAL_UNCERTAINTY,
        precision=2,
        expected_range="0-100",
        caveats="Minimum pLDDT value",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
    "plddt_median": RegistryEntry(
        metric_id="plddt_median",
        display_name="pLDDT Median",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.GLOBAL_CONFIDENCE,
        precision=2,
        expected_range="0-100",
        caveats="Median pLDDT value",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
    "contact_prob_mean": RegistryEntry(
        metric_id="contact_prob_mean",
        display_name="Contact Probability Mean",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="unitless",
        direction=MetricDirection.HIGHER_BETTER,
        portfolio_category=MetricPortfolioCategory.INTERFACE_CONFIDENCE,
        precision=4,
        expected_range="0-1",
        caveats="Mean contact probability, higher = more confident contacts",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
    "pae_mean": RegistryEntry(
        metric_id="pae_mean",
        display_name="PAE Mean",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="Angstroms",
        direction=MetricDirection.LOWER_BETTER,
        portfolio_category=MetricPortfolioCategory.STRUCTURAL_UNCERTAINTY,
        precision=2,
        expected_range="0-100",
        caveats="Mean Predicted Aligned Error, lower = more confident alignment",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
    "pde_mean": RegistryEntry(
        metric_id="pde_mean",
        display_name="PDE Mean",
        source="confidence_json",
        scope=MetricScope.GLOBAL,
        units="Angstroms",
        direction=MetricDirection.LOWER_BETTER,
        portfolio_category=MetricPortfolioCategory.STRUCTURAL_UNCERTAINTY,
        precision=2,
        expected_range="0-100",
        caveats="Mean Predicted Distance Error, lower = more confident distances",
        resolution_status="resolved_available",
        resolution_reason="Calculated from confidence JSON",
    ),
}


def create_default_registry_table() -> RegistryTable:
    """Create a registry table from the default AF3 metric registry."""
    return RegistryTable(
        records=list(AF3_METRIC_REGISTRY.values())
    )


__all__ = [
    "ValidationError",
    "TableSchema",
    "ExperimentTable",
    "ConditionTable",
    "ReplicateTable",
    "MeasurementTable",
    "RegistryTable",
    "QCFindingsTable",
    "ExclusionTable",
    "MappingTable",
    "ChainPairTable",
    "AF3_METRIC_REGISTRY",
    "create_default_registry_table",
]
