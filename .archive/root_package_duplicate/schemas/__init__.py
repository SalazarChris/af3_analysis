"""
Schemas module for AF3 Confidence Analysis Pipeline.

Contains data model definitions, table schemas, and validation logic.
"""

from af3_analysis.schemas.enums import (
    AnalysisMode,
    Definedness,
    ConditionStatus,
    ReplicateStatus,
    QCSeverity,
    ComparabilityStatus,
    MetricResolutionStatus,
)

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
)

from af3_analysis.schemas.tables import (
    ExperimentTable,
    ConditionTable,
    ReplicateTable,
    MeasurementTable,
    RegistryTable,
    QCFindingsTable,
    ExclusionTable,
    MappingTable,
)

from af3_analysis.schemas.validation import (
    validate_record,
    validate_table,
    ValidationError,
)

__all__ = [
    "AnalysisMode",
    "Definedness",
    "ConditionStatus",
    "ReplicateStatus",
    "QCSeverity",
    "ComparabilityStatus",
    "MetricResolutionStatus",
    "Experiment",
    "Condition",
    "Replicate",
    "BestModelReference",
    "Measurement",
    "RegistryEntry",
    "QCFinding",
    "Exclusion",
    "Mapping",
    "ExperimentTable",
    "ConditionTable",
    "ReplicateTable",
    "MeasurementTable",
    "RegistryTable",
    "QCFindingsTable",
    "ExclusionTable",
    "MappingTable",
    "validate_record",
    "validate_table",
    "ValidationError",
]
