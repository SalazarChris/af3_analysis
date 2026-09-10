"""
Enumerated types for AF3 Confidence Analysis Pipeline.

Defines all valid values for analysis modes, statuses, and classifications.
"""

from enum import Enum
from typing import get_args


class AnalysisMode(Enum):
    """Analysis mode configuration."""
    CANONICAL_ANALYSIS = "canonical_analysis"
    LEGACY_SUMMARY_DESCRIPTIVE = "legacy_summary_descriptive"


class Definedness(Enum):
    """
    Status of a metric value.
    
    - present: Valid numeric value exists.
    - undefined_by_composition: Cannot exist for this condition (e.g., inter-chain metric in monomer).
    - missing_technical: Should exist but extraction failed.
    - not_collected: Phase 1 schema did not extract it.
    """
    PRESENT = "present"
    UNDEFINED_BY_COMPOSITION = "undefined_by_composition"
    MISSING_TECHNICAL = "missing_technical"
    NOT_COLLECTED = "not_collected"


class ConditionStatus(Enum):
    """Status of a condition."""
    PLANNED = "planned"
    COMPLETE = "complete"
    INCOMPLETE = "incomplete"
    DUPLICATE_CANDIDATE = "duplicate_candidate"
    EXCLUDED_TECHNICAL = "excluded_technical"


class ReplicateStatus(Enum):
    """Status of a replicate/prediction."""
    VALID = "valid"
    FAILED = "failed"
    PARTIAL = "partial"
    DUPLICATE_REFERENCE = "duplicate_reference"
    EXCLUDED_TECHNICAL = "excluded_technical"


class QCSeverity(Enum):
    """Severity level for QC findings."""
    ERROR = "ERROR"
    WARNING = "WARNING"
    INFO = "INFO"


class ComparabilityStatus(Enum):
    """
    Comparability status for metric comparisons.
    
    - full: Identical mapped composition and denominator.
    - shared_scope_only: Global values are composition-dependent, but named shared scope is comparable.
    - not_comparable: No defensible common scope.
    """
    FULL = "full"
    SHARED_SCOPE_ONLY = "shared_scope_only"
    NOT_COMPARABLE = "not_comparable"


class MetricResolutionStatus(Enum):
    """Status of metric resolution in the registry."""
    RESOLVED_AVAILABLE = "resolved_available"
    RESOLVED_EXCLUDED = "resolved_excluded"
    UNAVAILABLE = "unavailable"
    EXPLORATORY = "exploratory"


class MetricScope(Enum):
    """Scope of a metric measurement."""
    GLOBAL = "global"
    CHAIN = "chain"
    CHAIN_PAIR = "chain_pair"
    RESIDUE = "residue"
    ATOM = "atom"
    MATRIX_SUMMARY = "matrix_summary"


class MetricPortfolioCategory(Enum):
    """Portfolio category for metrics."""
    GLOBAL_CONFIDENCE = "global_confidence"
    STRUCTURAL_UNCERTAINTY = "structural_uncertainty"
    INTERFACE_CONFIDENCE = "interface_confidence"
    LOCAL_CONFIDENCE = "local_confidence"
    DISORDER_FLEXIBILITY = "disorder_flexibility"
    PREDICTION_CONSISTENCY = "prediction_consistency"


class MetricDirection(Enum):
    """Direction of metric change (higher/lower is better)."""
    HIGHER_BETTER = "higher_better"
    LOWER_BETTER = "lower_better"
    NEUTRAL = "neutral"


class AnalysisLevel(Enum):
    """Level at which a metric is analyzed."""
    SAMPLE = "sample"
    SEED = "seed"


__all__ = [
    "AnalysisMode",
    "Definedness",
    "ConditionStatus",
    "ReplicateStatus",
    "QCSeverity",
    "ComparabilityStatus",
    "MetricResolutionStatus",
    "MetricScope",
    "MetricPortfolioCategory",
    "MetricDirection",
    "AnalysisLevel",
]
