"""
Enumerated types for structural analysis subsystem.

Defines comparison statuses, metric directions, and structural quality codes.
"""

from enum import Enum


class ComparisonStatus(Enum):
    """Status of a pairwise structural comparison."""
    COMPARABLE = "comparable"
    PARTIALLY_COMPARABLE = "partially_comparable"
    NOT_COMPARABLE = "not_comparable"


class ComparisonReason(Enum):
    """Detailed reason for a comparison status."""
    VALID = "valid"
    INSUFFICIENT_ATOMS = "insufficient_atoms"
    MISSING_CHAIN = "missing_chain"
    SEQUENCE_MISMATCH = "sequence_mismatch"
    NO_COMMON_RESIDUES = "no_common_residues"
    MISSING_REFERENCE = "missing_reference"
    UNSUPPORTED_ENTITY_TYPE = "unsupported_entity_type"
    INSUFFICIENT_REPLICATES = "insufficient_replicates"
    PARSE_ERROR = "parse_error"
    MODEL_NUMBER_MISMATCH = "model_number_mismatch"
    EMPTY_STRUCTURE = "empty_structure"


class AtomSelection(Enum):
    """Atom selection mode for alignment and metrics."""
    CA = "ca"
    BACKBONE = "backbone"
    ALL_HEAVY = "all_heavy"
    ALL = "all"


class AlignmentMethod(Enum):
    """Structural alignment method."""
    KABSCH = "kabsch"


class ComparisonMode(Enum):
    """How comparisons are structured across conditions and seeds."""
    MATCHED_SEED = "matched_seed"
    ALL_VS_ALL = "all_vs_all"
    SEED_PAIRED = "seed_paired"


class MetricDirection(Enum):
    """Direction of metric change (lower/higher is geometrically meaningful)."""
    LOWER_BETTER = "lower_better"
    HIGHER_BETTER = "higher_better"


class MetricScope(Enum):
    """Scope of a structural metric."""
    GLOBAL = "global"
    CHAIN = "chain"
    INTERFACE = "interface"
    REGION = "region"


class MetricApplicability(Enum):
    """Whether a metric applies to a given structural pair."""
    APPLICABLE = "applicable"
    NOT_APPLICABLE = "not_applicable"
    PARTIALLY_APPLICABLE = "partially_applicable"


class ParseStatus(Enum):
    """Status of CIF file parsing."""
    SUCCESS = "success"
    PARSE_ERROR = "parse_error"
    EMPTY_STRUCTURE = "empty_structure"
    MODEL_NUMBER_MISMATCH = "model_number_mismatch"
