"""
Structural/Geometric Analysis Subsystem
========================================

Generic 3D structural analysis layer for AF3 predicted structures.

Provides mmCIF parsing, normalized structural representation,
alignment, pluggable metrics, cross-condition comparison, and
normalized output tables.

This subsystem is project-agnostic and must not contain
hard-coded biological assumptions.
"""

from af3_analysis.structural.config import StructuralConfig
from af3_analysis.structural.representation import (
    NormalisedStructure,
    AtomRecord,
    ResidueRecord,
    EntityInfo,
    ChainInfo,
)
from af3_analysis.structural.alignment import (
    AlignmentResult,
    align_structures,
)
from af3_analysis.structural.comparison import (
    StructuralComparison,
    compare_conditions,
)
from af3_analysis.structural.tables import write_structural_tables

__all__ = [
    "StructuralConfig",
    "NormalisedStructure",
    "AtomRecord",
    "ResidueRecord",
    "EntityInfo",
    "ChainInfo",
    "AlignmentResult",
    "align_structures",
    "StructuralComparison",
    "compare_conditions",
    "write_structural_tables",
]
