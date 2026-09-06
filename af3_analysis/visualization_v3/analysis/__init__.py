"""
V3 Analysis Submodules.

Contains confidence-geometry integration, effect sizes,
factorial analysis, and seed reproducibility.
"""

from .confidence import (
    integrate_confidence_geometry,
    ConfidenceGeometryAnalysis,
)
from .effects import (
    calculate_structural_effect_sizes,
    StructuralEffectSize,
)
from .factorial import (
    factorial_structural_analysis,
    FactorialContrast,
)
from .reproducibility import (
    calculate_seed_reproducibility,
    SeedReproducibility,
)

__all__ = [
    "integrate_confidence_geometry",
    "ConfidenceGeometryAnalysis",
    "calculate_structural_effect_sizes",
    "StructuralEffectSize",
    "factorial_structural_analysis",
    "FactorialContrast",
    "calculate_seed_reproducibility",
    "SeedReproducibility",
]
