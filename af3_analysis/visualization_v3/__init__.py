"""
V3 Visualization Pipeline for AF3 Structural Analysis.

This module provides an additive visualization layer that extends the existing
AF3 analysis pipeline with enhanced structural/geometric analysis capabilities.

The V3 pipeline:
- Consumes existing pipeline outputs (tables, structures)
- Does NOT modify, rewrite, or replace any existing module
- Writes exclusively to run_dir/v3/
- Implements 20 figures (V3-F01 through V3-F20)

Design principles:
- Scientific terminology: structural difference, geometric displacement,
  condition-associated change, predicted structural cluster
- Configuration-driven: no hard-coded biological assumptions
- Failure-safe: single figure failure does not terminate pipeline
- Protein-agnostic: works for any AF3 project with structural data
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [
    "V3Config",
    "V3StructuralConfig",
    "V3ClusteringConfig",
]
