"""
Preprocessing module for AF3 Confidence Analysis Pipeline.

Provides metadata loading, canonicalization, QC, aggregation, and mappings.
"""

from af3_analysis.preprocessing.metadata import load_condition_metadata

__all__ = [
    "load_condition_metadata",
]
