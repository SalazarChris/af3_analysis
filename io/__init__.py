"""
I/O module for AF3 Confidence Analysis Pipeline.

Provides data loading, provenance tracking, and file I/O operations.
"""

from af3_analysis.io.phase1_loader import Phase1Loader
from af3_analysis.io.artifact_inventory import ArtifactInventory
from af3_analysis.io.provenance import Provenance, Checksum
from af3_analysis.io.wide_to_long import (
    seed_aggregated_to_long,
    long_to_seed_aggregated,
    load_seed_aggregated_long,
)

# Lazy import: ParquetStore requires pyarrow, which is optional.
# Importing it eagerly breaks modules that only need provenance/checksum
# when pyarrow is not installed.
try:
    from af3_analysis.io.parquet_store import ParquetStore
except ImportError:
    ParquetStore = None  # type: ignore[assignment,misc]

__all__ = [
    "Phase1Loader",
    "ArtifactInventory",
    "Provenance",
    "Checksum",
    "ParquetStore",
    "seed_aggregated_to_long",
    "long_to_seed_aggregated",
    "load_seed_aggregated_long",
]
