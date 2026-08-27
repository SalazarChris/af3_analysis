"""
Provenance tracking for AF3 Confidence Analysis Pipeline.

Tracks source files, checksums, and audit trail for reproducibility.
"""

import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
from datetime import datetime


@dataclass(frozen=True)
class Checksum:
    """File checksum for provenance."""
    algorithm: str
    value: str
    
    @classmethod
    def sha256(cls, value: str) -> "Checksum":
        return cls(algorithm="sha256", value=value)
    
    @classmethod
    def from_file(cls, path: Path, algorithm: str = "sha256") -> "Checksum":
        """Compute checksum of a file."""
        with open(path, "rb") as f:
            data = f.read()
            value = hashlib.sha256(data).hexdigest()
        return cls(algorithm=algorithm, value=value)


@dataclass(frozen=True)
class SourceFile:
    """Source file with provenance."""
    path: Path
    checksum: Checksum
    file_size: int
    modification_time: str
    file_type: str  # e.g., "csv", "json", "parquet"


@dataclass(frozen=True)
class RowProvenance:
    """Provenance for a single row."""
    source_file: SourceFile
    source_row_id: int
    source_row_checksum: Optional[str] = None


@dataclass
class Provenance:
    """
    Provenance tracking for data transformations.
    
    Maintains an audit trail of source files, transformations, and data lineage.
    """
    
    _sources: Dict[str, SourceFile] = field(default_factory=dict)
    _transformations: List[Dict[str, Any]] = field(default_factory=list)
    _created_at: str = field(default_factory=lambda: datetime.utcnow().isoformat() + "Z")
    
    def register_source(self, name: str, path: Path) -> SourceFile:
        """Register a source file and return its provenance."""
        checksum = Checksum.from_file(path)
        source = SourceFile(
            path=path,
            checksum=checksum,
            file_size=path.stat().st_size,
            modification_time=datetime.fromtimestamp(
                path.stat().st_mtime
            ).isoformat() + "Z",
            file_type=path.suffix[1:],  # Remove the dot
        )
        self._sources[name] = source
        return source
    
    def register_transformation(
        self,
        transformation_type: str,
        input_sources: List[str],
        output_name: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Record a data transformation."""
        self._transformations.append({
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "transformation_type": transformation_type,
            "input_sources": input_sources,
            "output_name": output_name,
            "metadata": metadata or {},
        })
    
    def get_source(self, name: str) -> Optional[SourceFile]:
        """Get source file provenance."""
        return self._sources.get(name)
    
    def get_all_sources(self) -> Dict[str, SourceFile]:
        """Get all registered sources."""
        return self._sources.copy()
    
    def get_transformations(self) -> List[Dict[str, Any]]:
        """Get transformation history."""
        return self._transformations.copy()
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert provenance to dictionary."""
        return {
            "created_at": self._created_at,
            "sources": {
                name: {
                    "path": str(src.path),
                    "checksum": src.checksum.value,
                    "file_size": src.file_size,
                    "modification_time": src.modification_time,
                    "file_type": src.file_type,
                }
                for name, src in self._sources.items()
            },
            "transformations": self._transformations,
        }


__all__ = [
    "Checksum",
    "SourceFile",
    "RowProvenance",
    "Provenance",
]
