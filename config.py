"""
Configuration management for AF3 Confidence Analysis Pipeline.

Parses and validates a JSON configuration file, then creates
an immutable run context with resolved paths and metadata.
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set
from dataclasses import dataclass, field, asdict
import hashlib


class ConfigError(Exception):
    """Base configuration error."""
    pass


class ConfigValidationError(ConfigError):
    """Configuration validation error."""
    pass


class ConfigPathError(ConfigError):
    """Configuration path resolution error."""
    pass


@dataclass(frozen=True)
class AnalysisMode:
    """Analysis mode configuration."""
    name: str
    canonical_analysis: bool = True
    requires_raw_af3_root: bool = True
    
    @classmethod
    def from_name(cls, name: str) -> "AnalysisMode":
        if name == "canonical_analysis":
            return cls(name=name, canonical_analysis=True, requires_raw_af3_root=True)
        elif name == "legacy_summary_descriptive":
            return cls(name=name, canonical_analysis=False, requires_raw_af3_root=False)
        else:
            raise ConfigValidationError(f"Unknown analysis mode: {name}")


@dataclass(frozen=True)
class ResamplingConfig:
    """Bootstrap and permutation resampling configuration."""
    bootstrap_iterations: int = 2000
    permutation_iterations: int = 2000
    confidence_level: float = 0.95
    analysis_seed: Optional[int] = None
    
    def __post_init__(self):
        if not 0 < self.confidence_level < 1:
            raise ConfigValidationError(
                f"confidence_level must be in (0, 1), got {self.confidence_level}"
            )
        if self.bootstrap_iterations < 100:
            raise ConfigValidationError(
                f"bootstrap_iterations must be >= 100, got {self.bootstrap_iterations}"
            )
        if self.permutation_iterations < 100:
            raise ConfigValidationError(
                f"permutation_iterations must be >= 100, got {self.permutation_iterations}"
            )


@dataclass(frozen=True)
class InputPaths:
    """Input directory paths."""
    raw_af3_root: Optional[Path] = None
    legacy_csv_root: Optional[Path] = None
    mapping_files: Optional[Dict[str, Path]] = None


@dataclass(frozen=True)
class OutputPaths:
    """Output directory paths."""
    output_root: Path
    run_id: str
    
    @property
    def run_dir(self) -> Path:
        return self.output_root / self.run_id
    
    @property
    def manifest_dir(self) -> Path:
        return self.run_dir / "manifest"
    
    @property
    def tables_dir(self) -> Path:
        return self.run_dir / "tables"
    
    @property
    def figures_dir(self) -> Path:
        return self.run_dir / "figures"
    
    @property
    def reports_dir(self) -> Path:
        return self.run_dir / "reports"
    
    @property
    def logs_dir(self) -> Path:
        return self.run_dir / "logs"
    
    @property
    def parquet_dir(self) -> Path:
        return self.run_dir / "parquet"


@dataclass
class AnalysisConfig:
    """Main analysis configuration."""
    run_id: str
    output_root: Path
    random_seed: int
    
    # Resampling configuration
    bootstrap_iterations: int = 2000
    permutation_iterations: int = 2000
    confidence_level: float = 0.95
    
    # Analysis mode
    analysis_mode: str = "canonical_analysis"
    
    # Output switches
    export_parquet: bool = True
    export_csv: bool = True
    generate_figures: bool = True
    generate_reports: bool = True
    
    # Input paths
    raw_af3_root: Optional[Path] = None
    legacy_csv_root: Optional[Path] = None
    
    # Optional settings
    coordinate_analysis_enabled: bool = False
    mapping_id: Optional[str] = None
    reference_condition: Optional[str] = None
    
    @property
    def resampling(self) -> ResamplingConfig:
        return ResamplingConfig(
            bootstrap_iterations=self.bootstrap_iterations,
            permutation_iterations=self.permutation_iterations,
            confidence_level=self.confidence_level,
            analysis_seed=self.random_seed,
        )
    
    @property
    def input_paths(self) -> InputPaths:
        return InputPaths(
            raw_af3_root=self.raw_af3_root,
            legacy_csv_root=self.legacy_csv_root,
        )
    
    @property
    def output_paths(self) -> OutputPaths:
        return OutputPaths(output_root=self.output_root, run_id=self.run_id)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AnalysisConfig":
        """Create config from dictionary (typically JSON)."""
        # Required fields
        run_id = data.get("run_id")
        if not run_id:
            raise ConfigValidationError("Missing required field: run_id")
        
        output_root = data.get("output_root")
        if not output_root:
            raise ConfigValidationError("Missing required field: output_root")
        
        random_seed = data.get("random_seed")
        if random_seed is None:
            raise ConfigValidationError("Missing required field: random_seed")
        
        # Parse paths
        output_root_path = Path(output_root)
        raw_af3_root = Path(data["raw_af3_root"]) if data.get("raw_af3_root") else None
        legacy_csv_root = Path(data["legacy_csv_root"]) if data.get("legacy_csv_root") else None
        
        # Optional with defaults
        bootstrap_iterations = data.get("bootstrap_iterations", 2000)
        permutation_iterations = data.get("permutation_iterations", 2000)
        confidence_level = data.get("confidence_level", 0.95)
        analysis_mode = data.get("analysis_mode", "canonical_analysis")
        coordinate_analysis_enabled = data.get("coordinate_analysis_enabled", False)
        mapping_id = data.get("mapping_id")
        reference_condition = data.get("reference_condition")
        
        return cls(
            run_id=str(run_id),
            output_root=output_root_path,
            random_seed=int(random_seed),
            bootstrap_iterations=int(bootstrap_iterations),
            permutation_iterations=int(permutation_iterations),
            confidence_level=float(confidence_level),
            analysis_mode=str(analysis_mode),
            coordinate_analysis_enabled=bool(coordinate_analysis_enabled),
            mapping_id=mapping_id,
            reference_condition=reference_condition,
            raw_af3_root=raw_af3_root,
            legacy_csv_root=legacy_csv_root,
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary for JSON serialization."""
        return {
            "run_id": self.run_id,
            "output_root": str(self.output_root),
            "random_seed": self.random_seed,
            "bootstrap_iterations": self.bootstrap_iterations,
            "permutation_iterations": self.permutation_iterations,
            "confidence_level": self.confidence_level,
            "analysis_mode": self.analysis_mode,
            "raw_af3_root": str(self.raw_af3_root) if self.raw_af3_root else None,
            "legacy_csv_root": str(self.legacy_csv_root) if self.legacy_csv_root else None,
            "coordinate_analysis_enabled": self.coordinate_analysis_enabled,
            "mapping_id": self.mapping_id,
            "reference_condition": self.reference_condition,
        }
    
    def to_manifest(self) -> Dict[str, Any]:
        """Create run manifest data."""
        from . import __version__
        return {
            "run_id": self.run_id,
            "analysis_mode": self.analysis_mode,
            "created_at": datetime.utcnow().isoformat() + "Z",
            "version": __version__,
            "config_hash": self._compute_config_hash(),
            "inputs": self.input_paths.to_dict() if hasattr(self.input_paths, "to_dict") else self.to_dict(),
            "outputs": {
                "output_root": str(self.output_root),
                "run_dir": str(self.output_paths.run_dir),
            },
        }
    
    def _compute_config_hash(self) -> str:
        """Compute SHA256 hash of config for provenance."""
        config_str = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.sha256(config_str.encode()).hexdigest()


def load_config(config_path: str) -> AnalysisConfig:
    """Load and validate configuration from JSON file."""
    path = Path(config_path)
    if not path.exists():
        raise ConfigPathError(f"Configuration file not found: {path}")
    
    try:
        with open(path, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise ConfigValidationError(f"Invalid JSON in configuration file: {e}")
    
    return AnalysisConfig.from_dict(data)


def create_run_directory(config: AnalysisConfig) -> Path:
    """Create run directory structure from configuration."""
    run_dir = config.output_paths.run_dir
    run_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    config.output_paths.manifest_dir.mkdir(exist_ok=True)
    config.output_paths.tables_dir.mkdir(exist_ok=True)
    config.output_paths.figures_dir.mkdir(exist_ok=True)
    config.output_paths.reports_dir.mkdir(exist_ok=True)
    config.output_paths.logs_dir.mkdir(exist_ok=True)
    config.output_paths.parquet_dir.mkdir(exist_ok=True)
    
    return run_dir


def create_config_interactive(
    raw_af3_root: str,
    output_dir: str = "af3_analysis_output",
    run_id: Optional[str] = None,
    bootstrap_iterations: int = 2000,
    permutation_iterations: int = 2000,
    confidence_level: float = 0.95,
    random_seed: int = 42,
    reference_condition: Optional[str] = None,
) -> AnalysisConfig:
    """Create a configuration interactively, e.g., from UI prompts.
    
    If run_id is not provided, it is auto-generated based on the current timestamp.
    """
    if not run_id:
        from datetime import datetime
        run_id = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    return AnalysisConfig(
        run_id=run_id,
        output_root=Path(output_dir),
        random_seed=random_seed,
        bootstrap_iterations=bootstrap_iterations,
        permutation_iterations=permutation_iterations,
        confidence_level=confidence_level,
        raw_af3_root=Path(raw_af3_root),
        reference_condition=reference_condition,
    )


__all__ = [
    "AnalysisConfig",
    "ConfigError",
    "ConfigValidationError",
    "ConfigPathError",
    "load_config",
    "create_run_directory",
    "create_config_interactive",
]
