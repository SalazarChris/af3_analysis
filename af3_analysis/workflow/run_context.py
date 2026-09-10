"""
Run context management for AF3 Confidence Analysis Pipeline.

Manages the mutable state during analysis while preserving immutability
of configuration and final outputs.
"""

import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field


@dataclass(frozen=True)
class AnalysisStage:
    """Analysis stage definition."""
    name: str
    order: int
    required: bool = True
    depends_on: Optional[str] = None


@dataclass
class RunContext:
    """
    Mutable run context tracking analysis progress.
    
    This is the primary interface for modules to access:
    - Configuration (immutable)
    - I/O paths
    - Logging
    - Data tables (mutable during analysis)
    - Stage tracking
    """
    config: Any  # AnalysisConfig
    logger: Any  # RunLogger
    
    # Stage tracking
    _stage_order: List[str] = field(default_factory=list)
    _stage_status: Dict[str, str] = field(default_factory=dict)
    _stage_data: Dict[str, Any] = field(default_factory=dict)
    
    # Data tables (created during analysis)
    tables: Dict[str, Any] = field(default_factory=dict)
    
    # Metrics registry
    resolved_registry: Any = None
    
    # QC findings and exclusions
    qc_findings: List[Dict[str, Any]] = field(default_factory=list)
    exclusion_log: List[Dict[str, Any]] = field(default_factory=list)
    
    def __post_init__(self):
        # Initialize stage tracking
        self._stage_status = {
            "config_loaded": "completed",
            "run_directory_created": "completed",
            "input_inventory": "pending",
            "canonicalization": "pending",
            "qc_gates": "pending",
            "seed_aggregation": "pending",
            "statistical_analysis": "pending",
            "reporting": "pending",
        }
    
    def start_stage(self, stage_name: str) -> None:
        """Mark a stage as starting."""
        if stage_name not in self._stage_status:
            self._stage_order.append(stage_name)
            self._stage_status[stage_name] = "in_progress"
        self.logger.info(f"Starting stage: {stage_name}")
    
    def complete_stage(self, stage_name: str, **metadata: Any) -> None:
        """Mark a stage as completed."""
        if stage_name in self._stage_status:
            self._stage_status[stage_name] = "completed"
        self.logger.info(f"Completed stage: {stage_name}", **metadata)
    
    def fail_stage(self, stage_name: str, error: str, **metadata: Any) -> None:
        """Mark a stage as failed."""
        if stage_name in self._stage_status:
            self._stage_status[stage_name] = "failed"
        self.logger.error(f"Failed stage: {stage_name}: {error}", **metadata)
        raise RuntimeError(f"Stage {stage_name} failed: {error}")
    
    def get_stage_status(self, stage_name: str) -> str:
        """Get the status of a stage."""
        return self._stage_status.get(stage_name, "unknown")
    
    def is_stage_complete(self, stage_name: str) -> bool:
        """Check if a stage has been completed."""
        return self._stage_status.get(stage_name) == "completed"
    
    def set_table(self, name: str, table: Any) -> None:
        """Store a data table in the context."""
        self.tables[name] = table
    
    def get_table(self, name: str) -> Any:
        """Retrieve a data table from the context."""
        return self.tables.get(name)
    
    def has_table(self, name: str) -> bool:
        """Check if a table exists."""
        return name in self.tables
    
    def add_qc_finding(self, finding: Dict[str, Any]) -> None:
        """Add a QC finding."""
        self.qc_findings.append(finding)
    
    def add_exclusion(self, exclusion: Dict[str, Any]) -> None:
        """Add an exclusion record."""
        self.exclusion_log.append(exclusion)
    
    def get_run_time(self) -> str:
        """Get current timestamp."""
        return datetime.utcnow().isoformat() + "Z"
    
    def to_manifest_entry(self) -> Dict[str, Any]:
        """Create a manifest entry for the current run context."""
        return {
            "run_id": self.config.run_id,
            "status": "completed" if all(
                s == "completed" 
                for s in self._stage_status.values() 
                if s in self._stage_order
            ) else "in_progress",
            "stages": self._stage_status,
            "tables_created": list(self.tables.keys()),
            "qc_findings_count": len(self.qc_findings),
            "exclusions_count": len(self.exclusion_log),
        }


def create_run_context(config: Any, logger: Any) -> RunContext:
    """Create a new run context."""
    return RunContext(config=config, logger=logger)


__all__ = [
    "RunContext",
    "AnalysisStage",
    "create_run_context",
]
