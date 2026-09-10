"""
Structured logging utilities for AF3 Confidence Analysis Pipeline.

Creates run.log (human-readable) and run.jsonl (machine-readable) files.
"""

import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional


class JSONFormatter(logging.Formatter):
    """Format log records as JSON lines."""
    
    def format(self, record: logging.LogRecord) -> str:
        log_data = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }
        
        # Add extra fields if present
        if hasattr(record, "extra_fields"):
            log_data.update(record.extra_fields)
        
        return json.dumps(log_data, default=str)


class RunLogger:
    """
    Dual-output logger for run.log and run.jsonl.
    
    run.log: Human-readable log with timestamps
    run.jsonl: Machine-readable JSON lines for processing
    """
    
    def __init__(self, log_dir: Path, run_id: str, level: int = logging.INFO):
        self.log_dir = log_dir
        self.run_id = run_id
        self._setup_logging(level)
    
    def _setup_logging(self, level: int) -> None:
        """Configure logging handlers."""
        self.logger = logging.getLogger(f"af3_analysis.{self.run_id}")
        self.logger.setLevel(level)
        self.logger.propagate = False
        
        # Clear existing handlers
        self.logger.handlers = []
        
        # Human-readable log (run.log)
        log_path = self.log_dir / "run.log"
        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setLevel(level)
        file_format = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler.setFormatter(file_format)
        self.logger.addHandler(file_handler)
        
        # Machine-readable log (run.jsonl)
        jsonl_path = self.log_dir / "run.jsonl"
        jsonl_handler = logging.FileHandler(jsonl_path, mode="w")
        jsonl_handler.setLevel(level)
        jsonl_handler.setFormatter(JSONFormatter())
        self.logger.addHandler(jsonl_handler)
        
        # Also output to console
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(file_format)
        self.logger.addHandler(console_handler)
    
    def log(self, level: int, message: str, **extra: Any) -> None:
        """Log a message with optional extra fields."""
        extra_fields = getattr(self, "_extra_fields", {})
        extra_fields.update(extra)
        
        # Create custom record with extra fields
        record = self.logger.makeRecord(
            self.logger.name, level, "", 0, message, (), None
        )
        record.extra_fields = extra_fields
        
        self.logger.handle(record)
    
    def info(self, message: str, **extra: Any) -> None:
        self.log(logging.INFO, message, **extra)
    
    def warning(self, message: str, **extra: Any) -> None:
        self.log(logging.WARNING, message, **extra)
    
    def error(self, message: str, **extra: Any) -> None:
        self.log(logging.ERROR, message, **extra)
    
    def debug(self, message: str, **extra: Any) -> None:
        self.log(logging.DEBUG, message, **extra)
    
    def start_event(self, event: str, **metadata: Any) -> None:
        """Log a structured start event."""
        self.info(f"Starting: {event}", event_type="start", **metadata)
    
    def end_event(self, event: str, **metadata: Any) -> None:
        """Log a structured end event."""
        self.info(f"Completed: {event}", event_type="end", **metadata)
    
    def qc_event(self, rule: str, severity: str, affected_keys: list, **metadata: Any) -> None:
        """Log a QC finding event."""
        self.info(
            f"QC: {rule} [{severity}]",
            event_type="qc_finding",
            rule_id=rule,
            severity=severity,
            affected_keys=affected_keys,
            **metadata,
        )
    
    def exclusion_event(self, primary_key: str, rule: str, source_path: str, **metadata: Any) -> None:
        """Log an exclusion event."""
        self.info(
            f"Exclusion: {primary_key} - {rule}",
            event_type="exclusion",
            primary_key=primary_key,
            rule=rule,
            source_path=str(source_path),
            **metadata,
        )


def create_logger(log_dir: Path, run_id: str, level: int = logging.INFO) -> RunLogger:
    """Create a RunLogger instance."""
    return RunLogger(log_dir, run_id, level)


__all__ = [
    "RunLogger",
    "create_logger",
]
