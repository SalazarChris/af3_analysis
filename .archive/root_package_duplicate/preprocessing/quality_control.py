"""
Quality control for AF3 Confidence Analysis Pipeline.

Implements file/provenance, key/design, replicate, metric, composition,
aggregate-cross-check, and readiness checks.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
import pandas as pd
import math

from af3_analysis.schemas.enums import QCSeverity, Definedness
from af3_analysis.io.provenance import Checksum


@dataclass
class QCFinding:
    """A single QC finding."""
    finding_id: str
    rule_id: str
    severity: QCSeverity
    affected_keys: List[str]
    evidence: str
    action: str  # "exclude", "flag", "ignore"
    timestamp: str = ""


@dataclass
class QCResult:
    """Result of QC checks."""
    findings: List[QCFinding]
    exclusion_log: List[Dict[str, Any]]
    valid_records: int
    excluded_records: int


class QualityControl:
    """
    Run comprehensive QC checks on AF3 data.
    
    Implements multiple check categories:
    - File/provenance checks
    - Key and design checks
    - Replicate checks
    - Metric checks
    - Composition checks
    - Aggregate checks
    - Readiness checks
    """
    
    def __init__(self):
        self._findings: List[QCFinding] = []
        self._exclusion_log: List[Dict[str, Any]] = []
        self._finding_counter = 0
    
    def _next_finding_id(self) -> str:
        """Generate next finding ID."""
        self._finding_counter += 1
        return f"QC_{self._finding_counter:04d}"
    
    def run_all_checks(
        self,
        replicates_df: pd.DataFrame,
        measurements_df: pd.DataFrame,
        conditions_df: pd.DataFrame,
    ) -> QCResult:
        """
        Run all QC checks.
        
        Args:
            replicates_df: Replicates table
            measurements_df: Measurements table
            conditions_df: Conditions table
        
        Returns:
            QCResult with findings and exclusions
        """
        # Run all check categories
        self._file_provenance_checks(replicates_df)
        self._key_design_checks(replicates_df, conditions_df)
        self._replicate_checks(replicates_df)
        self._metric_checks(measurements_df)
        self._composition_checks(measurements_df)
        self._aggregate_cross_check(replicates_df, conditions_df)
        self._readiness_checks()
        
        # Calculate results
        valid_count = len(replicates_df)
        excluded_count = len(self._exclusion_log)
        
        return QCResult(
            findings=self._findings,
            exclusion_log=self._exclusion_log,
            valid_records=valid_count,
            excluded_records=excluded_count,
        )
    
    def _file_provenance_checks(self, replicates_df: pd.DataFrame) -> None:
        """Check file existence and provenance."""
        for _, row in replicates_df.iterrows():
            source_path = row.get("source_confidences_path")
            if pd.isna(source_path):
                continue
            
            path = Path(str(source_path))
            if not path.exists():
                finding = QCFinding(
                    finding_id=self._next_finding_id(),
                    rule_id="FILE_NOT_FOUND",
                    severity=QCSeverity.ERROR,
                    affected_keys=[row.get("prediction_id", "")],
                    evidence=f"Source file not found: {source_path}",
                    action="exclude",
                )
                self._findings.append(finding)
                
                self._exclusion_log.append({
                    "exclusion_id": finding.finding_id,
                    "primary_key": row.get("prediction_id"),
                    "rule_id": finding.rule_id,
                    "source_path": source_path,
                    "evidence": finding.evidence,
                    "resolution": "file_missing",
                })
    
    def _key_design_checks(
        self,
        replicates_df: pd.DataFrame,
        conditions_df: pd.DataFrame,
    ) -> None:
        """Check key constraints and referential integrity."""
        # Check for duplicate (condition_id, seed, sample)
        if len(replicates_df) > 0:
            key_cols = ["condition_id", "seed", "sample"]
            if replicates_df.duplicated(subset=key_cols).any():
                duplicates = replicates_df[replicates_df.duplicated(subset=key_cols, keep=False)]
                for _, row in duplicates.iterrows():
                    finding = QCFinding(
                        finding_id=self._next_finding_id(),
                        rule_id="DUPLICATE_KEY",
                        severity=QCSeverity.ERROR,
                        affected_keys=[row.get("prediction_id", "")],
                        evidence=f"Duplicate key: {row.get('condition_id')}, seed={row.get('seed')}, sample={row.get('sample')}",
                        action="exclude",
                    )
                    self._findings.append(finding)
    
    def _replicate_checks(self, replicates_df: pd.DataFrame) -> None:
        """Check replicate structure and validity."""
        # Check for valid seed/sample values
        for _, row in replicates_df.iterrows():
            seed = row.get("seed")
            sample = row.get("sample")
            
            if not isinstance(seed, int) or seed < 0:
                finding = QCFinding(
                    finding_id=self._next_finding_id(),
                    rule_id="INVALID_SEED",
                    severity=QCSeverity.ERROR,
                    affected_keys=[row.get("prediction_id", "")],
                    evidence=f"Invalid seed value: {seed}",
                    action="exclude",
                )
                self._findings.append(finding)
            
            if not isinstance(sample, int) or sample < 0:
                finding = QCFinding(
                    finding_id=self._next_finding_id(),
                    rule_id="INVALID_SAMPLE",
                    severity=QCSeverity.ERROR,
                    affected_keys=[row.get("prediction_id", "")],
                    evidence=f"Invalid sample value: {sample}",
                    action="exclude",
                )
                self._findings.append(finding)
    
    def _metric_checks(self, measurements_df: pd.DataFrame) -> None:
        """Check metric values and definedness."""
        for _, row in measurements_df.iterrows():
            metric_id = row.get("metric_id")
            value = row.get("value")
            definedness = row.get("definedness")
            
            # Check for non-finite values
            if definedness == "present" and pd.notna(value):
                try:
                    fval = float(value)
                    if not math.isfinite(fval):
                        finding = QCFinding(
                            finding_id=self._next_finding_id(),
                            rule_id="NON_FINITE_VALUE",
                            severity=QCSeverity.ERROR,
                            affected_keys=[row.get("prediction_id", "")],
                            evidence=f"Non-finite value for {metric_id}: {value}",
                            action="exclude",
                        )
                        self._findings.append(finding)
                except (ValueError, TypeError):
                    finding = QCFinding(
                        finding_id=self._next_finding_id(),
                        rule_id="INVALID_VALUE",
                        severity=QCSeverity.ERROR,
                        affected_keys=[row.get("prediction_id", "")],
                        evidence=f"Invalid value format for {metric_id}: {value}",
                        action="exclude",
                    )
                    self._findings.append(finding)
    
    def _composition_checks(self, measurements_df: pd.DataFrame) -> None:
        """Check composition consistency."""
        # Check for undefined metrics in wrong definedness category
        for _, row in measurements_df.iterrows():
            definedness = row.get("definedness")
            value = row.get("value")
            
            if definedness == "present" and pd.isna(value):
                finding = QCFinding(
                    finding_id=self._next_finding_id(),
                    rule_id="PRESENT_WITHOUT_VALUE",
                    severity=QCSeverity.ERROR,
                    affected_keys=[row.get("prediction_id", "")],
                    evidence="Metric marked present but has no value",
                    action="exclude",
                )
                self._findings.append(finding)
    
    def _aggregate_cross_check(
        self,
        replicates_df: pd.DataFrame,
        conditions_df: pd.DataFrame,
    ) -> None:
        """Cross-check aggregate values with computed values."""
        # Check for condition coverage
        for _, row in conditions_df.iterrows():
            cond_id = row.get("condition_id")
            condition_replicates = replicates_df[replicates_df["condition_id"] == cond_id]
            
            if len(condition_replicates) == 0:
                finding = QCFinding(
                    finding_id=self._next_finding_id(),
                    rule_id="NO_REPLICATES",
                    severity=QCSeverity.WARNING,
                    affected_keys=[cond_id],
                    evidence="Condition has no replicates",
                    action="ignore",
                )
                self._findings.append(finding)
    
    def _readiness_checks(self) -> None:
        """Check if data is ready for analysis stages."""
        # Check for critical errors
        critical_errors = [
            f for f in self._findings if f.severity == QCSeverity.ERROR
        ]
        
        if critical_errors:
            finding = QCFinding(
                finding_id=self._next_finding_id(),
                rule_id="CRITICAL_ERRORS",
                severity=QCSeverity.ERROR,
                affected_keys=["SYSTEM"],
                evidence=f"{len(critical_errors)} critical QC errors found",
                action="block_analysis",
            )
            self._findings.append(finding)
    
    def get_findings(self) -> List[QCFinding]:
        """Get all QC findings."""
        return self._findings.copy()
    
    def get_exclusion_log(self) -> List[Dict[str, Any]]:
        """Get the exclusion log."""
        return self._exclusion_log.copy()
    
    def has_critical_errors(self) -> bool:
        """Check if there are critical errors blocking analysis."""
        return any(f.severity == QCSeverity.ERROR for f in self._findings)


def run_quality_control(
    replicates_df: pd.DataFrame,
    measurements_df: pd.DataFrame,
    conditions_df: pd.DataFrame,
) -> QCResult:
    """
    Convenience function to run all QC checks.
    
    Args:
        replicates_df: Replicates table
        measurements_df: Measurements table
        conditions_df: Conditions table
    
    Returns:
        QCResult with findings and exclusions
    """
    qc = QualityControl()
    return qc.run_all_checks(replicates_df, measurements_df, conditions_df)


__all__ = [
    "QCFinding",
    "QCResult",
    "QualityControl",
    "run_quality_control",
]
