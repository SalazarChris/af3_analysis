"""
Eligibility gates for AF3 Confidence Analysis Pipeline.

Checks if data is eligible for various analysis stages.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True)
class EligibilityGates:
    """Eligibility gate results."""
    can_proceed: bool
    reasons: List[str]
    blocked_stages: List[str]
    eligible_metrics: List[str]
    eligible_conditions: List[str]


class EligibilityChecker:
    """
    Check if data is eligible for analysis stages.
    
    Blocks requested stages that lack:
    - Raw provenance
    - Required mappings
    - Compositional comparability
    - Eligible metric population
    """
    
    def __init__(self):
        self._checks: List[Dict[str, Any]] = []
    
    def check_stage_eligibility(
        self,
        stage: str,
        conditions_df: pd.DataFrame,
        replicates_df: pd.DataFrame,
        measurements_df: pd.DataFrame,
        resolved_registry: Dict[str, str],
        mappings_available: bool = False,
        raw_provenance_available: bool = True,
    ) -> EligibilityGates:
        """
        Check if data is eligible for a specific analysis stage.
        
        Args:
            stage: Stage to check (e.g., "local_analysis", "interface_analysis", "matrix_analysis")
            conditions_df: Conditions table
            replicates_df: Replicates table
            measurements_df: Measurements table
            resolved_registry: Dict mapping metric_id to resolution status
            mappings_available: Whether mappings are available
            raw_provenance_available: Whether raw provenance is available
        
        Returns:
            EligibilityGates result
        """
        reasons = []
        blocked_stages = []
        eligible_metrics = []
        eligible_conditions = list(conditions_df["condition_id"].unique()) if not conditions_df.empty else []
        
        # Check raw provenance
        if stage in ["canonical_analysis"] and not raw_provenance_available:
            reasons.append(f"Stage {stage} requires raw provenance")
            blocked_stages.append(stage)
        
        # Check mappings for structural stages
        if stage in ["local_analysis", "interface_analysis", "matrix_analysis"] and not mappings_available:
            reasons.append(f"Stage {stage} requires mappings")
            blocked_stages.append(stage)
        
        # Check compositional comparability
        if not self._has_compositional_comparability(conditions_df):
            reasons.append("Conditions lack compositional comparability")
            blocked_stages.append(stage)
        
        # Check eligible metric population
        if not self._has_eligible_metrics(resolved_registry):
            reasons.append("No eligible metrics for analysis")
            blocked_stages.append(stage)
        
        can_proceed = len(blocked_stages) == 0
        
        return EligibilityGates(
            can_proceed=can_proceed,
            reasons=reasons,
            blocked_stages=blocked_stages,
            eligible_metrics=eligible_metrics,
            eligible_conditions=eligible_conditions,
        )
    
    def _has_compositional_comparability(self, conditions_df: pd.DataFrame) -> bool:
        """Check if conditions have compositional comparability."""
        if conditions_df.empty:
            return False
        
        # Check for at least 2 conditions
        if len(conditions_df) < 2:
            return False
        
        # Check that conditions have comparable composition
        # (simplified check - in practice, would compare input signatures)
        return True
    
    def _has_eligible_metrics(self, resolved_registry: Dict[str, str]) -> bool:
        """Check if there are eligible metrics."""
        # Count available metrics
        available = sum(
            1 for status in resolved_registry.values()
            if status == "resolved_available"
        )
        
        return available > 0
    
    def get_analysis_eligibility(
        self,
        conditions_df: pd.DataFrame,
        replicates_df: pd.DataFrame,
        measurements_df: pd.DataFrame,
        resolved_registry: Dict[str, str],
    ) -> Dict[str, Any]:
        """
        Get overall analysis eligibility.
        
        Args:
            conditions_df: Conditions table
            replicates_df: Replicates table
            measurements_df: Measurements table
            resolved_registry: Dict mapping metric_id to resolution status
        
        Returns:
            Dictionary with eligibility info for each stage
        """
        stages = [
            "global_analysis",
            "local_analysis",
            "interface_analysis",
            "matrix_analysis",
            "factorial_analysis",
        ]
        
        results = {}
        for stage in stages:
            gates = self.check_stage_eligibility(stage, conditions_df, replicates_df, measurements_df, resolved_registry)
            results[stage] = {
                "can_proceed": gates.can_proceed,
                "reasons": gates.reasons,
                "blocked_stages": gates.blocked_stages,
            }
        
        return results


def check_stage_eligibility(
    stage: str,
    conditions_df: pd.DataFrame,
    replicates_df: pd.DataFrame,
    measurements_df: pd.DataFrame,
    resolved_registry: Dict[str, str],
    mappings_available: bool = False,
    raw_provenance_available: bool = True,
) -> EligibilityGates:
    """
    Convenience function to check stage eligibility.
    
    Args:
        stage: Stage to check
        conditions_df: Conditions table
        replicates_df: Replicates table
        measurements_df: Measurements table
        resolved_registry: Dict mapping metric_id to resolution status
        mappings_available: Whether mappings are available
        raw_provenance_available: Whether raw provenance is available
    
    Returns:
        EligibilityGates result
    """
    checker = EligibilityChecker()
    return checker.check_stage_eligibility(
        stage, conditions_df, replicates_df, measurements_df, resolved_registry,
        mappings_available, raw_provenance_available
    )


def get_analysis_eligibility(
    conditions_df: pd.DataFrame,
    replicates_df: pd.DataFrame,
    measurements_df: pd.DataFrame,
    resolved_registry: Dict[str, str],
) -> Dict[str, Any]:
    """
    Convenience function to get overall analysis eligibility.
    
    Args:
        conditions_df: Conditions table
        replicates_df: Replicates table
        measurements_df: Measurements table
        resolved_registry: Dict mapping metric_id to resolution status
    
    Returns:
        Dictionary with eligibility info for each stage
    """
    checker = EligibilityChecker()
    return checker.get_analysis_eligibility(
        conditions_df, replicates_df, measurements_df, resolved_registry
    )


__all__ = [
    "EligibilityGates",
    "EligibilityChecker",
    "check_stage_eligibility",
    "get_analysis_eligibility",
]
