"""
Definedness classification for AF3 Confidence Analysis Pipeline.

Classifies missingness as present, undefined_by_composition, missing_technical, or not_collected.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
import pandas as pd

from af3_analysis.schemas.enums import Definedness


class DefinednessClassifier:
    """
    Classify the definedness of metrics.
    
    Categories:
    - present: Valid numeric value exists
    - undefined_by_composition: Cannot exist (e.g., interface metric in monomer)
    - missing_technical: Should exist but extraction failed
    - not_collected: Phase 1 schema did not extract it
    """
    
    def __init__(self):
        self._definedness_matrix: Dict[str, Dict[str, str]] = {}
    
    def classify_metric(
        self,
        metric_id: str,
        value: Optional[float],
        condition_composition: Optional[Dict[str, Any]] = None,
        phase1_schema: Optional[List[str]] = None,
    ) -> str:
        """
        Classify the definedness of a metric.
        
        Args:
            metric_id: The metric identifier
            value: The metric value (or None if missing)
            condition_composition: Condition-specific composition info
            phase1_schema: List of metrics collected in Phase 1
        
        Returns:
            Definedness category
        """
        # First check if undefined by composition (regardless of value)
        if condition_composition and self._is_undefined_by_composition(metric_id, condition_composition):
            return Definedness.UNDEFINED_BY_COMPOSITION.value
        
        # Check if value is present and finite
        if value is not None and self._is_finite(value):
            return Definedness.PRESENT.value
        
        # Check if not collected (Phase 1)
        if phase1_schema and metric_id not in phase1_schema:
            return Definedness.NOT_COLLECTED.value
        
        # Otherwise, missing_technical
        return Definedness.MISSING_TECHNICAL.value
    
    def _is_finite(self, value: float) -> bool:
        """Check if a value is finite."""
        import math
        return math.isfinite(value)
    
    def _is_undefined_by_composition(
        self,
        metric_id: str,
        condition_composition: Dict[str, Any],
    ) -> bool:
        """Check if a metric is undefined by composition."""
        # Interface metrics undefined for monomers
        interface_metrics = ["iptm", "chain_pair_iptm", "interface_pae", "interface_pde", "contact_prob"]
        
        if metric_id in interface_metrics and not condition_composition.get("is_complex", False):
            return True
        
        # Add more composition rules as needed
        return False
    
    def classify_dataframe(
        self,
        df: pd.DataFrame,
        condition_composition: Optional[Dict[str, Any]] = None,
        phase1_schema: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Classify definedness for all metrics in a DataFrame.
        
        Args:
            df: DataFrame with metric_id and value columns
            condition_composition: Condition composition info
            phase1_schema: Phase 1 schema metrics
        
        Returns:
            DataFrame with definedness column added
        """
        df = df.copy()
        df["definedness"] = df.apply(
            lambda row: self.classify_metric(
                row.get("metric_id"),
                row.get("value"),
                condition_composition,
                phase1_schema,
            ),
            axis=1,
        )
        return df
    
    def create_definedness_matrix(
        self,
        measurements_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Create a definedness matrix showing which metrics are defined for each prediction.
        
        Args:
            measurements_df: Measurements table
        
        Returns:
            DataFrame with predictions as rows and metrics as columns
        """
        if measurements_df.empty:
            return pd.DataFrame()
        
        # Pivot to wide format
        matrix = measurements_df.pivot_table(
            index="prediction_id",
            columns="metric_id",
            values="definedness",
            aggfunc="first",
        )
        
        return matrix


def classify_definedness(
    measurements_df: pd.DataFrame,
    condition_composition: Optional[Dict[str, Any]] = None,
    phase1_schema: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Convenience function to classify definedness.
    
    Args:
        measurements_df: Measurements table
        condition_composition: Condition composition info
        phase1_schema: Phase 1 schema metrics
    
    Returns:
        DataFrame with definedness column added
    """
    classifier = DefinednessClassifier()
    return classifier.classify_dataframe(measurements_df, condition_composition, phase1_schema)


__all__ = [
    "DefinednessClassifier",
    "classify_definedness",
]
