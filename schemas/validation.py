"""
Validation utilities for AF3 Confidence Analysis Pipeline.

Provides validation functions for records and tables.
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional
import pandas as pd

from af3_analysis.schemas.tables import (
    TableSchema,
    ValidationError,
)
from af3_analysis.schemas.enums import Definedness, MetricScope


def validate_record(table_schema: TableSchema, record: Dict[str, Any]) -> bool:
    """
    Validate a single record against a table schema.
    
    Returns True if valid, raises ValidationError if invalid.
    """
    table_schema.validate_record(record)
    return True


def validate_table(table_schema: TableSchema, df: pd.DataFrame) -> bool:
    """
    Validate a DataFrame table against a table schema.
    
    Returns True if valid, raises ValidationError if invalid.
    """
    table_schema.validate_table(df)
    return True


def validate_definedness_consistency(
    definedness: str,
    value: Optional[float],
    metric_id: str,
) -> bool:
    """
    Validate consistency between definedness and value.
    
    Rules:
    - 'present' definedness requires a non-null value
    - Other definedness requires null value
    """
    definedness_values = [d.value for d in Definedness]
    if definedness not in definedness_values:
        raise ValidationError(
            f"Invalid definedness '{definedness}' for metric '{metric_id}'",
            {"definedness": definedness, "metric_id": metric_id},
        )
    
    if definedness == "present" and value is None:
        raise ValidationError(
            f"Value must be present when definedness is 'present' for metric '{metric_id}'",
            {"metric_id": metric_id, "definedness": definedness},
        )
    
    if definedness != "present" and value is not None:
        raise ValidationError(
            f"Value must be None when definedness is '{definedness}' for metric '{metric_id}'",
            {"metric_id": metric_id, "definedness": definedness, "value": value},
        )
    
    return True


def validate_scope_type(metric_id: str, scope_type: str) -> MetricScope:
    """
    Validate and convert scope_type string to MetricScope enum.
    """
    scope_values = [s.value for s in MetricScope]
    if scope_type not in scope_values:
        raise ValidationError(
            f"Invalid scope_type '{scope_type}' for metric '{metric_id}'",
            {"metric_id": metric_id, "scope_type": scope_type},
        )
    return MetricScope(scope_type)


def validate_primary_key(record: Dict[str, Any], primary_keys: List[str]) -> bool:
    """Validate that all primary key fields are present and non-null."""
    for pk in primary_keys:
        if pk not in record:
            raise ValidationError(f"Missing primary key field '{pk}'", record)
        if record[pk] is None:
            raise ValidationError(f"Primary key field '{pk}' is null", record)
    return True


def validate_foreign_key(
    record: Dict[str, Any],
    foreign_keys: Dict[str, str],
    reference_tables: Dict[str, pd.DataFrame],
) -> bool:
    """Validate that foreign keys reference valid records."""
    for fk_col, ref_table_col in foreign_keys.items():
        if fk_col not in record:
            continue  # Foreign keys are optional
        ref_table_name, ref_col = ref_table_col.split(".")
        if ref_table_name not in reference_tables:
            continue  # Skip if reference table not available
        ref_df = reference_tables[ref_table_name]
        ref_values = set(ref_df[ref_col].dropna())
        if record[fk_col] not in ref_values:
            raise ValidationError(
                f"Foreign key '{fk_col}'={record[fk_col]} not found in {ref_table_name}.{ref_col}",
                record,
            )
    return True


def validate_checksum(checksum: Optional[str], path: Optional[Path]) -> bool:
    """
    Validate checksum format.
    
    Expected format: SHA256 hex digest (64 characters).
    """
    if checksum is None:
        return True  # Optional field
    if not re.match(r"^[a-fA-F0-9]{64}$", checksum):
        raise ValidationError(
            f"Invalid checksum format: expected SHA256 hex digest (64 chars), got '{checksum[:20]}...'",
            {"checksum": checksum, "path": str(path) if path else None},
        )
    return True


def validate_path_exists(path: Optional[Path], field_name: str) -> bool:
    """Validate that a path exists (if provided)."""
    if path is None:
        return True  # Optional field
    if not path.exists():
        raise ValidationError(
            f"Path for '{field_name}' does not exist: {path}",
            {field_name: str(path)},
        )
    return True


def validate_numeric_range(
    value: Optional[float],
    metric_id: str,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
) -> bool:
    """Validate that a numeric value is within expected range."""
    if value is None:
        return True  # Optional field
    
    if min_val is not None and value < min_val:
        raise ValidationError(
            f"Value {value} for metric '{metric_id}' is below minimum {min_val}",
            {"metric_id": metric_id, "value": value, "min": min_val},
        )
    
    if max_val is not None and value > max_val:
        raise ValidationError(
            f"Value {value} for metric '{metric_id}' is above maximum {max_val}",
            {"metric_id": metric_id, "value": value, "max": max_val},
        )
    
    return True


def validate_finite_value(value: Optional[float], metric_id: str) -> bool:
    """Validate that a numeric value is finite (not NaN or inf)."""
    if value is None:
        return True  # Optional field
    
    import math
    if not math.isfinite(value):
        raise ValidationError(
            f"Value for metric '{metric_id}' must be finite, got {value}",
            {"metric_id": metric_id, "value": value},
        )
    return True


__all__ = [
    "validate_record",
    "validate_table",
    "validate_definedness_consistency",
    "validate_scope_type",
    "validate_primary_key",
    "validate_foreign_key",
    "validate_checksum",
    "validate_path_exists",
    "validate_numeric_range",
    "validate_finite_value",
]
