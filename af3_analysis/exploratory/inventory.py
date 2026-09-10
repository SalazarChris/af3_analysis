"""Design inventory and exploratory analysis for AF3 Confidence Analysis Pipeline."""

import pandas as pd
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass


@dataclass
class DesignInventory:
    """Summary of experimental design."""
    condition_count: int
    factor_levels: Dict[str, int]
    replicates_per_condition: Dict[str, int]
    is_balanced: bool
    seed_sets: Dict[str, set]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "condition_count": self.condition_count,
            "factor_levels": self.factor_levels,
            "replicates_per_condition": self.replicates_per_condition,
            "is_balanced": self.is_balanced,
            "seed_sets": {k: list(v) for k, v in self.seed_sets.items()},
        }


def is_fully_crossed(factors_df: pd.DataFrame) -> Tuple[bool, List[Dict[str, Any]]]:
    """
    Check if the factorial design is fully crossed.
    
    Args:
        factors_df: Wide factor matrix with conditions as rows
    
    Returns:
        Tuple of (is_fully_crossed, missing_cells)
    """
    if factors_df.empty:
        return True, []
    
    # Get factor columns (exclude n_replicates and other metadata)
    factor_cols = [c for c in factors_df.columns if c.startswith("factor_")]
    if not factor_cols:
        return True, []
    
    # Count unique level combinations
    unique_combinations = factors_df[factor_cols].drop_duplicates()
    expected_count = 1
    for col in factor_cols:
        expected_count *= factors_df[col].nunique()
    
    missing = []
    if unique_combinations.shape[0] < expected_count:
        # Find missing cells (would require full factorial generation)
        missing.append({
            "type": "missing_combinations",
            "missing_count": expected_count - unique_combinations.shape[0],
        })
    
    return len(missing) == 0, missing


def inventory_markdown(
    inv: DesignInventory,
    constants_report: Optional[pd.DataFrame] = None,
    collinear_pairs: Optional[List[Tuple[str, str, float]]] = None,
) -> str:
    """
    Create human-readable design inventory markdown.
    
    Args:
        inv: DesignInventory
        constants_report: DataFrame of constant columns
        collinear_pairs: List of (col1, col2, correlation) tuples
    
    Returns:
        Markdown string
    """
    lines = ["# Design Inventory", ""]
    
    lines.extend([
        f"**Condition Count:** {inv.condition_count}",
        "",
        "## Factor Levels",
        "",
    ])
    
    for factor, levels in inv.factor_levels.items():
        lines.append(f"- {factor}: {levels} levels")
    
    lines.extend([
        "",
        "## Replicates per Condition",
        "",
    ])
    
    for cond, count in inv.replicates_per_condition.items():
        lines.append(f"- {cond}: {count} replicates")
    
    lines.extend([
        "",
        f"**Balanced Design:** {'Yes' if inv.is_balanced else 'No'}",
        "",
    ])
    
    if constants_report is not None and not constants_report.empty:
        lines.extend([
            "## Constant Factors",
            "",
        ])
        for _, row in constants_report.iterrows():
            lines.append(f"- {row['factor']}: {row['value']}")
        lines.append("")
    
    if collinear_pairs is not None and len(collinear_pairs) > 0:
        lines.extend([
            "## Collinear Factor Pairs",
            "",
        ])
        for col1, col2, corr in collinear_pairs:
            lines.append(f"- {col1} ↔ {col2}: r = {corr:.3f}")
        lines.append("")
    
    return "\n".join(lines)


__all__ = [
    "DesignInventory",
    "is_fully_crossed",
    "inventory_markdown",
]
