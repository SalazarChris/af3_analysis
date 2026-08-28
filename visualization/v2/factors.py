"""
V2 Factor Parsing
=================
Decompose raw condition names into PTM state, DNA presence, and
environment — without hardcoding biological specifics into the
visualization layer.

Strategy
--------
The condition names follow the pattern:

    pou_<ptm>[_dna][_NA<val>_HOH<val>_CL<val>]

We parse these generically using the experiment metadata when
available, and fall back to regex decomposition otherwise.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Regex patterns for environment parsing
# ---------------------------------------------------------------------------

_ENV_RE = re.compile(
    r"NA(\d+)_HOH(\d+)_CL(\d+)", re.IGNORECASE
)


# ---------------------------------------------------------------------------
# Core decomposition
# ---------------------------------------------------------------------------

def parse_condition_name(
    name: str,
    design: Any = None,
) -> Dict[str, Any]:
    """Parse a condition name into its factor components.

    Parameters
    ----------
    name : str
        Raw condition name, e.g. ``"pou_tpo101_dna_NA100_HOH1000_CL100"``.
    design : ExperimentDesign, optional
        If provided, metadata is used for labels.

    Returns
    -------
    dict with keys: ``ptm_state``, ``has_dna``, ``environment``, ``env_label``,
    ``ptm_label``, ``display_label``.
    """
    lower = name.lower()

    # --- PTM state detection (order matters: multi-PTM before single) ---
    if "tpo101" in lower and "sep102" in lower:
        ptm_state = "T101 + S102"
    elif "tpo101" in lower:
        ptm_state = "T101"
    elif "sep102" in lower:
        ptm_state = "S102"
    else:
        ptm_state = "Baseline"

    # --- DNA detection ---
    has_dna = "_dna" in lower or "_dna_" in lower

    # --- Environment detection ---
    env_match = _ENV_RE.search(name)
    if env_match:
        na, hoh, cl = env_match.groups()
        env_key = f"NA{na}_HOH{hoh}_CL{cl}"
        env_label = f"Na{na} / Hoh{hoh} / Cl{cl}"
    else:
        env_key = "baseline"
        env_label = "Baseline"

    # --- Display labels ---
    if design is not None and name in design.conditions:
        ptm_label = design.conditions[name].label
    else:
        ptm_label = _ptm_display(ptm_state, has_dna)

    return {
        "condition_name": name,
        "ptm_state": ptm_state,
        "has_dna": has_dna,
        "environment": env_key,
        "env_label": env_label,
        "ptm_label": ptm_label,
        "display_label": _build_display(ptm_state, has_dna, env_key),
    }


def _ptm_display(ptm_state: str, has_dna: bool) -> str:
    """Short human-readable label (no environment)."""
    parts: List[str] = []
    if ptm_state != "Baseline":
        parts.append(ptm_state)
    if has_dna:
        parts.append("+ DNA")
    return "Baseline " + " ".join(parts) if not parts else " ".join(parts)


def _build_display(ptm_state: str, has_dna: bool, env_key: str) -> str:
    """Full display label used for tick labels etc."""
    dna = " + DNA" if has_dna else ""
    env = "" if env_key == "baseline" else f" ({env_key})"
    return f"{ptm_state}{dna}{env}"


# ---------------------------------------------------------------------------
# DataFrame augmentation
# ---------------------------------------------------------------------------

def add_factor_columns(
    df: pd.DataFrame,
    design: Any = None,
) -> pd.DataFrame:
    """Add ``ptm_state``, ``has_dna``, ``environment``, ``env_label`` columns.

    Works on *seed_aggregated* or *descriptive_stats* DataFrames that
    contain a ``condition_name`` column.
    """
    if "condition_name" not in df.columns:
        raise ValueError("DataFrame must contain a 'condition_name' column.")

    parsed = df["condition_name"].apply(
        lambda c: pd.Series(parse_condition_name(c, design=design))
    )

    # Only add columns that don't already exist
    for col in ["ptm_state", "has_dna", "environment", "env_label"]:
        if col not in df.columns:
            df[col] = parsed[col].values

    return df


# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------

def get_unique_ptm_states(df: pd.DataFrame) -> List[str]:
    """Ordered unique PTM states present in the DataFrame."""
    from .config import PTM_ORDER
    present = set(df["ptm_state"].unique()) if "ptm_state" in df.columns else set()
    return [p for p in PTM_ORDER if p in present]


def get_unique_environments(df: pd.DataFrame) -> List[str]:
    """Sorted unique environments present."""
    if "environment" not in df.columns:
        return []
    return sorted(df["environment"].unique())


def get_n_env_for_ptm(df: pd.DataFrame, ptm_state: str) -> Dict[str, int]:
    """Count how many environments exist for each PTM/DNA combination."""
    subset = df[df["ptm_state"] == ptm_state] if "ptm_state" in df.columns else df
    if subset.empty:
        return {}
    return (
        subset.groupby(["ptm_state", "has_dna"])["environment"]
        .nunique()
        .to_dict()
    )
