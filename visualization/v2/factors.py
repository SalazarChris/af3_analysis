"""
V2 Factor Parsing
=================
Decompose raw condition names into factor components (state,
ligand presence, environment) without hardcoding biological
specifics.

Strategy
--------
1. When experiment metadata (``ExperimentDesign``) is available,
   factor values are read directly from the authoritative attribute
   definitions — no string parsing needed.
2. When no metadata is available, a generic regex decomposition
   extracts environment patterns and known ligand tokens; the
   remainder of the condition name is treated as the "state" factor.
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

# Known ligand tokens (case-insensitive) — extend as needed
_LIGAND_TOKENS: Set[str] = {"dna"}


def parse_condition_name(
    name: str,
    design: Any = None,
) -> Dict[str, Any]:
    """Parse a condition name into its factor components.

    When *design* is provided, factor values come directly from the
    experiment metadata.  Otherwise a generic decomposition is used.

    Parameters
    ----------
    name : str
        Raw condition name.
    design : ExperimentDesign, optional
        Experiment metadata for authoritative factor definitions.

    Returns
    -------
    dict with keys: ``ptm_state``, ``has_dna``, ``environment``,
    ``env_label``, ``ptm_label``, ``display_label``.
    """
    lower = name.lower()

    # ------------------------------------------------------------------
    # Path A: metadata-driven decomposition (preferred)
    # ------------------------------------------------------------------
    if design is not None and name in design.conditions:
        cond_def = design.conditions[name]
        attrs = cond_def.attributes

        # Derive ptm_state from the first categorical attribute that
        # is NOT a ligand and NOT an environment attribute.
        ptm_state = None
        has_dna = False
        for attr_name, attr_def in design.attributes.items():
            val = attrs.get(attr_name)
            if val is None:
                continue
            if attr_def.attr_type == "binary":
                # Binary attrs are treated as ligand flags
                if val is True and attr_name.lower() in _LIGAND_TOKENS:
                    has_dna = True
                elif val is True:
                    # Generic binary attribute — include in state
                    label = attr_def.name if attr_def.name else attr_name
                    ptm_state = label if ptm_state is None else f"{ptm_state} + {label}"
            else:
                # Categorical attrs define the state
                label = str(val) if val != "baseline" and val != "Baseline" else "Baseline"
                ptm_state = label if ptm_state is None else f"{ptm_state} + {label}"

        if ptm_state is None:
            ptm_state = "Baseline"

        # Environment from metadata — check for NA/HOH/CL attributes
        env_key, env_label = _extract_env_from_name(name)

        return {
            "condition_name": name,
            "ptm_state": ptm_state,
            "has_dna": has_dna,
            "environment": env_key,
            "env_label": env_label,
            "ptm_label": cond_def.label,
            "display_label": _build_display(ptm_state, has_dna, env_key),
        }

    # ------------------------------------------------------------------
    # Path B: generic regex decomposition (no metadata)
    # ------------------------------------------------------------------
    tokens = lower.split("_")

    # Detect ligand tokens (e.g. _dna)
    ligand_hits = [t for t in tokens if t in _LIGAND_TOKENS]
    has_dna = "dna" in ligand_hits

    # Detect environment suffix
    env_key, env_label = _extract_env_from_name(name)
    env_match = _ENV_RE.search(name)
    env_suffix_start = env_match.start() if env_match else None

    # The "state" is everything that is not a ligand token and not
    # the environment suffix.  We also skip the first token (assumed
    # to be the protein/system name).
    state_tokens: List[str] = []
    for tok in tokens:
        if tok in _LIGAND_TOKENS:
            continue
        if env_suffix_start is not None and name.lower().find(tok) >= env_suffix_start:
            continue
        state_tokens.append(tok)
    # Drop the first token (protein prefix)
    if state_tokens:
        state_tokens = state_tokens[1:]

    ptm_state = " ".join(state_tokens).strip().title() if state_tokens else "Baseline"

    # --- Display labels ---
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


def _extract_env_from_name(name: str) -> Tuple[str, str]:
    """Extract environment key and label from a condition name."""
    env_match = _ENV_RE.search(name)
    if env_match:
        na, hoh, cl = env_match.groups()
        return f"NA{na}_HOH{hoh}_CL{cl}", f"Na{na} / Hoh{hoh} / Cl{cl}"
    return "baseline", "Baseline"


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
    from .config import derive_ptm_order
    if "ptm_state" not in df.columns:
        return []
    return derive_ptm_order(df["ptm_state"].unique().tolist())


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
