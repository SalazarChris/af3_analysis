"""
Wide-to-long adapter for AF3 seed-aggregated data.

Converts the wide-format ``seed_aggregated.csv`` produced by the analysis
pipeline into the long-format DataFrame expected by the ``statistics/`` and
``exploratory/`` subpackages.

Design principles
-----------------
* **Pure adapter** — no statistical calculations, no method choices.
* **Non-destructive** — never mutates the input DataFrame.
* **Explicit** — every metric is classified with a ``scope_type`` and
  ``scope_id``; non-metric columns are never silently included.
* **Reversible** — the long format can be pivoted back to the original
  wide format without loss.

Expected downstream schema
---------------------------
The long-format DataFrame returned by :func:`seed_aggregated_to_long`
contains exactly these columns:

.. code-block::

    condition_id   str   – condition identifier (e.g. ``cond_001``)
    condition_name str   – human-readable condition name
    seed           int   – seed identifier
    metric_id      str   – metric identifier (e.g. ``pLDDT_mean``)
    scope_type     str   – ``"global"`` or ``"chain"``
    scope_id       str   – empty string for global, chain letter for chain
    value          float – metric value (may be NaN)

Column conventions
------------------
* **Global metrics** have ``scope_type="global"`` and ``scope_id=""``.
* **Chain-level metrics** have ``scope_type="chain"`` and ``scope_id``
  set to the chain letter (e.g. ``"A"``).
* Columns whose name matches ``chain_<X>_residues`` are excluded because
  they represent fixed sequence lengths, not measured confidence values.
* The ``condition_name`` column is carried through for downstream
  labeling but is *not* treated as a metric.
"""

from __future__ import annotations

import re
from typing import List, Optional, Sequence

import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_ID_COLUMNS: List[str] = ["condition_id", "condition_name", "seed"]

# Columns that are structural metadata, not metrics.
_EXCLUDED_COLUMNS: List[str] = ["condition_id", "condition_name", "seed"]

# Pattern for chain-level metric columns: chain_<X>_<metric>
_CHAIN_METRIC_RE = re.compile(r"^chain_([A-Z])_(.+)$")

# Pattern for chain residue-count columns to exclude.
_CHAIN_RESIDUES_RE = re.compile(r"^chain_[A-Z]_residues$")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def seed_aggregated_to_long(
    df: pd.DataFrame,
    *,
    include_condition_name: bool = True,
    extra_id_columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Convert a wide-format seed-aggregated DataFrame to long format.

    Parameters
    ----------
    df : pandas.DataFrame
        Wide-format DataFrame as produced by the analysis pipeline's
        ``_stage_run_analysis``.  Expected columns include at least
        ``condition_id``, ``condition_name``, and ``seed``, plus one or
        more numeric metric columns.
    include_condition_name : bool, optional
        If *True* (default), carry ``condition_name`` through to the
        long-format output.  Set to *False* if only the identifier is
        needed.
    extra_id_columns : sequence of str, optional
        Additional columns to treat as identifier (non-metric) columns.
        These are preserved in the output but not melted into the
        ``metric_id`` dimension.

    Returns
    -------
    pandas.DataFrame
        Long-format DataFrame with columns:
        ``condition_id``, ``condition_name`` (if requested), ``seed``,
        ``metric_id``, ``scope_type``, ``scope_id``, ``value``.
    """
    id_cols = list(_ID_COLUMNS)
    if not include_condition_name and "condition_name" in id_cols:
        id_cols.remove("condition_name")
    if extra_id_columns:
        id_cols.extend(c for c in extra_id_columns if c not in id_cols)

    # Validate required columns
    missing = [c for c in _ID_COLUMNS[:2] + ["seed"] if c not in df.columns]  # condition_id, condition_name, seed
    if missing:
        raise ValueError(f"Input DataFrame is missing required columns: {missing}")

    # Classify each non-id column
    global_metrics: List[str] = []
    chain_metrics: List[str] = []  # (metric_name, chain_letter) stored temporarily

    for col in df.columns:
        if col in id_cols:
            continue
        if not pd.api.types.is_numeric_dtype(df[col]):
            continue
        # Exclude chain residue-count columns (fixed sequence lengths)
        if _CHAIN_RESIDUES_RE.match(col):
            continue

        m = _CHAIN_METRIC_RE.match(col)
        if m:
            chain_letter = m.group(1)
            metric_name = m.group(2)
            chain_metrics.append((metric_name, chain_letter, col))
        else:
            global_metrics.append(col)

    # Build long-format rows
    rows: List[dict] = []

    # -- Global metrics --
    for metric_col in global_metrics:
        for _, row in df.iterrows():
            entry: dict = {
                "condition_id": row["condition_id"],
                "seed": int(row["seed"]),
                "metric_id": metric_col,
                "scope_type": "global",
                "scope_id": "",
                "value": row[metric_col] if pd.notna(row[metric_col]) else None,
            }
            if "condition_name" in id_cols:
                entry["condition_name"] = row["condition_name"]
            rows.append(entry)

    # -- Chain-level metrics --
    for metric_name, chain_letter, col in chain_metrics:
        for _, row in df.iterrows():
            entry = {
                "condition_id": row["condition_id"],
                "seed": int(row["seed"]),
                "metric_id": metric_name,
                "scope_type": "chain",
                "scope_id": chain_letter,
                "value": row[col] if pd.notna(row[col]) else None,
            }
            if "condition_name" in id_cols:
                entry["condition_name"] = row["condition_name"]
            rows.append(entry)

    long_df = pd.DataFrame(rows, columns=(
        [c for c in ["condition_id", "condition_name", "seed"] if c in id_cols]
        + ["metric_id", "scope_type", "scope_id", "value"]
    ))

    return long_df


def long_to_seed_aggregated(df_long: pd.DataFrame) -> pd.DataFrame:
    """Pivot a long-format DataFrame back to wide seed-aggregated form.

    This is the inverse of :func:`seed_aggregated_to_long`.  It is
    provided for verification and round-trip testing.

    Parameters
    ----------
    df_long : pandas.DataFrame
        Long-format DataFrame with columns ``condition_id``,
        ``condition_name`` (optional), ``seed``, ``metric_id``,
        ``scope_type``, ``scope_id``, ``value``.

    Returns
    -------
    pandas.DataFrame
        Wide-format DataFrame with one row per (condition_id, seed).
    """
    # Reconstruct the original column name from metric_id + scope
    df = df_long.copy()
    df["wide_column"] = df.apply(_reconstruct_wide_column, axis=1)

    id_cols = ["condition_id", "seed"]
    if "condition_name" in df.columns:
        id_cols.insert(1, "condition_name")

    wide = df.pivot_table(
        index=id_cols,
        columns="wide_column",
        values="value",
        aggfunc="first",
    ).reset_index()

    wide.columns.name = None
    return wide


def _reconstruct_wide_column(row: pd.Series) -> str:
    """Reconstruct the original wide-format column name from a long row."""
    if row["scope_type"] == "chain" and row["scope_id"]:
        return f"chain_{row['scope_id']}_{row['metric_id']}"
    return row["metric_id"]


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------


def load_seed_aggregated_long(
    path: str,
    *,
    include_condition_name: bool = True,
) -> pd.DataFrame:
    """Load ``seed_aggregated.csv`` and return it in long format.

    Parameters
    ----------
    path : str
        Path to ``seed_aggregated.csv`` (or any wide-format CSV with
        the expected schema).
    include_condition_name : bool, optional
        Whether to include ``condition_name`` in the output (default *True*).

    Returns
    -------
    pandas.DataFrame
        Long-format DataFrame.
    """
    df = pd.read_csv(path)
    return seed_aggregated_to_long(df, include_condition_name=include_condition_name)
