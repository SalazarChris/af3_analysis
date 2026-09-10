"""
Factor analysis utilities for AF3 Confidence Analysis Pipeline.

Provides marginal analysis and monotonicity flagging for factorial designs.
"""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def factor_marginals(
    df: pd.DataFrame,
    metric_col: str,
    factor_cols: List[str],
) -> Dict[str, pd.DataFrame]:
    """Compute marginal means for each factor level.

    Parameters
    ----------
    df : pd.DataFrame
        Long-format table with one row per replicate.
    metric_col : str
        Column containing the metric values.
    factor_cols : list[str]
        Columns representing experimental factors.

    Returns
    -------
    dict[str, DataFrame]
        Mapping factor name -> DataFrame with columns
        ``[level, mean, sd, n]``.
    """
    results: Dict[str, pd.DataFrame] = {}
    for factor in factor_cols:
        if factor not in df.columns:
            continue
        grouped = df.groupby(factor)[metric_col].agg(["mean", "std", "count"])
        grouped = grouped.reset_index()
        grouped.columns = [factor, "mean", "sd", "n"]
        results[factor] = grouped
    return results


def monotonicity_flag(
    marginals: pd.DataFrame,
    level_col: str = None,
    mean_col: str = "mean",
) -> bool:
    """Check whether marginal means are monotonically ordered.

    Parameters
    ----------
    marginals : pd.DataFrame
        Output from :func:`factor_marginals` for a single factor.
    level_col : str, optional
        Column containing factor levels.  If *None*, use the first
        column that is not *mean_col*.
    mean_col : str
        Column containing mean values.

    Returns
    -------
    bool
        *True* if means are monotonically increasing or decreasing
        across levels (sorted alphabetically / numerically).
    """
    if level_col is None:
        level_col = [c for c in marginals.columns if c != mean_col][0]

    sorted_df = marginals.sort_values(level_col)
    means = sorted_df[mean_col].values

    if len(means) < 2:
        return True

    diffs = np.diff(means)
    return bool(np.all(diffs >= 0) or np.all(diffs <= 0))
