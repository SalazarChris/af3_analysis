"""
V3 Pairwise Structural Distances.

Calculates a pairwise structural distance matrix.
Uses pairwise RMSD.

Preserves metadata for every prediction:
- condition
- seed
- prediction ID

Caches this calculation because it can be expensive.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..model import StructureData
from .comparability import (
    find_common_space_prepared,
    kabsch_rmsd_prepared,
    prepare_structure,
    PreparedStructure,
)
from .rmsd import calculate_rmsd

logger = logging.getLogger(__name__)

# Progress log interval (pairs). Keeps long O(N^2) runs visibly alive
# without flooding logs for small datasets.
PROGRESS_LOG_EVERY = 500

# Minimum pair count before progress messages are raised to WARNING so they
# remain visible even if only warnings are being shown (e.g. quiet console).
PROGRESS_LOG_HEAVY_RUN_MIN_PAIRS = 20_000


@dataclass
class PairwiseDistanceMatrix:
    """Pairwise structural distance matrix with metadata."""

    # Distance matrix
    matrix: np.ndarray  # (N, N) array of RMSD values
    n_structures: int

    # Metadata for each structure
    predictions: List[str]  # prediction_id for each row/col
    conditions: List[str]  # condition_id for each row/col
    seeds: List[int]  # seed for each row/col
    samples: List[int]  # sample for each row/col

    # Validity mask (where RMSD could be calculated)
    valid: np.ndarray  # (N, N) boolean mask

    # Summary statistics
    mean_distance: float
    median_distance: float
    min_distance: float
    max_distance: float

    # Status
    status: str = "complete"
    reason: str = ""


def calculate_pairwise_rmsd_matrix(
    structures: List[StructureData],
    *,
    alignment_atom: str = "CA",
    min_common_atoms: int = 3,
    min_sequence_identity: float = 0.5,
    min_coverage: float = 0.80,
    cache_dir: Optional[Path] = None,
    cache_key: Optional[str] = None,
) -> PairwiseDistanceMatrix:
    """
    Calculate pairwise RMSD matrix for all structures.

    Parameters
    ----------
    structures : list of StructureData
    alignment_atom : str
    min_common_atoms : int
    min_sequence_identity : float
    min_coverage : float
    cache_dir : Path, optional
        Directory to cache the matrix.
    cache_key : str, optional
        Cache key (e.g., hash of configuration).

    Returns
    -------
    PairwiseDistanceMatrix
    """
    n = len(structures)

    if n == 0:
        return PairwiseDistanceMatrix(
            matrix=np.array([]).reshape(0, 0),
            n_structures=0,
            predictions=[],
            conditions=[],
            seeds=[],
            samples=[],
            valid=np.array([]).reshape(0, 0),
            mean_distance=np.nan,
            median_distance=np.nan,
            min_distance=np.nan,
            max_distance=np.nan,
        )

    # Extract metadata
    prediction_ids = [s.prediction_id for s in structures]
    conditions = [s.condition_id for s in structures]
    seeds = [s.seed for s in structures]
    samples = [s.sample for s in structures]

    # Initialize matrix
    matrix = np.full((n, n), np.nan)
    valid = np.zeros((n, n), dtype=bool)

    # Prepare per-structure inventories once (O(N) setup) so each of the
    # N*(N-1)/2 comparisons resolves the common space with set operations
    # instead of rebuilding residue dictionaries per pair. The prepared path
    # is only used for CA alignment; other alignment atoms fall back to the
    # general per-pair implementation (same results, validated by tests).
    prepared: List[PreparedStructure] = [
        prepare_structure(s) for s in structures
    ]
    use_fast_path = alignment_atom == "CA"

    n_pairs = n * (n - 1) // 2
    t_start = time.perf_counter()
    logger.info(
        "[V3 Pairwise] Computing %d pairwise RMSD values for %d structures",
        n_pairs, n,
    )

    # Calculate pairwise RMSD
    # Only calculate upper triangle, mirror to lower
    n_done = 0
    progress_level = (
        logging.WARNING if n_pairs >= PROGRESS_LOG_HEAVY_RUN_MIN_PAIRS else logging.INFO
    )
    for i in range(n):
        for j in range(i + 1, n):
            if use_fast_path:
                status, _coverage, common_residues, n_ca_pairs = (
                    find_common_space_prepared(
                        prepared[i],
                        prepared[j],
                        min_sequence_identity=min_sequence_identity,
                    )
                )
                if (
                    status != "not_comparable"
                    and n_ca_pairs >= min_common_atoms
                ):
                    rmsd_value = kabsch_rmsd_prepared(
                        prepared[i], prepared[j], common_residues
                    )
                else:
                    rmsd_value = None
            else:
                result = calculate_rmsd(
                    structures[i],
                    structures[j],
                    alignment_atom=alignment_atom,
                    min_common_atoms=min_common_atoms,
                    min_sequence_identity=min_sequence_identity,
                    min_coverage=min_coverage,
                )
                rmsd_value = result["rmsd"]

            if rmsd_value is not None:
                matrix[i, j] = rmsd_value
                matrix[j, i] = rmsd_value
                valid[i, j] = True
                valid[j, i] = True

            n_done += 1
            if n_done % PROGRESS_LOG_EVERY == 0:
                elapsed = time.perf_counter() - t_start
                rate = n_done / elapsed if elapsed > 0 else 0.0
                eta_s = (n_pairs - n_done) / rate if rate > 0 else float("nan")
                logger.log(
                    progress_level,
                    "[V3 Pairwise] %d/%d pairs (%.0f%%), %.1f pairs/s, ETA %.0f s",
                    n_done, n_pairs, 100.0 * n_done / n_pairs, rate, eta_s,
                )

    elapsed = time.perf_counter() - t_start
    logger.info("[V3 Pairwise] Done: %d pairs in %.1f s", n_pairs, elapsed)

    # Set diagonal to 0 (self-comparison)
    np.fill_diagonal(matrix, 0.0)
    np.fill_diagonal(valid, True)

    # Calculate summary statistics
    valid_distances = matrix[valid]
    valid_distances = valid_distances[valid_distances > 0]  # Exclude self

    if len(valid_distances) > 0:
        mean_dist = float(np.mean(valid_distances))
        median_dist = float(np.median(valid_distances))
        min_dist = float(np.min(valid_distances))
        max_dist = float(np.max(valid_distances))
    else:
        mean_dist = median_dist = min_dist = max_dist = np.nan

    return PairwiseDistanceMatrix(
        matrix=matrix,
        n_structures=n,
        predictions=prediction_ids,
        conditions=conditions,
        seeds=seeds,
        samples=samples,
        valid=valid,
        mean_distance=mean_dist,
        median_distance=median_dist,
        min_distance=min_dist,
        max_distance=max_dist,
    )


def save_pairwise_matrix(
    matrix: PairwiseDistanceMatrix,
    output_dir: Path,
    filename: str = "pairwise_rmsd_matrix.csv",
) -> Path:
    """Save pairwise matrix to CSV with metadata."""
    # Create DataFrame with prediction IDs as index/columns
    df = pd.DataFrame(
        matrix.matrix,
        index=matrix.predictions,
        columns=matrix.predictions,
    )

    # Add metadata columns
    metadata = pd.DataFrame({
        "prediction_id": matrix.predictions,
        "condition": matrix.conditions,
        "seed": matrix.seeds,
        "sample": matrix.samples,
    }, index=matrix.predictions)

    # Save matrix
    df.to_csv(output_dir / filename)

    # Save metadata next to the matrix (name derived from the matrix
    # filename so loader/save stay in sync when a cache key is used).
    matrix_path = output_dir / filename
    metadata_path = matrix_path.with_name(matrix_path.stem + "_metadata.csv")
    metadata.to_csv(metadata_path)

    return matrix_path


def load_pairwise_matrix(
    matrix_path: Path,
    metadata_path: Optional[Path] = None,
) -> PairwiseDistanceMatrix:
    """Load pairwise matrix from CSV."""
    df = pd.read_csv(matrix_path, index_col=0)

    predictions = list(df.index)
    matrix_array = df.values

    if metadata_path is not None:
        metadata = pd.read_csv(metadata_path, index_col=0)
        conditions = metadata.get("condition", [""] * len(predictions)).tolist()
        seeds = metadata.get("seed", [0] * len(predictions)).tolist()
        samples = metadata.get("sample", [0] * len(predictions)).tolist()
    else:
        conditions = [""] * len(predictions)
        seeds = [0] * len(predictions)
        samples = [0] * len(predictions)

    # Reconstruct validity mask: NaN cells were never computable; the
    # diagonal (self-comparison) is always valid. A zero RMSD between two
    # distinct structures is a legitimate measurement and must stay valid.
    valid_array = ~np.isnan(matrix_array)
    np.fill_diagonal(valid_array, True)

    # Statistics follow the compute-side convention: exclude self-comparisons
    # and zero/negative values.
    n_pred = len(predictions)
    off_diagonal = ~np.eye(n_pred, dtype=bool)
    distances = matrix_array[valid_array & off_diagonal]
    distances = distances[distances > 0]
    if len(distances) > 0:
        mean_dist = float(np.mean(distances))
        median_dist = float(np.median(distances))
        min_dist = float(np.min(distances))
        max_dist = float(np.max(distances))
    else:
        mean_dist = median_dist = min_dist = max_dist = np.nan

    return PairwiseDistanceMatrix(
        matrix=matrix_array,
        n_structures=len(predictions),
        predictions=predictions,
        conditions=conditions,
        seeds=seeds,
        samples=samples,
        valid=valid_array,
        mean_distance=mean_dist,
        median_distance=median_dist,
        min_distance=min_dist,
        max_distance=max_dist,
    )
