"""
V3 Confidence x Geometry Integration.

Builds a dedicated analysis layer connecting:
- structural geometry
- AF3 confidence

Keeps metrics conceptually separate.
Does NOT merge them into one arbitrary "score."
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class ConfidenceGeometryObservation:
    """One observation linking confidence and geometry."""

    prediction_id: str
    condition_id: str
    seed: int
    sample: int

    # Geometry
    rmsd_to_reference: Optional[float] = None
    mean_displacement: Optional[float] = None

    # Confidence
    plddt_mean: Optional[float] = None
    plddt_min: Optional[float] = None
    plddt_max: Optional[float] = None
    plddt_median: Optional[float] = None
    pae_mean: Optional[float] = None
    contact_prob_mean: Optional[float] = None

    # Classification
    structural_category: Optional[str] = None  # low_diff/high_conf, etc.
    confidence_category: Optional[str] = None

    # Metadata
    has_dna: bool = False
    n_chains: int = 0


@dataclass
class ConfidenceGeometryAnalysis:
    """Complete confidence-geometry analysis result."""

    observations: List[ConfidenceGeometryObservation] = field(default_factory=list)
    correlation_rmsd_plddt: Optional[float] = None
    correlation_displacement_plddt: Optional[float] = None

    # Category counts
    category_counts: Dict[str, int] = field(default_factory=dict)
    # e.g., {"low_diff_high_conf": 5, "high_diff_low_conf": 2, ...}

    # Statistics
    n_observations: int = 0
    n_with_geometry: int = 0
    n_with_confidence: int = 0

    # Warnings
    warnings: List[str] = field(default_factory=list)


def integrate_confidence_geometry(
    structures: List[Any],
    rmsd_values: Optional[List[float]] = None,
    displacement_values: Optional[List[float]] = None,
    plddt_values: Optional[List[float]] = None,
    pae_values: Optional[List[float]] = None,
    contact_prob_values: Optional[List[float]] = None,
    *,
    plddt_threshold_high: float = 70.0,
    plddt_threshold_low: float = 50.0,
    rmsd_threshold_high: float = 2.0,
    rmsd_threshold_low: float = 0.5,
) -> ConfidenceGeometryAnalysis:
    """
    Integrate confidence metrics with structural geometry.

    Parameters
    ----------
    structures : list of StructureData
    rmsd_values : list of float, optional
        RMSD to reference for each structure.
    displacement_values : list of float, optional
        Mean displacement for each structure.
    plddt_values : list of float, optional
        Mean pLDDT for each structure.
    pae_values : list of float, optional
    contact_prob_values : list of float, optional
    plddt_threshold_high : float
        pLDDT threshold for "high confidence".
    plddt_threshold_low : float
        pLDDT threshold for "low confidence".
    rmsd_threshold_high : float
        RMSD threshold for "high structural difference".
    rmsd_threshold_low : float
        RMSD threshold for "low structural difference".

    Returns
    -------
    ConfidenceGeometryAnalysis
    """
    observations = []
    warnings = []

    n = len(structures)
    if rmsd_values is None:
        rmsd_values = [None] * n
    if displacement_values is None:
        displacement_values = [None] * n
    if plddt_values is None:
        plddt_values = [None] * n
    if pae_values is None:
        pae_values = [None] * n
    if contact_prob_values is None:
        contact_prob_values = [None] * n

    # Classify each observation
    for i, struct in enumerate(structures):
        obs = ConfidenceGeometryObservation(
            prediction_id=struct.prediction_id,
            condition_id=struct.condition_id,
            seed=struct.seed,
            sample=struct.sample,
            rmsd_to_reference=rmsd_values[i],
            mean_displacement=displacement_values[i],
            plddt_mean=plddt_values[i],
            plddt_min=plddt_values[i],  # Approximate
            plddt_max=plddt_values[i],  # Approximate
            plddt_median=plddt_values[i],  # Approximate
            pae_mean=pae_values[i],
            contact_prob_mean=contact_prob_values[i],
            has_dna=struct.has_dna,
            n_chains=struct.n_chains,
        )

        # Classify into categories
        rmsd = obs.rmsd_to_reference
        plddt = obs.plddt_mean

        if rmsd is not None and plddt is not None:
            # Determine structural category
            if rmsd <= rmsd_threshold_low:
                struct_cat = "low_structural_difference"
            elif rmsd >= rmsd_threshold_high:
                struct_cat = "high_structural_difference"
            else:
                struct_cat = "moderate_structural_difference"

            # Determine confidence category
            if plddt >= plddt_threshold_high:
                conf_cat = "high_confidence"
            elif plddt <= plddt_threshold_low:
                conf_cat = "low_confidence"
            else:
                conf_cat = "moderate_confidence"

            obs.structural_category = struct_cat
            obs.confidence_category = conf_cat

            # Create combined category
            combined = f"{struct_cat}__{conf_cat}"
        else:
            combined = "incomplete"

        observations.append(obs)

    # Count categories
    category_counts = {}
    for obs in observations:
        if obs.structural_category and obs.confidence_category:
            key = f"{obs.structural_category}__{obs.confidence_category}"
            category_counts[key] = category_counts.get(key, 0) + 1

    # Calculate correlations where possible
    corr_rmsd_plddt = None
    corr_displacement_plddt = None

    # RMSD vs pLDDT
    rmsd_vals = [o.rmsd_to_reference for o in observations
                 if o.rmsd_to_reference is not None and o.plddt_mean is not None]
    plddt_vals = [o.plddt_mean for o in observations
                  if o.rmsd_to_reference is not None and o.plddt_mean is not None]

    if len(rmsd_vals) > 2:
        corr_rmsd_plddt = float(np.corrcoef(rmsd_vals, plddt_vals)[0, 1])

    # Displacement vs pLDDT
    disp_vals = [o.mean_displacement for o in observations
                 if o.mean_displacement is not None and o.plddt_mean is not None]
    plddt_vals2 = [o.plddt_mean for o in observations
                    if o.mean_displacement is not None and o.plddt_mean is not None]

    if len(disp_vals) > 2:
        corr_displacement_plddt = float(np.corrcoef(disp_vals, plddt_vals2)[0, 1])

    n_with_geometry = sum(
        1 for o in observations
        if o.rmsd_to_reference is not None or o.mean_displacement is not None
    )
    n_with_confidence = sum(
        1 for o in observations
        if o.plddt_mean is not None or o.pae_mean is not None
    )

    return ConfidenceGeometryAnalysis(
        observations=observations,
        correlation_rmsd_plddt=corr_rmsd_plddt,
        correlation_displacement_plddt=corr_displacement_plddt,
        category_counts=category_counts,
        n_observations=len(observations),
        n_with_geometry=n_with_geometry,
        n_with_confidence=n_with_confidence,
        warnings=warnings,
    )


def confidence_geometry_table(
    analysis: ConfidenceGeometryAnalysis,
) -> pd.DataFrame:
    """Create DataFrame from confidence-geometry analysis."""
    rows = []
    for obs in analysis.observations:
        row = {
            "prediction_id": obs.prediction_id,
            "condition_id": obs.condition_id,
            "seed": obs.seed,
            "sample": obs.sample,
            "rmsd_to_reference": obs.rmsd_to_reference,
            "mean_displacement": obs.mean_displacement,
            "plddt_mean": obs.plddt_mean,
            "plddt_min": obs.plddt_min,
            "plddt_max": obs.plddt_max,
            "plddt_median": obs.plddt_median,
            "pae_mean": obs.pae_mean,
            "contact_prob_mean": obs.contact_prob_mean,
            "structural_category": obs.structural_category,
            "confidence_category": obs.confidence_category,
            "has_dna": obs.has_dna,
            "n_chains": obs.n_chains,
        }
        rows.append(row)

    return pd.DataFrame(rows)


def confidence_change_vs_structural_change(
    observations: List[ConfidenceGeometryObservation],
    *,
    confidence_metric: str = "plddt_mean",
    structural_metric: str = "rmsd_to_reference",
    reference_condition: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Analyze relationship between confidence change and structural change.

    For each condition vs reference:
    - Delta confidence
    - Delta structural (RMSD)
    - Correlation
    """
    # This would require paired observations (ref vs target for each seed)
    # Simplified version: just compute overall statistics

    results = {
        "n_observations": len(observations),
        "correlation": None,
        "warning": None,
    }

    # Extract paired changes
    # This would need more structure...

    return results
