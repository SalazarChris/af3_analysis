"""
V2 Label Utilities
==================
Human-readable labels for metrics, axes, legends, and figure captions.
"""

from __future__ import annotations

from typing import Dict, List


# Metric column → axis label
METRIC_LABELS: Dict[str, str] = {
    "pLDDT_mean":          "Mean pLDDT",
    "pLDDT_min":           "Minimum pLDDT",
    "pae_mean":            "Mean PAE (Å)",
    "contact_prob_mean":   "Mean Contact Probability",
    "ranking_score":       "Ranking Score",
}

# Metric column → short panel letter
METRIC_PANEL_LETTERS: Dict[str, str] = {
    "pLDDT_mean":        "A",
    "pLDDT_min":         "B",
    "pae_mean":          "C",
    "contact_prob_mean": "D",
    "ranking_score":     "E",
}

# Metrics ordered as they appear in Figure 7
FIGURE7_METRICS: List[str] = [
    "pLDDT_mean",
    "pLDDT_min",
    "pae_mean",
    "contact_prob_mean",
    "ranking_score",
]

# Metrics for Figure 8 (pairplot)
FIGURE8_METRICS: List[str] = [
    "pLDDT_mean",
    "pae_mean",
    "contact_prob_mean",
    "ranking_score",
]


def metric_label(col: str) -> str:
    """Human-readable axis label for a metric column."""
    return METRIC_LABELS.get(col, col.replace("_", " ").title())


def panel_letter(col: str) -> str:
    """Single-letter panel identifier."""
    return METRIC_PANEL_LETTERS.get(col, "?")


def figure_title(fig_id: str) -> str:
    """Standardised figure title."""
    titles = {
        "fig7": "Distribution of AF3 Confidence Metrics Across Experimental Conditions",
        "fig8": "Relationships Among AF3 Confidence Metrics",
        "effects": "Factor Effects on AF3 Confidence Metrics",
    }
    return titles.get(fig_id, fig_id)
