"""
V2 Visualization Configuration
==============================
Centralized style settings for publication-quality AF3 figures.
All visual constants live here — individual plotting functions import from here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Output defaults
# ---------------------------------------------------------------------------

DPI: int = 300
OUTPUT_FORMAT: str = "png"  # "png" or "pdf"

# Single-column (7–9 in) and double-column (10–13 in)
SINGLE_COL_WIDTH: float = 7.5
DOUBLE_COL_WIDTH: float = 11.0

# ---------------------------------------------------------------------------
# Typography hierarchy
# ---------------------------------------------------------------------------

FONT_SIZES: Dict[str, int] = {
    "figure_title": 16,
    "panel_title": 12,
    "axis_label": 11,
    "tick_label": 9,
    "legend_title": 10,
    "legend_text": 9,
    "annotation": 8,
}

# ---------------------------------------------------------------------------
# Factor colour palette
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Factor colour palette (no hardcoded biological names)
# ---------------------------------------------------------------------------

# Fixed colour for Baseline; all other states get husl colours dynamically
_BASE_PTM_COLOUR: str = "#4C72B0"  # steel-blue for Baseline

# DNA → line-style / marker
DNA_STYLE: Dict[bool, Dict[str, str]] = {
    False: {"linestyle": "-",  "marker": "o", "label": "No DNA"},
    True:  {"linestyle": "--", "marker": "^", "label": "DNA"},
}

# Environment palette — one muted colour per env, used only in facets
ENV_PALETTE: str = "Set2"  # enough for ≤8 environments; fall back to "husl"

# Seed point style (used when individual seeds are overlaid)
SEED_STYLE: Dict[str, object] = {
    "alpha": 0.35,
    "s": 18,
    "linewidth": 0.3,
    "edgecolor": "white",
}

# ---------------------------------------------------------------------------
# Figure safety limits
# ---------------------------------------------------------------------------

MAX_FIGURE_HEIGHT_INCHES: float = 16.0
MAX_FIGURE_WIDTH_INCHES: float = 20.0
MAX_LEGEND_ENTRIES: int = 12
MAX_HEATMAP_CELLS: int = 500

# ---------------------------------------------------------------------------
# Derived helpers
# ---------------------------------------------------------------------------

DNA_ORDER: List[bool] = [False, True]


def derive_ptm_order(ptm_states: List[str]) -> List[str]:
    """Return stable ordering: Baseline first, then alphabetical."""
    others = sorted(s for s in ptm_states if s != "Baseline")
    return (["Baseline"] + others) if "Baseline" in ptm_states else others


def build_ptm_colours(ptm_states: List[str]) -> Dict[str, str]:
    """Generate a colour map for arbitrary PTM states.

    ``Baseline`` always gets the fixed blue; new states use husl.
    """
    palette: Dict[str, str] = {}
    husl = sns.color_palette("husl", n_colors=max(len(ptm_states), 8))
    husl_idx = 0
    for s in ptm_states:
        if s == "Baseline":
            palette[s] = _BASE_PTM_COLOUR
        else:
            palette[s] = mcolors.to_hex(husl[husl_idx])
            husl_idx += 1
    return palette


def apply_v2_style() -> None:
    """Apply a consistent matplotlib/seaborn style for V2 figures."""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
        "font.size": FONT_SIZES["tick_label"],
        "axes.titlesize": FONT_SIZES["panel_title"],
        "axes.labelsize": FONT_SIZES["axis_label"],
        "xtick.labelsize": FONT_SIZES["tick_label"],
        "ytick.labelsize": FONT_SIZES["tick_label"],
        "legend.fontsize": FONT_SIZES["legend_text"],
        "legend.title_fontsize": FONT_SIZES["legend_title"],
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.15,
    })


def get_ptm_color(ptm_state: str, palette: Optional[Dict[str, str]] = None) -> str:
    """Return colour for *ptm_state* using *palette* or fallback husl."""
    if palette and ptm_state in palette:
        return palette[ptm_state]
    if ptm_state == "Baseline":
        return _BASE_PTM_COLOUR
    idx = abs(hash(ptm_state)) % 8
    return mcolors.to_hex(sns.color_palette("husl", n_colors=8)[idx])


def get_dna_style(has_dna: bool) -> Dict[str, str]:
    """Return line-style dict for a DNA state."""
    return DNA_STYLE.get(has_dna, DNA_STYLE[False])


def make_ptm_palette(ptm_states: Optional[List[str]] = None) -> Dict[str, str]:
    """Return ordered colour mapping for the given *ptm_states*.

    If *ptm_states* is ``None``, falls back to ``["Baseline"]``.
    """
    states = ptm_states if ptm_states else ["Baseline"]
    return build_ptm_colours(derive_ptm_order(states))
