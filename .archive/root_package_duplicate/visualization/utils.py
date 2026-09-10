from typing import Any, Dict, List, Optional, Tuple

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def get_display_label(
    identifier: str,
    mapping_type: str = 'metric',
    schema: Optional[Dict[str, Any]] = None,
) -> str:
    """Returns a human-readable display label for a given identifier.

    Parameters
    ----------
    identifier : str
        The identifier to label (condition name, metric name, etc.).
    mapping_type : str
        Either ``'metric'`` or ``'condition'``.
    schema : dict, optional
        Schema dict that may contain an ``'experiment_design'`` key
        with an :class:`ExperimentDesign` instance.  Used for
        condition labels when available.
    """
    if not isinstance(identifier, str):
        return str(identifier)

    # --- Condition labels: prefer metadata if available ---
    if mapping_type == 'condition' and schema is not None:
        design = schema.get('experiment_design')
        if design is not None:
            label = design.get_condition_label(identifier)
            if label != identifier:  # metadata provided a real label
                return label

    labels = {
        'metric': {
            'plddt_mean': 'Mean pLDDT',
            'plddt_max': 'Max pLDDT',
            'plddt_min': 'Min pLDDT',
            'plddt_median': 'Median pLDDT',
            'pae_mean': 'Mean PAE',
            'pae_max': 'Max PAE',
            'pae_min': 'Min PAE',
            'pae_median': 'Median PAE',
            'contact_prob_mean': 'Mean Contact Probability',
            'contact_prob_max': 'Max Contact Probability',
            'contact_prob_min': 'Min Contact Probability',
            'contact_prob_median': 'Median Contact Probability',
            'ranking_score': 'Ranking Score',
            'fraction_disordered': 'Fraction Disordered',
            'ptm_score': 'PTM Score',
            'clash_score': 'Clash Score',
        },
    }
    
    ident_lower = identifier.lower()
    
    if mapping_type == 'metric':
        if ident_lower.startswith('chain_') and ident_lower.endswith('_plddt'):
            chain_letter = identifier.split('_')[1].upper()
            return f"Chain {chain_letter} pLDDT"
        if ident_lower.startswith('chain_') and ident_lower.endswith('_residues'):
            chain_letter = identifier.split('_')[1].upper()
            return f"Chain {chain_letter} Residue Count"
            
        mapping = labels.get('metric', {})
        if ident_lower in mapping:
            return mapping[ident_lower]
            
    elif mapping_type == 'condition':
        # Fallback for unknown conditions (no metadata)
        clean = identifier.replace('_', ' ').title()
        return clean

    # Generic fallback: replace underscores and title case
    return identifier.replace('_', ' ').title()


def get_condition_order(conditions: list, schema: dict = None) -> list:
    """Returns a deterministic ordering for conditions.

    Uses metadata-defined order when available; otherwise sorts alphabetically.
    """
    # Prefer schema-defined order (from metadata)
    if schema and 'condition_order' in schema and schema['condition_order']:
        order = schema['condition_order']
        return [c for c in order if c in conditions] + [c for c in conditions if c not in order]
    
    # Alphabetical fallback (no POU-specific hard-coding)
    return sorted(list(conditions))


def get_condition_style_map(conditions: list, schema: dict = None) -> dict:
    """Returns a consistent color palette mapping for conditions.

    Uses ``husl`` which provides perceptually distinct colors up to 256+
    categories.  ``Set2`` only has 8 colors and would silently repeat,
    making the plot unreadable for >8 conditions.
    """
    ordered_conditions = get_condition_order(conditions, schema)
    n = len(ordered_conditions)
    palette = sns.color_palette("husl", n_colors=max(n, 1))
    return dict(zip(ordered_conditions, palette))


# ---------------------------------------------------------------------------
# Preflight / safety-limit helpers
# ---------------------------------------------------------------------------

# Absolute upper bounds — exceed these and rendering either hangs or produces
# an unusable image.
MAX_FIGURE_HEIGHT_INCHES = 50.0
MAX_FIGURE_WIDTH_INCHES = 30.0
MAX_ANNO_COUNT = 500       # max text annotations in a single figure
MAX_TICK_LABELS = 150      # max tick labels on any single axis
MAX_HEATMAP_CELLS = 500    # max cells before annotations are suppressed
MAX_PAIRPLOT_GROUPS = 40   # max hue groups in a pairplot / ECDF
MAX_Y_LABELS_FOREST = 100  # max y-axis labels in the forest plot
MAX_ECDF_LINES = 30        # max simultaneous ECDF lines
DPI = 100


def preflight_check(
    fig_name: str,
    n_rows: int = 0,
    n_cols: int = 0,
    n_axes: int = 0,
    fig_width: float = 0.0,
    fig_height: float = 0.0,
    n_annotations: int = 0,
    n_tick_labels_x: int = 0,
    n_tick_labels_y: int = 0,
    n_heatmap_cells: int = 0,
    n_legend_entries: int = 0,
) -> Tuple[bool, List[str]]:
    """Check a figure against safety limits before rendering.

    Returns ``(ok, warnings)`` where *ok* is ``True`` if all hard limits
    are satisfied and *warnings* lists any soft-limit breaches.
    """
    warnings_list: List[str] = []

    width_px = fig_width * DPI
    height_px = fig_height * DPI

    if fig_height > MAX_FIGURE_HEIGHT_INCHES:
        return False, [
            f"{fig_name}: figure height {fig_height:.1f}in "
            f"exceeds max {MAX_FIGURE_HEIGHT_INCHES}in "
            f"({height_px:.0f}px)"
        ]
    if fig_width > MAX_FIGURE_WIDTH_INCHES:
        return False, [
            f"{fig_name}: figure width {fig_width:.1f}in "
            f"exceeds max {MAX_FIGURE_WIDTH_INCHES}in "
            f"({width_px:.0f}px)"
        ]

    if n_annotations > MAX_ANNO_COUNT:
        warnings_list.append(
            f"{fig_name}: {n_annotations} annotations exceeds soft limit "
            f"{MAX_ANNO_COUNT} — annotations may be suppressed"
        )
    if n_tick_labels_x > MAX_TICK_LABELS:
        warnings_list.append(
            f"{fig_name}: {n_tick_labels_x} x-tick labels exceeds "
            f"{MAX_TICK_LABELS}"
        )
    if n_tick_labels_y > MAX_TICK_LABELS:
        warnings_list.append(
            f"{fig_name}: {n_tick_labels_y} y-tick labels exceeds "
            f"{MAX_TICK_LABELS}"
        )
    if n_heatmap_cells > MAX_HEATMAP_CELLS:
        warnings_list.append(
            f"{fig_name}: {n_heatmap_cells} heatmap cells exceeds "
            f"{MAX_HEATMAP_CELLS} — annotations will be suppressed"
        )
    if n_legend_entries > MAX_PAIRPLOT_GROUPS:
        warnings_list.append(
            f"{fig_name}: {n_legend_entries} legend entries exceeds "
            f"{MAX_PAIRPLOT_GROUPS}"
        )

    return True, warnings_list


def format_figure_layout(fig, title: str):
    """Applies standardized formatting to a figure."""
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    

def adjust_legend(ax, **kwargs):
    """Adjusts legend to avoid covering data."""
    if ax.get_legend() is not None:
        sns.move_legend(ax, "center left", bbox_to_anchor=(1.02, 0.5), **kwargs)

