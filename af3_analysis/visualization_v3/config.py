"""
V3 Visualization Pipeline Configuration.

All fields have project-agnostic defaults. No biological assumptions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Output defaults
# ---------------------------------------------------------------------------

DPI: int = 300
OUTPUT_FORMAT: str = "png"  # "png" or "pdf"

# Single-column (7-9 in) and double-column (10-13 in)
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
# Figure safety limits
# ---------------------------------------------------------------------------

MAX_FIGURE_HEIGHT_INCHES: float = 16.0
MAX_FIGURE_WIDTH_INCHES: float = 20.0
MAX_LEGEND_ENTRIES: int = 12
MAX_HEATMAP_CELLS: int = 500


# ---------------------------------------------------------------------------
# V3 configuration
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class V3FigureConfig:
    """Configuration for a single V3 figure."""
    figure_id: str
    enabled: bool = True
    title: Optional[str] = None
    # Figure-specific overrides can be added here


@dataclass(frozen=True)
class V3StructuralConfig:
    """Configuration for structural calculations."""
    # Alignment
    alignment_atom: str = "CA"  # CA, backbone, all_heavy
    min_common_atoms: int = 3
    min_sequence_identity: float = 0.5

    # Contacts
    contact_distance: float = 8.0  # Angstrom

    # Coverage
    minimum_coverage: float = 0.80

    # Output
    output_dpi: int = 300
    output_format: str = "png"

    # Seed handling
    seed_comparison_mode: str = "matched_seed"  # matched_seed, all_vs_all


@dataclass(frozen=True)
class V3ClusteringConfig:
    """Configuration for structural clustering."""
    method: str = "hierarchical"  # hierarchical, k_medoids
    n_clusters: Optional[int] = None  # None = auto/determine from data
    linkage: str = "average"  # single, complete, average, ward
    distance_metric: str = "rmsd_ca"  # rmsd_ca, rmsd_backbone, contact_overlap


@dataclass(frozen=True)
class V3Config:
    """Main V3 visualization configuration."""

    # Enable/disable
    enabled: bool = True

    # Reference resolution
    reference: Optional[Dict[str, Any]] = field(default_factory=dict)
    # Example: {"condition": "<condition_id_or_name>"} or
    # {"condition_a": "X", "condition_b": "Y"}

    # Figure toggles
    figures: Dict[str, bool] = field(default_factory=dict)
    # Default: enable all figures
    # Can be overridden: {"F01": True, "F02": False, ...}

    # Structural configuration
    structure: V3StructuralConfig = field(default_factory=V3StructuralConfig)

    # Clustering configuration
    clustering: V3ClusteringConfig = field(default_factory=V3ClusteringConfig)

    # Region/site definitions (for local geometry)
    sites: List[Dict[str, Any]] = field(default_factory=list)
    # Example: [{"label": "site1", "residue": 123, "radius": 5}]

    # Region/domain definitions (for domain motion)
    regions: List[Dict[str, Any]] = field(default_factory=list)
    # Example: [{"label": "domain1", "chain": "A", "start": 1, "end": 50}]

    # Data-driven detection of sites/regions when no manual definitions
    # are supplied. Manual definitions always take priority.
    auto_sites: bool = False
    auto_regions: bool = False

    # Output settings
    output_dpi: int = 300
    output_format: str = "png"

    # Logging
    verbose: bool = True

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "V3Config":
        """Create from dictionary (typically JSON)."""
        # Parse structural config
        struct_data = data.get("structure", {})
        structure = V3StructuralConfig(
            alignment_atom=struct_data.get("alignment_atom", "CA"),
            min_common_atoms=struct_data.get("min_common_atoms", 3),
            min_sequence_identity=struct_data.get("min_sequence_identity", 0.5),
            contact_distance=struct_data.get("contact_distance", 8.0),
            minimum_coverage=struct_data.get("minimum_coverage", 0.80),
            output_dpi=struct_data.get("output_dpi", 300),
            output_format=struct_data.get("output_format", "png"),
            seed_comparison_mode=struct_data.get("seed_comparison_mode", "matched_seed"),
        )

        # Parse clustering config
        cluster_data = data.get("clustering", {})
        clustering = V3ClusteringConfig(
            method=cluster_data.get("method", "hierarchical"),
            n_clusters=cluster_data.get("n_clusters"),
            linkage=cluster_data.get("linkage", "average"),
            distance_metric=cluster_data.get("distance_metric", "rmsd_ca"),
        )

        # Parse reference
        reference = data.get("reference", {})

        # Parse figures
        figures_data = data.get("figures", {})
        figures = {k: v for k, v in figures_data.items()}

        # Parse sites
        sites = data.get("sites", [])

        # Parse regions
        regions = data.get("regions", [])

        return cls(
            enabled=data.get("enabled", True),
            reference=reference,
            figures=figures,
            structure=structure,
            clustering=clustering,
            sites=sites,
            regions=regions,
            output_dpi=data.get("output_dpi", 300),
            output_format=data.get("output_format", "png"),
            verbose=data.get("verbose", True),
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "enabled": self.enabled,
            "reference": self.reference,
            "figures": self.figures,
            "structure": {
                "alignment_atom": self.structure.alignment_atom,
                "min_common_atoms": self.structure.min_common_atoms,
                "min_sequence_identity": self.structure.min_sequence_identity,
                "contact_distance": self.structure.contact_distance,
                "minimum_coverage": self.structure.minimum_coverage,
                "output_dpi": self.structure.output_dpi,
                "output_format": self.structure.output_format,
                "seed_comparison_mode": self.structure.seed_comparison_mode,
            },
            "clustering": {
                "method": self.clustering.method,
                "n_clusters": self.clustering.n_clusters,
                "linkage": self.clustering.linkage,
                "distance_metric": self.clustering.distance_metric,
            },
            "sites": self.sites,
            "regions": self.regions,
            "output_dpi": self.output_dpi,
            "output_format": self.output_format,
            "verbose": self.verbose,
        }


def create_v3_config_interactive(
    enabled: bool = True,
    reference_condition: Optional[str] = None,
    output_dpi: int = 300,
    output_format: str = "png",
) -> V3Config:
    """Create a V3 configuration with sensible defaults."""
    # Leave the figures dict empty so get_enabled_figures applies the
    # default suite (all figures except V3_DEFAULT_OFF_FIGURES).
    figures: Dict[str, bool] = {}

    reference = {}
    if reference_condition:
        reference["condition"] = reference_condition

    return V3Config(
        enabled=enabled,
        reference=reference,
        figures=figures,
        output_dpi=output_dpi,
        output_format=output_format,
    )


# Default figure list (all 20 figures)
V3_ALL_FIGURES = [f"F{i:02d}" for i in range(1, 21)]

# Figures excluded from the default suite. They remain fully available via
# explicit opt-in (config figures dict or the CLI --figures list). Rationale
# (docs/V3_STRUCTURAL_VISUALIZATION_REVIEW.md):
#   F01 — QC is reported by validation/tables; the figure adds no information
#         when all structures parse.
#   F06 — per-pair contact detail is consolidated into F07 (which renders
#         the recurring-pair panel itself when F06 is not part of the run).
#   F08 / F09 — interface figures are data-dependent and often reflect
#         sparse baselines rather than comparable changes.
#   F11 — domain motion is only meaningful with meaningful region
#         definitions; it skips without them.
#   F12 — the cluster-assignment presentation is hard to read/interpret
#         (dense distance matrix + sparse categorical strip); the cluster
#         information remains available via F13B's ordering and the
#         structural_clusters.csv table. Opt in via figures config.
V3_DEFAULT_OFF_FIGURES = ("F01", "F06", "F08", "F09", "F11", "F12")

# F13 renders as separate views (condition similarity / within-condition
# reproducibility / prediction-level matrix). The prediction-level matrix
# views are hidden from the default output: at hundreds of predictions the
# matrices are not human-interpretable. F13B (within-condition
# reproducibility) remains the default F13 view. Hidden views can be
# re-enabled with figures entries {"F13A": True, "F13D": True}.
V3_F13_HIDDEN_VIEWS = ("F13A", "F13D")


def get_enabled_figures(config: V3Config) -> List[str]:
    """Return list of enabled figure IDs from configuration.

    With no explicit figure toggles, the default suite is all figures
    except ``V3_DEFAULT_OFF_FIGURES``. With a partial figures dict, entries
    are overrides: an explicit True enables even a default-off figure, and
    an explicit False disables even a default-on one.
    """
    if not config.figures:
        return [f for f in V3_ALL_FIGURES if f not in V3_DEFAULT_OFF_FIGURES]

    enabled = []
    for fig_id in V3_ALL_FIGURES:
        default_on = fig_id not in V3_DEFAULT_OFF_FIGURES
        if config.figures.get(fig_id, default_on):
            enabled.append(fig_id)
    return enabled


__all__ = [
    "V3Config",
    "V3StructuralConfig",
    "V3ClusteringConfig",
    "V3FigureConfig",
    "load_v3_config",
    "create_v3_config_interactive",
    "get_enabled_figures",
    "V3_ALL_FIGURES",
    "V3_DEFAULT_OFF_FIGURES",
    "V3_F13_HIDDEN_VIEWS",
    "DPI",
    "OUTPUT_FORMAT",
    "SINGLE_COL_WIDTH",
    "DOUBLE_COL_WIDTH",
    "FONT_SIZES",
    "apply_v3_style",
]


def apply_v3_style() -> None:
    """Apply a consistent matplotlib/seaborn style for V3 figures."""
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
