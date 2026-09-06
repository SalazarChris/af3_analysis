"""
V3 Figure Generators.

Contains all 20 figure generators (V3-F01 through V3-F20).
Each figure has:
- unique ID
- clear title
- axis labels
- documented unit of analysis
- documented metric
- valid observation count
- reference information
- warnings where applicable
"""

from .f01_structural_qc import generate_f01_structural_qc
from .f02_global_structural_difference import generate_f02_global_structural_difference
from .f03_matched_seed_structural_difference import generate_f03_matched_seed_structural_difference
from .f04_per_residue_displacement import generate_f04_per_residue_displacement
from .f05_displacement_heatmap import generate_f05_displacement_heatmap
from .f06_contact_map_difference import generate_f06_contact_map_difference
from .f07_contact_change_summary import generate_f07_contact_change_summary
from .f08_interface_analysis import generate_f08_interface_analysis
from .f09_interface_change_map import generate_f09_interface_change_map
from .f10_local_geometry import generate_f10_local_geometry
from .f11_domain_motion import generate_f11_domain_motion
from .f12_structural_clustering import generate_f12_structural_clustering
from .f13_similarity_matrix import generate_f13_similarity_matrix
from .f14_mds_embedding import generate_f14_mds_embedding
from .f15_confidence_geometry import generate_f15_confidence_geometry
from .f16_confidence_change_vs_structural_change import generate_f16_confidence_change_vs_structural_change
from .f17_seed_reproducibility import generate_f17_seed_reproducibility
from .f18_effect_sizes import generate_f18_effect_sizes
from .f19_factorial_effects import generate_f19_factorial_effects
from .f20_structure_confidence_matrix import generate_f20_structure_confidence_matrix

__all__ = [
    "generate_f01_structural_qc",
    "generate_f02_global_structural_difference",
    "generate_f03_matched_seed_structural_difference",
    "generate_f04_per_residue_displacement",
    "generate_f05_displacement_heatmap",
    "generate_f06_contact_map_difference",
    "generate_f07_contact_change_summary",
    "generate_f08_interface_analysis",
    "generate_f09_interface_change_map",
    "generate_f10_local_geometry",
    "generate_f11_domain_motion",
    "generate_f12_structural_clustering",
    "generate_f13_similarity_matrix",
    "generate_f14_mds_embedding",
    "generate_f15_confidence_geometry",
    "generate_f16_confidence_change_vs_structural_change",
    "generate_f17_seed_reproducibility",
    "generate_f18_effect_sizes",
    "generate_f19_factorial_effects",
    "generate_f20_structure_confidence_matrix",
]
