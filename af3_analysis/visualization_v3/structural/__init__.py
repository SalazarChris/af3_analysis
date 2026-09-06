"""
V3 Structural Analysis Submodules.

Contains structural comparison, RMSD, displacement, contacts, interfaces,
clustering, regions, and pairwise analysis.
"""

from .comparability import (
    find_common_structural_space,
    CommonStructuralSpace,
    has_sufficient_coverage,
    get_common_ca_coords,
    get_common_backbone_coords,
)
from .rmsd import calculate_rmsd, calculate_rmsd_matrix
from .displacement import (
    calculate_per_residue_displacement,
    calculate_displacement_vector,
)
from .contacts import (
    calculate_contact_map,
    calculate_contact_difference,
    ContactMap,
)
from .interfaces import (
    find_interface_contacts,
    InterfaceContacts,
)
from .pairwise import calculate_pairwise_rmsd_matrix
from .clustering import (
    hierarchical_clustering,
    kmeans_clustering,
   mds_embedding,
)
from .regions import (
    calculate_local_geometry,
    RegionResult,
)

__all__ = [
    "find_common_structural_space",
    "CommonStructuralSpace",
    "has_sufficient_coverage",
    "get_common_ca_coords",
    "get_common_backbone_coords",
    "calculate_rmsd",
    "calculate_rmsd_matrix",
    "calculate_per_residue_displacement",
    "calculate_displacement_vector",
    "calculate_contact_map",
    "calculate_contact_difference",
    "ContactMap",
    "find_interface_contacts",
    "InterfaceContacts",
    "calculate_pairwise_rmsd_matrix",
    "hierarchical_clustering",
    "kmeans_clustering",
    "mds_embedding",
    "calculate_local_geometry",
    "RegionResult",
]
