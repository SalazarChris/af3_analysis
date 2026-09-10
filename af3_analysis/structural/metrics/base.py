"""
Abstract base class for structural metrics.

Every metric declares what it requires (applicability) and computes
a value for a single pair of aligned structures.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Optional

from af3_analysis.structural.alignment import AlignmentResult
from af3_analysis.structural.enums import MetricDirection, MetricScope
from af3_analysis.structural.representation import NormalisedStructure


class StructuralMetric(ABC):
    """
    Abstract base for a pluggable structural metric.

    Subclass and implement all abstract methods. Register instances
    with the MetricRegistry.
    """

    @property
    @abstractmethod
    def metric_id(self) -> str:
        """Unique identifier, e.g. 'rmsd_global_ca'."""
        ...

    @property
    @abstractmethod
    def display_name(self) -> str:
        """Human-readable name."""
        ...

    @property
    @abstractmethod
    def units(self) -> str:
        """Measurement units, e.g. 'Å', 'Å²', 'count'."""
        ...

    @property
    @abstractmethod
    def scope(self) -> MetricScope:
        """Scope of the metric."""
        ...

    @property
    @abstractmethod
    def direction(self) -> MetricDirection:
        """Whether lower or higher values are geometrically meaningful."""
        ...

    @abstractmethod
    def requires(self) -> List[str]:
        """
        Declares entity types needed.

        Examples:
            []                           — always applicable
            ["protein"]                  — requires at least one protein chain
            ["protein", "protein"]       — requires at least two protein chains
            ["protein", "dna"]           — requires protein + nucleic acid

        The comparison layer checks applicability before invoking the metric.
        """
        ...

    @abstractmethod
    def compute(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
        alignment: AlignmentResult,
        **kwargs,
    ) -> Optional[float]:
        """
        Compute the metric for one pair of aligned structures.

        Returns None if the metric is undefined for this pair.
        """
        ...

    def is_applicable(
        self,
        structure_a: NormalisedStructure,
        structure_b: NormalisedStructure,
    ) -> bool:
        """Check whether this metric's requirements are met."""
        required = self.requires()
        if not required:
            return True

        types_a = structure_a.get_entity_types()
        types_b = structure_b.get_entity_types()
        all_types = types_a | types_b

        # Check that all required entity types are present in at least one structure
        # For multi-chain requirements (e.g., ["protein", "protein"]),
        # we need at least that many chains of the required type
        from collections import Counter
        type_counts = Counter()
        for t in required:
            type_counts[t] += 1

        # Count available chains of each type
        available = Counter()
        for chain in structure_a.chains:
            if chain.polymer_type and "polypeptide" in chain.polymer_type:
                available["protein"] += 1
            elif chain.polymer_type and "polydeoxyribonucleotide" in chain.polymer_type:
                available["dna"] += 1
            elif chain.polymer_type and "polyribonucleotide" in chain.polymer_type:
                available["rna"] += 1
            elif chain.entity_type == "non-polymer":
                available["ligand"] += 1
            elif chain.entity_type == "water":
                available["water"] += 1
            else:
                available[chain.entity_type] += 1
        for chain in structure_b.chains:
            if chain.polymer_type and "polypeptide" in chain.polymer_type:
                available["protein"] += 1
            elif chain.polymer_type and "polydeoxyribonucleotide" in chain.polymer_type:
                available["dna"] += 1
            elif chain.polymer_type and "polyribonucleotide" in chain.polymer_type:
                available["rna"] += 1
            elif chain.entity_type == "non-polymer":
                available["ligand"] += 1
            elif chain.entity_type == "water":
                available["water"] += 1
            else:
                available[chain.entity_type] += 1

        for entity_type, count_needed in type_counts.items():
            if available.get(entity_type, 0) < count_needed:
                return False
        return True
