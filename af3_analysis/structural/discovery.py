"""
Structural discovery adapter.

Scans a raw AF3 output root for CIF structure files, maps them to
condition/seed/sample identifiers, validates against experiment metadata,
and produces a comprehensive structural dataset inventory.

Design principles:
- No biological assumptions.
- Condition identity comes from directory names or experiment metadata.
- Seed and sample are extracted from the directory/file naming convention.
- Validation is metadata-driven, not hardcoded.
"""

from __future__ import annotations

import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Data structures
# ------------------------------------------------------------------

@dataclass
class StructureRecord:
    """One discovered CIF file with its metadata."""
    path: Path
    condition_id: str
    seed: int
    sample: int
    checksum: str = ""

    @property
    def prediction_id(self) -> str:
        return f"{self.condition_id}_seed-{self.seed}_sample-{self.sample}"


@dataclass
class ConditionInventory:
    """Inventory of structural predictions for one condition."""
    condition_id: str
    n_predictions: int
    unique_seeds: List[int]
    unique_samples: List[int]
    records: List[StructureRecord] = field(default_factory=list)

    # Entity composition (filled after parsing)
    entity_compositions: List[Dict[str, Any]] = field(default_factory=list)
    chain_sets: List[Set[str]] = field(default_factory=list)
    polymer_type_sets: List[Set[str]] = field(default_factory=list)

    # Validation
    missing_seeds: List[int] = field(default_factory=list)
    missing_samples: List[Tuple[int, int]] = field(default_factory=list)
    duplicate_predictions: List[str] = field(default_factory=list)


@dataclass
class DiscoveryReport:
    """Complete discovery and validation report for the structural dataset."""
    root_path: Path
    n_total_cif_files: int
    n_discovered: int
    n_parse_success: int
    n_parse_failure: int
    n_orphan: int  # files not matching any known condition

    conditions: Dict[str, ConditionInventory] = field(default_factory=dict)
    parse_errors: List[Dict[str, Any]] = field(default_factory=list)
    orphan_files: List[Path] = field(default_factory=list)
    validation_warnings: List[str] = field(default_factory=list)
    validation_errors: List[str] = field(default_factory=list)

    @property
    def n_conditions(self) -> int:
        return len(self.conditions)

    @property
    def expected_seeds(self) -> Optional[List[int]]:
        """If all conditions share the same seeds, return them."""
        if not self.conditions:
            return None
        seed_sets = [set(inv.unique_seeds) for inv in self.conditions.values()]
        if all(s == seed_sets[0] for s in seed_sets):
            return sorted(seed_sets[0])
        return None

    @property
    def expected_samples_per_seed(self) -> Optional[int]:
        """If all conditions have the same samples per seed, return the count."""
        if not self.conditions:
            return None
        sample_counts = []
        for inv in self.conditions.values():
            if inv.unique_seeds:
                n = len(inv.unique_samples)
                sample_counts.append(n)
        if len(set(sample_counts)) == 1:
            return sample_counts[0]
        return None

    def summary(self) -> str:
        """Human-readable summary."""
        lines = [
            f"Structural Discovery Report",
            f"  Root: {self.root_path}",
            f"  CIF files found: {self.n_total_cif_files}",
            f"  Successfully discovered: {self.n_discovered}",
            f"  Parse success: {self.n_parse_success}",
            f"  Parse failures: {self.n_parse_failure}",
            f"  Orphan files: {self.n_orphan}",
            f"  Conditions: {self.n_conditions}",
            "",
        ]
        for cond_id in sorted(self.conditions):
            inv = self.conditions[cond_id]
            lines.append(f"  [{cond_id}]")
            lines.append(f"    Predictions: {inv.n_predictions}")
            lines.append(f"    Seeds: {sorted(inv.unique_seeds)}")
            lines.append(f"    Samples: {sorted(inv.unique_samples)}")
            if inv.missing_seeds:
                lines.append(f"    Missing seeds: {inv.missing_seeds}")
            if inv.duplicate_predictions:
                lines.append(f"    Duplicates: {len(inv.duplicate_predictions)}")
            if inv.chain_sets:
                unique_chains = set(frozenset(cs) for cs in inv.chain_sets)
                lines.append(f"    Chain sets: {[sorted(cs) for cs in unique_chains]}")
            if inv.polymer_type_sets:
                unique_pts = set(frozenset(pts) for pts in inv.polymer_type_sets)
                lines.append(f"    Polymer types: {[sorted(pts) for pts in unique_pts]}")
            lines.append("")

        if self.validation_warnings:
            lines.append("Warnings:")
            for w in self.validation_warnings:
                lines.append(f"  - {w}")
            lines.append("")

        if self.validation_errors:
            lines.append("Errors:")
            for e in self.validation_errors:
                lines.append(f"  - {e}")
            lines.append("")

        if self.parse_errors:
            lines.append(f"Parse errors ({len(self.parse_errors)}):")
            for err in self.parse_errors[:10]:
                lines.append(f"  - {err['condition_id']} s{err['seed']} m{err['sample']}: {err['reason']}")
            if len(self.parse_errors) > 10:
                lines.append(f"  ... and {len(self.parse_errors) - 10} more")
            lines.append("")

        return "\n".join(lines)


# ------------------------------------------------------------------
# Discovery
# ------------------------------------------------------------------

# Pattern for AF3 CIF files in seed-sample subdirectories
_CIF_NAME_PATTERN = re.compile(
    r"^(?P<cond>.+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model\.cif$"
)


def discover_structures(
    root: Path,
    *,
    expected_conditions: Optional[List[str]] = None,
) -> DiscoveryReport:
    """
    Discover AF3 structure files under a root directory.

    Parameters
    ----------
    root : Path
        Root directory containing condition subdirectories.
    expected_conditions : list of str, optional
        Expected condition IDs from experiment metadata.
        If provided, validates that discovered conditions match.

    Returns
    -------
    DiscoveryReport
        Complete discovery and validation report.
    """
    logger.info("[Structural] Discovering structures in %s", root)

    # Phase 1: Find all CIF files
    all_cif = list(root.rglob("*_model.cif"))
    all_cif = [p for p in all_cif if p.is_file()]
    logger.info("[Structural] Found %d CIF files", len(all_cif))

    # Phase 2: Parse condition/seed/sample from filenames
    records: List[StructureRecord] = []
    orphans: List[Path] = []
    seen: Set[Tuple[str, int, int]] = set()

    for cif_path in sorted(all_cif):
        match = _CIF_NAME_PATTERN.match(cif_path.name)
        if not match:
            orphans.append(cif_path)
            continue

        cond = match.group("cond")
        seed = int(match.group("seed"))
        sample = int(match.group("sample"))

        key = (cond, seed, sample)
        if key in seen:
            # Duplicate — keep first encountered
            continue
        seen.add(key)

        records.append(StructureRecord(
            path=cif_path,
            condition_id=cond,
            seed=seed,
            sample=sample,
        ))

    # Phase 3: Group by condition
    by_condition: Dict[str, List[StructureRecord]] = defaultdict(list)
    for rec in records:
        by_condition[rec.condition_id].append(rec)

    # Phase 4: Validate against expected conditions
    discovered_conditions = set(by_condition.keys())
    report = DiscoveryReport(
        root_path=root,
        n_total_cif_files=len(all_cif),
        n_discovered=len(records),
        n_parse_success=0,  # filled after parsing
        n_parse_failure=0,
        n_orphan=len(orphans),
        orphan_files=orphans,
    )

    if expected_conditions is not None:
        expected_set = set(expected_conditions)
        missing = expected_set - discovered_conditions
        extra = discovered_conditions - expected_set

        if missing:
            report.validation_errors.append(
                f"Conditions in metadata but not found on disk: {sorted(missing)}"
            )
        if extra:
            report.validation_warnings.append(
                f"Conditions found on disk but not in metadata: {sorted(extra)}"
            )

    # Phase 5: Build per-condition inventory
    # Determine expected seeds/samples from the first condition
    all_seeds: Set[int] = set()
    all_samples: Set[int] = set()
    for rec in records:
        all_seeds.add(rec.seed)
        all_samples.add(rec.sample)

    for cond_id in sorted(by_condition):
        cond_records = by_condition[cond_id]
        seeds = sorted(set(r.seed for r in cond_records))
        samples = sorted(set(r.sample for r in cond_records))

        # Check for completeness
        missing_seeds = sorted(set(all_seeds) - set(seeds))
        missing_samples = []
        for s in seeds:
            for m in all_samples:
                if not any(r.seed == s and r.sample == m for r in cond_records):
                    missing_samples.append((s, m))

        # Check for duplicates
        pred_keys = [(r.seed, r.sample) for r in cond_records]
        seen_keys = set()
        dupes = []
        for k in pred_keys:
            if k in seen_keys:
                dupes.append(f"s{k[0]}_m{k[1]}")
            seen_keys.add(k)

        inv = ConditionInventory(
            condition_id=cond_id,
            n_predictions=len(cond_records),
            unique_seeds=seeds,
            unique_samples=samples,
            records=cond_records,
            missing_seeds=missing_seeds,
            missing_samples=missing_samples,
            duplicate_predictions=dupes,
        )
        report.conditions[cond_id] = inv

        if missing_seeds:
            report.validation_warnings.append(
                f"Condition '{cond_id}': missing seeds {missing_seeds}"
            )
        if missing_samples:
            report.validation_warnings.append(
                f"Condition '{cond_id}': {len(missing_samples)} missing seed×sample combinations"
            )
        if dupes:
            report.validation_errors.append(
                f"Condition '{cond_id}': duplicate predictions: {dupes}"
            )

    logger.info(
        "[Structural] Discovered %d records across %d conditions",
        len(records), len(by_condition),
    )

    return report


def populate_composition(
    report: DiscoveryReport,
    parsed_structures: List[Any],  # NormalisedStructure objects
) -> None:
    """
    Populate entity composition and chain information from parsed structures.

    Modifies report.conditions in-place.

    Parameters
    ----------
    report : DiscoveryReport
        The discovery report to update.
    parsed_structures : list
        List of NormalisedStructure objects (already parsed).
    """
    # Group parsed structures by condition
    by_condition: Dict[str, List[Any]] = defaultdict(list)
    for s in parsed_structures:
        by_condition[s.condition_id].append(s)

    for cond_id, structs in by_condition.items():
        if cond_id not in report.conditions:
            continue
        inv = report.conditions[cond_id]
        for s in structs:
            inv.entity_compositions.append({
                "n_atoms": s.n_atoms,
                "n_chains": s.n_chains,
                "entity_types": sorted(s.get_entity_types()),
                "polymer_types": sorted(s.get_polymer_types()),
            })
            inv.chain_sets.append(set(s.chain_ids))
            inv.polymer_type_sets.append(s.get_polymer_types())

    # Report cross-condition composition differences
    all_compositions = {}
    for cond_id, inv in report.conditions.items():
        if inv.chain_sets:
            # Use the most common chain set as representative
            chain_set_counter = Counter(frozenset(cs) for cs in inv.chain_sets)
            most_common = chain_set_counter.most_common(1)[0][0]
            all_compositions[cond_id] = most_common

    # Check if all conditions have the same composition
    unique_compositions = set(frozenset(cs) for cs in all_compositions.values())
    if len(unique_compositions) > 1:
        report.validation_warnings.append(
            f"Entity composition varies across conditions: "
            f"{len(unique_compositions)} unique chain sets found"
        )


def get_all_records(report: DiscoveryReport) -> List[StructureRecord]:
    """Get all structure records from the report, sorted by condition/seed/sample."""
    records = []
    for inv in report.conditions.values():
        records.extend(inv.records)
    return sorted(records, key=lambda r: (r.condition_id, r.seed, r.sample))


__all__ = [
    "StructureRecord",
    "ConditionInventory",
    "DiscoveryReport",
    "discover_structures",
    "populate_composition",
    "get_all_records",
]
