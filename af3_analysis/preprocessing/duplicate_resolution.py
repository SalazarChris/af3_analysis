"""
Duplicate resolution for AF3 Confidence Analysis Pipeline.

Implements deterministic rules for resolving duplicate prediction candidates.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import pandas as pd


@dataclass(frozen=True)
class DuplicateGroup:
    """Group of records that are potential duplicates."""
    candidate_id: str
    members: List[Dict[str, Any]]
    resolution_rule: str
    resolution_reason: str


@dataclass(frozen=True)
class DuplicateResolution:
    """Resolution of duplicate records."""
    resolved_id: str
    kept_record: Dict[str, Any]
    excluded_records: List[Dict[str, Any]]
    resolution_rule: str
    resolution_reason: str


class DuplicateResolver:
    """
    Resolve duplicate prediction candidates deterministically.
    
    Rules:
    - Records with identical provenance (checksum, path) are duplicates
    - Records with identical (condition_id, seed, sample) are candidates
    - Resolve deterministically using a documented rule
    - Never silently merge - record all decisions
    """
    
    def __init__(self):
        self._duplicate_groups: List[DuplicateGroup] = []
        self._resolutions: List[DuplicateResolution] = []
    
    def find_duplicates(self, records: List[Dict[str, Any]]) -> List[DuplicateGroup]:
        """
        Find duplicate candidate groups.
        
        Args:
            records: List of record dicts with at least "condition_id", "seed", "sample"
        
        Returns:
            List of duplicate groups
        """
        # Group by (condition_id, seed, sample)
        groups: Dict[str, List[Dict[str, Any]]] = {}
        
        for record in records:
            key = f"{record.get('condition_id')}_seed-{record.get('seed')}_sample-{record.get('sample')}"
            if key not in groups:
                groups[key] = []
            groups[key].append(record)
        
        # Filter to groups with more than one member
        duplicates = []
        for key, members in groups.items():
            if len(members) > 1:
                duplicates.append(DuplicateGroup(
                    candidate_id=key,
                    members=members,
                    resolution_rule="",
                    resolution_reason="",
                ))
        
        self._duplicate_groups = duplicates
        return duplicates
    
    def resolve_duplicates(
        self,
        duplicate_groups: List[DuplicateGroup],
    ) -> List[DuplicateResolution]:
        """
        Resolve duplicate groups using deterministic rules.
        
        Rules (in order):
        1. If one record has provenance checksum, keep that one
        2. If both have checksums but different, keep the first
        3. If neither has checksum, keep the first
        
        Args:
            duplicate_groups: List of duplicate groups
        
        Returns:
            List of resolutions
        """
        resolutions = []
        
        for group in duplicate_groups:
            members = group.members
            
            # Try to find record with checksum
            with_checksum = [m for m in members if m.get("checksum")]
            
            if len(with_checksum) == 1:
                # Keep the one with checksum
                kept = with_checksum[0]
                excluded = [m for m in members if m != kept]
                resolution = DuplicateResolution(
                    resolved_id=group.candidate_id,
                    kept_record=kept,
                    excluded_records=excluded,
                    resolution_rule="keep_record_with_checksum",
                    resolution_reason="Record has file checksum for provenance",
                )
            else:
                # Keep first record
                kept = members[0]
                excluded = members[1:]
                resolution = DuplicateResolution(
                    resolved_id=group.candidate_id,
                    kept_record=kept,
                    excluded_records=excluded,
                    resolution_rule="keep_first_record",
                    resolution_reason="No distinguishing provenance; keeping first record",
                )
            
            resolutions.append(resolution)
        
        self._resolutions = resolutions
        return resolutions
    
    def apply_resolutions(
        self,
        records: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """
        Apply resolutions to filter out duplicates.
        
        Args:
            records: Original list of records
        
        Returns:
            Filtered list with duplicates removed
        """
        excluded_ids = set()
        
        for resolution in self._resolutions:
            excluded_ids.add(id(resolution))
        
        # Filter out excluded records
        kept = []
        for record in records:
            record_id = f"{record.get('condition_id')}_seed-{record.get('seed')}_sample-{record.get('sample')}"
            
            # Check if this record was excluded
            is_excluded = False
            for resolution in self._resolutions:
                if resolution.resolved_id == record_id:
                    for excluded in resolution.excluded_records:
                        if excluded is record:
                            is_excluded = True
                            break
            
            if not is_excluded:
                kept.append(record)
        
        return kept
    
    def get_resolution_log(self) -> List[Dict[str, Any]]:
        """Get a log of all resolutions for audit."""
        return [
            {
                "resolved_id": r.resolved_id,
                "kept_record": r.kept_record.get("prediction_id"),
                "excluded_count": len(r.excluded_records),
                "resolution_rule": r.resolution_rule,
                "resolution_reason": r.resolution_reason,
            }
            for r in self._resolutions
        ]


def resolve_duplicates(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Convenience function to resolve duplicates.
    
    Args:
        records: List of record dicts
    
    Returns:
        List of deduplicated records
    """
    resolver = DuplicateResolver()
    groups = resolver.find_duplicates(records)
    resolver.resolve_duplicates(groups)
    return resolver.apply_resolutions(records)


__all__ = [
    "DuplicateGroup",
    "DuplicateResolution",
    "DuplicateResolver",
    "resolve_duplicates",
]
