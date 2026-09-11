#!/usr/bin/env python3
"""
matcher.py - Condition matching for AF3 Server predictions.
"""

from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import re


class ConditionMatcher:
    """
    Match AF3 prediction folders to conditions in the registry.
    
    Supports multiple matching strategies:
    - Strategy A: Explicit mapping (prediction_folder field in metadata)
    - Strategy B: Folder name matching (exact match)
    - Strategy C: Pattern matching (extract identifiers)
    - Strategy D: Fallback to first condition
    """
    
    def __init__(self, registry):
        self.registry = registry
        self.matches: Dict[str, Dict[str, Any]] = {}  # folder_name -> condition
        self.unmatched: List[str] = []
        
    def match_folder_to_condition(self, folder_name: str) -> Optional[Dict[str, Any]]:
        """
        Match a folder name to a condition using multiple strategies.
        
        Args:
            folder_name: Name of AF3 prediction folder
            
        Returns:
            Matched condition or None
        """
        # Strategy A: Check for explicit prediction_folder mapping
        condition = self._match_explicit(folder_name)
        if condition:
            return condition
            
        # Strategy B: Direct name match (exact match)
        condition = self._match_exact(folder_name)
        if condition:
            return condition
            
        # Strategy C: Pattern matching
        condition = self._match_pattern(folder_name)
        if condition:
            return condition
            
        # Strategy D: Return first condition as fallback
        conditions = self.registry.get_all_conditions()
        if conditions:
            return conditions[0]
            
        return None
    
    def _match_explicit(self, folder_name: str) -> Optional[Dict[str, Any]]:
        """Strategy A: Match via explicit prediction_folder field."""
        for cond in self.registry.get_all_conditions():
            pred_folder = cond.get('prediction_folder')
            if pred_folder == folder_name:
                return cond
        return None
    
    def _match_exact(self, folder_name: str) -> Optional[Dict[str, Any]]:
        """Strategy B: Direct name match."""
        for cond in self.registry.get_all_conditions():
            if cond['condition_name'] == folder_name:
                return cond
        return None
    
    def _match_pattern(self, folder_name: str) -> Optional[Dict[str, Any]]:
        """Strategy C: Pattern-based matching."""
        conditions = self.registry.get_all_conditions()
        
        # Try to extract base name (e.g., '<base>' from '<base>_na100_hoh1000_cl100')
        base_patterns = self._extract_base_patterns(folder_name)
        
        for cond in conditions:
            cond_name = cond['condition_name']
            
            # Check if condition name is a prefix or base pattern
            if cond_name == folder_name:
                return cond
                
            for base_pattern in base_patterns:
                if cond_name == base_pattern:
                    return cond
                    
            # Try substring match (for cases like '<cond_a>' matching '<cond_a>_na100')
            if folder_name.startswith(cond_name + '_') or folder_name.startswith(cond_name):
                return cond
                
        return None
    
    def _extract_base_patterns(self, folder_name: str) -> List[str]:
        """Extract potential base patterns from folder name."""
        patterns = []
        parts = folder_name.split('_')
        
        # Try progressively shorter prefix combinations
        for i in range(1, min(len(parts) + 1, 4)):
            pattern = '_'.join(parts[:i])
            if pattern and pattern != folder_name:
                patterns.append(pattern)
                
        return patterns
    
    def match_all_folders(self, input_dir: str) -> Tuple[Dict[str, Dict[str, Any]], List[str]]:
        """
        Match all AF3 prediction folders to conditions.
        
        Args:
            input_dir: Root directory containing prediction folders
            
        Returns:
            Tuple of (matched_conditions, unmatched_folders)
        """
        input_path = Path(input_dir)
        matched = {}
        unmatched = []
        
        # Find all condition folders (excluding outputs, cache, etc.)
        condition_folders = []
        for item in input_path.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                # Skip outputs and other special directories
                if item.name in ['outputs', 'cache', '__pycache__', 'templates']:
                    continue
                # Should contain seed-XX directories
                if any(subdir.name.startswith('seed-') for subdir in item.iterdir()):
                    condition_folders.append(item)
        
        # Match each folder
        for folder in condition_folders:
            condition = self.match_folder_to_condition(folder.name)
            if condition:
                matched[folder.name] = condition
                print(f"  Matched: {folder.name} -> {condition['condition_name']} ({condition['internal_id']})")
            else:
                unmatched.append(folder.name)
                print(f"  UNMATCHED: {folder.name}")
                
        self.matches = matched
        self.unmatched = unmatched
        
        return matched, unmatched
    
    def get_match_summary(self) -> Dict[str, Any]:
        """Generate summary of matches."""
        return {
            'total_folders': len(self.matches) + len(self.unmatched),
            'matched': len(self.matches),
            'unmatched': len(self.unmatched),
            'unmatched_folders': self.unmatched,
            'match_details': self.matches
        }
