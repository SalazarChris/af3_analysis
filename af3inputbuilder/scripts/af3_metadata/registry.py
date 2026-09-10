#!/usr/bin/env python3
"""
registry.py - Condition registry for AF3 Server predictions.
"""

import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime


class ConditionRegistry:
    """
    Central registry for experimental conditions.
    
    Each condition receives a unique internal ID while preserving
    the original metadata. Supports multiple metadata sources.
    """
    
    def __init__(self):
        self.conditions: List[Dict[str, Any]] = []
        self.condition_ids: Dict[str, Dict[str, Any]] = {}  # condition_id -> condition
        self.source_map: Dict[str, List[str]] = {}  # source_file -> list of condition_ids
        self.next_internal_id = 1
        
    def add_condition(self, condition_data: Dict[str, Any], source_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Add a condition to the registry.
        
        Args:
            condition_data: Condition definition with name, factors, etc.
            source_file: Optional path to source metadata file
            
        Returns:
            Condition with internal ID added
        """
        # Generate unique internal ID
        internal_id = f"cond_{self.next_internal_id:04d}"
        self.next_internal_id += 1
        
        # Create condition record
        condition = {
            'internal_id': internal_id,
            'original_data': condition_data.copy(),
            'source_file': source_file,
            'created_at': datetime.now().isoformat()
        }
        
        # Extract name from original data
        condition_name = condition_data.get('condition_name') or condition_data.get('name') or internal_id
        condition['condition_name'] = condition_name
        
        # Extract factors
        factors = condition_data.get('factors', {})
        if not factors and 'factors' in condition_data:
            factors = condition_data['factors']
        condition['factors'] = factors
        
        # Extract condition_id if present
        condition_id = condition_data.get('condition_id') or condition_data.get('id')
        if condition_id:
            condition['condition_id'] = condition_id
            
        # Extract prediction_folder if present (explicit mapping)
        prediction_folder = condition_data.get('prediction_folder') or condition_data.get('folder')
        if prediction_folder:
            condition['prediction_folder'] = prediction_folder
        
        # Store condition
        self.conditions.append(condition)
        self.condition_ids[internal_id] = condition
        
        # Track source mapping
        if source_file:
            if source_file not in self.source_map:
                self.source_map[source_file] = []
            self.source_map[source_file].append(internal_id)
        
        return condition
    
    def add_conditions_from_metadata(self, metadata: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Add multiple conditions from parsed metadata.
        
        Args:
            metadata: Parsed metadata with 'conditions' list
            
        Returns:
            List of added conditions
        """
        added = []
        source_file = metadata.get('source_file')
        
        for condition_data in metadata.get('conditions', []):
            condition = self.add_condition(condition_data, source_file)
            added.append(condition)
            
        return added
    
    def get_condition_by_name(self, name: str) -> Optional[Dict[str, Any]]:
        """Find condition by name."""
        for cond in self.conditions:
            if cond['condition_name'] == name:
                return cond
        return None
    
    def get_condition_by_id(self, condition_id: str) -> Optional[Dict[str, Any]]:
        """Find condition by condition_id."""
        for cond in self.conditions:
            if cond.get('condition_id') == condition_id:
                return cond
        return None
    
    def get_condition_by_internal_id(self, internal_id: str) -> Optional[Dict[str, Any]]:
        """Find condition by internal ID."""
        return self.condition_ids.get(internal_id)
    
    def get_all_conditions(self) -> List[Dict[str, Any]]:
        """Return all conditions."""
        return self.conditions.copy()
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert registry to dictionary."""
        return {
            'metadata_version': '1.0',
            'created_at': datetime.now().isoformat(),
            'total_conditions': len(self.conditions),
            'source_files': list(self.source_map.keys()),
            'conditions': self.conditions
        }
    
    def save(self, output_path: str) -> Path:
        """Save registry to JSON file."""
        output_file = Path(output_path) / 'condition_registry.json'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
            
        return output_file
    
    def to_csv(self) -> List[Dict[str, Any]]:
        """
        Convert conditions to flat CSV-compatible format.
        
        Returns:
            List of dictionaries suitable for DataFrame conversion
        """
        rows = []
        for cond in self.conditions:
            row = {
                'internal_id': cond['internal_id'],
                'condition_name': cond['condition_name'],
                'condition_id': cond.get('condition_id', ''),
                'prediction_folder': cond.get('prediction_folder', ''),
                'source_file': cond.get('source_file', ''),
                'factors_json': json.dumps(cond.get('factors', {}))
            }
            # Add individual factor columns
            for key, value in cond.get('factors', {}).items():
                row[f'factor_{key}'] = value
            rows.append(row)
        return rows
