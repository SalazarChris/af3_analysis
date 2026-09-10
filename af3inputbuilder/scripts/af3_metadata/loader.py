#!/usr/bin/env python3
"""
loader.py - Metadata discovery and loading for AF3 Server predictions.
"""

import json
import os
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime


def discover_metadata_files(input_dir: str, depth: int = 3) -> List[Path]:
    """
    Search for JSON metadata files in input directory.
    
    Looks for files that might contain condition metadata:
    - *_metadata.json
    - *_conditions.json
    - *_experiment.json
    - *conditions*.json (case-insensitive)
    
    Args:
        input_dir: Root directory to search
        depth: Maximum directory depth to search
        
    Returns:
        List of paths to discovered metadata files
    """
    input_path = Path(input_dir)
    metadata_files = []
    
    search_extensions = ['.json']
    
    # Search recursively up to specified depth
    for root, dirs, files in os.walk(input_path):
        # Check current directory depth
        root_path = Path(root)
        rel_depth = len(root_path.relative_to(input_path).parts)
        if rel_depth > depth:
            continue
            
        for file in files:
            if not any(file.endswith(ext) for ext in search_extensions):
                continue
                
            # Check if filename suggests metadata
            lower_name = file.lower()
            if ('metadata' in lower_name or 
                'condition' in lower_name or 
                'experiment' in lower_name):
                metadata_files.append(root_path / file)
                
    return sorted(metadata_files)


def load_metadata(filepath: Path) -> Optional[Dict[str, Any]]:
    """
    Load and validate a metadata JSON file.
    
    Args:
        filepath: Path to metadata JSON file
        
    Returns:
        Parsed metadata dictionary or None if invalid
    """
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        # Basic validation
        if not isinstance(data, dict):
            print(f"  Warning: {filepath.name} - Root element is not a JSON object")
            return None
            
        if len(data) == 0:
            print(f"  Warning: {filepath.name} - Empty JSON object")
            return None
            
        return data
        
    except json.JSONDecodeError as e:
        print(f"  Error: {filepath.name} - Invalid JSON: {e}")
        return None
    except Exception as e:
        print(f"  Error reading {filepath.name}: {e}")
        return None


def parse_metadata_file(filepath: Path) -> Dict[str, Any]:
    """
    Parse a metadata file and extract condition information.
    
    This function handles different metadata file structures:
    - List of conditions
    - Single condition object
    - Nested structure with conditions under a key
    
    Returns:
        Dictionary with parsed conditions
    """
    data = load_metadata(filepath)
    if not data:
        return {'conditions': [], 'source_file': str(filepath)}
    
    conditions = []
    
    # Strategy 1: File is a list of conditions
    if isinstance(data, list):
        conditions = data
        
    # Strategy 2: File has a 'conditions' key
    elif isinstance(data, dict) and 'conditions' in data:
        conditions = data['conditions']
        
    # Strategy 3: File is a single condition object
    elif isinstance(data, dict) and ('condition_name' in data or 
                                      'condition_id' in data or
                                      'factors' in data):
        conditions = [data]
        
    # Strategy 4: Unknown structure - try to find condition-like objects
    else:
        # Look for any list that might contain conditions
        for key, value in data.items():
            if isinstance(value, list) and len(value) > 0:
                # Check if items look like conditions
                first_item = value[0]
                if isinstance(first_item, dict) and any(k in first_item for k in 
                                                       ['condition_name', 'condition_id', 'factors', 'prediction_folder']):
                    conditions = value
                    break
    
    return {
        'conditions': conditions,
        'source_file': str(filepath),
        'original_data': data,
        'load_timestamp': datetime.now().isoformat()
    }


def load_all_metadata(input_dir: str) -> List[Dict[str, Any]]:
    """
    Load all metadata files from input directory.
    
    Args:
        input_dir: Root directory containing metadata files
        
    Returns:
        List of parsed metadata dictionaries
    """
    metadata_files = discover_metadata_files(input_dir)
    
    all_metadata = []
    for filepath in metadata_files:
        parsed = parse_metadata_file(filepath)
        if parsed and parsed['conditions']:
            all_metadata.append(parsed)
            print(f"  Loaded: {filepath.name} ({len(parsed['conditions'])} condition(s))")
            
    return all_metadata
