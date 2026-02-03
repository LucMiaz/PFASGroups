#!/usr/bin/env python3
"""
Script to update PFAS group IDs in benchmark files and review-app lookup table
based on the ID mapping in PFAS_group_id_mapping.json
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List

def load_id_mapping(mapping_file: Path) -> Dict[int, int]:
    """Load the ID mapping from JSON file and convert keys to integers."""
    with open(mapping_file, 'r') as f:
        mapping = json.load(f)
    # Convert string keys to integers
    return {int(old_id): int(new_id) for old_id, new_id in mapping.items()}

def update_group_id(old_id: int, id_mapping: Dict[int, int]) -> int:
    """Update a single group ID using the mapping."""
    return id_mapping.get(old_id, old_id)

def update_group_ids_in_list(group_list: List[int], id_mapping: Dict[int, int]) -> List[int]:
    """Update all group IDs in a list."""
    return [update_group_id(gid, id_mapping) for gid in group_list]

def update_benchmark_result(result: Dict[str, Any], id_mapping: Dict[int, int]) -> Dict[str, Any]:
    """Update group IDs in a benchmark result object."""
    # Update detected_groups
    if 'detected_groups' in result:
        result['detected_groups'] = update_group_ids_in_list(result['detected_groups'], id_mapping)
    
    # Update matches
    if 'matches' in result:
        for match in result['matches']:
            if match.get('type') == 'group' and 'id' in match:
                old_id = match['id']
                new_id = update_group_id(old_id, id_mapping)
                match['id'] = new_id
                # Update match_id if it exists
                if 'match_id' in match:
                    match['match_id'] = f"G{new_id}"
    
    return result

def update_benchmark_file(file_path: Path, id_mapping: Dict[int, int]):
    """Update all group IDs in a benchmark JSON file."""
    print(f"Updating {file_path.name}...")
    
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Handle different structures
    if isinstance(data, list):
        for entry in data:
            if 'pfasgroups_result' in entry:
                update_benchmark_result(entry['pfasgroups_result'], id_mapping)
            elif 'expected_groups' in entry:
                entry['expected_groups'] = update_group_ids_in_list(entry['expected_groups'], id_mapping)
            if 'detected_groups' in entry:
                entry['detected_groups'] = update_group_ids_in_list(entry['detected_groups'], id_mapping)
            if 'pfasgroups_expected_groups' in entry:
                entry['pfasgroups_expected_groups'] = update_group_ids_in_list(entry['pfasgroups_expected_groups'], id_mapping)
            if 'pfasgroups_groups' in entry:
                entry['pfasgroups_groups'] = update_group_ids_in_list(entry['pfasgroups_groups'], id_mapping)
    elif isinstance(data, dict):
        if 'pfasgroups_result' in data:
            update_benchmark_result(data['pfasgroups_result'], id_mapping)
        if 'expected_groups' in data:
            data['expected_groups'] = update_group_ids_in_list(data['expected_groups'], id_mapping)
        if 'detected_groups' in data:
            data['detected_groups'] = update_group_ids_in_list(data['detected_groups'], id_mapping)
    
    # Write back
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"  ✓ Updated {file_path.name}")

def update_review_app_lookup(file_path: Path, id_mapping: Dict[int, int]):
    """Update the review-app pfas_groups_map.json lookup table."""
    print(f"Updating {file_path.name}...")
    
    with open(file_path, 'r', encoding='utf-8') as f:
        groups = json.load(f)
    
    # Create a reverse mapping to find which entries need updating
    # The lookup table has the old IDs, we need to update them to new IDs
    updated_groups = []
    for group in groups:
        old_id = group['id']
        new_id = update_group_id(old_id, id_mapping)
        group['id'] = new_id
        updated_groups.append(group)
    
    # Sort by new ID
    updated_groups.sort(key=lambda x: x['id'])
    
    # Write back
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(updated_groups, f, indent=2, ensure_ascii=False)
    
    print(f"  ✓ Updated {file_path.name}")

def main():
    # Set up paths
    benchmark_dir = Path(__file__).parent
    mapping_file = benchmark_dir / 'data' / 'PFAS_group_id_mapping.json'
    data_dir = benchmark_dir / 'data'
    scripts_data_dir = benchmark_dir / 'scripts' / 'data'
    review_app_dir = benchmark_dir / 'review-app'
    
    # Load ID mapping
    print("Loading ID mapping...")
    id_mapping = load_id_mapping(mapping_file)
    print(f"  Loaded {len(id_mapping)} ID mappings")
    
    # Update benchmark files in data/
    print("\nUpdating benchmark files in data/...")
    benchmark_files = [
        'pfas_oecd_benchmark_20260202_104629.json',
        'pfas_timing_benchmark_20260202_105102.json',
        'pfas_non_fluorinated_benchmark_20260202_105102.json',
        'pfas_highly_branched_benchmark_20260202_105121.json',
        'pfas_enhanced_benchmark_20260202_104204.json',
        'pfas_complex_branched_benchmark_20260202_105117.json',
        'telomer_validation_results.json'
    ]
    
    for filename in benchmark_files:
        file_path = data_dir / filename
        if file_path.exists():
            update_benchmark_file(file_path, id_mapping)
        else:
            print(f"  ⚠ File not found: {filename}")
    
    # Update benchmark files in scripts/data/
    print("\nUpdating benchmark files in scripts/data/...")
    if scripts_data_dir.exists():
        for file_path in scripts_data_dir.glob('pfas_definitions_benchmark_*.json'):
            update_benchmark_file(file_path, id_mapping)
    
    # Update review-app lookup table
    print("\nUpdating review-app lookup table...")
    review_app_lookup = review_app_dir / 'pfas_groups_map.json'
    if review_app_lookup.exists():
        update_review_app_lookup(review_app_lookup, id_mapping)
    else:
        print(f"  ⚠ File not found: pfas_groups_map.json")
    
    print("\n✓ All updates completed!")

if __name__ == '__main__':
    main()
