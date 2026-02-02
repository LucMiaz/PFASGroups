#!/usr/bin/env python3
"""
Add test metadata to PFAS groups and definitions JSON files.

This script adds "test" keys to each entry with:
- For PFAS groups: category (OECD/generic/telomer), examples, and generation info
- For PFAS definitions: category (definition), examples with TP/TN/FP/FN
"""

import json
import csv
import sys
import os
from pathlib import Path
from collections import defaultdict
import random

# Setup paths
script_dir = Path(__file__).parent
root_dir = script_dir.parent.parent
data_dir = script_dir.parent / 'data'
pfas_data_dir = root_dir / 'PFASgroups' / 'data'

# File paths
GROUPS_FILE = pfas_data_dir / 'PFAS_groups_smarts.json'
DEFINITIONS_FILE = pfas_data_dir / 'PFAS_definitions_smarts.json'
OECD_BENCHMARK = data_dir / 'pfas_oecd_benchmark_20260202_104629.json'
TELOMER_VALIDATION = data_dir / 'telomer_validation_results.json'
TEST_COMPOUNDS = script_dir / 'benchmark_test_compounds.csv'
ENHANCED_BENCHMARK = script_dir / 'enhanced_pfas_benchmark.py'


def load_smiles_map():
    """Extract smiles_map from enhanced_pfas_benchmark.py"""
    smiles_map = {}
    
    with open(ENHANCED_BENCHMARK, 'r') as f:
        content = f.read()
        
    # Find the smiles_map dictionary definition
    import re
    pattern = r'smiles_map\s*=\s*\{(.*?)\}'
    match = re.search(pattern, content, re.DOTALL)
    
    if match:
        # Parse the dictionary manually
        map_content = match.group(1)
        # Extract each group entry
        group_pattern = r'(\d+):\s*\{[^}]*\'smiles\':\s*[\'"]([^\'"]+)[\'"][^}]*\'mode\':\s*[\'"]([^\'"]+)[\'"][^}]*\}'
        
        for group_match in re.finditer(group_pattern, map_content):
            group_id = int(group_match.group(1))
            smiles = group_match.group(2)
            mode = group_match.group(3)
            smiles_map[group_id] = {'smiles': smiles, 'mode': mode}
    
    return smiles_map


def load_oecd_examples():
    """Load examples from OECD benchmark"""
    with open(OECD_BENCHMARK, 'r') as f:
        oecd_data = json.load(f)
    
    # Group by detected groups
    group_examples = defaultdict(list)
    
    for molecule in oecd_data:
        result = molecule.get('pfasgroups_result', {})
        if result.get('is_pfas'):
            groups = result.get('groups_detected', [])
            smiles = molecule.get('smiles', '')
            name = molecule.get('name', '')
            
            for group_id in groups:
                if len(group_examples[group_id]) < 5:  # Collect up to 5 examples
                    group_examples[group_id].append({
                        'smiles': smiles,
                        'name': name
                    })
    
    return group_examples


def load_telomer_examples():
    """Load telomer examples from validation results"""
    with open(TELOMER_VALIDATION, 'r') as f:
        telomer_data = json.load(f)
    
    telomer_examples = defaultdict(list)
    
    results = telomer_data.get('results', [])
    for result in results:
        if result.get('detected', False):
            groups = result.get('groups_detected', [])
            smiles = result.get('smiles', '')
            cid = result.get('cid', '')
            
            for group_id in groups:
                if len(telomer_examples[group_id]) < 5:
                    telomer_examples[group_id].append({
                        'smiles': smiles,
                        'cid': cid
                    })
    
    return telomer_examples


def load_test_compounds():
    """Load test compounds for definitions"""
    compounds = {
        'true_positives': [],
        'true_negatives': [],
        'false_positives': [],
        'false_negatives': []
    }
    
    with open(TEST_COMPOUNDS, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            smiles = row.get('smiles', '')
            expected = row.get('expected_pfas', '').lower() == 'true'
            category = row.get('category', '')
            
            # Store a subset for each category
            compound_info = {
                'smiles': smiles,
                'category': category
            }
            
            # This is a simplified version - would need actual benchmark results
            # to determine TP/TN/FP/FN for each definition
            if expected:
                if len(compounds['true_positives']) < 10:
                    compounds['true_positives'].append(compound_info)
            else:
                if len(compounds['true_negatives']) < 10:
                    compounds['true_negatives'].append(compound_info)
    
    return compounds


def categorize_groups(groups_data):
    """Categorize groups as OECD, generic, or telomer"""
    categorized = {}
    
    for group in groups_data:
        group_id = group['id']
        name = group['name'].lower()
        
        # Telomer groups (40-48)
        if 40 <= group_id <= 48:
            category = 'telomer'
        # Generic component groups (49-51)
        elif group_id in [49, 50, 51]:
            category = 'generic'
        # OECD groups (1-39)
        else:
            category = 'OECD'
        
        categorized[group_id] = category
    
    return categorized


def add_test_metadata_to_groups():
    """Add test metadata to PFAS groups"""
    print("📝 Adding test metadata to PFAS groups...")
    
    # Load data
    with open(GROUPS_FILE, 'r') as f:
        groups_data = json.load(f)
    
    smiles_map = load_smiles_map()
    oecd_examples = load_oecd_examples()
    telomer_examples = load_telomer_examples()
    categories = categorize_groups(groups_data)
    
    # Add test metadata
    for group in groups_data:
        group_id = group['id']
        category = categories.get(group_id, 'OECD')
        
        test_data = {'category': category}
        
        if category == 'telomer':
            # Telomer groups
            test_data['generate'] = {
                'smiles': smiles_map.get(group_id, {}).get('smiles', ''),
                'mode': smiles_map.get(group_id, {}).get('mode', 'insert'),
                'is_telomer': True
            }
            test_data['examples'] = telomer_examples.get(group_id, [])[:3]
            
        elif category == 'generic':
            # Generic groups
            test_data['generate'] = {
                'smiles': smiles_map.get(group_id, {}).get('smiles', ''),
                'mode': smiles_map.get(group_id, {}).get('mode', 'insert'),
                'is_telomer': False
            }
            test_data['examples'] = oecd_examples.get(group_id, [])[:3]
            
        else:  # OECD
            # OECD groups - just examples
            examples = oecd_examples.get(group_id, [])[:4]
            test_data['examples'] = examples
        
        group['test'] = test_data
    
    # Save
    with open(GROUPS_FILE, 'w') as f:
        json.dump(groups_data, f, indent=2)
    
    print(f"✅ Updated {len(groups_data)} groups in {GROUPS_FILE}")


def add_test_metadata_to_definitions():
    """Add test metadata to PFAS definitions"""
    print("\n📝 Adding test metadata to PFAS definitions...")
    
    # Load data
    with open(DEFINITIONS_FILE, 'r') as f:
        definitions_data = json.load(f)
    
    test_compounds = load_test_compounds()
    
    # Add test metadata
    for definition in definitions_data:
        test_data = {
            'category': 'definition',
            'examples': {
                'true_positives': test_compounds['true_positives'][:3],
                'true_negatives': test_compounds['true_negatives'][:3],
                'note': 'Examples from benchmark test suite'
            }
        }
        
        definition['test'] = test_data
    
    # Save
    with open(DEFINITIONS_FILE, 'w') as f:
        json.dump(definitions_data, f, indent=2)
    
    print(f"✅ Updated {len(definitions_data)} definitions in {DEFINITIONS_FILE}")


def main():
    """Main function"""
    print("="*60)
    print("ADDING TEST METADATA TO PFAS GROUPS AND DEFINITIONS")
    print("="*60)
    
    # Check files exist
    for filepath in [GROUPS_FILE, DEFINITIONS_FILE, OECD_BENCHMARK, 
                     TELOMER_VALIDATION, TEST_COMPOUNDS, ENHANCED_BENCHMARK]:
        if not Path(filepath).exists():
            print(f"❌ Error: File not found: {filepath}")
            return 1
    
    try:
        add_test_metadata_to_groups()
        add_test_metadata_to_definitions()
        
        print("\n" + "="*60)
        print("✅ TEST METADATA SUCCESSFULLY ADDED")
        print("="*60)
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
