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
    smiles_map = {
        # Keep existing mappings for groups 29-59
        29: {'smiles': 'O[H]', 'mode': 'attach'},
        30: {'smiles': 'C(=O)', 'mode': 'insert'},
        31: {'smiles': 'O', 'mode': 'insert'},
        32: {'smiles': 'C(=O)OC', 'mode': 'insert'},
        33: {'smiles': 'C(=O)O', 'mode': 'attach'},
        34: {'smiles': 'C(=O)N', 'mode': 'insert'},
        35: {'smiles': 'C(=O)Cl', 'mode': 'attach'},
        36: {'smiles': 'S(=O)(=O)O', 'mode': 'attach'},
        37: {'smiles': 'SO[H]', 'mode': 'attach'},
        38: {'smiles': 'S(=O)O[H]', 'mode': 'attach'},
        39: {'smiles': 'OS(=O)(=O)O', 'mode': 'attach'},
        40: {'smiles': 'P(=O)(O)O', 'mode': 'attach'},
        41: {'smiles': 'P(=O)O', 'mode': 'attach'},
        42: {'smiles': 'Cl', 'mode': 'attach'},
        43: {'smiles': 'Br', 'mode': 'attach'},
        44: {'smiles': 'I', 'mode': 'attach'},
        45: {'smiles': 'S(=O)(=O)N', 'mode': 'insert'},
        46: {'smiles': 'c1ncc[nH]1', 'mode': 'attach'},
        47: {'smiles': 'c1ncccc1', 'mode': 'attach'},
        48: {'smiles': 'c1ccc2OC(F)(F)Oc2c1', 'mode': 'attach'},
        49: {'smiles': 'N', 'mode': 'insert'},
        52: {'smiles': 'C(F)=C(F)(F)F', 'mode': 'attach'},
        53: {'smiles': 'C#C', 'mode': 'insert'},
        54: {'smiles': 'c1ccccc1', 'mode': 'attach'},
        55: {'smiles': 'C1(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)(F)', 'mode': 'attach'},
        56: {'smiles': 'C1(F)C(F)(F)(F)C(H)C(F)C(F)(F)C1(F)', 'mode': 'attach'},
        57: {'smiles': 'c1c(F)c(F)c(F)c(F)c1F', 'mode': 'attach'},
        58: {'smiles': 'c1c(F)c(F)ccc1F', 'mode': 'attach'},
        59: {'smiles': 'OO', 'mode': 'insert'},#peroxyde
        60: {'smiles': 'C(=O)OOC(=O)', 'mode': 'insert'},# Benzoyl peroxides
        61: {'smiles': '[SiH3]', 'mode': 'attach'},  # Silane - needs SiH3 not Si(C)(C)C
        62: {'smiles': '[Si](Cl)(Cl)Cl', 'mode': 'attach'},  # Trichlorosilane
        63: {'smiles': 'COC(=O)C=C', 'mode': 'attach'},  # Acrylate - needs CH2 anchor
        64: {'smiles': 'COC(=O)C(=C)C', 'mode': 'attach'},  # Methacrylate - needs CH2 anchor
        65: {'smiles': 'N(C)(C)CC(=O)O', 'mode': 'attach'},  # Betaine
        66: {'smiles': 'SN#C', 'mode': 'attach'},  # Thiocyanic acid
        67: {'smiles': 'C(=O)SC(=O)C(O)C(=O)O', 'mode': 'attach'},  # Thia keto propanoic acid - max_dist_from_CF:0, attaches directly
        68: {'smiles': 'OC1C(O)C(OC(C1O)C(=O)O)O', 'mode': 'attach'},  # Glucuronate
        # Telomer groups 69+ - use simplified patterns for generation
        69: {'smiles': 'CCCC[Si]([H])([H])[H]', 'mode': 'attach'},  # FT silane
        70: {'smiles': '[Si](Cl)(Cl)Cl', 'mode': 'attach'},  # FT trichlorosilane
        71: {'smiles': 'CCCCI', 'mode': 'attach'},  # FT iodide
        72: {'smiles': 'CCCCC(=O)[H]', 'mode': 'attach'},  # FT aldehyde
        73: {'smiles': 'CCCCC(=O)O', 'mode': 'attach'},  # FT carboxylic acids
        74: {'smiles': 'CCCCOCCCCO', 'mode': 'attach'},  # FT ethoxylates
        75: {'smiles': 'CCCCS(=O)(=O)O', 'mode': 'attach'},  # FT sulfonic acid
        76: {'smiles': 'CCCCOP(=O)(O)O', 'mode': 'attach'},  # FT monophosphate
        77: {'smiles': 'CCCCOP(=O)(O)OP(=O)(O)O', 'mode': 'attach'},  # FT diphosphate
        78: {'smiles': 'CCCCOP(=O)(O)OP(=O)(O)OP(=O)(O)O', 'mode': 'attach'},  # FT triphosphate
        79: {'smiles': 'CCCCOC(=O)C=C', 'mode': 'attach'},  # FT acrylate
        80: {'smiles': 'CCCCOC(=O)C(=C)C', 'mode': 'attach'},  # FT methacrylate
        81: {'smiles': 'CCCCO', 'mode': 'attach'},  # FT alcohol
        82: {'smiles': 'CCCCS', 'mode': 'insert'},  # FT sulfure
        83: {'smiles': 'CCCCS(=O)(=O)N', 'mode': 'insert'},  # FT sulfonamide
        84: {'smiles': 'CCCC[Si](OC)(OC)OC', 'mode': 'attach'},  # FT silyl
        85: {'smiles': 'CCCCN(C)(C)CC(=O)O', 'mode': 'attach'},  # FT betaine
        86: {'smiles': 'CCCCOC(=O)C', 'mode': 'attach'},  # FT ester
        87: {'smiles': 'CCCCSN#C', 'mode': 'attach'},  # FT thiocyanic acid
        88: {'smiles': 'CCCCC=C', 'mode': 'insert'},  # FT alkene
        89: {'smiles': 'CCCCN(C)(C)C', 'mode': 'attach'},  # FT trimethylamine
        90: {'smiles': 'CCCCS(=O)O', 'mode': 'attach'},  # FT sulfinic acid
        91: {'smiles': 'CCCCSO', 'mode': 'attach'},  # FT sulfenic acid
        92: {'smiles': 'OS(=O)(=O)O', 'mode': 'attach'},  # FT sulfuric acid
        93: {'smiles': 'CCCCC', 'mode': 'attach'},  # FT methyl
        94: {'smiles': 'CCCCOC1C(O)C(OC(C1O)C(=O)O)O', 'mode': 'attach'},  # FT glucuronic acid
        95: {'smiles': 'CCCCC(O)S(=O)(=O)O', 'mode': 'attach'},  # FT formaldehyde bisulfite
        96: {'smiles': 'CCCCC(=O)O', 'mode': 'attach'},  # FT unsaturated carboxylic acids (with alkene)
        97: {'smiles': 'CCCCC(=O)SC(=O)C(O)C(=O)O', 'mode': 'attach'},  # FT thia keto propanoic acid
        # Non-telomer groups 98+
        98: {'smiles': 'C[NH2+]CC(=O)O', 'mode': 'attach'},  # Glycine - needs charged nitrogen or N+ form
        # More telomers 99+
        99: {'smiles': 'CCCCNCC(=O)O', 'mode': 'attach'},  # FT glycine
        100: {'smiles': 'CS(=O)(=O)NCCO', 'mode': 'attach'},  # Sulfonamidoethanol - simpler anchor
        101: {'smiles': 'S(=O)(=O)NCCO', 'mode': 'attach'},  # FT sulfonamidoethanol
        102: {'smiles': 'CCCCOCC(=O)O', 'mode': 'attach'},  # FT ether carboxylic acids
        103: {'smiles': 'CCCCS(=O)(=O)CCC(=O)O', 'mode': 'attach'},  # FT sulfonyl propanoic acid
        104: {'smiles': 'S(=O)(=O)CCC(=O)O', 'mode': 'attach'},  # Sulfonyl propanoic acid - max_dist_from_CF:0, direct attachment
        105: {'smiles': 'CCCCS(=O)CCC(=O)NC(C)(C)CS(=O)(=O)O', 'mode': 'attach'},  # FT sulfinyl amido sulfonic acid
        106: {'smiles': 'S(=O)(CCC(=O)NC(C)(C)CS(=O)(=O)O)', 'mode': 'attach'},  # Sulfinyl amido sulfonic acid - max_dist_from_CF:0, direct attachment
        107: {'smiles': 'CCCCS(=O)', 'mode': 'attach'},  # FT sulfinyl
        108: {'smiles': 'CCCCS(=O)(=O)', 'mode': 'attach'},  # FT sulfone
        109: {'smiles': 'CCCCC(=O)SCC(O)C(=O)O', 'mode': 'attach'},  # FT acetylsulfanylhydroxypropanoic acid
        110: {'smiles': 'CCCCC(=O)SC(O)C(=O)O', 'mode': 'attach'},  # FT acetylhydroxyethanethioate
        111: {'smiles': 'CCCCNCCC(=O)O', 'mode': 'attach'},  # FT amino propanoic acid
        112: {'smiles': 'CCCCNCC[N+](C)(C)C', 'mode': 'attach'},  # FT amino ethyl trimethyl ammonium
    }
    return smiles_map


def load_oecd_examples():
    """Load examples from OECD benchmark"""
    with open(OECD_BENCHMARK, 'r') as f:
        oecd_data = json.load(f)
    
    # Group by detected groups
    group_examples = {}
    
    for molecule in oecd_data:
        result = molecule.get('pfasgroups_result', {})
        groups = result.get('detected_groups', [])
        smiles = molecule.get('molecule_data',{}).get('smiles', '')
        if smiles !='':
            if groups is not None:
                for group_id in groups:
                    group_examples.setdefault(group_id,[]).append(
                        smiles)
    return group_examples


def load_telomer_examples():
    """Load telomer examples from validation results"""
    with open(TELOMER_VALIDATION, 'r') as f:
        telomer_data = json.load(f)
    
    telomer_examples = {}
    results = telomer_data.get('results', [])
    for result in results:
        if result.get('detected', False):
            groups = result.get('telomer_groups', [])
            smiles = result.get('smiles', '')
            
            for group_id in groups:
                telomer_examples.setdefault(group_id.get('id'),[]).append(
                        smiles
                    )
    
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
        if 'telomer' in name.lower():
            category = 'telomer'
        # Generic component groups (49-51)
        elif group_id < 29:
            category = 'OECD'
        # OECD groups (1-39)
        else:
            category = 'generic'
        
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
        category = categories.get(group_id, 'generic')
        
        test_data = {'category': category}
        
        if category == 'telomer':
            # Telomer groups
            test_data['generate'] = {
                'smiles': smiles_map.get(group_id, {}).get('smiles', ''),
                'mode': smiles_map.get(group_id, {}).get('mode', 'insert'),
                'is_telomer': True
            }
            test_data['examples'] = telomer_examples.get(group_id, [])[:5]
            
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
        #add_test_metadata_to_definitions()
        
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
