#!/usr/bin/env python3
"""
Update PFAS definitions with comprehensive test compounds.

For PFASTRUCTv5: Uses test_set_for_PFASSTRUCTv5.tsv
For others: Infers appropriate test compounds based on definition criteria
"""

import json
import csv
from pathlib import Path

# Setup paths
script_dir = Path(__file__).parent
root_dir = script_dir.parent.parent
data_dir = root_dir / 'PFASgroups' / 'data'

DEFINITIONS_FILE = data_dir / 'PFAS_definitions_smarts.json'
PFASSTRUCT_TSV = script_dir.parent / 'data' / 'test_set_for_PFASSTRUCTv5.tsv'


def load_pfasstruct_tests():
    """Load test compounds from PFASTRUCTv5 TSV file."""
    positives = []
    negatives = []
    
    with open(PFASSTRUCT_TSV, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            smiles = row.get('SMILES', '').strip() if row.get('SMILES') else ''
            if not smiles:
                continue
            included = row.get('Included', '').strip().lower() == 'true'
            name = row.get('name', '').strip() if row.get('name') else ''
            
            if included:
                positives.append({'smiles': smiles, 'name': name})
            else:
                negatives.append({'smiles': smiles, 'name': name})
    
    return positives, negatives


def get_test_compounds():
    """Get test compounds for all definitions."""
    
    # OECD Definition (CF3 or CF2, not bonded to H, Cl, Br, I)
    oecd_positives = [
        {'smiles': 'C(F)(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluoropropane'},
        {'smiles': 'OC(=O)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluorobutanoic acid (PFBA)'},
        {'smiles': 'C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluorooctane (PFO)'},
        {'smiles': 'S(=O)(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluorobutane sulfonic acid (PFBS)'},
        {'smiles': 'C(F)(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'name': 'Perfluorodiethyl ether'},
        {'smiles': 'C(C(F)(F)F)(F)C(F)(F)C(F)(F)F', 'name': 'Polyfluorinated alkane'},
    ]
    oecd_negatives = [
        {'smiles': 'CCCCCC', 'name': 'Hexane (no fluorine)'},
        {'smiles': 'C(F)CC', 'name': 'Monofluoropropane (no CF2/CF3)'},
        {'smiles': 'FC(F)C', 'name': 'Contains CHF2 (excluded by OECD)'},
        {'smiles': 'C(F)(F)(Cl)C', 'name': 'Contains CFCl (excluded by OECD)'},
        {'smiles': 'CC(F)C', 'name': 'Only single CF (no CF2/CF3)'},
        {'smiles': 'C1=CC=CC=C1', 'name': 'Benzene (no fluorine)'},
    ]
    
    # EU PFAS Restriction (more strict exclusions)
    eu_positives = [
        {'smiles': 'C(F)(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluoropropane'},
        {'smiles': 'C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluorooctane'},
        {'smiles': 'S(=O)(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'PFBS'},
        {'smiles': 'C(F)(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'name': 'Perfluorodiethyl ether'},
        {'smiles': 'N(CCCC)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluoroalkylamine'},
    ]
    eu_negatives = [
        {'smiles': 'CCCCCC', 'name': 'Hexane'},
        {'smiles': 'FC(F)C', 'name': 'Contains CHF2'},
        {'smiles': 'C(F)(F)(Cl)C', 'name': 'Contains CFCl'},
        {'smiles': 'CC(F)C', 'name': 'Single CF only'},
        {'smiles': 'C1=CC=CC=C1', 'name': 'Benzene'},
        {'smiles': 'FCCC', 'name': 'Monofluoropropane'},
    ]
    
    # OPPT 2023 (Three structural patterns)
    oppt_positives = [
        {'smiles': 'C(F)(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluoropropane'},
        {'smiles': 'C(F)(F)C(F)(F)C(F)(F)C', 'name': 'Polyfluorobutane'},
        {'smiles': 'C(F)(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'name': 'Perfluoroether'},
        {'smiles': 'C(F)(F)(C(F)(F)F)C(F)(F)C(F)(F)F', 'name': 'Branched perfluoroalkane'},
        {'smiles': 'C(F)(F)(F)C(OC)(F)C(F)(F)F', 'name': 'Fluoroether'},
    ]
    oppt_negatives = [
        {'smiles': 'CCCCCC', 'name': 'Hexane'},
        {'smiles': 'C(F)CCC', 'name': 'Monofluorobutane'},
        {'smiles': 'CC(F)(F)C', 'name': 'Single CF2 (no pattern match)'},
        {'smiles': 'FC(F)C(F)F', 'name': 'Contains CHF2 (limited fluorination)'},
        {'smiles': 'C1=CC=CC=C1', 'name': 'Benzene'},
    ]
    
    # UK PFAS Definition (CF3 or CF2-CF, no Cl/Br/I)
    uk_positives = [
        {'smiles': 'C(F)(F)(F)C(F)(F)C(F)(F)F', 'name': 'Perfluoropropane'},
        {'smiles': 'C(F)(F)(F)C(F)(F)C', 'name': 'Polyfluoropropane'},
        {'smiles': 'C(F)(F)C(F)(F)C', 'name': 'Polyfluoropropane'},
        {'smiles': 'OC(=O)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'PFBA'},
        {'smiles': 'S(=O)(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'name': 'PFBS'},
        {'smiles': 'C(F)(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'name': 'Perfluoroether'},
    ]
    uk_negatives = [
        {'smiles': 'CCCCCC', 'name': 'Hexane'},
        {'smiles': 'C(F)CC', 'name': 'Monofluoropropane'},
        {'smiles': 'C(F)(F)(Cl)C', 'name': 'Contains CFCl (excluded)'},
        {'smiles': 'C(F)(F)(Br)C', 'name': 'Contains CFBr (excluded)'},
        {'smiles': 'CC(F)C', 'name': 'Single CF only'},
        {'smiles': 'C1=CC=CC=C1', 'name': 'Benzene'},
    ]
    
    # PFASTRUCTv5 (from TSV file)
    pfasstruct_pos, pfasstruct_neg = load_pfasstruct_tests()
    
    return {
        1: {'positives': oecd_positives, 'negatives': oecd_negatives},
        2: {'positives': eu_positives, 'negatives': eu_negatives},
        3: {'positives': oppt_positives, 'negatives': oppt_negatives},
        4: {'positives': uk_positives, 'negatives': uk_negatives},
        5: {'positives': pfasstruct_pos, 'negatives': pfasstruct_neg},
    }


def update_definitions():
    """Update definitions file with test compounds."""
    
    # Load definitions
    with open(DEFINITIONS_FILE, 'r') as f:
        definitions = json.load(f)
    
    # Get test compounds
    test_compounds = get_test_compounds()
    
    # Update each definition
    for definition in definitions:
        def_id = definition['id']
        if def_id in test_compounds:
            tests = test_compounds[def_id]
            
            # Update test section
            if 'test' not in definition:
                definition['test'] = {}
            
            definition['test']['category'] = 'definition'
            definition['test']['examples'] = {
                'positives': [item['smiles'] for item in tests['positives']],
                'negatives': [item['smiles'] for item in tests['negatives']]
            }
            
            print(f"Updated Definition {def_id} ({definition['name']}):")
            print(f"  True Positives: {len(tests['positives'])}")
            print(f"  True Negatives: {len(tests['negatives'])}")
    
    # Write back
    with open(DEFINITIONS_FILE, 'w') as f:
        json.dump(definitions, f, indent=2, ensure_ascii=False)
    
    print(f"\n✅ Successfully updated {DEFINITIONS_FILE}")


if __name__ == '__main__':
    update_definitions()
