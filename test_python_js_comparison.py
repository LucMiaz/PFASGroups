"""
Test script to compare Python and JavaScript PFAS group analysis results.
Tests both bycomponent=True and bycomponent=False flavors.
"""
import json
import requests
from rdkit import Chem
import sys
import os

# Add parent directory to path to import PFASgroups
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PFASgroups.core import parse_groups_in_mol

# Test molecules covering different PFAS group types
TEST_MOLECULES = [
    {
        'name': 'PFOA (Perfluorooctanoic acid)',
        'smiles': 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
        'expected_groups': [1]  # Perfluoroalkyl carboxylic acids
    },
    {
        'name': 'PFOS (Perfluorooctane sulfonic acid)',
        'smiles': 'C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
        'expected_groups': [4]  # Perfluoroalkane sulfonic acids
    },
    {
        'name': '6:2 FTOH (6:2 Fluorotelomer alcohol)',
        'smiles': 'C(CCO)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
        'expected_groups': [7]  # Fluorotelomer alcohols
    },
    {
        'name': 'N-EtFOSE (N-Ethyl perfluorooctane sulfonamidoethanol)',
        'smiles': 'C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)NCCO)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
        'expected_groups': [13]  # Perfluoroalkane sulfonamides
    },
    {
        'name': 'Simple PFAS with short chain',
        'smiles': 'FC(F)(F)C(F)(F)C(=O)O',
        'expected_groups': [1]  # Perfluoroalkyl carboxylic acids
    },
    {
        'name': 'PFBA (Perfluorobutanoic acid)',
        'smiles': 'C(=O)(O)C(C(C(F)(F)F)(F)F)(F)F',
        'expected_groups': [1]
    },
    {
        'name': 'PFBS (Perfluorobutane sulfonic acid)',
        'smiles': 'C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F',
        'expected_groups': [4]
    },
]

def normalize_group_result(group_data):
    """Normalize group result for comparison."""
    return {
        'id': group_data.get('id'),
        'name': group_data.get('name'),
        'matchCount': group_data.get('matchCount', group_data.get('match_count', 0)),
        'nCFchain': group_data.get('nCFchain', group_data.get('n_CFchain', 0)),
    }

def compare_results(python_result, js_result, molecule_name, flavor):
    """Compare Python and JavaScript results for a molecule."""
    print(f"\n{'='*80}")
    print(f"Testing: {molecule_name} (bycomponent={flavor})")
    print(f"{'='*80}")
    
    # Extract groups from Python result
    py_groups = {}
    if python_result and 'detected_groups' in python_result:
        for group in python_result['detected_groups']:
            py_groups[group['id']] = normalize_group_result(group)
    
    # Extract groups from JS result
    js_groups = {}
    if js_result and 'groups' in js_result:
        for group in js_result['groups']:
            # For JS, extract the appropriate flavor
            group_data = None
            if flavor:
                group_data = group.get('bycomponent')
            else:
                group_data = group.get('default')
            
            if group_data:
                js_groups[group['id']] = {
                    'id': group['id'],
                    'name': group['name'],
                    'matchCount': group_data.get('matchCount', 0),
                    'nCFchain': group_data.get('nCFchain', 0),
                }
    
    # Compare detected group IDs
    py_ids = set(py_groups.keys())
    js_ids = set(js_groups.keys())
    
    print(f"\nPython detected groups: {sorted(py_ids)}")
    print(f"JavaScript detected groups: {sorted(js_ids)}")
    
    # Check for differences
    only_in_python = py_ids - js_ids
    only_in_js = js_ids - py_ids
    in_both = py_ids & js_ids
    
    match = True
    
    if only_in_python:
        print(f"\n❌ Groups detected ONLY in Python: {sorted(only_in_python)}")
        match = False
    
    if only_in_js:
        print(f"\n❌ Groups detected ONLY in JavaScript: {sorted(only_in_js)}")
        match = False
    
    # Compare details for groups detected in both
    if in_both:
        print(f"\n✓ Groups detected in BOTH: {sorted(in_both)}")
        
        for group_id in sorted(in_both):
            py_g = py_groups[group_id]
            js_g = js_groups[group_id]
            
            print(f"\n  Group {group_id} ({py_g['name']}):")
            print(f"    Python    - matchCount: {py_g['matchCount']}, nCFchain: {py_g['nCFchain']}")
            print(f"    JavaScript - matchCount: {js_g['matchCount']}, nCFchain: {js_g['nCFchain']}")
            
            if py_g['matchCount'] != js_g['matchCount']:
                print(f"    ⚠️  matchCount differs!")
                match = False
            
            if py_g['nCFchain'] != js_g['nCFchain']:
                print(f"    ⚠️  nCFchain differs!")
                match = False
    
    if match and py_ids == js_ids:
        print(f"\n✅ PERFECT MATCH - Python and JavaScript produce identical results")
        return True
    else:
        print(f"\n❌ MISMATCH DETECTED")
        return False

def test_molecule(smiles, name, bycomponent=False):
    """Test a single molecule with both Python and JavaScript."""
    # Python analysis
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"⚠️  Failed to parse SMILES in Python: {smiles}")
            return None, None
        
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        py_result = parse_groups_in_mol(mol, formula, bycomponent=bycomponent)
    except Exception as e:
        print(f"⚠️  Python analysis failed: {e}")
        return None, None
    
    # JavaScript analysis
    try:
        response = requests.post(
            'http://localhost:3000/analyze',
            json={
                'input': smiles,
                'inputType': 'smiles',
                'formula': None
            },
            timeout=10
        )
        
        if response.status_code != 200:
            print(f"⚠️  JavaScript API returned status {response.status_code}")
            return py_result, None
        
        js_result = response.json()
    except requests.exceptions.ConnectionError:
        print(f"⚠️  Cannot connect to JavaScript server at http://localhost:3000")
        print(f"    Please start the server with: node server.js")
        return py_result, None
    except Exception as e:
        print(f"⚠️  JavaScript analysis failed: {e}")
        return py_result, None
    
    return py_result, js_result

def main():
    """Run comparison tests."""
    print("="*80)
    print("PFAS Groups Analysis Comparison: Python vs JavaScript")
    print("="*80)
    
    # Check if JS server is running
    try:
        response = requests.get('http://localhost:3000/health', timeout=5)
        if response.status_code == 200:
            health = response.json()
            print(f"\n✓ JavaScript server is running")
            print(f"  RDKit status: {health.get('rdkit', 'unknown')}")
            print(f"  RDKit version: {health.get('version', 'unknown')}")
        else:
            print(f"\n⚠️  JavaScript server returned unexpected status: {response.status_code}")
            return
    except requests.exceptions.ConnectionError:
        print(f"\n❌ Cannot connect to JavaScript server at http://localhost:3000")
        print(f"   Please start the server with: cd PFASgroupsJS && node server.js")
        return
    
    results_summary = {
        'total': 0,
        'perfect_match': 0,
        'mismatch': 0,
        'errors': 0
    }
    
    # Test both flavors
    for flavor in [False, True]:
        flavor_name = "bycomponent=True" if flavor else "bycomponent=False (default)"
        print(f"\n\n{'#'*80}")
        print(f"# Testing with {flavor_name}")
        print(f"{'#'*80}")
        
        for test_case in TEST_MOLECULES:
            results_summary['total'] += 1
            
            py_result, js_result = test_molecule(
                test_case['smiles'],
                test_case['name'],
                bycomponent=flavor
            )
            
            if py_result is None or js_result is None:
                results_summary['errors'] += 1
                continue
            
            match = compare_results(py_result, js_result, test_case['name'], flavor)
            
            if match:
                results_summary['perfect_match'] += 1
            else:
                results_summary['mismatch'] += 1
    
    # Print summary
    print(f"\n\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Total tests: {results_summary['total']}")
    print(f"✅ Perfect matches: {results_summary['perfect_match']}")
    print(f"❌ Mismatches: {results_summary['mismatch']}")
    print(f"⚠️  Errors: {results_summary['errors']}")
    
    if results_summary['mismatch'] == 0 and results_summary['errors'] == 0:
        print(f"\n🎉 ALL TESTS PASSED - Python and JavaScript implementations are identical!")
        return 0
    else:
        print(f"\n⚠️  Some tests failed - implementations may differ")
        return 1

if __name__ == '__main__':
    sys.exit(main())
