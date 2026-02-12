#!/usr/bin/env python3
"""
Test PFAS definitions against expected matches from comparison CSV.

This test reads the pfas_definition_comparison.csv file and validates that
each PFAS definition correctly identifies molecules as matching or not matching
according to the expected values in the CSV.
"""

import pytest
import csv
from pathlib import Path
from typing import Dict, List

from PFASgroups.PFASDefinitionModel import PFASDefinition


@pytest.fixture(scope='session')
def data_dir():
    """Get path to data directory."""
    script_dir = Path(__file__).parent
    root_dir = script_dir.parent
    return root_dir / 'PFASgroups' / 'data'


@pytest.fixture(scope='session')
def definitions(data_dir):
    """Load all PFAS definitions for testing."""
    import json
    definitions_file = data_dir / 'PFAS_definitions_smarts.json'
    
    with open(definitions_file, 'r') as f:
        definitions_data = json.load(f)
    
    definitions_dict = {}
    for data in definitions_data:
        try:
            definition = PFASDefinition(
                id=data['id'],
                name=data['name'],
                smarts=data.get('smarts', []),
                fluorineRatio=data.get('fluorineRatio'),
                description=data.get('description', ''),
                requireBoth=data.get('requireBoth', False),
                includeHydrogen=data.get('includeHydrogen', False)
            )
            definitions_dict[data['name']] = definition
        except Exception as e:
            print(f"Warning: Failed to load definition {data.get('id', '?')}: {e}")
    
    return definitions_dict


@pytest.fixture(scope='session')
def comparison_data():
    """Load comparison data from CSV file."""
    csv_file = Path(__file__).parent / 'pfas_definition_comparison.csv'
    
    test_cases = []
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            test_case = {
                'smiles': row['SMILES'],
                'description': row['Description'],
                'expected': {}
            }
            
            # Extract expected values for each definition
            for key, value in row.items():
                if key.startswith('Expected_'):
                    # Extract definition name (e.g., Expected_OECD -> OECD)
                    def_key = key.replace('Expected_', '')
                    # Convert string to boolean
                    test_case['expected'][def_key] = value.strip().lower() == 'true'
            
            test_cases.append(test_case)
    
    return test_cases


# Map CSV column names to definition names
DEFINITION_NAME_MAP = {
    'OECD': 'OECD Definition',
    'EU': 'EU PFAS Restriction',
    'OPPT': 'OPPT 2023',
    'UK': 'UK PFAS Definition',
    'PFASTRUCTv5': 'PFASTRUCTv5'
}


def test_definition_comparison_all(definitions, comparison_data):
    """Test all definitions against all molecules in the comparison CSV."""
    failures = []
    total_tests = 0
    
    for test_case in comparison_data:
        smiles = test_case['smiles']
        description = test_case['description']
        
        for def_key, expected_match in test_case['expected'].items():
            total_tests += 1
            
            # Get the full definition name
            def_name = DEFINITION_NAME_MAP.get(def_key, def_key)
            
            # Skip if definition not found
            if def_name not in definitions:
                failures.append({
                    'smiles': smiles,
                    'description': description,
                    'definition': def_name,
                    'expected': expected_match,
                    'actual': None,
                    'error': f'Definition "{def_name}" not found'
                })
                continue
            
            definition = definitions[def_name]
            
            try:
                actual_match = definition.applies_to_molecule(smiles)
                
                if actual_match != expected_match:
                    failures.append({
                        'smiles': smiles,
                        'description': description,
                        'definition': def_name,
                        'expected': expected_match,
                        'actual': actual_match,
                        'error': None
                    })
            except Exception as e:
                failures.append({
                    'smiles': smiles,
                    'description': description,
                    'definition': def_name,
                    'expected': expected_match,
                    'actual': None,
                    'error': str(e)
                })
    
    # Generate detailed failure report
    if failures:
        report = [f"\n{'='*80}"]
        report.append(f"PFAS Definition Comparison Test Failures: {len(failures)}/{total_tests}")
        report.append('='*80)
        
        for i, failure in enumerate(failures, 1):
            report.append(f"\nFailure {i}:")
            report.append(f"  SMILES: {failure['smiles']}")
            report.append(f"  Description: {failure['description']}")
            report.append(f"  Definition: {failure['definition']}")
            report.append(f"  Expected: {failure['expected']}")
            report.append(f"  Actual: {failure['actual']}")
            if failure['error']:
                report.append(f"  Error: {failure['error']}")
        
        report.append(f"\n{'='*80}\n")
        pytest.fail('\n'.join(report))


@pytest.mark.parametrize('def_key', ['OECD', 'EU', 'OPPT', 'UK', 'PFASTRUCTv5'])
def test_definition_comparison_by_definition(definitions, comparison_data, def_key):
    """Test each definition separately (allows for better test reporting)."""
    def_name = DEFINITION_NAME_MAP.get(def_key, def_key)
    
    if def_name not in definitions:
        pytest.skip(f'Definition "{def_name}" not found')
    
    definition = definitions[def_name]
    failures = []
    
    for test_case in comparison_data:
        smiles = test_case['smiles']
        description = test_case['description']
        
        # Skip if this definition is not expected in this test case
        if def_key not in test_case['expected']:
            continue
        
        expected_match = test_case['expected'][def_key]
        
        try:
            actual_match = definition.applies_to_molecule(smiles)
            
            if actual_match != expected_match:
                failures.append({
                    'smiles': smiles,
                    'description': description,
                    'expected': expected_match,
                    'actual': actual_match
                })
        except Exception as e:
            failures.append({
                'smiles': smiles,
                'description': description,
                'expected': expected_match,
                'actual': None,
                'error': str(e)
            })
    
    if failures:
        report = [f"\n{def_name} failures:"]
        for failure in failures:
            report.append(f"  {failure['smiles']} ({failure['description']})")
            report.append(f"    Expected: {failure['expected']}, Got: {failure['actual']}")
            if 'error' in failure:
                report.append(f"    Error: {failure['error']}")
        
        pytest.fail('\n'.join(report))


if __name__ == '__main__':
    # Allow running this test file directly
    pytest.main([__file__, '-v'])
