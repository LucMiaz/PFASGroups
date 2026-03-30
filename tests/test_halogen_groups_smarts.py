#!/usr/bin/env python3
"""
Pytest-compatible test suite for HalogenGroupoups and PFASDefinitions.

This module provides pytest test functions for validating PFAS groups and definitions.
Can be run with pytest or the standalone run_groups_definitions_tests.py script.

Usage:
    # Run with pytest
    pytest tests/test_pytest.py -v
    
    # Run specific test
    pytest tests/test_pytest.py::test_all_groups -v
    
    # Run with specific markers
    pytest tests/test_pytest.py -m groups
    pytest tests/test_pytest.py -m definitions
"""

import pytest
import json
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from HalogenGroups import HalogenGroup, PFASDefinition


@pytest.fixture(scope="module")
def groups_file():
    """Path to Halogen groups JSON file."""
    return Path(__file__).parent.parent / 'PFASGroups' / 'data' / 'Halogen_groups_smarts.json'


@pytest.fixture(scope="module")
def definitions_file():
    """Path to PFAS definitions JSON file."""
    return Path(__file__).parent.parent / 'PFASGroups' / 'data' / 'PFAS_definitions_smarts.json'


@pytest.fixture(scope="module")
def all_groups(groups_file):
    """Load all Halogen groups from JSON."""
    with open(groups_file, 'r') as f:
        groups_data = json.load(f)
    
    groups = []
    for data in groups_data:
        try:
            group = HalogenGroup(**data)
            groups.append(group)
        except Exception as e:
            pytest.fail(f"Failed to load group {data.get('id', '?')}: {e}")
    
    return groups


@pytest.fixture(scope="module")
def all_definitions(definitions_file):
    """Load all PFAS definitions from JSON."""
    with open(definitions_file, 'r') as f:
        definitions_data = json.load(f)
    
    definitions = []
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
            definitions.append(definition)
        except Exception as e:
            pytest.fail(f"Failed to load definition {data.get('id', '?')}: {e}")
    
    return definitions


@pytest.mark.groups
@pytest.mark.parametrize("group_id", range(1, 120))  # Test groups 1-119
def test_individual_group(all_groups, group_id):
    """Test individual Halogen group by ID."""
    # Find the group
    group = next((g for g in all_groups if g.id == group_id), None)
    
    if group is None:
        pytest.skip(f"Group {group_id} not found")
    
    # Run test
    result = group.test()
    
    # Skip if no test data
    if result.get('passed') is None:
        pytest.skip(f"No test data for group {group_id}: {group.name}")
    
    # Check for failures
    if not result['passed']:
        failure_msg = f"Group {group_id} ({group.name}) failed {len(result['failures'])}/{result['total_tests']} tests:\n"
        for failure in result['failures'][:5]:  # Show first 5 failures
            failure_msg += f"  - {failure['smiles'][:50]}...\n"
            failure_msg += f"    Expected: {failure['expected']}, Got: {failure['got']}\n"
            failure_msg += f"    Error: {failure['error']}\n"
        pytest.fail(failure_msg)
    
    # Assert passed
    assert result['passed'], f"Group {group_id} ({group.name}) failed"


@pytest.mark.groups
def test_all_groups(all_groups):
    """Test all Halogen groups at once."""
    failed_groups = []
    
    for group in all_groups:
        result = group.test()
        
        # Skip if no test data
        if result.get('passed') is None:
            continue
        
        if not result['passed']:
            failed_groups.append({
                'id': group.id,
                'name': group.name,
                'failures': result['failures']
            })
    
    if failed_groups:
        msg = f"{len(failed_groups)} groups failed:\n"
        for fg in failed_groups[:10]:  # Show first 10 failures
            msg += f"  - Group {fg['id']}: {fg['name']} ({len(fg['failures'])} failures)\n"
        pytest.fail(msg)
    
    assert len(failed_groups) == 0


@pytest.mark.definitions
@pytest.mark.parametrize("definition_id", range(1, 6))  # Test definitions 1-5
def test_individual_definition(all_definitions, definition_id):
    """Test individual PFAS definition by ID."""
    # Find the definition
    definition = next((d for d in all_definitions if d.id == definition_id), None)
    
    if definition is None:
        pytest.skip(f"Definition {definition_id} not found")
    
    # Run test
    result = definition.test()
    
    # Skip if no test data
    if result.get('passed') is None:
        pytest.skip(f"No test data for definition {definition_id}: {definition.name}")
    
    # Check for failures
    if not result['passed']:
        stats = result.get('stats', {})
        failure_msg = f"Definition {definition_id} ({definition.name}) failed:\n"
        failure_msg += f"  TP: {stats.get('true_positives', 0)}, "
        failure_msg += f"TN: {stats.get('true_negatives', 0)}, "
        failure_msg += f"FP: {stats.get('false_positives', 0)}, "
        failure_msg += f"FN: {stats.get('false_negatives', 0)}\n"
        
        for failure in result['failures'][:5]:  # Show first 5 failures
            failure_msg += f"  - [{failure.get('type', 'unknown')}] {failure['smiles'][:50]}...\n"
            failure_msg += f"    Expected: {failure['expected']}, Got: {failure['got']}\n"
            failure_msg += f"    Error: {failure['error']}\n"
        pytest.fail(failure_msg)
    
    # Assert passed
    assert result['passed'], f"Definition {definition_id} ({definition.name}) failed"


@pytest.mark.definitions
def test_all_definitions(all_definitions):
    """Test all PFAS definitions at once."""
    failed_definitions = []
    
    for definition in all_definitions:
        result = definition.test()
        
        # Skip if no test data
        if result.get('passed') is None:
            continue
        
        if not result['passed']:
            failed_definitions.append({
                'id': definition.id,
                'name': definition.name,
                'stats': result.get('stats', {}),
                'failures': result['failures']
            })
    
    if failed_definitions:
        msg = f"{len(failed_definitions)} definitions failed:\n"
        for fd in failed_definitions:
            stats = fd['stats']
            msg += f"  - Definition {fd['id']}: {fd['name']}\n"
            msg += f"    TP: {stats.get('true_positives', 0)}, "
            msg += f"TN: {stats.get('true_negatives', 0)}, "
            msg += f"FP: {stats.get('false_positives', 0)}, "
            msg += f"FN: {stats.get('false_negatives', 0)}\n"
        pytest.fail(msg)
    
    assert len(failed_definitions) == 0


@pytest.mark.smoke
def test_can_load_groups(groups_file):
    """Smoke test: Can load groups file."""
    assert groups_file.exists(), f"Groups file not found: {groups_file}"
    
    with open(groups_file, 'r') as f:
        data = json.load(f)
    
    assert isinstance(data, list), "Groups data should be a list"
    assert len(data) > 0, "Groups data should not be empty"


@pytest.mark.smoke
def test_can_load_definitions(definitions_file):
    """Smoke test: Can load definitions file."""
    assert definitions_file.exists(), f"Definitions file not found: {definitions_file}"
    
    with open(definitions_file, 'r') as f:
        data = json.load(f)
    
    assert isinstance(data, list), "Definitions data should be a list"
    assert len(data) > 0, "Definitions data should not be empty"


@pytest.mark.smoke
def test_group_model_has_test_method():
    """Smoke test: HalogenGroup has test() method."""
    assert hasattr(HalogenGroup, 'test'), "HalogenGroup should have test() method"


@pytest.mark.smoke
def test_definition_model_has_test_method():
    """Smoke test: PFASDefinition has test() method."""
    assert hasattr(PFASDefinition, 'test'), "PFASDefinition should have test() method"


if __name__ == '__main__':
    # Run with pytest if available, otherwise print message
    try:
        import pytest
        sys.exit(pytest.main([__file__, '-v']))
    except ImportError:
        print("pytest not installed. Install with: pip install pytest")
        print("Or use: python tests/run_groups_definitions_tests.py")
        sys.exit(1)
