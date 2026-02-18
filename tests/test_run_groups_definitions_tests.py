#!/usr/bin/env python3
"""
Comprehensive test script for HalogenGroups and PFASDefinitions.

This script validates all PFAS groups and definitions against their test metadata,
ensuring that patterns correctly identify positive examples and reject negative examples.

Usage:
    python run_groups_definitions_tests.py [options]
    
Options:
    --groups-only     Test only PFAS groups
    --definitions-only Test only PFAS definitions
    --verbose, -v     Show detailed output for each test
    --group-id ID     Test specific group by ID
    --definition-id ID Test specific definition by ID
    
Examples:
    # Test everything
    python run_groups_definitions_tests.py
    
    # Test only groups with verbose output
    python run_groups_definitions_tests.py --groups-only -v
    
    # Test specific group
    python run_groups_definitions_tests.py --group-id 1
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import pytest

# Add parent directory to path to import HalogenGroups
sys.path.insert(0, str(Path(__file__).parent.parent))

from HalogenGroups import HalogenGroup, PFASDefinition


def load_groups(groups_file: Path) -> List[HalogenGroup]:
    """Load all PFAS groups from JSON file.
    
    Parameters
    ----------
    groups_file : Path
        Path to Halogen_groups_smarts.json
    
    Returns
    -------
    List[HalogenGroup]
        List of HalogenGroup model instances
    """
    with open(groups_file, 'r') as f:
        groups_data = json.load(f)
    
    groups = []
    for data in groups_data:
        try:
            # Initialize HalogenGroup directly from JSON data - no ComponentsSolver needed
            group = HalogenGroup(**data)
            groups.append(group)
        except Exception as e:
            print(f"Warning: Failed to load group {data.get('id', '?')}: {e}")
    
    return groups


def load_definitions(definitions_file: Path) -> List[PFASDefinition]:
    """Load all PFAS definitions from JSON file.
    
    Parameters
    ----------
    definitions_file : Path
        Path to PFAS_definitions_smarts.json
    
    Returns
    -------
    List[PFASDefinition]
        List of PFASDefinition model instances
    """
    with open(definitions_file, 'r') as f:
        definitions_data = json.load(f)
    
    definitions = []
    for data in definitions_data:
        try:
            # PFASDefinition requires specific arguments: id, name, smarts, fluorineRatio, description
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
            print(f"Warning: Failed to load definition {data.get('id', '?')}: {e}")
    
    return definitions


# ============================================================================
# PYTEST FIXTURES
# ============================================================================

@pytest.fixture(scope='session')
def data_dir():
    """Get path to data directory."""
    script_dir = Path(__file__).parent
    root_dir = script_dir.parent
    return root_dir / 'HalogenGroups' / 'data'


@pytest.fixture(scope='session')
def groups(data_dir):
    """Load all PFAS groups for testing."""
    groups_file = data_dir / 'Halogen_groups_smarts.json'
    return load_groups(groups_file)


@pytest.fixture(scope='session')
def definitions(data_dir):
    """Load all PFAS definitions for testing."""
    definitions_file = data_dir / 'PFAS_definitions_smarts.json'
    return load_definitions(definitions_file)


# ============================================================================
# VALIDATION FUNCTIONS (used by both pytest and standalone)
# ============================================================================

def validate_groups(groups: List[HalogenGroup], verbose: bool = False, specific_id: int = None) -> Tuple[int, int, List[Dict]]:
    """Test all PFAS groups against their test metadata.
    
    Parameters
    ----------
    groups : List[HalogenGroup]
        List of PFAS groups to test
    verbose : bool
        If True, print detailed output for each test
    specific_id : int, optional
        If provided, only test the group with this ID
    
    Returns
    -------
    Tuple[int, int, List[Dict]]
        (passed_count, failed_count, failures_details)
    """
    print("\n" + "="*80)
    print("TESTING PFAS GROUPS")
    print("="*80)
    
    passed = 0
    failed = 0
    all_failures = []
    
    for group in groups:
        # Skip if testing specific ID
        if specific_id is not None and group.id != specific_id:
            continue
        
        if verbose:
            print(f"\nTesting Group {group.id}: {group.name}")
        
        try:
            result = group.test()
            
            if result.get('passed') is None:
                if verbose:
                    print(f"  ⚠️  No test data available")
                continue
            
            if result['passed']:
                passed += 1
                if verbose:
                    print(f"  ✅ PASSED ({result['total_tests']} tests)")
            else:
                failed += 1
                if verbose:
                    print(f"  ❌ FAILED ({len(result['failures'])}/{result['total_tests']} failures)")
                    for failure in result['failures']:
                        print(f"     - {failure['smiles'][:50]}...")
                        print(f"       Expected: {failure['expected']}, Got: {failure['got']}")
                        print(f"       Error: {failure['error']}")
                
                all_failures.append({
                    'type': 'group',
                    'id': group.id,
                    'name': group.name,
                    'result': result
                })
        except Exception as e:
            failed += 1
            if verbose:
                print(f"  ❌ EXCEPTION: {str(e)}")
            all_failures.append({
                'type': 'group',
                'id': group.id,
                'name': group.name,
                'exception': str(e)
            })
    
    print(f"\n{'='*80}")
    print(f"Groups Summary: {passed} passed, {failed} failed")
    print(f"{'='*80}\n")
    
    return passed, failed, all_failures


def validate_definitions(definitions: List[PFASDefinition], verbose: bool = False, specific_id: int = None) -> Tuple[int, int, List[Dict]]:
    """Test all PFAS definitions against their test metadata.
    
    Parameters
    ----------
    definitions : List[PFASDefinition]
        List of PFAS definitions to test
    verbose : bool
        If True, print detailed output for each test
    specific_id : int, optional
        If provided, only test the definition with this ID
    
    Returns
    -------
    Tuple[int, int, List[Dict]]
        (passed_count, failed_count, failures_details)
    """
    print("\n" + "="*80)
    print("TESTING PFAS DEFINITIONS")
    print("="*80)
    
    passed = 0
    failed = 0
    all_failures = []
    
    for definition in definitions:
        # Skip if testing specific ID
        if specific_id is not None and definition.id != specific_id:
            continue
        
        if verbose:
            print(f"\nTesting Definition {definition.id}: {definition.name}")
        
        try:
            result = definition.test()
            
            if result.get('passed') is None:
                if verbose:
                    print(f"  ⚠️  No test data available")
                continue
            
            stats = result.get('stats', {})
            total_correct = stats.get('true_positives', 0) + stats.get('true_negatives', 0)
            total_incorrect = stats.get('false_positives', 0) + stats.get('false_negatives', 0)
            
            if result['passed']:
                passed += 1
                if verbose:
                    print(f"  ✅ PASSED ({result['total_tests']} tests)")
                    print(f"     TP: {stats.get('true_positives', 0)}, "
                          f"TN: {stats.get('true_negatives', 0)}, "
                          f"FP: {stats.get('false_positives', 0)}, "
                          f"FN: {stats.get('false_negatives', 0)}")
            else:
                failed += 1
                if verbose:
                    print(f"  ❌ FAILED ({len(result['failures'])}/{result['total_tests']} failures)")
                    print(f"     TP: {stats.get('true_positives', 0)}, "
                          f"TN: {stats.get('true_negatives', 0)}, "
                          f"FP: {stats.get('false_positives', 0)}, "
                          f"FN: {stats.get('false_negatives', 0)}")
                    
                    for failure in result['failures'][:5]:  # Show first 5 failures
                        print(f"     - [{failure['type']}] {failure['smiles'][:50]}...")
                        print(f"       Expected: {failure['expected']}, Got: {failure['got']}")
                        print(f"       Error: {failure['error']}")
                
                all_failures.append({
                    'type': 'definition',
                    'id': definition.id,
                    'name': definition.name,
                    'result': result
                })
        except Exception as e:
            failed += 1
            if verbose:
                print(f"  ❌ EXCEPTION: {str(e)}")
            all_failures.append({
                'type': 'definition',
                'id': definition.id,
                'name': definition.name,
                'exception': str(e)
            })
    
    print(f"\n{'='*80}")
    print(f"Definitions Summary: {passed} passed, {failed} failed")
    print(f"{'='*80}\n")
    
    return passed, failed, all_failures


# ============================================================================
# PYTEST TEST FUNCTIONS
# ============================================================================

def test_all_pfas_groups(groups):
    """Test all PFAS groups against their test metadata."""
    failures = []
    
    for group in groups:
        result = group.test()
        
        # Skip if no test data
        if result.get('passed') is None:
            continue
        
        if not result['passed']:
            failures.append({
                'group': group,
                'result': result
            })
    
    # Build detailed error message if any failures
    if failures:
        error_lines = [f"\n{len(failures)} group(s) failed:"]
        for item in failures:
            group = item['group']
            result = item['result']
            error_lines.append(f"\n  Group {group.id}: {group.name}")
            error_lines.append(f"    Failed: {len(result['failures'])}/{result['total_tests']} tests")
            for failure in result['failures'][:3]:  # Show first 3 failures per group
                error_lines.append(f"      SMILES: {failure['smiles'][:60]}")
                error_lines.append(f"      Expected: {failure['expected']}, Got: {failure['got']}")
        pytest.fail('\n'.join(error_lines))


def test_all_pfas_definitions(definitions):
    """Test all PFAS definitions against their test metadata."""
    failures = []
    
    for definition in definitions:
        result = definition.test()
        
        # Skip if no test data
        if result.get('passed') is None:
            continue
        
        if not result['passed']:
            failures.append({
                'definition': definition,
                'result': result
            })
    
    # Build detailed error message if any failures
    if failures:
        error_lines = [f"\n{len(failures)} definition(s) failed:"]
        for item in failures:
            definition = item['definition']
            result = item['result']
            stats = result.get('stats', {})
            error_lines.append(f"\n  Definition {definition.id}: {definition.name}")
            error_lines.append(f"    Failed: {len(result['failures'])}/{result['total_tests']} tests")
            error_lines.append(f"    TP: {stats.get('true_positives', 0)}, TN: {stats.get('true_negatives', 0)}, "
                             f"FP: {stats.get('false_positives', 0)}, FN: {stats.get('false_negatives', 0)}")
            for failure in result['failures'][:3]:  # Show first 3 failures per definition
                error_lines.append(f"      [{failure['type']}] {failure['smiles'][:60]}")
                error_lines.append(f"      Expected: {failure['expected']}, Got: {failure['got']}")
        pytest.fail('\n'.join(error_lines))


# ============================================================================
# STANDALONE FUNCTIONS
# ============================================================================

def write_failure_report(failures: List[Dict], output_file: Path):
    """Write detailed failure report to file.
    
    Parameters
    ----------
    failures : List[Dict]
        List of failure dictionaries
    output_file : Path
        Path to output file
    """
    with open(output_file, 'w') as f:
        f.write("PFAS GROUPS AND DEFINITIONS TEST FAILURE REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        for failure in failures:
            f.write(f"\n{'='*80}\n")
            f.write(f"{failure['type'].upper()} {failure['id']}: {failure['name']}\n")
            f.write(f"{'='*80}\n\n")
            
            if 'exception' in failure:
                f.write(f"EXCEPTION: {failure['exception']}\n")
            else:
                result = failure['result']
                f.write(f"Total Tests: {result['total_tests']}\n")
                f.write(f"Failures: {len(result['failures'])}\n")
                f.write(f"Category: {result['category']}\n\n")
                
                if 'stats' in result:
                    stats = result['stats']
                    f.write(f"Statistics:\n")
                    f.write(f"  True Positives:  {stats.get('true_positives', 0)}\n")
                    f.write(f"  True Negatives:  {stats.get('true_negatives', 0)}\n")
                    f.write(f"  False Positives: {stats.get('false_positives', 0)}\n")
                    f.write(f"  False Negatives: {stats.get('false_negatives', 0)}\n\n")
                
                f.write("Failure Details:\n")
                for i, fail in enumerate(result['failures'], 1):
                    f.write(f"\n{i}. SMILES: {fail['smiles']}\n")
                    f.write(f"   Expected: {fail['expected']}\n")
                    f.write(f"   Got: {fail['got']}\n")
                    if 'type' in fail:
                        f.write(f"   Type: {fail['type']}\n")
                    f.write(f"   Error: {fail['error']}\n")
    
    print(f"Detailed failure report written to: {output_file}")


def main():
    """Main test execution."""
    parser = argparse.ArgumentParser(
        description='Test PFAS groups and definitions against test metadata',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--groups-only', action='store_true',
                        help='Test only PFAS groups')
    parser.add_argument('--definitions-only', action='store_true',
                        help='Test only PFAS definitions')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show detailed output for each test')
    parser.add_argument('--group-id', type=int,
                        help='Test specific group by ID')
    parser.add_argument('--definition-id', type=int,
                        help='Test specific definition by ID')
    parser.add_argument('--output', '-o', type=Path,
                        default=Path('test_failures_report.txt'),
                        help='Output file for failure report (default: test_failures_report.txt)')
    
    args = parser.parse_args()
    
    # Determine paths
    script_dir = Path(__file__).parent
    root_dir = script_dir.parent
    data_dir = root_dir / 'HalogenGroups' / 'data'
    
    groups_file = data_dir / 'Halogen_groups_smarts.json'
    definitions_file = data_dir / 'PFAS_definitions_smarts.json'
    
    # Check files exist
    if not groups_file.exists():
        print(f"Error: Groups file not found: {groups_file}")
        return 1
    if not definitions_file.exists():
        print(f"Error: Definitions file not found: {definitions_file}")
        return 1
    
    all_failures = []
    total_passed = 0
    total_failed = 0
    
    # Test groups
    if not args.definitions_only:
        groups = load_groups(groups_file)
        print(f"Loaded {len(groups)} PFAS groups")
        
        passed, failed, failures = validate_groups(
            groups, 
            verbose=args.verbose,
            specific_id=args.group_id
        )
        total_passed += passed
        total_failed += failed
        all_failures.extend(failures)
    
    # Test definitions
    if not args.groups_only:
        definitions = load_definitions(definitions_file)
        print(f"Loaded {len(definitions)} PFAS definitions")
        
        passed, failed, failures = validate_definitions(
            definitions,
            verbose=args.verbose,
            specific_id=args.definition_id
        )
        total_passed += passed
        total_failed += failed
        all_failures.extend(failures)
    
    # Print overall summary
    print("\n" + "="*80)
    print("OVERALL SUMMARY")
    print("="*80)
    print(f"Total Passed: {total_passed}")
    print(f"Total Failed: {total_failed}")
    print(f"Success Rate: {100*total_passed/(total_passed+total_failed) if total_passed+total_failed > 0 else 0:.1f}%")
    print("="*80 + "\n")
    
    # Write failure report if there are failures
    if all_failures:
        write_failure_report(all_failures, args.output)
        return 1
    else:
        print("✅ All tests passed!")
        return 0


if __name__ == '__main__':
    sys.exit(main())
