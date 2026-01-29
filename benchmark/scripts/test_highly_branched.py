#!/usr/bin/env python3
"""
Test suite for highly branched PFAS compounds
Tests functional groups 29-59 (excluding cyclic 54-57 and groups without SMARTS 49-50)
attached at different distances (0, 1, 2 bonds) from perfluorinated components of varying sizes
"""

import sys
import os
import json
from datetime import datetime

# Add parent directory to path to import pfasgroups
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

try:
    from PFASgroups.parser import parse_smiles
except ImportError:
    print("❌ Error: Could not import pfasgroups module")
    print("   Make sure you're running from the benchmark directory")
    print("   and pfasgroups is installed in your environment")
    sys.exit(1)

# Functional groups to test (29-59, excluding cyclic and no-SMARTS)
FUNCTIONAL_GROUPS = {
    29: { "name": 'alcohol', "smarts": 'O', "distance": [0, 1, 2] },
    30: { "name": 'ketone', "smarts": 'FC(F)(F)C(=O)', "distance": [0, 1, 2] },
    31: { "name": 'ether', "smarts": 'CO', "distance": [0, 1, 2] },
    32: { "name": 'ester', "smarts": 'FC(F)(F)OC(=O)', "distance": [0, 1, 2] },
    33: { "name": 'carboxylic acid', "smarts": 'OC(=O)', "distance": [0, 1, 2] },
    34: { "name": 'amide', "smarts": 'NC(=O)', "distance": [0, 1, 2] },
    35: { "name": 'acyl halide', "smarts": 'ClC(=O)', "distance": [0, 1, 2] },
    36: { "name": 'sulfonic acid', "smarts": 'OS(=O)(=O)', "distance": [0, 1, 2] },
    37: { "name": 'sulfenic acid', "smarts": 'OS', "distance": [0, 1, 2] },
    38: { "name": 'sulfinic acid', "smarts": 'OS(=O)', "distance": [0, 1, 2] },
    39: { "name": 'phosphonic acid', "smarts": 'OP(=O)(O)', "distance": [0, 1, 2] },
    40: { "name": 'phosphinic acid', "smarts": 'P(=O)(O)', "distance": [0, 1, 2] },
    41: { "name": 'chloride', "smarts": 'Cl', "distance": [0, 1, 2] },
    42: { "name": 'bromide', "smarts": 'Br', "distance": [0, 1, 2] },
    43: { "name": 'iodide', "smarts": 'I', "distance": [0, 1, 2] },
    44: { "name": 'sulfonamide', "smarts": 'CS(=O)(=O)N', "distance": [0, 1, 2] },
    45: { "name": 'heterocyclic azole', "smarts": 'C1NCNC1', "distance": [0, 1, 2] },
    46: { "name": 'heterocyclic azine', "smarts": 'c1cccnc1', "distance": [0, 1, 2] },
    48: { "name": 'amine', "smarts": 'N', "distance": [0, 1, 2] },
    51: { "name": 'alkene', "smarts": 'C=C', "distance": [0, 1, 2] },
    52: { "name": 'alkyne', "smarts": 'C#C', "distance": [0, 1, 2] },
    53: { "name": 'side-chain aromatics', "smarts": 'c1ccccc1', "distance": [0, 1, 2] },
    58: { "name": 'peroxydes', "smarts": 'OO', "distance": [0, 1, 2] },
    59: { "name": 'benzoyl peroxydes', "smarts": 'c1ccccc1C(=O)OOC(=O)', "distance": [0, 1, 2] }
}


def attach_functional_group(component_size, func_group, distance):
    """
    Attach functional group to perfluorinated component at specified distance
    
    Args:
        component_size: Number of CF2 units (1-5)
        func_group: Functional group SMILES fragment
        distance: Distance from component start (0, 1, or 2 CF2 units)
    
    Returns:
        Complete SMILES string
    """
    # Build: funcGroup + (CF2)distance + (CF2)component_size-distance + CF3
    smiles = func_group
    
    # Add CF2 units for distance
    for _ in range(distance):
        smiles += 'C(F)(F)'
    
    # Add remaining CF2 units
    for _ in range(distance, component_size):
        smiles += 'C(F)(F)'
    
    # Add terminal CF3
    smiles += 'C(F)(F)F'
    
    return smiles


def analyze_molecule(smiles):
    """Analyze a molecule using pfasgroups"""
    try:
        result = parse_smiles(smiles)
        return result
    except Exception as e:
        print(f"     Error analyzing: {str(e)}")
        return None


def run_benchmark():
    """Run the highly branched compounds benchmark"""
    print("🧪 PFAS Highly Branched Compounds Benchmark")
    print("=" * 80)
    print(f"Testing {len(FUNCTIONAL_GROUPS)} functional groups")
    print("Component sizes: 1-5 CF2 units")
    print("Distances: 0, 1, 2 bonds from component start")
    print("=" * 80)
    print()
    
    results = {
        'total': 0,
        'success': 0,
        'failed': 0,
        'details': []
    }
    
    # Test each component size
    for component_size in range(1, 6):
        print(f"\n📊 Testing component size {component_size}")
        print("-" * 80)
        
        for group_id, group_info in FUNCTIONAL_GROUPS.items():
            for distance in group_info['distance']:
                results['total'] += 1
                smiles = attach_functional_group(component_size, group_info['smarts'], distance)
                
                print(f"\nTest {results['total']}:")
                print(f"  Group {group_id}: {group_info['name']}")
                print(f"  Component Size: {component_size}, Distance: {distance}")
                print(f"  SMILES: {smiles}")
                
                analysis = analyze_molecule(smiles)
                
                if analysis and len(analysis)>0:
                    # Check if the expected group was detected
                    found_group = any(m['id'] == group_id for a in analysis for m in a['matches'])
                    found_perfluoro = any(m['id'] in [49, 50] for a in analysis for m in a['matches'])
                    
                    if found_group:
                        results['success'] += 1
                        print(f"  ✅ PASS - Group {group_id} ({group_info['name']}) detected")
                        
                        # Log component information if available
                        for m in analysis:
                            for match in m['matches']:
                                if match['id'] == group_id and 'components' in match and match['components']:
                                    comp = match['components'][0]
                                    size = comp.get('size', len(comp.get('component', [])))
                                    diameter = comp.get('diameter', 'N/A')
                                    radius = comp.get('radius', 'N/A')
                                    ecc = comp.get('eccentricity', 'N/A')
                                    print(f"     Component size: {size}")
                                    print(f"     Diameter/Radius: {diameter}/{radius}")
                                    if isinstance(ecc, float):
                                        print(f"     Eccentricity: {ecc:.3f}")
                                    else:
                                        print(f"     Eccentricity: {ecc}")
                    else:
                        results['failed'] += 1
                        print(f"  ❌ FAIL - Group {group_id} ({group_info['name']}) NOT detected")
                        detected = ', '.join([f"{m['id']}:{m['name']}" for a in analysis for m in a['matches'][:5]])
                        print(f"     Detected groups: {detected}")
                    
                    if found_perfluoro:
                        for m in analysis:
                            for match in m['matches']:
                                if match['id'] in [49, 50]:
                                    chain_size = match.get('nCFchain', [None])[0]
                                    print(f"     Perfluoroalkyl component: ✓ (size: {chain_size})")
                                    break
                else:
                    results['failed'] += 1
                    print(f"  ❌ FAIL - Analysis error or invalid SMILES")
                
                # Store detailed results
                results['details'].append({
                    'component_size': component_size,
                    'group_id': group_id,
                    'group_name': group_info['name'],
                    'distance': distance,
                    'smiles': smiles,
                    'passed': analysis and len(analysis) and any(m['id'] == group_id for a in analysis for m in a['matches'])
                })
    
    # Print summary
    print("\n" + "=" * 80)
    print("📈 TEST SUMMARY")
    print("=" * 80)
    print(f"Total tests: {results['total']}")
    success_pct = (results['success'] / results['total'] * 100) if results['total'] > 0 else 0
    failed_pct = (results['failed'] / results['total'] * 100) if results['total'] > 0 else 0
    print(f"✅ Passed: {results['success']} ({success_pct:.1f}%)")
    print(f"❌ Failed: {results['failed']} ({failed_pct:.1f}%)")
    
    # Failures by group
    failures_by_group = {}
    for detail in results['details']:
        if not detail['passed']:
            gid = detail['group_id']
            if gid not in failures_by_group:
                failures_by_group[gid] = {
                    'name': detail['group_name'],
                    'count': 0,
                    'distances': set(),
                    'component_sizes': set()
                }
            failures_by_group[gid]['count'] += 1
            failures_by_group[gid]['distances'].add(detail['distance'])
            failures_by_group[gid]['component_sizes'].add(detail['component_size'])
    
    if failures_by_group:
        print("\n❌ Failures by Functional Group:")
        print("-" * 80)
        for gid, info in sorted(failures_by_group.items()):
            print(f"  Group {gid} ({info['name']}): {info['count']} failures")
            print(f"    At distances: {', '.join(map(str, sorted(info['distances'])))}")
            print(f"    With component sizes: {', '.join(map(str, sorted(info['component_sizes'])))}")
    
    # Success rate by distance
    print("\n📊 Success Rate by Distance:")
    print("-" * 80)
    for dist in range(3):
        dist_tests = [d for d in results['details'] if d['distance'] == dist]
        dist_success = sum(1 for d in dist_tests if d['passed'])
        dist_pct = (dist_success / len(dist_tests) * 100) if dist_tests else 0
        print(f"  Distance {dist}: {dist_success}/{len(dist_tests)} ({dist_pct:.1f}%)")
    
    # Success rate by component size
    print("\n📊 Success Rate by Component Size:")
    print("-" * 80)
    for size in range(1, 6):
        size_tests = [d for d in results['details'] if d['component_size'] == size]
        size_success = sum(1 for d in size_tests if d['passed'])
        size_pct = (size_success / len(size_tests) * 100) if size_tests else 0
        print(f"  Size {size}: {size_success}/{len(size_tests)} ({size_pct:.1f}%)")
    
    print("\n" + "=" * 80)
    print("✨ Test suite complete!")
    print()
    
    # Save results to JSON
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"data/pfas_highly_branched_benchmark_{timestamp}.json"
    
    output_data = {
        'metadata': {
            'timestamp': timestamp,
            'test_type': 'highly_branched',
            'total_tests': results['total'],
            'success_rate': success_pct,
            'failure_rate': failed_pct
        },
        'summary': {
            'total': results['total'],
            'passed': results['success'],
            'failed': results['failed']
        },
        'details': results['details'],
        'failures_by_group': {
            str(gid): {
                'name': info['name'],
                'count': info['count'],
                'distances': sorted(list(info['distances'])),
                'component_sizes': sorted(list(info['component_sizes']))
            }
            for gid, info in failures_by_group.items()
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"📁 Results saved to: {output_file}")
    print()
    
    return results


if __name__ == "__main__":
    run_benchmark()
