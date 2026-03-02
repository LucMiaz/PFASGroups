#!/usr/bin/env python3
"""
Test suite for highly branched PFAS compounds
Tests functional groups 29-59 (excluding cyclic 54-57 and groups without SMARTS 49-50)
attached at different distances (0, 1, 2 bonds) from perfluorinated components of varying sizes
Compares HalogenGroups and PFAS-Atlas detection
"""

import sys
import os
import json
from pathlib import Path
from datetime import datetime

# Add parent directory to path to import HalogenGroups
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

try:
    from PFASGroups.parser import parse_smiles
    HalogenGroupS_AVAILABLE = True
except ImportError:
    print("❌ Error: Could not import HalogenGroups module")
    print("   Make sure you're running from the benchmark directory")
    print("   and HalogenGroups is installed in your environment")
    HalogenGroupS_AVAILABLE = False
    sys.exit(1)

# Try to import PFAS-Atlas
repo_root = Path(__file__).resolve().parents[2]
atlas_dir = repo_root.parent / 'PFAS-atlas'
try:
    sys.path.append(str(atlas_dir))
    from classification_helper.classify_pfas import classify_pfas_molecule
    ATLAS_AVAILABLE = True
    print("✅ PFAS-Atlas available")
except ImportError:
    try:
        sys.path.append(str(atlas_dir / 'classification_helper'))
        from classify_pfas import classify_pfas_molecule
        ATLAS_AVAILABLE = True
        print("✅ PFAS-Atlas available (fallback import)")
    except ImportError:
        print("⚠️  PFAS-Atlas not available - will only test HalogenGroups")
        ATLAS_AVAILABLE = False

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
    39: { "name": 'sulfuric acid', "smarts": 'OS(=O)(=O)O', "distance": [0, 1, 2] },
    40: { "name": 'phosphonic acid', "smarts": 'P(=O)(O)O', "distance": [0, 1, 2] },
    41: { "name": 'phosphinic acid', "smarts": 'P(=O)(O)', "distance": [0, 1, 2] },
    42: { "name": 'chloride', "smarts": 'Cl', "distance": [0, 1, 2] },
    43: { "name": 'bromide', "smarts": 'Br', "distance": [0, 1, 2] },
    44: { "name": 'iodide', "smarts": 'I', "distance": [0, 1, 2] },
    45: { "name": 'sulfonamide', "smarts": 'CS(=O)(=O)N', "distance": [0, 1, 2] },
    46: { "name": 'heterocyclic azole', "smarts": 'C1NCNC1', "distance": [0, 1, 2] },
    47: { "name": 'heterocyclic azine', "smarts": 'c1cccnc1', "distance": [0, 1, 2] },
    48: { "name": 'amine', "smarts": 'N', "distance": [0, 1, 2] },
    52: { "name": 'alkene', "smarts": 'C=C', "distance": [0, 1, 2] },
    53: { "name": 'alkyne', "smarts": 'C#C', "distance": [0, 1, 2] },
    54: { "name": 'side-chain aromatics', "smarts": 'c1ccccc1', "distance": [0, 1, 2] },
    59: { "name": 'peroxydes', "smarts": 'OO', "distance": [0, 1, 2] },
    60: { "name": 'benzoyl peroxydes', "smarts": 'c1ccccc1C(=O)OOC(=O)', "distance": [0, 1, 2] }
}

QUICK_RUN = os.getenv("PFAS_BENCH_QUICK") == "1"
BRANCHED_GROUP_LIMIT = int(os.getenv("PFAS_BENCH_BRANCHED_GROUPS", "0") or "0")
BRANCHED_COMPONENT_MAX = int(os.getenv("PFAS_BENCH_BRANCHED_COMPONENTS", "5") or "5")
BRANCHED_DISTANCE_MAX = int(os.getenv("PFAS_BENCH_BRANCHED_DISTANCE_MAX", "2") or "2")


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


def analyze_molecule_HalogenGroups(smiles):
    """Analyze a molecule using HalogenGroups"""
    try:
        result = parse_smiles(smiles)
        return result
    except Exception as e:
        print(f"     Error analyzing with HalogenGroups: {str(e)}")
        return None


def analyze_molecule_atlas(smiles):
    """Analyze a molecule using PFAS-Atlas"""
    if not ATLAS_AVAILABLE:
        return None
    try:
        predictions = classify_pfas_molecule(smiles)
        return predictions
    except Exception as e:
        print(f"     Error analyzing with PFAS-Atlas: {str(e)}")
        return None


def run_benchmark():
    """Run the highly branched compounds benchmark"""
    print("🧪 PFAS Highly Branched Compounds Benchmark")
    print("=" * 80)
    group_items = list(FUNCTIONAL_GROUPS.items())
    if BRANCHED_GROUP_LIMIT > 0:
        group_items = group_items[:BRANCHED_GROUP_LIMIT]
    component_max = max(1, BRANCHED_COMPONENT_MAX)
    distance_max = max(0, BRANCHED_DISTANCE_MAX)
    
    print(f"Testing {len(group_items)} functional groups")
    print(f"Component sizes: 1-{component_max} CF2 units")
    print(f"Distances: 0-{distance_max} bonds from component start")
    if ATLAS_AVAILABLE:
        print("Comparing HalogenGroups and PFAS-Atlas detection")
    else:
        print("Testing HalogenGroups only (PFAS-Atlas not available)")
    print("=" * 80)
    print()
    
    results = {
        'total': 0,
        'HalogenGroups_success': 0,
        'HalogenGroups_failed': 0,
        'atlas_success': 0,
        'atlas_failed': 0,
        'details': []
    }
    
    # Test each component size
    for component_size in range(1, component_max + 1):
        print(f"\n📊 Testing component size {component_size}")
        print("-" * 80)
        
        for group_id, group_info in group_items:
            distances = [d for d in group_info['distance'] if d <= distance_max]
            for distance in distances:
                results['total'] += 1
                smiles = attach_functional_group(component_size, group_info['smarts'], distance)
                
                print(f"\nTest {results['total']}:")
                print(f"  Group {group_id}: {group_info['name']}")
                print(f"  Component Size: {component_size}, Distance: {distance}")
                print(f"  SMILES: {smiles}")
                
                # Test with HalogenGroups
                pfas_analysis = analyze_molecule_HalogenGroups(smiles)
                pfas_passed = False
                
                if pfas_analysis and len(pfas_analysis) > 0:
                    # Check if the expected group was detected
                    found_group = any(m['id'] == group_id for a in pfas_analysis for m in a['matches'])
                    found_perfluoro = any(m['id'] in [49, 50] for a in pfas_analysis for m in a['matches'])
                    
                    if found_group:
                        results['HalogenGroups_success'] += 1
                        pfas_passed = True
                        print(f"  ✅ HalogenGroups PASS - Group {group_id} ({group_info['name']}) detected")
                        
                        # Log component information
                        for m in pfas_analysis:
                            for match in m['matches']:
                                if match['id'] == group_id and 'components' in match and match['components']:
                                    comp = match['components'][0]
                                    size = comp.get('size', len(comp.get('component', [])))
                                    print(f"       Component size: {size}")
                                    break
                    else:
                        results['HalogenGroups_failed'] += 1
                        print(f"  ❌ HalogenGroups FAIL - Group {group_id} ({group_info['name']}) NOT detected")
                        detected = ', '.join([f"{m['id']}:{m.get('group_name', m.get('definition_name', 'unknown'))}" for a in pfas_analysis for m in a['matches'][:5]])
                        print(f"     Detected groups: {detected}")
                    
                    if found_perfluoro:
                        for m in pfas_analysis:
                            for match in m['matches']:
                                if match['id'] in [49, 50]:
                                    chain_size = match.get('nCFchain', [None])[0]
                                    print(f"     Perfluoroalkyl component: ✓ (size: {chain_size})")
                                    break
                else:
                    results['HalogenGroups_failed'] += 1
                    print(f"  ❌ HalogenGroups FAIL - Analysis error or invalid SMILES")
                
                # Test with PFAS-Atlas
                atlas_passed = False
                atlas_predictions = None
                if ATLAS_AVAILABLE:
                    atlas_predictions = analyze_molecule_atlas(smiles)
                    if atlas_predictions and len(atlas_predictions) >= 1:
                        results['atlas_success'] += 1
                        atlas_passed = True
                        print(f"  ✅ PFAS-Atlas PASS - Classified as PFAS: {atlas_predictions[0]}")
                        if len(atlas_predictions) >= 2:
                            print(f"     Second class: {atlas_predictions[1]}")
                    else:
                        results['atlas_failed'] += 1
                        print(f"  ❌ PFAS-Atlas FAIL - Not classified as PFAS")
                
                # Store detailed results
                results['details'].append({
                    'component_size': component_size,
                    'group_id': group_id,
                    'group_name': group_info['name'],
                    'distance': distance,
                    'smiles': smiles,
                    'HalogenGroups_passed': pfas_passed,
                    'atlas_passed': atlas_passed if ATLAS_AVAILABLE else None,
                    'atlas_predictions': atlas_predictions if ATLAS_AVAILABLE else None
                })
    
    # Print summary
    print("\n" + "=" * 80)
    print("📈 TEST SUMMARY")
    print("=" * 80)
    print(f"Total tests: {results['total']}")
    print(f"\nHalogenGroups Results:")
    print(f"  Passed: {results['HalogenGroups_success']} ({results['HalogenGroups_success']/results['total']*100:.1f}%)")
    print(f"  Failed: {results['HalogenGroups_failed']} ({results['HalogenGroups_failed']/results['total']*100:.1f}%)")
    
    if ATLAS_AVAILABLE:
        print(f"\nPFAS-Atlas Results:")
        print(f"  Passed: {results['atlas_success']} ({results['atlas_success']/results['total']*100:.1f}%)")
        print(f"  Failed: {results['atlas_failed']} ({results['atlas_failed']/results['total']*100:.1f}%)")
    
    HalogenGroups_success_pct = results['HalogenGroups_success'] / results['total'] * 100
    
    # Failures by group for HalogenGroups
    HalogenGroups_failures_by_group = {}
    for detail in results['details']:
        if not detail['HalogenGroups_passed']:
            gid = detail['group_id']
            if gid not in HalogenGroups_failures_by_group:
                HalogenGroups_failures_by_group[gid] = {
                    'name': detail['group_name'],
                    'count': 0,
                    'distances': set(),
                    'component_sizes': set()
                }
            HalogenGroups_failures_by_group[gid]['count'] += 1
            HalogenGroups_failures_by_group[gid]['distances'].add(detail['distance'])
            HalogenGroups_failures_by_group[gid]['component_sizes'].add(detail['component_size'])
    
    if HalogenGroups_failures_by_group:
        print("\n❌ HalogenGroups Failures by Functional Group:")
        print("-" * 80)
        for gid, info in sorted(HalogenGroups_failures_by_group.items()):
            print(f"  Group {gid} ({info['name']}): {info['count']} failures")
            print(f"    At distances: {', '.join(map(str, sorted(info['distances'])))}")
            print(f"    With component sizes: {', '.join(map(str, sorted(info['component_sizes'])))}")
    
    # Failures by group for PFAS-Atlas
    atlas_failures_by_group = {}
    if ATLAS_AVAILABLE:
        for detail in results['details']:
            if detail['atlas_passed'] is not None and not detail['atlas_passed']:
                gid = detail['group_id']
                if gid not in atlas_failures_by_group:
                    atlas_failures_by_group[gid] = {
                        'name': detail['group_name'],
                        'count': 0,
                        'distances': set(),
                        'component_sizes': set()
                    }
                atlas_failures_by_group[gid]['count'] += 1
                atlas_failures_by_group[gid]['distances'].add(detail['distance'])
                atlas_failures_by_group[gid]['component_sizes'].add(detail['component_size'])
        
        if atlas_failures_by_group:
            print("\n❌ PFAS-Atlas Failures by Functional Group:")
            print("-" * 80)
            for gid, info in sorted(atlas_failures_by_group.items()):
                print(f"  Group {gid} ({info['name']}): {info['count']} failures")
                print(f"    At distances: {', '.join(map(str, sorted(info['distances'])))}")
                print(f"    With component sizes: {', '.join(map(str, sorted(info['component_sizes'])))}")
    
    # Success rate by distance for HalogenGroups
    print("\n📊 HalogenGroups Success Rate by Distance:")
    print("-" * 80)
    for dist in range(3):
        dist_tests = [d for d in results['details'] if d['distance'] == dist]
        dist_success = sum(1 for d in dist_tests if d['HalogenGroups_passed'])
        dist_pct = (dist_success / len(dist_tests) * 100) if dist_tests else 0
        print(f"  Distance {dist}: {dist_success}/{len(dist_tests)} ({dist_pct:.1f}%)")
    
    # Success rate by component size for HalogenGroups
    print("\n📊 HalogenGroups Success Rate by Component Size:")
    print("-" * 80)
    for size in range(1, 6):
        size_tests = [d for d in results['details'] if d['component_size'] == size]
        size_success = sum(1 for d in size_tests if d['HalogenGroups_passed'])
        size_pct = (size_success / len(size_tests) * 100) if size_tests else 0
        print(f"  Size {size}: {size_success}/{len(size_tests)} ({size_pct:.1f}%)")
    
    # PFAS-Atlas analysis
    if ATLAS_AVAILABLE:
        print("\n📊 PFAS-Atlas Success Rate by Distance:")
        print("-" * 80)
        for dist in range(3):
            dist_tests = [d for d in results['details'] if d['distance'] == dist and d['atlas_passed'] is not None]
            dist_success = sum(1 for d in dist_tests if d['atlas_passed'])
            dist_pct = (dist_success / len(dist_tests) * 100) if dist_tests else 0
            print(f"  Distance {dist}: {dist_success}/{len(dist_tests)} ({dist_pct:.1f}%)")
        
        print("\n📊 PFAS-Atlas Success Rate by Component Size:")
        print("-" * 80)
        for size in range(1, 6):
            size_tests = [d for d in results['details'] if d['component_size'] == size and d['atlas_passed'] is not None]
            size_success = sum(1 for d in size_tests if d['atlas_passed'])
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
            'HalogenGroups_success_rate': HalogenGroups_success_pct,
            'atlas_available': ATLAS_AVAILABLE,
            'atlas_success_rate': (results['atlas_success'] / results['total'] * 100) if ATLAS_AVAILABLE else None
        },
        'summary': {
            'total': results['total'],
            'HalogenGroups_passed': results['HalogenGroups_success'],
            'HalogenGroups_failed': results['HalogenGroups_failed'],
            'atlas_passed': results['atlas_success'] if ATLAS_AVAILABLE else None,
            'atlas_failed': results['atlas_failed'] if ATLAS_AVAILABLE else None
        },
        'details': results['details'],
        'HalogenGroups_failures_by_group': {
            str(gid): {
                'name': info['name'],
                'count': info['count'],
                'distances': sorted(list(info['distances'])),
                'component_sizes': sorted(list(info['component_sizes']))
            }
            for gid, info in HalogenGroups_failures_by_group.items()
        },
        'atlas_failures_by_group': {
            str(gid): {
                'name': info['name'],
                'count': info['count'],
                'distances': sorted(list(info['distances'])),
                'component_sizes': sorted(list(info['component_sizes']))
            }
            for gid, info in atlas_failures_by_group.items()
        } if ATLAS_AVAILABLE and atlas_failures_by_group else {}
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"📁 Results saved to: {output_file}")
    print()
    
    return results


if __name__ == "__main__":
    run_benchmark()
