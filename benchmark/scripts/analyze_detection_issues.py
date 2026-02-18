#!/usr/bin/env python3
"""Analyze detection issues in benchmark results"""

import json
import sys

# Load the most recent benchmark data
with open('/home/luc/git/HalogenGroups/benchmark/data/pfas_enhanced_benchmark_20260201_181935.json', 'r') as f:
    results = json.load(f)

# Analyze detection rates for groups 67, 104, 106 and others with issues
target_groups = [67, 104, 106]

print("🔍 ANALYZING DETECTION ISSUES\n")
print("=" * 80)

for group_id in target_groups:
    # Find all molecules for this group
    group_molecules = [r for r in results if r.get('molecule_data', {}).get('group_id') == group_id]
    
    if not group_molecules:
        print(f"\n❌ Group {group_id}: No molecules found in benchmark data")
        continue
    
    detected_count = 0
    not_detected = []
    
    for mol_result in group_molecules:
        mol_data = mol_result.get('molecule_data', {})
        pfas_result = mol_result.get('HalogenGroups_result', {})
        detected_groups = pfas_result.get('detected_groups', [])
        
        smiles = mol_data.get('smiles', 'N/A')
        
        if group_id in detected_groups:
            detected_count += 1
        else:
            not_detected.append({
                'smiles': smiles,
                'detected_groups': detected_groups,
                'generation_type': mol_data.get('generation_type', 'unknown')
            })
    
    total = len(group_molecules)
    detection_rate = (detected_count / total * 100) if total > 0 else 0
    
    print(f"\n📊 Group {group_id}:")
    print(f"   Total molecules: {total}")
    print(f"   Detected: {detected_count}/{total} ({detection_rate:.1f}%)")
    print(f"   Not detected: {len(not_detected)}")
    
    if not_detected:
        print(f"\n   ❌ Molecules NOT detected as group {group_id}:")
        for i, mol in enumerate(not_detected[:10], 1):  # Show first 10
            print(f"   {i}. {mol['smiles']}")
            print(f"      Detected as: {mol['detected_groups']}")
        
        if len(not_detected) > 10:
            print(f"   ... and {len(not_detected) - 10} more")

# Also check for all groups with <100% detection
print(f"\n\n🔍 ALL GROUPS WITH <100% DETECTION:\n")
print("=" * 80)

# Group all molecules by target group
groups_summary = {}
for r in results:
    mol_data = r.get('molecule_data', {})
    group_id = mol_data.get('group_id')
    
    if group_id is None:
        continue
    
    if group_id not in groups_summary:
        groups_summary[group_id] = {
            'total': 0,
            'detected': 0,
            'not_detected_molecules': []
        }
    
    groups_summary[group_id]['total'] += 1
    
    pfas_result = r.get('HalogenGroups_result', {})
    detected_groups = pfas_result.get('detected_groups', [])
    
    if group_id in detected_groups:
        groups_summary[group_id]['detected'] += 1
    else:
        groups_summary[group_id]['not_detected_molecules'].append({
            'smiles': mol_data.get('smiles', 'N/A'),
            'detected': detected_groups
        })

# Sort by detection rate
sorted_groups = sorted(groups_summary.items(), key=lambda x: x[1]['detected']/max(x[1]['total'], 1))

for group_id, stats in sorted_groups:
    if stats['total'] > 0:
        detection_rate = (stats['detected'] / stats['total']) * 100
        if detection_rate < 100:
            print(f"\nGroup {group_id}: {stats['detected']}/{stats['total']} ({detection_rate:.1f}% detection)")
            if stats['not_detected_molecules']:
                print(f"   Undetected molecules (showing first 3):")
                for i, mol in enumerate(stats['not_detected_molecules'][:3], 1):
                    print(f"   {i}. {mol['smiles']}")
                    print(f"      Detected as: {mol['detected']}")
