#!/usr/bin/env python3
"""Extract detection statistics from benchmark datasets."""

import json
from pathlib import Path
from collections import Counter

# Load all benchmark files
datasets = {
    'Enhanced': 'data/pfas_enhanced_benchmark_20260129_003129.json',
    'OECD': 'data/pfas_oecd_benchmark_20260129_003748.json',
    'Complex Branched': 'data/pfas_complex_branched_benchmark_20260129_004804.json',
    'Highly Branched': 'data/pfas_highly_branched_benchmark_20260129_004813.json',
    'Telomer Validation': 'data/telomer_validation_results.json'
}

print('PFAS Group Detection Statistics Across Benchmark Datasets')
print('='*80)
print()

for name, path in datasets.items():
    if not Path(path).exists():
        print(f'{name}: File not found')
        continue
    
    with open(path) as f:
        data = json.load(f)
    
    print(f'{name}:')
    
    if name == 'Telomer Validation':
        total = data.get('total_molecules', 0)
        detected = data.get('telomer_detected', 0)
        group_counts = data.get('group_counts', [])
        print(f'  Total molecules: {total}')
        print(f'  Telomers detected: {detected} ({100*detected/total:.1f}%)')
        print(f'  Unique telomer groups: {len(group_counts)}')
        print(f'  Top 5 telomer groups:')
        for group in group_counts[:5]:
            print(f'    Group {group["id"]}: {group["count"]} detections ({group["name"]})')
    elif name == 'Highly Branched':
        # Different structure - has 'passed' field, not 'detected'
        details = data.get('details', [])
        total = len(details)
        detected = sum(1 for d in details if d.get('passed', False))
        group_counts = Counter()
        for d in details:
            if d.get('passed'):
                # Each test has a target group_id
                gid = d.get('group_id')
                if gid:
                    group_counts[gid] += 1
        print(f'  Total molecules: {total}')
        print(f'  PFAS groups detected: {detected} ({100*detected/total:.1f}%)')
        print(f'  Unique PFAS groups: {len(group_counts)}')
        print(f'  Top 5 groups:')
        for gid, count in group_counts.most_common(5):
            print(f'    Group {gid}: {count} detections')
    elif name == 'Complex Branched':
        # List of test objects, each with molecules
        if not isinstance(data, list):
            print(f'  Unexpected data format')
            continue
        
        total_molecules = 0
        detected = 0
        group_counts = Counter()
        
        for test in data:
            molecules = test.get('molecules', [])
            total_molecules += len(molecules)
            
            for mol in molecules:
                if mol.get('pfasgroups_detected', False):
                    detected += 1
                    groups = mol.get('pfasgroups_groups', [])
                    for gid in groups:
                        group_counts[gid] += 1
        
        print(f'  Total molecules: {total_molecules}')
        print(f'  Molecules with PFAS groups: {detected} ({100*detected/total_molecules:.1f}%)')
        print(f'  Unique PFAS groups detected: {len(group_counts)}')
        print(f'  Top 5 groups:')
        for gid, count in group_counts.most_common(5):
            print(f'    Group {gid}: {count} detections')
    else:
        # Handle both dict and list formats
        if isinstance(data, list):
            molecules = data
        else:
            molecules = data.get('molecules', [])
        
        if not molecules:
            print(f'  No molecule data found')
            continue
        
        total = len(molecules)
        
        # Count molecules with PFASgroups detections
        detected = 0
        group_counts = Counter()
        
        for mol in molecules:
            # Handle different formats
            pfas_results = mol.get('pfasgroups_result', mol.get('pfasgroups_results', {}))
            if isinstance(pfas_results, dict):
                # Check for detected_groups (enhanced format) or matches
                detected_groups = pfas_results.get('detected_groups', [])
                matches = pfas_results.get('matches', [])
                
                if detected_groups or matches:
                    detected += 1
                    
                    # Count from detected_groups if available
                    for gid in detected_groups:
                        group_counts[gid] += 1
                    
                    # Or count from matches
                    for match in matches:
                        if match.get('type') in ['PFASgroup', 'group']:
                            group_counts[match['id']] += 1
        
        print(f'  Total molecules: {total}')
        print(f'  Molecules with PFAS groups: {detected} ({100*detected/total:.1f}%)')
        print(f'  Unique PFAS groups detected: {len(group_counts)}')
        print(f'  Top 5 groups:')
        for gid, count in group_counts.most_common(5):
            print(f'    Group {gid}: {count} detections')
    
    print()

print('='*80)
print('Summary table for LaTeX:')
print()
print('Dataset & Total & Detected & Detection Rate & Unique Groups \\\\')
print('\\hline')

for name, path in datasets.items():
    if not Path(path).exists():
        continue
    
    with open(path) as f:
        data = json.load(f)
    
    if name == 'Telomer Validation':
        total = data.get('total_molecules', 0)
        detected = data.get('telomer_detected', 0)
        unique = len(data.get('group_counts', []))
        rate = 100*detected/total if total > 0 else 0
        print(f'{name} & {total} & {detected} & {rate:.1f}\\% & {unique} \\\\')
    elif name == 'Highly Branched':
        details = data.get('details', [])
        total = len(details)
        detected = sum(1 for d in details if d.get('passed', False))
        group_counts = Counter()
        for d in details:
            if d.get('passed'):
                gid = d.get('group_id')
                if gid:
                    group_counts[gid] += 1
        unique = len(group_counts)
        rate = 100*detected/total if total > 0 else 0
        print(f'{name} & {total} & {detected} & {rate:.1f}\\% & {unique} \\\\')
    elif name == 'Complex Branched':
        if not isinstance(data, list):
            continue
        
        total_molecules = 0
        detected = 0
        group_counts = Counter()
        
        for test in data:
            molecules = test.get('molecules', [])
            total_molecules += len(molecules)
            
            for mol in molecules:
                if mol.get('pfasgroups_detected', False):
                    detected += 1
                    groups = mol.get('pfasgroups_groups', [])
                    for gid in groups:
                        group_counts[gid] += 1
        
        unique = len(group_counts)
        rate = 100*detected/total_molecules if total_molecules > 0 else 0
        print(f'{name} & {total_molecules} & {detected} & {rate:.1f}\\% & {unique} \\\\')
    else:
        # Handle both dict and list formats
        if isinstance(data, list):
            molecules = data
        else:
            molecules = data.get('molecules', [])
        
        total = len(molecules)
        detected = 0
        group_counts = Counter()
        
        for mol in molecules:
            pfas_results = mol.get('pfasgroups_result', mol.get('pfasgroups_results', {}))
            if isinstance(pfas_results, dict):
                detected_groups = pfas_results.get('detected_groups', [])
                matches = pfas_results.get('matches', [])
                
                if detected_groups or matches:
                    detected += 1
                    for gid in detected_groups:
                        group_counts[gid] += 1
                    for match in matches:
                        if match.get('type') in ['PFASgroup', 'group']:
                            group_counts[match['id']] += 1
        
        unique = len(group_counts)
        rate = 100*detected/total if total > 0 else 0
        print(f'{name} & {total} & {detected} & {rate:.1f}\\% & {unique} \\\\')
