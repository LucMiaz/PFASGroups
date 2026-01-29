#!/usr/bin/env python
"""Check telomer validation results summary."""
import json

with open('../data/telomer_validation_results.json', 'r') as f:
    data = json.load(f)

print(f"Detection rate: {data['detection_rate']:.2f}%")
print(f"Total: {data['total_molecules']}, Detected: {data['telomer_detected']}")
print("\nGroup frequencies:")
for g in sorted(data['group_counts'], key=lambda x: x['count'], reverse=True):
    pct = 100 * g['count'] / data['total_molecules']
    print(f"  ID {g['id']:2d} ({g['name']:40s}): {g['count']:3d} ({pct:.1f}%)")
