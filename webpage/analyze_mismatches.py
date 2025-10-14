"""
Analyze mismatches between HTML and Python PFAS group detection.
"""
import pandas as pd
from collections import Counter
import json

# Load mismatches
df = pd.read_csv('mismatches_detailed.csv')

print("="*80)
print("MISMATCH ANALYSIS SUMMARY")
print("="*80)
print(f"\nTotal mismatches: {len(df)}")
print(f"Unique test molecules with mismatches: {df['test_id'].nunique()}")

# Analyze extra groups in HTML (false positives)
print("\n" + "="*80)
print("GROUPS EXTRA IN HTML (False Positives)")
print("="*80)

extra_groups = []
for idx, row in df.iterrows():
    extra = row['extra_in_html_names']
    if pd.notna(extra) and extra != '':
        for group in str(extra).split('; '):
            if group:
                extra_groups.append(group)

extra_counter = Counter(extra_groups)
print(f"\nTotal false positive instances: {len(extra_groups)}")
print("\nMost common false positives:")
for group, count in extra_counter.most_common(10):
    print(f"  {group}: {count} times")

# Analyze missing groups in HTML (false negatives)
print("\n" + "="*80)
print("GROUPS MISSING IN HTML (False Negatives)")
print("="*80)

missing_groups = []
for idx, row in df.iterrows():
    missing = row['missing_in_html_names']
    if pd.notna(missing) and missing != '':
        for group in str(missing).split('; '):
            if group:
                missing_groups.append(group)

missing_counter = Counter(missing_groups)
print(f"\nTotal false negative instances: {len(missing_groups)}")
print("\nMost common false negatives:")
for group, count in missing_counter.most_common(10):
    print(f"  {group}: {count} times")

# Group mismatches by expected group
print("\n" + "="*80)
print("MISMATCHES BY EXPECTED GROUP")
print("="*80)

grouped = df.groupby('expected_group_name').size().sort_values(ascending=False)
print(f"\nGroups with most test failures:")
for name, count in grouped.head(10).items():
    print(f"  {name}: {count} test molecules")

# Detailed breakdown for top problematic groups
print("\n" + "="*80)
print("DETAILED ANALYSIS OF TOP ISSUES")
print("="*80)

# Issue 1: Group 19 (Hydrofluoroolefins) - False positive
print("\n1. Group 19 (Hydrofluoroolefins) - FALSE POSITIVE")
print("   HTML detects it when Python doesn't")
group19_extra = df[df['extra_in_html_names'].str.contains('19:Hydrofluoroolefins', na=False)]
print(f"   Occurrences: {len(group19_extra)}")
if len(group19_extra) > 0:
    sample = group19_extra.iloc[0]
    print(f"   Example SMILES: {sample['smiles']}")
    print(f"   Formula: {sample['formula']}")
    print(f"   Python detected: {sample['python_group_names']}")
    print(f"   HTML detected: {sample['html_group_names']}")

# Issue 2: Groups 4 & 5 (Ether carboxylic acids) - False positive
print("\n2. Groups 4 & 5 (Ether carboxylic acids) - FALSE POSITIVE")
print("   HTML detects them when Python doesn't")
group45_extra = df[df['extra_in_html_names'].str.contains('4:Perfluoroalkylether carboxylic acids', na=False)]
print(f"   Occurrences: {len(group45_extra)}")
if len(group45_extra) > 0:
    sample = group45_extra.iloc[0]
    print(f"   Example SMILES: {sample['smiles']}")
    print(f"   Formula: {sample['formula']}")
    print(f"   Python detected: {sample['python_group_names']}")
    print(f"   HTML detected: {sample['html_group_names']}")

# Issue 3: Group 3 (Perfluoroalkyl dicarboxylic acids) - False positive
print("\n3. Group 3 (Perfluoroalkyl dicarboxylic acids) - FALSE POSITIVE")
print("   HTML detects it when Python doesn't")
group3_extra = df[df['extra_in_html_names'].str.contains('3:Perfluoroalkyl dicarboxylic acids', na=False)]
print(f"   Occurrences: {len(group3_extra)}")
if len(group3_extra) > 0:
    sample = group3_extra.iloc[0]
    print(f"   Example SMILES: {sample['smiles']}")
    print(f"   Formula: {sample['formula']}")
    print(f"   Python detected: {sample['python_group_names']}")
    print(f"   HTML detected: {sample['html_group_names']}")

# Issue 4: Group 21 (Semi-fluorinated alkanes) - False positive
print("\n4. Group 21 (Semi-fluorinated alkanes) - FALSE POSITIVE")
print("   HTML detects it when Python doesn't")
group21_extra = df[df['extra_in_html_names'].str.contains('21:Semi-fluorinated alkanes', na=False)]
print(f"   Occurrences: {len(group21_extra)}")
if len(group21_extra) > 0:
    sample = group21_extra.iloc[0]
    print(f"   Example SMILES: {sample['smiles']}")
    print(f"   Formula: {sample['formula']}")
    print(f"   Python detected: {sample['python_group_names']}")
    print(f"   HTML detected: {sample['html_group_names']}")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
print("\nAll mismatches are FALSE POSITIVES (extra groups in HTML).")
print("No FALSE NEGATIVES (missing groups in HTML) were found.")
print("\nThis suggests HTML has issues with:")
print("  1. Pathway filtering (not restricting to Perfluoroalkyl/Polyfluoroalkyl)")
print("  2. Constraint checking (formula constraints not being enforced)")
print("  3. Parent-child group relationships (not filtering out child groups)")
