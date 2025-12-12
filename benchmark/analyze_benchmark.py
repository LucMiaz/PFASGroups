import pandas as pd
import json

# Load the benchmark results
df = pd.read_csv('direct_benchmark_results.csv')

print('PFASGroups vs PFAS-atlas Benchmark Analysis')
print('=' * 60)

print(f'Total molecules analyzed: {len(df)}')
pfas_groups_coverage = (df['pfasgroups_classified'] == True).sum()
atlas_coverage = (df['atlas_classified'] == True).sum()
print(f'PFASGroups coverage: {pfas_groups_coverage}/{len(df)} ({pfas_groups_coverage/len(df)*100:.1f}%)')
print(f'PFAS-atlas coverage: {atlas_coverage}/{len(df)} ({atlas_coverage/len(df)*100:.1f}%)')

print('\n📊 PFAS-atlas Classification Results:')
atlas_counts = df['atlas_first_class'].value_counts()
for class_name, count in atlas_counts.items():
    print(f'  {class_name}: {count} ({count/len(df)*100:.1f}%)')

# Check agreement between systems
print('\n🔍 System Agreement Analysis:')
# Note: This is a complex comparison since the systems use different classification schemes
# PFAS-atlas uses broad categories, PFASGroups uses specific functional groups

# Count molecules where PFASGroups found groups vs PFAS-atlas classification
pfas_groups_detected = df['n_detected_groups'] > 0
pfas_atlas_pfas = ~df['atlas_first_class'].isin(['Not PFAS'])

detected_count = pfas_groups_detected.sum()
atlas_pfas_count = pfas_atlas_pfas.sum()
print(f'PFASGroups detected PFAS groups: {detected_count}/{len(df)} ({detected_count/len(df)*100:.1f}%)')
print(f'PFAS-atlas classified as PFAS: {atlas_pfas_count}/{len(df)} ({atlas_pfas_count/len(df)*100:.1f}%)')

# Cross-tabulation
print('\n📈 Cross-tabulation (PFASGroups groups detected vs PFAS-atlas classification):')
crosstab = pd.crosstab(pfas_groups_detected, df['atlas_first_class'], margins=True)
print(crosstab)

# Look at specific examples where systems might disagree
print('\n🔬 Examples of Different Classifications:')
not_pfas_by_atlas = df[df['atlas_first_class'] == 'Not PFAS']
if len(not_pfas_by_atlas) > 0:
    print(f'\nMolecules classified as "Not PFAS" by PFAS-atlas ({len(not_pfas_by_atlas)} total):')
    for idx, row in not_pfas_by_atlas.head(5).iterrows():
        print(f'  SMILES: {row["smiles"]}')
        print(f'  PFASGroups: {row["detected_groups"]} groups detected')
        print(f'  Origin: {row["origin"]}')
        print()

# Look at specificity
print('📊 PFASGroups Specificity Analysis:')
specific_count = (df['is_specific'] == True).sum()
print(f'Molecules with specific group detection: {specific_count}/{len(df)} ({specific_count/len(df)*100:.1f}%)')

# Group distribution
print('\n📋 Most Common PFASGroups Detected:')
all_groups = []
for groups_str in df['detected_groups'].dropna():
    if isinstance(groups_str, str) and groups_str.startswith('['):
        try:
            groups = eval(groups_str)
            all_groups.extend(groups)
        except:
            pass

if all_groups:
    from collections import Counter
    group_counts = Counter(all_groups)
    print('Top 10 most detected groups:')
    for group_id, count in group_counts.most_common(10):
        print(f'  Group {group_id}: {count} molecules ({count/len(df)*100:.1f}%)')

print('\n✅ Analysis complete!')