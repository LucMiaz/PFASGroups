import pandas as pd

# Load the specificity test results
df = pd.read_csv('specificity_test_results.csv')

# Find cases with more than 5 detected groups (low specificity)
low_spec = df[(df['valid_smiles'] == True) & (df['n_detected_groups'] > 5)]

print('TOP 10 LOW SPECIFICITY CASES:')
print('=' * 80)

for i, (_, row) in enumerate(low_spec.nlargest(10, 'n_detected_groups').iterrows()):
    print(f'{i+1}. Expected Groups: {row["group_ids"]}')
    print(f'   Detected Groups: {row["detected_groups"]}')
    print(f'   Total Detected: {row["n_detected_groups"]}')
    print(f'   SMILES: {row["smiles"][:100]}{"..." if len(row["smiles"]) > 100 else ""}')
    print()

print(f'\nSUMMARY:')
print(f'Total low specificity cases (>5 groups): {len(low_spec)}')
print(f'Worst case detected {low_spec["n_detected_groups"].max()} groups')
print(f'Average groups in low specificity cases: {low_spec["n_detected_groups"].mean():.1f}')

# Find the single false negative case
false_negatives = df[(df['valid_smiles'] == True) & (df['expected_group_detected'] == False)]
print(f'\nFALSE NEGATIVE CASES: {len(false_negatives)}')
if len(false_negatives) > 0:
    for i, (_, row) in enumerate(false_negatives.iterrows()):
        print(f'{i+1}. Expected Groups: {row["group_ids"]}')
        print(f'   Detected Groups: {row["detected_groups"]}')
        print(f'   SMILES: {row["smiles"][:100]}{"..." if len(row["smiles"]) > 100 else ""}')
        print()