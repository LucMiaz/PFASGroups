"""
Test the new summary methods for ResultsModel and MoleculeResult
"""
from HalogenGroups import parse_smiles

# Test molecule - perfluorooctane sulfonic acid
test_smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"

print("Testing summary methods")
print("=" * 80)

# Parse the molecule
results = parse_smiles(test_smiles)

# Test MoleculeResult.summary()
print("\n1. Testing MoleculeResult.summary():")
print("-" * 80)
results[0].summary()

# Test ResultsModel.summary()
print("\n2. Testing ResultsModel.summary():")
print("-" * 80)
results.summary()

# Test with multiple molecules
print("\n3. Testing with multiple molecules:")
print("-" * 80)
multi_smiles = [
    "FC(F)(F)C(F)(F)C(=O)O",  # PFPA
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFBA
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
]

multi_results = parse_smiles(multi_smiles)
multi_results.summary()

print("\nDone!")
