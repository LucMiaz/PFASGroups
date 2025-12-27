"""
Test script for compile_smartsPath and compile_smartsPaths functions.
"""
from PFASgroups import (
    compile_smartsPath, 
    compile_smartsPaths, 
    get_smartsPaths, 
    get_PFASGroups,
    parse_smiles,
    PFASGroup
)

print("Testing compile_smartsPath...")
print("-" * 60)

# Test 1: Single path compilation
print("\n1. Testing single path compilation (Perchlorinated)")
perchlo_path = compile_smartsPath(
    "[C;X4;H0](Cl)(Cl)!@!=!#[C;X4;H0](Cl)(Cl)",
    "[C;X4;H0](Cl)(Cl)Cl"
)
print(f"   Result: {type(perchlo_path)} with {len(perchlo_path)} elements")
print(f"   Chain mol: {perchlo_path[0]}")
print(f"   End mol: {perchlo_path[1]}")

# Test 2: Multiple paths compilation
print("\n2. Testing multiple paths compilation")
custom_paths_dict = {
    'Perchlorinated': {
        'chain': '[C;X4;H0](Cl)(Cl)!@!=!#[C;X4;H0](Cl)(Cl)',
        'end': '[C;X4;H0](Cl)(Cl)Cl'
    },
    'Perbrominated': {
        'chain': '[C;X4;H0](Br)(Br)!@!=!#[C;X4;H0](Br)(Br)',
        'end': '[C;X4;H0](Br)(Br)Br'
    },
    'MixedHalo': {
        'chain': '[C;X4]([F,Cl])([F,Cl])!@!=!#[C;X4]([F,Cl])',
        'end': '[C;X4]([F,Cl])([F,Cl])[F,Cl]'
    }
}

compiled_paths = compile_smartsPaths(custom_paths_dict)
print(f"   Compiled {len(compiled_paths)} paths:")
for name in compiled_paths:
    print(f"   - {name}")

# Test 3: Integration with existing paths
print("\n3. Testing integration with default paths")
default_paths = get_smartsPaths()
print(f"   Default paths: {list(default_paths.keys())}")

default_paths.update(compiled_paths)
print(f"   After merge: {list(default_paths.keys())}")

# Test 4: Use custom paths in parsing
print("\n4. Testing custom paths in parsing")

# Create a custom group that uses the new path
groups = get_PFASGroups()
groups.append(PFASGroup(
    id=3000,
    name="Perchlorinated carboxylic acid",
    smarts1="[C](Cl)(Cl)Cl",
    smarts2="C(=O)O",
    smartsPath="Perchlorinated",
    constraints={"nCl": [3, None], "nC": [2, None]}
))

# Test with a chlorinated compound
test_smiles = ["ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"]
results = parse_smiles(test_smiles, smartsPaths=default_paths, pfas_groups=groups)

print(f"   Test SMILES: {test_smiles[0]}")
print(f"   Found {len(results[0])} matches:")
for group, match_count, chain_lengths, _ in results[0]:
    if group.id == 3000:
        print(f"   ✓ Matched: {group.name}")
        print(f"     Chain length: {chain_lengths}")

# Test 5: Error handling
print("\n5. Testing error handling")
try:
    bad_paths = compile_smartsPaths({
        'BadPath': {'chain': '[C](Cl)'}  # Missing 'end' key
    })
    print("   ERROR: Should have raised ValueError")
except ValueError as e:
    print(f"   ✓ Correctly raised ValueError: {e}")

print("\n" + "=" * 60)
print("All tests completed successfully!")
print("=" * 60)
