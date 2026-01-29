import sys
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')

print("Script started")

try:
    print("Importing PFASgroups...")
    from PFASgroups import parse_smiles
    print("✓ Import successful")
    
    smiles = 'FC(F)(F)C(=O)O'
    
    print(f"\nTest 1: Default parameters")
    result = parse_smiles(smiles)
    print(f"✓ Success: {len(result)} results")
    
    print(f"\nTest 2: limit_effective_graph_resistance=0")
    result = parse_smiles(smiles, limit_effective_graph_resistance=0)
    print(f"✓ Success: {len(result)} results")
    
    print(f"\nTest 3: limit_effective_graph_resistance=None")
    result = parse_smiles(smiles, limit_effective_graph_resistance=None)
    print(f"✓ Success: {len(result)} results")
    
    print(f"\nTest 4: compute_component_metrics=False")
    result = parse_smiles(smiles, compute_component_metrics=False)
    print(f"✓ Success: {len(result)} results")
    
    print("\n✓ All tests passed!")
    
except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
