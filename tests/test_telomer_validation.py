import sys
print("Starting telomer validation test...")

try:
    from pathlib import Path
    print("✓ Path imported")
    
    from rdkit import Chem
    print("✓ RDKit imported")
    
    # Test SDF reading
    sdf_path = Path('benchmark/data/PubChem_fluorotelomers.sdf')
    print(f"Looking for: {sdf_path.absolute()}")
    print(f"Exists: {sdf_path.exists()}")
    
    if sdf_path.exists():
        supplier = Chem.SDMolSupplier(str(sdf_path))
        mols = []
        for i, mol in enumerate(supplier):
            if mol is not None:
                mols.append(mol)
        print(f"✓ Found {len(mols)} molecules in SDF file")
        
        # Test PFASgroups import
        from PFASgroups import parse_smiles
        print("✓ PFASgroups imported")
        
        # Test on first molecule
        if mols:
            smiles = Chem.MolToSmiles(mols[0])
            print(f"\nTesting first molecule: {smiles[:50]}...")
            result = parse_smiles(smiles)
            print(f"✓ Parse successful, got {len(result) if result else 0} results")
        
        print("\n✅ All tests passed! Script should work.")
    else:
        print("✗ SDF file not found!")

except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
