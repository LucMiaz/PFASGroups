#!/usr/bin/env python3
"""
Debug script to isolate the uint32 overflow issue in PFAS-atlas
"""

import sys
sys.path.insert(0, '/home/luc/git/PFAS-atlas')

def test_rdkit_basic():
    """Test basic RDKit operations"""
    try:
        from rdkit import Chem
        print("✓ RDKit basic import successful")
        
        test_smiles = "FC(F)(F)C(F)(F)F"
        mol = Chem.MolFromSmiles(test_smiles)
        print(f"✓ Mol creation successful: {mol}")
        return True
    except Exception as e:
        print(f"✗ RDKit basic test failed: {e}")
        return False

def test_atom_count():
    """Test the atom counting function"""
    try:
        from classification_helper.atom_count import count_Atom
        test_smiles = "FC(F)(F)C(F)(F)F"
        result = count_Atom(test_smiles)
        print(f"✓ Atom count successful: {result}")
        return True
    except Exception as e:
        print(f"✗ Atom count failed: {e}")
        return False

def test_mol_converter():
    """Test the mol converter functions"""
    try:
        from rdkit_helper import mol_from_smiles, rdkit_smiles_from_mol
        test_smiles = "FC(F)(F)C(F)(F)F"
        
        print("Testing mol_from_smiles...")
        mol = mol_from_smiles(test_smiles)
        print(f"✓ mol_from_smiles successful: {mol}")
        
        print("Testing rdkit_smiles_from_mol...")
        result_smiles = rdkit_smiles_from_mol(mol)
        print(f"✓ rdkit_smiles_from_mol successful: {result_smiles}")
        return True
    except Exception as e:
        print(f"✗ Mol converter failed: {e}")
        return False

def test_mhfp_step_by_step():
    """Test MHFP encoding step by step"""
    try:
        print("Testing MHFP encoder...")
        from mhfp.encoder import MHFPEncoder
        print("✓ MHFP encoder import successful")
        
        encoder = MHFPEncoder()
        print("✓ MHFP encoder creation successful")
        
        from rdkit import Chem
        test_smiles = "FC(F)(F)C(F)(F)F"
        mol = Chem.MolFromSmiles(test_smiles)
        print(f"✓ Molecule created: {mol}")
        
        print("Attempting MHFP encoding...")
        mhfps = encoder.encode_mol(mol)
        print(f"✓ MHFP encoding successful: {type(mhfps)}, length: {len(mhfps) if hasattr(mhfps, '__len__') else 'N/A'}")
        return True
    except Exception as e:
        print(f"✗ MHFP encoding failed: {e}")
        return False

def test_calculate_mhfp():
    """Test the calculate_MHFP function directly"""
    try:
        from calculate_mhfp import calculate_MHFP
        test_smiles = "FC(F)(F)C(F)(F)F"
        
        print("Testing calculate_MHFP function...")
        result = calculate_MHFP(test_smiles)
        print(f"✓ calculate_MHFP successful: {type(result)}, length: {len(result) if hasattr(result, '__len__') else 'N/A'}")
        return True
    except Exception as e:
        print(f"✗ calculate_MHFP failed: {e}")
        return False

def test_full_classification():
    """Test the full classification pipeline"""
    try:
        from classification_helper import classify_pfas_molecule
        test_smiles = "FC(F)(F)C(F)(F)F"
        
        print("Testing full classification...")
        result = classify_pfas_molecule(test_smiles)
        print(f"✓ Full classification successful: {result}")
        return True
    except Exception as e:
        print(f"✗ Full classification failed: {e}")
        return False

def main():
    """Run all tests to isolate the issue"""
    print("🔍 Debugging PFAS-atlas uint32 overflow issue")
    print("=" * 50)
    
    tests = [
        ("RDKit Basic", test_rdkit_basic),
        ("Atom Count", test_atom_count),
        ("Mol Converter", test_mol_converter),
        ("MHFP Step-by-step", test_mhfp_step_by_step),
        ("Calculate MHFP", test_calculate_mhfp),
        ("Full Classification", test_full_classification)
    ]
    
    for test_name, test_func in tests:
        print(f"\n📋 Testing {test_name}...")
        try:
            success = test_func()
            if not success:
                print(f"❌ {test_name} failed - stopping here")
                break
        except Exception as e:
            print(f"❌ {test_name} crashed: {e}")
            break
        print(f"✅ {test_name} passed")
    
    print("\n🏁 Debug session complete")

if __name__ == "__main__":
    main()