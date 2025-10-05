"""
Test script to demonstrate the difference between breaking all non-fluorinated bonds
vs. only breaking bonds connected to (but not part of) the fluorinated chain.
"""

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import sys
import os

# Add the PFASgroups directory to the path  
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_bond_breaking_strategies():
    """Test the different bond breaking strategies"""
    
    # Test molecule: PFOA-like with additional side chains
    smiles = "C(C(C(C(F)(F)F)(F)F)(F)F)(C(=O)O)(F)F"  # PFOA-like molecule
    
    print(f"Testing molecule: {smiles}")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Could not parse SMILES")
        return
        
    print(f"Formula: {CalcMolFormula(mol)}")
    print(f"Number of atoms: {mol.GetNumAtoms()}")
    print(f"Number of bonds: {mol.GetNumBonds()}")
    
    # Import the functions from the modified file
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "fragmentation", 
        "c:\\Users\\luc\\git\\PFASgroups\\PFASgroups\\fragmentation.py"
    )
    fragmentation = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fragmentation)
    
    # Test both strategies
    print("\n" + "="*60)
    print("STRATEGY 1: Breaking bonds connected to fluorinated chain")
    print("="*60)
    
    # Find fluorinated chains
    fluorinated_atoms = fragmentation.find_fluorinated_chains(mol)
    print(f"Fluorinated chain atoms: {sorted(fluorinated_atoms)}")
    
    # Get bonds connected to fluorinated chain
    connected_bonds = fragmentation.get_bonds_connected_to_fluorinated_path(mol, fluorinated_atoms)
    print(f"Bonds connected to fluorinated chain: {connected_bonds}")
    
    # Show what these bonds connect
    print("Bond details:")
    for bond_idx in connected_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        begin_in_chain = bond.GetBeginAtomIdx() in fluorinated_atoms
        end_in_chain = bond.GetEndAtomIdx() in fluorinated_atoms
        print(f"  Bond {bond_idx}: {begin_atom.GetSymbol()}{bond.GetBeginAtomIdx()} ({'F-chain' if begin_in_chain else 'other'}) - "
              f"{end_atom.GetSymbol()}{bond.GetEndAtomIdx()} ({'F-chain' if end_in_chain else 'other'})")
    
    # Generate degradation products using new strategy
    products_connected = fragmentation.generate_degradation_products(mol, max_breaks=2)
    print(f"\nDegradation products (connected bonds): {len(products_connected)}")
    
    print("\n" + "="*60)
    print("STRATEGY 2: Breaking all non-fluorinated bonds")
    print("="*60)
    
    # Get all non-fluorinated bonds
    all_non_fluorinated = fragmentation.get_non_fluorinated_bonds(mol, fluorinated_atoms)
    print(f"All non-fluorinated bonds: {all_non_fluorinated}")
    
    # Show difference
    bonds_between_non_fluorinated = [b for b in all_non_fluorinated if b not in connected_bonds]
    print(f"Bonds between non-fluorinated atoms (excluded in new strategy): {bonds_between_non_fluorinated}")
    
    # Show what these excluded bonds connect
    if bonds_between_non_fluorinated:
        print("Excluded bond details:")
        for bond_idx in bonds_between_non_fluorinated:
            bond = mol.GetBondWithIdx(bond_idx)
            begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
            end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
            print(f"  Bond {bond_idx}: {begin_atom.GetSymbol()}{bond.GetBeginAtomIdx()} - "
                  f"{end_atom.GetSymbol()}{bond.GetEndAtomIdx()} (both non-fluorinated)")
    
    print(f"\nTotal breakable bonds:")
    print(f"  Connected to fluorinated chain: {len(connected_bonds)}")
    print(f"  All non-fluorinated: {len(all_non_fluorinated)}")
    print(f"  Difference: {len(all_non_fluorinated) - len(connected_bonds)}")
    
    print(f"\nThis means the new strategy is more selective and will generate fewer, ")
    print(f"more chemically relevant degradation products by only breaking bonds that ")
    print(f"disconnect functional groups from the fluorinated backbone.")

if __name__ == "__main__":
    test_bond_breaking_strategies()