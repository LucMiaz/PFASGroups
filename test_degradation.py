"""
Test script for the new degradation product functions in fragmentation.py

This script demonstrates how to use the new functions to generate degradation products
for PFAS molecules by breaking bonds outside the fluorinated chain.
"""

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import sys
import os

# Add the PFASgroups directory to the path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import the fragmentation module
try:
    from PFASgroups.fragmentation import (
        find_fluorinated_chains,
        get_non_fluorinated_bonds, 
        generate_degradation_products,
        analyze_degradation_pathways
    )
    print("Successfully imported degradation functions!")
except ImportError as e:
    print(f"Import error: {e}")
    # For testing, we can import directly from the file
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "fragmentation", 
        "c:\\Users\\luc\\git\\PFASgroups\\PFASgroups\\fragmentation.py"
    )
    fragmentation = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fragmentation)
    
    find_fluorinated_chains = fragmentation.find_fluorinated_chains
    get_non_fluorinated_bonds = fragmentation.get_non_fluorinated_bonds
    generate_degradation_products = fragmentation.generate_degradation_products
    analyze_degradation_pathways = fragmentation.analyze_degradation_pathways
    print("Successfully imported from file!")

def test_degradation_functions():
    """Test the degradation product generation functions"""
    
    # Test molecules (PFAS examples)
    test_smiles = [
        "C(C(C(C(F)(F)F)(F)F)(F)F)(C(=O)O)(F)F",  # PFOA-like molecule
        "C(C(F)(F)F)(F)(F)C(=O)O",  # Shorter perfluoroalkyl acid
        "C(C(C(F)(F)F)(F)F)(F)(F)N",  # Perfluoroalkyl amine
    ]
    
    for i, smiles in enumerate(test_smiles):
        print(f"\n{'='*60}")
        print(f"Testing molecule {i+1}: {smiles}")
        print(f"{'='*60}")
        
        # Create molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Could not parse SMILES: {smiles}")
            continue
            
        print(f"Formula: {CalcMolFormula(mol)}")
        print(f"Number of atoms: {mol.GetNumAtoms()}")
        print(f"Number of bonds: {mol.GetNumBonds()}")
        
        # Find fluorinated chains
        fluorinated_atoms = find_fluorinated_chains(mol)
        print(f"Fluorinated chain atoms: {len(fluorinated_atoms)} atoms")
        print(f"Atom indices: {sorted(fluorinated_atoms)}")
        
        # Get breakable bonds
        breakable_bonds = get_non_fluorinated_bonds(mol, fluorinated_atoms)
        print(f"Breakable bonds (not in fluorinated chain): {len(breakable_bonds)}")
        print(f"Bond indices: {breakable_bonds}")
        
        # Generate degradation products
        print("\nGenerating degradation products...")
        degradation_products = generate_degradation_products(mol, max_breaks=2)
        
        print(f"Total degradation products: {len(degradation_products)}")
        
        # Show some examples
        count = 0
        for inchi_key, formulas in degradation_products.items():
            if count >= 5:  # Limit output
                break
            for formula, data in formulas.items():
                frag_mol = data['mol']
                print(f"  - Formula: {formula}, SMILES: {Chem.MolToSmiles(frag_mol)}, "
                      f"Breaks: {data['n_breaks']}")
                count += 1
        
        if len(degradation_products) > 5:
            print(f"  ... and {len(degradation_products) - 5} more products")
        
        # Analyze degradation pathways
        print("\nAnalyzing degradation pathways...")
        analysis = analyze_degradation_pathways(mol, max_breaks=2)
        
        print(f"Original molecule weight: {analysis['original_molecule']['molecular_weight']:.2f}")
        print(f"Products by number of breaks:")
        for n_breaks, count in analysis['degradation_summary']['products_by_breaks'].items():
            print(f"  {n_breaks} breaks: {count} products")

if __name__ == "__main__":
    test_degradation_functions()