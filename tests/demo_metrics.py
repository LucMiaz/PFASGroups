"""
Demonstration of component eccentricity and SMARTS centrality metrics.

These metrics provide insights into the molecular structure:
1. Eccentricity: Measures how linear vs branched the fluorinated component is
2. SMARTS Centrality: Measures how central the functional group is within the component
"""

from rdkit import Chem
from PFASgroups.parser import parse_groups_in_mol

def analyze_molecule(smiles, description):
    """Analyze a molecule and print its metrics."""
    print(f"\n{description}")
    print(f"SMILES: {smiles}")
    print("-" * 70)
    
    mol = Chem.MolFromSmiles(smiles)
    matches = parse_groups_in_mol(mol)
    
    for pfas_group, count, comp_lengths, components in matches:
        if components:
            print(f"\nPFAS Group: {pfas_group.name}")
            for i, comp in enumerate(components[:3]):  # Show up to 3 components
                print(f"  Component {i+1}:")
                print(f"    Length: {comp['length']} atoms")
                print(f"    Type: {comp['SMARTS']}")
                print(f"    Eccentricity: {comp['eccentricity']:.3f} ", end="")
                if comp['eccentricity'] > 0.9:
                    print("(linear)")
                elif comp['eccentricity'] > 0.7:
                    print("(slightly branched)")
                elif comp['eccentricity'] > 0.5:
                    print("(moderately branched)")
                else:
                    print("(highly branched)")
                
                print(f"    SMARTS Centrality: {comp['smarts_centrality']:.3f} ", end="")
                if comp['smarts_centrality'] < 0.2:
                    print("(peripheral)")
                elif comp['smarts_centrality'] < 0.5:
                    print("(off-center)")
                elif comp['smarts_centrality'] < 0.8:
                    print("(near center)")
                else:
                    print("(central)")
            
            if len(components) > 3:
                print(f"  ... and {len(components) - 3} more components")
            break  # Only show first PFAS group

if __name__ == "__main__":
    print("=" * 70)
    print("COMPONENT METRICS DEMONSTRATION")
    print("=" * 70)
    
    # Test cases showing different structural features
    test_cases = [
        # Linear, terminal functional group
        ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
         "Test 1: Long linear chain, terminal COOH"),
        
        # Linear, central functional group
        ("FC(F)(F)C(F)(F)C(F)(F)C(=O)C(F)(F)C(F)(F)C(F)(F)F",
         "Test 2: Linear chain, central ketone"),
        
        # Branched, terminal functional group
        ("FC(F)(C(F)(F)F)C(F)(C(F)(F)F)C(F)(C(F)(F)F)C(=O)O",
         "Test 3: Highly branched, terminal COOH"),
        
        # Branched, more central functional group
        ("FC(F)(F)C(F)(C(F)(F)F)C(=O)C(F)(C(F)(F)F)C(F)(F)F",
         "Test 4: Branched with central ketone"),
        
        # PFOS - sulfonic acid
        ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",
         "Test 5: PFOS - linear with terminal sulfonic acid"),
        
        # Fluorotelomer alcohol
        ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)CCO",
         "Test 6: Fluorotelomer alcohol (6:2 FTOH)"),
    ]
    
    for smiles, description in test_cases:
        analyze_molecule(smiles, description)
    
    print("\n" + "=" * 70)
    print("METRICS GUIDE")
    print("=" * 70)
    print("\nEccentricity (measures branching):")
    print("  1.0 = Perfectly linear (no branch points)")
    print("  0.7-0.9 = Slightly branched")
    print("  0.5-0.7 = Moderately branched")
    print("  0.0-0.5 = Highly branched")
    print("\nSMARTS Centrality (functional group position):")
    print("  1.0 = At exact center of component")
    print("  0.5-1.0 = Near center")
    print("  0.2-0.5 = Off-center")
    print("  0.0-0.2 = At periphery (end of component)")
    print("\nApplications:")
    print("  - Filter for linear vs branched PFAS")
    print("  - Identify terminal vs internal functional groups")
    print("  - Structural complexity analysis")
    print("  - Database clustering by structural features")
    print("=" * 70)
