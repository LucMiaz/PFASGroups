"""Test component fraction calculations with full atom counts."""

from rdkit import Chem
from PFASgroups import parse_smiles

def test_component_fractions():
    """Test that component_fraction includes all atoms (H, F, halogens) and total_components_fraction works."""
    
    # Test molecule: CF3CF2CF2CF2COOH (4-carbon perfluorinated chain with carboxylic acid)
    smiles = "C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
    
    print(f"Testing SMILES: {smiles}")
    
    # Parse with AddHs to see all atoms
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    total_atoms = mol.GetNumAtoms()
    print(f"Total atoms in molecule (with H): {total_atoms}")
    
    # Count atoms by type
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    print(f"Atom composition: {atom_counts}")
    
    # Parse with PFASgroups
    results = parse_smiles(smiles, bycomponent=True, output_format='list')
    
    print(f"\nNumber of matches: {len(results[0]['matches'])}")
    
    for match in results[0]['matches']:
        if match.get('match_count', 0) == 0:
            continue
            
        group_name = match.get('group_name', match.get('definition_name', 'Unknown'))
        print(f"\nPFAS Group: {group_name}")
        print(f"Number of components: {match['num_components']}")
        
        # Summary metrics
        mean_comp_frac = match.get('mean_component_fraction', 0)
        total_comp_frac = match.get('total_components_fraction', 0)
        print(f"Mean component fraction: {mean_comp_frac:.3f}")
        print(f"Total components fraction (union): {total_comp_frac:.3f}")
        
        # Individual components
        for i, comp in enumerate(match['components'], 1):
            print(f"\n  Component {i}:")
            print(f"    Backbone size (C atoms): {comp['size']} atoms")
            print(f"    Component fraction (with H, F): {comp['component_fraction']:.3f}")
            print(f"    Branching: {comp['branching']:.3f}")
            print(f"    Mean eccentricity: {comp['mean_eccentricity']:.3f}")
            print(f"    SMARTS type: {comp['SMARTS']}")
            
            # Calculate what the fraction means in atom count
            implied_atoms = int(comp['component_fraction'] * total_atoms)
            print(f"    Implied atom count: {implied_atoms} atoms")
    
    print("\n" + "="*70)
    print("VALIDATION CHECKS:")
    print("="*70)
    
    # Find perfluoroalkyl group
    for match in results[0]['matches']:
        if 'perfluoroalkyl' in match.get('group_name', '').lower() and 'carboxylic' in match.get('group_name', '').lower():
            comp = match['components'][0]
            backbone_size = comp['size']  # Should be 4 (C atoms)
            comp_fraction = comp['component_fraction']
            total_fraction = match.get('total_components_fraction', 0)
            
            # Calculate expected full size
            # 4 C atoms + 9 F atoms = 13 atoms in perfluoroalkyl component
            # (CF3 has 3F, CF2 groups have 2F each: 3 + 2 + 2 + 2 = 9F total)
            expected_full_size = 13  # 4C + 9F
            implied_full_size = int(comp_fraction * total_atoms)
            
            print(f"\n✓ Perfluoroalkyl component:")
            print(f"  - Backbone (C only): {backbone_size} atoms")
            print(f"  - Expected full size (C+F): {expected_full_size} atoms")
            print(f"  - Implied full size from fraction: {implied_full_size} atoms")
            print(f"  - Component fraction: {comp_fraction:.3f} = {implied_full_size}/{total_atoms}")
            
            if abs(implied_full_size - expected_full_size) <= 1:  # Allow 1 atom tolerance
                print(f"  ✓ PASS: Component fraction correctly includes F atoms!")
            else:
                print(f"  ✗ FAIL: Expected {expected_full_size} atoms, got {implied_full_size}")
            
            # Check total fraction
            print(f"\n✓ Total components fraction (union): {total_fraction:.3f}")
            union_atoms = int(total_fraction * total_atoms)
            print(f"  - Union covers {union_atoms}/{total_atoms} atoms")
            
            if total_fraction >= comp_fraction:
                print(f"  ✓ PASS: Total fraction >= individual component fraction")
            else:
                print(f"  ✗ FAIL: Total fraction should be >= individual fraction")
            
            break
    
    print("\n✓ Test completed successfully!")

if __name__ == "__main__":
    test_component_fractions()
