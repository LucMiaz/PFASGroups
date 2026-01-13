"""Test component ratio calculations"""
from rdkit import Chem
from PFASgroups import parse_smiles, get_PFASGroups
import json

def count_carbons(mol):
    """Count the number of carbon atoms in a molecule"""
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

def test_component_ratios():
    """Test that component ratios are computed correctly based on carbon atoms"""
    
    test_molecules = [
        # Perfluoroalkyl carboxylic acid
        "C(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)C(=O)O",
        # PFOA
        "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",
        # Simple perfluorinated compound
        "FC(F)(F)C(F)(F)C(F)(F)C(=O)O",
    ]
    
    print("=" * 100)
    print("TESTING COMPONENT RATIO CALCULATIONS")
    print("=" * 100)
    
    for smiles in test_molecules:
        print(f"\nTest molecule: {smiles}")
        mol = Chem.MolFromSmiles(smiles)
        
        total_carbons = count_carbons(mol)
        total_atoms = mol.GetNumAtoms()
        
        print(f"  Total carbons: {total_carbons}")
        print(f"  Total atoms: {total_atoms}")
        
        # Parse with bycomponent=True
        results = parse_smiles([smiles], bycomponent=True, output_format='dict')
        
        if results and len(results) > 0:
            result = results[0]
            if 'matches' in result:
                for match in result['matches']:
                    if isinstance(match, dict) and 'group_name' in match:
                        group_name = match['group_name']
                        match_count = match.get('match_count', 0)
                        
                        if match_count > 0:
                            print(f"\n  Group: {group_name}")
                            print(f"    Match count: {match_count}")
                            
                            # Check component fractions
                            mean_component_fraction = match.get('mean_component_fraction', 0)
                            total_components_fraction = match.get('total_components_fraction', 0)
                            
                            print(f"    Mean component fraction: {mean_component_fraction:.4f}")
                            print(f"    Total components fraction: {total_components_fraction:.4f}")
                            
                            # Check individual components
                            components = match.get('components', [])
                            if components:
                                print(f"    Number of components: {len(components)}")
                                sum_individual_fractions = 0
                                
                                for i, comp in enumerate(components):
                                    comp_fraction = comp.get('component_fraction', 0)
                                    comp_carbons = comp.get('size', 0)  # Component size = carbon atoms
                                    smarts_matches = comp.get('smarts_matches', [])
                                    
                                    # Count carbons in SMARTS matches
                                    smarts_carbons = 0
                                    if smarts_matches:
                                        for atom_idx in smarts_matches:
                                            if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                                                smarts_carbons += 1
                                    
                                    sum_individual_fractions += comp_fraction
                                    
                                    print(f"      Component {i+1}:")
                                    print(f"        Component carbons: {comp_carbons}")
                                    print(f"        SMARTS match carbons: {smarts_carbons}")
                                    print(f"        Total carbons in component: {comp_carbons + smarts_carbons}")
                                    print(f"        Component fraction: {comp_fraction:.4f}")
                                    print(f"        Expected fraction (carbons/total_carbons): {(comp_carbons + smarts_carbons)/total_carbons:.4f}")
                                
                                print(f"    Sum of individual fractions: {sum_individual_fractions:.4f}")
                                print(f"    Mean component fraction: {mean_component_fraction:.4f}")
                                print(f"    Total components fraction: {total_components_fraction:.4f}")
                                
                                # Validation
                                if sum_individual_fractions / len(components) - mean_component_fraction > 0.001:
                                    print(f"    [ERROR] Mean doesn't match sum/count")
                                else:
                                    print(f"    [OK] Mean component fraction is correct")
                                
                                if sum_individual_fractions > total_components_fraction + 0.001:
                                    print(f"    [ERROR] Sum of individual fractions > total fraction!")
                                    print(f"       This suggests overlapping components or incorrect total calculation")
                                elif total_components_fraction > 1.001:
                                    print(f"    [ERROR] Total components fraction > 1.0!")
                                else:
                                    print(f"    [OK] Total components fraction is reasonable")

if __name__ == '__main__':
    test_component_ratios()
