"""Test component ratio calculations"""
from rdkit import Chem
from HalogenGroups import parse_smiles

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
    
    for smiles in test_molecules:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        
        total_carbons = count_carbons(mol)
        total_atoms = mol.GetNumAtoms()
        assert total_carbons > 0
        assert total_atoms > 0
        
        # Parse with bycomponent=True
        results = parse_smiles([smiles], bycomponent=True, output_format='list')
        assert results
        
        if results and len(results) > 0:
            result = results[0]
            if 'matches' in result:
                for match in result['matches']:
                    if isinstance(match, dict) and 'group_name' in match:
                        match_count = match.get('match_count', 0)
                        if match_count > 0:
                            mean_component_fraction = match.get('mean_component_fraction', 0)
                            total_components_fraction = match.get('total_components_fraction', 0)
                            assert mean_component_fraction >= 0
                            assert 0 <= total_components_fraction <= 1
                            components = match.get('components', [])
                            if components:
                                sum_individual_fractions = 0
                                for comp in components:
                                    comp_fraction = comp.get('component_fraction', 0)
                                    sum_individual_fractions += comp_fraction
                                    assert comp_fraction >= 0
                                if len(components) > 0:
                                    mean_from_sum = sum_individual_fractions / len(components)
                                    assert abs(mean_from_sum - mean_component_fraction) <= 0.01
                                    # total_components_fraction = union coverage; individual fractions
                                    # may overlap (same carbon in multiple augmented components),
                                    # so total can exceed the sum — only validate the range [0, 1].
                                    max_component_fraction = max(c.get('component_fraction', 0) for c in components)
                                    assert total_components_fraction >= min(1.0, max_component_fraction)

if __name__ == '__main__':
    test_component_ratios()
