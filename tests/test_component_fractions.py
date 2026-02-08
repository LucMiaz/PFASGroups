"""Test component fraction calculations based on carbon atoms."""

from rdkit import Chem
from PFASgroups import parse_smiles

def test_component_fractions():
    """Test that component_fraction uses carbon atoms and total_components_fraction is valid."""
    
    # Test molecule: CF3CF2CF2CF2COOH (4-carbon perfluorinated chain with carboxylic acid)
    smiles = "C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
    
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    assert total_carbons > 0
    
    # Count atoms by type
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    assert atom_counts
    
    results = parse_smiles(smiles, bycomponent=True, output_format='list')
    
    assert results[0]['matches']
    
    for match in results[0]['matches']:
        if match.get('match_count', 0) == 0:
            continue
            
        group_name = match.get('group_name', match.get('definition_name', 'Unknown'))
        assert match['num_components'] >= 1
        
        # Summary metrics
        mean_comp_frac = match.get('mean_component_fraction', 0)
        total_comp_frac = match.get('total_components_fraction', 0)
        assert mean_comp_frac >= 0
        assert 0 <= total_comp_frac <= 1
        
        # Individual components
        for i, comp in enumerate(match['components'], 1):
            assert comp['size'] > 0
            comp_fraction = comp.get('component_fraction', 0)
            assert comp_fraction >= 0

            component_atoms = comp.get('component', [])
            smarts_matches = comp.get('smarts_matches', [])
            smarts_extra_atoms = comp.get('smarts_extra_atoms', 0)

            component_carbons = sum(
                1
                for atom_idx in component_atoms
                if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C'
            )
            smarts_carbons_not_in_component = 0
            for atom_idx in smarts_matches or []:
                if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C' and atom_idx not in component_atoms:
                    smarts_carbons_not_in_component += 1

            expected_fraction = (component_carbons + smarts_carbons_not_in_component + smarts_extra_atoms) / total_carbons
            assert abs(comp_fraction - expected_fraction) <= 1e-6
    
    # Total fraction should be >= max component fraction
    for match in results[0]['matches']:
        if match.get('match_count', 0) == 0:
            continue
        components = match.get('components', [])
        if not components:
            continue
        max_component_fraction = max(c.get('component_fraction', 0) for c in components)
        total_fraction = match.get('total_components_fraction', 0)
        # Total fraction is capped at 1.0; components can exceed 1.0 with extra SMARTS carbons.
        assert total_fraction >= min(1.0, max_component_fraction)
        assert total_fraction <= 1
        break
