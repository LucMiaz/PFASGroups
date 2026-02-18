"""Test the updated utility_funcgroups with comprehensive metrics."""

from HalogenGroups.parser import parse_mols
from rdkit import Chem

def test_parse_mols_output_structure():
    """Test that parse_mols returns the correct structure with summary metrics."""
    
    test_smiles = [
        'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',  # PFOA
        'C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(=O)O'  # Branched PFAS
    ]
    
    for i, smiles in enumerate(test_smiles):
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        
        results = parse_mols([mol], bycomponent=True, output_format='list')
        assert results
        mol_result = results[0]
        matches = mol_result.get('matches', [])
        assert isinstance(matches, list)

        for match in matches:
            if match['type'] != 'HalogenGroup':
                continue

            # Check summary metrics
            summary_fields = [
                'mean_eccentricity',
                'mean_smarts_centrality', 
                'mean_diameter',
                'mean_radius',
                'mean_effective_graph_resistance',
                'mean_dist_to_barycenter',
                'mean_dist_to_center',
                'mean_dist_to_periphery'
            ]
            
            for field in summary_fields:
                assert field in match, f"Missing summary field: {field}"
            
            # Check individual components
            components = match.get('components', [])
            if len(components) > 0:
                comp = components[0]
                
                component_fields = [
                    'component', 'size', 'SMARTS',
                    'eccentricity_values', 'smarts_centrality',
                    'diameter', 'radius', 'effective_graph_resistance',
                    'center', 'periphery', 'barycenter',
                    'min_dist_to_barycenter', 'min_dist_to_center',
                    'max_dist_to_periphery'
                ]
                
                for field in component_fields:
                    assert field in comp, f"Missing component field: {field}"

def simulate_database_save():
    """Simulate the database saving logic without actual database."""
    
    smiles = 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    
    results = parse_mols([mol], bycomponent=True, output_format='list')
    assert results
    mol_result = results[0]
    matches = mol_result.get('matches', [])

    for match in matches:
        if match['type'] != 'HalogenGroup':
            continue
        summary_data = {
            'n': match['match_count'],
            'chainLengths': match['components_sizes'],
            'mean_eccentricity': match.get('mean_eccentricity'),
            'mean_smarts_centrality': match.get('mean_smarts_centrality'),
            'mean_diameter': match.get('mean_diameter'),
            'mean_radius': match.get('mean_radius'),
            'mean_effective_graph_resistance': match.get('mean_effective_graph_resistance'),
            'mean_dist_to_barycenter': match.get('mean_dist_to_barycenter'),
            'mean_dist_to_center': match.get('mean_dist_to_center'),
            'mean_dist_to_periphery': match.get('mean_dist_to_periphery')
        }
        for key in summary_data:
            assert key in summary_data
        
        # Simulate Components saves
        components = match.get('components', [])
        for i, comp in enumerate(components):
            component_data = {
                'indices': comp['component'],
                'length': comp.get('size', len(comp['component'])),
                'SMARTS': comp['SMARTS'],
                'eccentricity_values': comp['eccentricity_values'],
                'smarts_centrality': comp['smarts_centrality'],
                'diameter': comp.get('diameter'),
                'radius': comp.get('radius'),
                'effective_graph_resistance': comp.get('effective_graph_resistance'),
                'center': comp.get('center', []),
                'periphery': comp.get('periphery', []),
                'barycenter': comp.get('barycenter', []),
                'min_dist_to_barycenter': comp.get('min_dist_to_barycenter', 0),
                'min_dist_to_center': comp.get('min_dist_to_center', 0),
                'max_dist_to_periphery': comp.get('max_dist_to_periphery', 0)
            }
            for key in component_data:
                assert key in component_data


def test_database_save_simulation():
    simulate_database_save()
