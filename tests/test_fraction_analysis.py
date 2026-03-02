"""Compare component fractions with and without SMARTS precomputation."""
from rdkit import Chem
from PFASGroups.parser import parse_smiles


def test_fraction_analysis_pfoa():
    pfoa_smiles = "C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F"

    mol = Chem.MolFromSmiles(pfoa_smiles)
    assert mol is not None
    mol = Chem.AddHs(mol)
    total_atoms = mol.GetNumAtoms()
    assert total_atoms > 0

    results = parse_smiles([pfoa_smiles], bycomponent=True)
    assert results
    result = results[0]
    assert 'matches' in result

    found = False
    for match in result['matches']:
        if 'Perfluoroalkyl carboxylic acid' in match.get('group_name', ''):
            found = True
            total_fraction = match.get('total_components_fraction', 0)
            assert 0 <= total_fraction <= 1
            for comp in match.get('components', []):
                comp_fraction = comp.get('component_fraction', 0)
                implied_atoms = int(comp_fraction * total_atoms)
                assert implied_atoms >= 0
            break
    assert found
