from rdkit import Chem
from PFASGroups.parser import parse_groups_in_mol


def _first_component(mol):
    matches, *_ = parse_groups_in_mol(mol)
    for _, _, _, components in matches:
        if components:
            return components[0]
    return None


def test_basic_metrics_presence():
    mol1 = Chem.MolFromSmiles('C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O')
    mol2 = Chem.MolFromSmiles('FC(F)(C(F)(F)F)C(F)(C(F)(F)F)C(=O)O')
    mol3 = Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O')

    for mol in [mol1, mol2, mol3]:
        assert mol is not None
        comp = _first_component(mol)
        assert comp is not None
        assert 'size' in comp
        assert 'smarts_centrality' in comp
