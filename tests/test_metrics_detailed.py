from rdkit import Chem
from PFASGroups.parser import parse_groups_in_mol


def _first_component(mol):
    matches, _ = parse_groups_in_mol(mol)
    for _, _, _, components in matches[:1]:
        if components:
            return components[0]
    return None


def test_metrics_detailed_components():
    mols = [
        Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O'),
        Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(=O)C(F)(F)C(F)(F)F'),
        Chem.MolFromSmiles('FC(F)(C(F)(F)F)C(F)(C(F)(F)C(F)(F)F)C(F)(C(F)(F)F)C(=O)O'),
        Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O'),
    ]

    for mol in mols:
        assert mol is not None
        comp = _first_component(mol)
        assert comp is not None
        assert 'size' in comp
        assert 'smarts_centrality' in comp
