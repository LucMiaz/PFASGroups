"""
Unit tests for PFASgroups core functions.
"""

from rdkit import Chem
from PFASgroups.core import parse_PFAS_groups
from PFASgroups.generate_mol import generate_random_mol

def test_parse_PFAS_groups_basic():
    """Test parse_PFAS_groups on a simple perfluorinated carboxylic acid."""
    mol = generate_random_mol(6, [{'group_smiles':'C(=O)O','n':1,'mode':'attach','neighbours':['C']}], perfluorinated=True)
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    result = parse_PFAS_groups(mol, formula)
    assert isinstance(result, list)

def test_generate_random_mol():
    """Test that generate_random_mol returns a valid RDKit molecule."""
    mol = generate_random_mol(5, [{'group_smiles':'C(=O)O','n':1,'mode':'attach','neighbours':['C']}], perfluorinated=True)
    assert Chem.MolToSmiles(mol) is not None
