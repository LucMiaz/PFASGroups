"""
Test cases for PFASgroups based on example generation scripts.
"""

import pytest
from rdkit import Chem
from .core import parse_PFAS_groups
from .generate_mol import generate_random_mol

# Example test for a generic PFAS group

def test_generic_pfas_group():
    smiles = Chem.MolToSmiles(generate_random_mol(6, [{'group_smiles':'C(=O)O','n':1,'mode':'attach','neighbours':['C']}], perfluorinated=True))
    result = parse_PFAS_groups([smiles])
    assert result is not None

# More tests can be added based on generate_generic_pfas_examples.py and generate_OECD_pfas_examples.py
