from pathlib import Path

import pytest
from rdkit import Chem

from PFASgroups import parse_smiles


def test_telomer_validation_sdf():
    sdf_path = Path('benchmark/data/PubChem_fluorotelomers.sdf')
    if not sdf_path.exists():
        pytest.skip("Telomer SDF not available")

    supplier = Chem.SDMolSupplier(str(sdf_path))
    mols = [mol for mol in supplier if mol is not None]
    assert mols

    smiles = Chem.MolToSmiles(mols[0])
    result = parse_smiles(smiles)
    assert result
