"""
Pytest test suite for PFAS group identification.

This module contains focused tests for PFAS group detection using known SMILES
and their expected group identifications. Uses pytest assertions for clear pass/fail.
"""

import pytest
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger

# Silence RDKit warnings
RDLogger.DisableLog('rdApp.*')

# Import PFASgroups functionality
try:
    from PFASgroups.parser import parse_groups_in_mol, parse_mol
except ImportError:
    from ..core import parse_groups_in_mol, parse_mol


class TestPFASIdentification:
    """Test class for PFAS group identification."""
    
    def test_pfoa_detection(self):
        """Test detection of PFOA (perfluorooctanoic acid) - should detect group 1."""
        smiles = "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"  # PFOA
        expected_groups = [1]  # Perfluoroalkyl carboxylic acids
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected in PFOA. Detected: {detected_groups}"
    
    def test_pfos_detection(self):
        """Test detection of PFOS (perfluorooctane sulfonic acid) - should detect group 6."""
        smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"  # PFOS
        expected_groups = [6]  # Perfluoroalkyl sulfonic acids
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected in PFOS. Detected: {detected_groups}"
    
    def test_pfhxa_detection(self):
        """Test detection of PFHxA (perfluorohexanoic acid) - should detect group 1."""
        smiles = "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"  # PFHxA (C6)
        expected_groups = [1]  # Perfluoroalkyl carboxylic acids
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected in PFHxA. Detected: {detected_groups}"
    
    def test_pfhxs_detection(self):
        """Test detection of PFHxS (perfluorohexane sulfonic acid) - should detect group 6."""
        smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"  # PFHxS (C6)
        expected_groups = [6]  # Perfluoroalkyl sulfonic acids
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected in PFHxS. Detected: {detected_groups}"
    
    def test_pfbs_detection(self):
        """Test detection of PFBS (perfluorobutane sulfonic acid) - should detect group 6."""
        smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"  # PFBS (C4)
        expected_groups = [6]  # Perfluoroalkyl sulfonic acids
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected in PFBS. Detected: {detected_groups}"
    
    def test_ftoh_detection(self):
        """Test detection of 8:2 FTOH (fluorotelomer alcohol) - should detect groups for telomer alcohol."""
        smiles = "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"  # 8:2 FTOH
        expected_groups = [15]  # n:2 fluorotelomer alcohols
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected in 8:2 FTOH. Detected: {detected_groups}"
    
    def test_simple_carboxylic_acid(self):
        """Test detection of carboxylic acid functional group in PFAS - should detect group 33."""
        smiles = "OC(=O)C(F)(F)C(F)(F)F"  # Simple perfluorinated carboxylic acid
        expected_groups = [33]  # carboxylic acid functional group
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula, include_PFAS_definitions=True)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected. Detected: {detected_groups}"
    
    def test_sulfonic_acid_functional_group(self):
        """Test detection of sulfonic acid functional group in PFAS - should detect group 36."""
        smiles = "FC(F)(F)C(F)(F)S(=O)(=O)O"  # Simple perfluorinated sulfonic acid
        expected_groups = [36]  # sulfonic acid functional group
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula, include_PFAS_definitions=True)
        
        detected_groups = [match[0].id for match in matches]
        
        for group in expected_groups:
            assert group in detected_groups, f"Expected group {group} not detected. Detected: {detected_groups}"


@pytest.mark.parametrize("smiles,expected_groups,description", [
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", [1], "PFPA (perfluoropentanoic acid)"),
    ("FC(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O", [6], "PFPrS (perfluoropropane sulfonic acid)"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", [1], "PFDA (perfluorodecanoic acid)"),
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", [15, 29, 49, 50], "6:2 FTOH"),
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F", [15, 29, 49, 50], "4:2 FTOH"),
])
def test_parametrized_pfas_compounds(smiles, expected_groups, description):
    """Parametrized test for various PFAS compounds."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES for {description}: {smiles}"
    
    mol_with_h = Chem.AddHs(mol)
    formula = CalcMolFormula(mol_with_h)
    matches = parse_groups_in_mol(mol_with_h, formula=formula)
    
    detected_groups = [match[0].id for match in matches]
    
    for group in expected_groups:
        assert group in detected_groups, f"Expected group {group} not detected in {description}. Detected: {detected_groups}"


def test_non_pfas_exclusion():
    """Test that non-PFAS compounds are correctly excluded (should not detect PFAS groups)."""
    non_pfas_smiles = [
        "CCC(=O)O",  # Regular carboxylic acid (propanoic acid)
        "CCCCCCCCC(=O)O",  # Long chain fatty acid (nonanoic acid)  
        "CCCCS(=O)(=O)O",  # Regular sulfonic acid
    ]
    
    pfas_groups = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
    
    for smiles in non_pfas_smiles:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol_with_h)
        matches = parse_groups_in_mol(mol_with_h, formula=formula)
        
        detected_groups = [match[0].id for match in matches]
        detected_pfas_groups = [g for g in detected_groups if g in pfas_groups]
        
        assert len(detected_pfas_groups) == 0, f"Non-PFAS compound {smiles} incorrectly detected PFAS groups: {detected_pfas_groups}"


class TestPFASDefinitions:
    """Test class for PFAS definitions detection."""
    
    def test_pfas_definition_1_detection(self):
        """Test detection of PFAS Definition 1 (OECD Definition)."""
        smiles = "OC(=O)CC(F)(F)CC"  # Simple perfluorinated carboxylic acid
        expected_definitions = [1]  # OECD Definition
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        result = parse_mol(mol, include_PFAS_definitions=True)
        assert 'matches' in result, "No matches found in result"
        
        # Extract PFAS definitions from matches
        detected_definitions = []
        for match in result['matches']:
            if match.get('type') == 'PFASdefinition':
                detected_definitions.append(match['id'])
        
        for definition in expected_definitions:
            assert definition in detected_definitions, f"Expected definition {definition} not detected. Detected: {detected_definitions}"
    
    def test_pfas_definition_2_detection(self):
        """Test detection of PFAS Definition 2 (EU PFAS Restriction)."""
        smiles = "OC(=O)CC(F)(F)CC"  # Simple perfluoroalkyl ester
        expected_definitions = [2]  # EU PFAS Restriction
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        result = parse_mol(mol, include_PFAS_definitions=True)
        assert 'matches' in result, "No matches found in result"
        
        # Extract PFAS definitions from matches
        detected_definitions = []
        for match in result['matches']:
            if match.get('type') == 'PFASdefinition':
                detected_definitions.append(match['id'])
        
        for definition in expected_definitions:
            assert definition in detected_definitions, f"Expected definition {definition} not detected. Detected: {detected_definitions}"
    
    def test_pfas_definition_3_detection(self):
        """Test detection of PFAS Definition 3 (OPPT 2023)."""
        smiles = "OCCC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F"  # Short-chain telomer alcohol
        expected_definitions = [3]  # OPPT 2023
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        result = parse_mol(mol, include_PFAS_definitions=True)
        assert 'matches' in result, "No matches found in result"
        
        # Extract PFAS definitions from matches
        detected_definitions = []
        for match in result['matches']:
            if match.get('type') == 'PFASdefinition':
                detected_definitions.append(match['id'])
        
        for definition in expected_definitions:
            assert definition in detected_definitions, f"Expected definition {definition} not detected. Detected: {detected_definitions}"
    def test_pfas_definition_4_detection(self):
        """Test detection of PFAS Definition 4 (UK)."""
        smiles = "OCCC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F"  # Short-chain telomer alcohol
        expected_definitions = [4]  # UK
        
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        result = parse_mol(mol, include_PFAS_definitions=True)
        assert 'matches' in result, "No matches found in result"
        
        # Extract PFAS definitions from matches
        detected_definitions = []
        for match in result['matches']:
            if match.get('type') == 'PFASdefinition':
                detected_definitions.append(match['id'])
        
        for definition in expected_definitions:
            assert definition in detected_definitions, f"Expected definition {definition} not detected. Detected: {detected_definitions}"
    def test_pfas_definition_5_detection(self):
        """Test detection of PFAS Definition 5 (PFASStructv5)."""
        smiles = "OC(F)(F)C(F)(F)C(F)(F)"  # Simple perfluorinated compound
        expected_definitions = [5]  # PFASStructv5
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Invalid SMILES: {smiles}"
        
        result = parse_mol(mol, include_PFAS_definitions=True)
        assert 'matches' in result, "No matches found in result"
        
        # Extract PFAS definitions from matches
        detected_definitions = []
        for match in result['matches']:
            if match.get('type') == 'PFASdefinition':
                detected_definitions.append(match['id'])
        
        for definition in expected_definitions:
            assert definition in detected_definitions, f"Expected definition {definition} not detected. Detected: {detected_definitions}"
    


class TestNonPFASDefinitions:
    """Test class to verify non-PFAS compounds don't trigger PFAS definitions."""
    
    def test_non_pfas_no_definitions(self):
        """Test that non-PFAS compounds do not trigger PFAS definitions."""
        non_pfas_smiles = [
            "FC(F)CC(=O)O",  # Regular carboxylic acid
            "FC(F)CCCCCCCC(=O)O",  # Long chain fatty acid
            "FC(F)CCCS(=O)(=O)O",  # Regular sulfonic acid
            "FC(F)CC(CCCCCC(=O)O)",  # Branched fatty acid
            "FC(F)CCCCCCC(=O)O",  # Palmitic acid
        ]
        
        for smiles in non_pfas_smiles:
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"Invalid SMILES: {smiles}"
            
            result = parse_mol(mol, include_PFAS_definitions=True)
            assert 'matches' in result, "No matches found in result"
            
            # Extract PFAS definitions from matches
            detected_definitions = []
            for match in result['matches']:
                if match.get('type') == 'PFASdefinition':
                    detected_definitions.append(match['id'])
            
            assert len(detected_definitions) == 0, f"Non-PFAS compound {smiles} incorrectly triggered PFAS definitions: {detected_definitions}"
    
    def test_fluorinated_non_pfas_exclusion(self):
        """Test that non-fluorinated compounds don't trigger PFAS definitions."""
        fluorinated_non_pfas = [
            "CC(C)C",  # Simple alkane
            "CCCC",  # Butane
            "CCC(C)CC",  # Branched alkane
        ]
        
        for smiles in fluorinated_non_pfas:
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"Invalid SMILES: {smiles}"
            
            result = parse_mol(mol, include_PFAS_definitions=True)
            assert 'matches' in result, "No matches found in result"
            
            # Extract PFAS definitions from matches
            detected_definitions = []
            for match in result['matches']:
                if match.get('type') == 'PFASdefinition':
                    detected_definitions.append(match['id'])
            
            # These should not trigger PFAS definitions (or very few)
            major_pfas_definitions = [1, 2, 3]  # Main PFAS definitions
            detected_major = [d for d in detected_definitions if d in major_pfas_definitions]
            
            assert len(detected_major) == 0, f"Simple fluorinated compound {smiles} incorrectly triggered major PFAS definitions: {detected_major}"


@pytest.mark.parametrize("smiles,expected_definitions,description", [
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", [1, 3], "PFOA - should trigger definitions 1 and 3"),
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O", [1, 3], "PFOS - should trigger definitions 1 and 3"),
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", [3], "6:2 FTOH - polyfluorinated telomer"),
    ("OCCC(F)(F)C(F)(F)F", [3], "Short-chain telomer alcohol"),
])
def test_parametrized_pfas_definitions(smiles, expected_definitions, description):
    """Parametrized test for PFAS definitions detection."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES for {description}: {smiles}"
    
    result = parse_mol(mol, include_PFAS_definitions=True)
    assert 'matches' in result, "No matches found in result"
    
    # Extract PFAS definitions from matches
    detected_definitions = []
    for match in result['matches']:
        if match.get('type') == 'PFASdefinition':
            detected_definitions.append(match['id'])
    
    for definition in expected_definitions:
        assert definition in detected_definitions, f"Expected definition {definition} not detected in {description}. Detected: {detected_definitions}"


if __name__ == "__main__":
    # Run tests if called directly
    pytest.main([__file__])