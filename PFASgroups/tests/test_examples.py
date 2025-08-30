"""
Test cases for PFASgroups based on example generation scripts.

This module provides comprehensive test generation for both OECD-defined and generic PFAS groups,
combining functionality from generate_OECD_pfas_examples.py and generate_generic_pfas_examples.py.
"""

import sys
import os
import csv
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Handle imports - try relative first, then absolute
from core import parse_PFAS_groups
from generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_group,get_attachment


# Try to import pytest, but make it optional
try:
    import pytest
    PYTEST_AVAILABLE = True
except ImportError:
    PYTEST_AVAILABLE = False
    print("pytest not available - test classes will be defined but not runnable with pytest")

# Set random seed for reproducibility
np.random.seed(2025)

# Define OECD PFAS group templates
OECD_PFAS_GROUPS = [
    (1, "Perfluoroalkyl carboxylic acids", [{"group_smiles":"C(=O)O[H]",'n':1,'mode':'attach','neighbours':['C']}], 'Perfluoroalkyl'),
    (2, "Polyfluoroalkyl carboxylic acid", [{"group_smiles":"C(=O)O[H]",'n':1,'mode':'attach','neighbours':['C']}],  'Polyfluoroalkyl'),
    (3, "Perfluoroalkyl dicarboxylic acids", [{"group_smiles":"C(=O)O[H]",'n':2,'mode':'attach','neighbours':['C']}], 'Perfluoroalkyl'),
    (4, "Perfluoroalkylether carboxylic acids", [{"group_smiles":"C(=O)O[H]",'n':1,'mode':'attach','neighbours':['C']},
                                                 {"group_smiles":"O",'n':1,'mode':'insert'}], 'Perfluoroalkyl'),
    (5, "Polyfluoroalkylether carboxylic acid", [{"group_smiles":"C(=O)O[H]",'n':1,'mode':'attach','neighbours':['C']},
                                                 {"group_smiles":"O",'n':1,'mode':'insert'}], 'Polyfluoroalkyl'),
    (6, "Perfluoroalkyl sulfonic acids", [{"group_smiles":"S(=O)(=O)O[H]", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (7, "Polyfluoroalkyl sulfonic acid", [{"group_smiles":"S(=O)(=O)O[H]", 'n':1, 'mode':'attach'}], 'Polyfluoroalkyl'),
    (8, "Perfluoroalkyl disulfonic acids", [{"group_smiles":"S(=O)(=O)O[H]", 'n':2, 'mode':'attach'}], 'Perfluoroalkyl'),
    (9, "Perfluoroalkylether sulfonic acids", [
        {"group_smiles":"S(=O)(=O)O[H]", 'n':1, 'mode':'attach'},
        {"group_smiles":"O", 'n':1, 'mode':'insert'}
    ], 'Perfluoroalkyl'),
    (10, "Polyfluoroalkylether sulfonic acid", [
        {"group_smiles":"S(=O)(=O)O[H]", 'n':1, 'mode':'attach'},
        {"group_smiles":"O", 'n':1, 'mode':'insert'}
    ], 'Polyfluoroalkyl'),
    (11, "Perfluoroalkyl sulfinic acids", [{"group_smiles":"S(=O)O[H]", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (12, "Perfluoroalkyl phosphonic acids", [{"group_smiles":"P(=O)(O[H])O[H]", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (13, "Perfluoroalkyl phosphinic acids", [{"group_smiles":"P(=O)O[H]", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (14, "Perfluoroalkyl alcohols", [{"group_smiles":"O", 'n':"[1,3]", 'mode':'attach'}], 'Perfluoroalkyl'),
    (15, "fluorotelomer alcohols", [{"group_smiles":"CCCO[H]", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (16, "Perfluoropolyethers", [{"group_smiles":"O", 'n':"[2,4]", 'mode':'insert'}], 'Perfluoroalkyl'),
    (17, "Hydrofluoroethers", [{"group_smiles":"O", 'n':"[1,5]", 'mode':'insert'}], 'Polyfluoroalkyl'),
    (18, "Perfluoroalkene", [{"group_smiles":"C(F)=C(F)", 'n':1, 'mode':'insert'}], 'Perfluoroalkyl'),
    (19, "Hydrofluoroolefins", [{"group_smiles":"C(F)=C([H])", 'n':1, 'mode':'insert'}], 'Polyfluoroalkyl'),
    #skip hydrofluorocarbons (many tested are HFCs)
    #(20, "Hydrofluorocarbons", [{"group_smiles":"C(F)[H]", 'n':1, 'mode':'insert'}], 'Polyfluoroalkyl'),
    (21, "Semi-fluorinated alkanes", [{"group_smiles":"C(F)[H]", 'n':1, 'mode':'insert'}], 'Polyfluoroalkyl'),
    (22, "Side-chain fluorinated aromatics", [{"group_smiles":"c1ccccc1", 'n':1, 'mode':'attach'}], 'Polyfluoroalkyl'),
    (23, "Perfluoroalkane", [{"group_smiles":"C(F)(F)", 'n':1, 'mode':'insert'}], 'Perfluoroalkyl'),
    (24, "Perfluoroalkyl-tert-amines", [{"group_smiles":"N(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (25, "Perfluoroalkyl iodides", [{"group_smiles":"I", 'n':"[1,3]", 'mode':'attach'}], 'Perfluoroalkyl'),
    (26, "Perfluoroalkane sulfonyl fluorides", [{"group_smiles":"S(=O)(=O)F", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl'),
    (27, "Perfluoroalkyl ketones", [{"group_smiles":"C(=O)", 'n':1, 'mode':'insert'}], 'Perfluoroalkyl'),
    (28, "Semi-fluoroalkyl ketones", [{"group_smiles":"C(=O)", 'n':1, 'mode':'insert'}], 'Polyfluoroalkyl'),
]

# Define generic PFAS group templates
GENERIC_PFAS_GROUPS = [
    (29, 'alcohol', "O[H]", 'attach'),
    (30, 'ketone', "C(=O)", 'insert'),
    (31, 'ether', "O", 'insert'),
    (32, 'ester', "C(=O)OC", 'insert'),
    (33, 'carboxylic acid', "C(=O)O", 'attach'),
    (34, 'amide', "C(=O)N", 'insert'),
    (35, 'acyl halide', "C(=O)Cl", 'attach'),
    (36, 'sulfonic acid', "S(=O)(=O)O", 'attach'),
    (37, 'sulfenic acid', "SO[H]", 'attach'),
    (38, 'sulfinic acid', "S(=O)O[H]", 'attach'),
    (39, 'phosphonic acid', "P(=O)(O)O", 'attach'),
    (40, 'phosphinic acid', "P(=O)O", 'attach'),
    (41, 'ethene', "C(F)=C(F)", 'insert'),
    (42, 'iodide', "C(F)I", 'insert'),
    (43, 'sulfonamide', "S(=O)(=O)N", 'insert'),
    (44, 'azole', "c1ncc[nH]1", 'attach'),
    (45, 'azine', "c1ncccc1", 'attach'),
    (46, 'benzodioxole', "c1ccc2OC(F)(F)Oc2c1", 'attach'),
    (47, 'amine', "N", 'insert'),
    # skip alkane (most tested are alkanes)
    #(48, 'alkane', 'C(F)(F)', 'insert'),
    (49, 'alkene', 'C(F)=C(F)', 'insert'),
    (50, 'alkyne', "C#C", 'insert'),
    (51, 'Side-chain aromatics', "c1ccccc1", 'attach'),
]

IGNORE_GROUPS = [48,20]

class TestPFASGroups:
    """Test class for PFAS group classification using synthetic examples."""
    
    def test_oecd_pfas_groups(self):
        """Test OECD PFAS group detection with synthetic compounds."""
        print("Testing OECD PFAS groups...")
        test_results = []
        
        for group_id, group_name, template, pathtype in OECD_PFAS_GROUPS:  # Test first 5 groups
            for n in np.random.randint(5, 15, size=3):  # 3 different chain lengths
                try:
                    mol = generate_random_mol(n, template, 
                                            perfluorinated=(pathtype == 'Perfluoroalkyl'))
                    if mol is None:
                        continue
                        
                    smiles = Chem.MolToSmiles(mol)
                    
                    # Test classification
                    formula = CalcMolFormula(mol)
                    matches = parse_PFAS_groups(mol, formula)
                    
                    # Check if target group is detected
                    detected_groups = [match[0].id for match in matches if match[0].id not in IGNORE_GROUPS]
                    is_detected = group_id in detected_groups
                    
                    test_results.append({
                        'group_ids': group_id,
                        'group_name': group_name,
                        'chain_length': n,
                        'pathtype': pathtype,
                        'smiles': smiles,
                        'detected': is_detected,
                        'all_matches': detected_groups
                    })
                    
                except Exception as e:
                    print(f"Error testing group {group_id}, n={n}: {e}")
        
        # Verify that at least some groups are detected correctly
        if test_results:
            detection_rate = sum(result['detected'] for result in test_results) / len(test_results)
            print(f"OECD detection rate: {detection_rate:.2%}")
            if PYTEST_AVAILABLE:
                assert detection_rate > 0.5, f"Detection rate too low: {detection_rate:.2%}"
            return detection_rate > 0.5
        else:
            print("No valid test results generated")
            return False
    
    def test_generic_pfas_groups(self):
        """Test generic PFAS group detection with synthetic compounds."""
        print("Testing generic PFAS groups...")
        test_results = []
        
        for group_id, group_name, template, insertion_mode in GENERIC_PFAS_GROUPS:  # Test first 5 groups
            for pathtype in ['Perfluoroalkyl', 'Polyfluoroalkyl']:
                for n in np.random.randint(5, 15, size=2):  # 2 different chain lengths
                    try:
                        funcgroup_template = [{"group_smiles": template, 
                                             'n': 1, 
                                             'mode': insertion_mode, 
                                             'neighbours': ['C']}]
                        
                        mol = generate_random_mol(n, funcgroup_template, 
                                                perfluorinated=(pathtype == 'Perfluoroalkyl'))
                        if mol is None:
                            continue
                            
                        smiles = Chem.MolToSmiles(mol)
                        
                        # Test classification
                        formula = CalcMolFormula(mol)
                        matches = parse_PFAS_groups(mol, formula)
                        
                        # Check if target group is detected
                        detected_groups = [match[0].id for match in matches]
                        is_detected = group_id in detected_groups
                        
                        test_results.append({
                            'group_ids': group_id,
                            'group_name': group_name,
                            'chain_length': n,
                            'pathtype': pathtype,
                            'smiles': smiles,
                            'detected': is_detected,
                            'all_matches': detected_groups
                        })
                        
                    except Exception as e:
                        print(f"Error testing generic group {group_id}, n={n}: {e}")
        
        # Verify that at least some groups are detected correctly
        if test_results:
            detection_rate = sum(result['detected'] for result in test_results) / len(test_results)
            print(f"Generic detection rate: {detection_rate:.2%}")
            if PYTEST_AVAILABLE:
                assert detection_rate > 0.3, f"Generic detection rate too low: {detection_rate:.2%}"
            return detection_rate > 0.3
        else:
            print("No valid test results generated")
            return False

    def run_all_tests(self):
        """Run all tests manually (for when pytest is not available)."""
        print("Running all PFAS group tests...")
        
        try:
            oecd_result = self.test_oecd_pfas_groups()
            print(f"OECD tests: {'PASSED' if oecd_result else 'FAILED'}")
        except Exception as e:
            print(f"OECD tests: FAILED - {e}")
            oecd_result = False
        
        try:
            generic_result = self.test_generic_pfas_groups()
            print(f"Generic tests: {'PASSED' if generic_result else 'FAILED'}")
        except Exception as e:
            print(f"Generic tests: FAILED - {e}")
            generic_result = False
        
        try:
            specificity_result = self.test_specificity()
            print(f"Specificity tests: {'PASSED' if specificity_result else 'FAILED'}")
        except Exception as e:
            print(f"Specificity tests: FAILED - {e}")
            specificity_result = False
        
        overall_result = oecd_result and generic_result and specificity_result
        print(f"Overall test result: {'PASSED' if overall_result else 'FAILED'}")
        return overall_result

    def test_specificity(self):
        """Test PFAS group specificity using the dedicated test function."""
        print("Testing PFAS group specificity...")
        try:
            results_df = test_pfas_group_specificity(verbose=False)
            
            # Calculate success metrics
            valid_tests = results_df[results_df['valid_smiles'] == True]
            detection_rate = valid_tests['expected_group_detected'].mean()
            specificity_rate = valid_tests['is_specific'].mean()
            
            # Pass if detection rate > 80% and specificity rate > 60%
            success = detection_rate > 0.8 and specificity_rate > 0.6
            
            print(f"Specificity test - Detection: {detection_rate:.1%}, Specificity: {specificity_rate:.1%}")
            return success
            
        except Exception as e:
            print(f"Specificity test failed: {e}")
            return False

def create_specificity_test_molecules():
    """
    Create simple, unambiguous test molecules for each PFAS group to test specificity.
    Returns a DataFrame with expected group matches.
    """
    
    # Define simple test molecules for each PFAS group
    # These are designed to be unambiguous and should match only their target group(s)
    BASE_CHAINS = ['CCC','CCCC','CCCCC','CC(C)(C)C','CC(C)CC','CC(C)(C)CC','CC(C)(C)C(C)C','CCCCCC','CCC(C)(C)C(C)(C)CC']

    test_results = []
    F_CHAINS = [fluorinate_mol(Chem.MolFromSmiles(chain), perfluorinated=True, p=1.0) for chain in BASE_CHAINS]
    polyF_CHAINS = [fluorinate_mol(Chem.MolFromSmiles(chain), perfluorinated=False, p=0.8) for chain in BASE_CHAINS]
    for mF,mpF in zip(F_CHAINS,polyF_CHAINS):
        # prepare 
        atom_index_attach = get_attachment(mF,4, atom_symbols=['C'],neighbors_symbols = {'C':['F','H']})
        atom_index_insert= get_attachment(mpF,4, atom_symbols=['C'],neighbors_symbols = {'C':['C']})
        # filter attachment points to have two distinct C atoms
        a,b = atom_index_attach.pop()
        c,d = atom_index_insert.pop()
        atoms_indices = [{'attach': [a,b], 'insert':[c,d]}]
        e,g = a,c
        while e == a:
            e,f = atom_index_attach.pop()
        while g == c:
            g,h = atom_index_insert.pop()
        atoms_indices.append({'attach': [e,f], 'insert': [g,h]})
        for group_id, group_name, templates, pathtype in OECD_PFAS_GROUPS:
            mol = Chem.Mol(mF) if pathtype == 'Perfluoroalkyl' else Chem.Mol(mpF)
            # Add functional group
            for idx, template in enumerate(templates):
                mol = append_functional_group(mol,template['group_smiles'],
                                insertion = template.get('mode','attach'),
                                m = template.get('n',1),
                                atom_indices = [atoms_indices[idx].get(template.get('mode','attach'),[])],
                                neighbor_atoms=['F','H'] if template.get('mode','attach')=='attach' else ['C'],
                                sanitize = False)

            if mol is None:
                    continue   
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
            # Test classification
            formula = CalcMolFormula(mol)
            test_results.append((group_id, inchi,formula,inchikey,pathtype))
        for group_id, group_name, group_smiles, insertion_mode in GENERIC_PFAS_GROUPS:
            for pathtype in ['Perfluoroalkyl', 'Polyfluoroalkyl']:
                m = mF if pathtype == 'Perfluoroalkyl' else mpF
                mol = append_functional_group(m,group_smiles,
                            insertion = insertion_mode,
                            m = 1,
                            atom_indices = [atoms_indices[0][insertion_mode]],
                            neighbor_atoms=['F','H'] if insertion_mode=='attach' else ['C'],
                            sanitize = False)
                if mol is None:
                    continue
                inchi = Chem.MolToInchi(mol)
                inchikey = Chem.MolToInchiKey(mol)
                # Test classification
                formula = CalcMolFormula(mol)
                test_results.append((group_id, inchi,formula,inchikey,pathtype))
    # Add expected ids for equivalent groups:
    # [[x,y],...] means group x is also y
    per = 'Perfluoroalkyl'
    poly = 'Polyfluoroalkyl'
    equivalent_groups = [[41,49,[per,poly]],
                         [41,50,[per,poly]],
                         [49,50,[per,poly]],
                         [6,36,[per]],
                         [7,36,[per,poly]],
                         [1,33,[per]],
                         [25,42,[per]],
                         [13,40,[per]],
                         [22,51,[per,poly]],
                         [1,2,[per]],
                         [4,1,[per]],
                         [4,2,[per]],
                         [4,5,[poly]],
                         [5,2,[poly]],
                         [6,7,[per]],
                         [8,6,[per]],
                         [8,7,[per]],
                         [9,6,[per]],
                         [9,7,[per]],
                         [9,10,[per]],
                         [14,29,[per]],
                         [2,33,[per,poly]],
                         [7,36,[per,poly]]]]
    for x,y,types in equivalent_groups:
        for t in types:
            for group_id, inchi,formula,inchikey,pathtype in test_results:
                if group_id == x and pathtype == t:
                    test_results.append((y,inchi,formula,inchikey,t))
    specificity_test_molecules = pd.DataFrame(test_results, columns=['group_ids','inchi','formula','inchikey','pathtype'])
    specificity_test_molecules = specificity_test_molecules.groupby(['inchi','inchikey','formula','pathtype']).agg({'group_ids':lambda x: sorted(set([int(y) for y in x]))}).reset_index()
    return specificity_test_molecules

def test_pfas_group_specificity(test_molecules=None, output_file='tests/specificity_test_results.csv', verbose=True):
    """
    Test the specificity of each PFAS group by using simple, unambiguous test molecules.
    
    Args:
        output_file: Optional CSV file to save results
        verbose: Whether to print detailed results
    
    Returns:
        DataFrame with test results including expected vs actual matches
    """
    
    print("Testing PFAS group specificity...")
    if test_molecules is None:
        test_molecules = create_specificity_test_molecules()
    print(test_molecules.head())
    results = []
    for i,(inchi, inchikey, formula, pathtype, group_ids) in test_molecules.iterrows():
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                print(f"Warning: Invalid Inchi {inchi}")
                results.append({
                    'groups_id': group_ids,
                    'inchi':inchi,
                    'smiles': None,
                    'valid_smiles': False,
                    'expected_group_detected': False,
                    'detected_groups': [],
                    'n_detected_groups': 0,
                    'is_specific': False,
                    'error': 'Invalid SMILES'
                })
                continue
            smiles = Chem.MolToSmiles(mol)
            formula = CalcMolFormula(mol)
            matches = parse_PFAS_groups(mol, formula)
            # Extract detected group IDs
            # For specificity testing, we don't need to filter by pathtype since we're testing specific molecules
            # The pathtype filtering was causing issues with dual-SMARTS groups where no chains are found
            detected_groups = [match[0].id for match in matches]
            expected_group_detected = len(set(group_ids).intersection(detected_groups))==len(group_ids)
            
            # Check specificity: ideally, only the expected group should be detected
            # But some overlap is expected (e.g., generic groups matching specific OECD groups)
            n_detected_groups = len(detected_groups)
            is_specific = group_ids == detected_groups  # Allow some overlap
            
            results.append({
                'group_ids': group_ids,
                'smiles': smiles,
                'inchi':inchi,
                'valid_smiles': True,
                'expected_group_detected': expected_group_detected,
                'detected_groups': detected_groups,
                'n_detected_groups': n_detected_groups,
                'is_specific': is_specific,
                'error': None
            })
            
            if verbose:
                status = "✓" if expected_group_detected else "✗"
                specificity = f"({n_detected_groups} groups)" if n_detected_groups > 1 else ""
                print(f"{status} Group {group_ids}: {specificity}")
                if not expected_group_detected:
                    print(f"    Expected: {group_ids}, Detected: {detected_groups}")
                elif n_detected_groups > 3:
                    print(f"    Low specificity - detected {n_detected_groups} groups: {detected_groups}")
            
        except Exception as e:
            results.append({
                'group_ids': group_ids,
                'smiles': None,
                'inchi':inchi,
                'valid_smiles': True,
                'expected_group_detected': False,
                'detected_groups': [],
                'n_detected_groups': 0,
                'is_specific': False,
                'error': str(e)
            })
            
            if verbose:
                print(f"✗ Group {group_ids}: ERROR - {e}")
    
    # Convert to DataFrame for analysis
    df = pd.DataFrame(results)
    
    # Calculate summary statistics
    valid_tests = df[df['valid_smiles'] == True]
    n_valid = len(valid_tests)
    n_detected = len(valid_tests[valid_tests['expected_group_detected'] == True])
    n_specific = len(valid_tests[valid_tests['is_specific'] == True])
    
    detection_rate = n_detected / n_valid if n_valid > 0 else 0
    specificity_rate = n_specific / n_valid if n_valid > 0 else 0
    
    print(f"\nSpecificity Test Summary:")
    print(f"Valid test molecules: {n_valid}")
    print(f"Expected groups detected: {n_detected} ({detection_rate:.1%})")
    print(f"Specific detections: {n_specific} ({specificity_rate:.1%})")
    
    # Identify problematic groups
    problematic = valid_tests[valid_tests['expected_group_detected'] == False]
    if len(problematic) > 0:
        print(f"\nProblematic groups (not detected):")
        for _, row in problematic.iterrows():
            print(f"  Group {row['group_ids']}")
    
    low_specificity = valid_tests[valid_tests['n_detected_groups'] > 5]
    if len(low_specificity) > 0:
        print(f"\nLow specificity groups (>5 matches):")
        for _, row in low_specificity.iterrows():
            print(f"  Group {row['group_ids']}: ({row['n_detected_groups']} matches)")
    
    # Save results if requested
    print(f'{output_file=}')
    if output_file:
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
    
    return df

def analyze_group_overlap(specificity_results):
    """
    Analyze which groups tend to co-occur and identify potential overlap issues.
    
    Args:
        specificity_results: DataFrame from test_pfas_group_specificity()
    
    Returns:
        Dictionary with overlap analysis
    """
    
    print("\nAnalyzing group overlap patterns...")
    
    # Create a co-occurrence matrix
    valid_results = specificity_results[specificity_results['valid_smiles'] == True]
    
    # Count co-occurrences
    overlap_counts = {}
    for _, row in valid_results.iterrows():
        expected_groups = row['group_ids']
        detected_groups = row['detected_groups']
        
        for detected_group in detected_groups:
            for expected_group in expected_groups:
                if expected_group != detected_group:
                    pair = tuple(sorted([expected_group, detected_group]))
                    overlap_counts[pair] = overlap_counts.get(pair, 0) + 1
    
    # Sort by frequency
    common_overlaps = sorted(overlap_counts.items(), key=lambda x: x[1], reverse=True)
    with open('data/PFAS_groups_smarts.json', 'r') as f:
        import json
        group_smarts = json.load(f)
    group_smarts = {int(v['id']):v['name'] for v in group_smarts}
    print("Most common group overlaps:")
    for (group1, group2), count in common_overlaps[:10]:
        # Get group names
        #name1 = valid_results[valid_results['group_ids'] == group1]['group_name'].iloc[0]
        #name2 = valid_results[valid_results['group_ids'] == group2]['group_name'].iloc[0]
        print(f"  Groups {group1:2d} ({group_smarts[group1]}) - {group2:2d} ({group_smarts[group2]}) :({count} times)")
    
    return {
        'overlap_counts': overlap_counts,
        'common_overlaps': common_overlaps
    }

def generate_oecd_test_compounds(output_file='tests/oecd_test_compounds.csv', n_compounds_per_group=10):
    """Generate comprehensive test dataset for OECD PFAS groups."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['group_ids', 'group_name', 'replicate', 'chain_length', 'pathtype', 'smiles'])
        
        for group_id, group_name, template, pathtype in OECD_PFAS_GROUPS:
            for replicate in range(n_compounds_per_group):
                n = np.random.randint(5, 25)  # Random chain length 5-25
                try:
                    mol = generate_random_mol(n, template, 
                                            perfluorinated=(pathtype == 'Perfluoroalkyl'))
                    smiles = Chem.MolToSmiles(mol)
                    
                    if mol:
                        writer.writerow([group_id, group_name, replicate, n, pathtype, smiles])
                    else:
                        print(f"Invalid molecule for group {group_id}, n={n}")
                        
                except Exception as e:
                    print(f"Error generating compound for group {group_id}, n={n}: {e}")
    
    print(f'OECD test compounds generated in {output_file}')

def generate_generic_test_compounds(output_file='tests/generic_test_compounds.csv', n_compounds_per_group=10):
    """Generate comprehensive test dataset for generic PFAS groups."""
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['group_ids', 'group_name', 'replicate', 'chain_length', 'pathtype', 'smiles'])
        
        for group_id, group_name, template, insertion_mode in GENERIC_PFAS_GROUPS:
            for pathtype in ['Perfluoroalkyl', 'Polyfluoroalkyl']:
                for replicate in range(n_compounds_per_group):
                    n = np.random.randint(5, 20)  # Random chain length 5-20
                    try:
                        funcgroup_template = [{"group_smiles": template, 
                                             'n': 1, 
                                             'mode': insertion_mode, 
                                             'neighbours': ['C']}]
                        
                        mol = generate_random_mol(n, funcgroup_template, 
                                                perfluorinated=(pathtype == 'Perfluoroalkyl'))
                        smiles = Chem.MolToSmiles(mol)
                        
                        if mol:
                            writer.writerow([group_id, group_name, replicate, n, pathtype, smiles])
                        else:
                            print(f"Invalid molecule for generic group {group_id}, n={n}")
                            
                    except Exception as e:
                        print(f"Error generating compound for generic group {group_id}, n={n}: {e}")
    
    print(f'Generic test compounds generated in {output_file}')

def validate_test_compounds(input_file, output_file=None):
    """Validate generated test compounds using the PFAS classification algorithm."""
    df = pd.read_csv(input_file)
    
    print(f'Validating {len(df)} test compounds...')
    results = []
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        try:
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is None:
                results.append({'detected': False, 'error': 'Invalid SMILES'})
                continue
                
            formula = CalcMolFormula(mol)
            matches = parse_PFAS_groups(mol, formula)
            
            # Check if the expected group is detected
            detected_groups = [match[0].id for match in matches]
            is_detected = row['group_ids'] in detected_groups
            
            results.append({
                'detected': is_detected,
                'detected_groups': detected_groups,
                'n_matches': len(matches),
                'error': None
            })
            
        except Exception as e:
            results.append({'detected': False, 'error': str(e)})
    
    # Add results to dataframe
    for i, result in enumerate(results):
        for key, value in result.items():
            df.loc[i, key] = value
    
    if output_file:
        df.to_csv(output_file, index=False)
        print(f'Validation results saved to {output_file}')
    
    # Print summary statistics
    detection_rate = df['detected'].mean()
    error_rate = df['error'].notna().mean()
    
    print(f'Overall detection rate: {detection_rate:.2%}')
    print(f'Error rate: {error_rate:.2%}')
    
    # Group-wise statistics
    group_stats = df.groupby('group_name')['detected'].agg(['count', 'sum', 'mean'])
    group_stats.columns = ['total', 'detected', 'detection_rate']
    print('\nGroup-wise detection rates:')
    print(group_stats.sort_values('detection_rate', ascending=False))
    
    return df

def test_pfoa_like_compounds():
    """Test detection of PFOA-like compounds (perfluoroalkyl carboxylic acids)."""
    print("Testing PFOA-like compounds...")
    # Generate PFOA (C8) and similar compounds
    test_smiles = [
        "C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F",  # PFOA
        "C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)(F)F",  # PFHA (C6)
        "C(C(C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F",  # PFDA (C10)
    ]
    
    success_count = 0
    for i, smiles in enumerate(test_smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            Chem.AddHs(mol)
            formula = CalcMolFormula(mol)
            matches = parse_PFAS_groups(mol, formula)
            
            # Should detect group 1 (Perfluoroalkyl carboxylic acids)
            detected_groups = [match[0].id for match in matches]
            print(detected_groups)
            print(formula)
            if 1 in detected_groups:
                print(f"✓ PFCA detected in compound {i+1}")
                success_count += 1
            else:
                print(f"✗ PFCA not detected in compound {i+1}: {smiles}")
                
            if PYTEST_AVAILABLE:
                assert 1 in detected_groups, f"PFCA not detected in {smiles}"
                
        except Exception as e:
            print(f"✗ Error testing compound {i+1}: {e}")
    
    success_rate = success_count / len(test_smiles)
    print(f"PFOA-like compound detection rate: {success_rate:.1%}")
    return success_rate > 0.5

def run_quick_test():
    """Run a quick test to verify the module is working."""
    print("Running quick functionality test...")
    
    try:
        # Test basic functionality
        tester = TestPFASGroups()
        
        # Run a simplified version of the tests
        print("Testing basic PFAS group detection...")
        
        # Test PFOA detection
        pfoa_success = test_pfoa_like_compounds()
        
        if pfoa_success:
            print("✓ Quick test PASSED - Module is working correctly")
            return True
        else:
            print("✗ Quick test FAILED - Some issues detected")
            return False
            
    except Exception as e:
        print(f"✗ Quick test FAILED with error: {e}")
        return False

if __name__ == "__main__":
    print("PFASgroups Test Examples Module")
    print("=" * 40)
    
    # Check command line arguments
    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
        
        if command == "quick":
            # Run quick test
            run_quick_test()
        
        elif command == "full":
            # Run full test suite
            tester = TestPFASGroups()
            tester.run_all_tests()
        
        elif command == "generate":
            # Generate test datasets
            print("Generating test datasets...")
            n_compounds = int(sys.argv[2]) if len(sys.argv) > 2 else 5
            
            print(f"Generating OECD test compounds ({n_compounds} per group)...")
            generate_oecd_test_compounds('oecd_test_compounds.csv', n_compounds_per_group=n_compounds)
            
            print(f"Generating generic test compounds ({n_compounds} per group)...")
            generate_generic_test_compounds('generic_test_compounds.csv', n_compounds_per_group=n_compounds)
        
        elif command == "validate":
            # Validate existing test compounds
            oecd_file = sys.argv[2] if len(sys.argv) > 2 else 'oecd_test_compounds.csv'
            generic_file = sys.argv[3] if len(sys.argv) > 3 else 'generic_test_compounds.csv'
            
            if os.path.exists(oecd_file):
                print(f"Validating OECD compounds from {oecd_file}...")
                validate_test_compounds(oecd_file, 'oecd_validated.csv')
            
            if os.path.exists(generic_file):
                print(f"Validating generic compounds from {generic_file}...")
                validate_test_compounds(generic_file, 'generic_validated.csv')
        
        elif command == "specificity":
            # Run specificity tests
            print("Running PFAS group specificity tests...")
            print("This tests that each group only matches appropriate compounds\n")
            
            test_molecules = create_specificity_test_molecules()
            results = test_pfas_group_specificity(test_molecules)
            analyze_group_overlap(results)
        
        else:
            print(f"Unknown command: {command}")
            print("Available commands:")
            print("  quick       - Run quick functionality test")
            print("  full        - Run full test suite") 
            print("  specificity - Test group specificity (no false positives)")
            print("  generate [n] - Generate test datasets (n compounds per group)")
            print("  validate [oecd_file] [generic_file] - Validate test compounds")
    
    else:
        # Default: run quick test
        print("No command specified. Running quick test...")
        print("Use 'python test_examples.py --help' for more options")
        print()
        run_quick_test()
        
        print("\nAvailable commands:")
        print("  python test_examples.py quick")
        print("  python test_examples.py full") 
        print("  python test_examples.py specificity")
        print("  python test_examples.py generate [n]")
        print("  python test_examples.py validate [files]")
