"""
Test cases for PFASgroups based on example generation scripts.

This module provides comprehensive test generation for both OECD-defined and generic PFAS groups,
combining functionality from generate_OECD_pfas_examples.py and generate_generic_pfas_examples.py.
"""
from rich import print
import sys
import os
import csv
import json
from rdkit import RDLogger

# Silence RDKit warnings
RDLogger.DisableLog('rdApp.*')
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import rdBase
print(rdBase.rdkitVersion)
import networkx as nx
from datetime import datetime
import time

tests_folder = os.path.dirname(os.path.abspath(__file__))
data_folder =  os.path.join(os.path.dirname(tests_folder),'data')
# Handle imports - try relative first, then absolute
try:
    from ..core import parse_PFAS_groups
    from ..generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_group,get_attachment
except ImportError:
    try:
        from .core import parse_PFAS_groups
        from .generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_group,get_attachment
    except ImportError:
        try:
            from PFASgroups.core import parse_PFAS_groups
            from PFASgroups.generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_group,get_attachment
        except ImportError as e:
            print(f"Error importing PFASgroups modules: {e}")
            raise e

# Try to import pytest, but make it optional
try:
    import pytest
    PYTEST_AVAILABLE = True
    
    # Add pytest hooks for better integration
    def pytest_sessionstart(session):
        """Called after the Session object has been created."""
        global TEST_SUMMARY_DATA
        TEST_SUMMARY_DATA['test_start_time'] = datetime.now()
        print(f"Starting PFAS Groups test session at {TEST_SUMMARY_DATA['test_start_time']}")
    
    def pytest_sessionfinish(session, exitstatus):
        """Called after whole test run finished."""
        global TEST_SUMMARY_DATA
        try:
            summary = generate_test_summary()
            # Save detailed results
            if TEST_SUMMARY_DATA['oecd_results']:
                pd.DataFrame(TEST_SUMMARY_DATA['oecd_results']).to_csv('oecd_test_results.csv', index=False)
            if TEST_SUMMARY_DATA['generic_results']:
                pd.DataFrame(TEST_SUMMARY_DATA['generic_results']).to_csv('generic_test_results.csv', index=False)
            if TEST_SUMMARY_DATA['specificity_results'] is not None:
                TEST_SUMMARY_DATA['specificity_results'].to_csv('specificity_test_results.csv', index=False)
        except Exception as e:
            print(f"Warning: Could not generate test summary: {e}")
    
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
    (16, "Perfluoropolyethers", [{"group_smiles":"O", 'n':2, 'mode':'insert'}], 'Perfluoroalkyl'),
    (17, "Hydrofluoroethers", [{"group_smiles":"O", 'n':"[1,5]", 'mode':'insert'}], 'Polyfluoroalkyl'),
    (18, "Perfluoroalkene", [{"group_smiles":"C(F)=C(F)", 'n':1, 'mode':'insert'}], 'Perfluoroalkyl'),
    (19, "Hydrofluoroolefins", [{"group_smiles":"C(F)=C([H])", 'n':1, 'mode':'insert'}], 'Polyfluoroalkyl'),
    (20, "Hydrofluorocarbons", [{"group_smiles":"C(F)[H]", 'n':1, 'mode':'insert'}], 'Polyfluoroalkyl'),
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
    (48, 'alkane', 'C(F)(F)', 'insert'),
    (49, 'alkene', 'C(F)=C(F)', 'insert'),
    (50, 'alkyne', "C#C", 'insert'),
    (51, 'Side-chain aromatics', "c1ccccc1", 'attach'),
]

IGNORE_GROUPS = []  # Test all groups including broad ones

# Global variable to store test results for summary generation
TEST_SUMMARY_DATA = {
    'oecd_results': [],
    'generic_results': [],
    'specificity_results': None,
    'test_start_time': None,
    'test_end_time': None
}

def generate_test_summary(output_file=f'{tests_folder}/results/test_summary_report.json'):
    """
    Generate a comprehensive test summary report with accuracy and specificity metrics.
    
    Args:
        output_file: Path to save the summary report
    
    Returns:
        dict: Summary statistics
    """
    global TEST_SUMMARY_DATA
    
    if TEST_SUMMARY_DATA['test_start_time'] is None:
        TEST_SUMMARY_DATA['test_start_time'] = datetime.now()
    
    TEST_SUMMARY_DATA['test_end_time'] = datetime.now()
    test_duration = (TEST_SUMMARY_DATA['test_end_time'] - TEST_SUMMARY_DATA['test_start_time']).total_seconds()
    
    summary = {
        'test_metadata': {
            'timestamp': datetime.now().isoformat(),
            'test_start_time': TEST_SUMMARY_DATA['test_start_time'].isoformat(),
            'test_end_time': TEST_SUMMARY_DATA['test_end_time'].isoformat(),
            'test_duration_seconds': test_duration,
            'python_version': sys.version,
            'rdkit_available': True,  # If we got here, RDKit is available
            'pytest_available': PYTEST_AVAILABLE
        },
        'oecd_test_results': {},
        'generic_test_results': {},
        'specificity_test_results': {},
        'overall_summary': {}
    }
    
    # Analyze OECD test results
    if TEST_SUMMARY_DATA['oecd_results']:
        oecd_df = pd.DataFrame(TEST_SUMMARY_DATA['oecd_results'])
        total_oecd = len(oecd_df)
        detected_oecd = len(oecd_df[oecd_df['detected'] == True])
        oecd_detection_rate = detected_oecd / total_oecd if total_oecd > 0 else 0
        
        # Group-wise analysis
        group_stats = oecd_df.groupby('group_ids').agg({
            'detected': ['count', 'sum', 'mean'],
            'group_name': 'first'
        }).round(3)
        group_stats.columns = ['total_tests', 'successful_detections', 'detection_rate', 'group_name']
        
        summary['oecd_test_results'] = {
            'total_tests': int(total_oecd),
            'successful_detections': int(detected_oecd),
            'overall_detection_rate': round(oecd_detection_rate, 3),
            'group_wise_results': group_stats.reset_index().to_dict('records'),
            'worst_performing_groups': group_stats.nsmallest(5, 'detection_rate').reset_index().to_dict('records'),
            'best_performing_groups': group_stats.nlargest(5, 'detection_rate').reset_index().to_dict('records')
        }
    
    # Analyze generic test results
    if TEST_SUMMARY_DATA['generic_results']:
        generic_df = pd.DataFrame(TEST_SUMMARY_DATA['generic_results'])
        total_generic = len(generic_df)
        detected_generic = len(generic_df[generic_df['detected'] == True])
        generic_detection_rate = detected_generic / total_generic if total_generic > 0 else 0
        
        # Group-wise analysis
        group_stats = generic_df.groupby('group_ids').agg({
            'detected': ['count', 'sum', 'mean'],
            'group_name': 'first'
        }).round(3)
        group_stats.columns = ['total_tests', 'successful_detections', 'detection_rate', 'group_name']
        
        summary['generic_test_results'] = {
            'total_tests': int(total_generic),
            'successful_detections': int(detected_generic),
            'overall_detection_rate': round(generic_detection_rate, 3),
            'group_wise_results': group_stats.reset_index().to_dict('records'),
            'worst_performing_groups': group_stats.nsmallest(5, 'detection_rate').reset_index().to_dict('records'),
            'best_performing_groups': group_stats.nlargest(5, 'detection_rate').reset_index().to_dict('records')
        }
    
    # Analyze specificity test results
    if TEST_SUMMARY_DATA['specificity_results'] is not None:
        specificity_df = TEST_SUMMARY_DATA['specificity_results']
        valid_tests = specificity_df[specificity_df['valid_smiles'] == True]
        
        if len(valid_tests) > 0:
            detection_rate = valid_tests['expected_group_detected'].mean()
            specificity_rate = valid_tests['is_specific'].mean()
            
            # Analyze failed cases
            failed_detection = valid_tests[valid_tests['expected_group_detected'] == False]
            failed_specificity = valid_tests[valid_tests['is_specific'] == False]
            
            # Calculate average number of detected groups (measure of specificity)
            avg_detected_groups = valid_tests['n_detected_groups'].mean()
            
            summary['specificity_test_results'] = {
                'total_tests': int(len(valid_tests)),
                'detection_rate': round(detection_rate, 3),
                'specificity_rate': round(specificity_rate, 3),
                'average_detected_groups_per_test': round(avg_detected_groups, 2),
                'failed_detection_cases': int(len(failed_detection)),
                'failed_specificity_cases': int(len(failed_specificity)),
                'specificity_details': {
                    'tests_with_exact_match': int(len(valid_tests[valid_tests['n_detected_groups'] == 1])),
                    'tests_with_multiple_matches': int(len(valid_tests[valid_tests['n_detected_groups'] > 1])),
                    'tests_with_no_matches': int(len(valid_tests[valid_tests['n_detected_groups'] == 0]))
                }
            }
    
    # Calculate overall summary
    total_tests = 0
    total_successful = 0
    
    if 'total_tests' in summary['oecd_test_results']:
        total_tests += summary['oecd_test_results']['total_tests']
        total_successful += summary['oecd_test_results']['successful_detections']
    
    if 'total_tests' in summary['generic_test_results']:
        total_tests += summary['generic_test_results']['total_tests']
        total_successful += summary['generic_test_results']['successful_detections']
    
    overall_accuracy = total_successful / total_tests if total_tests > 0 else 0
    
    summary['overall_summary'] = {
        'total_tests_run': total_tests,
        'total_successful_detections': total_successful,
        'overall_accuracy': round(overall_accuracy, 3),
        'specificity_rate': summary['specificity_test_results'].get('specificity_rate', 0) if 'specificity_rate' in summary['specificity_test_results'] else 0,
        'detection_rate': summary['specificity_test_results'].get('detection_rate', 0) if 'detection_rate' in summary['specificity_test_results'] else 0,
        'test_status': 'PASSED' if (overall_accuracy > 0.5 and 
                                   summary['specificity_test_results'].get('detection_rate', 0) > 0.8 and 
                                   summary['specificity_test_results'].get('specificity_rate', 0) > 0.6) else 'FAILED'
    }
    
    # Save summary to file
    try:
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"\nTest summary saved to: {output_file}")
    except Exception as e:
        print(f"Warning: Could not save summary to {output_file}: {e}")
    
    # Print summary to console
    print("\n" + "="*60)
    print("PFAS GROUPS TEST SUMMARY")
    print("="*60)
    print(f"Test Duration: {test_duration:.1f} seconds")
    print(f"Overall Accuracy: {summary['overall_summary']['overall_accuracy']:.1%}")
    print(f"Overall Status: {summary['overall_summary']['test_status']}")
    
    if 'total_tests' in summary['oecd_test_results']:
        print(f"\nOECD Tests: {summary['oecd_test_results']['successful_detections']}/{summary['oecd_test_results']['total_tests']} "
              f"({summary['oecd_test_results']['overall_detection_rate']:.1%})")
    
    if 'total_tests' in summary['generic_test_results']:
        print(f"Generic Tests: {summary['generic_test_results']['successful_detections']}/{summary['generic_test_results']['total_tests']} "
              f"({summary['generic_test_results']['overall_detection_rate']:.1%})")
    
    if 'detection_rate' in summary['specificity_test_results']:
        print(f"Specificity Tests:")
        print(f"  Detection Rate: {summary['specificity_test_results']['detection_rate']:.1%}")
        print(f"  Specificity Rate: {summary['specificity_test_results']['specificity_rate']:.1%}")
        print(f"  Avg Groups per Test: {summary['specificity_test_results']['average_detected_groups_per_test']:.1f}")
    
    print("="*60)
    
    return summary

class TestPFASGroups:
    """Test class for PFAS group classification using synthetic examples."""
    
    @classmethod
    def setup_class(cls):
        """Set up test class - called once before any tests in this class."""
        global TEST_SUMMARY_DATA
        TEST_SUMMARY_DATA['test_start_time'] = datetime.now()
        print(f"Starting PFAS Groups tests at {TEST_SUMMARY_DATA['test_start_time']}")
    
    @classmethod
    def teardown_class(cls):
        """Clean up after all tests - called once after all tests in this class."""
        global TEST_SUMMARY_DATA
        if PYTEST_AVAILABLE:
            # Generate summary when running with pytest
            summary = generate_test_summary()
            # Save detailed results as well
            if TEST_SUMMARY_DATA['oecd_results']:
                pd.DataFrame(TEST_SUMMARY_DATA['oecd_results']).to_csv('oecd_test_results.csv', index=False)
            if TEST_SUMMARY_DATA['generic_results']:
                pd.DataFrame(TEST_SUMMARY_DATA['generic_results']).to_csv('generic_test_results.csv', index=False)
            if TEST_SUMMARY_DATA['specificity_results'] is not None:
                TEST_SUMMARY_DATA['specificity_results'].to_csv('specificity_test_results.csv', index=False)
    
    def test_oecd_pfas_groups(self):
        """Test OECD PFAS group detection with synthetic compounds."""
        global TEST_SUMMARY_DATA
        print("Testing OECD PFAS groups...")
        test_results = []
        
        for group_id, group_name, template, pathtype in OECD_PFAS_GROUPS:  # Test first 5 groups
            existing_test = []
            for n in range(5, 15,1):  # 3 different chain lengths
                r = 0
                k = 0
                while r < 3 and k < 5:  # 3 replicates per chain length
                    k+=1
                    try:
                        mol = generate_random_mol(n, template, 
                                                perfluorinated=(pathtype == 'Perfluoroalkyl'))
                    except Exception as e:
                        pass
                        #print(f"Error testing group {group_id}, n={n}: {e}")
                    else:
                        if mol is None:
                            continue
                        inchikey = Chem.MolToInchiKey(mol)
                        if inchikey in existing_test:
                            continue
                        else:
                            existing_test.append(inchikey)
                            r+=1
                        smiles = Chem.MolToSmiles(mol)
                        
                        # Test classification
                        formula = CalcMolFormula(mol)
                        matches = parse_PFAS_groups(mol, formula)
                        
                        # Check if target group is detected
                        detected_groups = [match[0].id for match in matches if match[0].id not in IGNORE_GROUPS]
                        is_detected = group_id in detected_groups
                        
                        result = {
                            'group_ids': group_id,
                            'group_name': group_name,
                            'chain_length': n,
                            'pathtype': pathtype,
                            'smiles': smiles,
                            'detected': is_detected,
                            'all_matches': detected_groups
                        }
                        test_results.append(result)
                        
                        # Store for summary
                        TEST_SUMMARY_DATA['oecd_results'].append(result)
                        
        
        # Verify that at least some groups are detected correctly
        if test_results:
            detection_rate = sum(result['detected'] for result in test_results) / len(test_results)
            print(f"OECD detection rate: {detection_rate:.2%}")
            if PYTEST_AVAILABLE:
                assert detection_rate > 0.5, f"Detection rate too low: {detection_rate:.2%}"
            assert detection_rate > 0.5, f"Detection rate too low: {detection_rate:.2%}"
        else:
            assert False, "No valid test results generated"
    
    def test_generic_pfas_groups(self):
        """Test generic PFAS group detection with synthetic compounds."""
        global TEST_SUMMARY_DATA
        print("Testing generic PFAS groups...")
        test_results = []
        
        for group_id, group_name, template, insertion_mode in GENERIC_PFAS_GROUPS:  # Test first 5 groups
            for pathtype in ['Perfluoroalkyl', 'Polyfluoroalkyl']:
                existing_test = []
                for n in range(5, 15, 1):  # 2 different chain lengths
                    r = 1
                    k = 0
                    while r<3 and k< 6:  # 3 replicates per chain length
                        try:
                            funcgroup_template = [{"group_smiles": template, 
                                                'n': 1, 
                                                'mode': insertion_mode, 
                                                'neighbours': ['C']}]
                            
                            mol = generate_random_mol(n, funcgroup_template, 
                                                    perfluorinated=(pathtype == 'Perfluoroalkyl'))
                        except Exception as e:
                            pass
                            #print(f"Error testing generic group {group_id}, n={n}: {e}")
                        else:
                            if mol is None:
                                continue
                            inchikey = Chem.MolToInchiKey(mol)
                            if inchikey in existing_test:
                                continue
                            else:
                                existing_test.append(inchikey)
                                r+=1
                            smiles = Chem.MolToSmiles(mol)
                            
                            # Test classification
                            formula = CalcMolFormula(mol)
                            matches = parse_PFAS_groups(mol, formula)
                            
                            # Check if target group is detected
                            detected_groups = [match[0].id for match in matches if match[0].id not in IGNORE_GROUPS]
                            is_detected = group_id in detected_groups
                            
                            result = {
                                'group_ids': group_id,
                                'group_name': group_name,
                                'chain_length': n,
                                'pathtype': pathtype,
                                'smiles': smiles,
                                'detected': is_detected,
                                'all_matches': detected_groups
                            }
                            test_results.append(result)
                            
                            # Store for summary
                            TEST_SUMMARY_DATA['generic_results'].append(result)
                            
                        
        
        # Verify that at least some groups are detected correctly
        if test_results:
            detection_rate = sum(result['detected'] for result in test_results) / len(test_results)
            print(f"Generic detection rate: {detection_rate:.2%}")
            assert detection_rate > 0.3, f"Generic detection rate too low: {detection_rate:.2%}"
        else:
            assert False, "No valid test results generated"

    def run_all_tests(self):
        """Run all tests manually (for when pytest is not available)."""
        global TEST_SUMMARY_DATA
        TEST_SUMMARY_DATA['test_start_time'] = datetime.now()
        print("Running all PFAS group tests...")
        
        oecd_result = True
        generic_result = True
        specificity_result = True
        
        try:
            self.test_oecd_pfas_groups()
            print(f"OECD tests: PASSED")
        except Exception as e:
            print(f"OECD tests: FAILED - {e}")
            oecd_result = False
        
        try:
            self.test_generic_pfas_groups()
            print(f"Generic tests: PASSED")
        except Exception as e:
            print(f"Generic tests: FAILED - {e}")
            generic_result = False
        
        try:
            self.test_specificity()
            print(f"Specificity tests: PASSED")
        except Exception as e:
            print(f"Specificity tests: FAILED - {e}")
            specificity_result = False
    
        overall_result = oecd_result and generic_result and specificity_result
        
        # Generate summary for manual runs
        summary = generate_test_summary()
        
        # Save detailed results
        if TEST_SUMMARY_DATA['oecd_results']:
            pd.DataFrame(TEST_SUMMARY_DATA['oecd_results']).to_csv('oecd_test_results.csv', index=False)
            print("OECD test details saved to: oecd_test_results.csv")
        if TEST_SUMMARY_DATA['generic_results']:
            pd.DataFrame(TEST_SUMMARY_DATA['generic_results']).to_csv('generic_test_results.csv', index=False)
            print("Generic test details saved to: generic_test_results.csv")
        if TEST_SUMMARY_DATA['specificity_results'] is not None:
            TEST_SUMMARY_DATA['specificity_results'].to_csv('specificity_test_results.csv', index=False)
            print("Specificity test details saved to: specificity_test_results.csv")
        
        return overall_result

    def test_specificity(self):
        """Test PFAS group specificity using the dedicated test function."""
        global TEST_SUMMARY_DATA
        print("Testing PFAS group specificity...")
        try:
            results_df = df_test_pfas_group_specificity(verbose=False)
            
            # Store results for summary
            TEST_SUMMARY_DATA['specificity_results'] = results_df
            
            # Calculate success metrics
            valid_tests = results_df[results_df['valid_smiles'] == True]
            detection_rate = valid_tests['expected_group_detected'].mean()
            specificity_rate = valid_tests['is_specific'].mean()
            
            # Print detailed information about failed tests
            failed_detection = valid_tests[valid_tests['expected_group_detected'] == False]
            failed_specificity = valid_tests[valid_tests['is_specific'] == False]
            
            print(f"\nDetection Rate: {detection_rate:.1%} ({len(valid_tests[valid_tests['expected_group_detected'] == True])}/{len(valid_tests)})")
            print(f"Specificity Rate: {specificity_rate:.1%} ({len(valid_tests[valid_tests['is_specific'] == True])}/{len(valid_tests)})")
            
            if len(failed_detection) > 0:
                print(f"\nFailed Detection ({len(failed_detection)} cases):")
                for _, row in failed_detection.head(5).iterrows():  # Show first 5 failures
                    print(f"  Expected: {row['group_ids']}, Detected: {row['detected_groups']}")
            
            if len(failed_specificity) > 0:
                print(f"\nFailed Specificity ({len(failed_specificity)} cases):")
                for _, row in failed_specificity.head(5).iterrows():  # Show first 5 failures
                    print(f"  Expected: {row['group_ids']}, Detected: {row['detected_groups']}")
            
            # Export misclassified molecules for visualization
            try:
                from export_misclassified import export_misclassified_molecules
                print("\nExporting misclassified molecules...")
                export_misclassified_molecules()
            except Exception as export_error:
                print(f"Warning: Could not export misclassified molecules: {export_error}")
            
            # Pass if detection rate > 80% and specificity rate > 60%
            success = detection_rate > 0.8 and specificity_rate > 0.6
            
            assert success, f"Specificity test - Detection: {detection_rate:.1%}, Specificity: {specificity_rate:.1%}"
            
        except Exception as e:
            assert False, f"Specificity test failed: {e}"


def load_graph_from_json(filename):
    """Load the list-based JSON structure and convert to NetworkX directed graph with edge types."""
    with open(filename, 'r') as f:
        data = json.load(f)
    
    # Load the PFAS groups mapping to convert names to IDs
    try:
        with open(f'{data_folder}/PFAS_groups_smarts.json', 'r') as f:
            group_data = json.load(f)
        name_to_id = {group['name']: group['id'] for group in group_data}
    except FileNotFoundError:
        # Fallback: if no mapping file, create a temporary mapping
        print("Warning: PFAS_groups_smarts.json not found, using name-based graph")
        name_to_id = {}
    
    G = nx.DiGraph()
    Gper = nx.DiGraph()
    Gpoly = nx.DiGraph()
    Gboth = nx.DiGraph()
    Goneway = nx.DiGraph()

    # Process each edge in the list
    for edge_data in data:
        source_name = edge_data['source']
        target_name = edge_data['target']
        edge_type = edge_data['edge_type']
        
        # Convert names to IDs if mapping exists, otherwise use names
        source = name_to_id.get(source_name, source_name)
        target = name_to_id.get(target_name, target_name)
        
        # Add nodes to all graphs
        for graph in [G, Gper, Gpoly, Gboth, Goneway]:
            graph.add_node(source)
            graph.add_node(target)
        
        # Add edge to main graph (all edges)
        G.add_edge(source, target)
        
        # Add to specific graph based on type
        if edge_type == 'per':
            Gper.add_edge(source, target)
        elif edge_type == 'poly':
            Gpoly.add_edge(source, target)
        elif edge_type == 'both':
            Gboth.add_edge(source, target)
        elif edge_type == 'oneway':
            Goneway.add_edge(source, target)
    
    return G, Gper, Gpoly, Gboth, Goneway
    

def load_equivalent_groups_from_json(json_file=f'{tests_folder}/specificity_test_groups.json'):
    """
    Load equivalent groups from the JSON file structure (now in list format).
    Explore the directed graph to find all descendants for each source node.
    
    Args:
        json_file: Path to the JSON file containing the edge list structure
    Returns:
        List of tuples [x, y, pathtype] where x is a child group and y is a parent group
        pathtype indicates when the relationship applies:
        - 'Perfluoroalkyl': only for perfluorinated molecules (per edges)
        - 'Polyfluoroalkyl': only for polyfluorinated molecules (poly edges)
        - 'Both': for both types (both + oneway edges)
    """
    G, Gper, Gpoly, Gboth, Goneway = load_graph_from_json(json_file)
    
    equivalent_groups = []
    
    # Process the list of edges
    for source in G.nodes:
        perdescendants = nx.descendants(Gper, source)
        polydescendants = nx.descendants(Gpoly, source)
        bothdescendants = nx.descendants(Gboth, source)
        onewaydescendants = nx.descendants(Goneway, source)
        
        # Add per edges (only for perfluorinated molecules)
        for target in perdescendants:
            equivalent_groups.append([target, source, 'Perfluoroalkyl'])
        
        # Add poly edges (only for polyfluorinated molecules)
        for target in polydescendants:
            equivalent_groups.append([target, source, 'Polyfluoroalkyl'])
        
        # Add both edges (for BOTH perfluorinated and polyfluorinated)
        for target in bothdescendants:
            equivalent_groups.append([target, source, 'Perfluoroalkyl'])
            equivalent_groups.append([target, source, 'Polyfluoroalkyl'])
        
        # Add oneway edges (for BOTH perfluorinated and polyfluorinated)
        for target in onewaydescendants:
            equivalent_groups.append([target, source, 'Perfluoroalkyl'])
            equivalent_groups.append([target, source, 'Polyfluoroalkyl'])
    
    return equivalent_groups

def check_groups_are_related(group1, group2, graph):
    """
    Check if two groups are related in the directed graph (either as ancestor/descendant).
    
    Args:
        group1, group2: Group IDs to check
        graph: NetworkX directed graph
    Returns:
        bool: True if groups are connected (either direction)
    """
    if group1 == group2:
        return True
    
    # Check if both groups exist in the graph
    if group1 not in graph.nodes or group2 not in graph.nodes:
        return False
    
    # Check if group1 is ancestor of group2 or vice versa
    return (group2 in nx.descendants(graph, group1) or 
            group1 in nx.descendants(graph, group2))

def are_detected_groups_acceptable(expected_groups, detected_groups, json_file=f'{tests_folder}/specificity_test_groups.json'):
    """
    Check if detected groups are acceptable given the expected groups and graph relationships.
    
    Args:
        expected_groups: List of expected group IDs
        detected_groups: List of detected group IDs  
        json_file: Path to the JSON file containing the graph structure
    Returns:
        bool: True if all detected groups are either expected or related to expected groups
    """
    try:
        G, Gper, Gpoly, Gboth, Goneway = load_graph_from_json(json_file)
    except (FileNotFoundError, json.JSONDecodeError):
        # If graph file doesn't exist or is invalid, fall back to exact matching
        return sorted(expected_groups) == sorted(detected_groups)
    
    # All detected groups must be either:
    # 1. In the expected groups, OR
    # 2. Connected to at least one expected group in the graph
    for detected_group in detected_groups:
        if detected_group in expected_groups:
            continue  # This group is expected
        
        # Check if this detected group is related to any expected group
        is_related = False
        for expected_group in expected_groups:
            if check_groups_are_related(detected_group, expected_group, G):
                is_related = True
                break
        
        if not is_related:
            return False  # Found an unrelated detected group
    
    return True

def create_specificity_test_molecules():
    """
    Create simple, unambiguous test molecules for each PFAS group to test specificity.
    Returns a DataFrame with expected group matches.
    """
    
    # Define test molecules for each PFAS group
    # These are designed to be unambiguous and should match only their target group(s)
    # Start with defining base carbon chains (no cycles, no double/triple bonds) with lengths 5-15 and up to 3 replicates for each length (randomly occuring duplicates are discarded)
    BASE_CHAINS = list(set([Chem.MolToSmiles(generate_random_carbon_chain(i,cycle = False, alkene=False, alkyne = False)) for i in [j for j in range(5,16) for _ in range(3)]]))

    test_results = []
    F_CHAINS = [fluorinate_mol(Chem.MolFromSmiles(chain), perfluorinated=True, p=1.0) for chain in BASE_CHAINS]
    polyF_CHAINS = [fluorinate_mol(Chem.MolFromSmiles(chain), perfluorinated=False, p=0.8) for chain in BASE_CHAINS]
    for mF,mpF in zip(F_CHAINS,polyF_CHAINS):
        # prepare 
        atom_index_attach = get_attachment(mF, 8, atom_symbols=['C'], neighbors_symbols={'C':['F','H']})  # Get more attachment points
        atom_index_insert = get_attachment(mpF, 8, atom_symbols=['C'], neighbors_symbols={'C':['C']})     # Get more attachment points for polyfluorinated
        atom_index_insert_per = get_attachment(mF, 8, atom_symbols=['C'], neighbors_symbols={'C':['C']})   # Get more attachment points for perfluorinated
        
        # filter attachment points to have multiple distinct C atoms for multi-functional groups
        attach_points = []
        insert_points = []
        insert_points_per = []
        
        # Collect multiple distinct attachment points
        used_atoms = set()
        for a, b in atom_index_attach:
            if a not in used_atoms and b not in used_atoms:
                attach_points.extend([a, b])
                used_atoms.update([a, b])
                if len(attach_points) >= 8:  # Collect up to 8 points
                    break
        
        used_atoms = set()
        for c, d in atom_index_insert:
            if c not in used_atoms and d not in used_atoms:
                insert_points.extend([c, d])
                used_atoms.update([c, d])
                if len(insert_points) >= 8:  # Collect up to 8 points
                    break
        
        atoms_indices = {'attach': attach_points, 'insert': insert_points}
        for group_id, group_name, templates, pathtype in OECD_PFAS_GROUPS:
            mol = Chem.Mol(mF) if pathtype == 'Perfluoroalkyl' else Chem.Mol(mpF)
            # Add functional group(s) - FIXED to handle multiple functional groups correctly
            for idx, template in enumerate(templates):
                n_groups = template.get('n', 1)
                mode = template.get('mode', 'attach')
                
                # Prepare atom indices based on the number of functional groups needed
                # This fixes the issue where dicarboxylic acids (n=2) were only getting 1 carboxylic acid group
                if isinstance(n_groups, int) and n_groups > 1:
                    # For multiple functional groups (like dicarboxylic acids), we need multiple attachment points
                    # Use appropriate insertion points based on pathtype for insert mode
                    available_indices = atoms_indices[mode]
                    atom_indices_for_template = []
                    
                    # Take n_groups pairs of indices for attachment
                    for i in range(min(n_groups, len(available_indices) // 2)):
                        atom_indices_for_template.append(available_indices[i*2:(i*2)+2])
                elif isinstance(n_groups, str) and '[' in n_groups and ']' in n_groups:
                    # Handle range cases like "[1,3]" or "[2,4]"
                    import re
                    range_match = re.findall(r'\[(\d+),(\d+)\]', n_groups)
                    if range_match:
                        min_n, max_n = map(int, range_match[0])
                        # Use the middle value or max for testing
                        actual_n = np.random.randint(min_n,max_n,1)[0]
                        # Use appropriate insertion points based on pathtype for insert mode
                        available_indices = atoms_indices[mode]
                        atom_indices_for_template = []
                        
                        # Take actual_n pairs of indices for attachment
                        for i in range(min(actual_n, len(available_indices) // 2)):
                            atom_indices_for_template.append(available_indices[i*2:(i*2)+2])
                        
                        # Update n_groups for the append_functional_group call
                        n_groups = actual_n
                    else:
                        # Fallback for malformed range strings
                        available_indices = atoms_indices.get(mode, [])
                        if len(available_indices) >= 2:
                            atom_indices_for_template = [available_indices[:2]]
                        else:
                            atom_indices_for_template = [available_indices]
                        n_groups = 1
                else:
                    # For single functional groups or special cases
                    # Use appropriate insertion points based on pathtype for insert mode
                    available_indices = atoms_indices.get(mode, [])
                    if len(available_indices) >= 2:
                        atom_indices_for_template = [available_indices[:2]]
                    else:
                        atom_indices_for_template = [available_indices]
                
                mol = append_functional_group(mol, template['group_smiles'],
                                insertion=mode,
                                m=n_groups,
                                atom_indices=atom_indices_for_template,
                                neighbor_atoms=['F','H'] if mode=='attach' else ['C'],
                                sanitize=False)
            if mol is None:
                continue   
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.MolToInchiKey(mol)
            smiles = Chem.MolToSmiles(mol)
            # Test classification
            formula = CalcMolFormula(mol)
            test_results.append((f"{group_name}-{'per' if pathtype=='Perfluoroalkyl' else 'poly'}",group_id, smiles,inchi,formula,inchikey,pathtype))
        for group_id, group_name, group_smiles, insertion_mode in GENERIC_PFAS_GROUPS:
            for pathtype in ['Perfluoroalkyl', 'Polyfluoroalkyl']:
                m = mF if pathtype == 'Perfluoroalkyl' else mpF
                available_indices = atoms_indices.get(insertion_mode, [])
                if len(available_indices) >= 2:
                    atom_indices_to_use = [available_indices[:2]]
                else:
                    atom_indices_to_use = [available_indices]
                    
                mol = append_functional_group(m, group_smiles,
                            insertion=insertion_mode,
                            m=1,
                            atom_indices=atom_indices_to_use,
                            neighbor_atoms=['F','H'] if insertion_mode=='attach' else ['C'],
                            sanitize=False)
                if mol is None:
                    continue
                inchi = Chem.MolToInchi(mol)
                inchikey = Chem.MolToInchiKey(mol)
                # Test classification
                formula = CalcMolFormula(mol)
                smiles = Chem.MolToSmiles(mol)
                test_results.append((f"{group_name}-{'per' if pathtype=='Perfluoroalkyl' else 'poly'}",group_id, smiles, inchi,formula,inchikey,pathtype))
    # Add expected ids for equivalent groups:
    # Load equivalent groups from JSON file structure
    per = 'Perfluoroalkyl'
    poly = 'Polyfluoroalkyl'
    equivalent_groups = load_equivalent_groups_from_json(f'{tests_folder}/specificity_test_groups.json')
    
    # Process equivalent groups and add corresponding test results
    for x_id, y_id, pathtype in equivalent_groups:
        ## Process both Polyfluoroalkyl and Perfluoroalkyl relationships
        ## For perfluoroalkyl groups (like diacids), we also need to add their parent groups
        if x_id is not None and y_id is not None:
            # Find all test results with group_id x and add them as group_id y
            for origin, group_id, smiles, inchi, formula, inchikey, result_pathtype in test_results:
                # Only add the parent group if the pathtype matches
                if group_id == x_id and pathtype == result_pathtype:
                    # Keep the original molecule's pathtype, not the relationship pathtype
                    test_results.append((origin,y_id, smiles, inchi, formula, inchikey, result_pathtype))
    specificity_test_molecules = pd.DataFrame(test_results, columns=['origin','group_ids','smiles','inchi','formula','inchikey','pathtype'])
    specificity_test_molecules = specificity_test_molecules.groupby(['inchi','inchikey','formula','pathtype']).agg({
        'group_ids': lambda x: sorted(set([int(y) for y in x])),
        'smiles': 'first',
        'origin':lambda x: ", ".join([str(y) for y in set(x)])
    }).reset_index()
    return specificity_test_molecules

def df_test_pfas_group_specificity(test_molecules=None, output_file=f'{tests_folder}/results/specificity_test_results.csv', verbose=True):
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
    results = []
    for i,(inchi, inchikey, formula, pathtype, group_ids,smiles,origin) in test_molecules.iterrows():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Warning: Invalid smiles {smiles}")
                results.append({
                    'origin':origin,
                    'groups_id': group_ids,
                    'inchi':inchi,
                    'smiles': None,
                    'valid_smiles': False,
                    'expected_group_detected': False,
                    'expected_groups': [],
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
            detected_groups = [match[0].id for match in matches if match[0].id not in IGNORE_GROUPS]
            expected_group_detected = len(set(group_ids).intersection(detected_groups))==len(group_ids)
            
            # Check specificity: ideally, only the expected group should be detected
            # But some overlap is expected (e.g., generic groups matching specific OECD groups)
            # Allow overlap when groups are connected in the directed graph
            n_detected_groups = len(detected_groups)
            is_specific = are_detected_groups_acceptable(group_ids, detected_groups)
            
            results.append({
                'origin':origin,
                'group_ids': group_ids,
                'smiles': smiles,
                'inchi':inchi,
                'valid_smiles': True,
                'expected_group_detected': expected_group_detected,
                'expected_groups': [group_ids],
                'detected_groups': detected_groups,
                'n_detected_groups': n_detected_groups,
                'is_specific': is_specific,
                'error': None
            })
            
            if verbose:
                status = "✓" if expected_group_detected else "✗"
                specificity_status = "✓" if is_specific else "✗"
                specificity = f"({n_detected_groups} groups)" if n_detected_groups > 1 else ""
                if expected_group_detected is False:
                    print(f"{status} Group {group_ids}: {specificity}")
                    print(f"    Expected: {group_ids}, Detected: {detected_groups}")
                elif not is_specific:
                    print(f"{status} Detection {specificity_status} Specificity Group {group_ids}: {specificity}")
                    print(f"    Expected: {group_ids}, Detected: {detected_groups}")
                    print(f"    Some detected groups not related to expected groups")
            
        except Exception as e:
            results.append({
                'origin':origin,
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
    # problematic = valid_tests[valid_tests['expected_group_detected'] == False]
    # if len(problematic) > 0:
    #     print(f"\nProblematic groups (not detected):")
    #     for _, row in problematic.iterrows():
    #         print(f"  Group {row['group_ids']}")
    
    # low_specificity = valid_tests[valid_tests['n_detected_groups'] > 5]
    # if len(low_specificity) > 0:
    #     print(f"\nLow specificity groups (>5 matches):")
    #     for _, row in low_specificity.iterrows():
    #         print(f"  Group {row['group_ids']}: ({row['n_detected_groups']} matches)")
    
    # Save results if requested
    print(f'{output_file=}')
    if output_file:
        try:
            df.to_csv(output_file, index=False)
        except:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
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
    
    # Count co-occurrences where detected groups are NOT expected
    unexpected_detections = {}
    for _, row in valid_results.iterrows():
        expected_groups = set(row['group_ids'])
        detected_groups = set(row['detected_groups'])
        
        # Find groups that were detected but not expected
        unexpected = detected_groups - expected_groups
        
        for unexpected_group in unexpected:
            for expected_group in expected_groups:
                pair = (expected_group, unexpected_group)  # (expected, unexpected)
                unexpected_detections[pair] = unexpected_detections.get(pair, 0) + 1
    
    # Sort by frequency
    common_unexpected = sorted(unexpected_detections.items(), key=lambda x: x[1], reverse=True)
    
    # Load group names for display
    with open(f'{data_folder}/PFAS_groups_smarts.json', 'r') as f:
        import json
        group_smarts = json.load(f)
    group_smarts = {int(v['id']):v['name'] for v in group_smarts}
    
    print("Most common unexpected group detections (Expected -> Unexpected):")
    for (expected_group, unexpected_group), count in common_unexpected[:15]:
        expected_name = group_smarts.get(expected_group, f"Unknown {expected_group}")
        unexpected_name = group_smarts.get(unexpected_group, f"Unknown {unexpected_group}")
        print(f"  Expected {expected_group:2d} ({expected_name})")
        print(f"    -> Also detected {unexpected_group:2d} ({unexpected_name}) ({count} times)")
        print()
    
    # Also show the original co-occurrence analysis for comparison
    print("\nMost common group co-occurrences (any overlap):")
    overlap_counts = {}
    for _, row in valid_results.iterrows():
        expected_groups = row['group_ids']
        detected_groups = row['detected_groups']
        
        for detected_group in detected_groups:
            for expected_group in expected_groups:
                if expected_group != detected_group:
                    pair = tuple(sorted([expected_group, detected_group]))
                    overlap_counts[pair] = overlap_counts.get(pair, 0) + 1
    
    common_overlaps = sorted(overlap_counts.items(), key=lambda x: x[1], reverse=True)
    for (group1, group2), count in common_overlaps[:10]:
        name1 = group_smarts.get(group1, f"Unknown {group1}")
        name2 = group_smarts.get(group2, f"Unknown {group2}")
        print(f"  Groups {group1:2d} ({name1}) - {group2:2d} ({name2}) ({count} times)")
    
    return {
        'unexpected_detections': unexpected_detections,
        'common_unexpected': common_unexpected,
        'overlap_counts': overlap_counts,
        'common_overlaps': common_overlaps
    }

def generate_oecd_test_compounds(output_file=f'{tests_folder}/results/oecd_test_compounds.csv', n_compounds_per_group=10):
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
                    pass
                    #print(f"Error generating compound for group {group_id}, n={n}: {e}")
    
    print(f'OECD test compounds generated in {output_file}')

def generate_generic_test_compounds(output_file=f'{tests_folder}/results/generic_test_compounds.csv', n_compounds_per_group=10):
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
                        pass
                        #print(f"Error generating compound for generic group {group_id}, n={n}: {e}")
    
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
            detected_groups = [match[0].id for match in matches if match[0].id not in IGNORE_GROUPS]
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
            detected_groups = [match[0].id for match in matches if match[0].id not in IGNORE_GROUPS]
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
    assert success_rate > 0.5
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
            assert True
        else:
            assert False,"✗ Quick test FAILED - Some issues detected"
            
    except Exception as e:
        assert False,f"✗ Quick test FAILED with error: {e}"

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
            results = df_test_pfas_group_specificity(test_molecules)
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
