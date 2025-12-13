"""
Complex Functional Groups Benchmark

This script generates molecules with 1-3 complex functional groups,
tests them with both PFAS-Atlas and PFASGroups, and benchmarks the results.

Based on the logic from test_examples.py but focused on multi-functional group compounds.
"""

import sys
import os
import pandas as pd
import numpy as np
import json
from datetime import datetime
from collections import defaultdict, Counter
import traceback
from tqdm import tqdm

# Add PFASGroups to path
sys.path.insert(0, '/home/luc/git/PFASGroups')

# RDKit imports
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# PFASGroups imports
try:
    from PFASgroups.core import parse_PFAS_groups
    from PFASgroups.generate_mol import (
        generate_random_mol, 
        generate_random_carbon_chain, 
        fluorinate_mol, 
        append_functional_group,
        get_attachment
    )
    print("✅ PFASGroups imported successfully")
except ImportError as e:
    print(f"❌ Error importing PFASGroups: {e}")
    sys.exit(1)

# PFAS-Atlas imports
try:
    sys.path.insert(0, '/home/luc/git/pfas-atlas')
    from pfas_atlas import predict_PFAS_class
    print("✅ PFAS-Atlas imported successfully")
    PFAS_ATLAS_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  PFAS-Atlas not available: {e}")
    PFAS_ATLAS_AVAILABLE = False

# Set random seed for reproducibility
np.random.seed(42)

# Define complex functional group combinations
COMPLEX_FUNCTIONAL_GROUPS = {
    'carboxylic_acid_sulfonic': {
        'groups': [
            {"group_smiles": "C(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "S(=O)(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Carboxylic acid + Sulfonic acid',
        'expected_pfasgroups': [2, 7],  # Polyfluoroalkyl carboxylic acid, sulfonic acid
        'expected_atlas': 'Complex structure'
    },
    'ether_carboxylic_ketone': {
        'groups': [
            {"group_smiles": "O", 'n': 1, 'mode': 'insert'},
            {"group_smiles": "C(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "C(=O)", 'n': 1, 'mode': 'insert'}
        ],
        'description': 'Ether + Carboxylic acid + Ketone',
        'expected_pfasgroups': [31, 2, 30],  # Ether, carboxylic acid, ketone
        'expected_atlas': 'Complex structure'
    },
    'sulfonic_ester_alcohol': {
        'groups': [
            {"group_smiles": "S(=O)(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "C(=O)OC", 'n': 1, 'mode': 'insert'},
            {"group_smiles": "O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Sulfonic acid + Ester + Alcohol',
        'expected_pfasgroups': [7, 33, 29],  # Sulfonic acid, ester, alcohol
        'expected_atlas': 'Complex structure'
    },
    'aromatic_side_chain': {
        'groups': [
            {"group_smiles": "c1ccccc1", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "C(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Aromatic ring + Carboxylic acid',
        'expected_pfasgroups': [51, 2],  # Side-chain aromatics, carboxylic acid
        'expected_atlas': 'Aromatic PFASs'
    },
    'dicarboxylic_ether': {
        'groups': [
            {"group_smiles": "C(=O)O[H]", 'n': 2, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "O", 'n': 1, 'mode': 'insert'}
        ],
        'description': 'Dicarboxylic acid + Ether',
        'expected_pfasgroups': [2, 31],  # Carboxylic acid (multiple), ether
        'expected_atlas': 'Complex structure'
    },
    'phosphonic_amide': {
        'groups': [
            {"group_smiles": "P(=O)(O[H])O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "C(=O)N", 'n': 1, 'mode': 'insert'}
        ],
        'description': 'Phosphonic acid + Amide',
        'expected_pfasgroups': [40, 34],  # Phosphonic acid, amide
        'expected_atlas': 'Complex structure'
    },
    'iodide_ketone_alcohol': {
        'groups': [
            {"group_smiles": "I", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "C(=O)", 'n': 1, 'mode': 'insert'},
            {"group_smiles": "O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Iodide + Ketone + Alcohol',
        'expected_pfasgroups': [42, 30, 29],  # Iodide, ketone, alcohol
        'expected_atlas': 'Complex structure'
    },
    'sulfonamide_ester': {
        'groups': [
            {"group_smiles": "S(=O)(=O)N", 'n': 1, 'mode': 'insert'},
            {"group_smiles": "C(=O)OC", 'n': 1, 'mode': 'insert'}
        ],
        'description': 'Sulfonamide + Ester',
        'expected_pfasgroups': [43, 33],  # Sulfonamide, ester
        'expected_atlas': 'Complex structure'
    }
}

def load_pfasgroups_definitions():
    """Load PFASGroups group definitions"""
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            groups_data = json.load(f)
        
        groups_lookup = {}
        for group in groups_data:
            groups_lookup[group['id']] = {
                'name': group['name'],
                'alias': group.get('alias', ''),
                'main_group': group.get('main_group', ''),
                'base_functional_groups': group.get('base_functional_groups', [])
            }
        
        return groups_lookup
    except Exception as e:
        print(f"Warning: Could not load PFASGroups definitions: {e}")
        return {}

def generate_complex_molecule(functional_group_config, chain_length=8, pathtype='Polyfluoroalkyl'):
    """
    Generate a molecule with multiple functional groups
    
    Args:
        functional_group_config: Configuration dict with 'groups' list
        chain_length: Length of carbon chain
        pathtype: 'Perfluoroalkyl' or 'Polyfluoroalkyl'
    
    Returns:
        RDKit Mol object or None if generation failed
    """
    try:
        # Generate base carbon chain
        base_mol = generate_random_carbon_chain(chain_length, cycle=False, alkene=False, alkyne=False)
        if base_mol is None:
            return None
        
        # Fluorinate the chain
        perfluorinated = (pathtype == 'Perfluoroalkyl')
        fluorination_prob = 1.0 if perfluorinated else 0.8
        mol = fluorinate_mol(base_mol, perfluorinated=perfluorinated, p=fluorination_prob)
        if mol is None:
            return None
        
        # Get attachment points for functional groups
        attach_points = get_attachment(mol, 10, atom_symbols=['C'], neighbors_symbols={'C': ['F', 'H']})
        insert_points = get_attachment(mol, 10, atom_symbols=['C'], neighbors_symbols={'C': ['C']})
        
        # Prepare atom indices
        attach_atoms = []
        insert_atoms = []
        
        # Collect distinct attachment points
        used_atoms = set()
        for a, b in attach_points:
            if a not in used_atoms and b not in used_atoms:
                attach_atoms.extend([a, b])
                used_atoms.update([a, b])
                if len(attach_atoms) >= 12:
                    break
        
        used_atoms = set()
        for c, d in insert_points:
            if c not in used_atoms and d not in used_atoms:
                insert_atoms.extend([c, d])
                used_atoms.update([c, d])
                if len(insert_atoms) >= 12:
                    break
        
        atoms_indices = {'attach': attach_atoms, 'insert': insert_atoms}
        
        # Add each functional group
        for idx, template in enumerate(functional_group_config['groups']):
            n_groups = template.get('n', 1)
            mode = template.get('mode', 'attach')
            
            # Get appropriate atom indices
            available_indices = atoms_indices.get(mode, [])
            if len(available_indices) < 2:
                continue  # Skip if no attachment points available
            
            # Prepare atom indices for this template
            if isinstance(n_groups, int) and n_groups > 1:
                # Multiple functional groups
                atom_indices_for_template = []
                for i in range(min(n_groups, len(available_indices) // 2)):
                    atom_indices_for_template.append(available_indices[i*2:(i*2)+2])
            else:
                # Single functional group
                atom_indices_for_template = [available_indices[:2]]
            
            # Add the functional group
            mol = append_functional_group(
                mol, 
                template['group_smiles'],
                insertion=mode,
                m=n_groups,
                atom_indices=atom_indices_for_template,
                neighbor_atoms=['F', 'H'] if mode == 'attach' else ['C'],
                sanitize=False
            )
            
            if mol is None:
                return None
        
        # Final sanitization
        try:
            Chem.SanitizeMol(mol)
            return mol
        except:
            return None
    
    except Exception as e:
        print(f"Error generating complex molecule: {e}")
        return None

def test_with_pfasgroups(mol, formula):
    """Test molecule with PFASGroups"""
    try:
        matches = parse_PFAS_groups(mol, formula)
        detected_groups = [match[0].id for match in matches]
        return detected_groups, None
    except Exception as e:
        return [], str(e)

def test_with_pfas_atlas(smiles):
    """Test molecule with PFAS-Atlas"""
    if not PFAS_ATLAS_AVAILABLE:
        return None, None, "PFAS-Atlas not available"
    
    try:
        result = predict_PFAS_class(smiles, mhfp_error_handling=True)
        if result is None:
            return None, None, "PFAS-Atlas returned None"
        
        first_class = result.get('first_class', 'Unknown')
        second_class = result.get('second_class', 'Unknown')
        return first_class, second_class, None
    except Exception as e:
        return None, None, str(e)

def generate_benchmark_dataset(n_molecules_per_config=20):
    """Generate benchmark dataset with complex functional groups"""
    
    print("🧪 Generating complex functional group benchmark dataset...")
    
    dataset = []
    groups_lookup = load_pfasgroups_definitions()
    
    total_configs = len(COMPLEX_FUNCTIONAL_GROUPS)
    
    for config_name, config in COMPLEX_FUNCTIONAL_GROUPS.items():
        print(f"\n📋 Generating molecules for: {config['description']}")
        
        successful_molecules = 0
        attempts = 0
        max_attempts = n_molecules_per_config * 5  # Allow more attempts
        
        existing_inchikeys = set()
        
        while successful_molecules < n_molecules_per_config and attempts < max_attempts:
            attempts += 1
            
            # Vary chain length and pathtype
            chain_length = np.random.randint(6, 16)
            pathtype = np.random.choice(['Perfluoroalkyl', 'Polyfluoroalkyl'])
            
            mol = generate_complex_molecule(config, chain_length, pathtype)
            
            if mol is None:
                continue
            
            try:
                # Generate identifiers
                smiles = Chem.MolToSmiles(mol)
                inchi = Chem.MolToInchi(mol)
                inchikey = Chem.MolToInchiKey(mol)
                formula = CalcMolFormula(mol)
                
                # Skip duplicates
                if inchikey in existing_inchikeys:
                    continue
                existing_inchikeys.add(inchikey)
                
                # Test with PFASGroups
                pfasgroups_detected, pfasgroups_error = test_with_pfasgroups(mol, formula)
                
                # Test with PFAS-Atlas
                atlas_first, atlas_second, atlas_error = test_with_pfas_atlas(smiles)
                
                # Record result
                result = {
                    'config_name': config_name,
                    'description': config['description'],
                    'expected_pfasgroups': config['expected_pfasgroups'],
                    'expected_atlas': config['expected_atlas'],
                    'chain_length': chain_length,
                    'pathtype': pathtype,
                    'smiles': smiles,
                    'inchi': inchi,
                    'inchikey': inchikey,
                    'formula': formula,
                    'pfasgroups_detected': pfasgroups_detected,
                    'pfasgroups_error': pfasgroups_error,
                    'atlas_first_class': atlas_first,
                    'atlas_second_class': atlas_second,
                    'atlas_error': atlas_error,
                    'molecule_index': successful_molecules + 1,
                    'generation_attempt': attempts
                }
                
                dataset.append(result)
                successful_molecules += 1
                
                if successful_molecules % 5 == 0:
                    print(f"  ✅ Generated {successful_molecules}/{n_molecules_per_config} molecules")
                
            except Exception as e:
                print(f"  ⚠️  Error processing molecule: {e}")
                continue
        
        if successful_molecules < n_molecules_per_config:
            print(f"  ⚠️  Only generated {successful_molecules}/{n_molecules_per_config} for {config_name}")
    
    return pd.DataFrame(dataset)

def analyze_pfasgroups_performance(df):
    """Analyze PFASGroups performance on complex molecules"""
    
    print("\n📊 PFASGROUPS PERFORMANCE ANALYSIS")
    print("=" * 60)
    
    analysis = {}
    groups_lookup = load_pfasgroups_definitions()
    
    for config_name in df['config_name'].unique():
        config_df = df[df['config_name'] == config_name]
        config_info = COMPLEX_FUNCTIONAL_GROUPS[config_name]
        
        print(f"\n🔬 {config_info['description']} ({len(config_df)} molecules)")
        
        # Calculate detection rates for expected groups
        expected_groups = config_info['expected_pfasgroups']
        detection_stats = []
        
        for expected_group in expected_groups:
            group_name = groups_lookup.get(expected_group, {}).get('name', f'Group {expected_group}')
            
            detected_count = 0
            for _, row in config_df.iterrows():
                if expected_group in row['pfasgroups_detected']:
                    detected_count += 1
            
            detection_rate = detected_count / len(config_df) if len(config_df) > 0 else 0
            detection_stats.append({
                'group_id': expected_group,
                'group_name': group_name,
                'detected_count': detected_count,
                'total_molecules': len(config_df),
                'detection_rate': detection_rate
            })
            
            print(f"  📋 Group {expected_group} ({group_name}): {detected_count}/{len(config_df)} ({detection_rate:.1%})")
        
        # Calculate overall performance
        perfect_matches = 0
        partial_matches = 0
        no_matches = 0
        
        for _, row in config_df.iterrows():
            detected_set = set(row['pfasgroups_detected'])
            expected_set = set(expected_groups)
            
            if expected_set.issubset(detected_set):
                if detected_set == expected_set:
                    perfect_matches += 1
                else:
                    partial_matches += 1
            else:
                no_matches += 1
        
        perfect_rate = perfect_matches / len(config_df) if len(config_df) > 0 else 0
        partial_rate = partial_matches / len(config_df) if len(config_df) > 0 else 0
        
        print(f"  🎯 Perfect matches: {perfect_matches}/{len(config_df)} ({perfect_rate:.1%})")
        print(f"  🔸 Partial matches: {partial_matches}/{len(config_df)} ({partial_rate:.1%})")
        print(f"  ❌ No matches: {no_matches}/{len(config_df)}")
        
        analysis[config_name] = {
            'config_info': config_info,
            'detection_stats': detection_stats,
            'perfect_matches': perfect_matches,
            'partial_matches': partial_matches,
            'no_matches': no_matches,
            'perfect_rate': perfect_rate,
            'partial_rate': partial_rate,
            'total_molecules': len(config_df)
        }
    
    return analysis

def analyze_pfas_atlas_performance(df):
    """Analyze PFAS-Atlas performance on complex molecules"""
    
    print("\n📊 PFAS-ATLAS PERFORMANCE ANALYSIS")
    print("=" * 60)
    
    if not PFAS_ATLAS_AVAILABLE:
        print("⚠️  PFAS-Atlas not available for analysis")
        return {}
    
    analysis = {}
    
    for config_name in df['config_name'].unique():
        config_df = df[df['config_name'] == config_name]
        config_info = COMPLEX_FUNCTIONAL_GROUPS[config_name]
        expected_atlas = config_info['expected_atlas']
        
        print(f"\n🔬 {config_info['description']} ({len(config_df)} molecules)")
        print(f"  Expected: {expected_atlas}")
        
        # Count predictions
        atlas_counts = config_df['atlas_second_class'].value_counts()
        total_valid = len(config_df[config_df['atlas_second_class'].notna()])
        
        for atlas_class, count in atlas_counts.items():
            if pd.notna(atlas_class):
                percentage = count / total_valid * 100 if total_valid > 0 else 0
                status = "✅" if atlas_class == expected_atlas else "❌"
                print(f"  {status} {atlas_class}: {count}/{total_valid} ({percentage:.1f}%)")
        
        # Calculate accuracy
        correct_predictions = len(config_df[config_df['atlas_second_class'] == expected_atlas])
        accuracy = correct_predictions / total_valid if total_valid > 0 else 0
        
        print(f"  🎯 Accuracy: {correct_predictions}/{total_valid} ({accuracy:.1%})")
        
        # Count errors
        errors = len(config_df[config_df['atlas_error'].notna()])
        if errors > 0:
            print(f"  ⚠️  Errors: {errors}")
        
        analysis[config_name] = {
            'config_info': config_info,
            'expected_atlas': expected_atlas,
            'atlas_counts': dict(atlas_counts),
            'accuracy': accuracy,
            'correct_predictions': correct_predictions,
            'total_valid': total_valid,
            'errors': errors
        }
    
    return analysis

def compare_systems_performance(df):
    """Compare PFASGroups vs PFAS-Atlas performance"""
    
    print("\n🔄 SYSTEMS COMPARISON")
    print("=" * 60)
    
    comparison = {}
    
    for config_name in df['config_name'].unique():
        config_df = df[df['config_name'] == config_name]
        config_info = COMPLEX_FUNCTIONAL_GROUPS[config_name]
        
        print(f"\n🆚 {config_info['description']}")
        
        # PFASGroups performance
        expected_groups = set(config_info['expected_pfasgroups'])
        pfas_perfect = 0
        pfas_partial = 0
        
        for _, row in config_df.iterrows():
            detected_set = set(row['pfasgroups_detected'])
            if expected_groups.issubset(detected_set):
                if detected_set == expected_groups:
                    pfas_perfect += 1
                else:
                    pfas_partial += 1
        
        pfas_success_rate = (pfas_perfect + pfas_partial) / len(config_df) if len(config_df) > 0 else 0
        
        # PFAS-Atlas performance
        atlas_correct = 0
        if PFAS_ATLAS_AVAILABLE:
            expected_atlas = config_info['expected_atlas']
            atlas_correct = len(config_df[config_df['atlas_second_class'] == expected_atlas])
            atlas_valid = len(config_df[config_df['atlas_second_class'].notna()])
            atlas_success_rate = atlas_correct / atlas_valid if atlas_valid > 0 else 0
        else:
            atlas_success_rate = 0
        
        print(f"  🔬 PFASGroups: {pfas_perfect + pfas_partial}/{len(config_df)} success ({pfas_success_rate:.1%})")
        print(f"     - Perfect: {pfas_perfect}, Partial: {pfas_partial}")
        
        if PFAS_ATLAS_AVAILABLE:
            print(f"  🤖 PFAS-Atlas: {atlas_correct}/{len(config_df)} correct ({atlas_success_rate:.1%})")
        else:
            print(f"  🤖 PFAS-Atlas: Not available")
        
        comparison[config_name] = {
            'pfas_success_rate': pfas_success_rate,
            'pfas_perfect': pfas_perfect,
            'pfas_partial': pfas_partial,
            'atlas_success_rate': atlas_success_rate,
            'atlas_correct': atlas_correct,
            'total_molecules': len(config_df)
        }
    
    return comparison

def generate_summary_report(df, pfas_analysis, atlas_analysis, comparison):
    """Generate comprehensive summary report"""
    
    summary = {
        'benchmark_info': {
            'timestamp': datetime.now().isoformat(),
            'total_molecules': len(df),
            'functional_group_configs': len(COMPLEX_FUNCTIONAL_GROUPS),
            'pfas_atlas_available': PFAS_ATLAS_AVAILABLE,
            'successful_generations': len(df[df['smiles'].notna()])
        },
        'pfasgroups_performance': {},
        'atlas_performance': {},
        'system_comparison': comparison,
        'overall_statistics': {}
    }
    
    # Aggregate PFASGroups performance
    total_perfect = sum(config['perfect_matches'] for config in pfas_analysis.values())
    total_partial = sum(config['partial_matches'] for config in pfas_analysis.values())
    total_molecules = sum(config['total_molecules'] for config in pfas_analysis.values())
    
    overall_pfas_success = (total_perfect + total_partial) / total_molecules if total_molecules > 0 else 0
    
    summary['pfasgroups_performance'] = {
        'overall_success_rate': overall_pfas_success,
        'perfect_matches': total_perfect,
        'partial_matches': total_partial,
        'total_molecules': total_molecules,
        'by_config': pfas_analysis
    }
    
    # Aggregate PFAS-Atlas performance
    if PFAS_ATLAS_AVAILABLE:
        total_atlas_correct = sum(config['correct_predictions'] for config in atlas_analysis.values())
        total_atlas_valid = sum(config['total_valid'] for config in atlas_analysis.values())
        overall_atlas_accuracy = total_atlas_correct / total_atlas_valid if total_atlas_valid > 0 else 0
        
        summary['atlas_performance'] = {
            'overall_accuracy': overall_atlas_accuracy,
            'correct_predictions': total_atlas_correct,
            'total_valid': total_atlas_valid,
            'by_config': atlas_analysis
        }
    
    # Overall statistics
    summary['overall_statistics'] = {
        'pfasgroups_success_rate': overall_pfas_success,
        'atlas_accuracy': summary.get('atlas_performance', {}).get('overall_accuracy', 0),
        'molecules_per_config': total_molecules / len(COMPLEX_FUNCTIONAL_GROUPS),
        'successful_generation_rate': len(df[df['smiles'].notna()]) / len(df) if len(df) > 0 else 0
    }
    
    return summary

def save_results(df, summary):
    """Save benchmark results and summary"""
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save detailed results
    results_file = f'complex_functional_groups_benchmark_{timestamp}.csv'
    df.to_csv(results_file, index=False)
    print(f"📄 Detailed results saved: {results_file}")
    
    # Save summary
    summary_file = f'complex_benchmark_summary_{timestamp}.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"📄 Summary saved: {summary_file}")
    
    # Save latest results (overwrite)
    df.to_csv('complex_functional_groups_benchmark_latest.csv', index=False)
    with open('complex_benchmark_summary_latest.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"📄 Latest results saved: complex_functional_groups_benchmark_latest.csv")
    
    return results_file, summary_file

def main():
    """Main benchmark execution"""
    
    print("🚀 COMPLEX FUNCTIONAL GROUPS BENCHMARK")
    print("=" * 60)
    print(f"PFAS-Atlas Available: {PFAS_ATLAS_AVAILABLE}")
    print(f"Functional Group Configurations: {len(COMPLEX_FUNCTIONAL_GROUPS)}")
    
    try:
        # Generate benchmark dataset
        df = generate_benchmark_dataset(n_molecules_per_config=15)
        
        print(f"\n✅ Generated {len(df)} molecules successfully")
        
        # Analyze performance
        pfas_analysis = analyze_pfasgroups_performance(df)
        atlas_analysis = analyze_pfas_atlas_performance(df)
        comparison = compare_systems_performance(df)
        
        # Generate summary
        summary = generate_summary_report(df, pfas_analysis, atlas_analysis, comparison)
        
        # Save results
        results_file, summary_file = save_results(df, summary)
        
        # Print final summary
        print(f"\n🏆 BENCHMARK COMPLETE")
        print("=" * 60)
        print(f"Total molecules tested: {len(df)}")
        print(f"PFASGroups overall success: {summary['overall_statistics']['pfasgroups_success_rate']:.1%}")
        
        if PFAS_ATLAS_AVAILABLE:
            print(f"PFAS-Atlas overall accuracy: {summary['overall_statistics']['atlas_accuracy']:.1%}")
        
        print(f"Results saved: {results_file}")
        print(f"Summary saved: {summary_file}")
        
        return True
        
    except Exception as e:
        print(f"❌ Benchmark failed: {e}")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)