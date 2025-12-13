"""
Simplified Complex Functional Groups Benchmark with PFAS-Atlas Integration

This creates a focused benchmark with reliable functional group combinations
and attempts to integrate PFAS-Atlas properly.
"""

import sys
import os
import pandas as pd
import numpy as np
import json
from datetime import datetime
from collections import defaultdict, Counter
import traceback

# Add PFASGroups to path
sys.path.insert(0, '/home/luc/git/PFASGroups')

# RDKit imports
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# PFASGroups imports
from PFASgroups.core import parse_PFAS_groups
from PFASgroups.generate_mol import (
    generate_random_mol, 
    generate_random_carbon_chain, 
    fluorinate_mol, 
    append_functional_group,
    get_attachment
)

# Try to import PFAS-Atlas in multiple ways
PFAS_ATLAS_AVAILABLE = False
pfas_atlas_predict = None

# Method 1: Direct import from benchmark directory
try:
    from pfas_atlas_direct import predict_PFAS_class_safe
    pfas_atlas_predict = predict_PFAS_class_safe
    PFAS_ATLAS_AVAILABLE = True
    print("✅ PFAS-Atlas imported from benchmark directory")
except ImportError:
    pass

# Method 2: Import from pfas-atlas repo
if not PFAS_ATLAS_AVAILABLE:
    try:
        sys.path.insert(0, '/home/luc/git/pfas-atlas')
        from pfas_atlas import predict_PFAS_class
        pfas_atlas_predict = lambda x: predict_PFAS_class(x, mhfp_error_handling=True)
        PFAS_ATLAS_AVAILABLE = True
        print("✅ PFAS-Atlas imported from repository")
    except ImportError:
        print("⚠️  PFAS-Atlas not available")

# Set random seed
np.random.seed(42)

# Simplified functional group configurations that are more likely to work
SIMPLE_COMPLEX_GROUPS = {
    'carboxylic_sulfonic': {
        'groups': [
            {"group_smiles": "C(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "S(=O)(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Carboxylic acid + Sulfonic acid',
        'expected_pfasgroups': [2, 7],
        'expected_atlas': 'Complex structure'
    },
    'aromatic_carboxylic': {
        'groups': [
            {"group_smiles": "c1ccccc1", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "C(=O)O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Aromatic + Carboxylic acid',
        'expected_pfasgroups': [51, 2],
        'expected_atlas': 'Aromatic PFASs'
    },
    'alcohol_ether': {
        'groups': [
            {"group_smiles": "O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "O", 'n': 1, 'mode': 'insert'}
        ],
        'description': 'Alcohol + Ether',
        'expected_pfasgroups': [29, 31],
        'expected_atlas': 'Complex structure'
    },
    'ketone_alcohol': {
        'groups': [
            {"group_smiles": "C(=O)", 'n': 1, 'mode': 'insert'},
            {"group_smiles": "O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Ketone + Alcohol',
        'expected_pfasgroups': [30, 29],
        'expected_atlas': 'Complex structure'
    },
    'dicarboxylic': {
        'groups': [
            {"group_smiles": "C(=O)O[H]", 'n': 2, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Dicarboxylic acid',
        'expected_pfasgroups': [2],
        'expected_atlas': 'Complex structure'
    },
    'iodide_alcohol': {
        'groups': [
            {"group_smiles": "I", 'n': 1, 'mode': 'attach', 'neighbours': ['C']},
            {"group_smiles": "O[H]", 'n': 1, 'mode': 'attach', 'neighbours': ['C']}
        ],
        'description': 'Iodide + Alcohol',
        'expected_pfasgroups': [42, 29],
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

def generate_simple_complex_molecule(functional_group_config, chain_length=8, pathtype='Polyfluoroalkyl', max_attempts=10):
    """
    Generate a molecule with functional groups using a more reliable approach
    """
    for attempt in range(max_attempts):
        try:
            # Start with a simple linear chain
            base_mol = Chem.MolFromSmiles('C' * chain_length)
            if base_mol is None:
                continue
            
            # Fluorinate the chain
            perfluorinated = (pathtype == 'Perfluoroalkyl')
            fluorination_prob = 1.0 if perfluorinated else 0.7
            mol = fluorinate_mol(base_mol, perfluorinated=perfluorinated, p=fluorination_prob)
            if mol is None:
                continue
            
            # Add functional groups one by one
            for group_config in functional_group_config['groups']:
                # Find suitable attachment points
                mol_copy = Chem.Mol(mol)
                
                n_groups = group_config.get('n', 1)
                mode = group_config.get('mode', 'attach')
                
                if mode == 'attach':
                    # Find carbon atoms suitable for attachment
                    attach_points = []
                    for atom_idx in range(mol_copy.GetNumAtoms()):
                        atom = mol_copy.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() < 4:
                            attach_points.append(atom_idx)
                    
                    # Select attachment points
                    if len(attach_points) < n_groups:
                        continue  # Not enough attachment points
                    
                    selected_points = np.random.choice(attach_points, n_groups, replace=False)
                    
                    # Add functional groups
                    for point in selected_points:
                        mol_copy = append_functional_group(
                            mol_copy,
                            group_config['group_smiles'],
                            insertion='attach',
                            m=1,
                            atom_indices=[[point, point]],  # Using same atom for both indices
                            neighbor_atoms=['F', 'H'],
                            sanitize=False
                        )
                        if mol_copy is None:
                            break
                
                elif mode == 'insert':
                    # Find bonds suitable for insertion
                    insert_bonds = []
                    for bond_idx in range(mol_copy.GetNumBonds()):
                        bond = mol_copy.GetBondWithIdx(bond_idx)
                        if (bond.GetBeginAtom().GetSymbol() == 'C' and 
                            bond.GetEndAtom().GetSymbol() == 'C'):
                            insert_bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                    
                    if len(insert_bonds) < n_groups:
                        continue  # Not enough insertion points
                    
                    # Select insertion bonds
                    selected_bonds = []
                    for _ in range(n_groups):
                        if insert_bonds:
                            bond = insert_bonds.pop(np.random.randint(len(insert_bonds)))
                            selected_bonds.append(bond)
                    
                    # Add functional groups
                    for bond in selected_bonds:
                        mol_copy = append_functional_group(
                            mol_copy,
                            group_config['group_smiles'],
                            insertion='insert',
                            m=1,
                            atom_indices=[list(bond)],
                            neighbor_atoms=['C'],
                            sanitize=False
                        )
                        if mol_copy is None:
                            break
                
                if mol_copy is None:
                    break  # Failed to add this functional group
                
                mol = mol_copy
            
            if mol is not None:
                # Final sanitization
                try:
                    Chem.SanitizeMol(mol)
                    return mol
                except:
                    continue
        
        except Exception as e:
            continue
    
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
        result = pfas_atlas_predict(smiles)
        if result is None:
            return None, None, "PFAS-Atlas returned None"
        
        first_class = result.get('first_class', 'Unknown')
        second_class = result.get('second_class', 'Unknown')
        return first_class, second_class, None
    except Exception as e:
        return None, None, str(e)

def generate_benchmark_dataset(n_molecules_per_config=20):
    """Generate benchmark dataset with simple complex functional groups"""
    
    print("🧪 Generating simplified complex functional group benchmark...")
    
    dataset = []
    groups_lookup = load_pfasgroups_definitions()
    
    for config_name, config in SIMPLE_COMPLEX_GROUPS.items():
        print(f"\n📋 Generating molecules for: {config['description']}")
        
        successful_molecules = 0
        attempts = 0
        max_attempts = n_molecules_per_config * 10
        
        existing_inchikeys = set()
        
        while successful_molecules < n_molecules_per_config and attempts < max_attempts:
            attempts += 1
            
            # Vary parameters
            chain_length = np.random.randint(6, 12)
            pathtype = np.random.choice(['Perfluoroalkyl', 'Polyfluoroalkyl'])
            
            mol = generate_simple_complex_molecule(config, chain_length, pathtype)
            
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
                    'generation_attempts': attempts
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

def analyze_performance(df):
    """Analyze performance of both systems"""
    
    print("\n📊 PERFORMANCE ANALYSIS")
    print("=" * 60)
    
    analysis = {}
    groups_lookup = load_pfasgroups_definitions()
    
    for config_name in df['config_name'].unique():
        config_df = df[df['config_name'] == config_name]
        config_info = SIMPLE_COMPLEX_GROUPS[config_name]
        
        print(f"\n🔬 {config_info['description']} ({len(config_df)} molecules)")
        
        # PFASGroups analysis
        expected_groups = config_info['expected_pfasgroups']
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
        
        pfas_success_rate = (perfect_matches + partial_matches) / len(config_df) if len(config_df) > 0 else 0
        
        print(f"  🔬 PFASGroups:")
        print(f"    Perfect: {perfect_matches}/{len(config_df)} ({perfect_matches/len(config_df)*100:.1f}%)")
        print(f"    Partial: {partial_matches}/{len(config_df)} ({partial_matches/len(config_df)*100:.1f}%)")
        print(f"    Success: {perfect_matches + partial_matches}/{len(config_df)} ({pfas_success_rate:.1%})")
        
        # PFAS-Atlas analysis
        if PFAS_ATLAS_AVAILABLE:
            expected_atlas = config_info['expected_atlas']
            atlas_correct = len(config_df[config_df['atlas_second_class'] == expected_atlas])
            atlas_valid = len(config_df[config_df['atlas_second_class'].notna()])
            atlas_success_rate = atlas_correct / atlas_valid if atlas_valid > 0 else 0
            
            print(f"  🤖 PFAS-Atlas:")
            print(f"    Correct: {atlas_correct}/{atlas_valid} ({atlas_success_rate:.1%})")
            
            # Show prediction distribution
            atlas_counts = config_df['atlas_second_class'].value_counts()
            for atlas_class, count in atlas_counts.head(3).items():
                if pd.notna(atlas_class):
                    status = "✅" if atlas_class == expected_atlas else "❌"
                    print(f"    {status} {atlas_class}: {count}")
        
        analysis[config_name] = {
            'pfas_success_rate': pfas_success_rate,
            'perfect_matches': perfect_matches,
            'partial_matches': partial_matches,
            'total_molecules': len(config_df)
        }
    
    return analysis

def create_comparison_visualizations(df):
    """Create comparison visualizations and save results"""
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save detailed results
    results_file = f'simplified_complex_benchmark_{timestamp}.csv'
    df.to_csv(results_file, index=False)
    print(f"\n📄 Detailed results saved: {results_file}")
    
    # Create summary
    summary = {
        'benchmark_info': {
            'timestamp': datetime.now().isoformat(),
            'total_molecules': len(df),
            'functional_group_configs': len(SIMPLE_COMPLEX_GROUPS),
            'pfas_atlas_available': PFAS_ATLAS_AVAILABLE,
        },
        'results_by_config': {}
    }
    
    # Analyze each configuration
    for config_name in df['config_name'].unique():
        config_df = df[df['config_name'] == config_name]
        config_info = SIMPLE_COMPLEX_GROUPS[config_name]
        
        # PFASGroups performance
        expected_groups = set(config_info['expected_pfasgroups'])
        perfect = partial = 0
        
        for _, row in config_df.iterrows():
            detected_set = set(row['pfasgroups_detected'])
            if expected_groups.issubset(detected_set):
                if detected_set == expected_groups:
                    perfect += 1
                else:
                    partial += 1
        
        pfas_success = (perfect + partial) / len(config_df) if len(config_df) > 0 else 0
        
        # PFAS-Atlas performance
        atlas_success = 0
        if PFAS_ATLAS_AVAILABLE:
            expected_atlas = config_info['expected_atlas']
            correct = len(config_df[config_df['atlas_second_class'] == expected_atlas])
            valid = len(config_df[config_df['atlas_second_class'].notna()])
            atlas_success = correct / valid if valid > 0 else 0
        
        summary['results_by_config'][config_name] = {
            'description': config_info['description'],
            'expected_pfasgroups': config_info['expected_pfasgroups'],
            'expected_atlas': config_info['expected_atlas'],
            'total_molecules': len(config_df),
            'pfas_success_rate': pfas_success,
            'pfas_perfect': perfect,
            'pfas_partial': partial,
            'atlas_success_rate': atlas_success
        }
    
    # Save summary
    summary_file = f'simplified_complex_summary_{timestamp}.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"📄 Summary saved: {summary_file}")
    
    # Save latest versions
    df.to_csv('simplified_complex_benchmark_latest.csv', index=False)
    with open('simplified_complex_summary_latest.json', 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    return summary

def main():
    """Main benchmark execution"""
    
    print("🚀 SIMPLIFIED COMPLEX FUNCTIONAL GROUPS BENCHMARK")
    print("=" * 60)
    print(f"PFAS-Atlas Available: {PFAS_ATLAS_AVAILABLE}")
    print(f"Functional Group Configurations: {len(SIMPLE_COMPLEX_GROUPS)}")
    
    try:
        # Generate benchmark dataset
        df = generate_benchmark_dataset(n_molecules_per_config=20)
        
        print(f"\n✅ Generated {len(df)} molecules successfully")
        
        if len(df) == 0:
            print("❌ No molecules generated - exiting")
            return False
        
        # Analyze performance
        analysis = analyze_performance(df)
        
        # Create visualizations and save results
        summary = create_comparison_visualizations(df)
        
        # Print final summary
        print(f"\n🏆 BENCHMARK COMPLETE")
        print("=" * 60)
        print(f"Total molecules tested: {len(df)}")
        
        overall_pfas_success = sum(config['pfas_success_rate'] * config['total_molecules'] 
                                 for config in summary['results_by_config'].values()) / len(df)
        print(f"PFASGroups overall success: {overall_pfas_success:.1%}")
        
        if PFAS_ATLAS_AVAILABLE:
            overall_atlas_success = sum(config['atlas_success_rate'] * config['total_molecules'] 
                                      for config in summary['results_by_config'].values()) / len(df)
            print(f"PFAS-Atlas overall accuracy: {overall_atlas_success:.1%}")
        
        print("\n🎯 Best performing configurations:")
        sorted_configs = sorted(summary['results_by_config'].items(), 
                              key=lambda x: x[1]['pfas_success_rate'], reverse=True)
        for config_name, config_data in sorted_configs[:3]:
            print(f"  {config_data['description']}: {config_data['pfas_success_rate']:.1%}")
        
        return True
        
    except Exception as e:
        print(f"❌ Benchmark failed: {e}")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)