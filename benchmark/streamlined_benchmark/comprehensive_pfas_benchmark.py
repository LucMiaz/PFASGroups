"""
Streamlined PFAS Benchmarking System
Generate molecules for each PFAS group and test both PFASGroups and PFAS-Atlas detection
"""

import pandas as pd
import numpy as np
import json
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import random
import sys
import os
from datetime import datetime
from collections import defaultdict

# Add PFASGroups to path
sys.path.insert(0, '/home/luc/git/PFASGroups')
from PFASgroups.core import parse_PFAS_groups

# Try to import PFAS-Atlas
PFAS_ATLAS_AVAILABLE = False
try:
    sys.path.insert(0, '/home/luc/git/PFAS-atlas')
    from pfas_atlas import predict_PFAS_class
    PFAS_ATLAS_AVAILABLE = True
    print("✅ PFAS-Atlas imported successfully")
except ImportError:
    print("❌ PFAS-Atlas not available")

class PFASBenchmarkGenerator:
    """Generate molecules for specific PFAS functional groups"""
    
    def __init__(self):
        # Load PFAS groups definitions
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            self.pfas_groups = json.load(f)
        
        # Base fluorinated scaffolds
        self.base_scaffolds = {
            'short_perfluoro': 'C(F)(F)C(F)(F)C(F)(F)F',
            'medium_perfluoro': 'C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 
            'long_perfluoro': 'C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
            'polyfluoro': 'CCC(F)(F)C(F)(F)C(F)(F)F',
            'aromatic_side': 'c1ccc(C(F)(F)C(F)(F)F)cc1'
        }
        
        # Functional group modifications
        self.functional_modifications = {
            'alcohol': lambda smiles: smiles.replace('C(F)(F)F', 'C(F)(F)O', 1),
            'carboxylic_acid': lambda smiles: smiles.replace('C(F)(F)F', 'C(=O)O', 1),
            'sulfonic_acid': lambda smiles: smiles.replace('F', 'S(=O)(=O)O', 1),
            'ketone': lambda smiles: smiles.replace('C(F)(F)', 'C(=O)', 1),
            'ether': lambda smiles: smiles.replace('F', 'OC', 1),
            'ester': lambda smiles: smiles.replace('C(F)(F)F', 'C(=O)OC', 1),
            'amide': lambda smiles: smiles.replace('C(F)(F)F', 'C(=O)N', 1),
            'acyl_halide': lambda smiles: smiles.replace('C(F)(F)F', 'C(=O)Cl', 1),
            'sulfenic_acid': lambda smiles: smiles.replace('F', 'SO', 1),
            'sulfinic_acid': lambda smiles: smiles.replace('F', 'S(=O)O', 1),
            'phosphonic_acid': lambda smiles: smiles.replace('F', 'P(=O)(O)O', 1),
            'phosphinic_acid': lambda smiles: smiles.replace('F', 'P(=O)(O)', 1),
            'ethene': lambda smiles: smiles.replace('C(F)(F)C(F)(F)', 'C(F)=C(F)', 1),
            'iodide': lambda smiles: smiles.replace('F', 'I', 1),
            'sulfonamide': lambda smiles: smiles.replace('F', 'S(=O)(=O)N', 1),
            'amine': lambda smiles: smiles.replace('F', 'N', 1),
            'alkene': lambda smiles: smiles.replace('C(F)(F)C(F)(F)', 'C(F)=C', 1),
            'alkyne': lambda smiles: smiles.replace('C(F)(F)C(F)(F)', 'C#C', 1)
        }

    def generate_single_group_molecules(self, group_info, num_molecules=10):
        """Generate molecules for a specific PFAS group"""
        
        molecules = []
        group_id = group_info['id']
        group_name = group_info['name']
        base_groups = group_info.get('base_functional_groups', [])
        
        print(f"🧪 Generating molecules for Group {group_id}: {group_name}")
        
        attempts = 0
        max_attempts = num_molecules * 10
        
        while len(molecules) < num_molecules and attempts < max_attempts:
            attempts += 1
            
            try:
                # Choose base scaffold
                scaffold_name = random.choice(list(self.base_scaffolds.keys()))
                base_smiles = self.base_scaffolds[scaffold_name]
                
                # Apply functional group modifications
                current_smiles = base_smiles
                
                for func_group in base_groups:
                    if func_group in self.functional_modifications:
                        current_smiles = self.functional_modifications[func_group](current_smiles)
                
                # Validate molecule
                mol = Chem.MolFromSmiles(current_smiles)
                if mol is not None:
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                    
                    molecules.append({
                        'group_id': group_id,
                        'group_name': group_name,
                        'smiles': canonical_smiles,
                        'base_scaffold': scaffold_name,
                        'target_groups': base_groups,
                        'generation_type': 'single_group'
                    })
                    
            except Exception as e:
                continue
        
        success_rate = len(molecules) / attempts * 100 if attempts > 0 else 0
        print(f"  ✅ Generated {len(molecules)}/{num_molecules} molecules ({success_rate:.1f}% success)")
        
        return molecules

    def generate_multi_group_molecules(self, target_groups, num_molecules=10):
        """Generate molecules with multiple functional groups"""
        
        molecules = []
        print(f"🔬 Generating multi-group molecules: {', '.join([str(g) for g in target_groups])}")
        
        # Get group information
        group_infos = {g['id']: g for g in self.pfas_groups if g['id'] in target_groups}
        
        attempts = 0
        max_attempts = num_molecules * 15
        
        while len(molecules) < num_molecules and attempts < max_attempts:
            attempts += 1
            
            try:
                # Start with longer scaffold for multiple groups
                scaffold_name = random.choice(['medium_perfluoro', 'long_perfluoro'])
                current_smiles = self.base_scaffolds[scaffold_name]
                
                applied_groups = []
                
                # Apply modifications for each target group
                for group_id in target_groups:
                    if group_id not in group_infos:
                        continue
                        
                    group_info = group_infos[group_id]
                    base_groups = group_info.get('base_functional_groups', [])
                    
                    for func_group in base_groups:
                        if func_group in self.functional_modifications:
                            try:
                                new_smiles = self.functional_modifications[func_group](current_smiles)
                                # Check if modification was successful
                                if new_smiles != current_smiles:
                                    current_smiles = new_smiles
                                    applied_groups.append(group_id)
                                    break
                            except:
                                continue
                
                # Validate molecule
                mol = Chem.MolFromSmiles(current_smiles)
                if mol is not None and len(applied_groups) >= 2:
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                    
                    molecules.append({
                        'group_ids': applied_groups,
                        'target_groups': target_groups,
                        'smiles': canonical_smiles,
                        'base_scaffold': scaffold_name,
                        'generation_type': 'multi_group',
                        'num_applied_groups': len(applied_groups)
                    })
                    
            except Exception as e:
                continue
        
        success_rate = len(molecules) / attempts * 100 if attempts > 0 else 0
        print(f"  ✅ Generated {len(molecules)}/{num_molecules} multi-group molecules ({success_rate:.1f}% success)")
        
        return molecules

    def test_molecule_with_systems(self, molecule_data):
        """Test a molecule with both PFASGroups and PFAS-Atlas"""
        
        smiles = molecule_data['smiles']
        
        # Test with PFASGroups
        pfasgroups_result = {
            'detected_groups': [],
            'success': False,
            'error': None
        }
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                matches = parse_PFAS_groups(mol, formula)
                
                # Parse the complex tuple structure returned by parse_PFAS_groups
                group_ids = []
                if isinstance(matches, list) and matches:
                    for match_tuple in matches:
                        if isinstance(match_tuple, tuple) and len(match_tuple) > 0:
                            pfas_group = match_tuple[0]
                            if hasattr(pfas_group, 'id'):
                                group_ids.append(pfas_group.id)
                            elif hasattr(pfas_group, 'groupID'):
                                group_ids.append(pfas_group.groupID)
                
                pfasgroups_result['detected_groups'] = group_ids
                pfasgroups_result['success'] = len(group_ids) > 0
                pfasgroups_result['success'] = len(group_ids) > 0
                
        except Exception as e:
            pfasgroups_result['error'] = str(e)
        
        # Test with PFAS-Atlas
        atlas_result = {
            'first_class': None,
            'second_class': None,
            'success': False,
            'error': None
        }
        
        if PFAS_ATLAS_AVAILABLE:
            try:
                atlas_prediction = predict_PFAS_class(smiles, mhfp_error_handling=True)
                
                if atlas_prediction and isinstance(atlas_prediction, dict):
                    atlas_result['first_class'] = atlas_prediction.get('first_class')
                    atlas_result['second_class'] = atlas_prediction.get('second_class') 
                    atlas_result['success'] = atlas_result['first_class'] is not None
                else:
                    atlas_result['error'] = 'No prediction returned'
                    
            except Exception as e:
                atlas_result['error'] = str(e)
        else:
            atlas_result['error'] = 'PFAS-Atlas not available'
        
        return {
            'molecule_data': molecule_data,
            'pfasgroups_result': pfasgroups_result,
            'atlas_result': atlas_result
        }

def run_comprehensive_benchmark():
    """Run comprehensive PFAS benchmark"""
    
    print("🚀 COMPREHENSIVE PFAS BENCHMARK")
    print("=" * 50)
    
    generator = PFASBenchmarkGenerator()
    
    # Target groups (29-51 except 48, alkane is too generic)
    target_single_groups = [g for g in generator.pfas_groups if 29 <= g['id'] <= 51 and g['id'] != 48]
    
    all_results = []
    
    # Part 1: Single group molecules
    print(f"\n📋 PART 1: Single Group Molecules")
    print("=" * 40)
    
    for group_info in target_single_groups:
        molecules = generator.generate_single_group_molecules(group_info, num_molecules=8)
        
        for molecule in molecules:
            result = generator.test_molecule_with_systems(molecule)
            all_results.append(result)
    
    # Part 2: Multi-group molecules 
    print(f"\n📋 PART 2: Multi-Group Molecules")
    print("=" * 40)
    
    # Define interesting multi-group combinations
    multi_group_combinations = [
        [29, 33],  # alcohol + carboxylic acid
        [30, 47],  # ketone + amine  
        [31, 36],  # ether + sulfonic acid
        [32, 35],  # ester + acyl halide
        [34, 39],  # amide + phosphonic acid
        [37, 42],  # sulfenic acid + iodide
        [29, 30, 33],  # alcohol + ketone + carboxylic acid
        [31, 36, 47],  # ether + sulfonic acid + amine
        [32, 34, 39]   # ester + amide + phosphonic acid
    ]
    
    for combo in multi_group_combinations:
        molecules = generator.generate_multi_group_molecules(combo, num_molecules=6)
        
        for molecule in molecules:
            result = generator.test_molecule_with_systems(molecule)
            all_results.append(result)
    
    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f'pfas_comprehensive_benchmark_{timestamp}.json'
    
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n💾 Results saved: {output_file}")
    
    return all_results, output_file

if __name__ == "__main__":
    results, output_file = run_comprehensive_benchmark()