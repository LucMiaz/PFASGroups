"""
Enhanced PFAS Benchmark System
Comprehensive comparison between PFASGroups and PFAS-Atlas with larger datasets
"""

import json
import random
import time
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict
import sys
import os

# Add parent directory to path to import PFASGroups
sys.path.append('/home/luc/git/PFASGroups')

try:
    from PFASgroups.core import parse_PFAS_groups
    from PFASgroups.generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_groups
    PFASGROUPS_AVAILABLE = True
except ImportError:
    print("❌ PFASGroups not available")
    PFASGROUPS_AVAILABLE = False

# Try to import PFAS-Atlas
try:
    sys.path.append('/home/luc/git/PFAS-atlas')  
    from classification_helper.classify_pfas import classify_pfas_molecule
    ATLAS_AVAILABLE = True
    print("✅ PFAS-Atlas available")
except ImportError:
    try:
        sys.path.append('/home/luc/git/PFAS-atlas/classification_helper')
        from classify_pfas import classify_pfas_molecule
        ATLAS_AVAILABLE = True
        print("✅ PFAS-Atlas available (fallback import)")
    except ImportError:
        print("❌ PFAS-Atlas not available")
        ATLAS_AVAILABLE = False

class EnhancedPFASBenchmark:
    """Enhanced PFAS benchmark with comprehensive testing"""
    
    def __init__(self):
        # Load PFAS groups definitions
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            self.pfas_groups = json.load(f)
        
        # Target groups 29-51 (excluding 48)
        self.target_groups = [g for g in range(29, 52) if g != 48]
        
        # OECD target groups 1-28
        self.oecd_target_groups = list(range(1, 29))
        
        # Functional group definitions following test_examples.py format with proper modes
        self.functional_smarts = {
            29: {'name': 'alcohol', 'smiles': 'O[H]', 'mode': 'attach'},
            30: {'name': 'ketone', 'smiles': 'C(=O)', 'mode': 'insert'},
            31: {'name': 'ether', 'smiles': 'O', 'mode': 'insert'},
            32: {'name': 'ester', 'smiles': 'C(=O)OC', 'mode': 'insert'},
            33: {'name': 'carboxylic acid', 'smiles': 'C(=O)O', 'mode': 'attach'},
            34: {'name': 'amide', 'smiles': 'C(=O)N', 'mode': 'insert'},
            35: {'name': 'acyl halide', 'smiles': 'C(=O)Cl', 'mode': 'attach'},
            36: {'name': 'sulfonic acid', 'smiles': 'S(=O)(=O)O', 'mode': 'attach'},
            37: {'name': 'sulfenic acid', 'smiles': 'SO[H]', 'mode': 'attach'},
            38: {'name': 'sulfinic acid', 'smiles': 'S(=O)O[H]', 'mode': 'attach'},
            39: {'name': 'phosphonic acid', 'smiles': 'P(=O)(O)O', 'mode': 'attach'},
            40: {'name': 'phosphinic acid', 'smiles': 'P(=O)O', 'mode': 'attach'},
            41: {'name': 'ethene', 'smiles': 'C(F)=C(F)', 'mode': 'insert'},
            42: {'name': 'iodide', 'smiles': 'C(F)I', 'mode': 'insert'},
            43: {'name': 'sulfonamide', 'smiles': 'S(=O)(=O)N', 'mode': 'insert'},
            44: {'name': 'azole', 'smiles': 'c1ncc[nH]1', 'mode': 'attach'},
            45: {'name': 'azine', 'smiles': 'c1ncccc1', 'mode': 'attach'},
            46: {'name': 'benzodioxole', 'smiles': 'c1ccc2OC(F)(F)Oc2c1', 'mode': 'attach'},
            47: {'name': 'amine', 'smiles': 'N', 'mode': 'insert'},
            49: {'name': 'alkene', 'smiles': 'C(F)=C(F)', 'mode': 'insert'},
            50: {'name': 'alkyne', 'smiles': 'C#C', 'mode': 'insert'},
            51: {'name': 'Side-chain aromatics', 'smiles': 'c1ccccc1', 'mode': 'attach'}
        }
        
    def generate_single_group_molecules(self, group_id, count=15):
        """Generate molecules with a single functional group using proper generate_mol functions"""
        
        group_info = self.functional_smarts[group_id]
        molecules = []
        
        for i in range(count):
            try:
                # Generate base carbon chain with varying lengths
                chain_length = 3 + (i % 5)  # Chain lengths 3-7
                
                # Create functional group specification following test_examples.py format
                functional_group_spec = {
                    'group_smiles': group_info['smiles'],
                    'n': 1,
                    'mode': group_info['mode']  # Use the correct mode from functional_smarts
                }
                
                # Generate molecule using proper generate_mol function
                mol = generate_random_mol(
                    n=chain_length,
                    functional_groups=functional_group_spec,
                    perfluorinated=True,
                    cycle=(i % 7 == 0),  # Some cyclic molecules
                    alkene=(i % 5 == 0),  # Some alkenes
                    alkyne=(i % 6 == 0)   # Some alkynes
                )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    
                    molecule_data = {
                        'group_id': group_id,
                        'group_name': group_info['name'],
                        'smiles': smiles,
                        'chain_length': chain_length,
                        'target_groups': [group_id],
                        'generation_type': 'single_group'
                    }
                    molecules.append(molecule_data)
                    
            except Exception as e:
                print(f"Warning: Failed to generate molecule for group {group_id}: {e}")
                continue
        
        return molecules
    
    def generate_multi_group_molecules(self, group_ids, count=10):
        """Generate molecules with multiple functional groups using proper generate_mol functions"""
        
        molecules = []
        
        for i in range(count):
            try:
                # Generate base carbon chain with varying lengths for multi-group molecules
                chain_length = 5 + (i % 6)  # Chain lengths 5-10 for more space for multiple groups
                
                # Create functional group specifications for all target groups following test_examples.py format
                functional_groups_specs = []
                for group_id in group_ids:
                    group_info = self.functional_smarts[group_id]
                    functional_groups_specs.append({
                        'group_smiles': group_info['smiles'],
                        'n': 1,
                        'mode': group_info['mode']  # Use the correct mode for each group
                    })
                
                # Generate molecule using proper generate_mol function with multiple groups
                mol = generate_random_mol(
                    n=chain_length,
                    functional_groups=functional_groups_specs,
                    perfluorinated=True,
                    # cycle=(i % 8 == 0),  # Some cyclic molecules
                    # alkene=(i % 6 == 0),  # Some alkenes
                    # alkyne=(i % 7 == 0)   # Some alkynes
                )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    
                    molecule_data = {
                        'target_groups': group_ids,
                        'smiles': smiles,
                        'chain_length': chain_length,
                        'generation_type': 'multi_group'
                    }
                    molecules.append(molecule_data)
                    
            except Exception as e:
                print(f"Warning: Failed to generate multi-group molecule for groups {group_ids}: {e}")
                continue
        
        return molecules
    
    def test_with_pfasgroups(self, smiles):
        """Test molecule with PFASGroups detection"""
        
        start_time = time.perf_counter()
        pfasgroups_result = {
            'detected_groups': [],
            'success': False,
            'error': None,
            'execution_time': 0.0
        }
        
        if not PFASGROUPS_AVAILABLE:
            pfasgroups_result['error'] = 'PFASGroups not available'
            return pfasgroups_result
        
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
                
                pfasgroups_result['detected_groups'] = group_ids
                pfasgroups_result['success'] = len(group_ids) > 0
                
        except Exception as e:
            pfasgroups_result['error'] = str(e)
        finally:
            pfasgroups_result['execution_time'] = time.perf_counter() - start_time
        
        return pfasgroups_result
    
    def test_with_atlas(self, smiles):
        """Test molecule with PFAS-Atlas classification"""
        
        start_time = time.perf_counter()
        atlas_result = {
            'first_class': None,
            'second_class': None, 
            'success': False,
            'error': None,
            'execution_time': 0.0
        }
        
        if not ATLAS_AVAILABLE:
            atlas_result['error'] = 'PFAS-Atlas not available'
            return atlas_result
        
        try:
            # Use the correct PFAS-Atlas function classify_pfas_molecule
            predictions = classify_pfas_molecule(smiles)
            
            if predictions and len(predictions) >= 2:
                atlas_result['first_class'] = predictions[0]
                atlas_result['second_class'] = predictions[1]
                atlas_result['success'] = True
            elif predictions and len(predictions) >= 1:
                atlas_result['first_class'] = predictions[0]
                atlas_result['success'] = True
                    
        except Exception as e:
            atlas_result['error'] = str(e)
        finally:
            atlas_result['execution_time'] = time.perf_counter() - start_time
        
        return atlas_result
    
    def run_enhanced_benchmark(self, replicates = 40):
        """Run enhanced benchmark with larger datasets"""
        
        print("🚀 ENHANCED COMPREHENSIVE PFAS BENCHMARK")
        print("=" * 55)
        print(f"📊 Testing {len(self.target_groups)} functional groups with larger datasets")
        print(f"   • Single groups: {replicates} molecules each ({len(self.target_groups) * replicates} total)")
        print(f"   • Pairs: 7 combinations, 40 molecules each (280 total)")
        print(f"   • Triplets: 5 combinations, 40 molecules each (200 total)")
        
        all_results = []
        
        # Part 1: Enhanced Single group molecules (40 per group)
        print("\n📋 PART 1: Enhanced Single Group Molecules (40 each)")
        print("=" * 55)
        
        for group_id in self.target_groups:
            group_name = self.functional_smarts[group_id]['name']
            print(f"🧪 Generating molecules for Group {group_id}: {group_name}")
            
            molecules = self.generate_single_group_molecules(group_id, count=replicates)
            success_count = 0
            
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups
                    pfas_result = self.test_with_pfasgroups(mol_data['smiles'])
                    
                    # Test with PFAS-Atlas  
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'pfasgroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / replicates) * 100
            print(f"  ✅ Generated {success_count}/40 molecules ({success_rate:.1f}% success)")
        
        # Part 2: Enhanced Multi-group molecules  
        print("\n📋 PART 2: Enhanced Multi-Group Molecules")
        print("=" * 45)
        
        # 5 pairs with 10 molecules each
        multi_group_pairs = [
            [29, 33],  # alcohol + carboxylic acid
            [30, 47],  # ketone + amine
            [31, 36],  # ether + sulfonic acid  
            [32, 35],  # ester + acyl halide
            [34, 39],  # amide + phosphonic acid
            [33, 33],  # dicarboxylic acid
            [33, 36],  # carboxylic acid and sulfonic acid
        ]
        
        # 5 triplets with 10 molecules each
        multi_group_triplets = [
            [29, 30, 33],  # alcohol + ketone + carboxylic acid
            [31, 36, 47],  # ether + sulfonic acid + amine
            [32, 34, 39],  # ester + amide + phosphonic acid
            [30, 32, 47],  # ketone + ester + amine
            [29, 31, 43],  # alcohol + ether + sulfonamide
        ]
        
        # Test pairs
        print("🔬 Testing 5 functional group pairs (40 molecules each):")
        for combo in multi_group_pairs:
            combo_str = ' + '.join([self.functional_smarts[g]['name'] for g in combo])
            print(f"   • Groups {combo[0]}-{combo[1]}: {combo_str}")
            
            molecules = self.generate_multi_group_molecules(combo, count=replicates)
            success_count = 0
            
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups
                    pfas_result = self.test_with_pfasgroups(mol_data['smiles'])
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'pfasgroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / replicates) * 100
            print(f"     ✅ Generated {success_count}/40 molecules ({success_rate:.1f}% success)")
        
        # Test triplets
        print("\n🔬 Testing 5 functional group triplets (10 molecules each):")
        for combo in multi_group_triplets:
            combo_str = ' + '.join([self.functional_smarts[g]['name'] for g in combo])
            print(f"   • Groups {'-'.join(map(str, combo))}: {combo_str}")
            
            molecules = self.generate_multi_group_molecules(combo, count=replicates)
            success_count = 0
            
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups
                    pfas_result = self.test_with_pfasgroups(mol_data['smiles'])
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'pfasgroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / replicates) * 100
            print(f"     ✅ Generated {success_count}/40 molecules ({success_rate:.1f}% success)")
        
        # Save results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"pfas_enhanced_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        total_molecules = len(all_results)
        single_molecules = len(self.target_groups) * 15
        multi_molecules = total_molecules - single_molecules
        
        print(f"\n💾 Enhanced Benchmark Complete!")
        print(f"📊 Total molecules tested: {total_molecules}")
        print(f"   • Single group molecules: {single_molecules}")  
        print(f"   • Multi-group molecules: {multi_molecules}")
        print(f"💾 Results saved: {output_file}")
        
        return all_results, output_file
    
    def run_oecd_benchmark(self):
        """Run benchmark on OECD data using groups 1-28"""
        
        print("\n🚀 OECD PFAS BENCHMARK (Groups 1-28)")
        print("=" * 45)
        
        # Load OECD data
        oecd_file = '/home/luc/git/PFAS-atlas/input_data/OECD_4000/step3_OECD_Class_0812.csv'
        try:
            import pandas as pd
            oecd_data = pd.read_csv(oecd_file)
            print(f"📊 Loaded {len(oecd_data)} OECD molecules")
        except Exception as e:
            print(f"❌ Error loading OECD data: {e}")
            return [], None
            
        all_results = []
        
        # Process OECD molecules (limit to 1000 for reasonable runtime)
        max_molecules = min(1000, len(oecd_data))
        print(f"🧪 Testing {max_molecules} OECD molecules with PFASGroups (groups 1-28)")
        
        for idx, row in oecd_data.head(max_molecules).iterrows():
            smiles = row['SMILES']
            first_class = row.get('First_Class', 'Unknown')
            second_class = row.get('Second_Class', 'Unknown')
            
            if pd.isna(smiles) or smiles.strip() == '':
                continue
                
            # Test with PFASGroups
            pfas_result = self.test_with_pfasgroups(smiles)
            
            # Test with PFAS-Atlas
            atlas_result = self.test_with_atlas(smiles)
            
            molecule_data = {
                'smiles': smiles,
                'oecd_first_class': first_class,
                'oecd_second_class': second_class,
                'generation_type': 'oecd_molecule',
                'source': 'OECD_4000'
            }
            
            result = {
                'molecule_data': molecule_data,
                'pfasgroups_result': pfas_result,
                'atlas_result': atlas_result
            }
            
            all_results.append(result)
            
            if (idx + 1) % 100 == 0:
                print(f"  ✅ Processed {idx + 1}/{max_molecules} molecules")
        
        # Save OECD results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"pfas_oecd_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        print(f"\n💾 OECD Benchmark Complete!")
        print(f"📊 Total molecules tested: {len(all_results)}")
        print(f"💾 Results saved: {output_file}")
        
        return all_results, output_file

def main():
    """Main function to run enhanced benchmark"""
    
    benchmark = EnhancedPFASBenchmark()
    
    # Run both benchmarks
    print("Running Enhanced Functional Groups Benchmark...")
    enhanced_results, enhanced_file = benchmark.run_enhanced_benchmark()
    
    print("\nRunning OECD Benchmark...")
    oecd_results, oecd_file = benchmark.run_oecd_benchmark()
    
    print(f"\n🎯 Next steps:")
    print(f"   Enhanced analysis: python enhanced_analysis.py {enhanced_file} {oecd_file}")
    print(f"   Files generated:")
    print(f"     • {enhanced_file} (functional groups)")
    print(f"     • {oecd_file} (OECD data)")

if __name__ == "__main__":
    main()