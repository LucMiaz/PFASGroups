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
    
    def get_system_specifications(self):
        """Collect system specifications for reproducibility"""
        import platform
        import psutil
        import cpuinfo
        
        try:
            cpu_info = cpuinfo.get_cpu_info()
            cpu_name = cpu_info.get('brand_raw', 'Unknown CPU')
            cpu_count = psutil.cpu_count(logical=True)
            cpu_count_physical = psutil.cpu_count(logical=False)
        except:
            cpu_name = platform.processor()
            cpu_count = os.cpu_count()
            cpu_count_physical = os.cpu_count()
        
        memory = psutil.virtual_memory()
        
        specs = {
            'system': platform.system(),
            'platform': platform.platform(),
            'architecture': platform.architecture()[0],
            'python_version': platform.python_version(),
            'cpu_name': cpu_name,
            'cpu_cores_logical': cpu_count,
            'cpu_cores_physical': cpu_count_physical,
            'total_memory_gb': round(memory.total / (1024**3), 2),
            'available_memory_gb': round(memory.available / (1024**3), 2)
        }
        
        return specs
    
    def run_timing_benchmark(self, max_molecules=200, iterations=10):
        """Run timing benchmark comparing PFASGroups vs PFAS-Atlas on increasingly large molecules
        
        Args:
            max_molecules: Number of molecules to test
            iterations: Number of iterations for averaging (default 10)
        """
        
        # Get system specifications
        system_specs = self.get_system_specifications()
        
        print("\n🚀 TIMING PERFORMANCE BENCHMARK")
        print("=" * 55)
        print(f"📊 Testing {max_molecules} molecules with carboxylic acid functional groups")
        print(f"🔄 Running {iterations} iterations for statistical averaging")
        print("📏 Molecules will vary in size to test scaling performance (up to ~2000 atoms)")
        print(f"\n💻 System Specifications:")
        print(f"   • OS: {system_specs['system']} ({system_specs['architecture']})")
        print(f"   • CPU: {system_specs['cpu_name']}")
        print(f"   • Cores: {system_specs['cpu_cores_physical']} physical, {system_specs['cpu_cores_logical']} logical")
        print(f"   • Memory: {system_specs['total_memory_gb']} GB total, {system_specs['available_memory_gb']} GB available")
        print(f"   • Python: {system_specs['python_version']}")
        print(f"   • Platform: {system_specs['platform']}")
        
        timing_results = []
        
        for i in range(max_molecules):
            try:
                # Generate molecules of increasing size (chain lengths 3-50 to reach ~2000 atoms)
                # Use logarithmic scaling for better coverage of large molecules
                if i < max_molecules * 0.3:  # First 30%: small molecules (3-10 carbons)
                    chain_length = 3 + (i % 8)
                elif i < max_molecules * 0.6:  # Next 30%: medium molecules (10-25 carbons)
                    chain_length = 10 + (i % 16)
                else:  # Last 40%: large molecules (25-50 carbons for ~2000 atoms)
                    chain_length = 25 + (i % 26)
                
                # Create carboxylic acid functional group specification
                functional_group_spec = {
                    'group_smiles': 'C(=O)O',
                    'n': 1,
                    'mode': 'attach'
                }
                
                # Generate molecule using generate_random_mol with increased complexity for larger molecules
                complexity_factor = min(chain_length / 10, 3)  # Scale complexity with size
                mol = generate_random_mol(
                    n=chain_length,
                    functional_groups=functional_group_spec,
                    perfluorinated=True,
                    cycle=(i % max(5, int(15/complexity_factor)) == 0),  # More cycles for larger molecules
                    alkene=(i % max(4, int(12/complexity_factor)) == 0),   # More alkenes
                    alkyne=(i % max(6, int(18/complexity_factor)) == 0)   # More alkynes
                )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    
                    # Calculate molecular properties
                    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
                    num_atoms = mol.GetNumAtoms()
                    num_bonds = mol.GetNumBonds()
                    
                    # Run multiple timing iterations for statistical reliability
                    pfas_times = []
                    atlas_times = []
                    pfas_success_count = 0
                    atlas_success_count = 0
                    pfas_correct_count = 0  # Count correct functional group identification
                    atlas_correct_count = 0  # Count correct PFAS classification
                    detected_groups_list = []
                    atlas_classifications = []
                    
                    for iteration in range(iterations):
                        # Test with PFASGroups
                        pfas_result = self.test_with_pfasgroups(smiles)
                        
                        # Check if carboxylic acid (group 33) was correctly detected
                        correctly_detected = 33 in pfas_result.get('detected_groups', [])
                        
                        # Only include timing if functional group was correctly detected
                        if correctly_detected:
                            pfas_times.append(pfas_result['execution_time'])
                            pfas_correct_count += 1
                        
                        if pfas_result['success']:
                            pfas_success_count += 1
                        detected_groups_list.append(pfas_result['detected_groups'])
                        
                        # Test with PFAS-Atlas
                        atlas_result = self.test_with_atlas(smiles)
                        
                        # Check if Atlas correctly classified as PFAS (not "Not PFAS")
                        atlas_correct = atlas_result.get('first_class', '') != 'Not PFAS'
                        
                        # Only include timing if correctly classified as PFAS
                        if atlas_correct:
                            atlas_times.append(atlas_result['execution_time'])
                            atlas_correct_count += 1
                            
                        if atlas_result['success']:
                            atlas_success_count += 1
                        atlas_classifications.append((atlas_result['first_class'], atlas_result['second_class']))
                    
                    # Calculate statistics (only for correctly detected/classified molecules)
                    import numpy as np
                    pfas_time_avg = np.mean(pfas_times) if pfas_times else 0.0
                    pfas_time_std = np.std(pfas_times) if len(pfas_times) > 1 else 0.0
                    atlas_time_avg = np.mean(atlas_times) if atlas_times else 0.0
                    atlas_time_std = np.std(atlas_times) if len(atlas_times) > 1 else 0.0
                    
                    pfas_time_min = min(pfas_times) if pfas_times else 0.0
                    pfas_time_max = max(pfas_times) if pfas_times else 0.0
                    atlas_time_min = min(atlas_times) if atlas_times else 0.0
                    atlas_time_max = max(atlas_times) if atlas_times else 0.0
                    
                    # Most common classification results
                    pfas_success_rate = pfas_success_count / iterations
                    atlas_success_rate = atlas_success_count / iterations
                    most_common_groups = max(set(str(g) for g in detected_groups_list), key=lambda x: str(detected_groups_list).count(x)) if detected_groups_list else []
                    most_common_atlas = max(set(atlas_classifications), key=atlas_classifications.count) if any(atlas_classifications) else (None, None)
                    
                    timing_data = {
                        'molecule_id': i + 1,
                        'chain_length': chain_length,
                        'smiles': smiles,
                        'molecular_weight': mol_weight,
                        'num_atoms': num_atoms,
                        'num_bonds': num_bonds,
                        'iterations': iterations,
                        'pfasgroups_time_avg': pfas_time_avg,
                        'pfasgroups_time_std': pfas_time_std,
                        'pfasgroups_time_min': pfas_time_min,
                        'pfasgroups_time_max': pfas_time_max,
                        'atlas_time_avg': atlas_time_avg,
                        'atlas_time_std': atlas_time_std,
                        'atlas_time_min': atlas_time_min,
                        'atlas_time_max': atlas_time_max,
                        'pfasgroups_success_rate': pfas_success_rate,
                        'atlas_success_rate': atlas_success_rate,
                        'pfasgroups_correct_rate': pfas_correct_count / iterations,
                        'atlas_correct_rate': atlas_correct_count / iterations,
                        'pfasgroups_detected': eval(most_common_groups) if most_common_groups != '[]' else [],
                        'atlas_first_class': most_common_atlas[0],
                        'atlas_second_class': most_common_atlas[1],
                        'system_specs': system_specs
                    }
                    
                    timing_results.append(timing_data)
                    
                    if (i + 1) % 25 == 0:  # Report more frequently for longer tests
                        recent_results = timing_results[-25:]
                        avg_pfas_time = sum(r['pfasgroups_time_avg'] for r in recent_results) / len(recent_results) * 1000
                        avg_atlas_time = sum(r['atlas_time_avg'] for r in recent_results) / len(recent_results) * 1000
                        avg_atoms = sum(r['num_atoms'] for r in recent_results) / len(recent_results)
                        print(f"  ✅ Generated {i + 1}/{max_molecules} molecules | Avg times: PFASGroups {avg_pfas_time:.2f}ms, Atlas {avg_atlas_time:.2f}ms | Avg atoms: {avg_atoms:.1f}")
                        
            except Exception as e:
                print(f"Warning: Failed to generate molecule {i + 1}: {e}")
                continue
        
        # Save timing results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"pfas_timing_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(timing_results, f, indent=2, default=str)
        
        print(f"\n💾 Timing Benchmark Complete!")
        print(f"📊 Total molecules tested: {len(timing_results)}")
        print(f"💾 Results saved: {output_file}")
        
        # Quick statistics with improved analysis
        if timing_results:
            import numpy as np
            
            avg_pfas_times = [r['pfasgroups_time_avg'] * 1000 for r in timing_results]
            avg_atlas_times = [r['atlas_time_avg'] * 1000 for r in timing_results]
            
            pfas_overall_avg = np.mean(avg_pfas_times)
            pfas_overall_std = np.std(avg_pfas_times)
            atlas_overall_avg = np.mean(avg_atlas_times)
            atlas_overall_std = np.std(avg_atlas_times)
            
            print(f"\n📈 Comprehensive Timing Summary ({iterations} iterations per molecule):")
            print(f"   • PFASGroups: {pfas_overall_avg:.2f}±{pfas_overall_std:.2f}ms per molecule")
            print(f"   • PFAS-Atlas: {atlas_overall_avg:.2f}±{atlas_overall_std:.2f}ms per molecule")
            print(f"   • Speed ratio: {atlas_overall_avg/pfas_overall_avg:.1f}x (Atlas/PFASGroups)")
            
            # Enhanced scaling analysis with multiple size bins based on actual molecular size
            atom_ranges = [
                (0, 25, "Very Small"),
                (26, 50, "Small"), 
                (51, 100, "Medium"),
                (101, 200, "Large"),
                (201, 500, "Very Large"),
                (501, 2000, "Extremely Large")
            ]
            
            print(f"\n📏 Enhanced Size Scaling Analysis:")
            for min_atoms, max_atoms, label in atom_ranges:
                molecules_in_range = [r for r in timing_results if min_atoms <= r['num_atoms'] <= max_atoms]
                if molecules_in_range:
                    pfas_avg = np.mean([r['pfasgroups_time_avg'] * 1000 for r in molecules_in_range])
                    atlas_avg = np.mean([r['atlas_time_avg'] * 1000 for r in molecules_in_range])
                    atom_avg = np.mean([r['num_atoms'] for r in molecules_in_range])
                    count = len(molecules_in_range)
                    
                    print(f"   {label} ({min_atoms}-{max_atoms} atoms, n={count}, avg={atom_avg:.0f}): PFASGroups {pfas_avg:.2f}ms, Atlas {atlas_avg:.2f}ms")
            
            # Identify largest molecules tested
            largest_molecules = sorted(timing_results, key=lambda x: x['num_atoms'], reverse=True)[:5]
            print(f"\n🔬 Largest Molecules Tested:")
            for i, mol in enumerate(largest_molecules, 1):
                print(f"   {i}. {mol['num_atoms']} atoms, MW={mol['molecular_weight']:.1f}: PFASGroups {mol['pfasgroups_time_avg']*1000:.2f}ms, Atlas {mol['atlas_time_avg']*1000:.2f}ms")
        
        return timing_results, output_file
    
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
    
    # Ask user what benchmarks to run
    print("🚀 PFAS Benchmark Suite")
    print("Available benchmarks:")
    print("  1. Enhanced Functional Groups Benchmark (comprehensive)")
    print("  2. OECD Benchmark (real-world data)")
    print("  3. Timing Performance Benchmark (scaling analysis)")
    print("  4. All benchmarks")
    
    choice = input("\nChoose benchmark (1-4) or press Enter for all: ").strip()
    
    if choice in ['1', '4', '']:
        print("\nRunning Enhanced Functional Groups Benchmark...")
        enhanced_results, enhanced_file = benchmark.run_enhanced_benchmark()
    else:
        enhanced_results, enhanced_file = None, None
    
    if choice in ['2', '4', '']:
        print("\nRunning OECD Benchmark...")
        oecd_results, oecd_file = benchmark.run_oecd_benchmark()
    else:
        oecd_results, oecd_file = None, None
    
    if choice in ['3', '4', '']:
        print("\nRunning Timing Performance Benchmark...")
        timing_results, timing_file = benchmark.run_timing_benchmark(200, 10)  # 200 molecules, 10 iterations
    else:
        timing_results, timing_file = None, None
    
    print(f"\n🎯 Benchmark Complete! Next steps:")
    if enhanced_file and oecd_file:
        print(f"   Enhanced analysis: python enhanced_analysis.py {enhanced_file} {oecd_file}")
    if timing_file:
        print(f"   Timing analysis: python analyze_timing.py {timing_file}")
    
    print(f"\n📁 Files generated:")
    if enhanced_file:
        print(f"     • {enhanced_file} (functional groups)")
    if oecd_file:
        print(f"     • {oecd_file} (OECD data)")
    if timing_file:
        print(f"     • {timing_file} (timing performance)")

if __name__ == "__main__":
    main()