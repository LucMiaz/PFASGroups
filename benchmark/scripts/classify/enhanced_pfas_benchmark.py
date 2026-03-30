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
import networkx as nx
from rdkit.Chem import Descriptors


# Resolve repo root and prefer local PFASGroups
script_dir = os.path.dirname(os.path.abspath(__file__))
benchmark_dir = os.path.dirname(os.path.dirname(script_dir))
repo_root = os.path.dirname(benchmark_dir)

if repo_root not in sys.path:
    # Prefer local PFASGroups over any installed package
    sys.path.insert(0, repo_root)

try:
    from PFASGroups.parser import parse_mol
    from PFASGroups.generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_groups
    from PFASGroups.core import mol_to_nx, rdkit_disable_log
    PFASGroups_AVAILABLE = True
except ImportError:
    print("❌ PFASGroups not available")
    PFASGroups_AVAILABLE = False

# Try to import PFAS-Atlas
atlas_dir = os.path.join(os.path.dirname(repo_root), 'PFAS-atlas')
try:
    sys.path.append(atlas_dir)  
    from classification_helper.classify_pfas import classify_pfas_molecule
    ATLAS_AVAILABLE = True
    print("✅ PFAS-Atlas available")
except ImportError:
    try:
        sys.path.append(os.path.join(atlas_dir, 'classification_helper'))
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
        pfas_groups_path = os.path.join(repo_root, 'PFASGroups', 'data', 'Halogen_groups_smarts.json')
        with open(pfas_groups_path, 'r') as f:
            self.pfas_groups = json.load(f)
        
        # Load specificity test groups for OECD connections
        specificity_path = os.path.join(repo_root, 'tests', 'test_data','specificity_test_groups.json')
        with open(specificity_path, 'r') as f:
            self.specificity_groups = json.load(f)
        
        # Automatically determine max group ID from loaded data
        max_group_id = max(g['id'] for g in self.pfas_groups)
        
        # Build telomer group set by checking for "telomer" in group names
        self.telomer_groups = set()
        for group in self.pfas_groups:
            if 'telomer' in group.get('name', '').lower():
                self.telomer_groups.add(group['id'])
        
        print(f"✅ Loaded {len(self.pfas_groups)} PFAS groups (max ID: {max_group_id})")
        print(f"✅ Identified {len(self.telomer_groups)} telomer groups: {sorted(self.telomer_groups)}")
        
        # Target groups 29-max_id (excluding 49, 50, 51 which are componentSmarts-only groups)
        # Automatically includes all functional groups and telomers
        self.target_groups = [g for g in range(29, max_group_id + 1)]
        
        # OECD target groups 1-28
        self.oecd_target_groups = list(range(1, 29))
        
        # Functional group definitions following test_examples.py format with proper modes
        # Load group names and properties from JSON automatically
        self.functional_smarts = {}
        
        # Load from PFAS_groups_smarts.json for groups in target range
        # Read test generation data directly from the JSON file
        for group in self.pfas_groups:
            group_id = group['id']
            if group_id in self.target_groups:
                group_name = group.get('name', f'Group {group_id}')
                
                # Check if group has test.generate field
                if 'test' in group and 'generate' in group['test']:
                    generate_info = group['test']['generate']
                    smiles = generate_info.get('smiles', '')
                    mode = generate_info.get('mode', 'attach')
                    is_telomer = generate_info.get('is_telomer', False)
                    
                    # Only add groups that have valid SMILES for generation
                    if smiles:
                        # For telomer groups, strip any hydrocarbon linker prefix from the functional group
                        # Telomers often have SMILES like "CCCC[SiH3]" but the linker is built separately
                        # So we need to extract just the functional group part
                        if is_telomer:
                            import re
                            # Match leading carbons followed by heteroatom or bracket
                            match = re.search(r'^(C+)(\[|[NOSPB])', smiles)
                            if match:
                                # Remove the leading carbon chain
                                carbon_chain = match.group(1)
                                smiles = smiles[len(carbon_chain):]  # Strip the linker
                                print(f"ℹ️  Group {group_id} ({group_name}): Stripped linker '{carbon_chain}' from telomer SMILES")
                                print(f"    Functional group only: {smiles}")
                        
                        entry = {
                            'name': group_name,
                            'smiles': smiles,
                            'mode': mode
                        }
                        
                        if is_telomer:
                            entry['telomer'] = True
                            entry['ch2_range'] = generate_info.get('ch2_range', (2, 8))  # Default CH2 linker range for telomers
                            entry['linker_extra'] = generate_info.get('linker_extra', ['C'])  # Default linker is CH2
                        
                        self.functional_smarts[group_id] = entry
                    else:
                        print(f"⚠️  Group {group_id} ({group_name}) has empty SMILES in test.generate, skipping")
                else:
                    print(f"⚠️  Group {group_id} ({group_name}) missing test.generate field, skipping")
        
        print(f"✅ Loaded {len(self.functional_smarts)} functional group definitions for testing")
        
        # Update target_groups to only include groups that were successfully loaded into functional_smarts
        self.target_groups = [g for g in self.target_groups if g in self.functional_smarts]
        print(f"✅ Target groups updated to {len(self.target_groups)} groups with valid SMARTS patterns")
        
        # Build OECD group mappings after functional_smarts is defined
        self.build_oecd_mappings()
    
    def build_oecd_mappings(self):
        """Build mappings between OECD groups and enhanced functional groups using specificity connections"""
        
        # Map OECD groups (1-28) to their names and base functional groups
        self.oecd_group_info = {}
        self.base_to_oecd = {}
        
        for i, group in enumerate(self.pfas_groups[:28]):
            group_id = i + 1
            group_name = group['name'].lower()
            base_groups = group.get('base_functional_groups', [])
            
            self.oecd_group_info[group_id] = {
                'name': group_name,
                'base_groups': base_groups,
                'smarts1': group.get('smarts1'),
                'smarts2': group.get('smarts2')
            }
            
            # Map base functional groups to OECD group IDs
            for base_group in base_groups:
                base_key = base_group.lower()
                if base_key not in self.base_to_oecd:
                    self.base_to_oecd[base_key] = []
                self.base_to_oecd[base_key].append(group_id)
        
        # Build specificity connections from enhanced groups to OECD groups
        self.enhanced_to_oecd_connections = {}
        
        for conn in self.specificity_groups:
            source = conn['source'].lower()
            target = conn['target'].lower()
            edge_type = conn['edge_type']
            
            # Find source in base functional groups (OECD groups 1-28)
            source_oecd_groups = []
            
            # Direct base group match
            if source in self.base_to_oecd:
                source_oecd_groups.extend(self.base_to_oecd[source])
            
            # Direct OECD group name match
            for group_id, info in self.oecd_group_info.items():
                if source == info['name']:
                    source_oecd_groups.append(group_id)
            
            # Find target in enhanced functional groups (29-51, except 48)
            target_enhanced_groups = []
            for group_id in self.target_groups:
                if group_id in self.functional_smarts and target == self.functional_smarts[group_id]['name'].lower():
                    target_enhanced_groups.append(group_id)
            
            # Also check if target is an OECD group name
            for group_id, info in self.oecd_group_info.items():
                if target == info['name']:
                    target_enhanced_groups.append(group_id)
            
            # Create connections
            for src_group in source_oecd_groups:
                for tgt_group in target_enhanced_groups:
                    if src_group not in self.enhanced_to_oecd_connections:
                        self.enhanced_to_oecd_connections[src_group] = []
                    
                    self.enhanced_to_oecd_connections[src_group].append({
                        'target_group': tgt_group,
                        'edge_type': edge_type,
                        'connection_type': 'specificity'
                    })
        
        print(f"✅ Built OECD mappings: {len(self.oecd_group_info)} OECD groups, {len(self.enhanced_to_oecd_connections)} connections")
        
    def generate_single_group_molecules(self, group_id, count=15):
        """Generate molecules with a single functional group using proper generate_mol functions"""
        
        group_info = self.functional_smarts[group_id]
        molecules = []
        excluded_molecules = []  # Track molecules that fail round-trip SMILES validation
        
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
                    # cycle=(i % 7 == 0),  # Some cyclic molecules
                    # alkene=(i % 5 == 0),  # Some alkenes
                    # alkyne=(i % 6 == 0)   # Some alkynes
                )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)

                    # Validate round-trip SMILES to avoid downstream parse errors
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        excluded_molecules.append({
                            'smiles': smiles,
                            'chain_length': chain_length,
                            'target_group_id': group_id,
                            'target_group_name': group_info['name'],
                            'reason': 'invalid_smiles_roundtrip'
                        })
                        continue
                    
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
    
    def generate_fluorotelomer_molecules(self, group_id, count=15):
        """Generate fluorotelomer molecules with varying perfluorinated chain and CH2 linker lengths
        
        Fluorotelomers have structure: (CF2)n-CF3 + (CH2)m + functional_group
        where n varies the perfluorinated length and m varies the CH2 linker length
        
        Args:
            group_id: Fluorotelomer group ID (60-73)
            count: Number of molecules to generate
        """
        if group_id not in self.functional_smarts:
            print(f"Warning: Group {group_id} not found in functional_smarts")
            return []
        
        group_info = self.functional_smarts[group_id]
        molecules = []
        
        # Check if this is a telomer group
        is_telomer = group_info.get('telomer', False)
        ch2_range = group_info.get('ch2_range', (2, 8))
        linker_extra = group_info.get('linker_extra', [])
        if not is_telomer:
            # For non-telomer groups like 60, 61, generate normally
            return self.generate_single_group_molecules(group_id, count)
        
        # Generate telomer molecules with varying perfluorinated and CH2 lengths
        for i in range(count):
            try:
                # Vary perfluorinated chain length: CF3-(CF2)n where n = 1-5
                n_perfluoro = 2 + (i % 4)  # CF3-CF2 to CF3-CF2-CF2-CF2-CF2 (4:x to 10:x telomers)
                
                # Vary CH2 linker length within specified range
                ch2_min, ch2_max = ch2_range
                ch2_min = max(2, ch2_min)
                ch2_max = max(ch2_min, ch2_max)
                n_ch2 = ch2_min + (i % (ch2_max - ch2_min + 1))
                
                # Build perfluorinated chain: FC(F)(F)C(F)(F)...C(F)(F)
                # Start with FC(F)(F) then add C(F)(F) units
                if n_perfluoro == 1:
                    perfluoro_smiles = "FC(F)(F)"
                else:
                    perfluoro_smiles = "FC(C(F)(F)F)(F)" + "C(F)(F)" * (n_perfluoro - 2)
                
                # Build CH2 linker chain
                # For ethoxylates, add oxygen atoms between some CH2 groups
                # Pattern: -CH2-O-CH2-O-CH2- (ethylene oxide units)
                linker_parts = ["C"]
                j = 0
                for part in linker_extra:
                    linker_parts.append(part)
                    linker_parts.append("C")
                    j += 1
                ch2_remaining = n_ch2 - (1 + j)
                if ch2_remaining > 0:
                    linker_parts.extend(["C"] * ch2_remaining)
                linker_smiles = "".join(linker_parts)
                
                # Add functional group
                func_smiles = group_info['smiles']
                
                # Construct full SMILES by direct concatenation (atoms bond left to right in SMILES)
                # For most groups, structure is: perfluoro + linker + functional_group
                if group_info['mode'] == 'attach':
                    # Functional group attaches at end
                    full_smiles = f"{perfluoro_smiles}{linker_smiles}{func_smiles}"
                else:
                    # For insert mode, functional group goes in the middle (less common for telomers)
                    full_smiles = f"{perfluoro_smiles}{linker_smiles[:n_ch2//2]}{func_smiles}{linker_smiles[n_ch2//2:]}"
                
                # Validate and canonicalize SMILES
                mol = Chem.MolFromSmiles(full_smiles)
                if mol is not None:
                    canonical_smiles = Chem.MolToSmiles(mol)
                    
                    molecule_data = {
                        'group_id': group_id,
                        'group_name': group_info['name'],
                        'smiles': canonical_smiles,
                        'perfluoro_length': n_perfluoro,
                        'ch2_length': n_ch2,
                        'target_groups': [group_id],
                        'generation_type': 'fluorotelomer',
                        'telomer_notation': f"{n_perfluoro*2}:{n_ch2}"  # e.g., "6:2", "8:3"
                    }
                    molecules.append(molecule_data)
                    
            except Exception as e:
                print(f"Warning: Failed to generate fluorotelomer for group {group_id}: {e}")
                continue
        
        return molecules
    
    def generate_oecd_group_molecules(self, group_id, count=10):
        """Generate molecules representing OECD group patterns using existing PFAS structures"""
        
        if group_id not in self.oecd_group_info:
            print(f"Warning: OECD group {group_id} not found")
            return []
        
        group_info = self.oecd_group_info[group_id]
        molecules = []
        
        # For OECD groups, we don't generate new molecules but instead select 
        # existing PFAS structures that would match this group pattern
        # This represents real-world PFAS that fall into these OECD categories
        
        # Use the base functional groups to generate representative molecules
        base_groups = group_info['base_groups']
        
        if not base_groups:
            print(f"Warning: OECD group {group_id} has no base functional groups")
            return []
        
        for i in range(count):
            try:
                # Create molecule representing this OECD category
                chain_length = 4 + (i % 4)  # Varying chain lengths 4-7
                
                # Use the first base functional group for generation
                primary_base = base_groups[0].lower()
                
                # Map base group to enhanced functional group for generation
                enhanced_group_id = None
                for enh_id in self.target_groups:
                    if self.functional_smarts[enh_id]['name'].lower() == primary_base:
                        enhanced_group_id = enh_id
                        break
                
                if enhanced_group_id:
                    # Generate using enhanced functional group system
                    enh_group_info = self.functional_smarts[enhanced_group_id]
                    
                    functional_group_spec = {
                        'group_smiles': enh_group_info['smiles'],
                        'n': 1,
                        'mode': enh_group_info['mode']
                    }
                    
                    # Generate molecule using proper generate_mol function
                    mol = generate_random_mol(
                        n=chain_length,
                        functional_groups=functional_group_spec,
                        perfluorinated=True
                    )
                    
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        
                        molecule_data = {
                            'oecd_group_id': group_id,
                            'group_id': group_id,
                            'group_name': group_info['name'],
                            'smiles': smiles,
                            'chain_length': chain_length,
                            'target_groups': [group_id],
                            'generation_type': 'oecd_group',
                            'base_functional_groups': base_groups
                        }
                        molecules.append(molecule_data)
                
            except Exception as e:
                print(f"Warning: Failed to generate OECD molecule for group {group_id}: {e}")
                continue
        
        return molecules
    
    def test_with_PFASGroups(self, smiles, include_PFAS_definitions=True,
                             limit_effective_graph_resistance=None,
                             compute_component_metrics=True):
        """Test molecule with PFASGroups detection
        
        Args:
            smiles: SMILES string
            include_PFAS_definitions: Whether to include PFAS definitions (True for accuracy, False for specificity)
            limit_effective_graph_resistance: Limit or disable graph resistance computation (None=all, False/0=skip)
            compute_component_metrics: Whether to compute component graph metrics
        """
        
        start_time = time.perf_counter()
        PFASGroups_result = {
            'detected_groups': [],
            'detected_definitions': [],
            'matches': [],
            'success': False,
            'error': None,
            'execution_time': 0.0,
            'include_definitions': include_PFAS_definitions
        }
        
        if not PFASGroups_AVAILABLE:
            PFASGroups_result['error'] = 'PFASGroups not available'
            return PFASGroups_result
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Use parse_mol which returns dict with new format
                results = parse_mol(
                    mol,
                    include_PFAS_definitions=include_PFAS_definitions,
                    limit_effective_graph_resistance=limit_effective_graph_resistance,
                    compute_component_metrics=compute_component_metrics
                )
                
                # Extract groups and definitions from the new dictionary format
                # results is a dict with 'matches' key containing list of match dicts
                group_ids = []
                definition_ids = []
                all_matches = []
                
                if isinstance(results, dict) and 'matches' in results:
                    for match in results['matches']:
                        if match.get('type') == 'HalogenGroup':
                            group_ids.append(match['id'])
                            # Forward all summary metrics from parse_mol
                            _SUMMARY_KEYS = (
                                'mean_branching', 'total_branching', 'sum_component_branching_ratio',
                                'mean_smarts_centrality', 'mean_component_fraction', 'total_components_fraction',
                                'mean_eccentricity', 'median_eccentricity',
                                'mean_diameter', 'mean_radius',
                                'mean_effective_graph_resistance', 'mean_effective_graph_resistance_BDE',
                                'mean_dist_to_barycentre', 'mean_dist_to_centre', 'mean_dist_to_periphery',
                            )
                            all_matches.append({
                                'type': 'group',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('group_name'),
                                'match_count': match.get('match_count'),
                                'components_sizes': match.get('components_sizes', []),
                                'num_components': match.get('num_components', 0),
                                'components_types': match.get('components_types', []),
                                **{k: match.get(k) for k in _SUMMARY_KEYS},
                            })
                        elif match.get('type') == 'PFASdefinition':
                            definition_ids.append(match['id'])
                            all_matches.append({
                                'type': 'definition',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('definition_name')
                            })
                
                PFASGroups_result['detected_groups'] = group_ids
                PFASGroups_result['detected_definitions'] = definition_ids
                PFASGroups_result['matches'] = all_matches
                PFASGroups_result['success'] = len(group_ids) > 0 or len(definition_ids) > 0
                
        except Exception as e:
            PFASGroups_result['error'] = str(e)
        finally:
            PFASGroups_result['execution_time'] = time.perf_counter() - start_time
        
        return PFASGroups_result
    
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
            
            # Use specialized fluorotelomer generation for telomer groups (detected by name)
            if group_id in self.telomer_groups:
                molecules = self.generate_fluorotelomer_molecules(group_id, count=replicates)
            else:
                molecules = self.generate_single_group_molecules(group_id, count=replicates)
            success_count = 0
            total_molecules = len(molecules)
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups
                    pfas_result = self.test_with_PFASGroups(mol_data['smiles'])
                    
                    # Test with PFAS-Atlas  
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'PFASGroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / replicates) * 100
            print(f"  ✅ Generated {success_count}/{total_molecules} molecules ({success_rate:.1f}% success) for group {group_id}")
        
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
        
        # Filter combos to only include groups that exist in functional_smarts
        multi_group_pairs = [c for c in multi_group_pairs if all(g in self.functional_smarts for g in c)]
        multi_group_triplets = [c for c in multi_group_triplets if all(g in self.functional_smarts for g in c)]

        # Test pairs
        print("🔬 Testing 5 functional group pairs (40 molecules each):")
        for combo in multi_group_pairs:
            combo_str = ' + '.join([self.functional_smarts[g]['name'] for g in combo])
            print(f"   • Groups {combo[0]}-{combo[1]}: {combo_str}")
            
            molecules = self.generate_multi_group_molecules(combo, count=replicates)
            success_count = 0
            total_molecules = len(molecules)
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups - accuracy test
                    pfas_result = self.test_with_PFASGroups(mol_data['smiles'], include_PFAS_definitions=True)
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'PFASGroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / replicates) * 100
            print(f"     ✅ Generated {success_count}/{total_molecules} molecules ({success_rate:.1f}% success) for groups {combo[0]}-{combo[1]}")
        
        # Test triplets
        print("\nTesting 5 functional group triplets (10 molecules each):")
        for combo in multi_group_triplets:
            combo_str = ' + '.join([self.functional_smarts[g]['name'] for g in combo])
            print(f"   • Groups {'-'.join(map(str, combo))}: {combo_str}")
            
            molecules = self.generate_multi_group_molecules(combo, count=replicates)
            success_count = 0
            total_molecules = len(molecules)
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups - accuracy test
                    pfas_result = self.test_with_PFASGroups(mol_data['smiles'], include_PFAS_definitions=True)
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'PFASGroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / replicates) * 100
            print(f"     ✅ Generated {success_count}/{total_molecules} molecules ({success_rate:.1f}% success) for groups {'-'.join(map(str, combo))}")
        
        # Part 3: OECD groups (1-28) molecules
        print("\n📋 PART 3: OECD Group Molecules (1-28)")
        print("=" * 45)
        
        oecd_replicates = 15  # Fewer replicates for OECD groups
        
        for group_id in self.oecd_target_groups:
            group_name = self.oecd_group_info[group_id]['name']
            print(f"🧪 Generating molecules for OECD Group {group_id}: {group_name}")
            
            molecules = self.generate_oecd_group_molecules(group_id, count=oecd_replicates)
            success_count = 0
            
            for mol_data in molecules:
                if mol_data:
                    success_count += 1
                    
                    # Test with PFASGroups - accuracy test
                    pfas_result = self.test_with_PFASGroups(mol_data['smiles'], include_PFAS_definitions=True)
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'PFASGroups_result': pfas_result,
                        'atlas_result': atlas_result
                    }
                    
                    all_results.append(result)
            
            success_rate = (success_count / oecd_replicates) * 100
            print(f"  ✅ Generated {success_count}/{oecd_replicates} molecules ({success_rate:.1f}% success)")
        
        # Save results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_enhanced_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        total_molecules = len(all_results)
        enhanced_single_molecules = sum(1 for r in all_results if r['molecule_data'].get('generation_type') == 'single_group')
        oecd_molecules = sum(1 for r in all_results if r['molecule_data'].get('generation_type') == 'oecd_molecule')
        multi_molecules = total_molecules - enhanced_single_molecules - oecd_molecules
        
        print(f"\n💾 Enhanced Benchmark Complete!")
        print(f"📊 Total molecules tested: {total_molecules}")
        print(f"   • Enhanced single group molecules: {enhanced_single_molecules}")
        print(f"   • OECD group molecules: {oecd_molecules}")
        print(f"   • Multi-group molecules: {multi_molecules}")
        print(f"💾 Results saved: {output_file}")
        
        return all_results, output_file
    
    def get_system_specifications(self):
        """Collect system specifications for reproducibility"""
        import platform
        try:
            import psutil
        except ImportError:
            psutil = None
        try:
            import cpuinfo
        except ImportError:
            cpuinfo = None
        
        try:
            cpu_info = cpuinfo.get_cpu_info() if cpuinfo else {}
            cpu_name = cpu_info.get('brand_raw', 'Unknown CPU')
            if psutil:
                cpu_count = psutil.cpu_count(logical=True)
                cpu_count_physical = psutil.cpu_count(logical=False)
            else:
                cpu_count = os.cpu_count()
                cpu_count_physical = os.cpu_count()
        except Exception:
            cpu_name = platform.processor()
            cpu_count = os.cpu_count()
            cpu_count_physical = os.cpu_count()
        
        if psutil:
            memory = psutil.virtual_memory()
            total_memory_gb = round(memory.total / (1024**3), 2)
            available_memory_gb = round(memory.available / (1024**3), 2)
        else:
            total_memory_gb = 0
            available_memory_gb = 0
        
        specs = {
            'system': platform.system(),
            'platform': platform.platform(),
            'architecture': platform.architecture()[0],
            'python_version': platform.python_version(),
            'cpu_name': cpu_name,
            'cpu_cores_logical': cpu_count,
            'cpu_cores_physical': cpu_count_physical,
            'total_memory_gb': total_memory_gb,
            'available_memory_gb': available_memory_gb
        }
        
        return specs
    
    def calculate_molecule_complexity(self, mol):
        """Calculate molecule complexity using graph component metrics
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            dict: Dictionary of complexity metrics
        """
        try:
            # Convert RDKit molecule to NetworkX graph using PFASGroups.core function
            G = mol_to_nx(mol)
            
            # Calculate graph metrics
            metrics = {}
            
            # Eccentricity metrics (if graph is connected)
            if nx.is_connected(G):
                eccentricities = nx.eccentricity(G)
                metrics['diameter'] = nx.diameter(G)
                metrics['radius'] = nx.radius(G)
                metrics['avg_eccentricity'] = sum(eccentricities.values()) / len(eccentricities)
                metrics['max_eccentricity'] = max(eccentricities.values())
            else:
                # For disconnected graphs, use largest component
                largest_cc = max(nx.connected_components(G), key=len)
                subgraph = G.subgraph(largest_cc)
                eccentricities = nx.eccentricity(subgraph)
                metrics['diameter'] = nx.diameter(subgraph)
                metrics['radius'] = nx.radius(subgraph)
                metrics['avg_eccentricity'] = sum(eccentricities.values()) / len(eccentricities)
                metrics['max_eccentricity'] = max(eccentricities.values())
                metrics['num_components'] = nx.number_connected_components(G)
            
            # Additional complexity metrics
            metrics['avg_degree'] = sum(dict(G.degree()).values()) / G.number_of_nodes()
            metrics['density'] = nx.density(G)
            metrics['num_cycles'] = len(nx.cycle_basis(G))
            
            # Betweenness centrality (complexity of pathways)
            betweenness = nx.betweenness_centrality(G)
            metrics['avg_betweenness'] = sum(betweenness.values()) / len(betweenness)
            metrics['max_betweenness'] = max(betweenness.values())
            
            # Clustering coefficient (local connectivity)
            metrics['avg_clustering'] = nx.average_clustering(G)
            
            # Overall complexity score (weighted combination)
            metrics['complexity_score'] = (
                metrics['diameter'] * 0.2 +
                metrics['avg_eccentricity'] * 0.2 +
                metrics['avg_degree'] * 0.15 +
                metrics['density'] * 10 * 0.15 +
                metrics['num_cycles'] * 0.15 +
                metrics['avg_betweenness'] * 10 * 0.15
            )
            
            return metrics
            
        except Exception as e:
            print(f"Warning: Failed to calculate complexity metrics: {e}")
            return {
                'diameter': 0,
                'radius': 0,
                'avg_eccentricity': 0,
                'max_eccentricity': 0,
                'avg_degree': 0,
                'density': 0,
                'num_cycles': 0,
                'avg_betweenness': 0,
                'max_betweenness': 0,
                'avg_clustering': 0,
                'complexity_score': 0
            }
    
    def load_latest_timing_results(self, profile_label=None):
        """Load the most recent timing benchmark results from data directory
        
        Returns:
            List of timing results, or empty list if no previous results found
        """
        import glob
        
        # Find all timing benchmark JSON files
        if profile_label:
            pattern = os.path.join(script_dir, '..', 'data', f'pfas_timing_benchmark_{profile_label}_*.json')
        else:
            pattern = os.path.join(script_dir, '..', 'data', 'pfas_timing_benchmark_*.json')
        files = glob.glob(pattern)
        
        if not files:
            print("ℹ️  No previous timing results found")
            return []
        
        # Get the most recent file
        latest_file = max(files, key=os.path.getmtime)
        
        try:
            with open(latest_file, 'r') as f:
                results = json.load(f)
            
            print(f"✅ Loaded {len(results)} previous timing results from {os.path.basename(latest_file)}")
            return results
        except Exception as e:
            print(f"⚠️  Failed to load previous results: {e}")
            return []
    
    def run_timing_benchmark(self, max_molecules=4000, iterations=5, reuse_previous=False,
                             limit_effective_graph_resistance=None,
                             compute_component_metrics=True,
                             profile_label="full"):
        """Run timing benchmark comparing PFASGroups vs PFAS-Atlas with diverse functional groups
        
        Args:
            max_molecules: Number of molecules to test (default 2500)
            iterations: Number of iterations for averaging (default 5)
            reuse_previous: If True, load previous results and add new ones (default False)
            limit_effective_graph_resistance: Limit or disable graph resistance computation (None=all, False/0=skip)
            compute_component_metrics: Whether to compute component graph metrics
            profile_label: Label for this timing profile (used in output file name)
        """
        
        # Get system specifications
        system_specs = self.get_system_specifications()
        
        # Get available functional groups (IDs 29-114)
        excluded_timing_groups = []#{62, 103}
        available_groups = []
        for group_id in range(29, 115):
            if group_id in self.functional_smarts:
                group_info = self.functional_smarts[group_id]
                # Check if group has test generation data
                if 'smiles' in group_info and 'mode' in group_info:
                    if group_id in excluded_timing_groups:
                        continue
                    available_groups.append(group_id)
        
        print("\n🚀 TIMING PERFORMANCE BENCHMARK WITH GRAPH COMPLEXITY METRICS")
        print("=" * 70)
        print(f"🧪 Timing profile: {profile_label}")
        print(f"📊 Testing {max_molecules} molecules using grid-based sampling (IDs 29-114)")
        print(f"🔄 Running {iterations} iterations per molecule for statistical averaging")
        print(f"📏 Chain length range: 5-200 carbon atoms (binned for systematic coverage)")
        print(f"🧬 Available functional groups: {len(available_groups)} groups")
        print(f"📈 Using graph metrics for complexity: eccentricity, diameter, betweenness, etc.")
        if not compute_component_metrics:
            print("⚙️  Component graph metrics disabled")
        if limit_effective_graph_resistance is False or limit_effective_graph_resistance == 0:
            print("⚙️  Effective graph resistance disabled")
        print(f"\n💻 System Specifications:")
        print(f"   • OS: {system_specs['system']} ({system_specs['architecture']})")
        print(f"   • CPU: {system_specs['cpu_name']}")
        print(f"   • Cores: {system_specs['cpu_cores_physical']} physical, {system_specs['cpu_cores_logical']} logical")
        print(f"   • Memory: {system_specs['total_memory_gb']} GB total, {system_specs['available_memory_gb']} GB available")
        print(f"   • Python: {system_specs['python_version']}")
        print(f"   • Platform: {system_specs['platform']}")
        
        # Create grid of (chain_length, group_id) combinations
        # Define chain length bins for better coverage
        chain_length_bins = list(range(5, 50, 5)) + list(range(50, 100, 10)) + list(range(100, 201, 20))
        
        # Create all possible grid points
        grid_points = []
        for chain_length in chain_length_bins:
            for group_id in available_groups:
                grid_points.append((chain_length, group_id))
        
        # Shuffle grid randomly
        random.shuffle(grid_points)
        
        # If max_molecules is larger than grid, allow repetition by cycling through shuffled grid
        grid_cycle = grid_points * ((max_molecules // len(grid_points)) + 1)
        
        print(f"🔲 Created grid: {len(chain_length_bins)} chain lengths × {len(available_groups)} groups = {len(grid_points)} unique combinations")
        print(f"📍 Chain length bins: {chain_length_bins[:5]}... to ...{chain_length_bins[-3:]}")
        
        # Load previous results if requested
        if reuse_previous:
            timing_results = self.load_latest_timing_results(profile_label)
            if timing_results:
                print(f"🔄 Starting from {len(timing_results)} previous results, will add {max_molecules} more")
                print(f"📊 Target total: {len(timing_results) + max_molecules} molecules")
                max_molecules_total = len(timing_results) + max_molecules
            else:
                print(f"⚠️  No previous results found, starting fresh")
                max_molecules_total = max_molecules
        else:
            timing_results = []
            max_molecules_total = max_molecules
        
        excluded_molecules = []  # Track excluded molecules for review
        excluded_wrong_detection = 0
        failed_generations = 0
        grid_index = 0
        attempts = 0
        
        while len(timing_results) < max_molecules_total and attempts < max_molecules_total * 3:
            attempts += 1
            try:
                # Get next grid point (chain_length, group_id)
                chain_length, group_id = grid_cycle[grid_index]
                grid_index += 1
                group_info = self.functional_smarts[group_id]
                
                # Create functional group specification from JSON data
                functional_group_spec = {
                    'group_smiles': group_info['smiles'],
                    'n': 1,
                    'mode': group_info['mode']
                }
                
                # Generate molecule with random structural features
                if group_id in self.telomer_groups:
                    telomer_mols = self.generate_fluorotelomer_molecules(group_id, count=1)
                    if telomer_mols:
                        smiles = telomer_mols[0]['smiles']
                        mol = Chem.MolFromSmiles(smiles)
                    else:
                        mol = None
                else:
                    mol = generate_random_mol(
                        n=chain_length,
                        functional_groups=functional_group_spec,
                        perfluorinated=True,
                        cycle=(random.random() < 0.2),  # 20% chance of cycle
                        alkene=(random.random() < 0.15),  # 15% chance of alkene
                        alkyne=(random.random() < 0.1)   # 10% chance of alkyne
                    )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    
                    # Pre-validate: Check that PFASGroups detects the molecule
                    validation_result = self.test_with_PFASGroups(
                        smiles,
                        include_PFAS_definitions=True,
                        limit_effective_graph_resistance=limit_effective_graph_resistance,
                        compute_component_metrics=compute_component_metrics
                    )
                    if not validation_result['success'] or group_id not in validation_result['detected_groups']:
                        excluded_wrong_detection += 1
                        # Save excluded molecule for review
                        excluded_molecules.append({
                            'smiles': smiles,
                            'chain_length': chain_length,
                            'target_group_id': group_id,
                            'target_group_name': group_info['name'],
                            'detected_groups': validation_result['detected_groups'],
                            'success': validation_result['success'],
                            'reason': 'target_group_not_detected' if validation_result['success'] else 'validation_failed'
                        })
                        continue  # Skip molecules where the target group is not correctly identified
                    
                    # Calculate molecular properties
                    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
                    num_atoms = mol.GetNumAtoms()
                    num_bonds = mol.GetNumBonds()
                    
                    # Calculate graph complexity metrics
                    complexity_metrics = self.calculate_molecule_complexity(mol)
                    
                    # Calculate RDKit descriptors
                    num_rings = Descriptors.RingCount(mol)
                    num_aromatic_rings = Descriptors.NumAromaticRings(mol)
                    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
                    tpsa = Descriptors.TPSA(mol)
                    
                    # Run multiple timing iterations for statistical reliability
                    pfas_times = []
                    atlas_times = []
                    pfas_success_count = 0
                    atlas_success_count = 0
                    detected_groups_list = []
                    atlas_classifications = []
                    
                    for iteration in range(iterations):
                        # Test with PFASGroups - accuracy test
                        pfas_result = self.test_with_PFASGroups(
                            smiles,
                            include_PFAS_definitions=True,
                            limit_effective_graph_resistance=limit_effective_graph_resistance,
                            compute_component_metrics=compute_component_metrics
                        )
                        pfas_times.append(pfas_result['execution_time'])
                        if pfas_result['success'] and group_id in pfas_result['detected_groups']:
                            pfas_success_count += 1
                        detected_groups_list.append(pfas_result['detected_groups'])
                        
                        # Test with PFAS-Atlas
                        atlas_result = self.test_with_atlas(smiles)
                        atlas_times.append(atlas_result['execution_time'])
                        if atlas_result['success']:
                            atlas_success_count += 1
                        atlas_classifications.append((atlas_result['first_class'], atlas_result['second_class']))
                    
                    # Calculate statistics
                    import numpy as np
                    pfas_time_avg = np.mean(pfas_times)
                    pfas_time_std = np.std(pfas_times)
                    atlas_time_avg = np.mean(atlas_times)
                    atlas_time_std = np.std(atlas_times)
                    
                    # Most common classification results
                    pfas_success_rate = pfas_success_count / iterations
                    atlas_success_rate = atlas_success_count / iterations
                    most_common_groups = max(set(str(g) for g in detected_groups_list), key=lambda x: str(detected_groups_list).count(x)) if detected_groups_list else []
                    most_common_atlas = max(set(atlas_classifications), key=atlas_classifications.count) if any(atlas_classifications) else (None, None)
                    
                    timing_data = {
                        'molecule_id': len(timing_results) + 1,
                        'timing_profile': profile_label,
                        'compute_component_metrics': compute_component_metrics,
                        'limit_effective_graph_resistance': limit_effective_graph_resistance,
                        'chain_length': chain_length,
                        'target_group_id': group_id,
                        'target_group_name': group_info['name'],
                        'smiles': smiles,
                        'molecular_weight': mol_weight,
                        'num_atoms': num_atoms,
                        'num_bonds': num_bonds,
                        'num_rings': num_rings,
                        'num_aromatic_rings': num_aromatic_rings,
                        'num_rotatable_bonds': num_rotatable_bonds,
                        'tpsa': tpsa,
                        'complexity_diameter': complexity_metrics['diameter'],
                        'complexity_radius': complexity_metrics['radius'],
                        'complexity_avg_eccentricity': complexity_metrics['avg_eccentricity'],
                        'complexity_max_eccentricity': complexity_metrics['max_eccentricity'],
                        'complexity_avg_degree': complexity_metrics['avg_degree'],
                        'complexity_density': complexity_metrics['density'],
                        'complexity_num_cycles': complexity_metrics['num_cycles'],
                        'complexity_avg_betweenness': complexity_metrics['avg_betweenness'],
                        'complexity_max_betweenness': complexity_metrics['max_betweenness'],
                        'complexity_avg_clustering': complexity_metrics['avg_clustering'],
                        'complexity_score': complexity_metrics['complexity_score'],
                        'iterations': iterations,
                        'PFASGroups_time_avg': pfas_time_avg,
                        'PFASGroups_time_std': pfas_time_std,
                        'PFASGroups_time_min': min(pfas_times),
                        'PFASGroups_time_max': max(pfas_times),
                        'atlas_time_avg': atlas_time_avg,
                        'atlas_time_std': atlas_time_std,
                        'atlas_time_min': min(atlas_times),
                        'atlas_time_max': max(atlas_times),
                        'PFASGroups_success_rate': pfas_success_rate,
                        'atlas_success_rate': atlas_success_rate,
                        'PFASGroups_detected': eval(most_common_groups) if most_common_groups != '[]' else [],
                        'atlas_first_class': most_common_atlas[0],
                        'atlas_second_class': most_common_atlas[1],
                        'system_specs': system_specs
                    }
                    
                    timing_results.append(timing_data)
                    
                    if len(timing_results) % 50 == 0:  # Report every 50 molecules
                        recent_results = timing_results[-50:]
                        avg_pfas_time = sum(r['PFASGroups_time_avg'] for r in recent_results) / len(recent_results) * 1000
                        avg_atlas_time = sum(r['atlas_time_avg'] for r in recent_results) / len(recent_results) * 1000
                        avg_atoms = sum(r['num_atoms'] for r in recent_results) / len(recent_results)
                        avg_complexity = sum(r['complexity_score'] for r in recent_results) / len(recent_results)
                        print(f"  ✅ Generated {len(timing_results)}/{max_molecules_total} molecules | "
                              f"Avg times: PFASGroups {avg_pfas_time:.2f}ms, Atlas {avg_atlas_time:.2f}ms | "
                              f"Avg atoms: {avg_atoms:.1f} | Avg complexity: {avg_complexity:.2f}")
                else:
                    failed_generations += 1
                    
            except Exception as e:
                print(f"Warning: Failed to generate molecule (attempt {attempts}): {e}")
                failed_generations += 1
                continue
        
        # Save timing results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_timing_benchmark_{profile_label}_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(timing_results, f, indent=2, default=str)
        
        # Save excluded molecules for review
        if excluded_molecules:
            excluded_file = f"data/pfas_timing_excluded_{timestamp}.json"
            with open(excluded_file, 'w') as f:
                json.dump(excluded_molecules, f, indent=2, default=str)
            print(f"\n⚠️  EXCLUDED MOLECULES SAVED FOR REVIEW")
            print(f"📋 {len(excluded_molecules)} molecules excluded (target group not detected)")
            print(f"💾 Excluded molecules saved: {excluded_file}")
            print(f"\n📌 Please review excluded molecules to verify:")
            print(f"   • No bugs in molecular generation")
            print(f"   • SMARTS patterns are correctly defined")
            print(f"   • Expected behavior for edge cases")
        
        print(f"\n💾 Timing Benchmark Complete!")
        print(f"📊 Total molecules tested: {len(timing_results)}")
        print(f"⚠️ Molecules excluded (wrong detection): {excluded_wrong_detection}")
        print(f"❌ Failed generations: {failed_generations}")
        print(f"✅ Success rate: {len(timing_results)/(len(timing_results)+excluded_wrong_detection+failed_generations)*100:.1f}%")
        print(f"💾 Results saved: {output_file}")
        
        # Quick statistics with improved analysis
        if timing_results:
            import numpy as np
            
            avg_pfas_times = [r['PFASGroups_time_avg'] * 1000 for r in timing_results]
            avg_atlas_times = [r['atlas_time_avg'] * 1000 for r in timing_results]
            
            pfas_overall_avg = np.mean(avg_pfas_times)
            pfas_overall_std = np.std(avg_pfas_times)
            atlas_overall_avg = np.mean(avg_atlas_times)
            atlas_overall_std = np.std(avg_atlas_times)
            
            print(f"\n📈 Comprehensive Timing Summary ({iterations} iterations per molecule):")
            print(f"   • PFASGroups: {pfas_overall_avg:.2f}±{pfas_overall_std:.2f}ms per molecule")
            print(f"   • PFAS-Atlas: {atlas_overall_avg:.2f}±{atlas_overall_std:.2f}ms per molecule")
            if pfas_overall_avg > atlas_overall_avg:
                print(f"   • Speed: PFASGroups is {pfas_overall_avg/atlas_overall_avg:.1f}x slower than Atlas")
            else:
                print(f"   • Speed: PFASGroups is {atlas_overall_avg/pfas_overall_avg:.1f}x faster than Atlas")
            
            # Enhanced scaling analysis with multiple size bins
            atom_ranges = [
                (0, 50, "Very Small"),
                (51, 100, "Small"), 
                (101, 200, "Medium"),
                (201, 500, "Large"),
                (501, 1000, "Very Large"),
                (1001, 2000, "Extremely Large")
            ]
            
            print(f"\n📏 Enhanced Size Scaling Analysis:")
            for min_atoms, max_atoms, label in atom_ranges:
                molecules_in_range = [r for r in timing_results if min_atoms <= r['num_atoms'] <= max_atoms]
                if molecules_in_range:
                    pfas_avg = np.mean([r['PFASGroups_time_avg'] * 1000 for r in molecules_in_range])
                    atlas_avg = np.mean([r['atlas_time_avg'] * 1000 for r in molecules_in_range])
                    atom_avg = np.mean([r['num_atoms'] for r in molecules_in_range])
                    count = len(molecules_in_range)
                    
                    print(f"   {label} ({min_atoms}-{max_atoms} atoms, n={count}, avg={atom_avg:.0f}): PFASGroups {pfas_avg:.2f}ms, Atlas {atlas_avg:.2f}ms")
            
            # Identify largest molecules tested
            largest_molecules = sorted(timing_results, key=lambda x: x['num_atoms'], reverse=True)[:5]
            print(f"\n🔬 Largest Molecules Tested:")
            for i, mol in enumerate(largest_molecules, 1):
                print(f"   {i}. {mol['num_atoms']} atoms, MW={mol['molecular_weight']:.1f}: PFASGroups {mol['PFASGroups_time_avg']*1000:.2f}ms, Atlas {mol['atlas_time_avg']*1000:.2f}ms")
            
            # Complexity metrics analysis
            complexity_scores = [r['complexity_score'] for r in timing_results]
            diameters = [r['complexity_diameter'] for r in timing_results]
            eccentricities = [r['complexity_avg_eccentricity'] for r in timing_results]
            
            print(f"\n📊 Graph Complexity Metrics Summary:")
            print(f"   • Complexity Score: {np.mean(complexity_scores):.2f}±{np.std(complexity_scores):.2f} (range: {min(complexity_scores):.2f}-{max(complexity_scores):.2f})")
            print(f"   • Diameter: {np.mean(diameters):.2f}±{np.std(diameters):.2f} (range: {min(diameters)}-{max(diameters)})")
            print(f"   • Avg Eccentricity: {np.mean(eccentricities):.2f}±{np.std(eccentricities):.2f}")
            
            # Most complex molecules
            most_complex = sorted(timing_results, key=lambda x: x['complexity_score'], reverse=True)[:5]
            print(f"\n🧬 Most Complex Molecules:")
            for i, mol in enumerate(most_complex, 1):
                print(f"   {i}. Complexity={mol['complexity_score']:.2f}, {mol['num_atoms']} atoms, diameter={mol['complexity_diameter']}, "
                      f"cycles={mol['complexity_num_cycles']}: {mol['target_group_name']}")
        
        return timing_results, output_file
    
    def run_non_fluorinated_benchmark(self, molecules_per_group=50):
        """Run benchmark to test if systems correctly exclude non-fluorinated functional groups
        
        Tests carboxylic acid, sulfonic acid, and ether functional groups without fluorination
        to verify that PFAS detection systems properly exclude non-PFAS molecules.
        
        Args:
            molecules_per_group: Number of molecules to test per functional group
        """
        
        print("\n🚀 NON-FLUORINATED FUNCTIONAL GROUPS BENCHMARK")
        print("=" * 55)
        print(f"🧪 Testing specificity of functional group detection systems")
        print(f"📊 Testing 3 functional groups with {molecules_per_group} molecules each")
        print("🎯 Expected: Target functional groups should NOT be detected in non-fluorinated molecules")
        
        # Test functional groups: carboxylic acid, sulfonic acid, ether
        test_groups = {
            33: {'name': 'carboxylic acid', 'smiles': 'C(=O)O', 'mode': 'attach'},
            36: {'name': 'sulfonic acid', 'smiles': 'S(=O)(=O)O', 'mode': 'attach'},
            31: {'name': 'ether', 'smiles': 'O', 'mode': 'insert'}
        }
        
        all_results = []
        
        for group_id, group_info in test_groups.items():
            print(f"\n🧪 Testing non-fluorinated {group_info['name']} (Group {group_id})")
            
            group_results = {
                'group_id': group_id,
                'group_name': group_info['name'],
                'molecules_tested': 0,
                'PFASGroups_detections': 0,
                'atlas_detections': 0,
                'molecules': []
            }
            
            for i in range(molecules_per_group):
                try:
                    # Generate NON-fluorinated molecules manually using RDKit
                    chain_length = 3 + (i % 8)  # Chain lengths 3-10
                    
                    # Create base hydrocarbon chain without fluorination
                    if group_id == 33:  # carboxylic acid
                        # Create simple carboxylic acids without fluorine
                        base_smiles_options = [
                            "CCCC(=O)O",  # butanoic acid
                            "CCC(=O)O",   # propanoic acid
                            "CCCCC(=O)O", # pentanoic acid
                            "CCCCCC(=O)O", # hexanoic acid
                            "CC(C)C(=O)O", # isobutyric acid
                            "CC(C)CC(=O)O", # 3-methylbutanoic acid
                            "CCCCCCC(=O)O", # heptanoic acid
                            "CCCCCCCC(=O)O", # octanoic acid
                            "CC(=O)O", # acetic acid
                            "CCCCCCCCC(=O)O", # nonanoic acid
                            "CCCCCCCCCC(=O)O", # decanoic acid
                            "c1ccccc1C(=O)O", # benzoic acid
                        ]
                        smiles = base_smiles_options[i % len(base_smiles_options)]
                        
                    elif group_id == 36:  # sulfonic acid
                        # Create simple sulfonic acids without fluorine
                        base_smiles_options = [
                            "CCS(=O)(=O)O",     # ethanesulfonic acid
                            "CCCS(=O)(=O)O",    # propanesulfonic acid
                            "CCCCS(=O)(=O)O",   # butanesulfonic acid
                            "CS(=O)(=O)O",      # methanesulfonic acid
                            "CCCCCS(=O)(=O)O",  # pentanesulfonic acid
                            "CC(C)CS(=O)(=O)O", # 2-methylpropanesulfonic acid
                            "CCCCCCS(=O)(=O)O", # hexanesulfonic acid
                            "c1ccc(S(=O)(=O)O)cc1", # benzenesulfonic acid
                            "CCCCCCCS(=O)(=O)O", # heptanesulfonic acid
                            "CCCCCCCCS(=O)(=O)O", # octanesulfonic acid
                            "c1ccc(C)c(S(=O)(=O)O)c1", # toluenesulfonic acid
                            "c1cc(S(=O)(=O)O)ccc1C", # p-toluenesulfonic acid
                        ]
                        smiles = base_smiles_options[i % len(base_smiles_options)]
                        
                    elif group_id == 31:  # ether
                        # Create simple ethers without fluorine
                        base_smiles_options = [
                            "CCOC",        # ethyl methyl ether
                            "CCOCC",       # diethyl ether
                            "CCCOCC",      # ethyl propyl ether
                            "CCCOCCCC",    # dipropyl ether
                            "CC(C)OC",     # isopropyl methyl ether
                            "CCCCOC",      # butyl methyl ether
                            "CCOCCC",      # ethyl propyl ether
                            "CCCCOCC",    # butyl ethyl ether
                            "CCCCOCCCC", # dibutyl ether
                            "CCCCCOC", # pentyl methyl ether
                            "c1ccccc1OC", # anisole (methoxybenzene)
                            "CCCCCCOC", # hexyl methyl ether
                        ]
                        smiles = base_smiles_options[i % len(base_smiles_options)]
                    
                    # Verify molecule is valid and non-fluorinated
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None and 'F' not in smiles:
                        
                        group_results['molecules_tested'] += 1
                        
                        # Test with PFASGroups - specificity test (should NOT detect target group)
                        pfas_result = self.test_with_PFASGroups(smiles, include_PFAS_definitions=False)
                        
                        # Check if target functional group was incorrectly detected
                        if pfas_result['success'] and group_id in pfas_result['detected_groups']:
                            group_results['PFASGroups_detections'] += 1
                        
                        # Test with PFAS-Atlas - should NOT detect as PFAS
                        atlas_result = self.test_with_atlas(smiles)
                        
                        if atlas_result['success']:
                            group_results['atlas_detections'] += 1
                        
                        group_results['molecules'].append({
                            'smiles': smiles,
                            'PFASGroups_detected': pfas_result['success'],
                            'PFASGroups_detected_groups': pfas_result['detected_groups'],
                            'PFASGroups_target_group_detected': group_id in pfas_result['detected_groups'],
                            'atlas_detected': atlas_result['success'],
                            'atlas_first_class': atlas_result.get('first_class'),
                            'atlas_second_class': atlas_result.get('second_class')
                        })
                        group_results['molecules'].append({
                            'smiles': smiles,
                            'chain_length': len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']),
                            'contains_fluorine': False,
                            'PFASGroups_detected': pfas_result['success'],
                            'PFASGroups_groups': pfas_result['detected_groups'],
                            'PFASGroups_target_group_detected': group_id in pfas_result['detected_groups'] if pfas_result['success'] else False,
                            'atlas_detected': atlas_result['success'] and atlas_result['first_class'] != 'Not PFAS',
                            'atlas_first_class': atlas_result['first_class'],
                            'atlas_second_class': atlas_result['second_class']
                        })
                        
                except Exception as e:
                    print(f"Warning: Failed to generate non-fluorinated molecule {i+1} for group {group_id}: {e}")
                    continue
            
            pfas_false_positive_rate = (group_results['PFASGroups_detections'] / max(group_results['molecules_tested'], 1)) * 100
            atlas_false_positive_rate = (group_results['atlas_detections'] / max(group_results['molecules_tested'], 1)) * 100
            
            print(f"   ✅ Generated {group_results['molecules_tested']} non-fluorinated molecules")
            print(f"   ⚠️  PFASGroups false positives (target group {group_id}): {group_results['PFASGroups_detections']}/{group_results['molecules_tested']} ({pfas_false_positive_rate:.1f}%)")
            print(f"   ⚠️  Atlas false positives (any PFAS class): {group_results['atlas_detections']}/{group_results['molecules_tested']} ({atlas_false_positive_rate:.1f}%)")
            
            all_results.append(group_results)
        
        # Save non-fluorinated benchmark results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_non_fluorinated_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        # Summary statistics
        total_molecules = sum(r['molecules_tested'] for r in all_results)
        total_pfas_false_positives = sum(r['PFASGroups_detections'] for r in all_results)
        total_atlas_false_positives = sum(r['atlas_detections'] for r in all_results)
        
        overall_pfas_fpr = (total_pfas_false_positives / max(total_molecules, 1)) * 100
        overall_atlas_fpr = (total_atlas_false_positives / max(total_molecules, 1)) * 100
        
        print(f"\n💾 Non-Fluorinated Benchmark Complete!")
        print(f"📊 Total molecules tested: {total_molecules}")
        print(f"🎯 Expected result: 0% detection of target functional groups (these are NOT PFAS molecules)")
        print(f"📈 PFASGroups false positive rate (target groups): {overall_pfas_fpr:.1f}% ({total_pfas_false_positives}/{total_molecules})")
        print(f"📈 PFAS-Atlas false positive rate (any PFAS): {overall_atlas_fpr:.1f}% ({total_atlas_false_positives}/{total_molecules})")
        print(f"💾 Results saved: {output_file}")
        
        if overall_pfas_fpr == 0 and overall_atlas_fpr == 0:
            print(f"✅ PERFECT: Both systems correctly excluded all non-fluorinated molecules!")
        elif overall_pfas_fpr < 5 and overall_atlas_fpr < 5:
            print(f"✅ GOOD: Both systems have low false positive rates (<5%)")
        else:
            print(f"⚠️  WARNING: High false positive rates detected - systems may need calibration")
        
        return all_results, output_file
    
    def run_complex_branched_benchmark(self, molecules_per_test=20):
        """Run benchmark to test detection of complex branched PFAS molecules
        
        Tests both systems' ability to correctly identify highly branched, complex PFAS structures
        that might challenge pattern recognition algorithms.
        
        Args:
            molecules_per_test: Number of molecules to test per complexity category
        """
        
        print("\nCOMPLEX BRANCHED PFAS BENCHMARK")
        print("=" * 55)
        print(f"Testing detection of complex branched PFAS molecules")
        print(f"Testing {molecules_per_test} molecules per complexity category")
        print("Expected: Both systems should correctly identify these as PFAS")
        
        # Define complex branched PFAS molecules with different functional groups and complexity levels
        test_molecules = {
            'highly_branched_carboxylic_acid': {
                'smiles': 'C(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)C(=O)O',
                'description': 'Highly branched perfluorocarboxylic acid',
                'expected_PFASGroups': [34],  # carboxylic acid (group 34)
                'complexity': 'very_high'
            },
            'branched_sulfonic_acid': {
                'smiles': 'C(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)S(=O)(=O)O',
                'description': 'Branched perfluorosulfonic acid',
                'expected_PFASGroups': [37],  # sulfonic acid (group 37)
                'complexity': 'high'
            },
            'branched_ether_chain': {
                'smiles': 'C(C(F)(F)C(F)(F)F)(C(F)(F)F)OC(C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)F',
                'description': 'Branched perfluoroether',
                'expected_PFASGroups': [32],  # ether (group 32)
                'complexity': 'high'
            },
            'multi_functional_branched': {
                'smiles': 'C(C(F)(F)C(=O)O)(C(F)(F)C(F)(F)F)OC(F)(F)C(F)(F)S(=O)(=O)O',
                'description': 'Multi-functional branched PFAS (carboxylic acid + ether + sulfonic acid)',
                'expected_PFASGroups': [32, 34, 37],  # ether (32), carboxylic acid (34), sulfonic acid (37)
                'complexity': 'very_high'
            },
            'cyclic_branched': {
                'smiles': 'C1(C(F)(F)F)(C(F)(F)F)C(C(F)(F)C(F)(F)F)C(C(F)(F)F)C(C(=O)O)C1(F)F',
                'description': 'Cyclic branched perfluorocarboxylic acid',
                'expected_PFASGroups': [56, 57],  # Perfluoro cyclic (56), Polyfluoro cyclic (57)
                'complexity': 'high'
            },
            'aromatic_branched': {
                'smiles': 'c1(F)c(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)c(F)c(C(F)(F)C(=O)O)c(F)c1F',
                'description': 'Aromatic branched PFAS with carboxylic acid',
                'expected_PFASGroups': [22, 34, 51, 52, 55, 58, 59],  # side-chain fluorinated aromatics (22), carboxylic acid (34), perfluoroalkyl (51), polyfluoroalkyl (52), side-chain aromatics (55), perfluoroaryl (58), polyfluoroaryl (59)
                'complexity': 'very_high'
            }
        }
        
        all_results = []
        total_tests = 0
        successful_detections = {'PFASGroups': 0, 'atlas': 0}
        
        for test_name, mol_info in test_molecules.items():
            print(f"\nTesting {test_name}: {mol_info['description']}")
            print(f"   SMILES: {mol_info['smiles']}")
            print(f"   Complexity: {mol_info['complexity']}")
            print(f"   Expected PFASGroups: {mol_info['expected_PFASGroups']}")
            
            test_results = {
                'test_name': test_name,
                'description': mol_info['description'],
                'smiles': mol_info['smiles'],
                'complexity': mol_info['complexity'],
                'expected_PFASGroups': mol_info['expected_PFASGroups'],
                'molecules_tested': 0,
                'PFASGroups_detections': 0,
                'atlas_detections': 0,
                'PFASGroups_correct_detections': 0,
                'molecules': []
            }
            
            # Test the molecule multiple times to check consistency
            for i in range(molecules_per_test):
                try:
                    # Verify molecule is valid
                    mol = Chem.MolFromSmiles(mol_info['smiles'])
                    if mol is None:
                        print(f"   ❌ Invalid SMILES: {mol_info['smiles']}")
                        continue
                    
                    # Calculate molecule properties
                    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
                    num_atoms = mol.GetNumAtoms()
                    num_heavy_atoms = mol.GetNumHeavyAtoms()
                    num_fluorines = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F')
                    
                    test_results['molecules_tested'] += 1
                    total_tests += 1
                    
                    # Test with PFASGroups - accuracy test
                    pfas_result = self.test_with_PFASGroups(mol_info['smiles'], include_PFAS_definitions=True)
                    PFASGroups_detected = pfas_result['success']
                    detected_groups = pfas_result['detected_groups']
                    
                    if PFASGroups_detected:
                        test_results['PFASGroups_detections'] += 1
                        successful_detections['PFASGroups'] += 1
                        
                        # Check if expected groups were detected
                        expected_groups_detected = all(group in detected_groups for group in mol_info['expected_PFASGroups'])
                        if expected_groups_detected:
                            test_results['PFASGroups_correct_detections'] += 1
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_info['smiles'])
                    atlas_detected = atlas_result['success']
                    
                    if atlas_detected:
                        test_results['atlas_detections'] += 1
                        successful_detections['atlas'] += 1
                    
                    molecule_result = {
                        'iteration': i + 1,
                        'smiles': mol_info['smiles'],
                        'molecular_weight': mol_weight,
                        'num_atoms': num_atoms,
                        'num_heavy_atoms': num_heavy_atoms,
                        'num_fluorines': num_fluorines,
                        'PFASGroups_detected': PFASGroups_detected,
                        'PFASGroups_groups': detected_groups,
                        'PFASGroups_expected_groups': mol_info['expected_PFASGroups'],
                        'PFASGroups_correct': all(group in detected_groups for group in mol_info['expected_PFASGroups']) if PFASGroups_detected else False,
                        'atlas_detected': atlas_detected,
                        'atlas_first_class': atlas_result['first_class'],
                        'atlas_second_class': atlas_result['second_class'],
                        'PFASGroups_execution_time': pfas_result['execution_time'],
                        'atlas_execution_time': atlas_result['execution_time']
                    }
                    
                    test_results['molecules'].append(molecule_result)
                    
                except Exception as e:
                    print(f"   Warning: Failed to test iteration {i+1}: {e}")
                    continue
            
            # Calculate success rates for this test
            PFASGroups_detection_rate = (test_results['PFASGroups_detections'] / max(test_results['molecules_tested'], 1)) * 100
            PFASGroups_accuracy_rate = (test_results['PFASGroups_correct_detections'] / max(test_results['molecules_tested'], 1)) * 100
            atlas_detection_rate = (test_results['atlas_detections'] / max(test_results['molecules_tested'], 1)) * 100
            
            print(f"   Tested {test_results['molecules_tested']} iterations")
            print(f"   PFASGroups detection: {test_results['PFASGroups_detections']}/{test_results['molecules_tested']} ({PFASGroups_detection_rate:.1f}%)")
            print(f"   PFASGroups correct groups: {test_results['PFASGroups_correct_detections']}/{test_results['molecules_tested']} ({PFASGroups_accuracy_rate:.1f}%)")
            print(f"   Atlas detection: {test_results['atlas_detections']}/{test_results['molecules_tested']} ({atlas_detection_rate:.1f}%)")
            
            all_results.append(test_results)
        
        # Save complex branched benchmark results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_complex_branched_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        # Summary statistics
        overall_PFASGroups_rate = (successful_detections['PFASGroups'] / max(total_tests, 1)) * 100
        overall_atlas_rate = (successful_detections['atlas'] / max(total_tests, 1)) * 100
        
        print(f"\nComplex Branched Benchmark Complete!")
        print(f"Total tests: {total_tests} ({len(test_molecules)} molecule types × {molecules_per_test} iterations)")
        print(f"Expected result: High detection rates for complex PFAS structures")
        print(f"PFASGroups overall detection rate: {overall_PFASGroups_rate:.1f}% ({successful_detections['PFASGroups']}/{total_tests})")
        print(f"PFAS-Atlas overall detection rate: {overall_atlas_rate:.1f}% ({successful_detections['atlas']}/{total_tests})")
        print(f"Results saved: {output_file}")
        
        if overall_PFASGroups_rate >= 90 and overall_atlas_rate >= 90:
            print(f"EXCELLENT: Both systems show excellent detection of complex PFAS structures!")
        elif overall_PFASGroups_rate >= 75 and overall_atlas_rate >= 75:
            print(f"GOOD: Both systems show good detection rates for complex structures")
        else:
            print(f"CONCERNING: Lower than expected detection rates for complex PFAS structures")
        
        return all_results, output_file
    
    def run_oecd_benchmark(self, max_molecules=None):
        """Run benchmark on OECD data using groups 1-28"""
        
        print("\nOECD PFAS BENCHMARK (Groups 1-28)")
        print("=" * 45)
        
        # Load OECD data
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Use relative path: ../../../../PFAS-atlas from scripts/classify folder
        oecd_file = os.path.join(script_dir, '..', '..', '..', '..', 'PFAS-atlas', 'input_data', 'OECD_4000', 'step3_OECD_Class_0812.csv')
        oecd_file = os.path.normpath(oecd_file)  # Normalize the path
        try:
            import pandas as pd
            oecd_data = pd.read_csv(oecd_file)
            if max_molecules:
                oecd_data = oecd_data.head(max_molecules)
            print(f"Loaded {len(oecd_data)} OECD molecules")
        except Exception as e:
            print(f"❌ Error loading OECD data: {e}")
            print(f"   Expected path: {oecd_file}")
            import sys
            sys.exit(1)
            
        all_results = []
        
        # Process OECD molecules
        print(f"Testing {len(oecd_data)} OECD molecules with PFASGroups (groups 1-28)")
        
        for idx, row in oecd_data.iterrows():
            smiles = row['SMILES']
            first_class = row.get('First_Class', 'Unknown')
            second_class = row.get('Second_Class', 'Unknown')
            
            if pd.isna(smiles) or smiles.strip() == '':
                continue
                
            # Test with PFASGroups - accuracy test
            pfas_result = self.test_with_PFASGroups(smiles, include_PFAS_definitions=True)
            
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
                'PFASGroups_result': pfas_result,
                'atlas_result': atlas_result
            }
            
            all_results.append(result)
            
            if (idx + 1) % 100 == 0:
                print(f"  ✅ Processed {idx + 1}/{len(oecd_data)} molecules")
        
        # Save OECD results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_oecd_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        print(f"\nOECD Benchmark Complete!")
        print(f"Total molecules tested: {len(all_results)}")
        print(f"Results saved: {output_file}")
        
        return all_results, output_file
@rdkit_disable_log()
def main():
    """Main function to run enhanced benchmark"""
    
    def env_int(name, default):
        value = os.getenv(name)
        if value is None or value == "":
            return default
        try:
            return int(value)
        except ValueError:
            return default

    def timing_profile_options(profile_label):
        profile_key = (profile_label or "full").strip().lower()
        profiles = {
            "full": {
                "limit_effective_graph_resistance": None,
                "compute_component_metrics": True
            },
            "no_resistance": {
                "limit_effective_graph_resistance": False,
                "compute_component_metrics": True
            },
            "no_metrics": {
                "limit_effective_graph_resistance": False,
                "compute_component_metrics": False
            }
        }
        if profile_key not in profiles:
            print(f"⚠️  Unknown timing profile '{profile_label}', defaulting to 'full'")
            profile_key = "full"
        return profile_key, profiles[profile_key]
    
    benchmark = EnhancedPFASBenchmark()
    
    # Environment overrides for quick runs
    replicates = env_int("PFAS_BENCH_REPLICATES", 40)
    timing_molecules = env_int("PFAS_BENCH_TIMING_MOLECULES", 2500)
    timing_iterations = env_int("PFAS_BENCH_TIMING_ITERATIONS", 5)
    nonfluor_count = env_int("PFAS_BENCH_NONFLUOR", 50)
    complex_count = env_int("PFAS_BENCH_COMPLEX", 50)
    oecd_limit = env_int("PFAS_BENCH_OECD_LIMIT", 0)
    if oecd_limit == 0:
        oecd_limit = None
    
    # Ask user what benchmarks to run
    print("🚀 PFAS Benchmark Suite")
    print("Available benchmarks:")
    print("  1. Enhanced Functional Groups Benchmark (comprehensive)")
    print("  2. OECD Benchmark (real-world data)")
    print("  3. Timing Performance Benchmark (scaling analysis)")
    print("  4. Non-Fluorinated Exclusion Benchmark (false positive test)")
    print("  5. Complex Branched PFAS Benchmark (positive control validation)")
    print("  6. All benchmarks")
    
    choice = input("\nChoose benchmark (1-6) or press Enter for all: ").strip()
    
    if choice in ['1', '6', '']:
        print("\nRunning Enhanced Functional Groups Benchmark...")
        enhanced_results, enhanced_file = benchmark.run_enhanced_benchmark(replicates)
    else:
        enhanced_results, enhanced_file = None, None
    
    if choice in ['2', '6', '']:
        print("\nRunning OECD Benchmark...")
        oecd_results, oecd_file = benchmark.run_oecd_benchmark(oecd_limit)
    else:
        oecd_results, oecd_file = None, None
    
    if choice in ['3', '6', '']:
        print("\nRunning Timing Performance Benchmark...")
        
        # Ask if user wants to reuse previous results
        reuse = input("⏪ Reuse previous timing results and add more? (y/N): ").strip().lower()
        reuse_previous = reuse in ['y', 'yes']

        timing_profile_env = os.getenv("PFAS_BENCH_TIMING_PROFILE", "full")
        timing_profile, timing_options = timing_profile_options(timing_profile_env)
        print(f"🧪 Using timing profile: {timing_profile}")
        
        timing_results, timing_file = benchmark.run_timing_benchmark(
            timing_molecules,
            timing_iterations,
            reuse_previous=reuse_previous,
            profile_label=timing_profile,
            **timing_options
        )
    else:
        timing_results, timing_file = None, None
    
    if choice in ['4', '6', '']:
        print("\nRunning Non-Fluorinated Exclusion Benchmark...")
        nonfluor_results, nonfluor_file = benchmark.run_non_fluorinated_benchmark(nonfluor_count)
    else:
        nonfluor_results, nonfluor_file = None, None
    
    if choice in ['5', '6', '']:
        print("\nRunning Complex Branched PFAS Benchmark...")
        complex_results, complex_file = benchmark.run_complex_branched_benchmark(complex_count)
    else:
        complex_results, complex_file = None, None
    
    print(f"\nBenchmark Complete! Next steps:")
    if enhanced_file and oecd_file:
        print(f"   Enhanced analysis: python enhanced_analysis.py {enhanced_file} {oecd_file}")
    if timing_file:
        print(f"   Timing analysis: python analyze_timing.py {timing_file}")
    if nonfluor_file:
        print(f"   Non-fluorinated analysis: python analyze_nonfluorinated.py {nonfluor_file}")
    if complex_file:
        print(f"   Complex molecules analysis: python analyze_complex.py {complex_file}")
    if complex_file:
        print(f"   Complex molecules analysis: python analyze_complex.py {complex_file}")
    
    print(f"\n📁 Files generated:")
    if enhanced_file:
        print(f"     • {enhanced_file} (functional groups)")
    if oecd_file:
        print(f"     • {oecd_file} (OECD data)")
    if timing_file:
        print(f"     • {timing_file} (timing performance)")
    if nonfluor_file:
        print(f"     • {nonfluor_file} (non-fluorinated exclusion)")
    if complex_file:
        print(f"     • {complex_file} (complex branched PFAS)")

if __name__ == "__main__":
    main()
