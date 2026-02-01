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
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(os.path.dirname(script_dir))

sys.path.append(parent_dir)

try:
    from PFASgroups.parser import parse_mol
    from PFASgroups.generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol, append_functional_groups
    PFASGROUPS_AVAILABLE = True
except ImportError:
    print("❌ PFASGroups not available")
    PFASGROUPS_AVAILABLE = False

# Try to import PFAS-Atlas
atlas_dir = os.path.join(os.path.dirname(parent_dir), 'PFAS-atlas')
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
        pfas_groups_path = os.path.join(parent_dir, 'PFASgroups', 'data', 'PFAS_groups_smarts.json')
        with open(pfas_groups_path, 'r') as f:
            self.pfas_groups = json.load(f)
        
        # Load specificity test groups for OECD connections
        specificity_path = os.path.join(parent_dir, 'tests', 'specificity_test_groups.json')
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
        
        # Target groups 29-max_id (excluding 49, 50 which are smartsPath-only groups)
        # Automatically includes all functional groups and telomers
        self.target_groups = [g for g in range(29, max_group_id + 1) if g not in [49, 50]]
        
        # OECD target groups 1-28
        self.oecd_target_groups = list(range(1, 29))
        
        # Functional group definitions following test_examples.py format with proper modes
        # Load group names and properties from JSON automatically
        self.functional_smarts = {}
        
        # Load from PFAS_groups_smarts.json for groups in target range
        for group in self.pfas_groups:
            group_id = group['id']
            if group_id in self.target_groups:
                group_name = group.get('name', f'Group {group_id}')
                is_telomer = group_id in self.telomer_groups
                
                # Try to determine appropriate SMILES and mode from group data
                # For telomers, use basic functional group patterns
                # This is a simplified version - actual molecule generation will use the full SMARTS
                smiles_map = {
                    # Keep existing mappings for groups 29-59
                    29: {'smiles': 'O[H]', 'mode': 'attach'},
                    30: {'smiles': 'C(=O)', 'mode': 'insert'},
                    31: {'smiles': 'O', 'mode': 'insert'},
                    32: {'smiles': 'C(=O)OC', 'mode': 'insert'},
                    33: {'smiles': 'C(=O)O', 'mode': 'attach'},
                    34: {'smiles': 'C(=O)N', 'mode': 'insert'},
                    35: {'smiles': 'C(=O)Cl', 'mode': 'attach'},
                    36: {'smiles': 'S(=O)(=O)O', 'mode': 'attach'},
                    37: {'smiles': 'SO[H]', 'mode': 'attach'},
                    38: {'smiles': 'S(=O)O[H]', 'mode': 'attach'},
                    39: {'smiles': 'P(=O)(O)O', 'mode': 'attach'},
                    40: {'smiles': 'P(=O)O', 'mode': 'attach'},
                    41: {'smiles': 'Cl', 'mode': 'attach'},
                    42: {'smiles': 'Br', 'mode': 'attach'},
                    43: {'smiles': 'I', 'mode': 'attach'},
                    44: {'smiles': 'S(=O)(=O)N', 'mode': 'insert'},
                    45: {'smiles': 'c1ncc[nH]1', 'mode': 'attach'},
                    46: {'smiles': 'c1ncccc1', 'mode': 'attach'},
                    47: {'smiles': 'c1ccc2OC(F)(F)Oc2c1', 'mode': 'attach'},
                    48: {'smiles': 'N', 'mode': 'insert'},
                    51: {'smiles': 'C(F)=C(F)', 'mode': 'insert'},
                    52: {'smiles': 'C#C', 'mode': 'insert'},
                    53: {'smiles': 'c1ccccc1', 'mode': 'attach'},
                    54: {'smiles': 'C1(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)(F)', 'mode': 'attach'},
                    55: {'smiles': 'C1(F)C(F)(F)(F)C(H)C(F)C(F)(F)C1(F)', 'mode': 'attach'},
                    56: {'smiles': 'c1c(F)c(F)c(F)c(F)c1F', 'mode': 'attach'},
                    57: {'smiles': 'c1c(F)c(F)ccc1F', 'mode': 'attach'},
                    58: {'smiles': 'OO', 'mode': 'insert'},
                    59: {'smiles': 'C(=O)OOC(=O)', 'mode': 'insert'},
                    # Groups 60-66 non-telomers
                    60: {'smiles': 'C(=O)OOC(=O)', 'mode': 'insert'},  # Benzoyl peroxides
                    61: {'smiles': '[Si](C)(C)C', 'mode': 'attach'},  # Silane
                    62: {'smiles': '[Si](Cl)(Cl)Cl', 'mode': 'attach'},  # Trichlorosilane
                    63: {'smiles': 'OC(=O)C=C', 'mode': 'attach'},  # Acrylate
                    64: {'smiles': 'OC(=O)C(=C)C', 'mode': 'attach'},  # Methacrylate
                    65: {'smiles': 'N(C)(C)CC(=O)O', 'mode': 'attach'},  # Betaine
                    66: {'smiles': 'SN#C', 'mode': 'attach'},  # Thiocyanic acid
                    67: {'smiles': 'C(=O)SC(=O)C(O)C(=O)O', 'mode': 'attach'},  # Thia keto propanoic acid
                    68: {'smiles': 'OC1C(O)C(OC(C1O)C(=O)O)O', 'mode': 'attach'},  # Glucuronate
                    # Telomer groups 69+ - use simplified patterns for generation
                    69: {'smiles': '[Si]([H])([H])[H]', 'mode': 'attach'},  # FT silane
                    70: {'smiles': '[Si](Cl)(Cl)Cl', 'mode': 'attach'},  # FT trichlorosilane
                    71: {'smiles': 'I', 'mode': 'attach'},  # FT iodide
                    72: {'smiles': 'C(=O)[H]', 'mode': 'attach'},  # FT aldehyde
                    73: {'smiles': 'C(=O)O', 'mode': 'attach'},  # FT carboxylic acids
                    74: {'smiles': 'O', 'mode': 'attach'},  # FT ethoxylates
                    75: {'smiles': 'S(=O)(=O)O', 'mode': 'attach'},  # FT sulfonic acid
                    76: {'smiles': 'OP(=O)(O)O', 'mode': 'attach'},  # FT monophosphate
                    77: {'smiles': 'OP(=O)(O)OP(=O)(O)O', 'mode': 'attach'},  # FT diphosphate
                    78: {'smiles': 'OP(=O)(O)OP(=O)(O)OP(=O)(O)O', 'mode': 'attach'},  # FT triphosphate
                    79: {'smiles': 'OC(=O)C=C', 'mode': 'attach'},  # FT acrylate
                    80: {'smiles': 'OC(=O)C(=C)C', 'mode': 'attach'},  # FT methacrylate
                    81: {'smiles': 'O', 'mode': 'attach'},  # FT alcohol
                    82: {'smiles': 'S', 'mode': 'insert'},  # FT sulfure
                    83: {'smiles': 'S(=O)(=O)N', 'mode': 'insert'},  # FT sulfonamide
                    84: {'smiles': '[Si](OC)(OC)OC', 'mode': 'attach'},  # FT silyl
                    85: {'smiles': 'N(C)(C)CC(=O)O', 'mode': 'attach'},  # FT betaine
                    86: {'smiles': 'OC(=O)C', 'mode': 'attach'},  # FT ester
                    87: {'smiles': 'SN#C', 'mode': 'attach'},  # FT thiocyanic acid
                    88: {'smiles': 'C=C', 'mode': 'insert'},  # FT alkene
                    89: {'smiles': 'N(C)(C)C', 'mode': 'attach'},  # FT trimethylamine
                    90: {'smiles': 'S(=O)O', 'mode': 'attach'},  # FT sulfinic acid
                    91: {'smiles': 'SO', 'mode': 'attach'},  # FT sulfenic acid
                    92: {'smiles': 'OS(=O)(=O)O', 'mode': 'attach'},  # FT sulfuric acid
                    93: {'smiles': 'C', 'mode': 'attach'},  # FT methyl
                    94: {'smiles': 'OC1C(O)C(OC(C1O)C(=O)O)O', 'mode': 'attach'},  # FT glucuronic acid
                    95: {'smiles': 'C(O)S(=O)(=O)O', 'mode': 'attach'},  # FT formaldehyde bisulfite
                    96: {'smiles': 'C(=O)O', 'mode': 'attach'},  # FT unsaturated carboxylic acids (with alkene)
                    97: {'smiles': 'C(=O)SC(=O)C(O)C(=O)O', 'mode': 'attach'},  # FT thia keto propanoic acid
                    # Non-telomer groups 98+
                    98: {'smiles': 'NCC(=O)O', 'mode': 'attach'},  # Glycine
                    # More telomers 99+
                    99: {'smiles': 'NCC(=O)O', 'mode': 'attach'},  # FT glycine
                    100: {'smiles': 'S(=O)(=O)NCCO', 'mode': 'attach'},  # Sulfonamidoethanol
                    101: {'smiles': 'S(=O)(=O)NCCO', 'mode': 'attach'},  # FT sulfonamidoethanol
                    102: {'smiles': 'OCC(=O)O', 'mode': 'attach'},  # FT ether carboxylic acids
                    103: {'smiles': 'S(=O)(=O)CCC(=O)O', 'mode': 'attach'},  # FT sulfonyl propanoic acid
                    104: {'smiles': 'S(=O)(=O)CCC(=O)O', 'mode': 'attach'},  # Sulfonyl propanoic acid
                    105: {'smiles': 'S(=O)CCC(=O)NCC(C)(C)CS(=O)(=O)O', 'mode': 'attach'},  # FT sulfinyl amido sulfonic acid
                    106: {'smiles': 'S(=O)CCC(=O)NCC(C)(C)CS(=O)(=O)O', 'mode': 'attach'},  # Sulfinyl amido sulfonic acid
                    107: {'smiles': 'S(=O)', 'mode': 'attach'},  # FT sulfinyl
                    108: {'smiles': 'S(=O)(=O)', 'mode': 'attach'},  # FT sulfone
                    109: {'smiles': 'C(=O)SCC(O)C(=O)O', 'mode': 'attach'},  # FT acetylsulfanylhydroxypropanoic acid
                    110: {'smiles': 'C(=O)SC(O)C(=O)O', 'mode': 'attach'},  # FT acetylhydroxyethanethioate
                    111: {'smiles': 'NCCC(=O)O', 'mode': 'attach'},  # FT amino propanoic acid
                    112: {'smiles': 'NCC[N+](C)(C)C', 'mode': 'attach'},  # FT amino ethyl trimethyl ammonium
                    113: {'smiles': 'C', 'mode': 'attach'},  # Telomers (generic)
                }
                
                if group_id in smiles_map:
                    entry = {
                        'name': group_name,
                        'smiles': smiles_map[group_id]['smiles'],
                        'mode': smiles_map[group_id]['mode']
                    }
                    
                    if is_telomer:
                        entry['telomer'] = True
                        entry['ch2_range'] = (2, 8)
                        if 'ethoxylate' in group_name.lower():
                            entry['ethoxylate'] = True
                    
                    self.functional_smarts[group_id] = entry
        
        print(f"✅ Loaded {len(self.functional_smarts)} functional group definitions for testing")
        
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
                if target == self.functional_smarts[group_id]['name'].lower():
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
        is_ethoxylate = group_info.get('ethoxylate', False)
        ch2_range = group_info.get('ch2_range', (2, 8))
        
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
                n_ch2 = ch2_min + (i % (ch2_max - ch2_min + 1))
                
                # Build perfluorinated chain: FC(F)(F)C(F)(F)...C(F)(F)
                # Start with FC(F)(F) then add C(F)(F) units
                if n_perfluoro == 1:
                    perfluoro_smiles = "FC(F)(F)"
                else:
                    perfluoro_smiles = "FC(F)(F)" + "C(F)(F)" * (n_perfluoro - 1)
                
                # Build CH2 linker chain
                if is_ethoxylate:
                    # For ethoxylates, add oxygen atoms between some CH2 groups
                    # Pattern: -CH2-O-CH2-O-CH2- (ethylene oxide units)
                    linker_parts = ["C"]
                    for j in range(1, n_ch2):
                        if (j % 2 == 1):  # Add O every other position
                            linker_parts.append("OC")
                        else:
                            linker_parts.append("C")
                    linker_smiles = "".join(linker_parts)
                else:
                    # Regular CH2 chain
                    linker_smiles = "C" * n_ch2
                
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
    
    def test_with_pfasgroups(self, smiles, include_PFAS_definitions=True):
        """Test molecule with PFASGroups detection
        
        Args:
            smiles: SMILES string
            include_PFAS_definitions: Whether to include PFAS definitions (True for accuracy, False for specificity)
        """
        
        start_time = time.perf_counter()
        pfasgroups_result = {
            'detected_groups': [],
            'detected_definitions': [],
            'matches': [],
            'success': False,
            'error': None,
            'execution_time': 0.0,
            'include_definitions': include_PFAS_definitions
        }
        
        if not PFASGROUPS_AVAILABLE:
            pfasgroups_result['error'] = 'PFASGroups not available'
            return pfasgroups_result
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Use parse_mol which returns dict with new format
                results = parse_mol(mol, include_PFAS_definitions=include_PFAS_definitions)
                
                # Extract groups and definitions from the new dictionary format
                # results is a dict with 'matches' key containing list of match dicts
                group_ids = []
                definition_ids = []
                all_matches = []
                
                if isinstance(results, dict) and 'matches' in results:
                    for match in results['matches']:
                        if match.get('type') == 'PFASgroup':
                            group_ids.append(match['id'])
                            all_matches.append({
                                'type': 'group',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('group_name'),
                                'match_count': match.get('match_count'),
                                'components_sizes': match.get('components_sizes', []),
                                'num_components': match.get('num_components', 0),
                                'components_types': match.get('components_types', []),
                                # Summary metrics
                                'mean_eccentricity': match.get('mean_eccentricity', 0.0),
                                'mean_diameter': match.get('mean_diameter', float('nan')),
                                'mean_radius': match.get('mean_radius', float('nan'))
                            })
                        elif match.get('type') == 'PFASdefinition':
                            definition_ids.append(match['id'])
                            all_matches.append({
                                'type': 'definition',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('definition_name')
                            })
                
                pfasgroups_result['detected_groups'] = group_ids
                pfasgroups_result['detected_definitions'] = definition_ids
                pfasgroups_result['matches'] = all_matches
                pfasgroups_result['success'] = len(group_ids) > 0 or len(definition_ids) > 0
                
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
            
            # Use specialized fluorotelomer generation for telomer groups (detected by name)
            if group_id in self.telomer_groups:
                molecules = self.generate_fluorotelomer_molecules(group_id, count=replicates)
            else:
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
                    
                    # Test with PFASGroups - accuracy test
                    pfas_result = self.test_with_pfasgroups(mol_data['smiles'], include_PFAS_definitions=True)
                    
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
                    
                    # Test with PFASGroups - accuracy test
                    pfas_result = self.test_with_pfasgroups(mol_data['smiles'], include_PFAS_definitions=True)
                    
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
                    pfas_result = self.test_with_pfasgroups(mol_data['smiles'], include_PFAS_definitions=True)
                    
                    # Test with PFAS-Atlas
                    atlas_result = self.test_with_atlas(mol_data['smiles'])
                    
                    result = {
                        'molecule_data': mol_data,
                        'pfasgroups_result': pfas_result,
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
        enhanced_single_molecules = len(self.target_groups) * replicates
        oecd_molecules = len(self.oecd_target_groups) * oecd_replicates
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
        excluded_wrong_detection = 0
        failed_generations = 0
        
        for i in range(max_molecules):
            try:
                # Generate molecules of increasing size (chain lengths 3-50 to reach ~2000 atoms)
                # Use logarithmic scaling for better coverage of large molecules
                current_target = len(timing_results)
                if current_target < max_molecules * 0.3:  # First 30%: small molecules (3-10 carbons)
                    chain_length = 3 + (i % 8)
                elif current_target < max_molecules * 0.6:  # Next 30%: medium molecules (10-25 carbons)
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
                    # cycle=(i % max(5, int(15/complexity_factor)) == 0),  # More cycles for larger molecules
                    # alkene=(i % max(4, int(12/complexity_factor)) == 0),   # More alkenes
                    #alkyne=(i % max(6, int(18/complexity_factor)) == 0)   # More alkynes
                )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    
                    # Pre-validate: Check that PFASGroups correctly identifies carboxylic acid (group 33)
                    validation_result = self.test_with_pfasgroups(smiles, include_PFAS_definitions=True)
                    if not validation_result['success'] or 33 not in validation_result['detected_groups']:
                        excluded_wrong_detection += 1
                        continue  # Skip molecules where carboxylic acid is not correctly identified
                    
                    # Calculate molecular properties
                    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
                    num_atoms = mol.GetNumAtoms()
                    num_bonds = mol.GetNumBonds()
                    
                    # Run multiple timing iterations for statistical reliability
                    pfas_times = []
                    atlas_times = []
                    pfas_success_count = 0
                    atlas_success_count = 0
                    detected_groups_list = []
                    atlas_classifications = []
                    
                    for iteration in range(iterations):
                        # Test with PFASGroups - accuracy test
                        pfas_result = self.test_with_pfasgroups(smiles, include_PFAS_definitions=True)
                        pfas_times.append(pfas_result['execution_time'])
                        if pfas_result['success'] and 33 in pfas_result['detected_groups']:  # Ensure carboxylic acid is detected
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
                        'molecule_id': i + 1,
                        'chain_length': chain_length,
                        'smiles': smiles,
                        'molecular_weight': mol_weight,
                        'num_atoms': num_atoms,
                        'num_bonds': num_bonds,
                        'iterations': iterations,
                        'pfasgroups_time_avg': pfas_time_avg,
                        'pfasgroups_time_std': pfas_time_std,
                        'pfasgroups_time_min': min(pfas_times),
                        'pfasgroups_time_max': max(pfas_times),
                        'atlas_time_avg': atlas_time_avg,
                        'atlas_time_std': atlas_time_std,
                        'atlas_time_min': min(atlas_times),
                        'atlas_time_max': max(atlas_times),
                        'pfasgroups_success_rate': pfas_success_rate,
                        'atlas_success_rate': atlas_success_rate,
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
                else:
                    failed_generations += 1
                    
            except Exception as e:
                print(f"Warning: Failed to generate molecule {i + 1}: {e}")
                failed_generations += 1
                continue
            
            i += 1
        
        # Save timing results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_timing_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(timing_results, f, indent=2, default=str)
        
        print(f"\n💾 Timing Benchmark Complete!")
        print(f"📊 Total molecules tested: {len(timing_results)}")
        print(f"⚠️ Molecules excluded (wrong detection): {excluded_wrong_detection}")
        print(f"❌ Failed generations: {failed_generations}")
        print(f"✅ Success rate: {len(timing_results)/(len(timing_results)+excluded_wrong_detection+failed_generations)*100:.1f}%")
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
                'pfasgroups_detections': 0,
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
                        
                        test_results['molecules_tested'] += 1
                        
                        # Test with PFASGroups - specificity test (should NOT detect target group)
                        pfas_result = self.test_with_pfasgroups(smiles, include_PFAS_definitions=False)
                        
                        # Check if target functional group was incorrectly detected
                        if pfas_result['success'] and group_id in pfas_result['detected_groups']:
                            test_results['pfasgroups_detections'] += 1
                        
                        # Test with PFAS-Atlas - should NOT detect as PFAS
                        atlas_result = self.test_with_atlas(smiles)
                        
                        if atlas_result['success']:
                            test_results['atlas_detections'] += 1
                        
                        test_results['molecules'].append({
                            'smiles': smiles,
                            'pfasgroups_detected': pfas_result['success'],
                            'pfasgroups_detected_groups': pfas_result['detected_groups'],
                            'pfasgroups_target_group_detected': group_id in pfas_result['detected_groups'],
                            'atlas_detected': atlas_result['success'],
                            'atlas_first_class': atlas_result.get('first_class'),
                            'atlas_second_class': atlas_result.get('second_class')
                        })
                        group_results['molecules'].append({
                            'smiles': smiles,
                            'chain_length': len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']),
                            'contains_fluorine': False,
                            'pfasgroups_detected': pfas_result['success'],
                            'pfasgroups_groups': pfas_result['detected_groups'],
                            'pfasgroups_target_group_detected': group_id in pfas_result['detected_groups'] if pfas_result['success'] else False,
                            'atlas_detected': atlas_result['success'] and atlas_result['first_class'] != 'Not PFAS',
                            'atlas_first_class': atlas_result['first_class'],
                            'atlas_second_class': atlas_result['second_class']
                        })
                        
                except Exception as e:
                    print(f"Warning: Failed to generate non-fluorinated molecule {i+1} for group {group_id}: {e}")
                    continue
            
            pfas_false_positive_rate = (group_results['pfasgroups_detections'] / max(group_results['molecules_tested'], 1)) * 100
            atlas_false_positive_rate = (group_results['atlas_detections'] / max(group_results['molecules_tested'], 1)) * 100
            
            print(f"   ✅ Generated {group_results['molecules_tested']} non-fluorinated molecules")
            print(f"   ⚠️  PFASGroups false positives (target group {group_id}): {group_results['pfasgroups_detections']}/{group_results['molecules_tested']} ({pfas_false_positive_rate:.1f}%)")
            print(f"   ⚠️  Atlas false positives (any PFAS class): {group_results['atlas_detections']}/{group_results['molecules_tested']} ({atlas_false_positive_rate:.1f}%)")
            
            all_results.append(group_results)
        
        # Save non-fluorinated benchmark results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_non_fluorinated_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        # Summary statistics
        total_molecules = sum(r['molecules_tested'] for r in all_results)
        total_pfas_false_positives = sum(r['pfasgroups_detections'] for r in all_results)
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
        
        print("\n🚀 COMPLEX BRANCHED PFAS BENCHMARK")
        print("=" * 55)
        print(f"🧪 Testing detection of complex branched PFAS molecules")
        print(f"📊 Testing {molecules_per_test} molecules per complexity category")
        print("🎯 Expected: Both systems should correctly identify these as PFAS")
        
        # Define complex branched PFAS molecules with different functional groups and complexity levels
        test_molecules = {
            'highly_branched_carboxylic_acid': {
                'smiles': 'C(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)C(=O)O',
                'description': 'Highly branched perfluorocarboxylic acid',
                'expected_pfasgroups': [33],  # carboxylic acid
                'complexity': 'very_high'
            },
            'branched_sulfonic_acid': {
                'smiles': 'C(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)S(=O)(=O)O',
                'description': 'Branched perfluorosulfonic acid',
                'expected_pfasgroups': [36],  # sulfonic acid
                'complexity': 'high'
            },
            'branched_ether_chain': {
                'smiles': 'C(C(F)(F)C(F)(F)F)(C(F)(F)F)OC(C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)F',
                'description': 'Branched perfluoroether',
                'expected_pfasgroups': [31],  # ether
                'complexity': 'high'
            },
            'multi_functional_branched': {
                'smiles': 'C(C(F)(F)C(=O)O)(C(F)(F)C(F)(F)F)OC(F)(F)C(F)(F)S(=O)(=O)O',
                'description': 'Multi-functional branched PFAS (carboxylic acid + ether + sulfonic acid)',
                'expected_pfasgroups': [31, 33, 36],  # ether, carboxylic acid, sulfonic acid
                'complexity': 'very_high'
            },
            'cyclic_branched': {
                'smiles': 'C1(C(F)(F)F)(C(F)(F)F)C(C(F)(F)C(F)(F)F)C(C(F)(F)F)C(C(=O)O)C1(F)F',
                'description': 'Cyclic branched perfluorocarboxylic acid',
                'expected_pfasgroups': [50, 51, 55, 56],  # perfluoroalkyl, polyfluoroalkyl, perfluoro cyclic, polyfluoro cyclic (no aromatics - molecule is aliphatic)
                'complexity': 'high'
            },
            'aromatic_branched': {
                'smiles': 'c1c(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)c(F)c(C(F)(F)C(=O)O)c(F)c1F',
                'description': 'Aromatic branched PFAS with carboxylic acid',
                'expected_pfasgroups': [22, 33, 50, 51, 54, 57, 58],  # side-chain fluorinated aromatics, carboxylic acid, perfluoroalkyl, polyfluoroalkyl, side-chain aromatics, perfluoroaryl, polyfluoroaryl
                'complexity': 'very_high'
            }
        }
        
        all_results = []
        total_tests = 0
        successful_detections = {'pfasgroups': 0, 'atlas': 0}
        
        for test_name, mol_info in test_molecules.items():
            print(f"\n🔬 Testing {test_name}: {mol_info['description']}")
            print(f"   SMILES: {mol_info['smiles']}")
            print(f"   Complexity: {mol_info['complexity']}")
            print(f"   Expected PFASGroups: {mol_info['expected_pfasgroups']}")
            
            test_results = {
                'test_name': test_name,
                'description': mol_info['description'],
                'smiles': mol_info['smiles'],
                'complexity': mol_info['complexity'],
                'expected_pfasgroups': mol_info['expected_pfasgroups'],
                'molecules_tested': 0,
                'pfasgroups_detections': 0,
                'atlas_detections': 0,
                'pfasgroups_correct_detections': 0,
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
                    pfas_result = self.test_with_pfasgroups(mol_info['smiles'], include_PFAS_definitions=True)
                    pfasgroups_detected = pfas_result['success']
                    detected_groups = pfas_result['detected_groups']
                    
                    if pfasgroups_detected:
                        test_results['pfasgroups_detections'] += 1
                        successful_detections['pfasgroups'] += 1
                        
                        # Check if expected groups were detected
                        expected_groups_detected = all(group in detected_groups for group in mol_info['expected_pfasgroups'])
                        if expected_groups_detected:
                            test_results['pfasgroups_correct_detections'] += 1
                    
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
                        'pfasgroups_detected': pfasgroups_detected,
                        'pfasgroups_groups': detected_groups,
                        'pfasgroups_expected_groups': mol_info['expected_pfasgroups'],
                        'pfasgroups_correct': all(group in detected_groups for group in mol_info['expected_pfasgroups']) if pfasgroups_detected else False,
                        'atlas_detected': atlas_detected,
                        'atlas_first_class': atlas_result['first_class'],
                        'atlas_second_class': atlas_result['second_class'],
                        'pfasgroups_execution_time': pfas_result['execution_time'],
                        'atlas_execution_time': atlas_result['execution_time']
                    }
                    
                    test_results['molecules'].append(molecule_result)
                    
                except Exception as e:
                    print(f"   Warning: Failed to test iteration {i+1}: {e}")
                    continue
            
            # Calculate success rates for this test
            pfasgroups_detection_rate = (test_results['pfasgroups_detections'] / max(test_results['molecules_tested'], 1)) * 100
            pfasgroups_accuracy_rate = (test_results['pfasgroups_correct_detections'] / max(test_results['molecules_tested'], 1)) * 100
            atlas_detection_rate = (test_results['atlas_detections'] / max(test_results['molecules_tested'], 1)) * 100
            
            print(f"   ✅ Tested {test_results['molecules_tested']} iterations")
            print(f"   📊 PFASGroups detection: {test_results['pfasgroups_detections']}/{test_results['molecules_tested']} ({pfasgroups_detection_rate:.1f}%)")
            print(f"   🎯 PFASGroups correct groups: {test_results['pfasgroups_correct_detections']}/{test_results['molecules_tested']} ({pfasgroups_accuracy_rate:.1f}%)")
            print(f"   📊 Atlas detection: {test_results['atlas_detections']}/{test_results['molecules_tested']} ({atlas_detection_rate:.1f}%)")
            
            all_results.append(test_results)
        
        # Save complex branched benchmark results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_complex_branched_benchmark_{timestamp}.json"
        
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        # Summary statistics
        overall_pfasgroups_rate = (successful_detections['pfasgroups'] / max(total_tests, 1)) * 100
        overall_atlas_rate = (successful_detections['atlas'] / max(total_tests, 1)) * 100
        
        print(f"\n💾 Complex Branched Benchmark Complete!")
        print(f"📊 Total tests: {total_tests} ({len(test_molecules)} molecule types × {molecules_per_test} iterations)")
        print(f"🎯 Expected result: High detection rates for complex PFAS structures")
        print(f"📈 PFASGroups overall detection rate: {overall_pfasgroups_rate:.1f}% ({successful_detections['pfasgroups']}/{total_tests})")
        print(f"📈 PFAS-Atlas overall detection rate: {overall_atlas_rate:.1f}% ({successful_detections['atlas']}/{total_tests})")
        print(f"💾 Results saved: {output_file}")
        
        if overall_pfasgroups_rate >= 90 and overall_atlas_rate >= 90:
            print(f"🎉 EXCELLENT: Both systems show excellent detection of complex PFAS structures!")
        elif overall_pfasgroups_rate >= 75 and overall_atlas_rate >= 75:
            print(f"✅ GOOD: Both systems show good detection rates for complex structures")
        else:
            print(f"⚠️  CONCERNING: Lower than expected detection rates for complex PFAS structures")
        
        return all_results, output_file
    
    def run_oecd_benchmark(self):
        """Run benchmark on OECD data using groups 1-28"""
        
        print("\n🚀 OECD PFAS BENCHMARK (Groups 1-28)")
        print("=" * 45)
        
        # Load OECD data
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Use relative path: ../../../PFAS-atlas from scripts folder
        oecd_file = os.path.join(script_dir, '..', '..', '..', 'PFAS-atlas', 'input_data', 'OECD_4000', 'step3_OECD_Class_0812.csv')
        oecd_file = os.path.normpath(oecd_file)  # Normalize the path
        try:
            import pandas as pd
            oecd_data = pd.read_csv(oecd_file)
            print(f"📊 Loaded {len(oecd_data)} OECD molecules")
        except Exception as e:
            print(f"❌ Error loading OECD data: {e}")
            print(f"   Expected path: {oecd_file}")
            import sys
            sys.exit(1)
            
        all_results = []
        
        # Process OECD molecules
        print(f"🧪 Testing {len(oecd_data)} OECD molecules with PFASGroups (groups 1-28)")
        
        for idx, row in oecd_data.iterrows():
            smiles = row['SMILES']
            first_class = row.get('First_Class', 'Unknown')
            second_class = row.get('Second_Class', 'Unknown')
            
            if pd.isna(smiles) or smiles.strip() == '':
                continue
                
            # Test with PFASGroups - accuracy test
            pfas_result = self.test_with_pfasgroups(smiles, include_PFAS_definitions=True)
            
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
                print(f"  ✅ Processed {idx + 1}/{len(oecd_data)} molecules")
        
        # Save OECD results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"data/pfas_oecd_benchmark_{timestamp}.json"
        
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
    print("  4. Non-Fluorinated Exclusion Benchmark (false positive test)")
    print("  5. Complex Branched PFAS Benchmark (positive control validation)")
    print("  6. All benchmarks")
    
    choice = input("\nChoose benchmark (1-6) or press Enter for all: ").strip()
    
    if choice in ['1', '6', '']:
        print("\nRunning Enhanced Functional Groups Benchmark...")
        enhanced_results, enhanced_file = benchmark.run_enhanced_benchmark()
    else:
        enhanced_results, enhanced_file = None, None
    
    if choice in ['2', '6', '']:
        print("\nRunning OECD Benchmark...")
        oecd_results, oecd_file = benchmark.run_oecd_benchmark()
    else:
        oecd_results, oecd_file = None, None
    
    if choice in ['3', '6', '']:
        print("\nRunning Timing Performance Benchmark...")
        timing_results, timing_file = benchmark.run_timing_benchmark(200, 10)  # 200 molecules, 10 iterations
    else:
        timing_results, timing_file = None, None
    
    if choice in ['4', '6', '']:
        print("\nRunning Non-Fluorinated Exclusion Benchmark...")
        nonfluor_results, nonfluor_file = benchmark.run_non_fluorinated_benchmark(50)
    else:
        nonfluor_results, nonfluor_file = None, None
    
    if choice in ['5', '6', '']:
        print("\nRunning Complex Branched PFAS Benchmark...")
        complex_results, complex_file = benchmark.run_complex_branched_benchmark(50)
    else:
        complex_results, complex_file = None, None
    
    print(f"\n🎯 Benchmark Complete! Next steps:")
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
