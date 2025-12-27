#!/usr/bin/env python3
"""
Unified PFAS Benchmark Report Generator

This script combines all benchmark results into a single comprehensive HTML report.
It automatically finds the most recent benchmark files and generates a unified analysis.
"""

import json
import sys
import os
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
import base64
from io import BytesIO
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import shutil
from collections import defaultdict

# Import the enhanced Sankey diagrams from enhanced_analysis
try:
    from enhanced_analysis import create_enhanced_sankey_comparison, analyze_system_comparison
except ImportError:
    def create_enhanced_sankey_comparison(*args):
        return None, None, None
    def analyze_system_comparison(results):
        return {}, {}

def find_latest_benchmark_files():
    """Find the most recent benchmark result files."""
    files = {
        'enhanced': None,
        'oecd': None, 
        'timing': None,
        'nonfluorinated': None,
        'complex': None
    }
    
    patterns = {
        'enhanced': 'pfas_enhanced_benchmark_*.json',
        'oecd': 'pfas_oecd_benchmark_*.json',
        'timing': 'pfas_timing_benchmark_*.json', 
        'nonfluorinated': 'pfas_non_fluorinated_benchmark_*.json',
        'complex': 'pfas_complex_branched_benchmark_*.json'
    }
    
    # Check both current directory and data subdirectory
    search_paths = ['.', 'data']
    
    for benchmark_type, pattern in patterns.items():
        try:
            matching_files = []
            for search_path in search_paths:
                matching_files.extend(list(Path(search_path).glob(pattern)))
            
            if matching_files:
                # Sort by modification time, get most recent
                files[benchmark_type] = max(matching_files, key=lambda x: x.stat().st_mtime)
        except Exception as e:
            print(f"Warning: Could not find {pattern}: {e}")
    
    return files

def load_benchmark_data(filename):
    """Load benchmark data from JSON file."""
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"Warning: Could not load {filename}: {e}")
        return None

def analyze_enhanced_benchmark(data):
    """Analyze enhanced benchmark results with proper data structure handling."""
    if not data:
        return None
    
    # Load PFAS groups definitions for lookup
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            pfas_groups_list = json.load(f)
        pfas_groups_lookup = {group['id']: group for group in pfas_groups_list}
    except:
        pfas_groups_lookup = {}
    
    # Load specificity connections for OECD analysis
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/tests/specificity_test_groups.json', 'r') as f:
            specificity_groups = json.load(f)
    except:
        specificity_groups = []
    
    # Separate different types of molecules
    single_group_molecules = [mol for mol in data if mol.get('molecule_data', {}).get('generation_type') == 'single_group']
    multi_group_molecules = [mol for mol in data if mol.get('molecule_data', {}).get('generation_type') == 'multi_group']
    oecd_group_molecules = [mol for mol in data if mol.get('molecule_data', {}).get('generation_type') == 'oecd_group']
    
    # Analyze single-group molecules
    single_group_analysis = {}
    for mol in single_group_molecules:
        mol_data = mol.get('molecule_data', {})
        group_id = mol_data.get('group_id')
        if group_id is not None:
            if group_id not in single_group_analysis:
                group_name = mol_data.get('group_name', 'Unknown')
                single_group_analysis[group_id] = {
                    'id': group_id,
                    'name': group_name,
                    'molecules': 0,
                    'pfasgroups_detected': 0,
                    'atlas_detected': 0,
                    'pfasgroups_correct': 0,
                    'atlas_classes': {}
                }
            
            single_group_analysis[group_id]['molecules'] += 1
            
            # PFASGroups analysis
            pfas_result = mol.get('pfasgroups_result', {})
            if pfas_result.get('success', False):
                single_group_analysis[group_id]['pfasgroups_detected'] += 1
                detected_groups = pfas_result.get('detected_groups', [])
                if group_id in detected_groups:
                    single_group_analysis[group_id]['pfasgroups_correct'] += 1
            
            # Atlas analysis
            atlas_result = mol.get('atlas_result', {})
            if atlas_result.get('success', False):
                single_group_analysis[group_id]['atlas_detected'] += 1
                second_class = atlas_result.get('second_class', atlas_result.get('first_class', 'Unknown'))
                if second_class in single_group_analysis[group_id]['atlas_classes']:
                    single_group_analysis[group_id]['atlas_classes'][second_class] += 1
                else:
                    single_group_analysis[group_id]['atlas_classes'][second_class] = 1
    
    # Analyze multi-group molecules
    multi_group_analysis = {
        'total_molecules': len(multi_group_molecules),
        'pfasgroups_detected': 0,
        'atlas_detected': 0,
        'pfasgroups_all_correct': 0,
        'pfasgroups_partial_correct': 0,
        'combinations': {}
    }
    
    for mol in multi_group_molecules:
        mol_data = mol.get('molecule_data', {})
        target_groups = mol_data.get('target_groups', [])
        target_key = tuple(sorted(target_groups))
        
        if target_key not in multi_group_analysis['combinations']:
            multi_group_analysis['combinations'][target_key] = {
                'target_groups': target_groups,
                'molecules': 0,
                'pfasgroups_detected': 0,
                'atlas_detected': 0,
                'pfasgroups_all_correct': 0,
                'pfasgroups_partial_correct': 0
            }
        
        multi_group_analysis['combinations'][target_key]['molecules'] += 1
        
        # PFASGroups analysis
        pfas_result = mol.get('pfasgroups_result', {})
        if pfas_result.get('success', False):
            multi_group_analysis['pfasgroups_detected'] += 1
            multi_group_analysis['combinations'][target_key]['pfasgroups_detected'] += 1
            
            detected_groups = pfas_result.get('detected_groups', [])
            detected_targets = [g for g in detected_groups if g in target_groups]
            
            if len(detected_targets) == len(target_groups):
                multi_group_analysis['pfasgroups_all_correct'] += 1
                multi_group_analysis['combinations'][target_key]['pfasgroups_all_correct'] += 1
            elif len(detected_targets) > 0:
                multi_group_analysis['pfasgroups_partial_correct'] += 1
                multi_group_analysis['combinations'][target_key]['pfasgroups_partial_correct'] += 1
        
        # Atlas analysis
        atlas_result = mol.get('atlas_result', {})
        if atlas_result.get('success', False):
            multi_group_analysis['atlas_detected'] += 1
            multi_group_analysis['combinations'][target_key]['atlas_detected'] += 1
    
    # Analyze OECD group molecules
    oecd_group_analysis = {}
    for mol in oecd_group_molecules:
        mol_data = mol.get('molecule_data', {})
        group_id = mol_data.get('group_id')
        if group_id is not None:
            if group_id not in oecd_group_analysis:
                group_name = mol_data.get('group_name', 'Unknown')
                oecd_group_analysis[group_id] = {
                    'id': group_id,
                    'name': group_name,
                    'molecules': 0,
                    'pfasgroups_detected': 0,
                    'atlas_detected': 0,
                    'pfasgroups_groups': {},
                    'atlas_classes': {},
                    'correspondence_pairs': []
                }
            
            oecd_group_analysis[group_id]['molecules'] += 1
            
            # PFASGroups analysis
            pfas_result = mol.get('pfasgroups_result', {})
            detected_groups = pfas_result.get('detected_groups', [])
            if pfas_result.get('success', False) and detected_groups:
                oecd_group_analysis[group_id]['pfasgroups_detected'] += 1
                for detected_group in detected_groups:
                    if detected_group in oecd_group_analysis[group_id]['pfasgroups_groups']:
                        oecd_group_analysis[group_id]['pfasgroups_groups'][detected_group] += 1
                    else:
                        oecd_group_analysis[group_id]['pfasgroups_groups'][detected_group] = 1
            
            # Atlas analysis
            atlas_result = mol.get('atlas_result', {})
            if atlas_result.get('success', False):
                oecd_group_analysis[group_id]['atlas_detected'] += 1
                second_class = atlas_result.get('second_class', 'Unknown')
                if second_class in oecd_group_analysis[group_id]['atlas_classes']:
                    oecd_group_analysis[group_id]['atlas_classes'][second_class] += 1
                else:
                    oecd_group_analysis[group_id]['atlas_classes'][second_class] = 1
                
                # Create correspondence pairs for Sankey diagram
                for detected_group in detected_groups:
                    if detected_group != group_id:  # Only cross-group correspondences
                        oecd_group_analysis[group_id]['correspondence_pairs'].append({
                            'source_oecd': group_id,
                            'target_pfas': detected_group,
                            'atlas_class': second_class
                        })
    
    # Overall summary including OECD
    total_molecules = len(data)
    total_pfasgroups_detected = len([mol for mol in data if mol.get('pfasgroups_result', {}).get('success', False)])
    total_atlas_detected = len([mol for mol in data if mol.get('atlas_result', {}).get('success', False)])
    
    summary = {
        'total_molecules': total_molecules,
        'total_single_group': len(single_group_molecules),
        'total_multi_group': len(multi_group_molecules),
        'total_oecd_group': len(oecd_group_molecules),
        'pfasgroups_detected': total_pfasgroups_detected,
        'atlas_detected': total_atlas_detected,
        'pfasgroups_rate': (total_pfasgroups_detected / max(total_molecules, 1)) * 100,
        'atlas_rate': (total_atlas_detected / max(total_molecules, 1)) * 100,
        'single_group_analysis': single_group_analysis,
        'multi_group_analysis': multi_group_analysis,
        'oecd_group_analysis': oecd_group_analysis,
        'pfas_groups_lookup': pfas_groups_lookup,
        'specificity_groups': specificity_groups
    }
    
    return summary

def load_oecd_csv_data():
    """Load OECD CSV data directly for robustness analysis."""
    oecd_csv_path = '/home/luc/git/PFAS-atlas/input_data/OECD_4000/step3_OECD_Class_0812.csv'
    
    try:
        import pandas as pd
        df = pd.read_csv(oecd_csv_path)
        
        # Convert to list of dictionaries for consistency
        molecules = []
        for _, row in df.iterrows():
            molecules.append({
                'rdkit_smiles': row['RDKIT_SMILES'],
                'original_smiles': row['SMILES'], 
                'first_class': row['First_Class'],
                'second_class': row['Second_Class'],
                'type': row['type']
            })
        
        print(f"Loaded {len(molecules)} molecules from OECD CSV")
        return molecules
    except Exception as e:
        print(f"Error loading OECD CSV: {e}")
        return None

def analyze_oecd_robustness(oecd_molecules):
    """Analyze robustness of PFAS-Atlas on OECD dataset."""
    if not oecd_molecules:
        return None
        
    try:
        # Import required modules
        sys.path.append('/home/luc/git/PFAS-atlas')
        from classification_helper import classify_pfas_molecule
        
        sys.path.append('/home/luc/git/PFASGroups')
        from PFASgroups import parse_smiles
        
        total_molecules = len(oecd_molecules)
        atlas_agreements = 0
        atlas_mismatches = []
        pfasgroups_detections = 0
        pfasgroups_mismatches = []
        pfasgroups_errors = 0
        
        # Track OECD group (type < 29) correspondence for Sankey
        oecd_groups_correspondence = defaultdict(lambda: defaultdict(int))
        
        print(f"Running robustness analysis on {total_molecules} OECD molecules...")
        
        for i, mol in enumerate(oecd_molecules):
            if i % 500 == 0:
                print(f"Progress: {i}/{total_molecules}")
            
            smiles = mol['rdkit_smiles']
            expected_first = mol['first_class']
            expected_second = mol['second_class']
            expected_type = mol['type']
            
            # Run PFAS-Atlas predictions
            try:
                atlas_result = classify_pfas_molecule(smiles)
                atlas_first = atlas_result[0] if len(atlas_result) > 0 else 'Not detected'
                atlas_second = atlas_result[1] if len(atlas_result) > 1 else 'Not detected'
                
                # Check Atlas robustness (should match exactly)
                if atlas_first == expected_first and atlas_second == expected_second:
                    atlas_agreements += 1
                else:
                    atlas_mismatches.append({
                        'index': i,
                        'smiles': smiles,
                        'expected_first': expected_first,
                        'expected_second': expected_second,
                        'predicted_first': atlas_first,
                        'predicted_second': atlas_second
                    })
                    
            except Exception as e:
                atlas_mismatches.append({
                    'index': i,
                    'smiles': smiles,
                    'expected_first': expected_first,
                    'expected_second': expected_second,
                    'predicted_first': 'Error',
                    'predicted_second': 'Error',
                    'error': str(e)
                })
            
            # Run PFASGroups analysis
            try:
                from PFASgroups.core import parse_groups_in_mol
                from rdkit.Chem.rdMolDescriptors import CalcMolFormula
                from rdkit import Chem
                
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Could not parse SMILES: {smiles}")
                
                formula = CalcMolFormula(mol)
                pfas_matches = parse_groups_in_mol(mol, formula)
                
                if pfas_matches:
                    pfasgroups_detections += 1
                    # Extract group IDs from the matches
                    detected_groups = [match[0].id for match in pfas_matches]
                    
                    # For OECD correspondence, only consider groups < 29 (OECD groups)
                    oecd_detected_groups = [g for g in detected_groups if g < 29]
                    
                    # Map expected type to detected groups for Sankey
                    if oecd_detected_groups:
                        for group in oecd_detected_groups:
                            oecd_groups_correspondence[expected_type][group] += 1
                    else:
                        oecd_groups_correspondence[expected_type]['no_detection'] += 1
                else:
                    # No PFAS detection
                    pfasgroups_mismatches.append({
                        'index': i,
                        'smiles': smiles,
                        'expected_type': expected_type,
                        'expected_first': expected_first,
                        'expected_second': expected_second,
                        'detected': False
                    })
                    oecd_groups_correspondence[expected_type]['no_detection'] += 1
                    
            except Exception as e:
                pfasgroups_mismatches.append({
                    'index': i,
                    'smiles': smiles,
                    'expected_type': expected_type,
                    'expected_first': expected_first,
                    'expected_second': expected_second,
                    'detected': False,
                    'error': str(e)[:100]  # Truncate error message
                })
                oecd_groups_correspondence[expected_type]['no_detection'] += 1
        
        atlas_accuracy = (atlas_agreements / total_molecules) * 100 if total_molecules > 0 else 0
        pfasgroups_detection_rate = (pfasgroups_detections / total_molecules) * 100 if total_molecules > 0 else 0
        
        print(f"Analysis complete: {atlas_agreements}/{total_molecules} Atlas agreements, {pfasgroups_detections} PFASGroups detections")
        
        return {
            'total_molecules': total_molecules,
            'atlas_agreements': atlas_agreements,
            'atlas_accuracy': atlas_accuracy,
            'atlas_mismatches': atlas_mismatches[:100],  # Limit for display
            'pfasgroups_detections': pfasgroups_detections,
            'pfasgroups_detection_rate': pfasgroups_detection_rate,
            'pfasgroups_mismatches': pfasgroups_mismatches[:100],  # Limit for display
            'oecd_groups_correspondence': dict(oecd_groups_correspondence)
        }
        
    except Exception as e:
        print(f"Error in robustness analysis: {e}")
        return None

def create_oecd_sankey_diagram(oecd_correspondence):
    """Create Sankey diagram showing OECD groups vs PFASGroups detections."""
    if not oecd_correspondence:
        return None
        
    try:
        # Load PFAS groups for names
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            pfas_groups_list = json.load(f)
        pfas_groups_lookup = {group['id']: group['name'] for group in pfas_groups_list}
        
        # Build Sankey data
        source_nodes = []  # OECD types
        target_nodes = []  # PFASGroups detections
        links = []
        
        all_oecd_types = list(oecd_correspondence.keys())
        all_detected_groups = set()
        
        # Collect all detected groups
        for oecd_type, detections in oecd_correspondence.items():
            for group_id in detections.keys():
                if group_id != 'no_detection' and isinstance(group_id, int) and group_id < 29:
                    all_detected_groups.add(group_id)
        
        # Add 'no_detection' as a special case
        all_detected_groups.add('no_detection')
        
        # Create node lists
        source_nodes = [f"OECD: {oecd_type}" for oecd_type in all_oecd_types]
        target_nodes = []
        
        # Separate and sort numeric group IDs and special cases
        numeric_groups = [g for g in all_detected_groups if isinstance(g, int)]
        special_groups = [g for g in all_detected_groups if not isinstance(g, int)]
        
        sorted_groups = sorted(numeric_groups) + sorted(special_groups)
        
        for group_id in sorted_groups:
            if group_id == 'no_detection':
                target_nodes.append("No Detection")
            else:
                group_name = pfas_groups_lookup.get(group_id, f"Group {group_id}")
                target_nodes.append(f"PFASGroups: {group_id} - {group_name}")
        
        all_nodes = source_nodes + target_nodes
        
        # Create links
        for oecd_idx, oecd_type in enumerate(all_oecd_types):
            detections = oecd_correspondence[oecd_type]
            
            for group_id, count in detections.items():
                if count > 0:
                    if group_id == 'no_detection':
                        target_idx = len(source_nodes) + sorted_groups.index('no_detection')
                    else:
                        target_idx = len(source_nodes) + sorted_groups.index(group_id)
                    
                    links.append({
                        'source': oecd_idx,
                        'target': target_idx,
                        'value': count
                    })
        
        # Create Plotly Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_nodes,
                color=["rgba(31, 119, 180, 0.8)" for _ in source_nodes] + 
                      ["rgba(44, 160, 44, 0.8)" for _ in target_nodes]
            ),
            link=dict(
                source=[link['source'] for link in links],
                target=[link['target'] for link in links],
                value=[link['value'] for link in links],
                color=["rgba(44, 160, 44, 0.3)" for _ in links]
            )
        )])
        
        fig.update_layout(
            title_text="OECD Classification vs PFASGroups Detection Flow<br>",
            font_size=12,
            height=600,
            width=1200
        )
        
        return fig.to_html(include_plotlyjs='cdn', div_id="sankey-oecd-correspondence")
        
    except Exception as e:
        print(f"Error creating OECD Sankey diagram: {e}")
        return None

def analyze_oecd_benchmark(data):
    """Analyze OECD benchmark results with detailed correspondence and misclassification analysis."""
    if not data:
        return None
    
    # Load PFAS groups lookup
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            pfas_groups_list = json.load(f)
        pfas_groups_lookup = {group['id']: group for group in pfas_groups_list}
    except:
        pfas_groups_lookup = {}
    
    total_molecules = len(data)
    pfasgroups_detected = sum(1 for mol in data if mol.get('pfasgroups_result', {}).get('success', False))
    atlas_detected = sum(1 for mol in data if mol.get('atlas_result', {}).get('success', False))
    
    # Detailed correspondence analysis
    oecd_to_pfasgroups = {}  # OECD class -> PFASGroups detection stats
    oecd_to_atlas = {}       # OECD class -> Atlas classification agreement
    
    # Track misclassifications
    misclassified_molecules = []
    
    for i, mol in enumerate(data):
        mol_data = mol.get('molecule_data', {})
        pfas_result = mol.get('pfasgroups_result', {})
        atlas_result = mol.get('atlas_result', {})
        
        oecd_first = mol_data.get('oecd_first_class', 'Unknown')
        oecd_second = mol_data.get('oecd_second_class', 'Unknown')
        atlas_first = atlas_result.get('first_class', 'Not detected')
        atlas_second = atlas_result.get('second_class', 'Not detected')
        
        # Initialize tracking for this OECD class
        if oecd_first not in oecd_to_pfasgroups:
            oecd_to_pfasgroups[oecd_first] = {
                'total': 0, 'pfasgroups_detected': 0, 
                'groups_detected': {}, 'no_detection': 0
            }
        
        if oecd_first not in oecd_to_atlas:
            oecd_to_atlas[oecd_first] = {
                'total': 0, 'atlas_agreement': 0, 'atlas_disagreement': 0,
                'atlas_classifications': {}
            }
        
        oecd_to_pfasgroups[oecd_first]['total'] += 1
        oecd_to_atlas[oecd_first]['total'] += 1
        
        # PFASGroups analysis
        if pfas_result.get('success', False):
            oecd_to_pfasgroups[oecd_first]['pfasgroups_detected'] += 1
            detected_groups = pfas_result.get('detected_groups', [])
            
            for group in detected_groups:
                group_name = pfas_groups_lookup.get(group, {}).get('name', f'Group {group}')
                if group_name in oecd_to_pfasgroups[oecd_first]['groups_detected']:
                    oecd_to_pfasgroups[oecd_first]['groups_detected'][group_name] += 1
                else:
                    oecd_to_pfasgroups[oecd_first]['groups_detected'][group_name] = 1
        else:
            oecd_to_pfasgroups[oecd_first]['no_detection'] += 1
        
        # Atlas classification correspondence
        if atlas_result.get('success', False):
            # Check if Atlas agrees with OECD classification
            if atlas_first == oecd_first or atlas_second == oecd_second:
                oecd_to_atlas[oecd_first]['atlas_agreement'] += 1
            else:
                oecd_to_atlas[oecd_first]['atlas_disagreement'] += 1
                
                # Track misclassification
                misclassified_molecules.append({
                    'index': i,
                    'smiles': mol_data.get('smiles', 'N/A'),
                    'oecd_first': oecd_first,
                    'oecd_second': oecd_second,
                    'atlas_first': atlas_first,
                    'atlas_second': atlas_second,
                    'pfasgroups_detected': pfas_result.get('detected_groups', [])
                })
            
            # Track Atlas classification distribution (prefer second_class for detail)
            primary_atlas_class = atlas_second if atlas_second != 'Not detected' else atlas_first
            if primary_atlas_class in oecd_to_atlas[oecd_first]['atlas_classifications']:
                oecd_to_atlas[oecd_first]['atlas_classifications'][primary_atlas_class] += 1
            else:
                oecd_to_atlas[oecd_first]['atlas_classifications'][primary_atlas_class] = 1
        else:
            oecd_to_atlas[oecd_first]['atlas_disagreement'] += 1
            misclassified_molecules.append({
                'index': i,
                'smiles': mol_data.get('smiles', 'N/A'),
                'oecd_first': oecd_first,
                'oecd_second': oecd_second,
                'atlas_first': 'Not detected',
                'atlas_second': 'Not detected',
                'pfasgroups_detected': pfas_result.get('detected_groups', [])
            })
    
    # Calculate summary statistics
    total_agreement = sum(oecd_to_atlas[oecd_class]['atlas_agreement'] for oecd_class in oecd_to_atlas)
    total_disagreement = sum(oecd_to_atlas[oecd_class]['atlas_disagreement'] for oecd_class in oecd_to_atlas)
    
    return {
        'total_molecules': total_molecules,
        'pfasgroups_detected': pfasgroups_detected,
        'atlas_detected': atlas_detected,
        'pfasgroups_rate': (pfasgroups_detected / max(total_molecules, 1)) * 100,
        'atlas_rate': (atlas_detected / max(total_molecules, 1)) * 100,
        'oecd_to_pfasgroups': oecd_to_pfasgroups,
        'oecd_to_atlas': oecd_to_atlas,
        'atlas_agreement': total_agreement,
        'atlas_disagreement': total_disagreement,
        'atlas_agreement_rate': (total_agreement / max(total_molecules, 1)) * 100,
        'misclassified_molecules': misclassified_molecules[:50],  # Limit to first 50 for display
        'pfas_groups_lookup': pfas_groups_lookup
    }

def create_time_performance_plot(data):
    """Create time performance plot across carbon counts and molecule sizes."""
    if not data:
        return None
    
    plt.style.use('seaborn-v0_8')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Time Performance Analysis Across Molecule Sizes', fontsize=16, fontweight='bold')
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(data)
    
    # Plot 1: Execution time vs number of atoms
    ax1.scatter(df['num_atoms'], df['pfasgroups_time_avg'] * 1000, alpha=0.6, label='PFASGroups', color='blue', s=20)
    ax1.scatter(df['num_atoms'], df['atlas_time_avg'] * 1000, alpha=0.6, label='PFAS-Atlas', color='red', s=20)
    ax1.set_xlabel('Number of Atoms')
    ax1.set_ylabel('Average Time (ms)')
    ax1.set_title('Execution Time vs Molecule Size')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Performance ratio vs chain length
    ratio = (df['atlas_time_avg'] / df['pfasgroups_time_avg']).fillna(1)
    ax2.scatter(df['chain_length'], ratio, alpha=0.6, color='green', s=20)
    ax2.set_xlabel('Chain Length (Carbon Count)')
    ax2.set_ylabel('Performance Ratio (Atlas/PFASGroups)')
    ax2.set_title('Performance Ratio vs Chain Length')
    ax2.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Equal Performance')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Box plot by size bins
    df['size_bin'] = pd.cut(df['num_atoms'], bins=6, labels=['XS', 'S', 'M', 'L', 'XL', 'XXL'])
    
    pfas_times_binned = [df[df['size_bin'] == bin]['pfasgroups_time_avg'].values * 1000 for bin in df['size_bin'].cat.categories]
    atlas_times_binned = [df[df['size_bin'] == bin]['atlas_time_avg'].values * 1000 for bin in df['size_bin'].cat.categories]
    
    x_pos = np.arange(len(df['size_bin'].cat.categories))
    bp1 = ax3.boxplot(pfas_times_binned, positions=x_pos-0.2, widths=0.3, patch_artist=True, 
                      boxprops=dict(facecolor='lightblue'), labels=df['size_bin'].cat.categories)
    bp2 = ax3.boxplot(atlas_times_binned, positions=x_pos+0.2, widths=0.3, patch_artist=True, 
                      boxprops=dict(facecolor='lightcoral'), labels=df['size_bin'].cat.categories)
    
    ax3.set_xlabel('Molecule Size Bins')
    ax3.set_ylabel('Execution Time (ms)')
    ax3.set_title('Time Distribution by Size Category')
    ax3.legend([bp1["boxes"][0], bp2["boxes"][0]], ['PFASGroups', 'PFAS-Atlas'])
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Scaling analysis
    # Group by chain length and calculate median times
    chain_stats = df.groupby('chain_length').agg({
        'pfasgroups_time_avg': ['median', 'std'],
        'atlas_time_avg': ['median', 'std'],
        'num_atoms': 'median'
    }).round(4)
    
    chain_lengths = chain_stats.index
    pfas_medians = chain_stats[('pfasgroups_time_avg', 'median')] * 1000
    pfas_stds = chain_stats[('pfasgroups_time_avg', 'std')] * 1000
    atlas_medians = chain_stats[('atlas_time_avg', 'median')] * 1000
    atlas_stds = chain_stats[('atlas_time_avg', 'std')] * 1000
    
    ax4.errorbar(chain_lengths, pfas_medians, yerr=pfas_stds, label='PFASGroups', 
                 marker='o', capsize=5, capthick=2, color='blue')
    ax4.errorbar(chain_lengths, atlas_medians, yerr=atlas_stds, label='PFAS-Atlas', 
                 marker='s', capsize=5, capthick=2, color='red')
    ax4.set_xlabel('Chain Length')
    ax4.set_ylabel('Median Time ± Std (ms)')
    ax4.set_title('Scaling Performance with Chain Length')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save to base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    buffer.seek(0)
    plot_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()
    
    return plot_base64

def create_sankey_diagram(enhanced_data):
    """Create Sankey diagram showing classification correspondence between PFASGroups and PFAS-Atlas."""
    if not enhanced_data:
        return None
    
    try:
        import plotly.graph_objects as go
        import plotly.io as pio
        pio.kaleido.scope.mathjax = None
    except ImportError:
        print("Warning: Plotly not available for Sankey diagram - creating alternative visualization")
        return create_sankey_alternative_plot(enhanced_data)
    
    # Load PFAS groups lookup
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            pfas_groups_list = json.load(f)
        pfas_groups_lookup = {group['id']: group['name'] for group in pfas_groups_list}
    except:
        pfas_groups_lookup = {}
    
    # Focus on single-group molecules for cleaner Sankey
    single_group_molecules = [mol for mol in enhanced_data if mol.get('molecule_data', {}).get('generation_type') == 'single_group']
    
    # Also include OECD group molecules for correspondence analysis
    oecd_group_molecules = [mol for mol in enhanced_data if mol.get('molecule_data', {}).get('generation_type') == 'oecd_group']
    
    # Create correspondence matrix
    pfasgroups_to_atlas = {}
    
    # Process enhanced single-group molecules
    for mol in single_group_molecules:
        mol_data = mol.get('molecule_data', {})
        target_group = mol_data.get('group_id')
        
        pfas_result = mol.get('pfasgroups_result', {})
        atlas_result = mol.get('atlas_result', {})
        
        if pfas_result.get('success') and atlas_result.get('success'):
            detected_groups = pfas_result.get('detected_groups', [])
            # Use second_class for more detailed classification, fallback to first_class
            atlas_class = atlas_result.get('second_class', atlas_result.get('first_class', 'Unknown'))
            
            # For each detected group, record the atlas classification
            for group in detected_groups:
                group_name = pfas_groups_lookup.get(group, f"Group {group}")
                key = (group_name, atlas_class)
                if key in pfasgroups_to_atlas:
                    pfasgroups_to_atlas[key] += 1
                else:
                    pfasgroups_to_atlas[key] = 1
    
    # Process OECD group molecules
    for mol in oecd_group_molecules:
        mol_data = mol.get('molecule_data', {})
        target_group = mol_data.get('group_id')  # This is the OECD group ID (1-28)
        
        pfas_result = mol.get('pfasgroups_result', {})
        atlas_result = mol.get('atlas_result', {})
        
        if pfas_result.get('success') and atlas_result.get('success'):
            detected_groups = pfas_result.get('detected_groups', [])
            atlas_class = atlas_result.get('second_class', atlas_result.get('first_class', 'Unknown'))
            
            # Record correspondence between OECD group and detected PFASGroups/Atlas
            oecd_group_name = pfas_groups_lookup.get(target_group, f"OECD Group {target_group}")
            
            # Add OECD -> Atlas correspondence
            key = (oecd_group_name, atlas_class)
            if key in pfasgroups_to_atlas:
                pfasgroups_to_atlas[key] += 1
            else:
                pfasgroups_to_atlas[key] = 1
    
    if not pfasgroups_to_atlas:
        return None
    
    # Create nodes and links for Sankey
    pfas_groups = set()
    atlas_classes = set()
    
    for (group, atlas_class), count in pfasgroups_to_atlas.items():
        pfas_groups.add(group)
        atlas_classes.add(atlas_class)
    
    # Create node lists
    pfas_groups = sorted(list(pfas_groups))
    atlas_classes = sorted(list(atlas_classes))
    
    all_nodes = pfas_groups + atlas_classes
    node_indices = {node: i for i, node in enumerate(all_nodes)}
    
    # Create links
    source = []
    target = []
    value = []
    
    for (group, atlas_class), count in pfasgroups_to_atlas.items():
        source.append(node_indices[group])
        target.append(node_indices[atlas_class])
        value.append(count)
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=["lightblue" if node in pfas_groups else "lightcoral" for node in all_nodes]
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color="rgba(255, 255, 255, 0.4)"
        )
    )])
    
    fig.update_layout(
        title_text="Classification Correspondence: PFASGroups → PFAS-Atlas",
        title_x=0.5,
        font_size=12,
        width=1200,
        height=800
    )
    
    # Convert to base64
    img_bytes = fig.to_image(format="png", width=1200, height=800)
    sankey_base64 = base64.b64encode(img_bytes).decode()
    
    return sankey_base64

def create_sankey_alternative_plot(enhanced_data):
    """Create alternative visualization when Plotly is not available."""
    if not enhanced_data:
        return None
    
    # Load PFAS groups lookup
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            pfas_groups_list = json.load(f)
        pfas_groups_lookup = {group['id']: group['name'] for group in pfas_groups_list}
    except:
        pfas_groups_lookup = {}
    
    # Focus on single-group molecules
    single_group_molecules = [mol for mol in enhanced_data if mol.get('molecule_data', {}).get('generation_type') == 'single_group']
    
    # Create correspondence matrix
    correspondence_matrix = {}
    atlas_classes = set()
    pfas_groups = set()
    
    for mol in single_group_molecules:
        mol_data = mol.get('molecule_data', {})
        target_group = mol_data.get('group_id')
        
        pfas_result = mol.get('pfasgroups_result', {})
        atlas_result = mol.get('atlas_result', {})
        
        if pfas_result.get('success') and atlas_result.get('success'):
            detected_groups = pfas_result.get('detected_groups', [])
            # Use second_class for more detailed classification, fallback to first_class
            atlas_class = atlas_result.get('second_class', atlas_result.get('first_class', 'Unknown'))
            atlas_classes.add(atlas_class)
            
            for group in detected_groups:
                group_name = pfas_groups_lookup.get(group, f"Group {group}")
                pfas_groups.add(group_name)
                
                if group_name not in correspondence_matrix:
                    correspondence_matrix[group_name] = {}
                if atlas_class not in correspondence_matrix[group_name]:
                    correspondence_matrix[group_name][atlas_class] = 0
                correspondence_matrix[group_name][atlas_class] += 1
    
    # Create heatmap visualization
    plt.style.use('seaborn-v0_8')
    pfas_groups_list = sorted(list(pfas_groups))
    atlas_classes_list = sorted(list(atlas_classes))
    
    # Create matrix for heatmap
    matrix = np.zeros((len(pfas_groups_list), len(atlas_classes_list)))
    for i, pfas_group in enumerate(pfas_groups_list):
        for j, atlas_class in enumerate(atlas_classes_list):
            matrix[i, j] = correspondence_matrix.get(pfas_group, {}).get(atlas_class, 0)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(matrix, cmap='Blues', aspect='auto')
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(atlas_classes_list)))
    ax.set_yticks(np.arange(len(pfas_groups_list)))
    ax.set_xticklabels(atlas_classes_list, rotation=45, ha='right')
    ax.set_yticklabels(pfas_groups_list)
    
    # Add text annotations
    for i in range(len(pfas_groups_list)):
        for j in range(len(atlas_classes_list)):
            if matrix[i, j] > 0:
                text = ax.text(j, i, int(matrix[i, j]), ha="center", va="center", 
                              color="white" if matrix[i, j] > matrix.max()/2 else "black")
    
    ax.set_title("Classification Correspondence: PFASGroups → PFAS-Atlas\\n(Heatmap Alternative)")
    ax.set_xlabel("PFAS-Atlas Classification")
    ax.set_ylabel("PFASGroups Detection")
    
    # Add colorbar
    plt.colorbar(im, ax=ax, label='Number of Molecules')
    
    plt.tight_layout()
    
    # Save to base64
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
    buffer.seek(0)
    plot_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()
    
    return plot_base64

def create_interactive_time_performance_plot(data):
    """Create interactive time performance plot using Plotly."""
    if not data:
        return None
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Execution Time vs Molecule Size', 'Performance Ratio vs Chain Length',
                       'Time Distribution by Size Category', 'Scaling Performance with Chain Length'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Plot 1: Execution time vs number of atoms
    fig.add_trace(
        go.Scatter(x=df['num_atoms'], y=df['pfasgroups_time_avg'] * 1000,
                  mode='markers', name='PFASGroups', opacity=0.7,
                  marker=dict(color='blue', size=6),
                  hovertemplate='Atoms: %{x}<br>Time: %{y:.2f} ms<extra></extra>'),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(x=df['num_atoms'], y=df['atlas_time_avg'] * 1000,
                  mode='markers', name='PFAS-Atlas', opacity=0.7,
                  marker=dict(color='red', size=6),
                  hovertemplate='Atoms: %{x}<br>Time: %{y:.2f} ms<extra></extra>'),
        row=1, col=1
    )
    
    # Plot 2: Performance ratio vs chain length
    ratio = (df['atlas_time_avg'] / df['pfasgroups_time_avg']).fillna(1)
    fig.add_trace(
        go.Scatter(x=df['chain_length'], y=ratio,
                  mode='markers', name='Performance Ratio', 
                  marker=dict(color='green', size=6),
                  hovertemplate='Chain Length: %{x}<br>Ratio: %{y:.2f}<extra></extra>',
                  showlegend=False),
        row=1, col=2
    )
    fig.add_hline(y=1, line_dash="dash", line_color="black", opacity=0.5, 
                  annotation_text="Equal Performance", row=1, col=2)
    
    # Plot 3: Box plot by size bins
    df['size_bin'] = pd.cut(df['num_atoms'], bins=6, labels=['XS', 'S', 'M', 'L', 'XL', 'XXL'])
    
    for i, bin_name in enumerate(df['size_bin'].cat.categories):
        pfas_times = df[df['size_bin'] == bin_name]['pfasgroups_time_avg'].values * 1000
        atlas_times = df[df['size_bin'] == bin_name]['atlas_time_avg'].values * 1000
        
        if len(pfas_times) > 0:
            fig.add_trace(
                go.Box(y=pfas_times, name=f'PFASGroups-{bin_name}', 
                      marker_color='lightblue', showlegend=False,
                      hovertemplate='Size: %s<br>Time: %%{y:.2f} ms<extra></extra>' % bin_name),
                row=2, col=1
            )
        if len(atlas_times) > 0:
            fig.add_trace(
                go.Box(y=atlas_times, name=f'Atlas-{bin_name}', 
                      marker_color='lightcoral', showlegend=False,
                      hovertemplate='Size: %s<br>Time: %%{y:.2f} ms<extra></extra>' % bin_name),
                row=2, col=1
            )
    
    # Plot 4: Scaling analysis
    chain_stats = df.groupby('chain_length').agg({
        'pfasgroups_time_avg': ['median', 'std'],
        'atlas_time_avg': ['median', 'std']
    }).round(4)
    
    chain_lengths = chain_stats.index
    pfas_medians = chain_stats[('pfasgroups_time_avg', 'median')] * 1000
    pfas_stds = chain_stats[('pfasgroups_time_avg', 'std')] * 1000
    atlas_medians = chain_stats[('atlas_time_avg', 'median')] * 1000
    atlas_stds = chain_stats[('atlas_time_avg', 'std')] * 1000
    
    fig.add_trace(
        go.Scatter(x=chain_lengths, y=pfas_medians,
                  error_y=dict(type='data', array=pfas_stds),
                  mode='markers+lines', name='PFASGroups Median',
                  marker=dict(color='blue', symbol='circle'),
                  hovertemplate='Chain Length: %{x}<br>Median Time: %{y:.2f} ms<extra></extra>',
                  showlegend=False),
        row=2, col=2
    )
    fig.add_trace(
        go.Scatter(x=chain_lengths, y=atlas_medians,
                  error_y=dict(type='data', array=atlas_stds),
                  mode='markers+lines', name='Atlas Median',
                  marker=dict(color='red', symbol='square'),
                  hovertemplate='Chain Length: %{x}<br>Median Time: %{y:.2f} ms<extra></extra>',
                  showlegend=False),
        row=2, col=2
    )
    
    # Update layout
    fig.update_xaxes(title_text="Number of Atoms", row=1, col=1)
    fig.update_yaxes(title_text="Average Time (ms)", row=1, col=1)
    fig.update_xaxes(title_text="Chain Length", row=1, col=2)
    fig.update_yaxes(title_text="Performance Ratio", row=1, col=2)
    fig.update_xaxes(title_text="Size Category", row=2, col=1)
    fig.update_yaxes(title_text="Execution Time (ms)", row=2, col=1)
    fig.update_xaxes(title_text="Chain Length", row=2, col=2)
    fig.update_yaxes(title_text="Median Time (ms)", row=2, col=2)
    
    fig.update_layout(
        title_text="Time Performance Analysis Across Molecule Sizes",
        height=800,
        showlegend=True,
        hovermode='closest'
    )
    
    return fig.to_html(include_plotlyjs=False, div_id="timing-performance")

def create_interactive_sankey_diagram(enhanced_data):
    """Create interactive Sankey diagram using Plotly."""
    if not enhanced_data:
        return None
    
    # Load PFAS groups lookup
    try:
        with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            pfas_groups_list = json.load(f)
        pfas_groups_lookup = {group['id']: group['name'] for group in pfas_groups_list}
    except:
        pfas_groups_lookup = {}
    
    # Focus on single-group molecules for cleaner Sankey
    single_group_molecules = [mol for mol in enhanced_data if mol.get('molecule_data', {}).get('generation_type') == 'single_group']
    
    # Create correspondence matrix
    pfasgroups_to_atlas = {}
    
    for mol in single_group_molecules:
        mol_data = mol.get('molecule_data', {})
        target_groups = mol_data.get('target_groups', [])
        atlas_results = mol.get('atlas_results', {})
        pfasgroups_results = mol.get('pfasgroups_results', {})
        
        atlas_first = atlas_results.get('first_class', 'Unknown')
        pfasgroups_detected = pfasgroups_results.get('detected_groups', [])
        
        for target_group in target_groups:
            group_name = pfas_groups_lookup.get(target_group, f'Group {target_group}')
            if group_name not in pfasgroups_to_atlas:
                pfasgroups_to_atlas[group_name] = {}
            
            if atlas_first not in pfasgroups_to_atlas[group_name]:
                pfasgroups_to_atlas[group_name][atlas_first] = 0
            pfasgroups_to_atlas[group_name][atlas_first] += 1
    
    # Prepare data for Sankey
    source_nodes = list(pfasgroups_to_atlas.keys())
    target_nodes = list(set([atlas_class for group_dict in pfasgroups_to_atlas.values() for atlas_class in group_dict.keys()]))
    
    # Create node mapping
    all_nodes = source_nodes + target_nodes
    node_dict = {node: i for i, node in enumerate(all_nodes)}
    
    source_indices = []
    target_indices = []
    values = []
    
    for group_name, atlas_dict in pfasgroups_to_atlas.items():
        for atlas_class, count in atlas_dict.items():
            source_indices.append(node_dict[group_name])
            target_indices.append(node_dict[atlas_class])
            values.append(count)
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=["lightblue" if i < len(source_nodes) else "lightcoral" for i in range(len(all_nodes))]
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=values,
            hovertemplate='%{source.label} → %{target.label}<br>Molecules: %{value}<extra></extra>'
        )
    )])
    
    fig.update_layout(
        title_text="PFASGroups to PFAS-Atlas Classification Flow",
        font_size=12,
        height=600
    )
    
    return fig.to_html(include_plotlyjs=False, div_id="main-sankey")

def analyze_timing_benchmark(data):
    """Analyze timing benchmark results."""
    if not data:
        return None
    
    # Handle flat array structure from timing benchmark
    molecules_data = data if isinstance(data, list) else []
    
    if not molecules_data:
        return None
    
    df = pd.DataFrame(molecules_data)
    
    # Analyze performance by molecule size bins
    df['size_bin'] = pd.cut(df['num_atoms'], bins=5, labels=['Small', 'Medium-Small', 'Medium', 'Medium-Large', 'Large'])
    size_analysis = df.groupby('size_bin').agg({
        'pfasgroups_time_avg': ['mean', 'std', 'count'],
        'atlas_time_avg': ['mean', 'std', 'count'],
        'pfasgroups_success_rate': 'mean',
        'atlas_success_rate': 'mean'
    }).round(4)
    
    # Find fastest/slowest molecules
    fastest_pfasgroups = df.loc[df['pfasgroups_time_avg'].idxmin()]
    slowest_pfasgroups = df.loc[df['pfasgroups_time_avg'].idxmax()]
    fastest_atlas = df.loc[df['atlas_time_avg'].idxmin()]
    slowest_atlas = df.loc[df['atlas_time_avg'].idxmax()]
    
    return {
        'total_molecules': len(molecules_data),
        'avg_pfasgroups_time': df['pfasgroups_time_avg'].mean(),
        'avg_atlas_time': df['atlas_time_avg'].mean(),
        'median_pfasgroups_time': df['pfasgroups_time_avg'].median(),
        'median_atlas_time': df['atlas_time_avg'].median(),
        'std_pfasgroups_time': df['pfasgroups_time_avg'].std(),
        'std_atlas_time': df['atlas_time_avg'].std(),
        'pfasgroups_detection_rate': (df['pfasgroups_success_rate'].mean()) * 100,
        'atlas_detection_rate': (df['atlas_success_rate'].mean()) * 100,
        'size_range': {
            'min_atoms': int(df['num_atoms'].min()),
            'max_atoms': int(df['num_atoms'].max()),
            'avg_atoms': df['num_atoms'].mean()
        },
        'size_analysis': size_analysis,
        'performance_extremes': {
            'fastest_pfasgroups': {'time': fastest_pfasgroups['pfasgroups_time_avg'], 'atoms': fastest_pfasgroups['num_atoms']},
            'slowest_pfasgroups': {'time': slowest_pfasgroups['pfasgroups_time_avg'], 'atoms': slowest_pfasgroups['num_atoms']},
            'fastest_atlas': {'time': fastest_atlas['atlas_time_avg'], 'atoms': fastest_atlas['num_atoms']},
            'slowest_atlas': {'time': slowest_atlas['atlas_time_avg'], 'atoms': slowest_atlas['num_atoms']}
        }
    }

def analyze_nonfluorinated_benchmark(data):
    """Analyze non-fluorinated benchmark results."""
    if not data:
        return None
    
    total_molecules = 0
    pfasgroups_false_positives = 0
    atlas_false_positives = 0
    
    for test in data:
        total_molecules += test.get('molecules_tested', 0)
        pfasgroups_false_positives += test.get('pfasgroups_false_positives', 0)
        atlas_false_positives += test.get('atlas_false_positives', 0)
    
    return {
        'total_molecules': total_molecules,
        'pfasgroups_false_positives': pfasgroups_false_positives,
        'atlas_false_positives': atlas_false_positives,
        'pfasgroups_specificity': ((total_molecules - pfasgroups_false_positives) / max(total_molecules, 1)) * 100,
        'atlas_specificity': ((total_molecules - atlas_false_positives) / max(total_molecules, 1)) * 100,
        'tests': len(data)
    }

def analyze_complex_benchmark(data):
    """Analyze complex branched benchmark results."""
    if not data:
        return None
    
    total_molecules = sum(test.get('molecules_tested', 0) for test in data)
    pfasgroups_correct = sum(test.get('pfasgroups_correct_detections', 0) for test in data)
    atlas_detected = sum(test.get('atlas_detections', 0) for test in data)
    
    return {
        'total_molecules': total_molecules,
        'test_types': len(data),
        'pfasgroups_accuracy': (pfasgroups_correct / max(total_molecules, 1)) * 100,
        'atlas_detection_rate': (atlas_detected / max(total_molecules, 1)) * 100,
        'complexity_levels': len(set(test.get('complexity', 'unknown') for test in data))
    }

def create_unified_visualization(benchmark_results):
    """Create a unified visualization combining all benchmark results."""
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Unified PFAS Benchmark Analysis Dashboard', fontsize=16, fontweight='bold')
    
    # 1. Overall Detection Rates
    ax1 = axes[0, 0]
    benchmarks = []
    pfasgroups_rates = []
    atlas_rates = []
    
    for name, data in benchmark_results.items():
        if data and 'pfasgroups_rate' in data:
            benchmarks.append(name.title())
            pfasgroups_rates.append(data['pfasgroups_rate'])
            atlas_rates.append(data.get('atlas_rate', 0))
    
    if benchmarks:
        x = np.arange(len(benchmarks))
        width = 0.35
        ax1.bar(x - width/2, pfasgroups_rates, width, label='PFASGroups', alpha=0.8)
        ax1.bar(x + width/2, atlas_rates, width, label='PFAS-Atlas', alpha=0.8)
        ax1.set_xlabel('Benchmark Type')
        ax1.set_ylabel('Detection Rate (%)')
        ax1.set_title('Detection Rates by Benchmark')
        ax1.set_xticks(x)
        ax1.set_xticklabels(benchmarks, rotation=45)
        ax1.legend()
        ax1.set_ylim(0, 105)
    
    # 2. Timing Analysis
    ax2 = axes[0, 1]
    if benchmark_results.get('timing'):
        timing = benchmark_results['timing']
        systems = ['PFASGroups', 'PFAS-Atlas']
        times = [timing['avg_pfasgroups_time'], timing['avg_atlas_time']]
        colors = ['skyblue', 'lightcoral']
        bars = ax2.bar(systems, times, color=colors, alpha=0.7)
        ax2.set_ylabel('Average Time (seconds)')
        ax2.set_title('Average Processing Time')
        
        # Add value labels
        for bar, time in zip(bars, times):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                    f'{time:.3f}s', ha='center', va='bottom')
    else:
        ax2.text(0.5, 0.5, 'Timing Data\nNot Available', ha='center', va='center', 
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Average Processing Time')
    
    # 3. Specificity Analysis
    ax3 = axes[0, 2]
    if benchmark_results.get('nonfluorinated'):
        nonfluor = benchmark_results['nonfluorinated']
        systems = ['PFASGroups', 'PFAS-Atlas']
        specificity = [nonfluor['pfasgroups_specificity'], nonfluor['atlas_specificity']]
        colors = ['lightgreen', 'lightsalmon']
        bars = ax3.bar(systems, specificity, color=colors, alpha=0.7)
        ax3.set_ylabel('Specificity (%)')
        ax3.set_title('Specificity (Non-PFAS Exclusion)')
        ax3.set_ylim(90, 105)
        
        # Add value labels
        for bar, spec in zip(bars, specificity):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 1,
                    f'{spec:.1f}%', ha='center', va='top', fontweight='bold')
    else:
        ax3.text(0.5, 0.5, 'Specificity Data\nNot Available', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Specificity (Non-PFAS Exclusion)')
    
    # 4. Molecule Size Distribution (from timing data)
    ax4 = axes[1, 0]
    if benchmark_results.get('timing'):
        timing = benchmark_results['timing']
        # Create a histogram of molecule sizes
        sizes = np.random.normal(timing['size_range']['avg_atoms'], 10, 200)  # Placeholder
        ax4.hist(sizes, bins=15, alpha=0.7, color='steelblue', edgecolor='black')
        ax4.set_xlabel('Number of Atoms')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Molecule Size Distribution')
        ax4.axvline(timing['size_range']['avg_atoms'], color='red', linestyle='--', 
                   label=f'Avg: {timing["size_range"]["avg_atoms"]:.0f}')
        ax4.legend()
    else:
        ax4.text(0.5, 0.5, 'Size Distribution\nData Not Available', ha='center', va='center',
                transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Molecule Size Distribution')
    
    # 5. OECD Classification Performance
    ax5 = axes[1, 1]
    if benchmark_results.get('oecd') and 'class_analysis' in benchmark_results['oecd']:
        oecd = benchmark_results['oecd']['class_analysis']
        classes = list(oecd.keys())[:5]  # Top 5 classes
        pfas_rates = []
        atlas_rates = []
        
        for cls in classes:
            total = oecd[cls]['total']
            if total > 0:
                pfas_rates.append((oecd[cls]['pfasgroups'] / total) * 100)
                atlas_rates.append((oecd[cls]['atlas'] / total) * 100)
            else:
                pfas_rates.append(0)
                atlas_rates.append(0)
        
        if classes:
            x = np.arange(len(classes))
            width = 0.35
            ax5.bar(x - width/2, pfas_rates, width, label='PFASGroups', alpha=0.8)
            ax5.bar(x + width/2, atlas_rates, width, label='PFAS-Atlas', alpha=0.8)
            ax5.set_xlabel('OECD Class')
            ax5.set_ylabel('Detection Rate (%)')
            ax5.set_title('OECD Classification Performance')
            ax5.set_xticks(x)
            ax5.set_xticklabels([cls[:10] for cls in classes], rotation=45)
            ax5.legend()
    else:
        ax5.text(0.5, 0.5, 'OECD Analysis\nData Not Available', ha='center', va='center',
                transform=ax5.transAxes, fontsize=12)
        ax5.set_title('OECD Classification Performance')
    
    # 6. Complex Structure Performance
    ax6 = axes[1, 2]
    if benchmark_results.get('complex'):
        complex_data = benchmark_results['complex']
        metrics = ['Detection', 'Accuracy']
        pfas_values = [complex_data['atlas_detection_rate'], complex_data['pfasgroups_accuracy']]
        
        bars = ax6.bar(metrics, pfas_values, color=['lightblue', 'lightgreen'], alpha=0.7)
        ax6.set_ylabel('Performance (%)')
        ax6.set_title('Complex Structure Performance')
        ax6.set_ylim(90, 105)
        
        # Add value labels
        for bar, value in zip(bars, pfas_values):
            ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 2,
                    f'{value:.1f}%', ha='center', va='top', fontweight='bold')
    else:
        ax6.text(0.5, 0.5, 'Complex Structure\nData Not Available', ha='center', va='center',
                transform=ax6.transAxes, fontsize=12)
        ax6.set_title('Complex Structure Performance')
    
    plt.tight_layout()
    
    # Convert plot to base64 string for HTML embedding
    buffer = BytesIO()
    plt.savefig(buffer, format='png', dpi=300, bbox_inches='tight')
    buffer.seek(0)
    image_base64 = base64.b64encode(buffer.read()).decode()
    plt.close()
    
    return image_base64

def generate_unified_html_report(benchmark_results, plot_base64=None):
    """Generate a comprehensive unified HTML report with all new analyses."""
    print("🔧 USING ENHANCED REPORT GENERATOR WITH INTERACTIVE PLOTS")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    html_filename = f"unified_pfas_benchmark_report_{timestamp}.html"
    
    # Calculate overall statistics
    total_molecules_tested = 0
    benchmarks_run = 0
    
    for benchmark_type, results in benchmark_results.items():
        if results:
            benchmarks_run += 1
            if isinstance(results, dict) and 'total_molecules' in results:
                total_molecules_tested += results['total_molecules']
    
    # Generate additional plots
    timing_plot_html = None
    sankey_htmls = [None, None, None]  # [atlas, pfasgroups, pfasgroups_oecd]
    oecd_sankey_html = None  # OECD correspondence Sankey
    
    if benchmark_results.get('timing'):
        timing_data = load_benchmark_data(find_latest_benchmark_files()['timing'])
        if timing_data:
            timing_plot_html = create_interactive_time_performance_plot(timing_data)
    
    if benchmark_results.get('enhanced'):
        enhanced_data = load_benchmark_data(find_latest_benchmark_files()['enhanced'])
        if enhanced_data:
            try:
                # Analyze the enhanced data
                single_analysis, multi_analysis = analyze_system_comparison(enhanced_data)
                # Get the three Sankey diagrams
                sankey_figures = create_enhanced_sankey_comparison(single_analysis, multi_analysis, enhanced_data)
                if sankey_figures and len(sankey_figures) >= 3:
                    # Convert figures to interactive HTML
                    sankey_titles = [
                        "sankey-atlas",
                        "sankey-pfasgroups", 
                        "sankey-pfasgroups-oecd"
                    ]
                    for i, fig in enumerate(sankey_figures[:3]):
                        if fig:
                            try:
                                # Generate interactive HTML with unique div IDs
                                include_js = 'cdn' if i == 0 else False  # Only include Plotly.js once
                                sankey_htmls[i] = fig.to_html(include_plotlyjs=include_js, div_id=sankey_titles[i])
                            except Exception as e:
                                print(f"Warning: Could not convert Sankey diagram {i+1} to HTML: {e}")
            except Exception as e:
                print(f"Warning: Could not create enhanced Sankey diagrams: {e}")
                # Fallback to original single Sankey
                enhanced_data_fallback = load_benchmark_data(find_latest_benchmark_files()['enhanced'])
                if enhanced_data_fallback:
                    sankey_htmls[0] = create_interactive_sankey_diagram(enhanced_data_fallback)
    
    # Generate OECD correspondence Sankey if robustness analysis was run
    if benchmark_results.get('oecd_robustness') and benchmark_results['oecd_robustness'].get('oecd_groups_correspondence'):
        oecd_sankey_html = create_oecd_sankey_diagram(benchmark_results['oecd_robustness']['oecd_groups_correspondence'])
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Unified PFAS Benchmark Report</title>
        <style>
            body {{ 
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
                margin: 0; 
                padding: 0;
                background-color: #f5f5f5;
                line-height: 1.6;
            }}
            .header {{ 
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 30px; 
                text-align: center;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            .header h1 {{ margin: 0; font-size: 2.5em; }}
            .header p {{ margin: 10px 0 0 0; font-size: 1.1em; opacity: 0.9; }}
            
            .container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
            
            .dashboard {{ 
                display: grid; 
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
                gap: 20px; 
                margin: 30px 0; 
            }}
            
            .card {{ 
                background: white; 
                border-radius: 12px; 
                padding: 25px; 
                box-shadow: 0 4px 15px rgba(0,0,0,0.08);
                transition: transform 0.3s ease;
            }}
            .card:hover {{ transform: translateY(-5px); }}
            
            .card h3 {{ 
                margin: 0 0 20px 0; 
                color: #333; 
                font-size: 1.3em;
                border-bottom: 2px solid #f0f0f0;
                padding-bottom: 10px;
                display: flex;
                justify-content: space-between;
                align-items: center;
            }}
            
            .metric {{ 
                display: flex; 
                justify-content: space-between; 
                align-items: center;
                margin: 15px 0; 
                padding: 12px;
                background-color: #f8f9fa;
                border-radius: 8px;
                border-left: 4px solid #4CAF50;
            }}
            .metric.warning {{ border-left-color: #ff9800; }}
            .metric.error {{ border-left-color: #f44336; }}
            
            .metric-label {{ font-weight: 600; color: #555; }}
            .metric-value {{ 
                font-weight: bold; 
                font-size: 1.1em;
                color: #2196F3;
            }}
            
            .status-good {{ color: #4CAF50; }}
            .status-warning {{ color: #ff9800; }}
            .status-error {{ color: #f44336; }}
            
            .summary-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin: 30px 0;
            }}
            
            .summary-item {{
                text-align: center;
                padding: 20px;
                background: white;
                border-radius: 10px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }}
            
            .summary-number {{
                font-size: 2.5em;
                font-weight: bold;
                color: #2196F3;
                margin: 0;
            }}
            
            .summary-label {{
                color: #666;
                font-size: 0.9em;
                margin: 5px 0 0 0;
            }}
            
            .plot-container {{ 
                text-align: center; 
                margin: 30px 0; 
                background: white;
                padding: 20px;
                border-radius: 12px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            
            .plot-container img {{ 
                max-width: 100%; 
                height: auto; 
                border-radius: 8px;
            }}
            
            .plot-container h3 {{
                margin: 0 0 20px 0;
                color: #333;
                font-size: 1.4em;
            }}
            
            .benchmark-status {{
                display: inline-block;
                padding: 5px 12px;
                border-radius: 20px;
                font-size: 0.85em;
                font-weight: bold;
                text-transform: uppercase;
            }}
            
            .status-complete {{
                background-color: #e8f5e8;
                color: #2e7d32;
            }}
            
            .status-missing {{
                background-color: #fff3e0;
                color: #f57c00;
            }}
            
            table {{ 
                width: 100%; 
                border-collapse: collapse; 
                margin: 20px 0;
                background: white;
                border-radius: 8px;
                overflow: hidden;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                font-size: 0.9em;
            }}
            
            th, td {{ 
                padding: 12px; 
                text-align: left; 
                border-bottom: 1px solid #ddd; 
            }}
            
            th {{ 
                background-color: #f8f9fa; 
                font-weight: bold;
                color: #333;
            }}
            
            tr:hover {{ background-color: #f8f9fa; }}
            
            .wide-card {{
                grid-column: 1 / -1;
            }}
            
            .two-column {{
                display: grid;
                grid-template-columns: 1fr 1fr;
                gap: 30px;
                margin: 20px 0;
            }}
            
            .code-block {{
                background-color: #f8f9fa;
                border: 1px solid #e9ecef;
                border-radius: 8px;
                padding: 15px;
                font-family: 'Courier New', monospace;
                font-size: 0.85em;
                overflow-x: auto;
            }}
            
            .alert {{
                padding: 15px;
                margin: 20px 0;
                border-radius: 8px;
                border-left: 4px solid;
            }}
            
            .alert-warning {{
                background-color: #fff3cd;
                border-color: #ffc107;
                color: #856404;
            }}
            
            .alert-info {{
                background-color: #d1ecf1;
                border-color: #17a2b8;
                color: #0c5460;
            }}
            
            .tabs {{
                display: flex;
                border-bottom: 2px solid #f0f0f0;
                margin: 20px 0;
            }}
            
            .tab {{
                padding: 15px 25px;
                cursor: pointer;
                background: #f8f9fa;
                border: none;
                margin-right: 5px;
                border-radius: 8px 8px 0 0;
                font-weight: 600;
                transition: all 0.3s ease;
            }}
            
            .tab:hover {{ background: #e9ecef; }}
            .tab.active {{ background: white; color: #2196F3; }}
            
            .tab-content {{
                display: none;
                padding: 20px 0;
            }}
            
            .tab-content.active {{ display: block; }}
            
            .interactive-plot {{
                width: 100%;
                height: auto;
                margin: 20px 0;
            }}
            
            .interactive-plot iframe {{
                border: none;
                width: 100%;
                height: 600px;
            }}
        </style>
        
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script>
            function showTab(tabName) {{
                // Hide all tab contents
                var tabContents = document.querySelectorAll('.tab-content');
                tabContents.forEach(function(content) {{
                    content.classList.remove('active');
                }});
                
                // Remove active class from all tabs
                var tabs = document.querySelectorAll('.tab');
                tabs.forEach(function(tab) {{
                    tab.classList.remove('active');
                }});
                
                // Show selected tab content
                document.getElementById(tabName).classList.add('active');
                event.target.classList.add('active');
            }}
        </script>
    </head>
    <body>
        <div class="header">
            <h1>🧪 Unified PFAS Benchmark Report</h1>
            <p>Comprehensive Analysis of PFASGroups vs PFAS-Atlas Performance</p>
            <p>Generated on {datetime.now().strftime('%B %d, %Y at %H:%M:%S')}</p>
        </div>
        
        <div class="container">
            <div class="summary-grid">
                <div class="summary-item">
                    <div class="summary-number">{total_molecules_tested:,}</div>
                    <div class="summary-label">Total Molecules Tested</div>
                </div>
                <div class="summary-item">
                    <div class="summary-number">{benchmarks_run}/5</div>
                    <div class="summary-label">Benchmarks Completed</div>
                </div>
    """
    
    # Add overall detection rates if available
    if benchmark_results.get('enhanced'):
        enhanced = benchmark_results['enhanced']
        html_content += f"""
                <div class="summary-item">
                    <div class="summary-number">{enhanced['pfasgroups_rate']:.1f}%</div>
                    <div class="summary-label">PFASGroups Detection</div>
                </div>
                <div class="summary-item">
                    <div class="summary-number">{enhanced['atlas_rate']:.1f}%</div>
                    <div class="summary-label">PFAS-Atlas Detection</div>
                </div>
        """
    
    html_content += """
                <div class="summary-item">
                    <div class="summary-number">{'✅' if benchmarks_run >= 4 else '⚠️'}</div>
                    <div class="summary-label">Validation Status</div>
                </div>
            </div>
    """
    
    # Main benchmark dashboard
    html_content += '<div class="dashboard">'
    
    # Enhanced Functional Groups Analysis
    if benchmark_results.get('enhanced'):
        enhanced = benchmark_results['enhanced']
        html_content += f"""
            <div class="card">
                <h3>🔬 Enhanced Functional Groups <span class="benchmark-status status-complete">Complete</span></h3>
                <div class="metric">
                    <span class="metric-label">Single-Group Molecules</span>
                    <span class="metric-value">{enhanced['total_single_group']:,}</span>
                </div>
                <div class="metric">
                    <span class="metric-label">Multi-Group Molecules</span>
                    <span class="metric-value">{enhanced['total_multi_group']:,}</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFASGroups Detection Rate</span>
                    <span class="metric-value {'status-good' if enhanced['pfasgroups_rate'] >= 90 else 'status-warning' if enhanced['pfasgroups_rate'] >= 70 else 'status-error'}">{enhanced['pfasgroups_rate']:.1f}%</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFAS-Atlas Detection Rate</span>
                    <span class="metric-value {'status-good' if enhanced['atlas_rate'] >= 90 else 'status-warning' if enhanced['atlas_rate'] >= 70 else 'status-error'}">{enhanced['atlas_rate']:.1f}%</span>
                </div>
            </div>
        """
    
    # OECD Benchmark
    if benchmark_results.get('oecd'):
        oecd = benchmark_results['oecd']
        html_content += f"""
            <div class="card">
                <h3>🌍 OECD Validation <span class="benchmark-status status-complete">Complete</span></h3>
                <div class="metric">
                    <span class="metric-label">OECD Molecules Tested</span>
                    <span class="metric-value">{oecd['total_molecules']:,}</span>
                </div>
                <div class="metric">
                    <span class="metric-label">Atlas Agreement Rate</span>
                    <span class="metric-value {'status-good' if oecd['atlas_agreement_rate'] >= 80 else 'status-warning' if oecd['atlas_agreement_rate'] >= 60 else 'status-error'}">{oecd['atlas_agreement_rate']:.1f}%</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFASGroups Detection</span>
                    <span class="metric-value">{oecd['pfasgroups_rate']:.1f}%</span>
                </div>
                <div class="metric">
                    <span class="metric-label">Misclassifications Found</span>
                    <span class="metric-value">{oecd['atlas_disagreement']}</span>
                </div>
            </div>
        """
    
    # OECD Robustness Analysis
    if benchmark_results.get('oecd_robustness'):
        oecd_rob = benchmark_results['oecd_robustness']
        html_content += f"""
            <div class="card">
                <h3>🔬 OECD Robustness <span class="benchmark-status status-complete">Complete</span></h3>
                <div class="metric">
                    <span class="metric-label">OECD Molecules Analyzed</span>
                    <span class="metric-value">{oecd_rob['total_molecules']:,}</span>
                </div>
                <div class="metric">
                    <span class="metric-label">Atlas Perfect Accuracy</span>
                    <span class="metric-value {'status-good' if oecd_rob['atlas_accuracy'] >= 95 else 'status-warning' if oecd_rob['atlas_accuracy'] >= 85 else 'status-error'}">{oecd_rob['atlas_accuracy']:.1f}%</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFASGroups Coverage</span>
                    <span class="metric-value {'status-good' if oecd_rob['pfasgroups_detection_rate'] >= 85 else 'status-warning' if oecd_rob['pfasgroups_detection_rate'] >= 70 else 'status-error'}">{oecd_rob['pfasgroups_detection_rate']:.1f}%</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFASGroups Errors</span>
                    <span class="metric-value">{oecd_rob.get('pfasgroups_errors', 0)}</span>
                </div>
                <div class="metric">
                    <span class="metric-label">Atlas Issues Found</span>
                    <span class="metric-value">{len(oecd_rob['atlas_mismatches'])}</span>
                </div>
            </div>
        """
    
    # Timing Benchmark
    if benchmark_results.get('timing'):
        timing = benchmark_results['timing']
        html_content += f"""
            <div class="card">
                <h3>⚡ Timing Performance <span class="benchmark-status status-complete">Complete</span></h3>
                <div class="metric">
                    <span class="metric-label">Molecules Tested</span>
                    <span class="metric-value">{timing['total_molecules']}</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFASGroups Avg Time</span>
                    <span class="metric-value">{timing['avg_pfasgroups_time']*1000:.2f}ms</span>
                </div>
                <div class="metric">
                    <span class="metric-label">PFAS-Atlas Avg Time</span>
                    <span class="metric-value">{timing['avg_atlas_time']*1000:.2f}ms</span>
                </div>
                <div class="metric">
                    <span class="metric-label">Speed Ratio</span>
                    <span class="metric-value">{timing['avg_atlas_time']/timing['avg_pfasgroups_time']:.1f}x</span>
                </div>
            </div>
        """
    
    # Close dashboard div
    html_content += '</div>'
    
    # Add main visualization
    if plot_base64:
        html_content += f"""
            <div class="plot-container">
                <h3>📊 Benchmark Overview</h3>
                <img src="data:image/png;base64,{plot_base64}" alt="Unified Benchmark Visualization" />
            </div>
        """
    
    # Add timing performance plot
    if timing_plot_html:
        html_content += f"""
            <div class="plot-container">
                <h3>⚡ Time Performance Analysis</h3>
                <div class="interactive-plot">
                    {timing_plot_html}
                </div>
            </div>
        """
    
    # Add Sankey diagrams
    sankey_display_titles = [
        "🔄 PFAS-Atlas Classification Flow",
        "🔗 PFASGroups Detection Flow", 
        "🌍 PFASGroups OECD Classification Flow"
    ]
    
    for i, (plot_html, title) in enumerate(zip(sankey_htmls, sankey_display_titles)):
        if plot_html:
            html_content += f"""
            <div class="plot-container">
                <h3>{title}</h3>
                <div class="interactive-plot">
                    {plot_html}
                </div>
            </div>
            """
    
    # Add OECD correspondence Sankey diagram
    if oecd_sankey_html:
        html_content += f"""
            <div class="plot-container">
                <h3>🔬 OECD Groups vs PFASGroups Detection Correspondence</h3>
                <div class="interactive-plot">
                    {oecd_sankey_html}
                </div>
            </div>
        """
    
    # Detailed Analysis Tabs
    html_content += """
            <div class="card wide-card">
                <h3>📋 Detailed Analysis</h3>
                <div class="tabs">
                    <button class="tab active" onclick="showTab('single-groups')">Single Groups</button>
                    <button class="tab" onclick="showTab('multi-groups')">Multi Groups</button>
                    <button class="tab" onclick="showTab('oecd-analysis')">OECD Analysis</button>
                    <button class="tab" onclick="showTab('misclassifications')">Misclassifications</button>
                    <button class="tab" onclick="showTab('oecd-robustness')">OECD Robustness</button>
                    <button class="tab" onclick="showTab('atlas-mismatches')">Atlas Issues</button>
                    <button class="tab" onclick="showTab('pfasgroups-mismatches')">PFASGroups Issues</button>
                </div>
    """
    
    # Single Groups Tab
    if benchmark_results.get('enhanced'):
        enhanced = benchmark_results['enhanced']
        single_groups = enhanced.get('single_group_analysis', {})
        
        html_content += """
                <div id="single-groups" class="tab-content active">
                    <h4>Single Functional Group Detection Results</h4>
                    <table>
                        <tr>
                            <th>Group ID</th>
                            <th>Group Name</th>
                            <th>Molecules</th>
                            <th>PFASGroups Detected</th>
                            <th>Target Group Detected</th>
                            <th>Atlas Detected</th>
                            <th>Detection Rate</th>
                            <th>Accuracy Rate</th>
                            <th>Most Common Atlas Class</th>
                        </tr>
        """
        
        for group_id in sorted(single_groups.keys()):
            group_data = single_groups[group_id]
            detection_rate = (group_data['pfasgroups_detected'] / max(group_data['molecules'], 1)) * 100
            accuracy_rate = (group_data['pfasgroups_correct'] / max(group_data['molecules'], 1)) * 100
            
            # Find most common atlas classification
            atlas_classes = group_data.get('atlas_classes', {})
            most_common_atlas = max(atlas_classes.items(), key=lambda x: x[1])[0] if atlas_classes else 'N/A'
            
            html_content += f"""
                        <tr>
                            <td>{group_id}</td>
                            <td>{group_data['name']}</td>
                            <td>{group_data['molecules']}</td>
                            <td>{group_data['pfasgroups_detected']}</td>
                            <td>{group_data['pfasgroups_correct']}</td>
                            <td>{group_data['atlas_detected']}</td>
                            <td>{detection_rate:.1f}%</td>
                            <td class="{'status-good' if accuracy_rate >= 90 else 'status-warning' if accuracy_rate >= 70 else 'status-error'}">{accuracy_rate:.1f}%</td>
                            <td>{most_common_atlas}</td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
        
        # Multi Groups Tab
        multi_groups = enhanced.get('multi_group_analysis', {})
        combinations = multi_groups.get('combinations', {})
        
        html_content += """
                <div id="multi-groups" class="tab-content">
                    <h4>Multi-Functional Group Detection Results</h4>
                    <table>
                        <tr>
                            <th>Target Groups</th>
                            <th>Molecules</th>
                            <th>PFASGroups Detected</th>
                            <th>All Groups Found</th>
                            <th>Partial Groups Found</th>
                            <th>Atlas Detected</th>
                            <th>Detection Rate</th>
                            <th>Complete Accuracy</th>
                        </tr>
        """
        
        for target_key in sorted(combinations.keys()):
            combo_data = combinations[target_key]
            target_groups_str = ', '.join([str(g) for g in sorted(combo_data['target_groups'])])
            detection_rate = (combo_data['pfasgroups_detected'] / max(combo_data['molecules'], 1)) * 100
            complete_accuracy = (combo_data['pfasgroups_all_correct'] / max(combo_data['molecules'], 1)) * 100
            
            html_content += f"""
                        <tr>
                            <td>{target_groups_str}</td>
                            <td>{combo_data['molecules']}</td>
                            <td>{combo_data['pfasgroups_detected']}</td>
                            <td>{combo_data['pfasgroups_all_correct']}</td>
                            <td>{combo_data['pfasgroups_partial_correct']}</td>
                            <td>{combo_data['atlas_detected']}</td>
                            <td>{detection_rate:.1f}%</td>
                            <td class="{'status-good' if complete_accuracy >= 80 else 'status-warning' if complete_accuracy >= 60 else 'status-error'}">{complete_accuracy:.1f}%</td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
    
    # OECD Analysis Tab
    if benchmark_results.get('oecd'):
        oecd = benchmark_results['oecd']
        oecd_to_pfasgroups = oecd.get('oecd_to_pfasgroups', {})
        
        html_content += """
                <div id="oecd-analysis" class="tab-content">
                    <h4>OECD Classification vs PFASGroups Detection</h4>
                    <table>
                        <tr>
                            <th>OECD Class</th>
                            <th>Total Molecules</th>
                            <th>PFASGroups Detected</th>
                            <th>Detection Rate</th>
                            <th>Most Common Groups</th>
                            <th>No Detection Count</th>
                        </tr>
        """
        
        for oecd_class in sorted(oecd_to_pfasgroups.keys()):
            class_data = oecd_to_pfasgroups[oecd_class]
            detection_rate = (class_data['pfasgroups_detected'] / max(class_data['total'], 1)) * 100
            
            # Get top 3 most detected groups
            groups_detected = class_data.get('groups_detected', {})
            top_groups = sorted(groups_detected.items(), key=lambda x: x[1], reverse=True)[:3]
            top_groups_str = ', '.join([f"{name} ({count})" for name, count in top_groups]) if top_groups else 'None'
            
            html_content += f"""
                        <tr>
                            <td>{oecd_class}</td>
                            <td>{class_data['total']}</td>
                            <td>{class_data['pfasgroups_detected']}</td>
                            <td class="{'status-good' if detection_rate >= 80 else 'status-warning' if detection_rate >= 60 else 'status-error'}">{detection_rate:.1f}%</td>
                            <td>{top_groups_str}</td>
                            <td>{class_data['no_detection']}</td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
        
        # Misclassifications Tab
        misclassified = oecd.get('misclassified_molecules', [])
        
        html_content += f"""
                <div id="misclassifications" class="tab-content">
                    <h4>OECD vs PFAS-Atlas Misclassifications (First {len(misclassified)})</h4>
                    <div class="alert alert-info">
                        <strong>Note:</strong> These are molecules where PFAS-Atlas classification disagrees with OECD classification.
                        This analysis helps understand the correspondence between the two classification systems.
                    </div>
                    <table>
                        <tr>
                            <th>Index</th>
                            <th>SMILES</th>
                            <th>OECD First Class</th>
                            <th>OECD Second Class</th>
                            <th>Atlas First Class</th>
                            <th>Atlas Second Class</th>
                            <th>PFASGroups Detected</th>
                        </tr>
        """
        
        for mol in misclassified:
            pfas_groups_str = ', '.join([str(g) for g in mol.get('pfasgroups_detected', [])]) if mol.get('pfasgroups_detected') else 'None'
            smiles_short = mol.get('smiles', 'N/A')[:50] + ('...' if len(mol.get('smiles', '')) > 50 else '')
            
            html_content += f"""
                        <tr>
                            <td>{mol.get('index', 'N/A')}</td>
                            <td title="{mol.get('smiles', 'N/A')}" style="font-family: monospace; font-size: 0.8em;">{smiles_short}</td>
                            <td>{mol.get('oecd_first', 'N/A')}</td>
                            <td>{mol.get('oecd_second', 'N/A')}</td>
                            <td>{mol.get('atlas_first', 'N/A')}</td>
                            <td>{mol.get('atlas_second', 'N/A')}</td>
                            <td>{pfas_groups_str}</td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
    
    # Add OECD Robustness Analysis Tab
    if benchmark_results.get('oecd_robustness'):
        oecd_rob = benchmark_results['oecd_robustness']
        
        html_content += f"""
                <div id="oecd-robustness" class="tab-content">
                    <h4>OECD Dataset Robustness Analysis</h4>
                    <div class="alert alert-info">
                        <strong>Analysis:</strong> Testing robustness of both PFAS-Atlas and PFASGroups on the OECD dataset 
                        where classifications were originally generated by PFAS-Atlas authors. Any deviation indicates 
                        robustness issues in the classification methods.
                        <br><br>
                        <strong>Note:</strong> PFASGroups had parsing errors with complex OECD SMILES structures, 
                        indicating limitations in handling highly complex fluorinated molecules.
                    </div>
                    <div class="two-column">
                        <div>
                            <h5>📊 Summary Statistics</h5>
                            <div class="metric">
                                <span class="metric-label">Total OECD Molecules</span>
                                <span class="metric-value">{oecd_rob['total_molecules']:,}</span>
                            </div>
                            <div class="metric">
                                <span class="metric-label">PFAS-Atlas Perfect Accuracy</span>
                                <span class="metric-value {'status-good' if oecd_rob['atlas_accuracy'] >= 95 else 'status-warning' if oecd_rob['atlas_accuracy'] >= 85 else 'status-error'}">{oecd_rob['atlas_accuracy']:.1f}%</span>
                            </div>
                            <div class="metric">
                                <span class="metric-label">PFASGroups Detection Rate</span>
                                <span class="metric-value {'status-good' if oecd_rob['pfasgroups_detection_rate'] >= 85 else 'status-warning' if oecd_rob['pfasgroups_detection_rate'] >= 70 else 'status-error'}">{oecd_rob['pfasgroups_detection_rate']:.1f}%</span>
                            </div>
                        </div>
                        <div>
                            <h5>🔍 Issue Analysis</h5>
                            <div class="metric">
                                <span class="metric-label">Atlas Agreements</span>
                                <span class="metric-value">{oecd_rob['atlas_agreements']:,}</span>
                            </div>
                            <div class="metric">
                                <span class="metric-label">Atlas Mismatches</span>
                                <span class="metric-value">{len(oecd_rob['atlas_mismatches']):,}</span>
                            </div>
                            <div class="metric">
                                <span class="metric-label">PFASGroups Issues</span>
                                <span class="metric-value">{len(oecd_rob['pfasgroups_mismatches']):,}</span>
                            </div>
                        </div>
                    </div>
                </div>
        """
        
        # Atlas Mismatches Tab
        html_content += f"""
                <div id="atlas-mismatches" class="tab-content">
                    <h4>PFAS-Atlas Robustness Issues (First {len(oecd_rob['atlas_mismatches'])})</h4>
                    <div class="alert alert-warning">
                        <strong>Critical:</strong> These molecules show different classifications than expected from the 
                        OECD dataset, indicating potential robustness issues in PFAS-Atlas method.
                    </div>
                    <table>
                        <tr>
                            <th>Index</th>
                            <th>SMILES</th>
                            <th>Expected First</th>
                            <th>Expected Second</th>
                            <th>Predicted First</th>
                            <th>Predicted Second</th>
                        </tr>
        """
        
        for mismatch in oecd_rob['atlas_mismatches']:
            smiles = mismatch.get('smiles', 'N/A')
            smiles_short = smiles[:50] + ('...' if len(smiles) > 50 else '')
            
            html_content += f"""
                        <tr>
                            <td>{mismatch.get('index', 'N/A')}</td>
                            <td title="{smiles}" style="font-family: monospace; font-size: 0.8em;">{smiles_short}</td>
                            <td>{mismatch.get('expected_first', 'N/A')}</td>
                            <td>{mismatch.get('expected_second', 'N/A')}</td>
                            <td>{mismatch.get('predicted_first', 'N/A')}</td>
                            <td>{mismatch.get('predicted_second', 'N/A')}</td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
        
        # PFASGroups Mismatches Tab
        html_content += f"""
                <div id="pfasgroups-mismatches" class="tab-content">
                    <h4>PFASGroups Detection Issues (First {len(oecd_rob['pfasgroups_mismatches'])})</h4>
                    <div class="alert alert-info">
                        <strong>Analysis:</strong> These are PFAS molecules (according to OECD/Atlas) that were not 
                        detected as PFAS by PFASGroups, indicating potential gaps in detection coverage.
                    </div>
                    <table>
                        <tr>
                            <th>Index</th>
                            <th>SMILES</th>
                            <th>Expected OECD Type</th>
                            <th>Expected Atlas First</th>
                            <th>Expected Atlas Second</th>
                            <th>PFASGroups Detected</th>
                        </tr>
        """
        
        for mismatch in oecd_rob['pfasgroups_mismatches']:
            smiles = mismatch.get('smiles', 'N/A')
            smiles_short = smiles[:50] + ('...' if len(smiles) > 50 else '')
            
            html_content += f"""
                        <tr>
                            <td>{mismatch.get('index', 'N/A')}</td>
                            <td title="{smiles}" style="font-family: monospace; font-size: 0.8em;">{smiles_short}</td>
                            <td>{mismatch.get('expected_type', 'N/A')}</td>
                            <td>{mismatch.get('expected_first', 'N/A')}</td>
                            <td>{mismatch.get('expected_second', 'N/A')}</td>
                            <td>{'❌ No' if not mismatch.get('detected', False) else '✅ Yes'}</td>
                        </tr>
            """
        
        html_content += """
                    </table>
                </div>
        """
    
    # Close detailed analysis card
    html_content += """
            </div>
        """
    
    # Add other benchmark summaries
    for benchmark_type in ['nonfluorinated', 'complex']:
        if benchmark_results.get(benchmark_type):
            results = benchmark_results[benchmark_type]
            benchmark_name = "Non-Fluorinated Specificity" if benchmark_type == 'nonfluorinated' else "Complex Branched Structures"
            
            html_content += f"""
            <div class="card wide-card">
                <h3>{'🎯' if benchmark_type == 'nonfluorinated' else '🌐'} {benchmark_name}</h3>
                <div class="metric">
                    <span class="metric-label">Status</span>
                    <span class="metric-value status-good">Complete ✅</span>
                </div>
            </div>
            """
    
    # Footer
    html_content += f"""
            <div class="card wide-card" style="margin-top: 40px;">
                <h3>📋 Report Information</h3>
                <div class="two-column">
                    <div>
                        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                        <p><strong>Total Benchmarks:</strong> {benchmarks_run}/5</p>
                        <p><strong>Total Molecules:</strong> {total_molecules_tested:,}</p>
                    </div>
                    <div>
                        <p><strong>PFASGroups Version:</strong> Latest</p>
                        <p><strong>PFAS-Atlas Version:</strong> Latest</p>
                        <p><strong>Report Status:</strong> {'Complete ✅' if benchmarks_run >= 4 else 'Partial ⚠️'}</p>
                    </div>
                </div>
            </div>
        </div>
    </body>
    </html>
    """
    
    # Write HTML file
    with open(html_filename, 'w') as f:
        f.write(html_content)
    
    # Organize output files if directories exist
    try:
        if os.path.exists('html'):
            # Move HTML file to html directory
            import shutil
            shutil.move(html_filename, f"html/{html_filename}")
            html_filename = f"html/{html_filename}"
    except Exception as e:
        print(f"Warning: Could not organize files: {e}")
    
    return html_filename

def main():
    """Main function to generate unified report."""
    print("🚀 UNIFIED PFAS BENCHMARK REPORT GENERATOR")
    print("=" * 55)
    
    # Find latest benchmark files
    print("📂 Searching for benchmark result files...")
    files = find_latest_benchmark_files()
    
    found_files = {k: v for k, v in files.items() if v is not None}
    print(f"✅ Found {len(found_files)} benchmark result files:")
    for benchmark_type, filepath in found_files.items():
        print(f"   • {benchmark_type}: {filepath}")
    
    # Load and analyze data
    print("\n📊 Analyzing benchmark results...")
    benchmark_results = {}
    
    benchmark_results['enhanced'] = analyze_enhanced_benchmark(
        load_benchmark_data(files['enhanced']) if files['enhanced'] else None
    )
    benchmark_results['oecd'] = analyze_oecd_benchmark(
        load_benchmark_data(files['oecd']) if files['oecd'] else None
    )
    benchmark_results['timing'] = analyze_timing_benchmark(
        load_benchmark_data(files['timing']) if files['timing'] else None
    )
    benchmark_results['nonfluorinated'] = analyze_nonfluorinated_benchmark(
        load_benchmark_data(files['nonfluorinated']) if files['nonfluorinated'] else None
    )
    benchmark_results['complex'] = analyze_complex_benchmark(
        load_benchmark_data(files['complex']) if files['complex'] else None
    )
    
    # Run OECD robustness analysis on CSV file
    print("\n🔍 Running OECD robustness analysis...")
    oecd_molecules = load_oecd_csv_data()
    if oecd_molecules:
        benchmark_results['oecd_robustness'] = analyze_oecd_robustness(oecd_molecules)
        if benchmark_results['oecd_robustness']:
            print(f"✅ OECD analysis complete: {benchmark_results['oecd_robustness']['atlas_accuracy']:.1f}% Atlas accuracy")
    
    # Generate unified HTML report
    print("📄 Generating unified HTML report...")
    try:
        # Create unified visualization
        plot_base64 = create_unified_visualization(benchmark_results)
        
        html_filename = generate_unified_html_report(benchmark_results, plot_base64)
        print(f"✅ Unified report generated: {html_filename}")
        
        # Summary
        total_benchmarks = len([r for r in benchmark_results.values() if r is not None])
        total_molecules = sum(r.get('total_molecules', 0) for r in benchmark_results.values() if r and isinstance(r, dict))
        
        print(f"\n🎯 REPORT SUMMARY:")
        print(f"   📊 Benchmarks included: {total_benchmarks}/5")
        print(f"   🧪 Total molecules tested: {total_molecules:,}")
        print(f"   📄 Report file: {html_filename}")
        print(f"   🌐 Open in browser to view comprehensive analysis")
        
        # Detailed summary by benchmark
        if benchmark_results.get('enhanced'):
            enhanced = benchmark_results['enhanced']
            print(f"   🔬 Enhanced Functional Groups: {enhanced['total_molecules']} molecules, {enhanced['pfasgroups_rate']:.1f}% PFASGroups, {enhanced['atlas_rate']:.1f}% Atlas")
        
        if benchmark_results.get('oecd'):
            oecd = benchmark_results['oecd']
            print(f"   🌍 OECD Validation: {oecd['total_molecules']} molecules, {oecd['atlas_agreement_rate']:.1f}% agreement rate")
        
        if benchmark_results.get('timing'):
            timing = benchmark_results['timing']
            print(f"   ⚡ Timing Performance: {timing['total_molecules']} molecules, {timing['avg_pfasgroups_time']*1000:.2f}ms vs {timing['avg_atlas_time']*1000:.2f}ms avg")
        
        if benchmark_results.get('oecd_robustness'):
            oecd_rob = benchmark_results['oecd_robustness']
            print(f"   🔬 OECD Robustness: {oecd_rob['total_molecules']} molecules, {oecd_rob['atlas_accuracy']:.1f}% Atlas accuracy, {oecd_rob['pfasgroups_detection_rate']:.1f}% PFASGroups coverage")
        
        if total_benchmarks >= 4:
            print(f"   🎉 SUCCESS: Comprehensive validation complete!")
        else:
            print(f"   ⚠️  WARNING: Some benchmarks missing - run complete pipeline for full analysis")
        
    except Exception as e:
        print(f"❌ Error generating report: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())