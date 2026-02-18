#!/usr/bin/env python3
"""
Test fluorotelomer detection on PubChem fluorotelomer dataset.

This script:
1. Reads PubChem_fluorotelomers.sdf (molecules from PubChem search "fluorotelomer")
2. Tests detection with HalogenGroups
3. Calculates true/false positives/negatives
4. Saves results to review-app database
5. Generates summary statistics

Author: HalogenGroups Benchmark Team
Date: 2026-01-29
"""

import sys
from pathlib import Path
import json
import sqlite3
from datetime import datetime

from rdkit import Chem
from HalogenGroups import parse_smiles


def load_telomer_groups():
    """Load telomer group IDs automatically from PFAS_groups_smarts.json.
    
    Returns dictionary mapping group ID to group name for all groups containing
    'telomer' in their name (case-insensitive).
    """
    script_dir = Path(__file__).parent.parent.parent
    json_path = script_dir / 'HalogenGroups' / 'data' / 'PFAS_groups_smarts.json'
    
    telomer_groups = {}
    
    try:
        with open(json_path, 'r') as f:
            groups = json.load(f)
        
        for group in groups:
            group_id = group.get('id')
            group_name = group.get('name', '')
            
            # Check if "telomer" appears in the group name (case-insensitive)
            if 'telomer' in group_name.lower():
                telomer_groups[group_id] = group_name
        
        print(f"✓ Loaded {len(telomer_groups)} telomer groups from PFAS_groups_smarts.json")
        print(f"  Telomer group IDs: {sorted(telomer_groups.keys())}")
        
    except Exception as e:
        print(f"⚠ Warning: Could not load telomer groups from JSON: {e}")
        print(f"  Using empty telomer group list")
    
    return telomer_groups


# Automatically load telomer group IDs from PFAS_groups_smarts.json
TELOMER_GROUP_IDS = load_telomer_groups()


def read_sdf_file(sdf_path):
    """Read molecules from SDF file."""
    print(f"Reading SDF file: {sdf_path}")
    
    supplier = Chem.SDMolSupplier(str(sdf_path))
    molecules = []
    
    for i, mol in enumerate(supplier):
        if mol is None:
            print(f"  Warning: Could not parse molecule {i+1}")
            continue
        
        # Get properties
        props = mol.GetPropsAsDict()
        
        # Get SMILES
        smiles = Chem.MolToSmiles(mol)
        
        molecules.append({
            'index': i,
            'mol': mol,
            'smiles': smiles,
            'properties': props
        })
    
    print(f"✓ Loaded {len(molecules)} molecules from SDF")
    return molecules


def test_telomer_detection(molecules):
    """Test telomer detection on all molecules."""
    print(f"\nTesting telomer detection on {len(molecules)} molecules...")
    
    results = []
    
    for i, mol_data in enumerate(molecules, 1):
        smiles = mol_data['smiles']
        props = mol_data['properties']
        
        print(f"[{i}/{len(molecules)}] Testing molecule...", end='\r')
        
        try:
            # Parse with HalogenGroups
            parse_results = parse_smiles(smiles)
            
            if parse_results and len(parse_results) > 0:
                mol_result = parse_results[0]
                
                # Extract detected groups
                detected_groups = []
                telomer_detected = False
                
                if isinstance(mol_result, dict) and 'matches' in mol_result:
                    for match in mol_result['matches']:
                        if match.get('type') == 'HalogenGroup':
                            group_id = match['id']
                            group_name = match['group_name']
                            
                            detected_groups.append({
                                'id': group_id,
                                'name': group_name,
                                'match_count': match.get('match_count', 0)
                            })
                            
                            # Check if it's a telomer group
                            if group_id in TELOMER_GROUP_IDS:
                                telomer_detected = True
                
                # Determine if this is a true/false positive
                # Since these come from PubChem search "fluorotelomer", we assume they SHOULD be telomers
                # But we need manual inspection to be sure
                expected_telomer = True  # From PubChem search
                
                result = {
                    'index': mol_data['index'],
                    'smiles': smiles,
                    'properties': props,
                    'detected_groups': detected_groups,
                    'telomer_detected': telomer_detected,
                    'expected_telomer': expected_telomer,
                    'true_positive': telomer_detected and expected_telomer,
                    'false_negative': not telomer_detected and expected_telomer,
                    'false_positive': telomer_detected and not expected_telomer,
                    'true_negative': not telomer_detected and not expected_telomer,
                    'error': None
                }
                
            else:
                result = {
                    'index': mol_data['index'],
                    'smiles': smiles,
                    'properties': props,
                    'detected_groups': [],
                    'telomer_detected': False,
                    'expected_telomer': True,
                    'true_positive': False,
                    'false_negative': True,
                    'false_positive': False,
                    'true_negative': False,
                    'error': 'No parse results'
                }
                
        except Exception as e:
            result = {
                'index': mol_data['index'],
                'smiles': smiles,
                'properties': props,
                'detected_groups': [],
                'telomer_detected': False,
                'expected_telomer': True,
                'true_positive': False,
                'false_negative': True,
                'false_positive': False,
                'true_negative': False,
                'error': str(e)
            }
        
        results.append(result)
    
    print(f"\n✓ Completed testing {len(results)} molecules")
    return results


def calculate_statistics(results):
    """Calculate summary statistics."""
    stats = {
        'total': len(results),
        'true_positive': sum(1 for r in results if r['true_positive']),
        'false_negative': sum(1 for r in results if r['false_negative']),
        'false_positive': sum(1 for r in results if r['false_positive']),
        'true_negative': sum(1 for r in results if r['true_negative']),
        'errors': sum(1 for r in results if r['error']),
        'telomer_detected': sum(1 for r in results if r['telomer_detected']),
    }
    
    # Calculate rates
    if stats['total'] > 0:
        stats['detection_rate'] = stats['telomer_detected'] / stats['total'] * 100
        stats['error_rate'] = stats['errors'] / stats['total'] * 100
    
    # Calculate precision, recall, F1
    tp = stats['true_positive']
    fp = stats['false_positive']
    fn = stats['false_negative']
    
    if tp + fp > 0:
        stats['precision'] = tp / (tp + fp)
    else:
        stats['precision'] = 0.0
    
    if tp + fn > 0:
        stats['recall'] = tp / (tp + fn)
    else:
        stats['recall'] = 0.0
    
    if stats['precision'] + stats['recall'] > 0:
        stats['f1_score'] = 2 * (stats['precision'] * stats['recall']) / (stats['precision'] + stats['recall'])
    else:
        stats['f1_score'] = 0.0
    
    # Group detection frequency
    group_counts = {}
    for result in results:
        if result['telomer_detected']:
            for group in result['detected_groups']:
                if group['id'] in TELOMER_GROUP_IDS:
                    group_id = group['id']
                    if group_id not in group_counts:
                        group_counts[group_id] = {
                            'id': group_id,
                            'name': TELOMER_GROUP_IDS[group_id],
                            'count': 0
                        }
                    group_counts[group_id]['count'] += 1
    
    stats['group_counts'] = sorted(group_counts.values(), key=lambda x: x['count'], reverse=True)
    
    return stats


def save_to_review_app(results, stats, db_path):
    """Save results to review-app SQLite database."""
    print(f"\nSaving results to review-app database: {db_path}")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Create table for telomer validation if it doesn't exist
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS telomer_validation (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            smiles TEXT NOT NULL,
            pubchem_cid TEXT,
            pubchem_name TEXT,
            telomer_detected BOOLEAN,
            detected_groups TEXT,
            expected_telomer BOOLEAN,
            true_positive BOOLEAN,
            false_negative BOOLEAN,
            error TEXT,
            test_date TEXT,
            UNIQUE(smiles)
        )
    ''')
    
    # Insert results
    for result in results:
        pubchem_cid = result['properties'].get('PUBCHEM_COMPOUND_CID', '')
        pubchem_name = result['properties'].get('PUBCHEM_IUPAC_NAME', '')
        
        detected_groups_json = json.dumps(result['detected_groups'])
        
        cursor.execute('''
            INSERT OR REPLACE INTO telomer_validation 
            (smiles, pubchem_cid, pubchem_name, telomer_detected, detected_groups, 
             expected_telomer, true_positive, false_negative, error, test_date)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            result['smiles'],
            pubchem_cid,
            pubchem_name,
            result['telomer_detected'],
            detected_groups_json,
            result['expected_telomer'],
            result['true_positive'],
            result['false_negative'],
            result['error'],
            datetime.now().isoformat()
        ))
    
    # Save statistics
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS telomer_validation_stats (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            test_date TEXT,
            total INTEGER,
            true_positive INTEGER,
            false_negative INTEGER,
            false_positive INTEGER,
            true_negative INTEGER,
            errors INTEGER,
            detection_rate REAL,
            precision REAL,
            recall REAL,
            f1_score REAL,
            group_counts TEXT
        )
    ''')
    
    cursor.execute('''
        INSERT INTO telomer_validation_stats 
        (test_date, total, true_positive, false_negative, false_positive, true_negative,
         errors, detection_rate, precision, recall, f1_score, group_counts)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        datetime.now().isoformat(),
        stats['total'],
        stats['true_positive'],
        stats['false_negative'],
        stats['false_positive'],
        stats['true_negative'],
        stats['errors'],
        stats['detection_rate'],
        stats['precision'],
        stats['recall'],
        stats['f1_score'],
        json.dumps(stats['group_counts'])
    ))
    
    conn.commit()
    conn.close()
    
    print(f"✓ Saved {len(results)} results and statistics to database")


def print_summary(stats):
    """Print summary statistics."""
    print("\n" + "="*60)
    print("FLUOROTELOMER DETECTION SUMMARY")
    print("="*60)
    print(f"\nDataset: PubChem fluorotelomers")
    print(f"Total molecules: {stats['total']}")
    print(f"\nConfusion Matrix:")
    print(f"  True Positives:  {stats['true_positive']:4d} (telomers correctly detected)")
    print(f"  False Negatives: {stats['false_negative']:4d} (telomers missed)")
    print(f"  False Positives: {stats['false_positive']:4d} (non-telomers incorrectly detected)")
    print(f"  True Negatives:  {stats['true_negative']:4d} (non-telomers correctly rejected)")
    print(f"  Errors:          {stats['errors']:4d} (parsing/detection errors)")
    
    print(f"\nPerformance Metrics:")
    print(f"  Detection Rate: {stats['detection_rate']:.1f}%")
    print(f"  Precision:      {stats['precision']:.3f}")
    print(f"  Recall:         {stats['recall']:.3f}")
    print(f"  F1 Score:       {stats['f1_score']:.3f}")
    
    print(f"\nMost Frequently Detected Telomer Groups:")
    for i, group in enumerate(stats['group_counts'][:10], 1):
        print(f"  {i:2d}. ID {group['id']:2d}: {group['name']:<45s} ({group['count']:3d} detections)")
    
    print("\n" + "="*60)


def main():
    """Main execution."""
    script_dir = Path(__file__).parent.parent
    
    print(f"Script directory: {script_dir}")
    
    # Input file
    sdf_file = script_dir / 'data' / 'PubChem_fluorotelomers.sdf'
    
    print(f"Looking for SDF file: {sdf_file}")
    print(f"File exists: {sdf_file.exists()}")
    
    if not sdf_file.exists():
        print(f"Error: SDF file not found: {sdf_file}")
        return
    
    # Output paths
    results_file = script_dir / 'data' / 'telomer_validation_results.json'
    db_path = script_dir / 'review-app' / 'pfas_benchmark.db'
    
    # Read molecules
    molecules = read_sdf_file(sdf_file)
    
    # Test detection
    results = test_telomer_detection(molecules)
    
    # Calculate statistics
    stats = calculate_statistics(results)
    
    # Print summary
    print_summary(stats)
    
    # Save results to JSON
    output_data = {
        'test_date': datetime.now().isoformat(),
        'dataset': 'PubChem_fluorotelomers.sdf',
        'statistics': stats,
        'results': results
    }
    
    with open(results_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n✓ Saved detailed results to {results_file}")
    
    # Save to review-app database
    if db_path.exists():
        save_to_review_app(results, stats, db_path)
    else:
        print(f"\n⚠ Review-app database not found: {db_path}")
        print(f"  Skipping database save. Run review-app setup first.")
    
    print("\n✅ Telomer validation complete!")


if __name__ == '__main__':
    main()
