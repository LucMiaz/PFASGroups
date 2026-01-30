#!/usr/bin/env python3
"""
Validate telomer detection on PubChem fluorotelomer dataset.
Run from benchmark/ directory.
"""

from pathlib import Path
import json
from datetime import datetime

from PFASgroups import parse_smiles

# Import Chem after PFASgroups to avoid conflicts
from rdkit import Chem

# Telomer group IDs (62-87)
TELOMER_GROUP_IDS = set(range(68, 93))

def main():
    """Main execution."""
    print("="*60)
    print("Telomer Detection Validation")
    print("="*60)
    
    # Input file
    sdf_file = Path('../data/PubChem_fluorotelomers.sdf')
    
    if not sdf_file.exists():
        print(f"Error: SDF file not found: {sdf_file}")
        return
    
    print(f"\nLoading molecules from: {sdf_file}")
    
    # Load SDF file
    supplier = Chem.SDMolSupplier(str(sdf_file))
    molecules = []
    
    for i, mol in enumerate(supplier):
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            molecules.append({
                'index': i,
                'smiles': smiles,
                'mol': mol
            })
    
    print(f"Loaded {len(molecules)} molecules")
    
    # Test detection
    print(f"\nTesting telomer detection...")
    results = []
    telomer_detected_count = 0
    
    for i, mol_data in enumerate(molecules, 1):
        smiles = mol_data['smiles']
        
        if i % 10 == 0:
            print(f"  [{i}/{len(molecules)}] processed...", end='\r')
        
        try:
            # Parse with PFASgroups
            parse_results = parse_smiles(smiles)
            
            if parse_results and len(parse_results) > 0:
                mol_result = parse_results[0]
                
                # Extract detected telomer groups
                detected_telomers = []
                
                if isinstance(mol_result, dict) and 'matches' in mol_result:
                    for match in mol_result['matches']:
                        if match.get('type') == 'PFASgroup':
                            group_id = match['id']
                            if group_id in TELOMER_GROUP_IDS:
                                detected_telomers.append({
                                    'id': group_id,
                                    'name': match['group_name']
                                })
                
                if detected_telomers:
                    telomer_detected_count += 1
                    results.append({
                        'index': mol_data['index'],
                        'smiles': smiles,
                        'detected': True,
                        'telomer_groups': detected_telomers
                    })
                else:
                    results.append({
                        'index': mol_data['index'],
                        'smiles': smiles,
                        'detected': False
                    })
            else:
                results.append({
                    'index': mol_data['index'],
                    'smiles': smiles,
                    'detected': False,
                    'error': 'Parse returned empty'
                })
        
        except Exception as e:
            results.append({
                'index': mol_data['index'],
                'smiles': smiles,
                'detected': False,
                'error': str(e)
            })
    
    print(f"\n\n{'='*60}")
    print("Results Summary")
    print("="*60)
    print(f"Total molecules: {len(molecules)}")
    print(f"Telomers detected: {telomer_detected_count}")
    print(f"Detection rate: {100*telomer_detected_count/len(molecules):.1f}%")
    print(f"Missed: {len(molecules) - telomer_detected_count}")
    
    # Count which groups were detected
    group_counts = {}
    for result in results:
        if result.get('detected') and 'telomer_groups' in result:
            for group in result['telomer_groups']:
                gid = group['id']
                if gid not in group_counts:
                    group_counts[gid] = {'id': gid, 'name': group['name'], 'count': 0}
                group_counts[gid]['count'] += 1
    
    if group_counts:
        print(f"\nMost frequently detected telomer groups:")
        for i, (gid, data) in enumerate(sorted(group_counts.items(), key=lambda x: x[1]['count'], reverse=True)[:10], 1):
            print(f"  {i:2d}. ID {data['id']:2d}: {data['name']:<45s} ({data['count']:3d} detections)")
    
    # Save results
    output_file = Path('../data/telomer_validation_results.json')
    output_data = {
        'test_date': datetime.now().isoformat(),
        'dataset': 'PubChem_fluorotelomers.sdf',
        'total_molecules': len(molecules),
        'telomer_detected': telomer_detected_count,
        'detection_rate': 100*telomer_detected_count/len(molecules),
        'group_counts': list(group_counts.values()),
        'results': results
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n✓ Saved results to {output_file}")
    print("\n✅ Validation complete!")

if __name__ == '__main__':
    main()
