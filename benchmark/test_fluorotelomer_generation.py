"""Quick test for fluorotelomer molecule generation"""

import sys
import os

# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

from enhanced_pfas_benchmark import EnhancedPFASBenchmark

# Create benchmark instance
benchmark = EnhancedPFASBenchmark()

print("Testing Fluorotelomer Molecule Generation")
print("=" * 60)

# Test a few fluorotelomer groups
test_groups = [62, 66, 68, 69]  # Silane, carboxylic acid, ethoxylate, sulfonic acid

for group_id in test_groups:
    group_info = benchmark.functional_smarts[group_id]
    print(f"\n📋 Group {group_id}: {group_info['name']}")
    print(f"   Telomer: {group_info.get('telomer', False)}")
    if group_info.get('telomer'):
        print(f"   CH2 range: {group_info.get('ch2_range', 'N/A')}")
        print(f"   Ethoxylate: {group_info.get('ethoxylate', False)}")
    
    # Generate 5 molecules
    molecules = benchmark.generate_fluorotelomer_molecules(group_id, count=5)
    
    print(f"\n   Generated {len(molecules)} molecules:")
    for i, mol_data in enumerate(molecules, 1):
        telomer_notation = mol_data.get('telomer_notation', 'N/A')
        smiles = mol_data['smiles']
        print(f"   {i}. {telomer_notation}: {smiles}")

print("\n" + "=" * 60)
print("✅ Test complete!")
