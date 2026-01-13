"""Test SMARTS precomputation in PFASGroup"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from PFASgroups.core import parse_smiles

# Test with PFOA (Perfluorooctanoic acid)
pfoa_smiles = "C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F"

print("Testing PFOA with precomputed SMARTS atom counts...")
print(f"SMILES: {pfoa_smiles}")
print()

results = parse_smiles([pfoa_smiles], bycomponent=True)

if results:
    result = results[0]
    print(f"Formula: {result['formula']}")
    print(f"Total atoms: {result.get('n_atoms', 'N/A')}")
    print()
    
    if 'matches' in result:
        for match in result['matches']:
            print(f"Group: {match.get('group_name', 'Unknown')} (ID: {match.get('group_id', 'N/A')})")
            print(f"  Mean component fraction: {match.get('mean_component_fraction', 0):.3f}")
            print(f"  Total components fraction: {match.get('total_components_fraction', 0):.3f}")
            
            if 'components' in match:
                for i, comp in enumerate(match['components']):
                    print(f"  Component {i+1}:")
                    print(f"    Size (carbon backbone): {comp['size']}")
                    print(f"    Component fraction: {comp['component_fraction']:.3f}")
                    if comp.get('smarts_matches'):
                        print(f"    SMARTS matches: {len(comp['smarts_matches'])} atoms")
            print()
else:
    print("No results returned")
