"""Test component_fraction metric"""
from rdkit import Chem
from PFASgroups.core import parse_mols

# Test with a perfluorinated carboxylic acid
smiles = 'C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O'
mol = Chem.MolFromSmiles(smiles)
results = parse_mols([mol], bycomponent=True)

print(f"Testing SMILES: {smiles}")
print(f"Number of matches: {len(results[0]['matches'])}")

for match in results[0]['matches']:
    group_name = match.get('group_name', match.get('definition_name', 'Unknown'))
    print(f"\nPFAS Group: {group_name}")
    
    # Only process if it has components
    if 'components' not in match:
        continue
        
    print(f"Number of components: {len(match['components'])}")
    
    for i, comp in enumerate(match['components']):
        print(f"\n  Component {i+1}:")
        print(f"    Size: {comp['size']} atoms")
        print(f"    Component fraction: {comp['component_fraction']:.3f}")
        print(f"    Branching: {comp['branching']:.3f}")
        print(f"    Mean eccentricity: {comp['mean_eccentricity']:.3f}")
        print(f"    Median eccentricity: {comp['median_eccentricity']:.3f}")
        print(f"    SMARTS type: {comp['SMARTS']}")

print("\n✓ Test completed successfully!")
