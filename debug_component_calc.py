"""Debug component fraction calculation."""
from rdkit import Chem
from PFASgroups.core import parse_groups_in_mol
from PFASgroups.PFASGroupModel import PFASGroup
import json

# Load PFAS groups
with open('c:/Users/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    groups_data = json.load(f)

pfas_groups = [PFASGroup(**g) for g in groups_data]
group1 = pfas_groups[0]  # Perfluoroalkyl carboxylic acids

print(f"Group: {group1.name}")
print(f"SMARTS1: {group1.smarts1}")
print(f"smarts1_size: {group1.smarts1_size}")
print(f"smarts1_extra_atoms: {group1.smarts1_extra_atoms}")
print()

# Test molecule
smiles = "FC(F)(F)C(F)(F)C(F)(F)C(=O)O"
mol = Chem.MolFromSmiles(smiles)
print(f"Test molecule: {smiles}")
print(f"Total atoms: {mol.GetNumAtoms()}")
print(f"Total carbons: {sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')}")
print()

# Show atom indices
for atom in mol.GetAtoms():
    print(f"Atom {atom.GetIdx()}: {atom.GetSymbol()}")
print()

# Match SMARTS (group1.smarts1 is already a Mol object)
matches = mol.GetSubstructMatches(group1.smarts1)
print(f"SMARTS matches: {matches}")
print(f"SMARTS matches carbon at index: {matches[0][0]}")
print()

# Get components
results = parse_groups_in_mol(mol, pfas_groups=[group1])
print(f"Results type: {type(results)}")
print(f"Number of results: {len(results)}")
if results:
    print(f"First result type: {type(results[0])}")
    print(f"First result length: {len(results[0])}")
    pfas_group, match_count, component_sizes, matched_components = results[0]
    print(f"\nResults for {pfas_group.name}:")
    print(f"  Match count: {match_count}")
    print(f"  Component sizes: {component_sizes}")
    print(f"  Number of components: {len(matched_components)}")
    
    for i, comp in enumerate(matched_components):
        print(f"\n  Component {i+1}:")
        print(f"    Component atoms: {comp['component']}")
        print(f"    SMARTS atoms: {comp.get('smarts_matches', 'None')}")
        print(f"    Component fraction: {comp.get('component_fraction', 'N/A')}")
        
        # Manual calculation
        component_set = set(comp['component'])
        smarts_set = set(comp.get('smarts_matches', [])) if comp.get('smarts_matches') else set()
        
        component_carbons = sum(1 for idx in component_set if mol.GetAtomWithIdx(idx).GetSymbol() == 'C')
        smarts_carbons_in_component = sum(1 for idx in smarts_set if idx in component_set and mol.GetAtomWithIdx(idx).GetSymbol() == 'C')
        smarts_carbons_not_in_component = sum(1 for idx in smarts_set if idx not in component_set and mol.GetAtomWithIdx(idx).GetSymbol() == 'C')
        
        print(f"    Component carbons: {component_carbons}")
        print(f"    SMARTS carbons in component: {smarts_carbons_in_component}")
        print(f"    SMARTS carbons not in component: {smarts_carbons_not_in_component}")
        print(f"    smarts1_extra_atoms: {group1.smarts1_extra_atoms}")
        print(f"    adjusted_extra_atoms: {max(0, group1.smarts1_extra_atoms - smarts_carbons_in_component)}")
        print(f"    Total: {component_carbons} + {smarts_carbons_not_in_component} + {max(0, group1.smarts1_extra_atoms - smarts_carbons_in_component)}")
