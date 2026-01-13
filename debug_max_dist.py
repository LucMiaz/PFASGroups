"""Debug max_dist_from_CF issues"""
from rdkit import Chem
from PFASgroups import parse_smiles

# Test molecules
test_cases = [
    {
        "name": "N-hydroxysuccinimide ester with perfluoro chain",
        "smiles": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)CCC(=O)ON1C(=O)CCC1=O",
        "issue": "Heterocyclic acid detected incorrectly"
    },
    {
        "name": "PEG-fluorotelomer",
        "smiles": "COCCOCCOCCOCCOCCOCCOCCOCCOCCOCC(O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F",
        "issue": "Alcohol and ether detected in non-fluorinated region"
    }
]

for test in test_cases:
    print(f"\n{'='*80}")
    print(f"Test: {test['name']}")
    print(f"Issue: {test['issue']}")
    print(f"SMILES: {test['smiles']}")
    print(f"{'='*80}")
    
    mol = Chem.MolFromSmiles(test['smiles'])
    if mol is None:
        print("ERROR: Could not parse SMILES")
        continue
    
    # Count carbons and fluorines
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    num_fluorines = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F')
    print(f"\nMolecule stats: {num_carbons} carbons, {num_fluorines} fluorines")
    
    # Parse with current settings
    results = parse_smiles([test['smiles']], bycomponent=True, output_format='dict')
    
    if results and len(results) > 0:
        result = results[0]
        matches = result.get('matches', [])
        
        print(f"\nFound {len(matches)} group matches:")
        for match in matches:
            if isinstance(match, dict):
                group_name = match.get('group_name', 'Unknown')
                match_count = match.get('match_count', 0)
                
                if match_count > 0:
                    print(f"\n  Group: {group_name}")
                    print(f"    Match count: {match_count}")
                    
                    # Check components
                    components = match.get('components', [])
                    if components:
                        for i, comp in enumerate(components):
                            comp_atoms = comp.get('component', [])
                            smarts_matches = comp.get('smarts_matches', [])
                            
                            print(f"    Component {i+1}:")
                            print(f"      Component atoms: {comp_atoms}")
                            print(f"      SMARTS matches: {smarts_matches}")
                            
                            # Show which atoms are matched
                            if smarts_matches:
                                for atom_idx in smarts_matches:
                                    atom = mol.GetAtomWithIdx(atom_idx)
                                    neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                                    print(f"        Atom {atom_idx}: {atom.GetSymbol()} (neighbors: {neighbors})")
    else:
        print("No matches found")

print("\n" + "="*80)
print("Checking PFAS group definitions for max_dist_from_CF settings...")
print("="*80)

import json
from PFASgroups import PFAS_GROUPS_FILE

with open(PFAS_GROUPS_FILE, 'r') as f:
    groups = json.load(f)

# Find relevant groups
relevant_groups = ['heterocyclic acid', 'alcohol', 'ether']
for group in groups:
    name = group.get('name', '').lower()
    if any(keyword in name for keyword in relevant_groups):
        print(f"\nGroup: {group.get('name')}")
        print(f"  ID: {group.get('id')}")
        print(f"  max_dist_from_CF: {group.get('max_dist_from_CF', 'not set')}")
        print(f"  smartsPath: {group.get('smartsPath')}")
        print(f"  smarts1: {group.get('smarts1', 'none')}")
