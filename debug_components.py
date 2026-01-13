"""Debug fluorinated components and max_dist_from_CF"""
from rdkit import Chem
from PFASgroups.core import ComponentsSolver, mol_to_nx
from PFASgroups import parse_smiles
from PFASgroups.PFASGroupModel import PFASGroup
import json

# Load PFAS groups
with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    groups_data = json.load(f)
PFAS_GROUPS = [PFASGroup(**g) for g in groups_data]

# Test molecule 1: N-hydroxysuccinimide ester
smiles1 = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)CCC(=O)ON1C(=O)CCC1=O"
mol1 = Chem.MolFromSmiles(smiles1)
mol1 = Chem.AddHs(mol1)

print("="*80)
print("Test 1: N-hydroxysuccinimide ester")
print(f"SMILES: {smiles1}")
print("="*80)

# Print atom indices with their symbols
print("\nAtom indices:")
for i, atom in enumerate(mol1.GetAtoms()):
    symbol = atom.GetSymbol()
    neighbors = [mol1.GetAtomWithIdx(n).GetSymbol() for n in [x.GetIdx() for x in atom.GetNeighbors()]]
    print(f"  Atom {i}: {symbol} (neighbors: {neighbors})")

# Create components solver
with ComponentsSolver(mol1) as solver:
    print("\nFluorinated components:")
    for path_type, components in solver.components.items():
        print(f"\n  {path_type}:")
        for i, comp in enumerate(components):
            comp_list = sorted(list(comp))
            print(f"    Component {i}: {comp_list}")
            # Show symbols
            symbols = [mol1.GetAtomWithIdx(idx).GetSymbol() for idx in comp_list]
            print(f"      Symbols: {symbols}")
    
    # Test with max_dist=0
    print("\nExtended components (max_dist=0):")
    for path_type in solver.components.keys():
        ext_comps = solver.get(path_type, max_dist=0)
        print(f"\n  {path_type}:")
        for i, comp in enumerate(ext_comps):
            comp_list = sorted(list(comp))
            print(f"    Component {i}: {comp_list}")
    
    # Test with max_dist=1
    print("\nExtended components (max_dist=1):")
    for path_type in solver.components.keys():
        ext_comps = solver.get(path_type, max_dist=1)
        print(f"\n  {path_type}:")
        for i, comp in enumerate(ext_comps):
            comp_list = sorted(list(comp))
            print(f"    Component {i}: {comp_list}")

# Now test with actual PFAS group matching
print("\n\n" + "="*80)
print("Testing PFAS group matching:")
print("="*80)
results = parse_smiles([smiles1], bycomponent=True, output_format='dict')
if results:
    result = results[0]
    matches = result.get('matches', [])
    for match in matches:
        if isinstance(match, dict):
            group_name = match.get('group_name')
            if group_name in ['Heterocyclic azole', 'heterocyclic azole']:
                print(f"\n{group_name} match:")
                components = match.get('components', [])
                for comp in components:
                    comp_atoms = comp.get('component', [])
                    smarts_matches = comp.get('smarts_matches', [])
                    print(f"  Component atoms: {sorted(comp_atoms)}")
                    print(f"  SMARTS matches: {sorted(smarts_matches) if smarts_matches else None}")
                    
                # Check max_dist_from_CF setting
                pfas_group = next((g for g in PFAS_GROUPS if g.name == group_name), None)
                if pfas_group:
                    print(f"  max_dist_from_CF: {pfas_group.max_dist_from_CF}")

print("\n\n" + "="*80)
print("Test 2: PEG-fluorotelomer")
print("="*80)

smiles2 = "COCCOCCOCCOCCOCCOCCOCCOCCOCCOCC(O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F"
mol2 = Chem.MolFromSmiles(smiles2)
mol2 = Chem.AddHs(mol2)

print(f"SMILES: {smiles2}")
print(f"\nTotal atoms: {mol2.GetNumAtoms()}")

# Show which atoms are in PEG vs fluorinated region
print("\nPEG region (should be atoms 0-30):")
peg_atoms = []
for i in range(min(31, mol2.GetNumAtoms())):
    atom = mol2.GetAtomWithIdx(i)
    if atom.GetSymbol() in ['C', 'O']:
        peg_atoms.append(i)
        print(f"  Atom {i}: {atom.GetSymbol()}")

print("\nFluorinated region (should start around atom 32):")
fluoro_atoms = []
for i, atom in enumerate(mol2.GetAtoms()):
    if atom.GetSymbol() == 'F':
        fluoro_atoms.append(i)
        # Find the carbon it's connected to
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                print(f"  Atom {i}: F bonded to C{neighbor.GetIdx()}")
                break

# Check components
with ComponentsSolver(mol2) as solver:
    print("\nFluorinated components:")
    for path_type, components in solver.components.items():
        print(f"\n  {path_type}:")
        for i, comp in enumerate(components):
            comp_list = sorted(list(comp))
            print(f"    Component {i}: atoms {min(comp_list)}-{max(comp_list)} (total: {len(comp_list)})")

# Check PFAS matches
results2 = parse_smiles([smiles2], bycomponent=True, output_format='dict')
if results2:
    result = results2[0]
    matches = result.get('matches', [])
    for match in matches:
        if isinstance(match, dict):
            group_name = match.get('group_name')
            if group_name in ['alcohol', 'ether']:
                print(f"\n{group_name} match:")
                components = match.get('components', [])
                print(f"  Number of components: {len(components)}")
                if len(components) > 0:
                    first_comp = components[0]
                    comp_atoms = first_comp.get('component', [])
                    smarts_matches = first_comp.get('smarts_matches', [])
                    print(f"  First component atoms: {sorted(comp_atoms)[:10]}..." if len(comp_atoms) > 10 else f"  First component atoms: {sorted(comp_atoms)}")
                    print(f"  SMARTS matches: {sorted(smarts_matches)[:10]}..." if smarts_matches and len(smarts_matches) > 10 else f"  SMARTS matches: {sorted(smarts_matches) if smarts_matches else None}")
