"""Trace through fluorotelomer acrylate detection step by step"""
from rdkit import Chem
import json
from PFASgroups.ComponentsSolverModel import ComponentsSolver
from PFASgroups.PFASGroupModel import PFASGroup
import networkx as nx

# Test molecule
smiles = "C=CC(=O)OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
mol = Chem.MolFromSmiles(smiles)

print("="*80)
print("FLUOROTELOMER ACRYLATE DETECTION TRACE")
print("="*80)
print(f"SMILES: {smiles}")
print(f"Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Total atoms: {mol.GetNumAtoms()}\n")

# Load groups directly from JSON
with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    groups_data = json.load(f)

pfas_groups = [PFASGroup(**x) for x in groups_data if x.get('compute', True)]
group_79 = [g for g in pfas_groups if g.id == 79][0]

print("GROUP 79: Fluorotelomer acrylate")
print(f"  smartsPath: {group_79.smartsPath}")
print(f"  linker_smarts: {group_79.linker_smarts}")
print(f"  max_dist_from_CF: {group_79.max_dist_from_CF}")
print(f"  smarts: {group_79.smarts_str}")
print(f"  smarts_count: {group_79.smarts_count}")
print(f"  constraints: {group_79.constraints}")

print("\n" + "="*80)
print("STEP 1: SMARTS MATCHING")
print("="*80)

# Find SMARTS matches
matched = group_79.find_matched_atoms(mol)
print(f"find_matched_atoms result: {matched}")
print(f"SMARTS subset atoms: {sorted(group_79.subset)}")
print(f"all_matches: {[list(m) for m in group_79.all_matches]}")

if group_79.subset:
    for atom_idx in sorted(group_79.subset):
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        print(f"  Atom {atom_idx}: {atom.GetSymbol()} (degree={atom.GetDegree()}, neighbors={neighbors})")

print("\n" + "="*80)
print("STEP 2: COMPONENT DETECTION")
print("="*80)

# Create component solver
component_solver = ComponentsSolver(mol)

# Original Perfluoroalkyl components
orig_comps = component_solver.components.get('Perfluoroalkyl', [])
print(f"Original Perfluoroalkyl components: {len(orig_comps)}")
for i, comp in enumerate(orig_comps):
    print(f"  Component {i}: {sorted(comp)}")
    atoms_str = ','.join([mol.GetAtomWithIdx(a).GetSymbol() for a in sorted(comp)])
    print(f"    Atoms: {atoms_str}")

# Extended components at max_dist=16
component_solver.extend_components(16)
ext_comps = component_solver.extended_components.get('Perfluoroalkyl', {}).get(16, [])
print(f"\nExtended Perfluoroalkyl components (max_dist=16): {len(ext_comps)}")
for i, comp in enumerate(ext_comps):
    print(f"  Extended component {i}: size={len(comp)}, atoms={sorted(comp)[:10]}{'...' if len(comp) > 10 else ''}")

print("\n" + "="*80)
print("STEP 3: COMPONENT SATISFIES SMARTS CHECK")
print("="*80)

# Reset the extra atoms list before checking
group_79.component_specific_extra_atoms = []

for i, comp in enumerate(ext_comps):
    print(f"\nExtended component {i}:")
    
    # Manually check without side effects
    satisfies = True
    for j, (matches, min_count) in enumerate(zip(group_79.all_matches, group_79.smarts_count)):
        component_set = set(comp)
        found = sum(1 for match_tuple in matches if any(atom_idx in component_set for atom_idx in match_tuple))
        if found < min_count:
            satisfies = False
            break
    
    print(f"  component_satisfies_all_smarts: {satisfies}")
    
    if satisfies:
        # Check which SMARTS atoms are in this component
        smarts_in_comp = [atom for atom in group_79.subset if atom in comp]
        print(f"  SMARTS atoms in component: {smarts_in_comp}")

print("\n" + "="*80)
print("STEP 4: LINKER VALIDATION & AUGMENTATION")
print("="*80)

# Check linker matches
linker_smarts = group_79.linker_smarts
linker_matches = set([y for x in mol.GetSubstructMatches(linker_smarts) for y in x])
print(f"Linker SMARTS '[CH2X4]' matches {len(linker_matches)} atoms: {sorted(linker_matches)}")

for i, comp in enumerate(ext_comps):
    # Manually check if component satisfies SMARTS
    satisfies = True
    for j, (matches, min_count) in enumerate(zip(group_79.all_matches, group_79.smarts_count)):
        component_set = set(comp)
        found = sum(1 for match_tuple in matches if any(atom_idx in component_set for atom_idx in match_tuple))
        if found < min_count:
            satisfies = False
            break
    
    if not satisfies:
        continue
        
    print(f"\n--- Processing extended component {i} ---")
    orig_index = component_solver.component_to_original_index.get(('Perfluoroalkyl', 16, i), i)
    orig_comp = component_solver.components['Perfluoroalkyl'][orig_index]
    
    print(f"Original component: {sorted(orig_comp)}")
    print(f"Extended component size: {len(comp)}")
    
    # Manually trace augmentation
    augmented = set(orig_comp)
    print(f"\nStarting augmented component: {sorted(augmented)}")
    
    for smarts_atom in group_79.subset:
        print(f"\n  Checking SMARTS atom {smarts_atom}:")
        print(f"    In extended component: {smarts_atom in comp}")
        print(f"    In original component: {smarts_atom in orig_comp}")
        
        if smarts_atom in comp and smarts_atom not in orig_comp:
            try:
                # Find path
                G = component_solver.G
                path = nx.shortest_path(G, smarts_atom, list(orig_comp)[0])
                print(f"    Raw path: {path}")
                
                # Process path like shortest_path_to_component does
                border = path.pop(-1)
                while len(path) > 0 and border not in orig_comp:
                    border = path.pop(-1)
                
                if len(path) == 0:
                    print(f"    -> No valid path found")
                    continue
                    
                base_atom = border
                linker_atoms = [x for x in path if x != smarts_atom]
                
                print(f"    base_atom (component border): {base_atom}")
                print(f"    linker_atoms (intermediate): {linker_atoms}")
                
                # Check linker validation
                if len(linker_atoms) > 0:
                    validation_results = [atom in linker_matches for atom in linker_atoms]
                    print(f"    Linker validation: {validation_results}")
                    all_valid = all(atom in linker_matches for atom in linker_atoms)
                    print(f"    All linker atoms valid: {all_valid}")
                    
                    if not all_valid:
                        print(f"    -> REJECTED by linker validation")
                        continue
                else:
                    print(f"    -> Direct connection (no linker atoms to validate)")
                
                # Would add to augmented
                to_add = [smarts_atom] + linker_atoms + [base_atom]
                print(f"    -> Would add to augmented: {to_add}")
                augmented.update(to_add)
                
            except nx.NetworkXNoPath:
                print(f"    -> No path exists")
                continue
                
        elif smarts_atom in orig_comp:
            print(f"    -> Already in original component, adding to augmented")
            augmented.add(smarts_atom)
    
    print(f"\n  Final augmented component: size={len(augmented)}, atoms={sorted(augmented)}")
    
    # Check if any SMARTS made it in
    smarts_in_augmented = sum(1 for atom in group_79.subset if atom in augmented)
    print(f"  SMARTS atoms in augmented: {smarts_in_augmented}/{len(group_79.subset)}")
    
    if smarts_in_augmented == 0 and len(group_79.subset) > 0:
        print(f"  -> REJECTED: No SMARTS atoms in augmented component")
    else:
        print(f"  -> ACCEPTED: Augmented component is valid")

print("\n" + "="*80)
print("STEP 5: FINAL RESULT")
print("="*80)

# Run actual detection
from PFASgroups import parse_smiles
result = parse_smiles(smiles)

if result:
    # Handle dict result format with 'matches' key
    if isinstance(result, dict) and 'matches' in result:
        matches = result['matches']
        # Filter for PFASgroup type matches
        pfas_matches = [m for m in matches if m.get('type') == 'PFASgroup']
        group_ids = [m['id'] for m in pfas_matches]
        print(f"Total matches: {len(matches)}")
        print(f"PFAS group matches: {len(pfas_matches)}")
        if pfas_matches:
            print(f"Groups found:")
            for m in pfas_matches:
                print(f"  - Group {m['id']}: {m.get('group_name', 'N/A')}")
    # Handle list result format
    elif isinstance(result, list):
        print(f"Result is a list with {len(result)} items")
        if len(result) > 0:
            print(f"First item type: {type(result[0])}")
            if isinstance(result[0], dict):
                print(f"First item keys: {result[0].keys()}")
                # Check if it's the old format with 'pfasgroup' key
                if 'pfasgroup' in result[0]:
                    group_ids = [g['pfasgroup'].id for g in result]
                    print(f"Groups found:")
                    for g in result:
                        print(f"  - Group {g['pfasgroup'].id}: {g['pfasgroup'].name}")
                else:
                    print(f"Unknown dict format in list")
                    group_ids = []
            else:
                print(f"Unknown item type in list")
                group_ids = []
        else:
            print("Empty list")
            group_ids = []
    else:
        print(f"Unknown result format: {type(result)}")
        group_ids = []
else:
    print("No result returned")
    group_ids = []

print(f"\nDetected group IDs: {group_ids}")
print(f"Group 79 detected: {79 in group_ids}")
