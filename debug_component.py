from rdkit import Chem
from PFASgroups.core import find_alkyl_components, get_fluorinated_subgraph, parse_groups_in_mol

# Test the specific molecule 
mol = Chem.MolFromSmiles('O=C(Cl)C(F)(F)C(F)(C(F)(F)C(F)(F)F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F')

print(f'Molecule: {Chem.MolToSmiles(mol)}')

# Get fluorinated components (this should use the decorator)
try:
    fluorinated_components_dict = get_fluorinated_subgraph(mol)
    print('\nFluorinated components found:')
    for key, components in fluorinated_components_dict.items():
        print(f'  {key}: {len(components)} components')
        if components:
            print(f'    Component 0 length: {len(components[0])}')
        
    # Test find_alkyl_components with the polyfluoroalkyl components
    polyfluoro_components = fluorinated_components_dict.get('Polyfluoroalkyl', [])
    print(f'\nUsing Polyfluoroalkyl components: {len(polyfluoro_components)} components')
    
    if polyfluoro_components:
        # Test with acyl halide SMARTS (group 35)
        acyl_halide_smarts = Chem.MolFromSmarts('[#6$([#6][#6](=O)[#9,#17,#35,#53])]')
        print(f'\nAcyl halide SMARTS matches: {len(mol.GetSubstructMatches(acyl_halide_smarts))}')
        
        result = find_alkyl_components(mol, acyl_halide_smarts, polyfluoro_components)
        print(f'find_alkyl_components result: {result}')
    else:
        print('No Polyfluoroalkyl components found!')
        
except Exception as e:
    print(f'Error: {e}')
    import traceback
    traceback.print_exc()

print('\n' + '='*50)
print('Now testing full parse_groups_in_mol with bycomponent=True')
try:
    matches = parse_groups_in_mol(mol, bycomponent=True)
    detected = [match[0].id for match in matches]
    print(f'Bycomponent detected: {detected}')
except Exception as e:
    print(f'Error in parse_groups_in_mol: {e}')
    import traceback
    traceback.print_exc()