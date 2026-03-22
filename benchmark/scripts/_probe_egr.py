"""Probe PFASGroups parse_smiles batch mode."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')
from PFASGroups import parse_smiles

# Test batch mode with list
smis = ['FC(F)(F)C(=O)O', 'OCCO', 'CCCCCCCC(=O)O']
try:
    r = parse_smiles(smis, halogens='F', progress=False)
    arr = r.to_array(component_metrics=['effective_graph_resistance'])
    print('Batch mode OK: shape =', arr.shape, ' n_molecules =', r.n_molecules)
except Exception as e:
    print('Batch mode FAILED:', e)

# Test column_names (class-level or instance-level?)
r2 = parse_smiles('FC(F)(F)C(=O)O', halogens='F', progress=False)
names = r2.column_names(component_metrics=['binary'])
print('column_names(binary) count:', len(names), ' first:', names[0])
names_egr = r2.column_names(component_metrics=['effective_graph_resistance'])
print('column_names(EGR) count:', len(names_egr))

