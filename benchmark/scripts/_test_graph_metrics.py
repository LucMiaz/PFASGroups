"""Quick test of PFASGroups graph metrics API."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
print("Starting test...")
sys.stdout.flush()

from PFASGroups import generate_fingerprint, parse_smiles
import numpy as np

smiles_list = ['FC(F)(F)C(=O)O', 'C(C(F)(F)F)F', 'CC(=O)O']

GRAPH_METRICS = ['binary', 'total_component', 'effective_graph_resistance', 'min_dist_to_barycenter']
ALL_MOL_METRICS = ['n_components', 'total_size', 'mean_size', 'max_size',
                   'mean_branching', 'max_branching', 'mean_eccentricity',
                   'max_diameter', 'mean_component_fraction', 'max_component_fraction']

# 1. generate_fingerprint with graph metrics only (single SMILES → probe column names)
fp_probe, col_graph = generate_fingerprint(
    smiles_list[0],
    representation='vector',
    component_metrics=GRAPH_METRICS,
    halogens='F',
    saturation=None,
)
print(f"generate_fingerprint graph: shape={fp_probe.shape}, n_cols={len(col_graph)}")
print(f"  first 3 cols: {col_graph[:3]}")

# 2. parse_smiles + to_array — no mol metrics
r = parse_smiles(smiles_list)
arr_no_mol = r.to_array(component_metrics=GRAPH_METRICS)
print(f"parse_smiles to_array no mol: shape={arr_no_mol.shape}")
print(f"  as numpy dtype: {np.array(arr_no_mol).dtype}")

# 3. parse_smiles + to_array — with mol metrics
arr_mol = r.to_array(component_metrics=GRAPH_METRICS, molecule_metrics=ALL_MOL_METRICS)
print(f"parse_smiles to_array with mol: shape={arr_mol.shape}")
# extra columns
print(f"  extra cols: {arr_mol.shape[1] - arr_no_mol.shape[1]}")
# values of mol metrics columns
mol_matrix = np.array(arr_mol)
print(f"  mol-metric values[0]: {mol_matrix[0, -10:]}")

# 4. Check nan handling
print(f"  nan count (no mol): {np.isnan(np.array(arr_no_mol).astype(float)).sum()}")
print(f"  nan count (mol): {np.isnan(mol_matrix.astype(float)).sum()}")
