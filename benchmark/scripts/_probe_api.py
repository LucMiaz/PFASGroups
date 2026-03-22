"""Probe PFASGroups API to understand group count and available options."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')
import traceback

print("=== generate_fingerprint ===")
try:
    from PFASGroups import generate_fingerprint
    import inspect
    sig = inspect.signature(generate_fingerprint)
    print("Parameters:", list(sig.parameters.keys()))
    for name, p in sig.parameters.items():
        default = p.default if p.default is not inspect.Parameter.empty else "(required)"
        print(f"  {name}: default={default}")
except Exception as e:
    traceback.print_exc()

print()
print("=== generate_fingerprint binary ===")
try:
    fp, info = generate_fingerprint(
        'FC(F)(F)C(=O)O',
        representation='vector',
        count_mode='binary',
        halogens='F',
        saturation=None,
    )
    print("n_groups:", len(fp))
    print("info keys:", list(info.keys()))
    print("group_names:", info.get('group_names', '?'))
except Exception as e:
    traceback.print_exc()

print()
print("=== parse_smiles + to_array ===")
try:
    from PFASGroups import parse_smiles
    result = parse_smiles(['FC(F)(F)C(=O)O', 'OC(=O)C(F)(F)C(F)(F)F'])
    print("parse_smiles type:", type(result))
    print("parse_smiles methods:", [m for m in dir(result) if not m.startswith('_')])
except Exception as e:
    traceback.print_exc()

print()
print("=== to_array with component_metrics ===")
try:
    from PFASGroups import parse_smiles
    import inspect
    result = parse_smiles(['FC(F)(F)C(=O)O', 'OC(=O)C(F)(F)C(F)(F)F'])
    # Inspect to_array signature
    sig_arr = inspect.signature(result.to_array)
    print("to_array signature:", sig_arr)
    # Inspect EmbeddingArray methods
    arr_b = result.to_array(component_metrics=['binary'])
    print("EmbeddingArray type:", type(arr_b))
    print("EmbeddingArray methods:", [m for m in dir(arr_b) if not m.startswith('_')])
    # Try to get shape
    import numpy as np
    try:
        nparr = np.array(arr_b)
        print("np.array shape:", nparr.shape)
    except Exception as e2:
        print("np.array error:", e2)
    try:
        print("shape attr:", arr_b.shape)  # type: ignore
    except Exception as e2:
        print("no .shape attr:", e2)
    try:
        print("columns attr:", list(arr_b.columns)[:10])  # type: ignore
    except Exception as e2:
        print("no .columns attr:", e2)
    try:
        df = arr_b.as_dataframe()  # type: ignore
        print("as_dataframe shape:", df.shape)
        print("columns (first 5):", list(df.columns[:5]))
        print("columns (last 5):", list(df.columns[-5:]))
    except Exception as e2:
        print("as_dataframe error:", e2)
    try:
        df2 = arr_b.to_dataframe()  # type: ignore
        print("to_dataframe shape:", df2.shape)
    except Exception as e2:
        print("to_dataframe error:", e2)
except Exception as e:
    traceback.print_exc()

print()
print("=== to_array with all component_metrics ===")
try:
    result2 = parse_smiles(['FC(F)(F)C(=O)O', 'OC(=O)C(F)(F)C(F)(F)F'])
    arr_all = result2.to_array(component_metrics=['binary', 'total_component', 'effective_graph_resistance', 'min_dist_to_barycenter'])
    nparr2 = np.array(arr_all)
    print("binary+graph shape:", nparr2.shape)
except Exception as e:
    traceback.print_exc()

print()
print("=== column_names ===")
try:
    result3 = parse_smiles(['FC(F)(F)C(=O)O', 'OC(=O)C(F)(F)C(F)(F)F'])
    cols = result3.column_names(component_metrics=['binary'])
    print("column_names binary type:", type(cols))
    print("n cols binary:", len(cols))
    print("first 5:", cols[:5])
    print("last 5:", cols[-5:])
    cols_all = result3.column_names(component_metrics=['binary', 'total_component', 'effective_graph_resistance', 'min_dist_to_barycenter'])
    print("n cols binary+graph:", len(cols_all))
except Exception as e:
    traceback.print_exc()
