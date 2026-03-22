"""Quick import test for compare_fingerprints_toxcast."""
import importlib.util

spec = importlib.util.spec_from_file_location(
    'cmp', r'C:\Users\luc\git\PFASGroups\benchmark\scripts\compare_fingerprints_toxcast.py'
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
print('Import OK')
print('TXPP_TSV exists:', mod.TXPP_TSV.exists(), str(mod.TXPP_TSV))
print('TOXPRINT_TSV exists:', mod.TOXPRINT_TSV.exists())
print('DATASET_PKL exists:', mod.DATASET_PKL.exists())
print('EGR_CACHE path:', mod.EGR_CACHE)
print('FSET_COLORS_A:', list(mod.FSET_COLORS_A.keys()))
print('FSET_COLORS_B:', list(mod.FSET_COLORS_B.keys()))
