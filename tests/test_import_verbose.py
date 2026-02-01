#!/usr/bin/env python3
import sys
import traceback

print("=== Testing Import ===", flush=True)

try:
    print("Attempting: from PFASgroups import parse_smiles", flush=True)
    from PFASgroups import parse_smiles
    print(f"SUCCESS! parse_smiles = {parse_smiles}", flush=True)
except Exception as e:
    print(f"FAILED with exception: {type(e).__name__}: {e}", flush=True)
    print("\nFull traceback:", flush=True)
    traceback.print_exc()
    sys.exit(1)

print("\n=== Testing Parse ===", flush=True)
try:
    smiles = 'FC(F)(F)C(=O)O'
    result = parse_smiles(smiles)
    print(f"SUCCESS! Parsed {smiles}, got {len(result)} results", flush=True)
except Exception as e:
    print(f"FAILED: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

print("\n=== All Tests Passed ===", flush=True)
