import sys
import traceback

try:
    from PFASgroups.parser import parse_smiles
    print("SUCCESS: Import worked!")
except Exception as e:
    print(f"ERROR: {e}")
    traceback.print_exc()
