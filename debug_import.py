#!/usr/bin/env python3
import sys
import traceback

print("Step 1: Import networkx")
try:
    import networkx as nx
    print("  OK")
except Exception as e:
    print(f"  FAILED: {e}")
    traceback.print_exc()
    sys.exit(1)

print("Step 2: Import rdkit")
try:
    from rdkit import Chem
    print("  OK")
except Exception as e:
    print(f"  FAILED: {e}")
    traceback.print_exc()
    sys.exit(1)

print("Step 3: Import PFASgroups.core")
try:
    from PFASgroups import core
    print("  OK")
except Exception as e:
    print(f"  FAILED: {e}")
    traceback.print_exc()
    sys.exit(1)

print("Step 4: Import ComponentsSolverModel")
try:
    from PFASgroups import ComponentsSolverModel
    print("  OK")
except Exception as e:
    print(f"  FAILED: {e}")
    traceback.print_exc()
    sys.exit(1)

print("Step 5: Import parser")
try:
    from PFASgroups import parser
    print("  OK")
except Exception as e:
    print(f"  FAILED: {e}")
    traceback.print_exc()
    sys.exit(1)

print("Step 6: Import PFASgroups")
try:
    import PFASgroups
    print("  OK")
except Exception as e:
    print(f"  FAILED: {e}")
    traceback.print_exc()
    sys.exit(1)

print("\nAll imports successful!")
