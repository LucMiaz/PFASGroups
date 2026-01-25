import sys
print("Test 1: Basic imports")
import json
import time
import os
print("  Basic imports OK")

print("Test 2: Data imports")
import pandas as pd
import numpy as np
print("  Data imports OK")

print("Test 3: RDKit imports")
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
print("  RDKit imports OK")

print("Test 4: Path setup")
script_dir = 'c:/Users/luc/git/PFASGroups/benchmark'
parent_dir = 'c:/Users/luc/git/PFASGroups'
sys.path.append(parent_dir)
print(f"  Path: {parent_dir}")

print("Test 5: PFASgroups.core import")
try:
    from PFASgroups.parser import parse_mol
    print("  parse_mol imported OK")
except Exception as e:
    print(f"  FAILED: {e}")
    import traceback
    traceback.print_exc()

print("Test 6: PFASDefinitionModel import")
try:
    from PFASgroups.PFASDefinitionModel import PFASDefinition
    print("  PFASDefinition imported OK")
except Exception as e:
    print(f"  FAILED: {e}")
    import traceback
    traceback.print_exc()

print("\nAll tests complete!")
