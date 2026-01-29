#!/usr/bin/env python3
print("START")

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

print("About to import...")
from PFASgroups import parse_smiles
print("Import OK!")
