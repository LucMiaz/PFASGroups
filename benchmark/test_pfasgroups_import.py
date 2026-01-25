print("Testing PFASgroups imports...")
import sys
sys.path.append('c:/Users/luc/git/PFASGroups')

try:
    print("1. Importing parse_mol...")
    from PFASgroups.parser import parse_mol
    print("   Success!")
    
    print("2. Importing PFASDefinition...")
    from PFASgroups.PFASDefinitionModel import PFASDefinition  
    print("   Success!")
    
    print("All imports successful!")
    
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
