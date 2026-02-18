print("Testing HalogenGroups imports...")
import sys
sys.path.append('c:/Users/luc/git/HalogenGroups')

try:
    print("1. Importing parse_mol...")
    from HalogenGroups.parser import parse_mol
    print("   Success!")
    
    print("2. Importing PFASDefinition...")
    from HalogenGroups.PFASDefinitionModel import PFASDefinition  
    print("   Success!")
    
    print("All imports successful!")
    
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
