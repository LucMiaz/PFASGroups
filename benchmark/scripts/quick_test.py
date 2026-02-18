import sys
sys.path.append('c:/Users/luc/git/HalogenGroups')

print("Starting benchmark test...")

try:
    print("Importing modules...")
    from benchmark_pfas_definitions import PFASDefinitionBenchmark
    
    print("Creating benchmark instance...")
    benchmark = PFASDefinitionBenchmark()
    
    print("Running test...")
    result = benchmark.test_single_molecule('OC(=O)C(F)(F)C(F)(F)F')
    print(f"Test result: {result}")
    
    print("SUCCESS!")
    
except Exception as e:
    print(f"\nERROR: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
