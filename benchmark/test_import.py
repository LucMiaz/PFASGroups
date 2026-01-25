print("Starting test...")
import sys
sys.path.append('c:/Users/luc/git/PFASGroups')
print("Path appended")

try:
    print("Attempting import...")
    from benchmark_pfas_definitions import PFASDefinitionBenchmark
    print("Import successful")
    
    print("Creating benchmark instance...")
    benchmark = PFASDefinitionBenchmark()
    print("Benchmark initialized successfully")
    
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()

print("Test complete")
