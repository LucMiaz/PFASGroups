import sys
sys.path.append('c:/Users/luc/git/HalogenGroups')

from benchmark_pfas_definitions import PFASDefinitionBenchmark

print("Creating benchmark...")
benchmark = PFASDefinitionBenchmark()

print("\nTesting single molecule...")
result = benchmark.test_single_molecule('OC(=O)C(F)(F)C(F)(F)C(F)(F)F')
print(f"Valid: {result['valid']}")
print(f"Detected: {len(result['detected_ids'])} definitions")

print("\nTesting performance benchmark (5 molecules)...")
try:
    perf_result = benchmark.run_performance_benchmark(num_molecules=5)
    print(f"Performance test complete: {perf_result['summary']['total_molecules']} molecules tested")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
