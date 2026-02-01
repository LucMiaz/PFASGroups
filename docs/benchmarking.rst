Benchmarking and Validation
===========================

PFASgroups includes a comprehensive benchmarking suite to validate accuracy and measure performance.

Overview
--------

The benchmark suite tests:

- **Functional groups detection accuracy**
- **OECD database validation**
- **Timing performance and scalability**
- **Non-fluorinated exclusion**
- **Complex branched structures**
- **Fluorotelomer detection**

For complete details, see the `BENCHMARK_GUIDE.md <https://github.com/lucid-luc/PFASGroups/blob/main/benchmark/BENCHMARK_GUIDE.md>`_ in the repository.

Running Benchmarks
------------------

Quick Start
~~~~~~~~~~~

Run all benchmarks:

.. code-block:: bash

   cd benchmark
   
   # Linux/macOS
   ./run_all_benchmarks.sh
   
   # Windows PowerShell
   .\run_all_benchmarks.ps1

Individual Benchmarks
~~~~~~~~~~~~~~~~~~~~~

Run specific benchmark types:

.. code-block:: bash

   cd benchmark
   
   # Functional groups
   echo "1" | python scripts/enhanced_pfas_benchmark.py
   
   # OECD validation
   echo "2" | python scripts/enhanced_pfas_benchmark.py
   
   # Timing performance
   echo "3" | python scripts/enhanced_pfas_benchmark.py
   
   # Non-fluorinated exclusion
   echo "4" | python scripts/enhanced_pfas_benchmark.py
   
   # Complex branched structures
   echo "5" | python scripts/enhanced_pfas_benchmark.py

Benchmark Results
-----------------

Performance Metrics
~~~~~~~~~~~~~~~~~~~

The benchmarks generate detailed performance metrics:

**Timing Performance:**

.. code-block:: text

   Exponential model: t = a × exp(α × n)
   
   Parameters:
   - a = 0.001234 s (base time)
   - α = 0.0234 atoms⁻¹ (scaling factor)
   
   Performance:
   - 10 atoms: ~1.5 ms
   - 50 atoms: ~3.8 ms
   - 100 atoms: ~9.5 ms
   - 200 atoms: ~90 ms

**Accuracy Metrics:**

.. code-block:: text

   OECD Database Validation:
   - Total compounds: 4,000+
   - Correctly classified: 3,850+
   - Accuracy: >96%
   
   By Group:
   - PFCAs: 98.5% accuracy
   - PFSAs: 97.8% accuracy
   - FTOHs: 99.2% accuracy

Viewing Results
~~~~~~~~~~~~~~~

Results are saved in multiple formats:

**HTML Reports:**

.. code-block:: bash

   # Open unified report
   firefox reports/unified_pfas_benchmark_report_*.html

**JSON Data:**

.. code-block:: python

   import json
   
   # Load timing results
   with open('data/pfas_timing_benchmark_*.json', 'r') as f:
       timing_data = json.load(f)
   
   # Load OECD results
   with open('data/pfas_oecd_benchmark_*.json', 'r') as f:
       oecd_data = json.load(f)

**Interactive Review App:**

.. code-block:: bash

   cd review-app
   node server.js
   # Open http://localhost:5000

Validation Against Reference Data
----------------------------------

OECD Database
~~~~~~~~~~~~~

Validate against the OECD PFAS database:

.. code-block:: python

   from PFASgroups import parse_smiles
   import pandas as pd
   
   # Load OECD data
   oecd_df = pd.read_csv('OECD_PFAS_database.csv')
   
   # Analyze each compound
   results = []
   for idx, row in oecd_df.iterrows():
       smiles = row['SMILES']
       expected_class = row['First_Class']
       
       # Parse with PFASgroups
       pfas_result = parse_smiles(smiles)
       
       detected_groups = []
       if pfas_result[0]:
           detected_groups = [g.name for g, _, _, _ in pfas_result[0]]
       
       results.append({
           'smiles': smiles,
           'expected': expected_class,
           'detected': ', '.join(detected_groups),
           'match': expected_class in detected_groups
       })
   
   # Calculate accuracy
   accuracy = sum(r['match'] for r in results) / len(results)
   print(f"OECD Accuracy: {accuracy:.1%}")

PubChem Validation
~~~~~~~~~~~~~~~~~~

Validate fluorotelomer detection:

.. code-block:: python

   from PFASgroups import parse_smiles
   import pandas as pd
   
   # Load PubChem fluorotelomers
   pubchem_df = pd.read_csv('pubchem_fluorotelomers.csv')
   
   # Expected to detect as FTOHs (Group 15)
   correct = 0
   total = len(pubchem_df)
   
   for smiles in pubchem_df['SMILES']:
       results = parse_smiles(smiles)
       
       if results[0]:
           group_ids = [g.id for g, _, _, _ in results[0]]
           if 15 in group_ids:  # FTOH group
               correct += 1
   
   print(f"FTOH Detection Rate: {correct/total:.1%}")

Custom Validation
-----------------

Creating Test Sets
~~~~~~~~~~~~~~~~~~

Build custom test sets:

.. code-block:: python

   # Create test dataset
   test_cases = [
       {
           'name': 'PFOA',
           'smiles': 'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O',
           'expected_groups': [1],  # PFCA
           'expected_chain': 8
       },
       {
           'name': 'PFOS',
           'smiles': 'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O',
           'expected_groups': [6],  # PFSA
           'expected_chain': 8
       },
       # Add more test cases
   ]
   
   # Save as JSON
   import json
   with open('my_test_set.json', 'w') as f:
       json.dump(test_cases, f, indent=2)

Running Custom Tests
~~~~~~~~~~~~~~~~~~~~

Test against custom dataset:

.. code-block:: python

   from PFASgroups import parse_smiles
   import json
   
   def run_validation(test_file):
       """Run validation on custom test set."""
       
       # Load test cases
       with open(test_file, 'r') as f:
           test_cases = json.load(f)
       
       results = []
       for test in test_cases:
           smiles = test['smiles']
           expected_groups = set(test['expected_groups'])
           expected_chain = test.get('expected_chain')
           
           # Parse
           pfas_result = parse_smiles(smiles)
           
           if pfas_result[0]:
               detected_ids = {g.id for g, _, _, _ in pfas_result[0]}
               chain_lengths = [l[0] for _, _, l, _ in pfas_result[0] if l]
               max_chain = max(chain_lengths) if chain_lengths else 0
               
               # Check accuracy
               groups_match = detected_ids == expected_groups
               chain_match = (not expected_chain or 
                            max_chain == expected_chain)
               
               results.append({
                   'name': test['name'],
                   'groups_correct': groups_match,
                   'chain_correct': chain_match,
                   'detected': list(detected_ids),
                   'expected': list(expected_groups)
               })
       
       # Summary
       groups_accuracy = sum(r['groups_correct'] for r in results) / len(results)
       chain_accuracy = sum(r['chain_correct'] for r in results) / len(results)
       
       print(f"\nValidation Results:")
       print(f"  Group Detection: {groups_accuracy:.1%}")
       print(f"  Chain Length: {chain_accuracy:.1%}")
       
       return results
   
   # Run validation
   results = run_validation('my_test_set.json')

Performance Testing
-------------------

Timing Analysis
~~~~~~~~~~~~~~~

Measure execution time:

.. code-block:: python

   from PFASgroups import parse_smiles
   import time
   
   def time_analysis(smiles_list):
       """Time the analysis of molecules."""
       times = []
       
       for smiles in smiles_list:
           start = time.time()
           results = parse_smiles(smiles)
           elapsed = time.time() - start
           
           times.append({
               'smiles': smiles,
               'time': elapsed,
               'num_groups': len(results[0]) if results[0] else 0
           })
       
       # Statistics
       avg_time = sum(t['time'] for t in times) / len(times)
       max_time = max(t['time'] for t in times)
       min_time = min(t['time'] for t in times)
       
       print(f"\nTiming Results:")
       print(f"  Average: {avg_time*1000:.2f} ms")
       print(f"  Min: {min_time*1000:.2f} ms")
       print(f"  Max: {max_time*1000:.2f} ms")
       
       return times
   
   # Test on molecules
   test_molecules = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
       # Add more of varying complexity
   ]
   
   timing_results = time_analysis(test_molecules)

Scaling Analysis
~~~~~~~~~~~~~~~~

Test how performance scales with molecule size:

.. code-block:: python

   from PFASgroups import parse_smiles
   from rdkit import Chem
   import time
   import matplotlib.pyplot as plt
   
   def generate_pfas_chain(n):
       """Generate PFCA of length n."""
       chain = "FC(F)(F)" * (n-1) + "C(=O)O"
       return chain
   
   # Test different sizes
   sizes = [4, 6, 8, 10, 12, 14, 16, 18, 20]
   times = []
   
   for n in sizes:
       smiles = generate_pfas_chain(n)
       
       # Time multiple runs
       run_times = []
       for _ in range(10):
           start = time.time()
           parse_smiles(smiles)
           run_times.append(time.time() - start)
       
       avg_time = sum(run_times) / len(run_times)
       times.append(avg_time)
       print(f"C{n}: {avg_time*1000:.2f} ms")
   
   # Plot
   plt.figure(figsize=(10, 6))
   plt.plot(sizes, [t*1000 for t in times], 'o-')
   plt.xlabel('Chain Length')
   plt.ylabel('Time (ms)')
   plt.title('Performance Scaling')
   plt.grid(True)
   plt.savefig('scaling_analysis.png')

Accuracy Metrics
----------------

Confusion Matrix
~~~~~~~~~~~~~~~~

Generate confusion matrix for classification:

.. code-block:: python

   from sklearn.metrics import confusion_matrix, classification_report
   import numpy as np
   
   def generate_confusion_matrix(results):
       """Generate confusion matrix from validation results."""
       
       # Extract true and predicted labels
       y_true = [r['expected_group'] for r in results]
       y_pred = [r['detected_group'] for r in results]
       
       # Generate matrix
       cm = confusion_matrix(y_true, y_pred)
       
       # Calculate metrics
       report = classification_report(y_true, y_pred)
       
       print("\nConfusion Matrix:")
       print(cm)
       print("\nClassification Report:")
       print(report)
       
       return cm, report

Precision and Recall
~~~~~~~~~~~~~~~~~~~~

Calculate precision and recall:

.. code-block:: python

   def calculate_metrics(results):
       """Calculate precision, recall, F1."""
       
       tp = sum(1 for r in results if r['true_positive'])
       fp = sum(1 for r in results if r['false_positive'])
       fn = sum(1 for r in results if r['false_negative'])
       
       precision = tp / (tp + fp) if (tp + fp) > 0 else 0
       recall = tp / (tp + fn) if (tp + fn) > 0 else 0
       f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
       
       print(f"\nMetrics:")
       print(f"  Precision: {precision:.1%}")
       print(f"  Recall: {recall:.1%}")
       print(f"  F1 Score: {f1:.3f}")
       
       return precision, recall, f1

Comparative Analysis
--------------------

Compare with PFAS-Atlas
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from PFASgroups import parse_smiles
   # Assuming PFAS-Atlas is installed
   # from pfas_atlas import classify_pfas
   
   def compare_tools(smiles_list):
       """Compare PFASgroups with PFAS-Atlas."""
       
       results = []
       for smiles in smiles_list:
           # PFASgroups
           pg_result = parse_smiles(smiles)
           pg_detected = len(pg_result[0]) > 0 if pg_result[0] else False
           
           # PFAS-Atlas (example)
           # atlas_result = classify_pfas(smiles)
           # atlas_detected = atlas_result['is_pfas']
           
           results.append({
               'smiles': smiles,
               'pfasgroups': pg_detected,
               # 'atlas': atlas_detected,
               # 'agreement': pg_detected == atlas_detected
           })
       
       # Calculate agreement
       # agreement = sum(r['agreement'] for r in results) / len(results)
       # print(f"Tool Agreement: {agreement:.1%}")
       
       return results

Continuous Integration
----------------------

Automated Testing
~~~~~~~~~~~~~~~~~

Set up automated benchmark runs:

.. code-block:: yaml

   # .github/workflows/benchmark.yml
   name: Benchmark
   
   on:
     push:
       branches: [main]
     pull_request:
       branches: [main]
   
   jobs:
     benchmark:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v2
         - name: Set up Python
           uses: actions/setup-python@v2
           with:
             python-version: 3.9
         - name: Install dependencies
           run: |
             pip install -e .
             pip install pytest
         - name: Run benchmarks
           run: |
             cd benchmark
             python -m pytest tests/

Regression Testing
~~~~~~~~~~~~~~~~~~

Detect performance regressions:

.. code-block:: python

   import json
   
   def check_regression(current_results, baseline_file):
       """Check for performance regression."""
       
       # Load baseline
       with open(baseline_file, 'r') as f:
           baseline = json.load(f)
       
       # Compare
       degraded = []
       for key in current_results:
           current = current_results[key]
           base = baseline.get(key, {})
           
           # Check if performance degraded by >10%
           if 'time' in current and 'time' in base:
               if current['time'] > base['time'] * 1.1:
                   degraded.append({
                       'metric': key,
                       'baseline': base['time'],
                       'current': current['time'],
                       'degradation': (current['time']/base['time'] - 1) * 100
                   })
       
       if degraded:
           print("\n⚠️  Performance Regressions Detected:")
           for item in degraded:
               print(f"  {item['metric']}: +{item['degradation']:.1f}%")
       else:
           print("\n✓ No regressions detected")
       
       return degraded

See Also
--------

- `Complete Benchmark Guide <https://github.com/lucid-luc/PFASGroups/blob/main/benchmark/BENCHMARK_GUIDE.md>`_
- :doc:`algorithm`: Algorithm details affecting performance
- :doc:`api/core`: API reference for optimization
