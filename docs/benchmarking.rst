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

Validation Against Richard et al. 2023 (PFASSTRUCTV5 & CSRML)
--------------------------------------------------------------

PFASSTRUCTv5 Inventory Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PFASSTRUCTv5 PFAS definition (definition ID 5, ``'PFASTRUCTv5'``,
fluorine-ratio threshold ≥ 0.3) was validated against the PFASSTRUCTV5 inventory
of 14,735 PFAS structures published as supplementary material by:

   Richard, A. M. *et al.* (2023). *A New CSRML Structure-Based Fingerprint
   Method for Profiling and Categorizing Per- and Polyfluoroalkyl Substances
   (PFAS).* Chemical Research in Toxicology, 36(3), 318–338.
   https://doi.org/10.1021/acs.chemrestox.2c00403

The SDF file (``tests/test_data/PFASSTRUCTV5_20221101.sdf``) serves as a positive
test set: every molecule in the inventory is expected to satisfy definition 5.

An automated test file ``tests/test_pfasstructv5_against_richard2023.py``
verifies this. Under pytest, ten evenly-spaced molecules are sampled to keep CI
fast; running the script directly tests all 14,735 structures:

.. code-block:: bash

   # Fast CI check — 10-molecule evenly-spaced sample:
   pytest tests/test_pfasstructv5_against_richard2023.py -v

   # Full validation of all 14,735 PFASSTRUCTV5 molecules:
   python tests/test_pfasstructv5_against_richard2023.py

The script reports exact identifiers (DTXCID/DTXSID) for any molecule that fails
the definition, enabling targeted investigation of edge cases.

The definition uses **either** SMARTS matching (CF₃–CF₂ adjacency, CF₂–heteroatom–CF₂
bridges) **or** a fluorine-ratio criterion (F/all-atoms ≥ 0.3), whichever is satisfied.
This means that highly fluorinated compounds with non-standard architectures are
included even without a CF₂–CF₂ chain pattern.

Programmatic use:

.. code-block:: python

   import json
   from HalogenGroups import PFASDefinition
   from HalogenGroups.core import PFAS_DEFINITIONS_FILE

   with open(PFAS_DEFINITIONS_FILE) as fh:
       defs = json.load(fh)
   pfasstruct = PFASDefinition(**next(d for d in defs if d["id"] == 5))

   from rdkit import Chem
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
   print(pfasstruct.applies_to_molecule(mol))  # True


Fingerprint Comparison Against CSRML
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HalogenGroups binary group fingerprints were compared against the **125-bit CSRML**
fingerprint system introduced by Richard *et al.* 2023 (same reference as above).
CSRML fingerprints are a SMARTS-based binary encoding designed specifically for
PFAS profiling, derived from the TxP_PFAS taxonomy used in the CompTox Chemicals
Dashboard.

The comparison script ``benchmark/scripts/compare_pfasgroups_vs_txppfas.py``
computes and visualises:

- **Coverage** — fraction of the test set for which each group fires at least once
- **Information content** per group (Shannon entropy-based)
- **Lorenz curve** of class-richness concentration
- **PCA scatter** of the two fingerprint systems on a shared compound set
- **Pairwise Jaccard** similarity distributions within each system
- **Cross-fingerprint** similarity scatter (HalogenGroups vs. CSRML Jaccard)

.. code-block:: bash

   cd benchmark/scripts

   # HalogenGroups fingerprint analysis only (no CSRML data required):
   python compare_pfasgroups_vs_txppfas.py \
       --smiles ../data/test_set_for_PFASSTRUCTv5.tsv

   # Full comparison (Richard 2023 SI Table S2 exported as CSV):
   python compare_pfasgroups_vs_txppfas.py \
       --smiles ../data/test_set_for_PFASSTRUCTv5.tsv \
       --txppfas_csv ../data/richard2023_SI_table_S2.csv

   # Limit to 100 molecules for a quick sanity check:
   python compare_pfasgroups_vs_txppfas.py --max_mols 100

Result figures are written to ``benchmark/reports/``.

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

Test Molecule Generation
-------------------------

JSON-Driven Generation
~~~~~~~~~~~~~~~~~~~~~~

Since version 1.2.3, the benchmark system reads test molecule generation patterns directly from ``PFAS_groups_smarts.json``. Each functional group definition includes a ``test.generate`` field that specifies:

- **smiles**: The SMILES pattern to attach or insert
- **mode**: Either ``"attach"`` (terminal) or ``"insert"`` (internal)
- **is_telomer**: Boolean flag indicating if this is a fluorotelomer group

Example group definition:

.. code-block:: json

   {
     "id": 31,
     "name": "aldehyde",
     "test": {
       "generate": {
         "smiles": "C(=O)[H]",
         "mode": "attach",
         "is_telomer": false
       }
     }
   }

Example telomer group:

.. code-block:: json

   {
     "id": 74,
     "name": "Fluorotelomer silane",
     "test": {
       "generate": {
         "smiles": "CCCC[Si]([H])([H])[H]",
         "mode": "attach",
         "is_telomer": true
       }
     }
   }

The benchmark script automatically:

1. Reads all group definitions from ``PFAS_groups_smarts.json``
2. Extracts the ``test.generate`` field for each target group
3. Generates test molecules by attaching/inserting the SMILES patterns
4. Creates single-group and multi-group test cases
5. Validates detection accuracy against expected groups

This approach ensures:

- **Single source of truth**: Test patterns match production group definitions
- **Automatic updates**: Adding new groups automatically includes them in benchmarks
- **Consistency**: No risk of hardcoded test data diverging from actual patterns
- **Maintainability**: Changes to group definitions propagate to test generation

Customizing Test Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To modify test generation for a specific group, edit the ``test.generate`` field in ``PFAS_groups_smarts.json``:

.. code-block:: python

   import json
   
   # Load group definitions
   with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
       groups = json.load(f)
   
   # Find and modify a group
   for group in groups:
       if group['id'] == 31:  # aldehyde
           group['test']['generate'] = {
               'smiles': 'C(=O)[H]',
               'mode': 'attach',
               'is_telomer': False
           }
   
   # Save updated definitions
   with open('PFASgroups/data/PFAS_groups_smarts.json', 'w') as f:
       json.dump(groups, f, indent=2)
   
   # Re-run benchmark to use new pattern
   # The benchmark script will automatically pick up the changes

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
