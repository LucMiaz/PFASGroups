Benchmarking
============

.. code-block:: python

   from PFASGroups import parse_smiles, generate_fingerprint

   # Parse 10 k molecules
   smiles = [...]   # your list
   results = parse_smiles(smiles)
   fps, info = generate_fingerprint(smiles)   # shape (n, 116)

This page describes the validation studies performed to assess the accuracy
and completeness of the PFASGroups library.

.. contents:: Contents
   :local:
   :depth: 2

Validation against PFASSTRUCTv5
--------------------------------

**Dataset**: PFASSTRUCTv5 is a curated annotated PFAS structure database
containing ~14 000 structures
(`Schymanski et al. 2023 <https://zenodo.org/record/7370805>`_).

**Method**: Each structure was parsed with ``parse_smiles()`` using the
default fluorine-only mode and OECD 2021 group definitions.  The binary
PFAS/non-PFAS label from PFASSTRUCTv5 was used as the ground truth.

**Results**:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Metric
     - Value
     - Notes
   * - Sensitivity (recall)
     - > 0.97
     - Fraction of known PFAS correctly detected
   * - Specificity
     - > 0.90
     - Fraction of non-PFAS correctly excluded
   * - F1 score
     - > 0.96
     -

Failures are mostly due to:

* Highly functionalised PFAS with unusual connectivity not covered by
  any current SMARTS pattern
* Partially fluorinated structures at the boundary of the definition

Comparison with CSRML classifier
----------------------------------

The CSRML (Chemical Structure Rule Markup Language) classifier from the OECD
toolbox was used as an external reference.

**Key findings**:

* PFASGroups matches or exceeds CSRML accuracy on structures with at least
  one perfluoroalkyl chain.
* PFASGroups additionally detects fluorotelomer groups not covered by the
  CSRML rule set.
* For borderline polyfluoroalkyl structures sensitivity is comparable (both
  ~0.85).

Running the benchmark scripts
-------------------------------

Benchmark scripts are available in the ``benchmark/`` directory of the source
repository:

.. code-block:: bash

   cd benchmark
   # Reproduce PFASSTRUCTv5 validation
   python benchmark_pfasstructv5.py

   # Compare with CSRML classifier
   python benchmark_csrml.py

   # Accuracy report across all group categories
   python accuracy_report.py

Expected outputs are in ``benchmark/results/``.

Performance
-----------

Typical throughput on a modern laptop (single core):

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Task
     - Molecules
     - Time
   * - ``parse_smiles`` with defaults
     - 10 000
     - ~5 s
   * - ``generate_fingerprint``
     - 10 000
     - ~8 s
   * - ``parse_smiles`` with ``compute_component_metrics=True``
     - 10 000
     - ~12 s

For large datasets (>100 k molecules) consider using the
``compute_component_metrics=False`` flag to skip the effective graph
resistance computation:

.. code-block:: python

   results = parse_smiles(large_smiles_list, compute_component_metrics=False)