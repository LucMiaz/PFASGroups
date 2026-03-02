Molecule Prioritization
=======================

The ``prioritise_molecules`` function enables systematic ranking of PFAS compounds based on either similarity to reference compounds or intrinsic fluorination characteristics.

Overview
--------

Prioritization is useful for:

- **Regulatory compliance**: Identify molecules similar to regulated compounds
- **Risk assessment**: Rank compounds by environmental persistence or bioaccumulation potential
- **Screening**: Focus analytical resources on highest-priority molecules
- **Chemical inventory management**: Categorize compounds by hazard potential

Quick Start
-----------

.. code-block:: python

   from HalogenGroups import prioritise_molecules, get_priority_statistics
   
   # Prioritize by similarity to known compounds
   reference = ["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"]  # PFOA
   results, scores = prioritise_molecules(
       inventory_smiles,
       reference=reference
   )
   
   # Prioritize by fluorination characteristics
   results, scores = prioritise_molecules(
       inventory_smiles,
       a=1.0,  # Weight for total fluorination
       b=2.0,  # Weight for longest chains
       percentile=90
   )
   
   # Get statistics
   stats = get_priority_statistics(results, scores, top_n=10)
   print(f"Top compound: {stats['top_n_smiles'][0]}")

Functions
---------

prioritise_molecules()
^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: prioritise_molecules(molecules, reference=None, group_selection='all', count_mode='binary', a=1.0, b=1.0, percentile=90.0, return_scores=True, ascending=False)

   Prioritize PFAS molecules based on similarity to reference or intrinsic properties.
   
   :param molecules: Molecules to prioritize (SMILES list, Mol list, or ResultsModel)
   :type molecules: list or ResultsModel
   :param reference: Reference molecules for similarity comparison (optional)
   :type reference: list or ResultsModel, optional
   :param group_selection: PFAS group selection ('all', 'oecd', 'generic', etc.)
   :type group_selection: str
   :param count_mode: Fingerprint encoding ('binary', 'count', 'max_component')
   :type count_mode: str
   :param a: Weight for total fluorinated component size
   :type a: float
   :param b: Weight for component size percentile
   :type b: float
   :param percentile: Percentile for component sizes (0-100)
   :type percentile: float
   :param return_scores: Whether to return scores with results
   :type return_scores: bool
   :param ascending: Sort order (False = highest priority first)
   :type ascending: bool
   :returns: Prioritized results and optionally scores
   :rtype: ResultsModel or tuple(ResultsModel, np.ndarray)

get_priority_statistics()
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: get_priority_statistics(results, scores, top_n=10)

   Get statistics about prioritization results.
   
   :param results: Prioritized molecules
   :type results: ResultsModel
   :param scores: Priority scores
   :type scores: np.ndarray
   :param top_n: Number of top molecules to analyze
   :type top_n: int
   :returns: Statistics dictionary
   :rtype: dict

Prioritization Strategies
--------------------------

Reference-Based Prioritization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a reference list is provided, molecules are ranked by distributional similarity using KL divergence:

.. math::

   \text{Priority} = 1 - \text{KL}_{\text{minmax}}(P_{\text{mol}} \| P_{\text{ref}})

Where:
- Lower KL divergence = more similar to reference = higher priority
- Accounts for overall composition of PFAS groups
- Useful for finding analogs of regulated compounds

**Use cases:**

- Identify structural analogs of PFOA, PFOS, or other regulated PFASs
- Screen inventories for compounds similar to substances of concern
- Compare site-specific contamination to reference databases

**Example:**

.. code-block:: python

   # Find compounds similar to long-chain PFCAs
   reference = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHpA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFNA
   ]
   
   results, scores = prioritise_molecules(
       inventory,
       reference=reference,
       group_selection='all',
       count_mode='binary'
   )

Intrinsic Prioritization
^^^^^^^^^^^^^^^^^^^^^^^^^

Without a reference, molecules are ranked by fluorination characteristics:

.. math::

   \text{Priority} = a \cdot \sum(\text{component sizes}) + b \cdot P_p(\text{component sizes})

Where:
- :math:`\sum(\text{component sizes})`: Total fluorinated carbons
- :math:`P_p(\text{component sizes})`: pth percentile of component sizes
- :math:`a, b`: Tunable weights
- :math:`p`: Percentile level (0-100)

**Parameter Guidelines:**

+------------------------+--------+--------+-------------+----------------------------------+
| **Priority Criterion** | **a**  | **b**  | **p**       | **Rationale**                    |
+========================+========+========+=============+==================================+
| Long-chain compounds   | 0.5    | 2.0    | 90-95       | Emphasize longest chains         |
+------------------------+--------+--------+-------------+----------------------------------+
| Total fluorine burden  | 2.0    | 0.5    | 50-75       | Emphasize total fluorination     |
+------------------------+--------+--------+-------------+----------------------------------+
| Balanced approach      | 1.0    | 1.0    | 75-90       | Balance total and chain length   |
+------------------------+--------+--------+-------------+----------------------------------+
| Bioaccumulation        | 1.0    | 1.5    | 80-90       | Moderate chains, high total F    |
+------------------------+--------+--------+-------------+----------------------------------+

**Example:**

.. code-block:: python

   # Prioritize for environmental persistence (long chains)
   results, scores = prioritise_molecules(
       inventory,
       a=0.5,
       b=2.0,
       percentile=95
   )
   
   # Prioritize for total fluorination
   results, scores = prioritise_molecules(
       inventory,
       a=2.0,
       b=0.5,
       percentile=50
   )

Real-World Applications
-----------------------

Regulatory Screening
^^^^^^^^^^^^^^^^^^^^

Identify compounds that may fall under regulatory scrutiny:

.. code-block:: python

   # Screen for PFOA/PFOS-like compounds
   regulated = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
   ]
   
   results, scores = prioritise_molecules(
       site_inventory,
       reference=regulated,
       group_selection='oecd'
   )
   
   # Report high-priority compounds
   stats = get_priority_statistics(results, scores, top_n=20)
   print(f"Found {len([s for s in scores if s > 0.8])} compounds with >80% similarity")

Environmental Risk Assessment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rank compounds by persistence potential:

.. code-block:: python

   # Prioritize long-chain compounds (C8+)
   results, scores = prioritise_molecules(
       inventory,
       a=0.3,
       b=2.5,
       percentile=95
   )
   
   # Identify highest-risk compounds
   high_risk = []
   for i, (result, score) in enumerate(zip(results[:10], scores[:10])):
       # Check for long perfluorinated chains
       for match in result['matches']:
           if 'perfluoro' in match.get('group_name', '').lower():
               if max(match.get('components_sizes', [0])) >= 8:
                   high_risk.append((result['smiles'], score))
   
   print(f"Identified {len(high_risk)} high-risk compounds")

Analytical Method Development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Select representative compounds for method validation:

.. code-block:: python

   # Stratified sampling across priority ranges
   results, scores = prioritise_molecules(inventory, a=1.0, b=1.0, percentile=75)
   
   # Select compounds from each priority tier
   n_total = len(results)
   high_priority = results[:n_total//3]  # Top 33%
   medium_priority = results[n_total//3:2*n_total//3]  # Middle 33%
   low_priority = results[2*n_total//3:]  # Bottom 33%
   
   # Sample from each tier
   import random
   validation_set = (
       random.sample(list(high_priority), 5) +
       random.sample(list(medium_priority), 3) +
       random.sample(list(low_priority), 2)
   )

Output Interpretation
---------------------

Priority Scores
^^^^^^^^^^^^^^^

**Reference-based scores** (0-1 range):
- **0.8-1.0**: Very similar to reference, highest priority
- **0.6-0.8**: Moderately similar, medium-high priority
- **0.4-0.6**: Some similarity, medium priority
- **0.2-0.4**: Low similarity, low-medium priority
- **0.0-0.2**: Very different from reference, lowest priority

**Intrinsic scores** (unbounded, typically 0-100+):
- Relative ranking is more important than absolute values
- Compare within dataset, not across different analyses
- Higher scores indicate more/longer fluorinated chains

Statistics Dictionary
^^^^^^^^^^^^^^^^^^^^^^

The ``get_priority_statistics`` function returns:

.. code-block:: python

   {
       'n_molecules': 100,
       'score_mean': 45.2,
       'score_std': 12.3,
       'score_min': 8.5,
       'score_max': 78.9,
       'score_median': 42.1,
       'top_n_smiles': ['SMILES1', 'SMILES2', ...],
       'top_n_scores': [78.9, 67.4, ...],
       'top_n_groups': [('perfluoroalkyl', 10), ('carboxylic acid', 8), ...]
   }

Best Practices
--------------

1. **Choose appropriate reference compounds**
   - Use structurally relevant references
   - Include multiple reference compounds for robust comparison
   - Consider different functional group classes

2. **Tune parameters for your use case**
   - Test different a/b/percentile combinations
   - Validate against known priority compounds
   - Document your parameter choices

3. **Consider multiple prioritization strategies**
   - Compare reference-based and intrinsic approaches
   - Use domain knowledge to interpret scores
   - Combine with other data (toxicity, exposure, etc.)

4. **Handle edge cases**
   - Non-PFAS compounds will score low (as expected)
   - Single or few molecules may give unstable results
   - Verify results for unusual molecules

5. **Document and report**
   - Save priority scores with results
   - Report methodology and parameters
   - Include statistics in summaries

Limitations
-----------

- **Structural diversity**: Prioritization is based on fluorination patterns, not complete molecular structure
- **Biological activity**: Does not predict toxicity or biological effects
- **Data quality**: Results depend on accurate PFAS group identification
- **Parameter sensitivity**: Intrinsic scores can be sensitive to a/b/percentile choices
- **Reference dependency**: Reference-based scores depend heavily on reference selection

See Also
--------

- :doc:`fingerprint_analysis` - Fingerprint generation and analysis
- :doc:`api/core` - Core parsing functions
- :doc:`tutorial` - General PFASgroups tutorial
- ``examples/prioritization_examples.py`` - Working examples

References
----------

1. Kullback, S., & Leibler, R. A. (1951). On Information and Sufficiency. *The Annals of Mathematical Statistics*, 22(1):79-86.

2. Wang, Z., et al. (2017). A Never-Ending Story of Per- and Polyfluoroalkyl Substances (PFASs)? *Environmental Science & Technology*, 51(5):2508-2518.

3. Glüge, J., et al. (2020). An overview of the uses of per- and polyfluoroalkyl substances (PFAS). *Environmental Science: Processes & Impacts*, 22(12):2345-2373.
