Fragmentation Module
====================

This module provides functions for fragmenting molecules based on bond
dissociation energies (BDE).

.. module:: PFASgroups.fragmentation
   :synopsis: Molecular fragmentation utilities

Overview
--------

The fragmentation module enables breaking molecules at their weakest bonds,
useful for:

- Predicting degradation products
- Analyzing metabolic pathways
- Understanding molecular stability
- Generating molecular fragments for analysis

Main Functions
--------------

generate_fragments
^^^^^^^^^^^^^^^^^^

.. py:function:: generate_fragments(mol, bondWeights=None, nb_breakingpoints=5)

   Break a molecule on its weakest bonds using bond dissociation energies.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param bondWeights: Dictionary of bond dissociation energies (auto-loaded if None)
   :type bondWeights: dict, optional
   :param nb_breakingpoints: Number of bonds to break
   :type nb_breakingpoints: int, optional
   :returns: List of fragment molecules
   :rtype: list[rdkit.Chem.Mol]

   **Algorithm:**

   1. Calculate bond strengths from BDE lookup table
   2. Compute probability distribution (weaker bonds more likely to break)
   3. Randomly select bonds to break based on probabilities
   4. Generate fragments from broken bonds

   **Example:**

   .. code-block:: python

      from PFASgroups.fragmentation import generate_fragments
      from rdkit import Chem

      mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
      
      # Generate fragments
      fragments = generate_fragments(mol, nb_breakingpoints=3)

      for frag in fragments:
          print(Chem.MolToSmiles(frag))

generate_degradation_products
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: generate_degradation_products(mol, **kwargs)

   Generate potential degradation products by breaking weak bonds.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param kwargs: Additional parameters passed to generate_fragments
   :returns: List of degradation product molecules
   :rtype: list[rdkit.Chem.Mol]

get_fragments
^^^^^^^^^^^^^

.. py:function:: get_fragments(mol, idx, all=False)

   Fragment a molecule by breaking specified bonds.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param idx: List of bond indices to break
   :type idx: list[int]
   :param all: If True, break all bonds at once; if False, break one at a time
   :type all: bool, optional
   :returns: List of fragment molecules
   :rtype: list[rdkit.Chem.Mol]

   **Example:**

   .. code-block:: python

      from PFASgroups.fragmentation import get_fragments
      from rdkit import Chem

      mol = Chem.MolFromSmiles("CCCC")
      
      # Break bond at index 1
      fragments = get_fragments(mol, [1], all=True)

Bond Dissociation Energies
--------------------------

The module uses a lookup table of bond dissociation energies (in kJ/mol):

.. list-table:: Default BDE Values
   :header-rows: 1
   :widths: 20 20 20

   * - Bond
     - BDE (kJ/mol)
     - Notes
   * - C-C
     - 346
     - Single bond
   * - C-H
     - 414
     - Carbon-hydrogen
   * - C-F
     - 485
     - Strong C-F bond
   * - C-Cl
     - 339
     - 
   * - C-Br
     - 285
     - Weaker halogen
   * - F-F
     - 159
     - Weak
   * - H-H
     - 436
     - 
   * - H-F
     - 568
     - Very strong

Lower BDE values indicate weaker bonds that are more likely to break.

yield_scheme Decorator
^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: yield_scheme(name)

   Decorator that provides BDE data to fragmentation functions.

   :param name: Scheme name (currently only 'BDE' supported)
   :type name: str

Helper Functions
----------------

keysToInt
^^^^^^^^^

.. py:function:: keysToInt(x)

   Convert string dictionary keys to integers.

   :param x: Dictionary with string keys
   :type x: dict
   :returns: Dictionary with integer keys
   :rtype: dict

fragment
^^^^^^^^

.. py:function:: fragment(mol, bde_dict, max_depth=4, nb_breakingpoints=5)

   Fragment a molecule into smaller fragments iteratively.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param bde_dict: Bond dissociation energy dictionary
   :type bde_dict: dict
   :param max_depth: Maximum fragmentation depth
   :type max_depth: int, optional
   :param nb_breakingpoints: Bonds to break per iteration
   :type nb_breakingpoints: int, optional
   :returns: List of fragment molecules
   :rtype: list[rdkit.Chem.Mol]

Use Cases
---------

Environmental Degradation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups.fragmentation import generate_degradation_products
   from PFASgroups import parse_smiles
   from rdkit import Chem

   # Analyze potential environmental degradation
   pfos = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O")
   
   products = generate_degradation_products(pfos)
   
   for product in products:
       smiles = Chem.MolToSmiles(product)
       # Classify the degradation products
       results = parse_smiles([smiles])
       print(f"{smiles}: {[g.name for g, _, _, _ in results[0]]}")

Metabolite Prediction
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups.fragmentation import generate_fragments
   from rdkit import Chem

   def predict_metabolites(smiles):
       """Predict potential metabolites from weak bond cleavage."""
       mol = Chem.MolFromSmiles(smiles)
       fragments = generate_fragments(mol, nb_breakingpoints=2)
       
       metabolites = []
       for frag in fragments:
           # Add hydrogens to satisfy valence
           frag = Chem.AddHs(frag)
           try:
               Chem.SanitizeMol(frag)
               metabolites.append(Chem.MolToSmiles(frag))
           except:
               pass
       
       return metabolites

Notes and Limitations
---------------------

1. **Probabilistic**: Fragment selection is random based on BDE probabilities
2. **No radical chemistry**: Does not model radical intermediates
3. **Simplified BDE**: Uses average BDE values, not context-specific
4. **Valence issues**: Some fragments may have valence issues requiring post-processing
5. **No stereochemistry**: Stereochemistry is not preserved in fragments
