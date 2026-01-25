Homologue Generation
====================

This module provides functions for generating homologue series by iteratively
removing CF₂ moieties from fluorinated chains.

.. module:: PFASgroups.generate_homologues
   :synopsis: Homologue series generation

Overview
--------

Homologue generation creates shorter-chain analogues of PFAS compounds by
systematically removing repeating units (typically CF₂) from the fluorinated
chain. This is useful for:

- Exploring theoretical chemical space
- Generating training data for machine learning
- Understanding structure-activity relationships
- Predicting degradation products

Main Functions
--------------

generate_homologues
^^^^^^^^^^^^^^^^^^^

.. py:function:: generate_homologues(mol, smartsPathName='Perfluoroalkyl', smartsPaths=None, repeating='C(F)F', base_repeating=['C'])

   Generate homologue series by removing repeating units from fluorinated chains.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param smartsPathName: Pathway type ('Perfluoroalkyl' or 'Polyfluoroalkyl')
   :type smartsPathName: str, optional
   :param smartsPaths: Custom pathway definitions (optional)
   :type smartsPaths: dict, optional
   :param repeating: SMARTS pattern for the repeating unit
   :type repeating: str, optional
   :param base_repeating: Base atoms that should not be removed
   :type base_repeating: list[str], optional
   :returns: Dictionary mapping InChIKey → {formula: Mol}
   :rtype: dict

   **Algorithm:**

   1. Find all fluorinated chains in the molecule
   2. Identify consecutive repeating units (CF₂) in each chain
   3. Generate all combinations of unit removals
   4. Reconnect the molecular structure after removal
   5. Deduplicate results by InChIKey

   **Example - Perfluoroalkyl:**

   .. code-block:: python

      from PFASgroups import generate_homologues
      from rdkit import Chem

      # Start with PFOA (C8)
      mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
      
      # Generate homologues
      homologues = generate_homologues(mol, smartsPathName='Perfluoroalkyl')

      # Print results
      for inchikey, formulas in homologues.items():
          for formula, h_mol in formulas.items():
              smiles = Chem.MolToSmiles(h_mol)
              print(f"{formula}: {smiles}")

   Output::

      C7HF13O2: FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O
      C6HF11O2: FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O
      C5HF9O2: FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O
      C4HF7O2: FC(F)(F)C(F)(F)C(F)(F)C(=O)O
      ...

   **Example - Polyfluoroalkyl:**

   .. code-block:: python

      # For polyfluoroalkyl compounds
      homologues = generate_homologues(
          mol, 
          smartsPathName='Polyfluoroalkyl',
          repeating='C([F,H,I,Br,Cl])[F,H,I,Br,Cl]'
      )

find_chain
^^^^^^^^^^

.. py:function:: find_chain(mol, pathsmarts, endsmarts, repeating='C(F)(F)')

   Find fluorinated chains between terminal groups.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param pathsmarts: SMARTS for chain atoms
   :type pathsmarts: rdkit.Chem.Mol
   :param endsmarts: SMARTS for terminal groups
   :type endsmarts: rdkit.Chem.Mol
   :param repeating: SMARTS for repeating unit
   :type repeating: str, optional
   :returns: List of chain dictionaries with chain info
   :rtype: list[dict]

   Each chain dictionary contains:

   - ``chain``: Set of atom indices in the chain
   - ``start``: Start atom index
   - ``end``: End atom index
   - ``belly``: Set of repeating unit atom indices
   - ``consecutive_parts``: List of consecutive repeating unit segments

Use Cases
---------

Degradation Pathway Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups import generate_homologues
   from rdkit import Chem

   # Analyze potential degradation products
   pfoa = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
   
   homologues = generate_homologues(pfoa)
   
   print("Potential degradation products:")
   for inchikey, formulas in homologues.items():
       for formula, mol in formulas.items():
           n_carbons = formula.count('C')
           print(f"  C{n_carbons}: {Chem.MolToSmiles(mol)}")

Training Data Augmentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups import generate_homologues
   from rdkit import Chem

   def augment_with_homologues(smiles_list):
       """Expand dataset with homologues."""
       all_molecules = set()
       
       for smiles in smiles_list:
           mol = Chem.MolFromSmiles(smiles)
           if mol is None:
               continue
           
           all_molecules.add(smiles)
           
           homologues = generate_homologues(mol)
           for formulas in homologues.values():
               for h_mol in formulas.values():
                   all_molecules.add(Chem.MolToSmiles(h_mol))
       
       return list(all_molecules)

   # Original: 100 molecules → Augmented: potentially 500+
   augmented = augment_with_homologues(original_smiles)

Notes and Limitations
---------------------

1. **Linear chains only**: Works best with linear fluorinated chains
2. **Functional groups preserved**: Terminal functional groups are maintained
3. **Valence checking**: Resulting molecules are validated for correct valence
4. **Fragmentation**: Will raise error if removal causes fragmentation
5. **Minimum chain length**: At least one carbon must remain in the chain
