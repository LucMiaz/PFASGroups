Quick Start Guide
=================

This guide will get you up and running with PFASgroups in minutes.

Basic PFAS Detection
---------------------

Detecting PFAS groups in molecules:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Single SMILES string
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   results = parse_smiles(smiles)
   
   # Results format: [(PFASGroup, match_count, chain_lengths, matched_chains)]
   for group, count, lengths, chains in results[0]:
       print(f"{group.name}: {count} matches")
       print(f"  Chain lengths: {lengths}")
       print(f"  OECD Group ID: {group.id}")

**Output:**

.. code-block:: text

   Perfluoroalkyl carboxylic acids: 1 matches
     Chain lengths: [5]
     OECD Group ID: 1

Multiple Molecules
------------------

Process a list of SMILES:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",  # PFBA
       "FC(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
       "C(C(F)(F)F)F"  # Not a typical PFAS functional group
   ]
   
   results = parse_smiles(smiles_list)
   
   for i, smiles in enumerate(smiles_list):
       print(f"\n{smiles}:")
       if results[i]:
           for group, count, lengths, _ in results[i]:
               print(f"  ✓ {group.name} (chains: {lengths})")
       else:
           print("  ✗ No PFAS groups detected")

**Output:**

.. code-block:: text

   FC(F)(F)C(F)(F)C(=O)O:
     ✓ Perfluoroalkyl carboxylic acids (chains: [3])
   
   FC(F)(F)C(F)(F)S(=O)(=O)O:
     ✓ Perfluoroalkyl sulfonic acids (chains: [3])
   
   C(C(F)(F)F)F:
     ✗ No PFAS groups detected

Generating Fingerprints
------------------------

Create binary fingerprints for machine learning:

.. code-block:: python

   from PFASgroups import generate_fingerprint
   import numpy as np
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFBA
       "FC(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
       "C(C(F)(F)C(F)(F)C(F)(F)F)O"  # FTOH
   ]
   
   # Generate binary fingerprints
   fps, group_info = generate_fingerprint(smiles_list)
   
   print(f"Fingerprint shape: {fps.shape}")
   print(f"Number of groups: {len(group_info)}")
   
   # Use in ML model
   from sklearn.ensemble import RandomForestClassifier
   
   X = fps  # Features
   y = [1, 1, 0]  # Labels (1=PFAS, 0=not PFAS)
   
   model = RandomForestClassifier()
   model.fit(X, y)

**Different fingerprint representations:**

.. code-block:: python

   # Dictionary format (sparse)
   fps_dict, _ = generate_fingerprint(smiles_list, representation='dict')
   print(fps_dict[0])  # {group_id: 1, ...}
   
   # Count-based (instead of binary)
   fps_count, _ = generate_fingerprint(smiles_list, count_mode='count')
   print(fps_count[0])  # [0, 1, 0, 0, ...]  (counts, not binary)
   
   # Max chain length
   fps_chain, _ = generate_fingerprint(smiles_list, count_mode='max_chain')
   print(fps_chain[0])  # [0, 5, 0, 0, ...]  (max chain length per group)

Visualization
-------------

Visualize PFAS group assignments:

.. code-block:: python

   from PFASgroups import plot_pfasgroups
   
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   
   # Display in Jupyter notebook
   plot_pfasgroups(smiles, display=True)
   
   # Save to file
   plot_pfasgroups(smiles, path="pfba_annotated.png")
   
   # Generate SVG
   plot_pfasgroups(smiles, svg=True, path="pfba_annotated.svg")
   
   # Multiple molecules in grid
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "C(C(F)(F)C(F)(F)F)O"
   ]
   plot_pfasgroups(smiles_list, ncols=3, path="pfas_comparison.png")

Component-Based Analysis
------------------------

For molecules with disconnected or cyclic components:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Complex molecule with cycles
   smiles = "C1(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)(F)"
   
   # Standard analysis
   results_standard = parse_smiles(smiles, bycomponent=False)
   
   # Component-based analysis (better for cycles)
   results_component = parse_smiles(smiles, bycomponent=True)
   
   print("Standard:", len(results_standard[0]))
   print("Component-based:", len(results_component[0]))

Command Line Interface
----------------------

PFASgroups provides a CLI for quick analyses:

**Parse SMILES:**

.. code-block:: bash

   pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O"
   
   # Multiple SMILES
   pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O" "FC(F)(F)C(F)(F)S(=O)(=O)O"
   
   # From file
   pfasgroups parse --file molecules.smi

**Generate fingerprints:**

.. code-block:: bash

   pfasgroups fingerprint "FC(F)(F)C(F)(F)C(=O)O" --format dict
   
   # Save to file
   pfasgroups fingerprint --file molecules.smi --output fingerprints.json

**List available groups:**

.. code-block:: bash

   pfasgroups list-groups
   
   # Filter by main group
   pfasgroups list-groups --main "Perfluoroalkyl acids"

**List pathway types:**

.. code-block:: bash

   pfasgroups list-paths

Working with DataFrames
------------------------

Output results as pandas DataFrame:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "C(C(F)(F)C(F)(F)F)O"
   ]
   
   # Get DataFrame output
   df = parse_smiles(smiles_list, output_format='dataframe')
   
   print(df.head())
   print(f"\nShape: {df.shape}")
   
   # Filter for specific groups
   pfcas = df[df['group_name'] == 'Perfluoroalkyl carboxylic acids']
   print(f"\nPFCAs found: {len(pfcas)}")
   
   # Export to CSV
   df.to_csv('pfas_analysis.csv', index=False)

Reading from Files
------------------

Process SMILES from various file formats:

**From SMILES file (.smi):**

.. code-block:: python

   from PFASgroups import parse_smiles
   import pandas as pd
   
   # Read SMILES file
   with open('molecules.smi', 'r') as f:
       smiles_list = [line.strip() for line in f if line.strip()]
   
   results = parse_smiles(smiles_list)

**From CSV:**

.. code-block:: python

   import pandas as pd
   from PFASgroups import parse_smiles
   
   # Read CSV with SMILES column
   df = pd.read_csv('compounds.csv')
   smiles_list = df['SMILES'].tolist()
   
   results = parse_smiles(smiles_list)
   
   # Add results to DataFrame
   df['pfas_groups'] = [[g.name for g, _, _, _ in r] for r in results]

**From SDF:**

.. code-block:: python

   from rdkit import Chem
   from PFASgroups import parse_groups_in_mol
   
   # Read SDF file
   suppl = Chem.SDMolSupplier('compounds.sdf')
   
   results = []
   for mol in suppl:
       if mol is not None:
           result = parse_groups_in_mol(mol)
           results.append(result)

Error Handling
--------------

Robust error handling:

.. code-block:: python

   from PFASgroups import parse_smiles
   from rdkit import Chem
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",  # Valid
       "INVALID_SMILES",  # Invalid
       "CCO",  # Valid but not PFAS
   ]
   
   for smiles in smiles_list:
       try:
           # Check if valid SMILES first
           mol = Chem.MolFromSmiles(smiles)
           if mol is None:
               print(f"✗ Invalid SMILES: {smiles}")
               continue
           
           # Parse PFAS groups
           results = parse_smiles(smiles)
           
           if results and results[0]:
               print(f"✓ {smiles}: {len(results[0])} groups")
           else:
               print(f"○ {smiles}: No PFAS groups")
               
       except Exception as e:
           print(f"✗ Error processing {smiles}: {e}")

Next Steps
----------

- Read the :doc:`tutorial` for in-depth examples
- Learn about the :doc:`algorithm` behind the detection
- Explore :doc:`customization` for advanced use cases
- Check the :doc:`api/core` for complete API reference
- See :doc:`benchmarking` for validation and performance analysis
