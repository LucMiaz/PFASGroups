Tutorial
========

This comprehensive tutorial walks you through common use cases and advanced features of PFASgroups.

Part 1: Basic PFAS Detection
-----------------------------

Understanding PFAS Groups
~~~~~~~~~~~~~~~~~~~~~~~~~

PFAS groups represent functional groups attached to fluorinated chains. PFASgroups identifies 55+ different groups based on OECD terminology and extends them with generic classifications.

.. code-block:: python

   from PFASgroups import parse_smiles, get_PFASGroups
   
   # List all available groups
   groups = get_PFASGroups()
   
   # OECD groups (1-28)
   oecd_groups = groups[:28]
   for group in oecd_groups:
       print(f"{group.id}: {group.name} ({group.alias})")
   
   # Generic groups (29-55)
   generic_groups = groups[28:]
   for group in generic_groups:
       print(f"{group.id}: {group.name}")

Analyzing Common PFAS
~~~~~~~~~~~~~~~~~~~~~

Let's analyze several well-known PFAS compounds:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Common PFAS with their names
   pfas_compounds = {
       "PFOA": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
       "PFOS": "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",
       "PFBA": "FC(F)(F)C(F)(F)C(F)(F)C(=O)O",
       "6:2 FTOH": "C(C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)CO"
   }
   
   for name, smiles in pfas_compounds.items():
       print(f"\n{name}:")
       results = parse_smiles(smiles)
       
       if results[0]:
           for group, count, lengths, _ in results[0]:
               print(f"  ✓ {group.name}")
               print(f"    Chain length: {lengths}")
               print(f"    OECD Group: {group.id}")
       else:
           print("  ✗ No PFAS groups detected")

**Expected Output:**

.. code-block:: text

   PFOA:
     ✓ Perfluoroalkyl carboxylic acids
       Chain length: [8]
       OECD Group: 1
   
   PFOS:
     ✓ Perfluoroalkyl sulfonic acids
       Chain length: [8]
       OECD Group: 6
   
   PFBA:
     ✓ Perfluoroalkyl carboxylic acids
       Chain length: [4]
       OECD Group: 1
   
   6:2 FTOH:
     ✓ Fluorotelomer alcohols
       Chain length: [6]
       OECD Group: 15

Part 2: Working with Data
--------------------------

Batch Processing
~~~~~~~~~~~~~~~~

Process multiple molecules efficiently:

.. code-block:: python

   import pandas as pd
   from PFASgroups import parse_smiles
   
   # Load data from CSV
   df = pd.read_csv('chemicals.csv')
   smiles_list = df['SMILES'].tolist()
   
   # Parse all molecules
   results = parse_smiles(smiles_list)
   
   # Extract information
   pfas_detected = []
   for i, result in enumerate(results):
       if result:
           groups = [g.name for g, _, _, _ in result]
           chains = [lengths for _, _, lengths, _ in result]
           pfas_detected.append({
               'smiles': smiles_list[i],
               'num_groups': len(result),
               'groups': ', '.join(groups),
               'chains': str(chains)
           })
   
   # Create results DataFrame
   results_df = pd.DataFrame(pfas_detected)
   results_df.to_csv('pfas_results.csv', index=False)

DataFrame Output
~~~~~~~~~~~~~~~~

Use the built-in DataFrame output for easier analysis:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "CCO"  # Non-PFAS
   ]
   
   # Get DataFrame directly
   df = parse_smiles(smiles_list, output_format='dataframe')
   
   # Filter and analyze
   print("\nStatistics:")
   print(f"Total molecules: {len(smiles_list)}")
   print(f"PFAS detected: {df['smiles'].nunique()}")
   print(f"\nGroup distribution:")
   print(df['group_name'].value_counts())
   
   # Export to Excel
   df.to_excel('pfas_analysis.xlsx', index=False)

Part 3: Fingerprints for ML
----------------------------

Binary Fingerprints
~~~~~~~~~~~~~~~~~~~

Create binary fingerprints for classification:

.. code-block:: python

   from PFASgroups import generate_fingerprint
   from sklearn.ensemble import RandomForestClassifier
   from sklearn.model_selection import train_test_split
   import numpy as np
   
   # Training data
   pfas_smiles = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "C(C(F)(F)C(F)(F)F)CO"
   ]
   non_pfas_smiles = [
       "CCO",
       "CC(=O)O",
       "CCS(=O)(=O)O"
   ]
   
   # Generate fingerprints
   all_smiles = pfas_smiles + non_pfas_smiles
   fps, info = generate_fingerprint(all_smiles)
   
   # Labels (1=PFAS, 0=non-PFAS)
   y = np.array([1]*len(pfas_smiles) + [0]*len(non_pfas_smiles))
   
   # Train classifier
   X_train, X_test, y_train, y_test = train_test_split(fps, y, test_size=0.3)
   
   clf = RandomForestClassifier(n_estimators=100)
   clf.fit(X_train, y_train)
   
   # Predict
   accuracy = clf.score(X_test, y_test)
   print(f"Accuracy: {accuracy:.2%}")
   
   # Feature importance
   important_groups = np.argsort(clf.feature_importances_)[-5:]
   print("\nMost important groups:")
   for idx in important_groups:
       print(f"  Group {idx+1}: {info[idx+1]}")

Count-Based Fingerprints
~~~~~~~~~~~~~~~~~~~~~~~~~

Use count-based fingerprints for regression:

.. code-block:: python

   from PFASgroups import generate_fingerprint
   
   # Molecules with varying PFAS content
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",  # 1 PFCA group
       "FC(F)(F)C(F)(F)C(=O)OC(=O)C(F)(F)C(F)(F)F",  # 2 PFCA groups
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"  # 1 PFCA, longer chain
   ]
   
   # Count mode - counts occurrences
   fps_count, info = generate_fingerprint(smiles_list, count_mode='count')
   print("Count fingerprints:")
   print(fps_count)
   
   # Max chain mode - records maximum chain length
   fps_chain, info = generate_fingerprint(smiles_list, count_mode='max_chain')
   print("\nMax chain fingerprints:")
   print(fps_chain)

Part 4: Visualization
---------------------

Basic Visualization
~~~~~~~~~~~~~~~~~~~

Visualize detected groups:

.. code-block:: python

   from PFASgroups import plot_pfasgroups
   
   # Single molecule
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(=O)O"
   plot_pfasgroups(smiles, path="pfba_annotated.png")

Grid Layouts
~~~~~~~~~~~~

Compare multiple molecules:

.. code-block:: python

   from PFASgroups import plot_pfasgroups
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",  # PFBA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHxA
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"  # PFOA
   ]
   
   plot_pfasgroups(
       smiles_list,
       ncols=3,
       subwidth=400,
       subheight=400,
       path="pfca_series.png"
   )

SVG for Publications
~~~~~~~~~~~~~~~~~~~~

Generate vector graphics:

.. code-block:: python

   from PFASgroups import plot_pfasgroups
   
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   
   # SVG (scalable)
   plot_pfasgroups(smiles, svg=True, path="pfhxa.svg")
   
   # Can be imported into Inkscape, Illustrator, or LaTeX

Part 5: Advanced Features
--------------------------

Component-Based Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

For cyclic or complex structures:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Perfluorocyclic compound
   cyclic_smiles = "FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F"
   
   # Standard analysis
   results_std = parse_smiles(cyclic_smiles, bycomponent=False)
   print(f"Standard: {len(results_std[0])} groups")
   
   # Component-based (better for cycles)
   results_comp = parse_smiles(cyclic_smiles, bycomponent=True)
   print(f"Component: {len(results_comp[0])} groups")

Homologue Generation
~~~~~~~~~~~~~~~~~~~~

Generate shorter-chain analogues:

.. code-block:: python

   from PFASgroups import generate_homologues
   from rdkit import Chem
   
   # PFOA (C8)
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
   
   # Generate homologues (C7, C6, C5, C4)
   homologues = generate_homologues(mol, smartsPathName='Perfluoroalkyl')
   
   print(f"Generated {len(homologues)} homologues:")
   for inchikey, formulas in homologues.items():
       for formula, homologue_mol in formulas.items():
           smiles = Chem.MolToSmiles(homologue_mol)
           print(f"  {formula}: {smiles}")

Custom PFAS Groups
~~~~~~~~~~~~~~~~~~

Define custom groups:

.. code-block:: python

   from PFASgroups import PFASGroup, parse_smiles, get_PFASGroups
   from rdkit import Chem
   
   # Get default groups
   groups = get_PFASGroups()
   
   # Add custom group
   custom_group = PFASGroup(
       id=100,
       name="Perfluoroalkyl nitrates",
       alias="PFANs",
       smarts1=Chem.MolFromSmarts("[C]-[O]-[N+](=O)[O-]"),
       smartsPath="Perfluoroalkyl",
       constraints={
           "eq": {"N": 1, "O": 3},
           "gte": {"F": 1}
       }
   )
   
   groups.append(custom_group)
   
   # Use in parsing
   smiles = "FC(F)(F)C(F)(F)C(F)(F)ON(=O)=O"
   results = parse_smiles(smiles, pfas_groups=groups)
   
   if results[0]:
       for group, count, lengths, _ in results[0]:
           print(f"Detected: {group.name}")

Part 6: Real-World Applications
--------------------------------

Database Screening
~~~~~~~~~~~~~~~~~~

Screen a chemical database:

.. code-block:: python

   import pandas as pd
   from PFASgroups import parse_smiles
   from tqdm import tqdm
   
   # Load database
   db = pd.read_csv('chemical_database.csv')
   
   # Screen for PFAS
   pfas_results = []
   
   for idx, row in tqdm(db.iterrows(), total=len(db)):
       smiles = row['SMILES']
       results = parse_smiles(smiles)
       
       if results[0]:
           pfas_results.append({
               'id': row['compound_id'],
               'smiles': smiles,
               'name': row['name'],
               'num_groups': len(results[0]),
               'groups': [g.name for g, _, _, _ in results[0]],
               'chain_lengths': [l for _, _, l, _ in results[0]]
           })
   
   # Save PFAS hits
   pfas_df = pd.DataFrame(pfas_results)
   pfas_df.to_csv('pfas_hits.csv', index=False)
   
   print(f"\nFound {len(pfas_results)} PFAS compounds")
   print(f"Screening rate: {len(db)/len(pfas_results):.1f}:1")

Regulatory Compliance
~~~~~~~~~~~~~~~~~~~~~

Check against PFAS definitions:

.. code-block:: python

   from PFASgroups import parse_smiles, PFASDefinition
   from rdkit import Chem
   
   # Define EPA PFAS definition (example)
   epa_definition = PFASDefinition(
       id=1,
       name="EPA PFAS",
       smarts=["[C](F)(F)F"],
       fluorineRatio=0.0,
       requireBoth=False
   )
   
   # Test compounds
   compounds = [
       ("PFOA", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"),
       ("TFA", "FC(F)(F)C(=O)O"),
       ("Ethanol", "CCO")
   ]
   
   for name, smiles in compounds:
       mol = Chem.MolFromSmiles(smiles)
       is_pfas = epa_definition.applies_to_molecule(mol)
       
       # Also check specific groups
       results = parse_smiles(smiles)
       groups = [g.name for g, _, _, _ in results[0]] if results[0] else []
       
       print(f"\n{name}:")
       print(f"  EPA PFAS: {'Yes' if is_pfas else 'No'}")
       print(f"  Groups: {', '.join(groups) if groups else 'None'}")

Error Handling
~~~~~~~~~~~~~~

Robust processing with error handling:

.. code-block:: python

   from PFASgroups import parse_smiles
   from rdkit import Chem
   import logging
   
   # Set up logging
   logging.basicConfig(level=logging.INFO)
   logger = logging.getLogger(__name__)
   
   def process_smiles_safe(smiles_list):
       """Safely process SMILES with error handling."""
       results = []
       errors = []
       
       for i, smiles in enumerate(smiles_list):
           try:
               # Validate SMILES
               mol = Chem.MolFromSmiles(smiles)
               if mol is None:
                   errors.append((i, smiles, "Invalid SMILES"))
                   continue
               
               # Parse PFAS groups
               result = parse_smiles(smiles)
               results.append({
                   'index': i,
                   'smiles': smiles,
                   'success': True,
                   'groups': result[0] if result else []
               })
               
           except Exception as e:
               logger.error(f"Error processing {smiles}: {e}")
               errors.append((i, smiles, str(e)))
               results.append({
                   'index': i,
                   'smiles': smiles,
                   'success': False,
                   'error': str(e)
               })
       
       return results, errors
   
   # Use it
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",  # Valid PFAS
       "INVALID",  # Invalid SMILES
       "CCO"  # Valid non-PFAS
   ]
   
   results, errors = process_smiles_safe(smiles_list)
   
   print(f"\nProcessed: {len(results)}")
   print(f"Errors: {len(errors)}")

Next Steps
----------

- Explore :doc:`customization` for creating custom groups
- See :doc:`benchmarking` for validation and performance
- Check :doc:`api/core` for complete API reference
- Review :doc:`algorithm` for technical details
