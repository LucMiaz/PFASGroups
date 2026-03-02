Quick Start Guide
=================

This guide will get you up and running with PFASgroups in minutes.

Basic PFAS Detection
---------------------

Detecting PFAS groups in molecules:

.. code-block:: python

   from HalogenGroups import parse_smiles
   
   # Single SMILES string
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   results = parse_smiles(smiles)
   
   # results is a ResultsModel (list of MoleculeResult)
   for match in results[0].matches:
       print(f"{match.group_name}: {len(match.components)} component(s)")
       sizes = [len(c.atoms) for c in match.components]
       print(f"  Component sizes: {sizes}")
       print(f"  OECD Group ID: {match.group_id}")

**Output:**

.. code-block:: text

   Perfluoroalkyl carboxylic acids: 1 matches
     Chain lengths: [5]
     OECD Group ID: 1

Multiple Molecules
------------------

Process a list of SMILES:

.. code-block:: python

   from HalogenGroups import parse_smiles
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",  # PFBA
       "FC(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
       "C(C(F)(F)F)F"  # Not a typical PFAS functional group
   ]
   
   results = parse_smiles(smiles_list)
   
   for mol_result in results:
       print(f"\n{mol_result.smiles}:")
       group_matches = [m for m in mol_result.matches if m.is_group]
       if group_matches:
           for match in group_matches:
               sizes = [len(c.atoms) for c in match.components]
               print(f"  ✓ {match.group_name} (component sizes: {sizes})")
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

   from HalogenGroups.fingerprints import generate_fingerprint
   import numpy as np
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFBA
       "FC(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
       "C(C(F)(F)C(F)(F)C(F)(F)F)O"  # FTOH
   ]
   
   # Default: all halogens, per-saturation → 116 groups × 4 = 464 columns
   fps, group_info = generate_fingerprint(smiles_list)
   print(f"Fingerprint shape: {fps.shape}")   # (3, 464)
   
   # Fluorine only → 116 columns
   fps_f, info = generate_fingerprint(smiles_list, halogens='F')
   print(f"Fluorine-only shape: {fps_f.shape}")  # (3, 116)
   print(info['group_names'][:3])              # ['... [F]', '... [F]', ...]
   
   # Polyfluorinated components only
   fps_poly, _ = generate_fingerprint(smiles_list, halogens='F', saturation='poly')
   
   # Use in ML model
   from sklearn.ensemble import RandomForestClassifier
   X = fps_f   # shape (3, 116)
   y = [1, 1, 0]
   model = RandomForestClassifier()
   model.fit(X, y)

**Different fingerprint representations:**

.. code-block:: python

   # Dictionary format (sparse)
   fps_dict, _ = generate_fingerprint(smiles_list, representation='dict')
   print(fps_dict[0])  # {'Perfluoroalkyl carboxylic acids': 1, ...}
   
   # Count-based (instead of binary)
   fps_count, _ = generate_fingerprint(smiles_list, count_mode='count')
   
   # Max chain length
   fps_chain, _ = generate_fingerprint(smiles_list, count_mode='max_component')

**Using ResultsModel.to_fingerprint():**

.. code-block:: python

   from HalogenGroups import parse_smiles
   
   results = parse_smiles(smiles_list)
   
   # Default fingerprint (all halogens, per-saturation, 116 groups × 4 = 464 columns)
   fp = results.to_fingerprint()
   
   # Fluorine only → shape (n, 116)
   fp_f = results.to_fingerprint(halogens='F')
   
   # OECD groups only, polyfluorinated
   fp_oecd = results.to_fingerprint(
       group_selection='oecd', halogens='F', saturation='poly')
   
   print(repr(fp))
   # ResultsFingerprint(n_molecules=3, n_groups=464,
   #   group_selection='all', halogens=['F', 'Cl', 'Br', 'I'], saturation='per', count_mode='binary')

Visualization
-------------

Visualize PFAS group assignments:

.. code-block:: python

   from HalogenGroups import parse_smiles
   
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   results = parse_smiles(smiles)
   
   # Display in Jupyter notebook
   results[0].show(display=True)
   
   # Save to file
   img = results[0].show(display=False)
   img.save("pfba_annotated.png")
   
   # Export to SVG
   results[0].svg("pfba_annotated.svg")
   
   # Multiple molecules in grid
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "C(C(F)(F)C(F)(F)F)O"
   ]
   results = parse_smiles(smiles_list)
   grid_img = results.show(display=False, ncols=3)
   grid_img.save("pfas_comparison.png")

Component-Based Analysis
------------------------

For molecules with disconnected or cyclic components:

.. code-block:: python

   from HalogenGroups import parse_smiles
   
   # Complex molecule with cycles
   smiles = "C1(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)(F)"
   
   # Analysis - all components are detected automatically
   results = parse_smiles(smiles)
   
   print("Matched groups:", len(results[0].matches))
   print("Total components:", sum(len(m.components) for m in results[0].matches))

Command Line Interface
----------------------

PFASgroups provides a CLI for quick analyses:

**Parse SMILES:**

.. code-block:: bash

   halogengroups parse "FC(F)(F)C(F)(F)C(=O)O"
   
   # Multiple SMILES
   halogengroups parse "FC(F)(F)C(F)(F)C(=O)O" "FC(F)(F)C(F)(F)S(=O)(=O)O"
   
   # From file
   halogengroups parse --file molecules.smi

**Generate fingerprints:**

.. code-block:: bash

   halogengroups fingerprint "FC(F)(F)C(F)(F)C(=O)O" --format dict
   
   # Save to file
   halogengroups fingerprint --file molecules.smi --output fingerprints.json

**List available groups:**

.. code-block:: bash

   halogengroups list-groups
   
   # Filter by main group
   halogengroups list-groups --main "Perfluoroalkyl acids"

**List pathway types:**

.. code-block:: bash

   halogengroups list-paths

Working with DataFrames
------------------------

Output results as pandas DataFrame:

.. code-block:: python

   from HalogenGroups import parse_smiles
   import pandas as pd
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "C(C(F)(F)C(F)(F)F)O"
   ]
   
   results = parse_smiles(smiles_list)
   
   # Build a DataFrame from results
   rows = []
   for mol_result in results:
       for match in mol_result.matches:
           if match.is_group:
               rows.append({
                   'smiles': mol_result.smiles,
                   'group_name': match.group_name,
                   'group_id': match.group_id,
                   'n_components': len(match.components),
               })
   df = pd.DataFrame(rows)
   print(df.head())
   
   # Export to CSV
   df.to_csv('pfas_analysis.csv', index=False)

Reading from Files
------------------

Process SMILES from various file formats:

**From SMILES file (.smi):**

.. code-block:: python

   from HalogenGroups import parse_smiles
   import pandas as pd
   
   # Read SMILES file
   with open('molecules.smi', 'r') as f:
       smiles_list = [line.strip() for line in f if line.strip()]
   
   results = parse_smiles(smiles_list)

**From CSV:**

.. code-block:: python

   import pandas as pd
   from HalogenGroups import parse_smiles
   
   # Read CSV with SMILES column
   df = pd.read_csv('compounds.csv')
   smiles_list = df['SMILES'].tolist()
   
   results = parse_smiles(smiles_list)
   
   # Add results to DataFrame
   df['pfas_groups'] = [
       [m.group_name for m in mol_result.matches if m.is_group]
       for mol_result in results
   ]

**From SDF:**

.. code-block:: python

   from rdkit import Chem
   from HalogenGroups import parse_groups_in_mol
   
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

   from HalogenGroups import parse_smiles
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
           mol_result = results[0]
           group_matches = [m for m in mol_result.matches if m.is_group]
           
           if group_matches:
               print(f"✓ {smiles}: {len(group_matches)} groups")
           else:
               print(f"○ {smiles}: No PFAS groups")
               
       except Exception as e:
           print(f"✗ Error processing {smiles}: {e}")

Filtering Components by Halogen, Form, and Saturation
------------------------------------------------------

You can filter component matches by specific halogens, chemical forms, or saturation levels:

.. code-block:: python

   from HalogenGroups import parse_smiles
   
   smiles_list = [
       "C(C(F)(F)F)F",  # Simple fluorinated
       "FC(F)(F)C(F)(F)C(=O)O",  # Perfluorinated carboxylic acid
       "C(C(Cl)(Cl)Cl)Cl"  # Chlorinated
   ]
   
   # Filter only fluorine components
   results_f = parse_smiles(smiles_list, halogens='F')
   
   # Filter perfluorinated alkyl compounds
   results_pfa = parse_smiles(
       smiles_list,
       halogens='F',
       saturation='per',
       form='alkyl'
   )
   
   # Filter polyfluorinated cyclic compounds
   results_polyf_cyclic = parse_smiles(
       smiles_list,
       halogens='F',
       saturation='poly',
       form='cyclic'
   )
   
   # Filter multiple halogens (F and Cl)
   results_multi = parse_smiles(smiles_list, halogens=['F', 'Cl'])

**Valid filter options:**

- **halogens**: ``'F'``, ``'Cl'``, ``'Br'``, ``'I'``, or list like ``['F', 'Cl']`` (``None`` for all)
- **saturation**: ``'per'`` or ``'poly'`` (or list like ``['per', 'poly']`` for both, ``None`` for all)
- **form**: ``'alkyl'`` or ``'cyclic'`` (or list like ``['alkyl', 'cyclic']`` for both, ``None`` for all)

Command-line equivalents:

.. code-block:: bash

   # Filter for fluorine components only
   halogengroups parse --halogens F "C(C(F)(F)F)F"
   
   # Filter for perfluorinated alkyl
   halogengroups parse --halogens F --saturation per --form alkyl "FC(F)(F)C(F)(F)C(=O)O"
   
   # Filter for multiple halogens
   halogengroups parse --halogens F Cl "C(C(F)(F)F)Cl"
   
   # Filter for cyclic forms only
   halogengroups parse --form cyclic "your_smiles_here"

Next Steps
----------

- Read the :doc:`tutorial` for in-depth examples
- Learn about the :doc:`algorithm` behind the detection
- Explore :doc:`customization` for advanced use cases
- Check the :doc:`api/core` for complete API reference
- See :doc:`benchmarking` for validation and performance analysis
