Getting Started
===============

This guide will help you get up and running with PFASgroups quickly.

Installation
------------

Requirements
^^^^^^^^^^^^

- Python >= 3.7
- RDKit
- NumPy
- Pandas
- NetworkX
- tqdm
- svgutils

From PyPI
^^^^^^^^^

The easiest way to install PFASgroups is via pip:

.. code-block:: bash

   pip install PFASgroups

From Conda-Forge
^^^^^^^^^^^^^^^^

If you use conda, you can install from conda-forge:

.. code-block:: bash

   conda install -c conda-forge pfasgroups

From Source
^^^^^^^^^^^

For the latest development version:

.. code-block:: bash

   git clone https://github.com/lucmiaz/PFASGroups.git
   cd PFASGroups
   pip install -e .

Verifying Installation
^^^^^^^^^^^^^^^^^^^^^^

After installation, verify that PFASgroups is correctly installed:

.. code-block:: python

   import PFASgroups
   print(PFASgroups.__version__)

You should also have access to the command-line interface:

.. code-block:: bash

   pfasgroups --help

Quick Start
-----------

Basic Parsing
^^^^^^^^^^^^^

Parse SMILES strings to identify PFAS groups:

.. code-block:: python

   from PFASgroups import parse_smiles

   # Parse a single SMILES
   results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O")

   # Parse multiple SMILES
   smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
   results = parse_smiles(smiles_list)

   # Each result contains: (PFASGroup, match_count, chain_lengths, matched_components)
   for i, smiles_str in enumerate(smiles_list):
       print(f"\nAnalyzing: {smiles_str}")
       for group, n_matches, chain_lengths, components in results[i]:
           print(f"  {group.name}: {n_matches} matches, chains: {chain_lengths}")
           # Access component-level metrics when bycomponent=True
           if components and isinstance(components, list):
               for comp in components:
                   if isinstance(comp, dict):
                       print(f"    Component fraction: {comp.get('component_fraction', 'N/A'):.3f}")
                       print(f"    Branching: {comp.get('branching', 'N/A'):.3f}")

Understanding Results
^^^^^^^^^^^^^^^^^^^^^

The ``parse_smiles`` function returns a list of results, one for each input SMILES.
Each result is a list of tuples containing:

1. **PFASGroup**: The matched PFAS group object with name, ID, and constraints
2. **match_count**: Number of times this group pattern was matched
3. **chain_lengths**: List of carbon chain lengths found
4. **matched_chains**: Detailed information about matched chains (atom indices, pathway type)

Fingerprint Generation
^^^^^^^^^^^^^^^^^^^^^^

Generate PFAS fingerprints for machine learning:

.. code-block:: python

   from PFASgroups import generate_fingerprint

   smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]

   # Binary vector fingerprint (default)
   fingerprints, group_info = generate_fingerprint(smiles_list)
   print(fingerprints.shape)  # (2, 55) - 2 molecules, 55 PFAS groups

   # Dictionary format
   fps_dict, info = generate_fingerprint(smiles_list, representation='dict')
   
   # Count-based fingerprint
   fps_count, info = generate_fingerprint(smiles_list, count_mode='count')

Command Line Interface
^^^^^^^^^^^^^^^^^^^^^^

PFASgroups also provides a convenient CLI:

.. code-block:: bash

   # Parse SMILES
   pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O"

   # Generate fingerprints
   pfasgroups fingerprint "C(C(F)(F)F)F" --format dict

   # List available groups
   pfasgroups list-groups

   # List pathway types
   pfasgroups list-paths

Visualization
^^^^^^^^^^^^^

Visualize PFAS group assignments:

.. code-block:: python

   from PFASgroups import plot_pfasgroups

   # Plot PFAS groups on a molecule
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O")

   # Save as SVG
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O", svg=True, path="output.svg")

Next Steps
----------

- Read the :doc:`user_guide` for comprehensive examples
- Explore the :doc:`api/core` for detailed API documentation
- Learn about :doc:`customization` to extend PFASgroups
- Check :doc:`pfas_groups/oecd_groups` for PFAS classifications
