Molecule Generation
===================

This module provides utilities for generating random PFAS-like molecules
for testing and validation purposes.

.. module:: PFASgroups.generate_mol
   :synopsis: Random molecule generation utilities

Overview
--------

The molecule generation module is primarily used for:

- Creating synthetic test molecules
- Generating training data for machine learning
- Validating algorithm behavior
- Exploring chemical space

Main Functions
--------------

generate_random_carbon_chain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: generate_random_carbon_chain(n, cycle=False, alkene=False, alkyne=False)

   Generate a random carbon chain molecule with n carbons.

   :param n: Number of carbon atoms
   :type n: int
   :param cycle: Generate cyclic (aromatic) structure
   :type cycle: bool, optional
   :param alkene: Include C=C double bonds
   :type alkene: bool, optional
   :param alkyne: Include C≡C triple bonds
   :type alkyne: bool, optional
   :returns: RDKit molecule object
   :rtype: rdkit.Chem.Mol

   **Example:**

   .. code-block:: python

      from PFASgroups.generate_mol import generate_random_carbon_chain
      from rdkit import Chem

      # Simple chain
      mol = generate_random_carbon_chain(6)
      print(Chem.MolToSmiles(mol))  # e.g., "CCCCCC"

      # Cyclic structure
      mol = generate_random_carbon_chain(6, cycle=True)
      print(Chem.MolToSmiles(mol))  # e.g., "c1ccccc1"

      # With double bonds
      mol = generate_random_carbon_chain(6, alkene=True)
      print(Chem.MolToSmiles(mol))  # e.g., "CC=CCCC"

append_functional_group
^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: append_functional_group(mol, group_smiles, insertion='attach', m=1, atom_indices=None, neighbor_atoms=['C'], sanitize=False)

   Append a functional group to a molecule.

   :param mol: Base molecule
   :type mol: rdkit.Chem.Mol
   :param group_smiles: Functional group as SMILES string
   :type group_smiles: str
   :param insertion: Mode - 'attach' (append) or 'insert' (between atoms)
   :type insertion: str, optional
   :param m: Number of groups to attach
   :type m: int, optional
   :param atom_indices: Specific atom indices for attachment
   :type atom_indices: list, optional
   :param neighbor_atoms: Allowed neighbor atom types
   :type neighbor_atoms: list[str], optional
   :param sanitize: Sanitize molecule after attachment
   :type sanitize: bool, optional
   :returns: Modified molecule
   :rtype: rdkit.Chem.Mol

   **Example:**

   .. code-block:: python

      from PFASgroups.generate_mol import generate_random_carbon_chain, append_functional_group
      from rdkit import Chem

      # Create base chain
      mol = generate_random_carbon_chain(4)

      # Add carboxylic acid
      mol = append_functional_group(mol, "C(=O)O", insertion='attach')
      print(Chem.MolToSmiles(mol))

      # Insert ether linkage
      mol = generate_random_carbon_chain(4)
      mol = append_functional_group(mol, "O", insertion='insert')
      print(Chem.MolToSmiles(mol))

get_attachment
^^^^^^^^^^^^^^

.. py:function:: get_attachment(mol, m, atom_symbols=['C'], neighbors_symbols={'C': ['F', 'H']})

   Find attachment points in a molecule for functional group addition.

   :param mol: Molecule to analyze
   :type mol: rdkit.Chem.Mol
   :param m: Maximum number of attachment points to return
   :type m: int
   :param atom_symbols: Allowed central atom types
   :type atom_symbols: list[str], optional
   :param neighbors_symbols: Allowed neighbor types for each atom type
   :type neighbors_symbols: dict, optional
   :returns: List of (atom_index, neighbor_index) tuples
   :rtype: list[tuple]

attach_mol
^^^^^^^^^^

.. py:function:: attach_mol(mol, submol, atom_index)

   Attach a submolecule to a specific atom in the main molecule.

   :param mol: Main molecule
   :type mol: rdkit.Chem.Mol
   :param submol: Submolecule to attach
   :type submol: rdkit.Chem.Mol
   :param atom_index: Atom index in main molecule for attachment
   :type atom_index: int
   :returns: Combined molecule
   :rtype: rdkit.Chem.Mol

   The submolecule is attached at its first atom (index 0).

insert_mol
^^^^^^^^^^

.. py:function:: insert_mol(mol, group_mol, atom_index, neighbor_index)

   Insert a submolecule between two bonded atoms.

   :param mol: Main molecule
   :type mol: rdkit.Chem.Mol
   :param group_mol: Group to insert
   :type group_mol: rdkit.Chem.Mol
   :param atom_index: First atom of the bond
   :type atom_index: int
   :param neighbor_index: Second atom of the bond
   :type neighbor_index: int
   :returns: Modified molecule
   :rtype: rdkit.Chem.Mol

remove_atoms
^^^^^^^^^^^^

.. py:function:: remove_atoms(mol, idxs, removable=['H', 'F', 'Cl', 'Br', 'I'], show_on_error=False)

   Remove atoms and reconnect the molecular structure.

   :param mol: Molecule to modify
   :type mol: rdkit.Chem.Mol
   :param idxs: Atom indices to remove
   :type idxs: list[int]
   :param removable: Elements that can be removed along with specified atoms
   :type removable: list[str], optional
   :param show_on_error: Display molecule on error
   :type show_on_error: bool, optional
   :returns: Modified molecule
   :rtype: rdkit.Chem.Mol

Use Cases
---------

Generating Test Molecules
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups.generate_mol import generate_random_carbon_chain, append_functional_group
   from rdkit import Chem
   import numpy as np

   def generate_random_pfas(n_carbons=6, n_fluorines=None):
       """Generate a random PFAS-like molecule."""
       # Create carbon chain
       mol = generate_random_carbon_chain(n_carbons)
       
       # Add fluorines (if not specified, use random amount)
       if n_fluorines is None:
           n_fluorines = np.random.randint(1, n_carbons * 2)
       
       for _ in range(n_fluorines):
           mol = append_functional_group(mol, "F", insertion='attach')
       
       # Add functional group
       functional_groups = ["C(=O)O", "S(=O)(=O)O", "O", "N"]
       fg = np.random.choice(functional_groups)
       mol = append_functional_group(mol, fg, insertion='attach')
       
       return mol

   # Generate 10 random PFAS molecules
   molecules = [generate_random_pfas() for _ in range(10)]
   for mol in molecules:
       print(Chem.MolToSmiles(mol))

Synthetic Validation Dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups.generate_mol import generate_random_carbon_chain, append_functional_group
   from PFASgroups import parse_smiles
   from rdkit import Chem

   def generate_validation_set():
       """Generate molecules for each PFAS group."""
       validation_molecules = []
       
       # PFCAs: perfluoroalkyl carboxylic acids
       for n in range(4, 10):
           mol = generate_random_carbon_chain(n)
           # Fluorinate all carbons except the terminal one
           for i in range(n - 1):
               mol = append_functional_group(mol, "F", m=2)
           mol = append_functional_group(mol, "C(=O)O")
           validation_molecules.append(("PFCA", Chem.MolToSmiles(mol)))
       
       return validation_molecules

   validation = generate_validation_set()
   
   # Verify classification
   for expected_group, smiles in validation:
       results = parse_smiles([smiles])
       detected = [g.name for g, _, _, _ in results[0]]
       print(f"{expected_group}: {smiles} -> {detected}")

Notes and Limitations
---------------------

1. **Random seed**: Set ``np.random.seed()`` for reproducibility
2. **Valence checking**: Generated molecules are sanitized but may have unusual structures
3. **Not exhaustive**: Does not generate all possible structures
4. **Test purposes**: Primarily intended for testing, not for generating real compounds
5. **Stereochemistry**: Does not handle stereochemistry
