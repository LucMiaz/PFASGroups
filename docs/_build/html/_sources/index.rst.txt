HalogenGroups / PFASGroups
==========================

.. image:: https://badge.fury.io/py/PFASgroups.svg
   :target: https://badge.fury.io/py/PFASgroups
   :alt: PyPI version

.. image:: https://img.shields.io/badge/python-3.7+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.7+

.. image:: https://img.shields.io/badge/License-CC%20BY--ND%204.0-lightgrey.svg
   :target: https://creativecommons.org/licenses/by-nd/4.0/
   :alt: License: CC BY-ND 4.0

.. image:: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg
   :target: https://www.rdkit.org/
   :alt: Powered by RDKit

A comprehensive Python cheminformatics package for automated detection,
classification, and analysis of halogenated substances — with a focus on
Per- and Polyfluoroalkyl Substances (PFAS).

The package can be imported as **HalogenGroups** (all halogens: F, Cl, Br, I)
or as **PFASGroups** (fluorine-only, for backward compatibility):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Import
     - Default behavior
   * - ``from HalogenGroups import parse_smiles``
     - Searches for all halogens (F, Cl, Br, I) across 116 groups
   * - ``from PFASGroups import parse_smiles``
     - Searches for fluorine only (backward-compatible with v1/v2)

Quick example
-------------

.. code-block:: python

   from HalogenGroups import parse_smiles, generate_fingerprint
   import numpy as np

   smiles = ["CCCC(F)(F)F", "ClCCCl", "FC(F)(F)C(=O)O"]
   results = parse_smiles(smiles)

   for mol_result in results:
       print(mol_result.smiles, "->", len(mol_result.matches), "group matches")

   fps, group_names = generate_fingerprint(smiles)
   print(fps.shape)   # (3, 464)

Key capabilities
----------------

- **116 halogen groups**: 28 OECD, 45 generic, 43 fluorotelomer groups
- **Dual-halogen API**: parse for any combination of F, Cl, Br, I in one call
- **PFAS definition screening**: classify molecules against 5 regulatory frameworks
- **Fingerprinting**: generate 464-column multi-halogen fingerprints for ML
- **Dimensionality reduction**: PCA, t-SNE, UMAP on fingerprint matrices
- **Molecule prioritization**: rank molecules by structural novelty
- **Command-line interface**: ``halogengroups parse``, ``fingerprint``, ``list-groups``

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Core Concepts

   halogengroups
   pfas_definitions
   algorithm

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics

   tutorial
   customization
   prioritization
   benchmarking
   cli

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/core
   api/models
   api/fingerprint_analysis

.. toctree::
   :maxdepth: 1
   :caption: Project Info

   changelog
   contributing
   license

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`