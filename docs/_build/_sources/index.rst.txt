PFASGroups Documentation
=========================

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

**PFASGroups** is a Python cheminformatics package for automated detection,
classification, and analysis of Per- and Polyfluoroalkyl Substances (PFAS)
using SMARTS-based structural group matching.

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OCCOCCO"])

   for mol in results:
       print(mol.smiles, "—", len(mol.matches), "group match(es)")
   # CCCC(F)(F)F — 1 group match(es)
   # FC(F)(F)C(=O)O — 1 group match(es)
   # OCCOCCO — 0 group match(es)

Key capabilities
----------------

- **116 halogen groups**: 28 OECD, 45 generic, 43 fluorotelomer groups
- **PFAS definition screening**: classify molecules against 5 regulatory frameworks
- **Fingerprinting**: generate 116-column fingerprints (binary, count  or max-component) for ML
- **Group selection**: use all groups, OECD only, generic, telomers, or a custom subset
- **Dimensionality reduction**: PCA, t-SNE, UMAP on fingerprint matrices
- **Molecule prioritization**: rank molecules by structural novelty
- **Command-line interface**: ``pfasgroups parse``, ``fingerprint``, ``list-groups``
- **Multi-halogen support**: extend detection to Cl, Br, I via :doc:`halogengroups`

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Core Concepts

   pfas_definitions
   algorithm

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics

   tutorial
   customization
   prioritization
   benchmarking
   halogengroups
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