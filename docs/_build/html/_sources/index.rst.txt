PFASgroups Documentation
========================

.. image:: https://badge.fury.io/py/PFASgroups.svg
   :target: https://badge.fury.io/py/PFASgroups
   :alt: PyPI version

.. image:: https://img.shields.io/badge/python-3.7+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.7+

.. image:: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
   :target: https://creativecommons.org/licenses/by-nc/4.0/
   :alt: License: CC BY-NC 4.0

.. image:: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg
   :target: https://www.rdkit.org/
   :alt: Powered by RDKit

A comprehensive Python cheminformatics package for automated detection, classification, and analysis of Per- and Polyfluoroalkyl Substances (PFAS) in chemical databases.

Overview
--------

**PFASgroups** provides three main capabilities:

1. **PFAS Group Identification and Chain Analysis**: Automated detection of 55 functional groups directly connected to fluorinated chains, based on OECD terminology
2. **Fluorinated Chain Length Quantification**: Determination of per- or polyfluorinated alkyl chain lengths connected to functional groups
3. **Customization**: Both functional groups and pathway patterns can be customized via JSON configuration files

The algorithm combines SMARTS pattern matching, molecular formula constraints, and graph-based pathfinding using RDKit and NetworkX.

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash

   pip install PFASgroups

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from PFASgroups import parse_smiles, generate_fingerprint

   # Parse PFAS structures
   smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
   results = parse_smiles(smiles_list)

   # Each result contains: (PFASGroup, match_count, chain_lengths, matched_chains)
   for i, smiles_str in enumerate(smiles_list):
       print(f"\nAnalyzing: {smiles_str}")
       for group, n_matches, chain_lengths, chains in results[i]:
           print(f"  {group.name}: {n_matches} matches, chains: {chain_lengths}")

   # Generate fingerprints for machine learning
   fingerprints, group_info = generate_fingerprint(smiles_list)

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   tutorial
   algorithm
   pfas_groups

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/core
   api/models
   api/homologues
   api/fragmentation
   api/visualization
   api/generation

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics

   customization
   benchmarking
   cli

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   contributing
   changelog
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
