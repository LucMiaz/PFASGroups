.. PFASgroups documentation master file

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
   :alt: License

.. image:: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg
   :target: https://www.rdkit.org/
   :alt: Powered by RDKit

**PFASgroups** is a comprehensive Python cheminformatics package for automated detection, 
classification, and analysis of Per- and Polyfluoroalkyl Substances (PFAS) in chemical databases.

Key Features
------------

- **PFAS Group Identification**: Automated detection of 55 functional groups based on OECD terminology
- **Chain Length Analysis**: Quantification of per- or polyfluorinated alkyl chain lengths
- **Comprehensive Metrics**: Graph-theoretic analysis (diameter, radius, center, periphery, branching)
- **Molecular Coverage Analysis**: Quantify fluorination extent with component fraction metrics
- **Fingerprint Generation**: Generate PFAS fingerprints for machine learning applications
- **Homologue Generation**: Create shorter-chain analogues by removing CF₂ moieties
- **Visualization**: Plot PFAS group assignments on molecular structures
- **Customization**: Extend with custom PFAS groups and pathway patterns via JSON

Quick Example
-------------

.. code-block:: python

   from PFASgroups import parse_smiles, generate_fingerprint

   # Parse PFAS structures with comprehensive metrics
   smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
   results = parse_smiles(smiles_list, bycomponent=True)
   
   # Access metrics
   for result in results:
       for match in result['matches']:
           print(f"Group: {match['group_name']}")
           print(f"  Mean branching: {match['mean_branching']}")
           print(f"  Component fraction: {match['mean_component_fraction']}")
           print(f"  Total coverage: {match['total_components_fraction']}")

   # Generate fingerprints for machine learning
   fingerprints, group_info = generate_fingerprint(smiles_list)

Installation
------------

From PyPI:

.. code-block:: bash

   pip install PFASgroups

From Conda-Forge:

.. code-block:: bash

   conda install -c conda-forge pfasgroups

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   getting_started
   user_guide
   cli
   customization

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
   :caption: PFAS Classifications

   pfas_definitions
   pfas_groups/oecd_groups
   pfas_groups/generic_groups
   pfas_groups/pathway_types

.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   algorithm
   testing
   changelog
   contributing
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Acknowledgments
---------------

This project is part of the `ZeroPM project <https://zeropm.eu/>`_ (WP2) and has received 
funding from the European Union's Horizon 2020 research and innovation programme under 
grant agreement No 101036756.

Developed at the `Department of Environmental Science <https://aces.su.se>`_ at Stockholm University.
