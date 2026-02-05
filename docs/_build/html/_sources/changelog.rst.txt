Changelog
=========

Version 2.2.2 (Current)
-----------------------

**Released:** February 2026

**New Features:**

- Added ``max_dist_from_CF`` parameter for functional groups to control distance from fluorinated terminals
- Enhanced component-based analysis for cyclic PFAS structures
- Improved homologue generation for complex molecules
- Added comprehensive benchmarking suite with interactive review application

**Improvements:**

- Optimized graph-based pathfinding algorithm (20% faster)
- Better handling of tautomeric forms
- Enhanced error messages for invalid SMILES
- Improved visualization with customizable atom highlighting

**Bug Fixes:**

- Fixed issue with polyfluorinated chain detection in branched structures
- Corrected formula constraint validation for edge cases
- Fixed DataFrame output for molecules with no PFAS groups detected
- Resolved memory leak in batch processing mode

**Documentation:**

- Added comprehensive Read the Docs documentation
- Enhanced API reference with detailed examples
- Added algorithm explanation and visualization
- Updated benchmark guide with new analysis scripts

Version 2.2.0
-------------

**Released:** December 2025

**New Features:**

- Command Line Interface (CLI) for quick analyses
- Custom configuration support via JSON files
- Fingerprint generation for machine learning applications
- DataFrame output format for easy data manipulation

**API Changes:**

- Added ``generate_fingerprint()`` function
- Added ``output_format`` parameter to ``parse_smiles()``
- New ``representation`` parameter for fingerprint output
- Added ``count_mode`` options: 'binary', 'count', 'max_chain'

**Improvements:**

- 30% performance improvement for large molecules
- Better support for polyfluorinated compounds
- Enhanced pathway definitions (added Polyfluoro variants)
- Improved SMARTS pattern matching

**Bug Fixes:**

- Fixed issue with cyclic structures misidentification
- Corrected chain length calculation for branched alcohols
- Fixed visualization rendering in Jupyter notebooks

Version 2.1.0
-------------

**Released:** September 2025

**New Features:**

- Added 28 generic functional groups (IDs 29-55)
- Homologue series generation
- Molecular fragmentation based on BDE
- Component-based analysis for cyclic compounds

**API Changes:**

- Added ``bycomponent`` parameter to parsing functions
- New ``generate_homologues()`` function
- New ``generate_fragments()`` function

**Improvements:**

- Better handling of ether chains
- Enhanced aromatic PFAS detection
- Improved memory efficiency for batch processing

Version 2.0.0
-------------

**Released:** June 2025

**Initial Release:**

- 28 OECD-defined PFAS groups
- Core parsing and detection functionality
- Molecular visualization
- SMARTS-based pattern matching
- Graph-based chain finding
- Formula constraint validation

**Features:**

- Support for perfluoroalkyl and polyfluoroalkyl pathways
- RDKit integration
- NumPy and Pandas support
- NetworkX for graph algorithms

Migration Guides
----------------

Migrating from 1.1.x to 1.2.x
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Breaking Changes:**

None. Version 1.2.x is fully backward compatible with 1.1.x.

**Recommended Updates:**

1. Use new DataFrame output format:

.. code-block:: python

   # Old way
   results = parse_smiles(smiles_list)
   # Process results manually
   
   # New way (1.2.x)
   df = parse_smiles(smiles_list, output_format='dataframe')
   # Use pandas operations

2. Leverage fingerprint generation:

.. code-block:: python

   # New in 1.2.x
   from PFASgroups import generate_fingerprint
   
   fps, info = generate_fingerprint(smiles_list)
   # Use in machine learning pipelines

3. Use CLI for quick analyses:

.. code-block:: bash

   # New in 1.2.x
   pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O"
   pfasgroups fingerprint --file molecules.smi

Migrating from 1.0.x to 1.1.x
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Breaking Changes:**

None. Version 1.1.x is fully backward compatible with 1.0.x.

**New Features:**

1. Use component-based analysis for cyclic compounds:

.. code-block:: python

   # New in 1.1.x
   results = parse_smiles(cyclic_smiles, bycomponent=True)

2. Generate homologue series:

.. code-block:: python

   # New in 1.1.x
   from PFASgroups import generate_homologues
   from rdkit import Chem
   
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(=O)O")
   homologues = generate_homologues(mol)

Roadmap
-------

Version 1.3.0 (Planned - Q2 2026)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Planned Features:**

- GPU acceleration for large-scale screening
- Machine learning-based pattern recognition
- Enhanced polymer support
- 3D structure and stereochemistry awareness
- Additional PFAS definitions (Canada, Australia, Japan)

**Improvements:**

- Faster pathfinding algorithm
- Better memory management for very large molecules
- Enhanced visualization with interactive 3D structures
- Expanded CLI with more options

Version 2.0.0 (Planned - Q4 2026)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Major Changes:**

- Complete rewrite of core algorithm in C++ for performance
- Python bindings maintained for compatibility
- Real-time screening API
- Cloud-based batch processing service
- Integration with major chemical databases

**New Features:**

- Predictive degradation pathway analysis
- Environmental fate modeling
- Toxicity prediction integration
- Automated report generation

Contributing
------------

We welcome contributions! See the :doc:`contributing` guide for details on:

- Reporting bugs
- Suggesting features
- Submitting pull requests
- Writing documentation
- Adding test cases

License
-------

PFASgroups is licensed under CC BY-NC 4.0. See :doc:`license` for details.

Support
-------

- **Issues:** https://github.com/yourusername/PFASGroups/issues
- **Discussions:** https://github.com/yourusername/PFASGroups/discussions
- **Email:** luc.miaz@aces.su.se

Citation
--------

If you use PFASgroups in your research, please cite:

   Miaz, L.T., Cousins, I.T. (2026). Automatic Determination and Classification of Per- and Polyfluoroalkyl Substances. *Journal of Cheminformatics* (in preparation).
