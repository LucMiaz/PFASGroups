Changelog
=========

Version History
---------------

Version 1.2.2 (Current)
^^^^^^^^^^^^^^^^^^^^^^^

*Release Date: 2024*

**New Features:**

- Added ``max_dist_from_CF`` parameter for controlling functional group distance
  from fluorinated carbons
- New PFAS group definitions added
- Improved documentation with ReadTheDocs support

**Bug Fixes:**

- Fixed edge case in cyclic structure detection
- Improved handling of multi-chain molecules
- Fixed fingerprint generation for empty results

**Documentation:**

- Comprehensive API documentation
- User guide with examples
- Algorithm overview section

Version 1.2.1
^^^^^^^^^^^^^

**Changes:**

- Performance optimizations for large datasets
- Updated OECD group definitions
- Bug fixes for constraint validation

Version 1.2.0
^^^^^^^^^^^^^

**New Features:**

- ``bycomponent`` mode for complex structure analysis
- Generic PFAS groups (IDs 29-57)
- Fingerprint generation for machine learning
- SVG output for visualization

**API Changes:**

- ``parse_smiles`` now returns chain length information
- New ``generate_fingerprint`` function
- Updated output format for ``parse_mols``

Version 1.1.0
^^^^^^^^^^^^^

**New Features:**

- PFAS definition matching (EPA, OECD criteria)
- Homologue generation module
- Fragmentation utilities
- Command-line interface

**Dependencies:**

- Added NetworkX for graph operations
- Added tqdm for progress bars
- Added svgutils for SVG manipulation

Version 1.0.0
^^^^^^^^^^^^^

*Initial Release*

**Features:**

- Core PFAS group classification (28 OECD groups)
- SMARTS-based pattern matching
- Four pathway types (Perfluoroalkyl, Polyfluoroalkyl, Polyfluoro, Polyfluorobr)
- Basic visualization

**Dependencies:**

- RDKit >= 2021.03
- NumPy >= 1.19
- Pandas >= 1.0

Upgrade Guide
-------------

From 1.1.x to 1.2.x
^^^^^^^^^^^^^^^^^^^

**API Changes:**

.. code-block:: python

   # Old (1.1.x)
   results = parse_smiles(smiles)
   for mol_result in results:
       for group, count in mol_result:
           print(group.name)

   # New (1.2.x)
   results = parse_smiles(smiles)
   for mol_result in results:
       for group, count, chains, matched in mol_result:
           print(f"{group.name}: chains={chains}")

**New Parameters:**

- Add ``bycomponent=True`` for complex molecules
- Use ``max_dist_from_CF`` in custom group definitions

From 1.0.x to 1.1.x
^^^^^^^^^^^^^^^^^^^

**New Dependencies:**

.. code-block:: bash

   pip install networkx tqdm svgutils

**New Features to Adopt:**

.. code-block:: python

   # Use PFAS definitions
   from PFASgroups import get_PFASDefinitions

   # Generate homologues
   from PFASgroups import generate_homologues

Deprecation Notices
-------------------

**Deprecated in 1.2.0:**

- ``output_format='csv'`` - Use pandas DataFrame instead
- Direct access to ``PFASGroup.smarts1_string`` - Use compiled pattern

**Planned for 2.0.0:**

- Remove legacy output formats
- Unified API for all parsing functions
- Async support for batch processing

Roadmap
-------

**Planned Features:**

- Stereochemistry support
- Polymer structure handling
- Machine learning integration
- Web API service
- Additional regulatory definitions

**Under Consideration:**

- GPU acceleration for large datasets
- Extended halogen support (iodine, astatine)
- Custom pathway definition GUI
