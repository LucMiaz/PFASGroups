Changelog
=========

Version 2.2.3 (Current)
-----------------------

**Released:** February 2026

**New Features:**

- Added a ``ResultsModel`` container around the default ``parse_mols`` / ``parse_smiles``
   output, providing convenient helpers for navigating PFAS group matches and
   components (e.g. ``show()``, ``summarise()``, ``table()``, and plotting
   functions) while remaining fully compatible with existing list-of-dicts
   consumers and JSON export.
- Extended visualisation utilities for PFAS components:
   - ``plot_pfasgroups`` now supports filtering by component path type
      (e.g. perfluoroalkyl vs polyfluoroalkyl), optional SMARTS filtering,
      panel labels, and more robust handling of invalid SMILES.
   - New figure-generation helpers in the codebase reproduce the main-text and
      supplementary figures used in the PFASgroups manuscript (overview example
      and component/path-type illustrations).
- Improved documentation and manuscript alignment:
   - Updated the main article and Additional File~1 descriptions of PFASgroups
      to match the current implementation (114 PFAS groups, component-based
      graph metrics, universal component merging, fluorotelomer linker
      validation, and regulatory definition support).
   - Clarified JSON-based configuration of PFAS groups and regulatory
      definitions and how component SMARTSs control perfluoroalkyl and
      polyfluoroalkyl detection.



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
