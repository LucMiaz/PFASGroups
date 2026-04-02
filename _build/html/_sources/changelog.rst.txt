Changelog
=========

Version 3.1.3
--------------

**Released:** 2026

- Documentation rewrite for ReadTheDocs publication
- Both ``HalogenGroups`` and ``PFASGroups`` entry points documented
- All code examples updated to current API (``ResultsModel`` pattern)
- CLI reference page added
- API reference pages updated with autodoc

Version 3.2.0
--------------

**New features:**

- **``get_compiled_HalogenGroups()``**: Returns compiled :class:`HalogenGroup`
  instances (``compute=True`` groups only), suitable for extending with custom
  groups and passing directly to :func:`parse_smiles` or
  :func:`generate_fingerprint`.  The older :func:`get_HalogenGroups` continues
  to return raw JSON dicts.

  .. code-block:: python

     from HalogenGroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles

     groups = get_compiled_HalogenGroups()
     groups.append(HalogenGroup(
         group_id=200, name="my_custom_group",
         smarts="[CX4](F)(F)(F)",
         category="Custom",
         is_PFAS=True,
     ))
     results = parse_smiles(["FC(F)(F)C(F)(F)F"], halogen_groups=groups)

**API changes:**

- ``generate_fingerprint`` return value is now a 2-D ``numpy.ndarray`` of
  shape ``(n_molecules, n_groups)`` instead of a Python list.

Version 3.1.0
--------------

**New features:**

- **Multi-halogen fingerprinting**: :func:`~HalogenGroups.generate_fingerprint`
  and :meth:`~HalogenGroups.ResultsModel.to_fingerprint` now accept a
  ``halogens`` parameter.  Multiple halogens produce horizontally-stacked
  vectors (116 × n_halogens columns).
- **Saturation filter**: ``saturation`` parameter (``'saturated'`` | ``'unsaturated'``
  | ``None``) for both parsing and fingerprinting.
- **Corrected group-selection indexing**: selections resolve by canonical group
  IDs rather than list positions.
- **``ResultsFingerprint`` metadata**: ``halogens`` and ``saturation``
  stored in the object and shown in ``__repr__`` / ``summary()``.

Version 2.2.4
--------------

**New features:**

- **ResultsFingerprint class** with group selection, multiple encoding modes
  (binary / count / max_component), PCA, t-SNE, UMAP, and KL divergence.
- **Database persistence**: ``to_sql()`` / ``from_sql()`` for ResultsFingerprint
  and ResultsModel (SQLite and PostgreSQL via SQLAlchemy).
- 100+ new unit tests.

Version 2.2.3
--------------

**New features:**

- **ResultsModel container**: wraps ``parse_smiles`` / ``parse_mols`` output
  with helpers ``show()``, ``summarise()``, ``table()``, and ``to_dataframe()``.
- Extended visualisation utilities (``plot_pfasgroups``).
- Updated PFAS group definitions to 116 groups.

Version 2.0
-----------

- Complete rewrite of the SMARTS-based matcher.
- New JSON-based group definition format.
- Regulatory PFAS definition support (OECD, REACH, OPPT, UK EA, PFASSTRUCTv5).
- CLI entry points ``halogengroups`` / ``pfasgroups`` added.