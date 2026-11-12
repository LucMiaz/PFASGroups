Core API
========

.. currentmodule:: PFASGroups

.. code-block:: python

   from PFASGroups import parse_smiles, generate_fingerprint

   # Parse — returns PFASEmbeddingSet
   results = parse_smiles(["CCCC(F)(F)F", "ClCCCl"])

   # Fingerprint — returns (EmbeddingArray, list of column names)
   fps, cols = generate_fingerprint(["CCCC(F)(F)F", "ClCCCl"])
   print(fps.shape)             # (2, 118)   — one binary column per group
   print(cols[:2])

.. contents:: Contents
   :local:
   :depth: 2

Parsing
-------

parse_smiles
~~~~~~~~~~~~

.. autofunction:: parse_smiles

**Parameters:**

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Description
   * - ``smiles``
     - *required*
     - Single SMILES string or list of SMILES strings
   * - ``halogens``
     - ``'F'``
     - Halogen(s) to detect; string or list of strings
   * - ``saturation``
     - ``None``
     - ``'per'``, ``'poly'``, or ``None`` (all)
   * - ``form``
     - ``None``
     - Structural form filter; ``None`` = all forms
   * - ``compute_component_metrics``
     - ``True``
     - Compute effective graph resistance and related per-component metrics
   * - ``limit_effective_graph_resistance``
     - ``None``
     - Max atoms for which EGR is computed; ``None`` = unlimited
   * - ``include_PFAS_definitions``
     - ``True``
     - Classify each molecule against the six PFAS definitions (fluorine only)
   * - ``pfas_groups``
     - all compiled groups
     - Custom list of :class:`HalogenGroup` instances.  Unknown keywords are
       ignored, so a misspelt name (e.g. ``halogen_groups``) silently falls
       back to the default groups

**Returns:** :class:`PFASEmbeddingSet`

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "ClCCCl"])
   for mol in results:
       print(mol.smiles, len(mol.matches))

parse_mols
~~~~~~~~~~

.. autofunction:: parse_mols

Like :func:`parse_smiles` but accepts RDKit ``Mol`` objects directly:

.. code-block:: python

   from rdkit import Chem
   from PFASGroups import parse_mols

   mols = [Chem.MolFromSmiles(s) for s in ["CCCC(F)(F)F", "ClCCCl"]]
   results = parse_mols(mols)

parse_mol
~~~~~~~~~

.. autofunction:: parse_mol

Parse a single RDKit molecule.  Returns a :class:`MoleculeResult`.

Fingerprinting
--------------

generate_fingerprint
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: generate_fingerprint

**Parameters:**

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default
     - Description
   * - ``smiles``
     - *required*
     - Single SMILES string or list of SMILES strings
   * - ``selected_groups``
     - ``None``
     - 0-based indices of groups to include (list, range, or ``None`` = all)
   * - ``representation``
     - ``'vector'``
     - ``'vector'``, ``'dict'``, ``'sparse'``, ``'detailed'``, or ``'int'``
   * - ``component_metrics``
     - ``['binary']``
     - List of metrics: count modes (``'binary'``, ``'count'``, ``'max_component'``, ``'total_component'``) and/or graph metrics (e.g. ``'effective_graph_resistance'``)
   * - ``halogens``
     - ``'F'``
     - Halogen(s) to include; string or list of strings
   * - ``saturation``
     - ``'per'``
     - ``'per'``, ``'poly'``, or ``None``

**Returns:** ``(EmbeddingArray of shape (n_mols, n_cols), list of column names)``

``generate_fingerprint`` parses the SMILES with ``halogens`` and
``saturation`` and returns ``results.to_array(component_metrics=...)``, one
column per group and metric.  For group selection or one block of columns per
halogen, parse once and call :meth:`PFASEmbeddingSet.to_array` (below);
``selected_groups`` and ``representation`` are currently accepted but not used.

.. code-block:: python

   from PFASGroups import generate_fingerprint

   from PFASGroups import parse_smiles

   # Default: binary, one column per group
   fps, cols = generate_fingerprint(["CCCC(F)(F)F", "ClCCCl"])
   print(fps.shape)                   # (2, 118)
   print(cols[:2])                    # ['OECD perfluoroalkyl alcohols [binary]', ...]

   # Count or max-component encoding
   fps_max, _ = generate_fingerprint(["CCCC(F)(F)F"], count_mode="max_component")

   # OECD groups only, or one block per halogen: use to_array
   results = parse_smiles(["CCCC(F)(F)F"], halogens=['F', 'Cl', 'Br', 'I'])
   print(results.to_array(group_selection='oecd').shape)            # (1, 28)
   print(results.to_array(halogens=['F', 'Cl', 'Br', 'I']).shape)   # (1, 468)

Group library
--------------

get_compiled_HalogenGroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_compiled_HalogenGroups

Returns a list of compiled :class:`HalogenGroup` instances
(118 groups with ``compute=True``):

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))      # 118
   print(groups[0].name)

get_HalogenGroups
~~~~~~~~~~~~~~~~~

.. autofunction:: get_HalogenGroups

Returns raw JSON-like dicts (internal format). Prefer
:func:`get_compiled_HalogenGroups` for most uses.

load_HalogenGroups
~~~~~~~~~~~~~~~~~~

.. autofunction:: load_HalogenGroups

Compile a list of raw group dicts into :class:`HalogenGroup` instances.

get_PFASDefinitions
~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_PFASDefinitions

Returns the list of :class:`PFASDefinition` objects.

get_componentSMARTSs
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_componentSMARTSs

Returns component-level SMARTS patterns used for saturation filtering.

Molecule prioritization
------------------------

prioritise_molecules
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: prioritise_molecules

See :doc:`../prioritization` for detailed documentation.