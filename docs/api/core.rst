Core API
========

.. currentmodule:: PFASGroups

.. code-block:: python

   from PFASGroups import parse_smiles, generate_fingerprint

   # Parse — returns PFASEmbeddingSet
   results = parse_smiles(["CCCC(F)(F)F", "ClCCCl"])

   # Fingerprint — returns (numpy.ndarray, dict)
   fps, info = generate_fingerprint(["CCCC(F)(F)F", "ClCCCl"])
   print(fps.shape)             # (2, 116)   — 116 groups, F only, binary
   print(info['group_names'][:2])

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
     - ``False``
     - Classify each molecule against the five PFAS definitions
   * - ``halogen_groups``
     - ``None``
     - Custom list of :class:`HalogenGroup` instances

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

**Returns:** ``(numpy.ndarray of shape (n_mols, n_cols), dict)``

The second return value is a dict with keys:
``group_names`` (list), ``group_ids`` (list), ``selected_indices`` (list),
``halogens`` (list), ``saturation`` (str or None).

.. code-block:: python

   from PFASGroups import generate_fingerprint

   # Default: 116-col binary, F only
   fps, info = generate_fingerprint(["CCCC(F)(F)F", "ClCCCl"])
   print(fps.shape)                   # (2, 116)
   print(info['group_names'][:2])     # list of group name strings
   print(info['halogens'])            # ['F']

   # OECD groups only (indices 0–27)
   fps_oecd, _ = generate_fingerprint(["CCCC(F)(F)F"], selected_groups=range(0, 28))
   print(fps_oecd.shape)              # (1, 28)

   # Multi-halogen: 4 × 116 = 464 columns
   fps_all, _ = generate_fingerprint(["CCCC(F)(F)F"],
                                      halogens=['F', 'Cl', 'Br', 'I'])
   print(fps_all.shape)               # (1, 464)

Group library
--------------

get_compiled_HalogenGroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_compiled_HalogenGroups

Returns a list of compiled :class:`HalogenGroup` instances
(116 groups with ``compute=True``):

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))      # 116
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