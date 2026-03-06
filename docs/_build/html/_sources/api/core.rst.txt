Core API
========

.. currentmodule:: HalogenGroups

This page covers the top-level functions exported by both the
``HalogenGroups`` and ``PFASGroups`` namespaces.

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
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``smiles``
     - *required*
     - Single SMILES string or list of SMILES strings
   * - ``halogens``
     - ``'F'`` (PFASGroups) / ``['F','Cl','Br','I']`` (HalogenGroups)
     - Halogen(s) to include; string or list of strings
   * - ``saturation``
     - ``None``
     - ``'saturated'``, ``'unsaturated'``, or ``None`` (all)
   * - ``form``
     - ``None``
     - Structural form filter; ``None`` = all forms
   * - ``compute_component_metrics``
     - ``True``
     - Compute effective graph resistance and related per-component metrics
   * - ``limit_effective_graph_resistance``
     - ``None``
     - Maximum number of atoms for which effective graph resistance is
       computed; ``None`` = unlimited
   * - ``include_PFAS_definitions``
     - ``False``
     - If ``True``, classify each molecule against the five PFAS definitions
   * - ``halogen_groups``
     - ``None``
     - Custom list of :class:`HalogenGroup` instances; ``None`` = use the
       built-in 116-group library
   * - ``n_jobs``
     - ``1``
     - Number of parallel workers (requires ``joblib``)

**Returns:** :class:`ResultsModel`

**Example:**

.. code-block:: python

   from HalogenGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "ClCCCl"])
   for mol in results:
       print(mol.smiles, len(mol.matches))

parse_mols
~~~~~~~~~~

.. autofunction:: parse_mols

Like :func:`parse_smiles` but accepts RDKit ``Mol`` objects directly:

.. code-block:: python

   from rdkit import Chem
   from HalogenGroups import parse_mols

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
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``smiles``
     - *required*
     - List of SMILES strings
   * - ``halogens``
     - ``'F'`` / ``['F','Cl','Br','I']``
     - Halogen(s) for fingerprint columns
   * - ``saturation``
     - ``None``
     - Saturation filter
   * - ``count_mode``
     - ``'max_component'``
     - ``'max_component'``, ``'all'``, or ``'binary'``
   * - ``representation``
     - ``'vector'``
     - Only ``'vector'`` is currently supported (returns numpy array)

**Returns:** ``(numpy.ndarray of shape (n_mols, n_groups × n_halogens), dict)``

.. code-block:: python

   from HalogenGroups import generate_fingerprint

   fps, info = generate_fingerprint(["CCCC(F)(F)F", "ClCCCl"])
   print(fps.shape)        # (2, 464)
   print(info[0])          # ('perfluoroalkyl', 'F')

Group library
--------------

get_compiled_HalogenGroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_compiled_HalogenGroups

Returns a list of compiled :class:`HalogenGroup` instances (only groups with
``compute=True``; 116 groups total).

.. code-block:: python

   from HalogenGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))      # 116
   print(groups[0].name)

get_HalogenGroups
~~~~~~~~~~~~~~~~~

.. autofunction:: get_HalogenGroups

Returns raw JSON-like dicts (internal format).  Prefer
:func:`get_compiled_HalogenGroups` for most uses.

load_HalogenGroups
~~~~~~~~~~~~~~~~~~

.. autofunction:: load_HalogenGroups

Compile a list of raw group dicts into :class:`HalogenGroup` instances.

get_PFASDefinitions
~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_PFASDefinitions

Returns the list of :class:`PFASDefinition` objects used by
``include_PFAS_definitions=True``.

get_componentSmartss
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: get_componentSmartss

Returns the component-level SMARTS patterns used for saturation filtering.

Molecule prioritization
------------------------

prioritise_molecules
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: prioritise_molecules

See :doc:`../prioritization` for detailed documentation.