Installation
============

Requirements
------------

- Python 3.7 or newer
- `RDKit <https://www.rdkit.org/>`_ 2020.09 or newer
- `NetworkX <https://networkx.org/>`_ 2.5 or newer

Install with pip
----------------

The recommended installation method is via pip.  RDKit must already be
available in the target environment (it is not listed as a pip dependency
because it is best installed via conda or a pre-built wheel).

.. code-block:: bash

   pip install PFASgroups

Install with conda
------------------

If you manage environments with conda you can install both RDKit and the
package in one step:

.. code-block:: bash

   conda install -c conda-forge rdkit
   pip install PFASgroups

Or create a fresh environment:

.. code-block:: bash

   conda create -n halogengroups -c conda-forge python=3.10 rdkit networkx
   conda activate halogengroups
   pip install PFASgroups

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/lucdecrem/PFASGroups.git
   cd PFASGroups
   pip install -e .

Optional dependencies
---------------------

Some features require additional packages that are not installed by default:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Feature
     - Package
     - Install
   * - PCA / t-SNE
     - scikit-learn
     - ``pip install scikit-learn``
   * - UMAP
     - umap-learn
     - ``pip install umap-learn``
   * - Database export
     - SQLAlchemy
     - ``pip install sqlalchemy``
   * - PostgreSQL export
     - psycopg2
     - ``pip install psycopg2-binary``

Verify installation
-------------------

.. code-block:: python

   from HalogenGroups import parse_smiles, get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(f"Loaded {len(groups)} halogen groups")
   # Loaded 116 halogen groups

   results = parse_smiles(["CCCC(F)(F)F"])
   print(results[0].matches[0].group_name)
   # perfluoroalkyl

Troubleshooting
---------------

**ImportError: No module named rdkit**

Install RDKit via conda-forge (recommended) or use a wheel from
`https://www.lfd.uci.edu/~gohlke/pythonlibs/ <https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.

**ModuleNotFoundError: No module named 'HalogenGroups'**

Make sure you installed the package as ``PFASgroups`` (on PyPI the distribution
is called *PFASgroups*; both ``import HalogenGroups`` and ``import PFASGroups``
work after installation).

**scikit-learn / umap-learn not found**

These are optional.  They are only needed when calling
:meth:`~HalogenGroups.ResultsFingerprint.perform_pca`,
:meth:`~HalogenGroups.ResultsFingerprint.perform_tsne`, or
:meth:`~HalogenGroups.ResultsFingerprint.perform_umap`.