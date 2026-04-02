Installation
============

Requirements
------------

- Python 3.7 or newer
- `RDKit <https://www.rdkit.org/>`_ 2020.09 or newer
- `NetworkX <https://networkx.org/>`_ 2.5 or newer

Install with pip
----------------

The recommended installation method is via pip.  RDKit must already be available in the target environment. It is recommended to use an environment manager (like Conda/Mamba, e.g [Miniforge](https://github.com/conda-forge/miniforge)) and install RDKit via 

.. code-block:: bash

   mamba install -y -c rdkit rdkit

Then install PFASGroups from PyPI:

.. code-block:: bash

   pip install PFASgroups



Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/LucMiaz/PFASGroups.git
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
:meth:`~PFASGroups.PFASEmbeddingSet.perform_pca`,
:meth:`~PFASGroups.PFASEmbeddingSet.perform_tsne`, or
:meth:`~PFASGroups.PFASEmbeddingSet.perform_umap`.