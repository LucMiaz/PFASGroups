Installation
============

PFASgroups can be installed in multiple ways depending on your needs and preferences.

From PyPI (Recommended)
-----------------------

The simplest way to install PFASgroups is via pip:

.. code-block:: bash

   pip install PFASgroups

This will install PFASgroups and all required dependencies.

From Conda-Forge
----------------

If you prefer conda package management:

.. code-block:: bash

   conda install -c conda-forge pfasgroups

From Source
-----------

For development or to get the latest features:

.. code-block:: bash

   git clone https://github.com/yourusername/PFASGroups.git
   cd PFASGroups
   pip install -e .

The ``-e`` flag installs in editable mode, allowing you to modify the source code.

Requirements
------------

PFASgroups requires Python 3.7 or later and the following packages:

Core Dependencies
~~~~~~~~~~~~~~~~~

- **RDKit** (>=2020.03): Chemistry toolkit for molecular operations
- **NumPy** (>=1.19.0): Numerical computing
- **Pandas** (>=1.1.0): Data manipulation and analysis
- **NetworkX** (>=2.5): Graph algorithms for pathfinding
- **tqdm**: Progress bars

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~

For visualization:

- **matplotlib** (>=3.3.0): Static plotting
- **svgutils**: SVG manipulation

For benchmarking:

- **seaborn**: Statistical visualization
- **plotly**: Interactive plots
- **scipy**: Statistical analysis

Verifying Installation
----------------------

To verify that PFASgroups is installed correctly:

.. code-block:: python

   import PFASgroups
   print(PFASgroups.__version__)

You should see the version number printed without errors.

Testing the Installation
------------------------

Run a quick test:

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Test with PFOA (Perfluorooctanoic acid)
   result = parse_smiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
   
   if result:
       print("✓ Installation successful!")
       print(f"Detected groups: {[g.name for g, _, _, _ in result[0]]}")
   else:
       print("✗ Installation issue detected")

Expected output should show detection of "Perfluoroalkyl carboxylic acids" (PFCAs).

Troubleshooting
---------------

RDKit Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~

RDKit can sometimes be challenging to install. If you encounter issues:

**Using conda (recommended for RDKit):**

.. code-block:: bash

   conda create -n pfas python=3.9
   conda activate pfas
   conda install -c conda-forge rdkit
   pip install PFASgroups

**On Windows:**

RDKit wheels are available for recent Python versions. If pip installation fails, use conda.

**On Linux:**

Most distributions have RDKit packages:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install python3-rdkit
   
   # Fedora
   sudo dnf install python3-rdkit

Import Errors
~~~~~~~~~~~~~

If you get import errors:

1. Ensure you're in the correct Python environment
2. Check that all dependencies are installed: ``pip list | grep -i rdkit``
3. Try importing RDKit directly: ``python -c "import rdkit; print(rdkit.__version__)"``

Performance Issues
~~~~~~~~~~~~~~~~~~

For better performance on large datasets:

- Ensure NumPy is linked to an optimized BLAS library (Intel MKL or OpenBLAS)
- Consider using PyPy for pure-Python operations (note: RDKit may not be available)
- Use parallel processing when analyzing multiple molecules

Upgrading
---------

To upgrade to the latest version:

.. code-block:: bash

   pip install --upgrade PFASgroups

To upgrade from conda:

.. code-block:: bash

   conda update -c conda-forge pfasgroups

Development Installation
------------------------

For contributors:

.. code-block:: bash

   git clone https://github.com/yourusername/PFASGroups.git
   cd PFASGroups
   pip install -e ".[dev]"

This installs additional development dependencies including:

- pytest: Testing framework
- black: Code formatter
- flake8: Linter
- sphinx: Documentation generator
