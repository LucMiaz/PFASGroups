Contributing
============

Thank you for considering a contribution to HalogenGroups / PFASGroups!

.. contents:: Contents
   :local:
   :depth: 2

How to contribute
-----------------

1. **Open an issue** to discuss bugs, feature requests, or new group definitions
   before submitting a pull request.
2. **Fork** the repository on GitHub.
3. Create a **feature branch**: ``git checkout -b my-feature``
4. Make your changes and add tests.
5. Run the test suite (see below).
6. Open a **pull request** against the ``main`` branch.

Setting up a development environment
--------------------------------------

.. code-block:: bash

   git clone https://github.com/lucdecrem/PFASGroups.git
   cd PFASGroups
   conda create -n halogengroups-dev -c conda-forge python=3.10 rdkit networkx
   conda activate halogengroups-dev
   pip install -e ".[dev]"

Running tests
-------------

.. code-block:: bash

   pytest

To run only unit tests (no slow benchmarks):

.. code-block:: bash

   pytest -m "not slow"

To check coverage:

.. code-block:: bash

   pytest --cov=HalogenGroups --cov-report=term-missing

Code style
----------

The project uses `Black <https://black.readthedocs.io/>`_ for formatting and
`Flake8 <https://flake8.pycqa.org/>`_ for linting:

.. code-block:: bash

   black HalogenGroups/ tests/
   flake8 HalogenGroups/ tests/

Adding new halogen groups
--------------------------

New group definitions live in
``HalogenGroups/data/halogen_groups.json``.  Each entry must have:

.. code-block:: json

   {
     "group_id": 117,
     "name": "my_new_group",
     "category": "Generic",
     "smarts": "[CX4](F)(F)(F)[NX3]",
     "is_PFAS": true,
     "compute": true,
     "description": "Perfluoroalkyl amine"
   }

After adding a group:

1. Run ``pytest tests/test_groups.py`` to check that the new group loads and
   matches the expected test molecules.
2. Add at least one positive and one negative test molecule to
   ``tests/_test_examples.py``.

Reporting issues
----------------

Please include:

- Python and RDKit versions (``python --version``, ``python -c "import rdkit; print(rdkit.__version__)"``).
- The minimal SMILES that triggers the issue.
- Expected vs actual output.

Contact
-------

- **GitHub Issues**: https://github.com/lucdecrem/PFASGroups/issues
- **Email**: luc.miaz@aces.su.se