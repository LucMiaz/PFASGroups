Contributing
============

We welcome contributions to PFASgroups! This guide explains how to contribute
to the project.

Getting Started
---------------

Development Setup
^^^^^^^^^^^^^^^^^

1. **Fork the repository** on GitHub

2. **Clone your fork:**

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/PFASgroups.git
      cd PFASgroups

3. **Create a virtual environment:**

   .. code-block:: bash

      conda create -n pfasgroups-dev python=3.10
      conda activate pfasgroups-dev

4. **Install development dependencies:**

   .. code-block:: bash

      pip install -e ".[dev]"

5. **Install RDKit:**

   .. code-block:: bash

      conda install -c conda-forge rdkit

6. **Run tests to verify setup:**

   .. code-block:: bash

      pytest tests/

Code Style
----------

We follow PEP 8 with some modifications:

- Line length: 100 characters
- Use double quotes for strings
- Use type hints for function signatures

**Formatting:**

.. code-block:: bash

   # Format with black
   black PFASgroups/

   # Check with flake8
   flake8 PFASgroups/

**Example code style:**

.. code-block:: python

   from typing import List, Optional, Tuple
   from rdkit import Chem

   def parse_molecule(
       smiles: str,
       add_hydrogens: bool = True,
       sanitize: bool = True
   ) -> Optional[Chem.Mol]:
       """
       Parse a SMILES string into an RDKit molecule.

       Parameters
       ----------
       smiles : str
           The SMILES string to parse.
       add_hydrogens : bool, optional
           Whether to add explicit hydrogens. Default is True.
       sanitize : bool, optional
           Whether to sanitize the molecule. Default is True.

       Returns
       -------
       Optional[Chem.Mol]
           The parsed molecule, or None if parsing failed.

       Examples
       --------
       >>> mol = parse_molecule("CCO")
       >>> mol.GetNumAtoms()
       3
       """
       mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
       if mol is None:
           return None
       
       if add_hydrogens:
           mol = Chem.AddHs(mol)
       
       return mol

Testing
-------

Writing Tests
^^^^^^^^^^^^^

All new features should include tests:

.. code-block:: python

   # tests/test_new_feature.py
   import pytest
   from PFASgroups import new_function

   class TestNewFeature:
       """Tests for the new feature."""

       def test_basic_functionality(self):
           """Test basic usage of new_function."""
           result = new_function("input")
           assert result == "expected_output"

       def test_edge_cases(self):
           """Test edge cases."""
           assert new_function("") is None
           assert new_function(None) is None

       def test_error_handling(self):
           """Test error handling."""
           with pytest.raises(ValueError):
               new_function("invalid_input")

Running Tests
^^^^^^^^^^^^^

.. code-block:: bash

   # Run all tests
   pytest tests/

   # Run with coverage
   pytest --cov=PFASgroups --cov-report=html

   # Run specific test
   pytest tests/test_core.py::test_pfca_detection

   # Run with verbose output
   pytest -v tests/

Documentation
-------------

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cd docs/
   make html

   # View in browser
   open _build/html/index.html

Documentation Style
^^^^^^^^^^^^^^^^^^^

We use NumPy-style docstrings:

.. code-block:: python

   def example_function(param1: str, param2: int = 10) -> List[str]:
       """
       Short description of the function.

       Longer description explaining the purpose and behavior
       of the function in more detail.

       Parameters
       ----------
       param1 : str
           Description of param1.
       param2 : int, optional
           Description of param2. Default is 10.

       Returns
       -------
       List[str]
           Description of return value.

       Raises
       ------
       ValueError
           If param1 is empty.

       Examples
       --------
       >>> result = example_function("test", 5)
       >>> len(result)
       5

       See Also
       --------
       related_function : Does something related.

       Notes
       -----
       Additional implementation notes.

       References
       ----------
       .. [1] Citation or reference.
       """
       pass

Pull Requests
-------------

Submitting a Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Create a feature branch:**

   .. code-block:: bash

      git checkout -b feature/my-new-feature

2. **Make your changes** and commit:

   .. code-block:: bash

      git add .
      git commit -m "Add new feature: description"

3. **Push to your fork:**

   .. code-block:: bash

      git push origin feature/my-new-feature

4. **Open a Pull Request** on GitHub

Pull Request Checklist
^^^^^^^^^^^^^^^^^^^^^^

Before submitting, ensure:

- [ ] All tests pass locally
- [ ] New code has test coverage
- [ ] Documentation is updated
- [ ] Code follows style guidelines
- [ ] Commit messages are clear
- [ ] Branch is up to date with main

Adding New PFAS Groups
----------------------

To add a new PFAS group:

1. **Edit the JSON definition file:**

   .. code-block:: json

      {
          "id": 58,
          "name": "NewGroup",
          "smarts1": "[C](=O)O",
          "smarts2": null,
          "smartsPath": "Perfluoroalkyl",
          "constraints": {
              "eq": {"O": 2},
              "only": ["C", "F", "H", "O"]
          },
          "max_dist_from_CF": 2
      }

2. **Add test cases:**

   .. code-block:: python

      def test_new_group():
          results = parse_smiles(["FC(F)(F)C(F)(F)C(=O)O"])
          groups = [g.name for g, _, _, _ in results[0]]
          assert "NewGroup" in groups

3. **Update documentation** with group description

Reporting Issues
----------------

Bug Reports
^^^^^^^^^^^

Include:

- PFASgroups version
- Python version
- RDKit version
- Minimal reproducible example
- Expected vs actual behavior

Feature Requests
^^^^^^^^^^^^^^^^

Describe:

- The problem you're trying to solve
- Your proposed solution
- Alternative approaches considered

Code of Conduct
---------------

We follow the `Contributor Covenant <https://www.contributor-covenant.org/>`_
code of conduct. Please be respectful and inclusive in all interactions.

Contact
-------

- **Issues**: GitHub Issues
- **Email**: [maintainer email]
- **Discussions**: GitHub Discussions

License
-------

By contributing, you agree that your contributions will be licensed under
the MIT License.
