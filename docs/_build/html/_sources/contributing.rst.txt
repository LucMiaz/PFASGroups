Contributing
============

We welcome contributions to PFASgroups! This guide will help you get started.

Getting Started
---------------

Setting Up Development Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Fork the repository** on GitHub

2. **Clone your fork:**

   .. code-block:: bash

      git clone https://github.com/yourusername/PFASGroups.git
      cd PFASGroups

3. **Create a development environment:**

   .. code-block:: bash

      conda create -n pfasgroups-dev python=3.9
      conda activate pfasgroups-dev
      conda install -c conda-forge rdkit

4. **Install in development mode:**

   .. code-block:: bash

      pip install -e ".[dev]"

5. **Install pre-commit hooks:**

   .. code-block:: bash

      pre-commit install

Development Workflow
~~~~~~~~~~~~~~~~~~~~

1. **Create a feature branch:**

   .. code-block:: bash

      git checkout -b feature/my-new-feature

2. **Make your changes**

3. **Run tests:**

   .. code-block:: bash

      pytest tests/

4. **Check code style:**

   .. code-block:: bash

      black PFASgroups/
      flake8 PFASgroups/

5. **Commit your changes:**

   .. code-block:: bash

      git add .
      git commit -m "Add my new feature"

6. **Push to your fork:**

   .. code-block:: bash

      git push origin feature/my-new-feature

7. **Create a Pull Request** on GitHub

Types of Contributions
----------------------

Bug Reports
~~~~~~~~~~~

When reporting bugs, please include:

- **Description** of the bug
- **Steps to reproduce**
- **Expected behavior**
- **Actual behavior**
- **System information** (OS, Python version, RDKit version)
- **SMILES string** that causes the issue (if applicable)

**Example bug report:**

.. code-block:: text

   **Bug**: Incorrect chain length for branched PFCA
   
   **Environment:**
   - OS: Ubuntu 22.04
   - Python: 3.9.7
   - RDKit: 2023.03.1
   - PFASgroups: 1.2.2
   
   **SMILES:** FC(F)(F)C(F)(F)C(F)(C(F)(F)F)C(=O)O
   
   **Expected:** Chain length should be 4
   **Actual:** Chain length reported as 5
   
   **Steps to reproduce:**
   ```python
   from PFASgroups import parse_smiles
   results = parse_smiles("FC(F)(F)C(F)(F)C(F)(C(F)(F)F)C(=O)O")
   print(results[0][0][2])  # Chain lengths
   ```

Feature Requests
~~~~~~~~~~~~~~~~

For feature requests, please describe:

- **What** you want to achieve
- **Why** it would be useful
- **How** you envision it working (if you have ideas)
- **Examples** of use cases

Documentation
~~~~~~~~~~~~~

Documentation improvements are always welcome:

- Fix typos or clarify existing docs
- Add examples
- Improve API documentation
- Translate documentation (future)

**Documentation style:**

- Use reStructuredText (.rst) format
- Include code examples
- Cross-reference related sections
- Test code examples

Code Contributions
------------------

Coding Standards
~~~~~~~~~~~~~~~~

We follow PEP 8 with some modifications:

- **Line length**: 100 characters (not 79)
- **Formatting**: Use Black formatter
- **Imports**: Organize with isort
- **Docstrings**: Google or NumPy style

**Example:**

.. code-block:: python

   from typing import List, Tuple, Optional
   from rdkit import Chem
   
   def my_function(
       smiles: str,
       option: Optional[bool] = None
   ) -> Tuple[int, List[str]]:
       """
       Brief description of function.
       
       Longer description with more details if needed.
       
       Args:
           smiles: SMILES string of molecule
           option: Optional parameter description
       
       Returns:
           Tuple of (count, list of names)
       
       Raises:
           ValueError: If SMILES is invalid
       
       Examples:
           >>> my_function("CCO")
           (1, ["ethanol"])
       """
       # Implementation
       pass

Testing
~~~~~~~

All new features must include tests:

.. code-block:: python

   # tests/test_my_feature.py
   import pytest
   from PFASgroups import parse_smiles
   
   def test_my_feature():
       """Test description."""
       smiles = "FC(F)(F)C(F)(F)C(=O)O"
       results = parse_smiles(smiles)
       
       assert len(results[0]) > 0
       assert results[0][0][0].id == 1  # PFCA group

**Run tests:**

.. code-block:: bash

   # All tests
   pytest
   
   # Specific test file
   pytest tests/test_my_feature.py
   
   # With coverage
   pytest --cov=PFASgroups tests/

Adding New PFAS Groups
~~~~~~~~~~~~~~~~~~~~~~

To add a new PFAS group to the default set:

1. **Add to JSON file** (``PFASgroups/data/PFAS_groups_smarts.json``):

   .. code-block:: json

      {
        "id": 56,
        "name": "Your New Group",
        "alias": "YNG",
        "smarts1": "[SMARTS_PATTERN]",
        "smartsPath": "Perfluoroalkyl",
        "constraints": {
          "gte": {"F": 1}
        }
      }

2. **Add tests** (``tests/test_pfas_groups.py``):

   .. code-block:: python

      def test_new_group():
          """Test the new PFAS group."""
          smiles = "FC(F)(F)C(F)(F)C(=O)NEWGROUP"
          results = parse_smiles(smiles)
          assert 56 in [g.id for g, _, _, _ in results[0]]

3. **Update documentation** (``docs/pfas_groups.rst``)

4. **Add benchmark test cases** (``benchmark/data/test_compounds.csv``)

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

When optimizing:

1. **Profile first:**

   .. code-block:: python

      import cProfile
      import pstats
      
      profiler = cProfile.Profile()
      profiler.enable()
      
      # Your code here
      parse_smiles(smiles_list)
      
      profiler.disable()
      stats = pstats.Stats(profiler)
      stats.sort_stats('cumulative')
      stats.print_stats(20)

2. **Benchmark before and after:**

   .. code-block:: bash

      cd benchmark
      python scripts/timing_benchmark.py

3. **Document performance impact** in PR

Pull Request Process
--------------------

PR Checklist
~~~~~~~~~~~~

Before submitting:

- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] New tests added for new features
- [ ] Documentation updated
- [ ] Changelog updated (``CHANGELOG.md``)
- [ ] No merge conflicts

PR Template
~~~~~~~~~~~

.. code-block:: markdown

   ## Description
   Brief description of changes
   
   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Documentation
   - [ ] Performance improvement
   
   ## Changes Made
   - Change 1
   - Change 2
   
   ## Testing
   How was this tested?
   
   ## Checklist
   - [ ] Tests pass
   - [ ] Documentation updated
   - [ ] Code follows style guide

Review Process
~~~~~~~~~~~~~~

1. **Automated checks** run on PR
2. **Maintainer reviews** code
3. **Discussion** if needed
4. **Approval** and merge

Communication
-------------

Channels
~~~~~~~~

- **GitHub Issues**: Bug reports, feature requests
- **GitHub Discussions**: Questions, ideas, general discussion
- **Email**: luc.miaz@aces.su.se (for private matters)

Communication Guidelines
~~~~~~~~~~~~~~~~~~~~~~~~

- Be respectful and constructive
- Provide context and examples
- Search existing issues before creating new ones
- Use clear, descriptive titles
- Follow up on your issues/PRs

License
-------

By contributing, you agree that your contributions will be licensed under the CC BY-NC 4.0 License.

For commercial use, contributors should be aware that commercial licensing is handled separately.

Code of Conduct
---------------

Pledge
~~~~~~

We pledge to make participation in our project a harassment-free experience for everyone.

Standards
~~~~~~~~~

**Positive behavior:**

- Using welcoming language
- Being respectful of differing viewpoints
- Gracefully accepting constructive criticism
- Focusing on what is best for the community

**Unacceptable behavior:**

- Harassment or discriminatory language
- Trolling or insulting comments
- Publishing others' private information
- Other conduct considered inappropriate

Enforcement
~~~~~~~~~~~

Report violations to: luc.miaz@aces.su.se

Recognition
-----------

Contributors will be:

- Listed in ``CONTRIBUTORS.md``
- Acknowledged in release notes
- Credited in publications (for significant contributions)

Questions?
----------

Don't hesitate to ask questions:

- Open a GitHub Discussion
- Email the maintainers
- Check existing documentation and issues

We appreciate your interest in contributing to PFASgroups!
