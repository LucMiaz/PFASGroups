# PFASgroups Documentation

This directory contains the Sphinx documentation for PFASgroups.

## Building the Documentation

### Prerequisites

Install Sphinx and required extensions:

```bash
pip install -r requirements.txt
```

### Build HTML Documentation

```bash
make html
```

The generated documentation will be in `_build/html/`. Open `_build/html/index.html` in your browser.

### Build PDF Documentation

```bash
make latexpdf
```

Requires LaTeX installation (e.g., TeX Live, MiKTeX).

### Clean Build Files

```bash
make clean
```

## Documentation Structure

- `index.rst` - Main documentation entry point
- `installation.rst` - Installation instructions
- `quickstart.rst` - Quick start guide
- `tutorial.rst` - Comprehensive tutorial
- `algorithm.rst` - Algorithm explanation
- `pfas_groups.rst` - PFAS group classifications
- `api/` - API reference documentation
  - `core.rst` - Core functions
  - `models.rst` - Data models
  - `homologues.rst` - Homologue generation
  - `fragmentation.rst` - Molecular fragmentation
  - `visualization.rst` - Visualization functions
  - `generation.rst` - Molecule generation
- `customization.rst` - Customization guide
- `benchmarking.rst` - Benchmarking and validation
- `cli.rst` - Command-line interface
- `contributing.rst` - Contributing guidelines
- `changelog.rst` - Version history
- `license.rst` - License information

## Read the Docs

This documentation is configured for Read the Docs hosting.

### Configuration Files

- `conf.py` - Sphinx configuration
- `requirements.txt` - Python dependencies for building docs
- `.readthedocs.yaml` - Read the Docs configuration (in project root)

### Local Preview

To preview the documentation as it will appear on Read the Docs:

```bash
# Install dependencies
pip install -r requirements.txt

# Build HTML
make html

# Serve locally
python -m http.server 8000 --directory _build/html
```

Then open http://localhost:8000 in your browser.

## Contributing to Documentation

When adding new features to PFASgroups:

1. Update relevant `.rst` files
2. Add docstrings to new functions/classes (Google or NumPy style)
3. Add examples to the tutorial or quickstart
4. Update the API reference if needed
5. Add to changelog
6. Rebuild documentation locally to check for errors

### Docstring Style

Use Google or NumPy style docstrings:

```python
def my_function(param1, param2):
    """
    Brief description of function.
    
    Longer description with more details if needed.
    
    Args:
        param1 (str): Description of param1
        param2 (int): Description of param2
    
    Returns:
        bool: Description of return value
    
    Examples:
        >>> my_function("test", 42)
        True
    """
    pass
```

### Adding New Pages

1. Create new `.rst` file in appropriate directory
2. Add to `toctree` directive in `index.rst` or parent page
3. Use appropriate section markers: `===`, `---`, `~~~`, `^^^`
4. Include code examples and cross-references

### Cross-References

Reference other documentation pages:

```rst
See :doc:`quickstart` for examples
See :doc:`api/core` for API details
```

Reference Python objects:

```rst
:func:`parse_smiles`
:class:`PFASGroup`
:meth:`PFASGroup.formula_dict_satisfies_constraints`
```

## Maintenance

- Update version numbers in `conf.py`
- Keep examples synchronized with code changes
- Test all code examples
- Update screenshots and diagrams as UI changes
- Review and update external links periodically
