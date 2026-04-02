# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'HalogenGroups / PFASGroups'
copyright = '2026, Luc T. Miaz'
author = 'Luc T. Miaz'
release = '3.2.2'
version = '3.2.2'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'myst_parser',
    'sphinx_copybutton',
]

# Copy button: strip prompt characters from copy
copybutton_prompt_text = r'>>> |\.\.\. |\$ |\# |In \[\d+\]: | {2,5}\.\.\.: | {5,8}: '  # noqa: E501
copybutton_prompt_is_regexp = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'README.md', 'ResultsFingerprint_Guide.md', 'fingerprint_options.md']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = None
html_favicon = None

html_theme_options = {
    'logo_only': False,
    # 'display_version': True,  # Not supported in current RTD theme version
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': '#2980B9',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# -- Extension configuration -------------------------------------------------

# Napoleon settings for Google and NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'rdkit': ('https://www.rdkit.org/docs/', None),
}

# MyST-Parser configuration
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "html_image",
]

# Source file suffix
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# The master toctree document
master_doc = 'index'

# Autodoc configuration
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# ---------------------------------------------------------------------------
# Nitpick / cross-reference suppression
# ---------------------------------------------------------------------------
# Suppress warnings for common false-positive cross-references produced by
# Napoleon converting Google/NumPy-style docstring type-annotations to RST.
nitpick_ignore = [
    ('py:class', "default 'F'"),
    ('py:class', "default 'per'"),
    ('py:class', "default 'all'"),
    ('py:class', "default 'binary'"),
    ('py:class', "default None"),
    ('py:class', "default True"),
    ('py:class', "default False"),
    ('py:class', 'optional'),
    ('py:class', 'callable'),
    ('py:class', 'array'),
    ('py:class', 'column_names'),
    ('py:class', 'Chem.Mol'),
    ('py:class', 'rdkit.Chem.Mol'),
    ('py:class', 'optional [ref.class]'),
    ('py:class', 'SQLAlchemy Engine'),
    ('py:class', 'ResultsFingerprint'),
    ('py:class', 'GroupMatch'),
    ('py:class', 'MatchComponent'),
    ('py:class', 'PFASComponent'),
    ('py:class', 'PFASEmbedding'),
    ('py:class', 'PIL.Image.Image'),
    ('py:class', '-> both'),
    ('py:func', 'parse_smiles'),
    ('py:func', 'parse_mols'),
    ('py:func', 'generate_fingerprint'),
    ('py:func', 'get_HalogenGroups'),
    ('py:func', 'HalogenGroups.generate_fingerprint'),
    ('py:func', 'PFASGroups.parser.parse_smiles'),
    ('py:func', 'PFASGroups.parser.parse_mols'),
    ('py:func', 'PFASGroups.core.fragment_until_valence_is_correct'),
    ('py:func', 'rdkit.Chem.rdchem.Mol.GetSubstructMatches'),
    ('py:meth', 'HalogenGroups.ResultsModel.to_fingerprint'),
    ('py:meth', 'MoleculeResult.classify'),
    ('py:meth', 'PFASEmbedding.to_array'),
    ('py:meth', 'PFASEmbedding.compare_kld'),
    ('py:class', 'PFASGroups.fingerprints.PFASFingerprint'),
    ('py:class', 'default=None'),
]

# Suppress all "default XYZ", "optional", "callable", etc. via regex
nitpick_ignore_regex = [
    ('py:class', r"^default .*"),
    ('py:class', r"^optional$"),
    ('py:class', r"^callable$"),
    ('py:class', r"^PFASGroups\.PFASEmbeddings\.\w+$"),
    ('py:class', r"^sqlalchemy\..*"),
    ('py:obj',   r"^PFASEmbedding\..*"),
]
