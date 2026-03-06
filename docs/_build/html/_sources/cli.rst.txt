Command-Line Interface
======================

PFASGroups ships a command-line tool available as both ``halogengroups`` and
``pfasgroups`` (aliases for the same entry point).

Installation makes the command available in your shell after:

.. code-block:: bash

   pip install PFASGroups
   # or: pip install -e .


Synopsis
--------

.. code-block:: text

   halogengroups <command> [options] [smiles ...]

   Commands:
     parse          Parse SMILES and identify halogen groups
     fingerprint    Generate group fingerprints
     list-groups    List all available halogen groups
     list-paths     List available component SMARTS types


parse
-----

Identify halogen groups in one or more SMILES strings.

.. code-block:: bash

   halogengroups parse [options] [smiles ...]

**Positional arguments**

.. code-block:: text

   smiles    One or more SMILES strings (quoted if they contain spaces)

**Options**

.. code-block:: text

   -i, --input FILE          Read SMILES from file (one per line)
   -o, --output FILE         Write results to FILE (default: stdout)
   --format {json,csv}       Output format (default: json)
   --pretty                  Pretty-print JSON
   --halogens H [H ...]      Filter components by halogen(s): F Cl Br I
   --saturation {per,poly}   Filter by saturation level
   --form {alkyl,cyclic}     Filter by molecular form
   --no-component-metrics    Skip all graph-theory metrics (fastest mode)
   --limit-effective-graph-resistance N
                             Only compute effective resistance for components
                             with fewer than N atoms (0 = disable; omit = always)
   --groups-file FILE        Custom halogen groups JSON file
   --component_smarts-file FILE
                             Custom component SMARTS JSON file

**Examples**

.. code-block:: bash

   # Single SMILES
   halogengroups parse "FC(F)(F)C(F)(F)C(=O)O"

   # Multiple SMILES, pretty-printed JSON
   halogengroups parse --pretty \
       "FC(F)(F)C(F)(F)C(=O)O" \
       "FC(F)(F)C(F)(F)S(=O)(=O)O"

   # Fluorine-only, perfluorinated alkyl components
   halogengroups parse --halogens F --saturation per --form alkyl \
       "FC(F)(F)C(F)(F)C(=O)O"

   # All four halogens
   halogengroups parse --halogens F Cl Br I \
       "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"

   # From file, CSV output
   halogengroups parse --input smiles.txt --output results.csv --format csv

   # Skip graph metrics (fastest)
   halogengroups parse --no-component-metrics "FC(F)(F)C(F)(F)C(=O)O"

   # Compute resistance only for small components
   halogengroups parse --limit-effective-graph-resistance 100 \
       "FC(F)(F)C(F)(F)C(=O)O"

   # Custom group file
   halogengroups parse --groups-file my_groups.json "FC(F)(F)C(F)(F)C(=O)O"

**JSON output structure**

.. code-block:: json

   [
     {
       "smiles": "FC(F)(F)C(F)(F)C(=O)O",
       "matches": [
         {
           "group_id": 1,
           "group_name": "Perfluoroalkyl carboxylic acids",
           "match_count": 1,
           "components": [
             {
               "size": 3,
               "SMARTS": "Perfluoroalkyl",
               "branching": 1.0,
               "smarts_centrality": 0.5
             }
           ]
         }
       ]
     }
   ]

**CSV output columns**

.. code-block:: text

   smiles, group_id, group_name, match_count, component_sizes


fingerprint
-----------

Generate halogen-group fingerprints suitable for machine learning.

.. code-block:: bash

   halogengroups fingerprint [options] [smiles ...]

**Positional arguments**

.. code-block:: text

   smiles    One or more SMILES strings

**Options**

.. code-block:: text

   -i, --input FILE               Read SMILES from file
   -o, --output FILE              Write output to FILE (default: stdout)
   -g, --groups SPEC              Group selection as range "1-28" or
                                  comma-separated "1,2,3" (default: all)
   -f, --format {vector,dict,sparse,detailed,int}
                                  Fingerprint representation (default: vector)
   --count-mode {binary,count,max_chain}
                                  Encoding mode (default: binary)
   --halogens H [H ...]           Halogens to include (default: F)
   --output-format {json,csv}     File format (default: json)
   --pretty                       Pretty-print JSON

**Examples**

.. code-block:: bash

   # Binary vector (default)
   halogengroups fingerprint "FC(F)(F)C(F)(F)C(=O)O"

   # Dictionary representation
   halogengroups fingerprint --format dict "FC(F)(F)C(F)(F)C(=O)O"

   # OECD groups only (IDs 1–28)
   halogengroups fingerprint --groups 1-28 "FC(F)(F)C(F)(F)C(=O)O"

   # Count mode
   halogengroups fingerprint --count-mode count "FC(F)(F)C(F)(F)C(=O)O"

   # Multi-halogen stacked fingerprint (116 × 4 = 464 columns)
   halogengroups fingerprint --halogens F Cl Br I "FC(F)(F)C(F)(F)C(=O)O"

   # From file, save to CSV
   halogengroups fingerprint --input smiles.txt \
       --output fps.csv --output-format csv


list-groups
-----------

List all available halogen groups and their definitions.

.. code-block:: bash

   halogengroups list-groups [options]

**Options**

.. code-block:: text

   -o, --output FILE   Write to FILE (default: stdout)
   --pretty            Pretty-print JSON (default: true)

**Example output (excerpt)**

.. code-block:: json

   [
     {
       "id": 1,
       "name": "Perfluoroalkyl carboxylic acids",
       "alias": "PFCA",
       "componentSmarts": "Perfluoroalkyl",
       "componentSaturation": "per",
       "componentHalogens": "F",
       "componentForm": "alkyl"
     },
     {
       "id": 6,
       "name": "Perfluoroalkyl sulfonic acids",
       "alias": "PFSA"
     }
   ]

**Examples**

.. code-block:: bash

   # Print to console
   halogengroups list-groups

   # Save to file
   halogengroups list-groups --output groups.json


list-paths
----------

List available component SMARTS types (path/component definitions).

.. code-block:: bash

   halogengroups list-paths [options]

**Options**

.. code-block:: text

   -o, --output FILE   Write to FILE (default: stdout)

**Example output (excerpt)**

.. code-block:: json

   {
     "Perfluoroalkyl": {
       "component": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
       "end": "[C;X4;H0](F)(F)F",
       "halogen": "F",
       "form": "alkyl",
       "saturation": "per"
     },
     "Polyfluoroalkyl": {
       "component": "[C;X4;H1](F)!@!=!#[C;X4](F)",
       "end": "[C;X4;H1](F)F",
       "halogen": "F",
       "form": "alkyl",
       "saturation": "poly"
     }
   }


Global Options
--------------

These options can be placed before any sub-command:

.. code-block:: text

   --groups-file FILE            Custom halogen groups JSON
   --component_smarts-file FILE  Custom component SMARTS JSON


Environment
-----------

The CLI uses the Python environment in which PFASGroups is installed.  To use a
specific environment:

.. code-block:: bash

   # conda
   conda activate chem
   halogengroups parse "FC(F)(F)C(F)(F)C(=O)O"

   # direct Python invocation
   python -m PFASGroups.cli parse "FC(F)(F)C(F)(F)C(=O)O"


Performance Tips
----------------

For large input files:

.. code-block:: bash

   # Skip graph metrics (5–10× faster for large molecules)
   halogengroups parse --no-component-metrics --input big_file.smi \
       --output results.json

   # Skip effective graph resistance for molecules with > 50 atoms
   halogengroups parse --limit-effective-graph-resistance 50 \
       --input big_file.smi --output results.json

See the :doc:`benchmarking` page for
timing data.


See Also
--------

* :doc:`quickstart` — Python API quick start
* :doc:`api/core` — ``parse_smiles``, ``generate_fingerprint`` reference
* :doc:`customization` — using custom groups and path files with the CLI
