Command Line Interface
======================

PFASgroups provides a convenient command-line interface for quick analysis.

Basic Usage
-----------

After installation, the ``pfasgroups`` command is available:

.. code-block:: bash

   pfasgroups --help

Parse Command
-------------

Identify PFAS groups in structures:

.. code-block:: bash

   # Parse SMILES from command line
   pfasgroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

   # Parse from file (one SMILES per line)
   pfasgroups parse --input smiles.txt --output results.json

   # Component-based analysis
   pfasgroups parse --bycomponent "C(C(F)(F)F)F"

   # Pretty-print JSON output
   pfasgroups parse "C(C(F)(F)F)F" --pretty

Input File Format
^^^^^^^^^^^^^^^^^

Create a text file with one SMILES per line:

.. code-block:: text

   C(C(F)(F)F)F
   FC(F)(F)C(F)(F)C(=O)O
   FC(F)C(F)(F)F

Output Format
^^^^^^^^^^^^^

The parse command outputs JSON with comprehensive metrics:

.. code-block:: json

   [
     {
       "smiles": "FC(F)(F)C(F)(F)C(=O)O",
       "inchikey": "...",
       "matches": [
         {
           "group_name": "Perfluoroalkyl carboxylic acids",
           "id": 1,
           "match_count": 1,
           "components_sizes": [3],
           "mean_branching": 1.0,
           "mean_component_fraction": 0.812,
           "total_components_fraction": 0.875,
           "mean_eccentricity": 2.5,
           "median_eccentricity": 2.5,
           "components": [
             {
               "component": [0, 1, 2, 3],
               "size": 4,
               "component_fraction": 0.812,
               "branching": 1.0,
               "mean_eccentricity": 2.5,
               "SMARTS": "alkyl"
             }
           ]
         }
       ]
     }
   ]

**Key Metrics:**

- ``mean_branching``: Average linearity (1.0=linear, 0.0=branched)
- ``mean_component_fraction``: Average molecular coverage per component
- ``total_components_fraction``: Total coverage by union of all components
- ``mean_eccentricity``: Average graph-theoretic eccentricity
- ``component_fraction``: Individual component coverage (includes all attached H, F, Cl, Br, I)

Fingerprint Command
-------------------

Generate PFAS fingerprints:

.. code-block:: bash

   # Basic fingerprint
   pfasgroups fingerprint "C(C(F)(F)F)F"

   # From file
   pfasgroups fingerprint --input smiles.txt --output fingerprints.json

   # Select specific groups (range)
   pfasgroups fingerprint "C(C(F)(F)F)F" --groups 28-52

   # Select specific groups (list)
   pfasgroups fingerprint "C(C(F)(F)F)F" --groups 28,29,30

   # Different output formats
   pfasgroups fingerprint "C(C(F)(F)F)F" --format dict
   pfasgroups fingerprint "C(C(F)(F)F)F" --format sparse
   pfasgroups fingerprint "C(C(F)(F)F)F" --format int

   # Count mode
   pfasgroups fingerprint "C(C(F)(F)F)F" --count-mode count
   pfasgroups fingerprint "C(C(F)(F)F)F" --count-mode max_chain

List Groups Command
-------------------

Display available PFAS groups:

.. code-block:: bash

   # List all default groups
   pfasgroups list-groups

   # List with details
   pfasgroups list-groups --verbose

Output example:

.. code-block:: text

   ID  Name                                  Alias
   1   Perfluoroalkyl carboxylic acids       PFCAs
   2   Polyfluoroalkyl carboxylic acid       PolyFCAs
   3   Perfluoroalkyl dicarboxylic acids     PFdiCAs
   ...

List Paths Command
------------------

Display available pathway types:

.. code-block:: bash

   pfasgroups list-paths

Output:

.. code-block:: text

   Available pathway types:
   - Perfluoroalkyl: Fully fluorinated carbon chains
   - Polyfluoroalkyl: Partially fluorinated carbon chains
   - Polyfluoro: Branched polyfluorinated chains
   - Polyfluorobr: Brominated polyfluorinated chains

Using Custom Configuration
--------------------------

Custom Groups File
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pfasgroups parse --groups-file custom_groups.json "C(C(F)(F)F)F"

Custom Paths File
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pfasgroups parse --fpaths-file custom_fpaths.json "C(C(F)(F)F)F"

Both Custom Files
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pfasgroups fingerprint \
       --fpaths-file custom_fpaths.json \
       --groups-file custom_groups.json \
       --input smiles.txt

Scripting Examples
------------------

Bash Script
^^^^^^^^^^^

.. code-block:: bash

   #!/bin/bash
   
   # Process multiple files
   for file in data/*.txt; do
       output="${file%.txt}_results.json"
       pfasgroups parse --input "$file" --output "$output"
   done

PowerShell Script
^^^^^^^^^^^^^^^^^

.. code-block:: powershell

   # Process multiple files
   Get-ChildItem data\*.txt | ForEach-Object {
       $output = $_.FullName -replace '\.txt$', '_results.json'
       pfasgroups parse --input $_.FullName --output $output
   }

Piping with Other Tools
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Pipe SMILES from another command
   cat smiles.txt | pfasgroups parse --input -

   # Process and filter with jq
   pfasgroups parse "FC(F)(F)C(F)(F)C(=O)O" | jq '.[] | .groups[] | .name'
