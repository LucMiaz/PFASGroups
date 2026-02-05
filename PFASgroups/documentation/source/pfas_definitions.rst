PFAS Definitions
================

PFASgroups supports multiple regulatory and academic PFAS definitions, enabling flexible classification
of per- and polyfluoroalkyl substances according to different frameworks.

Overview
--------

PFAS definitions vary across regulatory bodies and scientific communities. PFASgroups implements
five major PFAS definitions using structure-based SMARTS patterns:

1. **OECD Definition** (2021)
2. **EU PFAS Restriction** (REACH)
3. **OPPT 2023** (US EPA)
4. **UK PFAS Definition**
5. **PFASTRUCTv5** (structural patterns + fluorine ratio)

These definitions can be used to classify molecules and filter chemical inventories according
to specific regulatory requirements.

Supported Definitions
---------------------

1. OECD Definition (2021)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Contains at least one CF₃ or CF₂ group (not bonded to H, Cl, Br, or I).

**SMARTS Pattern**:

.. code-block:: text

   [#6X4!$([#6H1])!$([#6][#17,#35,#53])](F)F

**Description**: 
This definition captures fully or partially fluorinated alkyl moieties where carbon atoms
bear at least two fluorine atoms. The pattern explicitly excludes carbons bonded to hydrogen
(preventing inclusion of CHF groups) and halogens (Cl, Br, I).

**Key Features**:

- Requires CF₃ or CF₂ groups
- Excludes CHF groups (partially fluorinated)
- Excludes carbons bonded to Cl, Br, I

**Examples**:

Positive matches:

.. code-block:: python

   "C(F)(F)(F)C(F)(F)C(F)(F)F"          # Perfluoropropane
   "OC(=O)C(F)(F)C(F)(F)C(F)(F)F"       # PFBA
   "S(=O)(=O)(O)C(F)(F)C(F)(F)F"        # PFBS

Negative matches:

.. code-block:: python

   "CCCCCC"         # No fluorination
   "C(F)CC"         # Only monofluoro
   "FC(F)C"         # CHF groups
   "C(F)(F)(Cl)C"   # Halogen substitution

2. EU PFAS Restriction (REACH)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Similar to OECD but with additional exclusion criteria for hydroxyl,
alkoxy, amino, and other reactive functional groups in specific positions.

**SMARTS Patterns**:

.. code-block:: text

   Pattern 1 (CF₃):
   [#6X4!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])
         !$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),
              #7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])
                    [#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])
         !$([#6H1])!$([#6][#17,#35,#53])](F)(F)F

   Pattern 2 (CF₂):
   [#6X4!$([#6](F)(F)(F))
         !$([#6]([#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])
              [#6H3,#6H2X4,#6X3$([#6]=[#8]),a,#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),
               #16$([#16H1,#16$([#16][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]),
               #7$([#7H2,#7H1$([#7#6H3,#7#6H2X4,#7a,#6X3$([#6]=[#8])]),
                    #7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])])
         !$([#6]([#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),
                   #7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])
              [#6H3,#6H2X4,#6X3$([#6]=[#8]),a,#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),
               #16$([#16H1,#16$([#16][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]),
               #7$([#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),
                    #7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])])
         !$([#6H1])!$([#6][#17,#35,#53])](F)F

**Description**:
This definition extends the OECD framework by excluding fluorinated carbons directly bonded to
certain functional groups (hydroxyl, alkoxy, amino groups) to avoid capturing metabolites and
transformation products with reduced bioaccumulation potential.

**Key Features**:

- Similar core pattern to OECD
- Excludes carbons bonded to -OH, -OR, -NH₂, -NHR, -NR₂ groups
- More restrictive than OECD definition

**Examples**:

Positive matches (same as OECD positives)

Additional negatives:

.. code-block:: python

   "C(F)(F)(F)C(F)(F)OH"    # CF₂ bonded to OH
   "C(F)(F)(F)C(F)(F)OC"    # CF₂ bonded to OMe

3. OPPT 2023 (US EPA)
^^^^^^^^^^^^^^^^^^^^^

**Definition**: Molecules meeting at least one of three structural criteria related to
fully fluorinated carbon chains and ether linkages.

**SMARTS Patterns**:

.. code-block:: text

   Pattern 1: [#6H0X4](F)(F)[#6H0X4](F)
   Pattern 2: [#6](F)(F)([F,#8,#6H0X4])[#8][#6](F)(F)[F,#8,#6H0X4]
   Pattern 3: [#6](F)(F)(F)[#6]([#6](F)(F)F)([F,#6H0X4])[F,#6H0X4]

**Description**:
The OPPT definition uses three distinct structural motifs:

1. Adjacent fully fluorinated carbons (CF₂-CF or similar)
2. Perfluorinated ether linkages
3. Branched perfluorinated structures (e.g., perfluoroisopropyl groups)

**Key Features**:

- Three independent inclusion criteria
- Emphasizes fully fluorinated segments
- Captures perfluorinated ethers
- Includes branched perfluorinated moieties

**Examples**:

Positive matches:

.. code-block:: python

   "C(F)(F)(F)C(F)(F)C(F)(F)F"             # Perfluoropropane
   "C(F)(F)C(F)(F)C(F)(F)C"                # Mixed fluorination
   "C(F)(F)(F)C(F)(F)OC(F)(F)C(F)(F)F"     # Perfluoroether
   "C(F)(F)(C(F)(F)F)C(F)(F)C(F)(F)F"      # Branched perfluoro

Negative matches:

.. code-block:: python

   "CCCCCC"           # No fluorination
   "C(F)CCC"          # Insufficient fluorination
   "FC(F)C(F)F"       # CHF groups

4. UK PFAS Definition
^^^^^^^^^^^^^^^^^^^^^

**Definition**: Contains at least one CF₃ or CF₂-CF₂ group (not bonded to Cl, Br, I).

**SMARTS Patterns**:

.. code-block:: text

   Pattern 1 (CF₃): F[#6H0X4!$([#6][#17,#35,#53])](F)F
   Pattern 2 (CF₂-CF₂): F[#6](F)[#6X4](F)(F)

**Description**:
The UK definition is similar to the OECD framework but allows CHF groups and focuses
on perfluoroalkyl and polyfluoroalkyl segments.

**Key Features**:

- Allows CHF groups (unlike OECD)
- Requires CF₃ or CF₂-CF₂ segment
- Excludes halogenated carbons (Cl, Br, I)

**Examples**:

Positive matches:

.. code-block:: python

   "C(F)(F)(F)C(F)(F)C(F)(F)F"    # Perfluoropropane
   "C(F)(F)(F)C(F)(F)C"           # CF₃-CF₂
   "C(F)(F)C(F)(F)C"              # CF₂-CF₂
   "FC(F)(F)C(F)(F)C(F)(F)F"      # Same structures

Negative matches:

.. code-block:: python

   "CCCCCC"              # No fluorination
   "C(F)CC"              # Only monofluoro
   "C(F)(F)(Cl)C"        # Halogenated
   "C(F)(F)(Br)C"        # Halogenated

5. PFASTRUCTv5
^^^^^^^^^^^^^^

**Definition**: Contains specific structural patterns OR has a fluorine ratio ≥ 0.3.

**SMARTS Patterns**:

.. code-block:: text

   Pattern 1: F[#6](F)(F)[#6](F)
   Pattern 2: F[#6](F)[#6](F)(F)
   Pattern 3: F[#6](F)(F)[#6]~[#6](F)(F)
   Pattern 4: [#6](F)(F)[#5,#7,#8,#14,#15,#16][#6](F)(F)

**Fluorine Ratio**: ≥ 0.3 (number of fluorine atoms / total number of atoms)

**Description**:
PFASTRUCTv5 combines structural patterns with a quantitative fluorine content threshold,
providing a dual classification system.

**Key Features**:

- Four structural inclusion patterns
- Alternative fluorine ratio criterion (F/total atoms ≥ 0.3)
- Captures a broader range of fluorinated substances
- Includes heteroatom bridges (B, N, O, Si, P, S)

**Examples**:

Positive matches:

.. code-block:: python

   # Structural pattern matches
   "C(=O)(C(F)(F)Cl)C(F)(F)Cl"                   # Pattern 3
   "COC(C(OC(C(=O)C(F)(F)F)(F)F)F)(F)F"          # Pattern 2
   "C(C(CC(F)(F)F)(F)F)C(CC(F)(F)F)(F)F"         # Pattern 1
   
   # Fluorine ratio matches
   "C(C(=O)F)(F)F"                               # F ratio = 3/7 = 0.43
   "C(C(=O)CCSCCC(C(C(C(F)(F)F)(F)F)(F)F)(F)F"   # High F content

Negative matches:

.. code-block:: python

   "CCCCC(CC(C(F)(F)Br)(F)Cl)Br"    # Insufficient F, halogens
   "C(C(OC(Cl)(Cl)Cl)(F)F)(F)Cl"    # Pattern not matched
   "COC(=C(F)Cl)F"                   # Low F ratio, no pattern

Using PFAS Definitions in Code
-------------------------------

Checking PFAS Classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups import parse_smiles

   # Parse with PFAS definitions included
   results = parse_smiles(
       ["FC(F)(F)C(F)(F)C(=O)O"],
       include_PFAS_definitions=True
   )
   
   # Check which definitions the molecule satisfies
   for result in results:
       for group, count, chains, components in result:
           if group.id <= 5:  # PFAS definitions have IDs 1-5
               print(f"Matches {group.name}: {count > 0}")

Filtering by Definition
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups import parse_smiles
   import pandas as pd

   # Parse large dataset
   df = parse_smiles(
       smiles_list,
       output_format='dataframe',
       include_PFAS_definitions=True
   )
   
   # Filter molecules that match OECD definition (ID=1)
   oecd_pfas = df[df['group_id'] == 1]['smiles'].unique()
   
   # Filter molecules that match EU definition (ID=2)
   eu_pfas = df[df['group_id'] == 2]['smiles'].unique()
   
   # Count molecules matching each definition
   definition_counts = df[df['group_id'] <= 5].groupby('group_name')['smiles'].nunique()
   print(definition_counts)

Comparing Definitions
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups import parse_smiles
   import pandas as pd
   
   # Parse with all definitions
   df = parse_smiles(
       smiles_list,
       output_format='dataframe',
       include_PFAS_definitions=True
   )
   
   # Create a comparison matrix
   definitions = [
       ('OECD Definition', 1),
       ('EU PFAS Restriction', 2),
       ('OPPT 2023', 3),
       ('UK PFAS Definition', 4),
       ('PFASTRUCTv5', 5)
   ]
   
   comparison = {}
   for name, def_id in definitions:
       matches = set(df[df['group_id'] == def_id]['smiles'])
       comparison[name] = len(matches)
   
   print("PFAS counts by definition:")
   for name, count in comparison.items():
       print(f"  {name}: {count}")

Custom PFAS Definitions
-----------------------

You can extend PFASgroups with custom PFAS definitions by modifying the 
``PFAS_definitions_smarts.json`` file:

.. code-block:: json

   [
     {
       "id": 6,
       "name": "Custom PFAS Definition",
       "smarts": [
         "[your SMARTS pattern here]"
       ],
       "fluorineRatio": null,
       "description": "Description of your definition",
       "test": {
         "category": "definition",
         "examples": {
           "positives": ["SMILES1", "SMILES2"],
           "negatives": ["SMILES3", "SMILES4"]
         }
       }
     }
   ]

**Configuration File Location**:

.. code-block:: python

   import PFASgroups
   import os
   
   # Find the definitions file
   package_dir = os.path.dirname(PFASgroups.__file__)
   definitions_file = os.path.join(package_dir, 'data', 'PFAS_definitions_smarts.json')
   print(f"PFAS definitions: {definitions_file}")

Regulatory Implications
-----------------------

Understanding which PFAS definition is most appropriate depends on the regulatory context:

**OECD Definition**: Most widely adopted internationally, suitable for global assessments

**EU REACH**: Required for compliance with European Union chemical regulations

**OPPT 2023**: Relevant for US EPA reporting and regulatory compliance

**UK Definition**: Specific to UK chemical management frameworks

**PFASTRUCTv5**: Useful for comprehensive screening of fluorinated substances

**Recommendation**: When conducting PFAS inventories or assessments, analyze substances
against multiple definitions to understand the regulatory landscape across jurisdictions.

Performance Considerations
--------------------------

PFAS definition checking is computationally efficient but adds overhead to parsing:

- **Without definitions**: ~10-50 ms per molecule
- **With definitions**: ~15-70 ms per molecule (5-40% overhead)

For large datasets (>10,000 molecules), consider:

1. Pre-filtering obvious non-PFAS (molecules without fluorine)
2. Parallel processing using multiprocessing
3. Caching results for repeated analyses

.. code-block:: python

   # Pre-filter molecules without fluorine
   import pandas as pd
   
   df = pd.DataFrame({'smiles': smiles_list})
   df['has_fluorine'] = df['smiles'].str.contains('F')
   
   # Only check molecules with fluorine
   fluorinated = df[df['has_fluorine']]['smiles'].tolist()
   results = parse_smiles(fluorinated, include_PFAS_definitions=True)

See Also
--------

- :doc:`getting_started` - Basic usage guide
- :doc:`user_guide` - Comprehensive user documentation
- :doc:`api/core` - API reference for parse functions
- :doc:`pfas_groups/oecd_groups` - OECD PFAS group classifications

References
----------

1. OECD (2021). "Reconciling Terminology of the Universe of Per- and Polyfluoroalkyl Substances: Recommendations and Practical Guidance"
2. EU REACH Restriction Proposal (2023). "Perfluoroalkyl and Polyfluoroalkyl Substances (PFASs)"
3. US EPA OPPT (2023). "PFAS Structure Category"
4. UK Environment Agency (2022). "PFAS: Sources, Pathways and Environmental Data"
5. PFASTRUCTv5 (2020). Structural classification system for PFAS substances
