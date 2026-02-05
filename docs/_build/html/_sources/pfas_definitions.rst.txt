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

**SMARTS Pattern**::

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

Positive matches::

   "C(F)(F)(F)C(F)(F)C(F)(F)F"          # Perfluoropropane
   "OC(=O)C(F)(F)C(F)(F)C(F)(F)F"       # PFBA
   "S(=O)(=O)(O)C(F)(F)C(F)(F)F"        # PFBS

Negative matches::

   "CCCCCC"         # No fluorination
   "C(F)CC"         # Only monofluoro
   "FC(F)C"         # CHF groups
   "C(F)(F)(Cl)C"   # Halogen substitution

2. EU PFAS Restriction (REACH)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Definition**: Similar to OECD but with additional exclusion criteria for hydroxyl,
alkoxy, amino, and other reactive functional groups in specific positions.

**Description**:
This definition extends the OECD framework by excluding fluorinated carbons directly bonded to
certain functional groups (hydroxyl, alkoxy, amino groups) to avoid capturing metabolites and
transformation products with reduced bioaccumulation potential.

**Key Features**:

- Similar core pattern to OECD
- Excludes carbons bonded to -OH, -OR, -NH₂, -NHR, -NR₂ groups
- More restrictive than OECD definition

3. OPPT 2023 (US EPA)
^^^^^^^^^^^^^^^^^^^^^

**Definition**: Molecules meeting at least one of three structural criteria related to
fully fluorinated carbon chains and ether linkages.

**SMARTS Patterns**::

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

4. UK PFAS Definition
^^^^^^^^^^^^^^^^^^^^^

**Definition**: Contains at least one CF₃ or CF₂-CF₂ group (not bonded to Cl, Br, I).

**SMARTS Patterns**::

   Pattern 1 (CF₃): F[#6H0X4!$([#6][#17,#35,#53])](F)F
   Pattern 2 (CF₂-CF₂): F[#6](F)[#6X4](F)(F)

**Description**:
The UK definition is similar to the OECD framework but allows CHF groups and focuses
on perfluoroalkyl and polyfluoroalkyl segments.

**Key Features**:

- Allows CHF groups (unlike OECD)
- Requires CF₃ or CF₂-CF₂ segment
- Excludes halogenated carbons (Cl, Br, I)

5. PFASTRUCTv5
^^^^^^^^^^^^^^

**Definition**: Contains specific structural patterns OR has a fluorine ratio ≥ 0.3.

**SMARTS Patterns**::

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

See Also
--------

- :doc:`quickstart` - Basic usage guide
- :doc:`tutorial` - Comprehensive tutorial
- :doc:`api/core` - API reference for parse functions

References
----------

1. OECD (2021). "Reconciling Terminology of the Universe of Per- and Polyfluoroalkyl Substances: Recommendations and Practical Guidance"
2. EU REACH Restriction Proposal (2023). "Perfluoroalkyl and Polyfluoroalkyl Substances (PFASs)"
3. US EPA OPPT (2023). "PFAS Structure Category"
4. UK Environment Agency (2022). "PFAS: Sources, Pathways and Environmental Data"
5. PFASTRUCTv5 (2020). Structural classification system for PFAS substances
