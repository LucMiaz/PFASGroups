PFAS Definitions
================

PFASGroups can classify molecules against five regulatory and scientific
PFAS definitions.  Pass ``include_PFAS_definitions=True`` to
:func:`~PFASGroups.parse_smiles` to enable classification:

.. code-block:: python

   from PFASGroups import parse_smiles

   smiles = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OCCOCCO"]

   results = parse_smiles(smiles, include_PFAS_definitions=True)

   for mol in results:
       names = [d.definition_name for d in mol.pfas_definition_matches]
       print(mol.smiles, "->", names if names else "(none)")

.. contents:: Supported definitions
   :local:
   :depth: 1

OECD 2021
---------

**Full name**: OECD Reconciled PFAS Definition (2021)

**Reference**: `OECD, 2021. Reconciling Terminology of the Universe of
Per- and Polyfluoroalkyl Substances. <https://www.oecd.org/en/publications/reconciling-terminology-of-the-universe-of-per-and-polyfluoroalkyl-substances-pfas_d2a7ea98-en.html>`_

**Rule**: A substance is a PFAS if it contains at least one fully
fluorinated methyl (CF₃–) or methylene (–CF₂–) group that is not part of
a –CF₂–O– ether linkage.

This definition is implemented through the 28 OECD halogen groups bundled
in the library.  A molecule matches if it has at least one group with
``group_category == 'OECD'``.

EU REACH
--------

**Reference**: European Chemicals Agency (ECHA) REACH guidance.

**Rule**: Contains at least one perfluoroalkyl or polyfluoroalkyl moiety with
a carbon chain length ≥ 4, and does not belong to any exempted polymer
category.

OPPT 2023
---------

**Full name**: US EPA Office of Pollution Prevention and Toxics PFAS
Definition (2023)

**Reference**: `US EPA PFAS Definition <https://www.epa.gov/pfas>`_

**Rule**: Aliphatic fluorinated compounds containing ≥ 1 C–F bond that
meets structural criteria for persistence or the ability to degrade to
persistent PFAS moieties.

UK Environment Agency
---------------------

**Full name**: UK Environment Agency Fluorinated Polymer Definition

**Reference**: UK EA Guidance on Fluorinated Polymers (2022).

**Rule**: Polymers or oligomers with a perfluorinated backbone or with
perfluorinated side-chains.

PFASSTRUCTv5
------------

**Reference**: `Schymanski et al., PFASSTRUCTV5 annotated structure database
<https://zenodo.org/record/7370805>`_

PFASSTRUCTv5 is a curated database of PFAS structures.  PFASGroups ships
an integrated classifier that reproduces the binary PFAS/non-PFAS label from
this database for molecules whose structural motifs overlap with the
library groups.

Retrieving definition objects
------------------------------

.. code-block:: python

   from PFASGroups import get_PFASDefinitions

   definitions = get_PFASDefinitions()
   for d in definitions:
       print(d.name, d.short_name)

Checking a single molecule
---------------------------

.. code-block:: python

   from PFASGroups import parse_smiles

   result = parse_smiles(["CCCC(F)(F)F"], include_PFAS_definitions=True)[0]
   for match in result.pfas_definition_matches:
       print(match.definition_name)
       print(match.rule_description)