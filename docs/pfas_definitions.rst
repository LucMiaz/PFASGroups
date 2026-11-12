PFAS Definitions
================

PFASGroups classifies molecules against six regulatory and scientific PFAS
definitions.  The classification runs by default during parsing
(``include_PFAS_definitions=True``) for fluorine (``'F'`` in ``halogens``);
the definition matches are stored in ``mol.matches`` next to the group
matches, with ``type == 'PFASdefinition'``:

.. code-block:: python

   from PFASGroups import parse_smiles

   smiles = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OC(=O)c1ccc(C(F)(F)F)cc1", "OCCOCCO"]

   results = parse_smiles(smiles)

   for mol in results:
       names = [m["definition_name"] for m in mol.matches if m.is_definition]
       print(mol.smiles, "->", names if names else "(none)")

.. code-block:: text

   CCCC(F)(F)F -> ['OECD Definition', 'EU PFAS Restriction', 'UK PFAS Definition', 'PFASSTRUCTv5']
   O=C(O)C(F)(F)F -> ['OECD Definition', 'EU PFAS Restriction', 'UK PFAS Definition', 'PFASSTRUCTv5']
   O=C(O)c1ccc(C(F)(F)F)cc1 -> ['OECD Definition', 'EU PFAS Restriction', 'UK PFAS Definition']
   OCCOCCO -> (none)

Pass ``include_PFAS_definitions=False`` to skip the classification.

.. contents:: Supported definitions
   :local:
   :depth: 1

How a definition is evaluated
-----------------------------

Each definition (:class:`~PFASGroups.PFASDefinition`) has one or more SMARTS
patterns and, optionally, a minimum fluorine ratio.  A molecule matches when
any SMARTS pattern is found **or** the fluorine ratio reaches the threshold
(``requireBoth=True`` requires both).  The fluorine ratio is computed either

* by mass (``ratioBasis: "mass"``, the default): F mass / total mass, or
* by atom count (``ratioBasis: "count"``): F atoms / total atoms,

with hydrogen included in the total unless ``includeHydrogen`` is ``False``.
The two bases differ because fluorine (19.0) is heavier than C, N and O:
4-(trifluoromethyl)benzoic acid, C\ :sub:`8`\ H\ :sub:`5`\ F\ :sub:`3`\ O\ :sub:`2`,
is at 0.31 by mass but 0.23 by count without hydrogen.

Each definition match records what was tested:

.. code-block:: python

   mol = parse_smiles("OC(=O)c1ccc(C(F)(F)F)cc1")[0]
   d = next(m for m in mol.matches if m.is_definition)
   print(dict(d))
   # {'match_id': 'D1', 'id': 1, 'definition_name': 'OECD Definition',
   #  'type': 'PFASdefinition', 'halogen': 'F', 'smarts_match': True,
   #  'fluorine_ratio': None, 'fluorine_ratio_threshold': None,
   #  'fluorine_ratio_basis': 'mass', 'ratio_match': True}

OECD Definition (id 1)
----------------------

**Reference**: `OECD, 2021. Reconciling Terminology of the Universe of
Per- and Polyfluoroalkyl Substances. <https://www.oecd.org/en/publications/reconciling-terminology-of-the-universe-of-per-and-polyfluoroalkyl-substances-pfas_d2a7ea98-en.html>`_

**Rule**: at least one fully fluorinated methyl (CF\ :sub:`3`-) or methylene
(-CF\ :sub:`2`-) carbon, without any H, Cl, Br or I attached to it.

EU PFAS Restriction (id 2)
--------------------------

**Reference**: ECHA, proposal for a universal PFAS restriction (2023).

**Rule**: as the OECD definition, with the additional exclusions of the
restriction proposal: CF\ :sub:`3`/CF\ :sub:`2` groups bonded to the
non-fluorinated -OH, -O-R, -NH\ :sub:`2`, -NHR, -NR\ :sub:`2` or -SH/-SR
moieties listed in the proposal do not count.

OPPT 2023 (id 3)
----------------

**Reference**: US EPA, TSCA Section 8(a)(7) PFAS reporting rule,
40 CFR 705.3 (2023).

**Rule**: any of the three structural criteria of the final rule:

* R-(CF\ :sub:`2`)-CF(R')R'', both carbons saturated;
* R-CF\ :sub:`2`\ OCF\ :sub:`2`-R', both CF\ :sub:`2` saturated;
* CF\ :sub:`3`\ C(CF\ :sub:`3`)R'R''.

R groups are unconstrained (they may be H), as in the final 2023 rule.
Validated against EPA's PFAS8a7v3 CompTox list: 99.96 % coverage
(13 035 / 13 040).

UK PFAS Definition (id 4)
-------------------------

**Rule**: at least one CF\ :sub:`3` group, or a -CF\ :sub:`2`\ -CF\ :sub:`2`-
unit, on saturated carbons without H, Cl, Br or I attached.

PFASSTRUCTv5 (id 5)
-------------------

**Reference**: `PFASSTRUCTv5 annotated structure list
<https://zenodo.org/record/7370805>`_; fluorine-ratio criterion from Gaines
et al. (2023), `doi:10.1002/ieam.4735 <https://doi.org/10.1002/ieam.4735>`_.

**Rule**: any of the PFASSTRUCTv5 substructures (-CF\ :sub:`2`\ -CF-,
CF\ :sub:`3`\ -C~CF\ :sub:`2`, CF\ :sub:`2`-X-CF\ :sub:`2` with X = B, N, O,
Si, P or S) **or** at least 30 % of the non-hydrogen atoms are fluorine
(``ratioBasis: "count"``, ``includeHydrogen: false``).  Gaines et al. define
the threshold "strictly on atom counting", so since version 3.4.6 the ratio is
computed by count, not by mass.

OPPT 2021 (id 6)
----------------

**Reference**: US EPA, proposed TSCA Section 8(a)(7) rule (2021).

**Rule**: R-(CF\ :sub:`2`)-CF(R')R'', where the CF\ :sub:`2` and CF carbons
are saturated and none of the R groups can be hydrogen.  Superseded by
OPPT 2023, which drops the hydrogen restriction; kept for comparison.

Retrieving definitions
----------------------

:func:`~PFASGroups.get_PFASDefinitions` returns the raw definitions as
dictionaries (``id``, ``name``, ``smarts``, ``fluorineRatio``,
``description``, ...):

.. code-block:: python

   from PFASGroups import get_PFASDefinitions

   for d in get_PFASDefinitions():
       print(d["id"], d["name"], "-", d["description"][:60])

Checking a single molecule
--------------------------

:class:`~PFASGroups.PFASDefinition` can be used on its own;
``match_details`` reports which criterion was met, and ``ratio_basis``
overrides the definition's basis for one call:

.. code-block:: python

   from PFASGroups import PFASDefinition, get_PFASDefinitions

   pfasstruct = PFASDefinition(**get_PFASDefinitions()[4])   # PFASSTRUCTv5
   smi = "OC(=O)c1ccc(C(F)(F)F)cc1"

   print(pfasstruct.applies_to_molecule(smi))                     # False
   print(pfasstruct.match_details(smi)["fluorine_ratio"])         # 0.2308 (count)
   print(pfasstruct.match_details(smi, ratio_basis="mass")["fluorine_ratio"])  # 0.3079
