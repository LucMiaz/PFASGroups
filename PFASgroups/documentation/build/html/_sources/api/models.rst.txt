Data Models
===========

This section documents the data model classes used in PFASgroups.

.. module:: PFASgroups.PFASGroupModel
   :synopsis: PFAS group data model

PFASGroup Class
---------------

.. py:class:: PFASGroup(id, name, smarts1, smarts2, smartsPath, constraints, **kwargs)

   Represents a PFAS group definition with SMARTS patterns and constraints.

   :param id: Unique identifier for the group
   :type id: int
   :param name: Descriptive name of the group
   :type name: str
   :param smarts1: Primary SMARTS pattern for functional group detection
   :type smarts1: str or None
   :param smarts2: Secondary SMARTS pattern (optional, for path finding)
   :type smarts2: str or None
   :param smartsPath: Pathway type ('Perfluoroalkyl', 'Polyfluoroalkyl', 'cyclic', or None)
   :type smartsPath: str or None
   :param constraints: Molecular formula constraints dictionary
   :type constraints: dict
   :param max_dist_from_CF: Maximum bond distance from fluorinated carbon terminal (default: 2)
   :type max_dist_from_CF: int, optional

   **Attributes:**

   .. py:attribute:: id
      :type: int

      Unique identifier for the PFAS group.

   .. py:attribute:: name
      :type: str

      Descriptive name of the PFAS group.

   .. py:attribute:: smarts1
      :type: rdkit.Chem.Mol or None

      Preprocessed primary SMARTS pattern (RDKit Mol object).

   .. py:attribute:: smarts2
      :type: rdkit.Chem.Mol or None

      Preprocessed secondary SMARTS pattern (RDKit Mol object).

   .. py:attribute:: smartsPath
      :type: str or None

      Pathway type for chain finding.

   .. py:attribute:: constraints
      :type: dict

      Dictionary of molecular formula constraints.

   .. py:attribute:: max_dist_from_CF
      :type: int

      Maximum bond distance from fluorinated carbon for groups without formula constraints.

   **Methods:**

   .. py:method:: __str__()

      Returns the group name as string representation.

   .. py:method:: formula_dict_satisfies_constraints(formula_dict)

      Check if a molecular formula satisfies all constraints.

      :param formula_dict: Dictionary mapping element symbols to counts
      :type formula_dict: dict[str, int]
      :returns: True if all constraints are satisfied
      :rtype: bool

   .. py:method:: constraint_gte(formula_dict)

      Check greater-than-or-equal constraints.

      :param formula_dict: Element count dictionary
      :type formula_dict: dict
      :returns: True if all gte constraints satisfied
      :rtype: bool

   .. py:method:: constraint_lte(formula_dict)

      Check less-than-or-equal constraints.

      :param formula_dict: Element count dictionary
      :type formula_dict: dict
      :returns: True if all lte constraints satisfied
      :rtype: bool

   .. py:method:: constraint_eq(formula_dict)

      Check exact equality constraints.

      :param formula_dict: Element count dictionary
      :type formula_dict: dict
      :returns: True if all eq constraints satisfied
      :rtype: bool

   .. py:method:: constraint_only(formula_dict)

      Check element restriction constraints.

      :param formula_dict: Element count dictionary
      :type formula_dict: dict
      :returns: True if only allowed elements present
      :rtype: bool

   .. py:method:: constraint_rel(formula_dict)

      Check relative ratio constraints.

      :param formula_dict: Element count dictionary
      :type formula_dict: dict
      :returns: True if all relative constraints satisfied
      :rtype: bool

   **Example:**

   .. code-block:: python

      from PFASgroups import PFASGroup

      # Create a custom PFAS group
      group = PFASGroup(
          id=999,
          name="My Custom Group",
          smarts1="[C](F)(F)F",
          smarts2="C(=O)O",
          smartsPath="Perfluoroalkyl",
          constraints={
              "eq": {"O": 2},
              "gte": {"F": 3},
              "only": ["C", "F", "O", "H"]
          },
          max_dist_from_CF=2
      )

      # Check if a formula satisfies constraints
      formula = {"C": 4, "F": 7, "O": 2, "H": 1}
      is_valid = group.formula_dict_satisfies_constraints(formula)

Constraint Types
^^^^^^^^^^^^^^^^

The ``constraints`` dictionary supports the following types:

**eq (Exact Count)**

Requires exact element count:

.. code-block:: python

   {"eq": {"O": 2, "S": 1}}  # Exactly 2 O and 1 S

**gte (Greater Than or Equal)**

Requires minimum element count:

.. code-block:: python

   {"gte": {"F": 3, "O": 1}}  # At least 3 F and 1 O

**lte (Less Than or Equal)**

Requires maximum element count:

.. code-block:: python

   {"lte": {"Cl": 2}}  # At most 2 Cl

**only (Element Restriction)**

Restricts to specified elements only:

.. code-block:: python

   {"only": ["C", "F", "H", "O"]}  # Only these elements allowed

**rel (Relative Ratio)**

Requires specific element ratios:

.. code-block:: python

   # C = (F + H) / 2 - 1
   {
       "rel": {
           "C": {
               "atoms": ["F", "H"],
               "add": -1,
               "div": 2
           }
       }
   }

   # With additional atoms
   {
       "rel": {
           "C": {
               "atoms": ["F"],
               "add": 0.5,
               "div": 2,
               "add_atoms": ["O"]
           }
       }
   }

The relative constraint formula is:

.. math::

   \text{element} = \frac{\sum \text{atoms}}{div} + \text{add} + \sum \text{add\_atoms}

---

.. module:: PFASgroups.PFASDefinitionModel
   :synopsis: PFAS definition data model

PFASDefinition Class
--------------------

.. py:class:: PFASDefinition(id, name, smarts, fluorineRatio, description, **kwargs)

   Represents a broader PFAS definition (e.g., OECD, EPA definitions).

   :param id: Unique identifier
   :type id: int
   :param name: Definition name
   :type name: str
   :param smarts: List of SMARTS patterns
   :type smarts: list[str]
   :param fluorineRatio: Minimum fluorine ratio threshold (optional)
   :type fluorineRatio: float or None
   :param description: Human-readable description
   :type description: str
   :param includeHydrogen: Include H in fluorine ratio calculation (default: True)
   :type includeHydrogen: bool, optional
   :param requireBoth: Require both SMARTS match AND fluorine ratio (default: False)
   :type requireBoth: bool, optional

   **Attributes:**

   .. py:attribute:: id
      :type: int

      Unique identifier.

   .. py:attribute:: name
      :type: str

      Definition name.

   .. py:attribute:: description
      :type: str

      Human-readable description.

   .. py:attribute:: fluorineRatio
      :type: float or None

      Minimum fluorine ratio threshold.

   .. py:attribute:: smarts_strings
      :type: list[str]

      Original SMARTS pattern strings.

   .. py:attribute:: smarts_patterns
      :type: list[rdkit.Chem.Mol]

      Preprocessed SMARTS patterns.

   .. py:attribute:: includeHydrogen
      :type: bool

      Whether to include hydrogen in ratio calculations.

   .. py:attribute:: requireBoth
      :type: bool

      Whether both SMARTS and ratio must match.

   **Methods:**

   .. py:method:: applies_to_molecule(mol_or_smiles, formula=None, **kwargs)

      Check if this PFAS definition applies to a molecule.

      :param mol_or_smiles: RDKit molecule or SMILES string
      :type mol_or_smiles: rdkit.Chem.Mol or str
      :param formula: Pre-computed formula dictionary (optional)
      :type formula: dict, optional
      :returns: True if the definition applies
      :rtype: bool

   .. py:method:: _compute_formula(mol, include_hydrogen)

      Compute molecular formula as a dictionary.

      :param mol: RDKit molecule
      :type mol: rdkit.Chem.Mol
      :param include_hydrogen: Include hydrogen atoms
      :type include_hydrogen: bool
      :returns: Element count dictionary
      :rtype: dict

   .. py:method:: _check_fluorine_ratio(formula, include_hydrogen)

      Check if fluorine ratio meets the threshold.

      :param formula: Element count dictionary
      :type formula: dict
      :param include_hydrogen: Include hydrogen in total count
      :type include_hydrogen: bool
      :returns: True if ratio threshold met
      :rtype: bool

   **Example:**

   .. code-block:: python

      from PFASgroups import PFASDefinition, get_PFASDefinitions
      from rdkit import Chem

      # Get default definitions
      definitions = get_PFASDefinitions()

      # Check if a molecule matches any definition
      mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
      
      for defn in definitions:
          if defn.applies_to_molecule(mol):
              print(f"Matches: {defn.name}")
              print(f"  Description: {defn.description}")

      # Create a custom definition
      custom_def = PFASDefinition(
          id=100,
          name="Custom PFAS Definition",
          smarts=["[C](F)(F)F", "[C](F)(F)[C](F)(F)"],
          fluorineRatio=0.3,
          description="Custom definition requiring 30% fluorine",
          includeHydrogen=True,
          requireBoth=True
      )
