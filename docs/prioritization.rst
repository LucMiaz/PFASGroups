Prioritization
==============

:func:`~PFASGroups.prioritise_molecules` (alias
:func:`~PFASGroups.prioritize_molecules`) ranks a set of molecules for
screening, in one of two modes:

* **without a reference**, by their intrinsic fluorination: total size of the
  fluorinated components and size of the largest ones;
* **with a reference set**, by how similar their group profile is to the
  reference (for example a list of known priority PFAS).

.. contents:: Contents
   :local:
   :depth: 2

Basic usage
-----------

.. code-block:: python

   from PFASGroups import prioritise_molecules

   molecules = [
       "CCCC(F)(F)F",
       "FC(F)(F)C(=O)O",
       "FCCC(F)(F)F",
       "ClCCCl",
       "BrCCBr",
   ]

   ranked, scores = prioritise_molecules(molecules)

``ranked`` is a :class:`~PFASGroups.PFASEmbeddingSet` sorted from highest to
lowest priority, and ``scores`` a NumPy array of the scores in the same order:

.. code-block:: python

   for mol, score in zip(ranked, scores):
       print(f"{mol.smiles:40s}  {score:.1f}")

With ``return_scores=False`` only the sorted ``PFASEmbeddingSet`` is returned.
``molecules`` can also be a list of RDKit molecules or an already parsed
``PFASEmbeddingSet``, which avoids parsing twice.

Intrinsic prioritization (no reference)
---------------------------------------

Without a reference the score of molecule :math:`i` is

.. math::

   \text{score}(i) = a \cdot \sum_k s_{ik} + b \cdot P_p(s_{i\cdot})

where :math:`s_{ik}` are the sizes (number of halogenated carbons) of the
fluorinated components of molecule :math:`i` and :math:`P_p` is their
``percentile``-th percentile.  ``a`` weights the total fluorination and ``b``
the presence of long chains; molecules without a fluorinated component score 0.

.. code-block:: python

   # Emphasise long chains
   ranked, scores = prioritise_molecules(molecules, a=0.5, b=2.0, percentile=95)

   # Emphasise total fluorine content
   ranked = prioritise_molecules(molecules, a=2.0, b=0.5, return_scores=False)

Similarity to a reference set
-----------------------------

With ``reference``, each molecule is embedded with
:meth:`~PFASGroups.PFASEmbeddingSet.to_array` and scored by the cosine
similarity between its vector and the mean vector of the reference set:

.. math::

   \text{score}(i) = \frac{\mathbf{f}_i \cdot \bar{\mathbf{f}}_\text{ref}}
                          {\lVert \mathbf{f}_i \rVert \, \lVert \bar{\mathbf{f}}_\text{ref} \rVert}

Scores lie in :math:`[0, 1]`; molecules that activate the same groups as the
reference rank highest, and molecules without any group match score 0.

.. code-block:: python

   known_priority = ["FC(F)(F)C(F)(F)C(F)(F)C(=O)O"]
   candidates = ["FC(F)(F)C(F)(F)C(=O)O", "FC(F)(F)C(F)(F)S(=O)(=O)O", "ClCCCl"]

   ranked, scores = prioritise_molecules(candidates, reference=known_priority)
   print(ranked[0].smiles)   # the perfluoroalkyl carboxylic acid ranks first

``group_selection`` and ``count_mode`` choose the embedding used for the
comparison.

Parameters
----------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default
     - Description
   * - ``molecules``
     - *required*
     - SMILES strings, RDKit molecules or a ``PFASEmbeddingSet``
   * - ``reference``
     - ``None``
     - Reference molecules (same types as ``molecules``); switches to the
       similarity mode
   * - ``group_selection``
     - ``'all'``
     - Groups in the embedding (reference mode): ``'all'``, ``'oecd'``,
       ``'generic'``, ``'telomers'`` or ``'generic+telomers'``
   * - ``count_mode``
     - ``'max_component'``
     - Embedding encoding (reference mode): ``'binary'``, ``'count'`` or
       ``'max_component'``
   * - ``halogens``
     - ``'F'``
     - Halogens used when parsing the reference
   * - ``saturation``
     - ``None``
     - ``'per'``, ``'poly'`` or ``None`` (both) when parsing the reference
   * - ``a``
     - ``1.0``
     - Weight of the total component size (intrinsic mode)
   * - ``b``
     - ``1.0``
     - Weight of the component-size percentile (intrinsic mode)
   * - ``percentile``
     - ``90.0``
     - Percentile of the component sizes (intrinsic mode)
   * - ``return_scores``
     - ``True``
     - Return ``(ranked, scores)``; if ``False``, only ``ranked``
   * - ``ascending``
     - ``False``
     - If ``True``, lowest priority first
   * - ``progress``
     - ``False``
     - Show a progress bar while parsing

Summary statistics
------------------

:func:`~PFASGroups.get_priority_statistics` summarises a ranking (mean,
spread, top-*n* molecules and their most common groups):

.. code-block:: python

   from PFASGroups import get_priority_statistics

   stats = get_priority_statistics(ranked, scores, top_n=3)
   print(stats["score_mean"], stats["top_n_smiles"])

Example: shortlist from a CSV file
----------------------------------

.. code-block:: python

   import pandas as pd
   from PFASGroups import prioritise_molecules

   df = pd.read_csv("pfas_candidates.csv")
   ranked, scores = prioritise_molecules(df["SMILES"].tolist())

   # Top 50 by intrinsic fluorination
   top50 = [mol.smiles for mol in ranked[:50]]
