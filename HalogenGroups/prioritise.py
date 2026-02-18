"""
PFAS Molecule Prioritization Module
====================================

This module provides functions to prioritize PFAS molecules based on:
1. Distributional similarity to a reference list (using fingerprint KL divergence)
2. Fluorinated component characteristics (total fluorine content and component size distribution)

Author: HalogenGroups
Version: 2.2.4
"""

from typing import Union, List, Tuple, Optional, Dict
import numpy as np
from rdkit import Chem

from .parser import parse_smiles, parse_mols
from .results_model import ResultsModel


def prioritise_molecules(
    molecules: Union[List[str], List[Chem.Mol], ResultsModel],
    reference: Optional[Union[List[str], List[Chem.Mol], ResultsModel]] = None,
    group_selection: str = 'all',
    count_mode: str = 'binary',
    a: float = 1.0,
    b: float = 1.0,
    percentile: float = 90.0,
    return_scores: bool = True,
    ascending: bool = False
) -> Union[ResultsModel, Tuple[ResultsModel, np.ndarray]]:
    """
    Prioritize PFAS molecules based on similarity to a reference or intrinsic properties.

    Parameters
    ----------
    molecules : list of str, list of rdkit.Chem.Mol, or ResultsModel
        Molecules to prioritize. Can be:
        - List of SMILES strings
        - List of RDKit molecule objects
        - ResultsModel object (pre-computed results)

    reference : list of str, list of rdkit.Chem.Mol, ResultsModel, or None
        Reference molecules for similarity comparison. If provided, molecules are
        prioritized by distributional similarity (lower KL divergence = higher priority).
        If None, prioritization is based on intrinsic fluorinated component properties.

    group_selection : str, default 'all'
        PFAS group selection for fingerprint generation when using reference:
        - 'all': All 115 groups (OECD + generic)
        - 'oecd': OECD-defined groups (1-28)
        - 'generic': Generic functional groups (29-115)
        - 'telomers': Telomer-related groups
        - 'generic+telomers': Combined selection

    count_mode : str, default 'binary'
        Fingerprint encoding mode when using reference:
        - 'binary': 1 if present, 0 if absent
        - 'count': Number of matches
        - 'max_component': Maximum component size

    a : float, default 1.0
        Weight for total fluorinated component size (sum of all component sizes).
        Used when reference is None. Higher values prioritize molecules with more
        total fluorination.

    b : float, default 1.0
        Weight for component size percentile. Used when reference is None.
        Higher values prioritize molecules with larger individual fluorinated components.

    percentile : float, default 90.0
        Percentile value (0-100) for component size distribution. Used when reference
        is None. Common values:
        - 90.0: Focus on largest 10% of components
        - 75.0: Focus on largest 25% of components
        - 50.0: Median component size

    return_scores : bool, default True
        If True, returns tuple of (prioritized_results, scores).
        If False, returns only prioritized_results.

    ascending : bool, default False
        Sort order. If False (default), highest priority first.
        If True, lowest priority first.

    Returns
    -------
    ResultsModel or tuple
        If return_scores=True: (prioritized_results, scores)
        If return_scores=False: prioritized_results only

        prioritized_results : ResultsModel
            Molecules sorted by priority
        scores : np.ndarray
            Priority scores for each molecule

    Examples
    --------
    # Priority by similarity to reference list
    >>> from HalogenGroups import prioritise_molecules
    >>> inventory = ["FC(F)(F)C(F)(F)C(=O)O", "FC(F)(F)C(F)(F)C(F)(F)C(=O)O"]
    >>> reference = ["FC(F)(F)C(F)(F)C(=O)O"]  # Known priority compounds
    >>> results, scores = prioritise_molecules(inventory, reference=reference)
    >>> print(f"Most similar: {results[0]['smiles']}")

    # Priority by fluorination characteristics
    >>> results, scores = prioritise_molecules(
    ...     inventory,
    ...     a=1.0,  # Weight for total fluorination
    ...     b=2.0,  # Weight for largest components
    ...     percentile=90
    ... )
    >>> print(f"Highest priority: {results[0]['smiles']}")

    # Emphasize total fluorine content
    >>> results = prioritise_molecules(inventory, a=2.0, b=0.5, return_scores=False)

    # Focus on molecules with longest chains
    >>> results = prioritise_molecules(inventory, a=0.5, b=2.0, percentile=95)

    Notes
    -----
    **Reference-based prioritization:**

    Uses distributional similarity (KL divergence) between fingerprints:
    - Lower KL divergence = more similar to reference = higher priority
    - Useful for prioritizing molecules similar to known compounds of concern
    - Accounts for overall composition of PFAS groups

    **Intrinsic prioritization (no reference):**

    Score = a × Σ(component_sizes) + b × percentile(component_sizes, p)

    Where:
    - Σ(component_sizes): Total number of fluorinated carbons across all components
    - percentile(component_sizes, p): Size of fluorinated components at pth percentile

    This approach prioritizes molecules based on:
    - Total fluorination burden (a parameter)
    - Presence of long perfluorinated chains (b and percentile parameters)

    **Tuning guidelines:**

    For environmental persistence concerns:
    - High b, high percentile (e.g., b=2.0, p=90): Long-chain compounds

    For bioaccumulation potential:
    - Balanced a and b (e.g., a=1.0, b=1.0): Both total and chain length

    For screening priority:
    - High a, moderate b (e.g., a=2.0, b=1.0, p=75): Total fluorine load

    See Also
    --------
    ResultsModel.to_fingerprint : Convert results to fingerprints
    ResultsFingerprint.compare_kld : Compare fingerprint distributions
    """
    # Step 1: Convert input to ResultsModel if needed
    if isinstance(molecules, ResultsModel):
        results = molecules
    elif isinstance(molecules, list):
        if len(molecules) == 0:
            raise ValueError("molecules list is empty")

        # Check if list contains SMILES strings or Mol objects
        if isinstance(molecules[0], str):
            results = parse_smiles(molecules)
        elif isinstance(molecules[0], Chem.Mol):
            results = parse_mols(molecules)
        else:
            raise TypeError(
                f"molecules must be list of SMILES strings or RDKit Mol objects, "
                f"got {type(molecules[0])}"
            )
    else:
        raise TypeError(
            f"molecules must be list of SMILES/Mol objects or ResultsModel, "
            f"got {type(molecules)}"
        )

    # Step 2: Compute priority scores
    if reference is not None:
        # Reference-based: distributional similarity using KL divergence
        scores = _prioritise_by_reference(
            results, reference, group_selection, count_mode
        )
    else:
        # Intrinsic: fluorinated component characteristics
        scores = _prioritise_by_components(results, a, b, percentile)

    # Step 3: Sort by priority
    sorted_indices = np.argsort(scores)
    if not ascending:
        sorted_indices = sorted_indices[::-1]

    # Create prioritized results
    prioritized_results = ResultsModel([results[i] for i in sorted_indices])

    if return_scores:
        return prioritized_results, scores[sorted_indices]
    else:
        return prioritized_results


def _prioritise_by_reference(
    results: ResultsModel,
    reference: Union[List[str], List[Chem.Mol], ResultsModel],
    group_selection: str,
    count_mode: str
) -> np.ndarray:
    """
    Prioritize by distributional similarity to reference using KL divergence.

    Lower KL divergence = higher similarity = higher priority score (inverted).

    Parameters
    ----------
    results : ResultsModel
        Molecules to prioritize
    reference : list of str, list of rdkit.Chem.Mol, or ResultsModel
        Reference molecules
    group_selection : str
        PFAS group selection
    count_mode : str
        Fingerprint encoding mode

    Returns
    -------
    np.ndarray
        Priority scores (higher = higher priority)
    """
    # Convert reference to ResultsModel if needed
    if isinstance(reference, ResultsModel):
        ref_results = reference
    elif isinstance(reference, list):
        if len(reference) == 0:
            raise ValueError("reference list is empty")

        if isinstance(reference[0], str):
            ref_results = parse_smiles(reference)
        elif isinstance(reference[0], Chem.Mol):
            ref_results = parse_mols(reference)
        else:
            raise TypeError(
                f"reference must be list of SMILES strings or RDKit Mol objects, "
                f"got {type(reference[0])}"
            )
    else:
        raise TypeError(
            f"reference must be list of SMILES/Mol objects or ResultsModel, "
            f"got {type(reference)}"
        )

    # Convert to fingerprints
    ref_fp = ref_results.to_fingerprint(
        group_selection=group_selection,
        count_mode=count_mode
    )

    # Compute KL divergence for each molecule
    kl_divergences = []

    for i in range(len(results)):
        # Create single-molecule fingerprint
        single_result = ResultsModel([results[i]])
        single_fp = single_result.to_fingerprint(
            group_selection=group_selection,
            count_mode=count_mode
        )

        # Compute KL divergence (minmax method)
        try:
            kl_div = ref_fp.compare_kld(single_fp, method='minmax')
            kl_divergences.append(kl_div)
        except Exception as e:
            # If comparison fails, assign worst score
            kl_divergences.append(1.0)

    kl_divergences = np.array(kl_divergences)

    # Convert to priority scores: lower KL = higher priority
    # Invert so that higher score = higher priority
    priority_scores = 1.0 - kl_divergences

    return priority_scores


def _prioritise_by_components(
    results: ResultsModel,
    a: float,
    b: float,
    percentile: float
) -> np.ndarray:
    """
    Prioritize by fluorinated component characteristics.

    Score = a × sum(component_sizes) + b × percentile(component_sizes, p)

    Parameters
    ----------
    results : ResultsModel
        Molecules to prioritize
    a : float
        Weight for total component size
    b : float
        Weight for component size percentile
    percentile : float
        Percentile value (0-100)

    Returns
    -------
    np.ndarray
        Priority scores (higher = higher priority)
    """
    scores = []

    for mol_result in results:
        # Extract all fluorinated component sizes from all matched groups
        all_component_sizes = []

        for match in mol_result['matches']:
            if match.get('type') == 'HalogenGroup' and 'components_sizes' in match:
                component_sizes = match['components_sizes']
                if isinstance(component_sizes, list) and len(component_sizes) > 0:
                    all_component_sizes.extend(component_sizes)

        if len(all_component_sizes) == 0:
            # No fluorinated components found
            score = 0.0
        else:
            all_component_sizes = np.array(all_component_sizes)

            # Compute score components
            total_size = float(np.sum(all_component_sizes))
            percentile_size = float(np.percentile(all_component_sizes, percentile))

            # Combined score
            score = a * total_size + b * percentile_size

        scores.append(score)

    return np.array(scores)


def get_priority_statistics(
    results: ResultsModel,
    scores: np.ndarray,
    top_n: int = 10
) -> Dict:
    """
    Get statistics about prioritization results.

    Parameters
    ----------
    results : ResultsModel
        Prioritized molecules
    scores : np.ndarray
        Priority scores
    top_n : int, default 10
        Number of top molecules to analyze

    Returns
    -------
    dict
        Statistics including:
        - 'n_molecules': Total number of molecules
        - 'score_mean': Mean priority score
        - 'score_std': Standard deviation of scores
        - 'score_min': Minimum score
        - 'score_max': Maximum score
        - 'top_n_smiles': SMILES of top N molecules
        - 'top_n_scores': Scores of top N molecules
        - 'top_n_groups': Most common groups in top N molecules

    Examples
    --------
    >>> results, scores = prioritise_molecules(smiles_list, reference=ref_list)
    >>> stats = get_priority_statistics(results, scores, top_n=10)
    >>> print(f"Mean score: {stats['score_mean']:.3f}")
    >>> print(f"Top molecule: {stats['top_n_smiles'][0]}")
    """
    stats = {
        'n_molecules': len(results),
        'score_mean': float(np.mean(scores)),
        'score_std': float(np.std(scores)),
        'score_min': float(np.min(scores)),
        'score_max': float(np.max(scores)),
        'score_median': float(np.median(scores)),
    }

    # Top N molecules
    n = min(top_n, len(results))
    stats['top_n_smiles'] = [results[i]['smiles'] for i in range(n)]
    stats['top_n_scores'] = scores[:n].tolist()

    # Most common groups in top N
    group_counts = {}
    for i in range(n):
        for match in results[i]['matches']:
            if match.get('type') == 'HalogenGroup':
                group_name = match.get('group_name', 'Unknown')
                group_counts[group_name] = group_counts.get(group_name, 0) + 1

    # Sort by frequency
    sorted_groups = sorted(group_counts.items(), key=lambda x: x[1], reverse=True)
    stats['top_n_groups'] = sorted_groups[:10]  # Top 10 most common groups

    return stats


# Alias for American spelling
prioritize_molecules = prioritise_molecules
