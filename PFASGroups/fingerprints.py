from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from .parser import load_HalogenGroups, parse_groups_in_mol
from typing import Union, List, Optional
import numpy as np


@load_HalogenGroups()
def generate_fingerprint(smiles: Union[str, List[str]],
                             selected_groups: Union[List[int], range, None] = None,
                             representation: str = 'vector',
                             count_mode: str = 'binary',
                             halogens: Union[str, List[str]] = 'F',
                             saturation: Optional[str] = 'per',
                             **kwargs):
    """
    Generate PFAS group fingerprints from SMILES strings.

    Parameters:
    -----------
    smiles : str or list of str
        SMILES string(s) to generate fingerprints for
    selected_groups : list of int, range, or None
        Indices of PFAS groups to include in fingerprint.
        If None, all groups are used. Examples: [28, 29, 30], range(28, 52)
    representation : str, default 'vector'
        Type of fingerprint representation:
        - 'vector': Binary/count vector (numpy array)
        - 'dict': Dictionary mapping group names to counts
        - 'sparse': Dictionary with only non-zero group counts
        - 'detailed': Full match information including chain details
        - 'int': Convert binary to int using int(x,2)
    count_mode : str, default 'binary'
        How to count matches:
        - 'binary': 1 if group present, 0 if absent
        - 'count': Number of matches found
        - 'max_component': Maximum component size for groups with components
    halogens : str or list of str, default 'F'
        Which halogen(s) to use for component SMARTS matching.
        - Single halogen (e.g. 'F'): standard fingerprint of length n_groups
        - List of halogens (e.g. ['F', 'Cl', 'Br', 'I']): fingerprints are generated per
          halogen and stacked, yielding length n_groups * n_halogens.
          Group names are suffixed with '[F]', '[Cl]', etc.
        Available: 'F', 'Cl', 'Br', 'I'
    saturation : str or None, default 'per'
        Saturation filter for component SMARTS matching:
        - 'per': perfluorinated / perhalogenated components only
        - 'poly': polyfluorinated / polyhalogenated components only
        - None: no saturation filter (use all)
        Only applies to groups that use component SMARTS (e.g. OECD groups).
        Groups without component SMARTS (generic/telomer) are unaffected.
    pfas_groups : list, optional
        Custom list of PFAS groups. If None, uses default groups.

    Returns:
    --------
    fingerprints : numpy.ndarray, dict, or list
        Fingerprint representation(s) based on 'representation' parameter.
        For single SMILES input, returns single fingerprint.
        For list input, returns list of fingerprints.
        When multiple halogens are given and representation='vector', the
        vectors are stacked horizontally (shape: n_molecules x (n_groups*n_hal)).
    group_info : dict
        Information about the groups used in fingerprinting:
        - 'group_names': List of group names in fingerprint order
        - 'group_ids': List of group IDs in fingerprint order
        - 'selected_indices': Indices of selected groups
        - 'halogens': Halogen(s) used
        - 'saturation': Saturation filter used

    Examples:
    ---------
    >>> # Binary vector for all groups, fluorine only (default)
    >>> fp, info = generate_fingerprint('CC(F)(F)C(=O)O')

    >>> # Stacked fingerprint for F and Cl
    >>> fp, info = generate_fingerprint('CC(F)(F)C(=O)O', halogens=['F', 'Cl'])

    >>> # Count-based dictionary representation
    >>> fp, info = generate_fingerprint(['CCF', 'CCFF'],
    ...                                  representation='dict',
    ...                                  count_mode='count')
    """
    pfas_groups = kwargs.get('pfas_groups')
    # Normalise halogens to a list
    if isinstance(halogens, str):
        halogens_list = [halogens]
    else:
        halogens_list = list(halogens)

    # Handle single SMILES input
    if isinstance(smiles, str):
        smiles_list = [smiles]
        single_input = True
    else:
        smiles_list = smiles
        single_input = False

    # Determine which groups to use
    if selected_groups is None:
        # Use all groups
        selected_indices = list(range(len(pfas_groups)))
        selected_pfas_groups = pfas_groups
    else:
        # Convert range to list if needed
        if isinstance(selected_groups, range):
            selected_indices = list(selected_groups)
        else:
            selected_indices = selected_groups

        # Validate indices
        max_index = len(pfas_groups) - 1
        invalid_indices = [i for i in selected_indices if i < 0 or i > max_index]
        if invalid_indices:
            raise ValueError(f"Invalid group indices {invalid_indices}. "
                           f"Valid range is 0-{max_index}")

        selected_pfas_groups = [pfas_groups[i] for i in selected_indices]

    # Build group names/ids, suffixed with halogen tag when multiple halogens
    multi_halogen = len(halogens_list) > 1
    if multi_halogen:
        group_names = [f"{g.name} [{h}]" for h in halogens_list for g in selected_pfas_groups]
        group_ids = [g.id for h in halogens_list for g in selected_pfas_groups]
    else:
        group_names = [g.name for g in selected_pfas_groups]
        group_ids = [g.id for g in selected_pfas_groups]

    # Create group info
    group_info = {
        'group_names': group_names,
        'group_ids': group_ids,
        'selected_indices': selected_indices,
        'total_groups': len(pfas_groups),
        'halogens': halogens_list,
        'saturation': saturation,
    }

    # Base kwargs without injected keys
    kwargs_base = {k: v for k, v in kwargs.items() if k != 'pfas_groups'}

    def _parse_for_halogen(mol, formula, halogen):
        """Run parse_groups_in_mol for one halogen."""
        kw = dict(kwargs_base)
        kw['halogens'] = halogen
        if saturation is not None:
            kw['saturation'] = saturation
        all_matches, _ = parse_groups_in_mol(
            mol, formula=formula, pfas_groups=pfas_groups,
            include_PFAS_definitions=False, **kw
        )
        match_dict = {}
        for group, match_count, component_sizes, matched_components in all_matches:
            match_dict[group.id] = {
                'group': group,
                'match_count': match_count,
                'component_sizes': component_sizes,
                'matched_components': matched_components,
            }
        return match_dict

    fingerprints = []

    for smiles_str in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles_str}")

            formula = CalcMolFormula(mol)

            # Collect match dicts per halogen
            match_dicts = {h: _parse_for_halogen(mol, formula, h) for h in halogens_list}
            # First halogen's match_dict used for non-stacked repr (dict/sparse/detailed)
            match_dict = match_dicts[halogens_list[0]]

            def _build_vector(md, groups):
                """Build a count/binary vector for one match_dict + group list."""
                vec = np.zeros(len(groups), dtype=int)
                for i, group in enumerate(groups):
                    if group.id in md:
                        match_info = md[group.id]
                        if count_mode == 'binary' or representation == 'int':
                            vec[i] = 1
                        elif count_mode == 'count':
                            vec[i] = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                vec[i] = max([comp['size'] for comp in match_info['matched_components']])
                            else:
                                vec[i] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                        else:
                            raise ValueError(f"Unknown count_mode: {count_mode}")
                return vec

            if representation == 'vector' or representation == 'int':
                if multi_halogen:
                    # Stack one vector per halogen
                    parts = [_build_vector(match_dicts[h], selected_pfas_groups)
                             for h in halogens_list]
                    fingerprint = np.concatenate(parts)
                else:
                    fingerprint = _build_vector(match_dict, selected_pfas_groups)

                if representation == 'int':
                    fingerprints.append(int(''.join([str(x) for x in fingerprint]), 2))
                else:
                    fingerprints.append(fingerprint)

            elif representation == 'dict':
                # Create dictionary with all groups (including zeros)
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        if count_mode == 'binary':
                            fingerprint[group.name] = 1
                        elif count_mode == 'count':
                            fingerprint[group.name] = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                fingerprint[group.name] = max([comp['size'] for comp in match_info['matched_components']])
                            else:
                                fingerprint[group.name] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                    else:
                        fingerprint[group.name] = 0
                fingerprints.append(fingerprint)

            elif representation == 'sparse':
                # Create dictionary with only non-zero entries
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        value = 0
                        if count_mode == 'binary':
                            value = 1
                        elif count_mode == 'count':
                            value = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                value = max([comp['size'] for comp in match_info['matched_components']])
                            else:
                                value = match_info['match_count'] if match_info['match_count'] > 0 else 0

                        if value > 0:
                            fingerprint[group.name] = value
                fingerprints.append(fingerprint)

            elif representation == 'detailed':
                # Return full match information for selected groups
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        fingerprint[group.name] = {
                            'group_id': group.id,
                            'match_count': match_info['match_count'],
                            'component_sizes': match_info['component_sizes'],
                            'matched_components': match_info['matched_components'],
                            'smarts1': [Chem.MolToSmarts(x) for x in group.smarts if x is not None],
                        }
                fingerprints.append(fingerprint)
            else:
                raise ValueError(f"Unknown representation: {representation}")

        except Exception as e:
            raise ValueError(f"Error processing SMILES '{smiles_str}': {str(e)}")

    # Return single fingerprint for single input, list for multiple inputs
    if single_input:
        return fingerprints[0], group_info
    else:
        return fingerprints, group_info
