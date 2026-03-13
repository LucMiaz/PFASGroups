from .parser import add_componentSmarts, load_HalogenGroups, load_PFASDefinitions
from .core import HALOGEN_GROUPS_FILE, PFAS_DEFINITIONS_FILE
from .HalogenGroupModel import HalogenGroup
import json
@add_componentSmarts()
def get_componentSMARTSs(**kwargs):
    return kwargs.get('componentSmartss')

@load_HalogenGroups()
def get_HalogenGroups(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_groups',[]).extend(kwargs.get('agg_pfas_groups',[]))
    with open(HALOGEN_GROUPS_FILE, 'r') as f:
        data = json.load(f)
    return data

@load_HalogenGroups()
def get_compiled_HalogenGroups(**kwargs):
    """Return compiled HalogenGroup instances (compute=True groups only).

    Unlike :func:`get_HalogenGroups` which returns raw JSON dicts, this
    function returns ready-to-use :class:`HalogenGroup` instances that can
    be directly passed to :func:`parse_smiles` or extended with custom groups.

    Returns
    -------
    list of HalogenGroup
        All compiled groups, suitable for passing as *pfas_groups* to
        :func:`parse_smiles` or :class:`~PFASGroups.fingerprints.PFASFingerprint`.

    Examples
    --------
    >>> from PFASgroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles
    >>> groups = get_compiled_HalogenGroups()
    >>> groups.append(HalogenGroup(
    ...     id=200, name="Perfluoroalkyl nitrates",
    ...     smarts={"[C$(C[ON+](=O)[O-])]": 1},
    ...     componentSaturation="per", componentHalogens="F",
    ...     componentForm="alkyl",
    ...     constraints={"eq": {"N": 1}, "gte": {"F": 1}},
    ... ))
    >>> results = parse_smiles(["FC(F)(F)C(F)(F)ON(=O)=O"], pfas_groups=groups)
    """
    return list(kwargs.get('pfas_groups', []))

def get_compiled_PFASGroups() -> list:
    """Return compiled HalogenGroup instances restricted to fluorine (PFAS).

    Similar to :func:`get_compiled_HalogenGroups` but every group has
    ``componentHalogens`` forced to ``'F'``, so the compiled component-SMARTS
    patterns are built for fluorine only.  Suitable for extending with custom
    PFAS-specific groups and passing directly to :func:`parse_smiles` or
    :class:`~PFASGroups.fingerprints.PFASFingerprint` when analysing PFAS/fluorinated compounds.

    Returns
    -------
    list of HalogenGroup
        All compute=True groups compiled with ``componentHalogens='F'``,
        suitable for passing as *pfas_groups* to :func:`parse_smiles` or
        :class:`~PFASGroups.fingerprints.PFASFingerprint`.

    Examples
    --------
    >>> from PFASGroups import get_compiled_PFASGroups, HalogenGroup, parse_smiles
    >>> groups = get_compiled_PFASGroups()
    >>> groups.append(HalogenGroup(
    ...     id=200, name="Perfluoroalkyl nitrates",
    ...     smarts={"[C$(C[ON+](=O)[O-])]": 1},
    ...     componentSmarts="Perfluoroalkyl",
    ...     componentSaturation="per", componentHalogens="F",
    ...     componentForm="alkyl",
    ...     constraints={"eq": {"N": 1}, "gte": {"F": 1}},
    ... ))
    >>> results = parse_smiles(["FC(F)(F)C(F)(F)ON(=O)=O"], pfas_groups=groups)

    See Also
    --------
    get_compiled_HalogenGroups : same but supports all halogens (F, Cl, Br, I).
    """
    with open(HALOGEN_GROUPS_FILE, 'r') as f:
        _pfg = json.load(f)
    return [
        HalogenGroup(**{**entry, 'componentHalogens': 'F'})
        for entry in _pfg
        if entry.get('compute', True)
    ]


@load_HalogenGroups()
def get_PFASGroups(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_groups',[]).extend(kwargs.get('agg_pfas_groups',[]))
    with open(HALOGEN_GROUPS_FILE, 'r') as f:
        data = json.load(f)
    for entry in data:
        entry["componentHalogen"] = ['F']
    return data

@load_PFASDefinitions()
def get_PFASDefinitions(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_definitions',[])
    with open(PFAS_DEFINITIONS_FILE, 'r') as f:
        data = json.load(f)
    return data
