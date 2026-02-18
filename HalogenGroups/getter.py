from .parser import add_componentSmarts, load_HalogenGroups, load_PFASDefinitions
from .core import HALOGEN_GROUPS_FILE,PFAS_DEFINITIONS_FILE
import json
@add_componentSmarts()
def get_componentSmartss(**kwargs):
    return kwargs.get('componentSmartss')

@load_HalogenGroups()
def get_HalogenGroups(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_groups',[]).extend(kwargs.get('agg_pfas_groups',[]))
    with open(HALOGEN_GROUPS_FILE, 'r') as f:
        data = json.load(f)
    return data

@load_PFASDefinitions()
def get_PFASDefinitions(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_definitions',[])
    with open(PFAS_DEFINITIONS_FILE, 'r') as f:
        data = json.load(f)
    return data
