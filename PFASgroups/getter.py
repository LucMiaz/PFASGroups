from .parser import add_smartsPath, load_PFASGroups, load_PFASDefinitions
from .core import PFAS_GROUPS_FILE,PFAS_DEFINITIONS_FILE
import json
@add_smartsPath()
def get_smartsPaths(**kwargs):
    return kwargs.get('smartsPaths')

@load_PFASGroups()
def get_PFASGroups(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_groups',[]).extend(kwargs.get('agg_pfas_groups',[]))
    with open(PFAS_GROUPS_FILE, 'r') as f:
        data = json.load(f)
    return data

@load_PFASDefinitions()
def get_PFASDefinitions(**kwargs):
    if kwargs.get('json_format',False):
        return kwargs.get('pfas_definitions',[])
    with open(PFAS_DEFINITIONS_FILE, 'r') as f:
        data = json.load(f)
    return data