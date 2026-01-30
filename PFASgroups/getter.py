from .parser import add_smartsPath, load_PFASGroups, load_PFASDefinitions
@add_smartsPath()
def get_smartsPaths(**kwargs):
    return kwargs.get('smartsPaths')

@load_PFASGroups()
def get_PFASGroups(**kwargs):
    return kwargs.get('pfas_groups',[]).extend(kwargs.get('agg_pfas_groups',[]))

@load_PFASDefinitions()
def get_PFASDefinitions(**kwargs):
    return kwargs.get('pfas_definitions')