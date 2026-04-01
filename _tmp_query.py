from PFASGroups.getter import get_compiled_HalogenGroups, get_HalogenGroups
from PFASGroups.HalogenGroupModel import HalogenGroup
compiled = get_compiled_HalogenGroups()
raw = get_HalogenGroups()
compiled_ids = set(g.id for g in compiled)
excluded = [g for g in raw if g['id'] not in compiled_ids]
print('Compiled groups: %d' % len(compiled))
print('Total raw groups: %d' % len(raw))
print('Excluded (not compiled):')
for g in sorted(excluded, key=lambda x: x['id']):
    hal = g.get('componentHalogens') or g.get('componentHalogen', '')
    cmt = g.get('compute', True)
    print('  id=%3d  compute=%-5s  halogens=%-12s  name=%s' % (g['id'], str(cmt), str(hal), g['name']))
