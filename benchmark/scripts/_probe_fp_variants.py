import sys
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')
from HalogenGroups import generate_fingerprint

smi = 'FC(F)(F)C(F)(F)C(F)(F)C(=O)O'
fp_bin, info = generate_fingerprint(smi, representation='vector', count_mode='binary', halogens='F', saturation='per')
print('binary shape:', fp_bin.shape, 'active:', int(fp_bin.sum()))

fp_mh, info_mh = generate_fingerprint(smi, representation='vector', count_mode='binary', halogens=['F','Cl'], saturation='per')
print('multi-hal (F+Cl) shape:', fp_mh.shape)

fp_poly, _ = generate_fingerprint(smi, representation='vector', count_mode='binary', halogens='F', saturation='poly')
print('poly active:', int(fp_poly.sum()))

fp_none, _ = generate_fingerprint(smi, representation='vector', count_mode='binary', halogens='F', saturation=None)
print('saturation=None active:', int(fp_none.sum()))

fp_mx, _ = generate_fingerprint(smi, representation='vector', count_mode='max_component', halogens='F', saturation='per')
print('max_comp nonzero:', fp_mx[fp_mx>0].tolist())

# detailed repr
fp_det, _ = generate_fingerprint(smi, representation='detailed', halogens='F', saturation='per', count_mode='binary')
k0 = list(fp_det.keys())[0]
print('detailed[0] keys:', list(fp_det[k0].keys()))
print('component_sizes sample:', fp_det[k0].get('component_sizes'))
mc = fp_det[k0].get('matched_components')
if mc:
    print('matched_components[0] keys:', list(mc[0].keys()))
    print('matched_components[0]:', mc[0])

# test mixed-halogen molecule
smi_mixed = 'FC(F)(F)C(F)(F)C(=O)O.ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O'
fp_mixed, info_mixed = generate_fingerprint(smi_mixed, representation='vector', count_mode='binary', halogens=['F','Cl'], saturation='per')
print('\nmixed compound (F+Cl) shape:', fp_mixed.shape, 'active:', int(fp_mixed.sum()))
active_groups = [info_mixed['group_names'][i] for i in range(len(fp_mixed)) if fp_mixed[i]>0]
print('active groups:', active_groups[:6])
