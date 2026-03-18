from PFASGroups import PFASFingerprint, parse_smiles
fp = PFASFingerprint(['FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F'])
print('n_molecules:', fp.n_molecules)
print('has_cache:', fp.has_cache)
print('match_cache len:', len(fp.match_cache))
e = fp.get_embedding(component_metrics=['binary'])
import numpy as np
mat = np.atleast_2d(np.asarray(e, dtype=float))
print('embedding shape:', mat.shape)
e2 = fp.get_embedding(component_metrics=['min_dist_to_barycenter'])
mat2 = np.atleast_2d(np.asarray(e2, dtype=float))
print('min_dist_to_barycenter shape:', mat2.shape)
print('OK')

