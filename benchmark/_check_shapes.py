import numpy as np
import pathlib

caches = sorted(pathlib.Path(r'c:\Users\luc\git\PFASGroups\benchmark\data').glob('pfg_*_cache.npy'))
for p in caches:
    a = np.load(p)
    print(f'{p.name:<50} shape={a.shape}  n_cols={a.shape[1]}')
