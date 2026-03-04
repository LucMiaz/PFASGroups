import sys
print("Python:", sys.version, flush=True)
import numpy as np
print("numpy:", np.__version__, flush=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
print("matplotlib:", matplotlib.__version__, flush=True)
from scipy.stats import kruskal
print("scipy: OK", flush=True)
print("ALL OK", flush=True)
