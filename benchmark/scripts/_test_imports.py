import sys, traceback
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')

try:
    print("numpy", flush=True)
    import numpy as np
    print("numpy OK", flush=True)

    print("pandas", flush=True)
    import pandas as pd
    print("pandas OK", flush=True)

    print("matplotlib", flush=True)
    import matplotlib
    print("matplotlib imported", flush=True)
    matplotlib.use('Agg')
    print("matplotlib Agg OK", flush=True)

    import matplotlib.pyplot as plt
    print("plt OK", flush=True)

    import matplotlib.patches as mpatches
    print("patches OK", flush=True)

    import matplotlib.gridspec as gridspec
    print("gridspec OK", flush=True)

    from matplotlib.colors import ListedColormap
    print("ListedColormap OK", flush=True)

    import matplotlib.ticker as mticker
    print("mticker OK", flush=True)

    print("scipy", flush=True)
    from scipy.stats import kruskal, mannwhitneyu
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist
    print("scipy OK", flush=True)

    print("All imports OK", flush=True)

except Exception:
    traceback.print_exc()
