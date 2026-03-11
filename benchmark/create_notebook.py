"""Generate fingerprint_structure_analysis.ipynb programmatically."""
import json
from pathlib import Path

OUT = Path(__file__).parent / "fingerprint_structure_analysis.ipynb"


def code(src: str, id_: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": id_,
        "metadata": {},
        "outputs": [],
        "source": src,
    }


def md(src: str, id_: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": id_,
        "metadata": {},
        "source": src,
    }


cells = []

# ── Cell 0 ── title ──────────────────────────────────────────────────────────
cells.append(md(
    "# PFASGroups fingerprint — chain-length × branching-degree analysis\n"
    "\n"
    "Systematic study of how all 39 fingerprint configurations respond to:\n"
    "1. **Chain length** — total fluorinated carbons C2–C8\n"
    "2. **Branching degree** — linear / α-mono / β-mono / γ-mono / gem-di / α+β double\n"
    "\n"
    "Analyses: overall discrimination ranking, chain-length sensitivity,\n"
    "branching discrimination, MDS projections, per-metric impact, heatmaps.\n",
    "md-title",
))

# ── Cell 1 ── imports ────────────────────────────────────────────────────────
cells.append(code(
    "from __future__ import annotations\n"
    "import sys, warnings\n"
    "from pathlib import Path\n"
    "\n"
    "import numpy as np\n"
    "import pandas as pd\n"
    "import matplotlib\n"
    "import matplotlib.pyplot as plt\n"
    "import matplotlib.patches as mpatches\n"
    "from matplotlib.colors import LinearSegmentedColormap\n"
    "from scipy.cluster.hierarchy import linkage, leaves_list\n"
    "from scipy.spatial.distance import squareform\n"
    "\n"
    "# ── Add repo root to sys.path ─────────────────────────────────────────────\n"
    "NOTEBOOK_DIR = Path(globals().get('__vsc_ipynb_file__', __file__)).parent\n"
    "REPO_ROOT = NOTEBOOK_DIR.parent\n"
    "if str(REPO_ROOT) not in sys.path:\n"
    "    sys.path.insert(0, str(REPO_ROOT))\n"
    "\n"
    "from HalogenGroups import parse_smiles\n"
    "\n"
    "warnings.filterwarnings('ignore')\n"
    "%matplotlib inline\n"
    "print('Imports OK')\n",
    "cell-imports",
))

# ── Cell 2 ── molecule grid ──────────────────────────────────────────────────
cells.append(md(
    "## 1. Molecule grid\n\n"
    "All molecules are PFCA (perfluoroalkyl carboxylic acids, `OC(=O)` head group).\n"
    "The fluorinated tail varies on two axes: chain length (C2–C8) and branching type.\n",
    "md-molecules",
))

cells.append(code(
    "def _cf2(n): return 'C(F)(F)' * n\n"
    "\n"
    "def linear(n):     # n = total fluorinated C\n"
    "    assert n >= 1\n"
    "    return 'OC(=O)' + _cf2(n - 1) + 'C(F)(F)F'\n"
    "\n"
    "def alpha(n):      # CF3 branch at alpha-C; requires n >= 3\n"
    "    assert n >= 3\n"
    "    tail = _cf2(n - 3) + 'C(F)(F)F' if n > 3 else ''\n"
    "    return 'OC(=O)C(F)(C(F)(F)F)' + (tail if tail else 'F')\n"
    "\n"
    "def beta(n):       # CF3 branch at beta-C; requires n >= 4\n"
    "    assert n >= 4\n"
    "    tail = _cf2(n - 4) + 'C(F)(F)F' if n > 4 else ''\n"
    "    return 'OC(=O)C(F)(F)C(F)(C(F)(F)F)' + (tail if tail else 'F')\n"
    "\n"
    "def gamma(n):      # CF3 branch at gamma-C; requires n >= 5\n"
    "    assert n >= 5\n"
    "    tail = _cf2(n - 5) + 'C(F)(F)F' if n > 5 else ''\n"
    "    return 'OC(=O)C(F)(F)C(F)(F)C(F)(C(F)(F)F)' + (tail if tail else 'F')\n"
    "\n"
    "def gem(n):        # two CF3 at alpha-C; requires n >= 4\n"
    "    assert n >= 4\n"
    "    tail = _cf2(n - 4) + 'C(F)(F)F' if n > 4 else ''\n"
    "    return 'OC(=O)C(C(F)(F)F)(C(F)(F)F)' + (tail if tail else 'F')\n"
    "\n"
    "def double_ab(n):  # CF3 at alpha AND beta; requires n >= 5\n"
    "    assert n >= 5\n"
    "    tail = _cf2(n - 5) + 'C(F)(F)F' if n > 5 else ''\n"
    "    return 'OC(=O)C(F)(C(F)(F)F)C(F)(C(F)(F)F)' + (tail if tail else 'F')\n"
    "\n"
    "BRANCHING_MAKERS = {\n"
    "    'linear':    (linear,    range(2, 9)),\n"
    "    'alpha':     (alpha,     range(3, 9)),\n"
    "    'beta':      (beta,      range(4, 9)),\n"
    "    'gamma':     (gamma,     range(5, 9)),\n"
    "    'gem':       (gem,       range(4, 9)),\n"
    "    'double_ab': (double_ab, range(5, 9)),\n"
    "}\n"
    "\n"
    "records = []\n"
    "for btype, (maker, rng) in BRANCHING_MAKERS.items():\n"
    "    for n in rng:\n"
    "        records.append(dict(\n"
    "            label=f'C{n}-{btype}', smiles=maker(n),\n"
    "            n_fluoro_c=n, branching=btype,\n"
    "        ))\n"
    "MOL_DF = pd.DataFrame(records)\n"
    "print(f'{len(MOL_DF)} molecules')\n"
    "print(MOL_DF.groupby('branching')['label'].count())\n"
    "MOL_DF\n",
    "cell-molecules",
))

# ── Cell 3 ── draw grid ──────────────────────────────────────────────────────
cells.append(code(
    "from rdkit import Chem\n"
    "from rdkit.Chem import Draw\n"
    "mols   = [Chem.MolFromSmiles(s) for s in MOL_DF['smiles']]\n"
    "labels_list = list(MOL_DF['label'])\n"
    "img = Draw.MolsToGridImage(mols, molsPerRow=7, subImgSize=(200, 130),\n"
    "                           legends=labels_list, returnPNG=False)\n"
    "display(img)\n",
    "cell-draw",
))

# ── Cell 4 ── FP configs ─────────────────────────────────────────────────────
cells.append(md("## 2. Fingerprint configurations (all 39)\n", "md-fp-configs"))

cells.append(code(
    "ALL_GRAPH_METRICS = [\n"
    "    'branching', 'mean_eccentricity', 'median_eccentricity', 'diameter',\n"
    "    'radius', 'effective_graph_resistance', 'component_fraction',\n"
    "    'min_dist_to_center', 'max_dist_to_periphery', 'min_dist_to_barycenter',\n"
    "    'min_resistance_dist_to_barycenter', 'min_resistance_dist_to_center',\n"
    "    'max_resistance_dist_to_periphery', 'size',\n"
    "]\n"
    "ALL_MOL_METRICS = [\n"
    "    'n_components', 'total_size', 'mean_size', 'max_size',\n"
    "    'mean_branching', 'max_branching', 'mean_eccentricity',\n"
    "    'max_diameter', 'mean_component_fraction', 'max_component_fraction',\n"
    "]\n"
    "\n"
    "def build_configs():\n"
    "    cfgs = []\n"
    "    for mode in ('binary', 'count', 'max_component', 'total_component'):\n"
    "        cfgs.append(('Count modes', mode, dict(count_mode=mode)))\n"
    "    for m in ALL_GRAPH_METRICS:\n"
    "        cfgs.append(('Graph metrics', f'b+{m}',\n"
    "                     dict(count_mode='binary', graph_metrics=[m])))\n"
    "    cfgs.append(('Combos', 'b+ALL_graph',\n"
    "                 dict(count_mode='binary', graph_metrics=list(ALL_GRAPH_METRICS))))\n"
    "    for combo in [\n"
    "        ['branching', 'diameter'],\n"
    "        ['branching', 'effective_graph_resistance'],\n"
    "        ['branching', 'mean_eccentricity', 'diameter'],\n"
    "        ['branching', 'mean_eccentricity', 'effective_graph_resistance'],\n"
    "        ['branching', 'mean_eccentricity', 'diameter', 'radius', 'effective_graph_resistance'],\n"
    "    ]:\n"
    "        cfgs.append(('Combos', 'b+[' + ','.join(combo) + ']',\n"
    "                     dict(count_mode='binary', graph_metrics=combo)))\n"
    "    for combo in [\n"
    "        ['branching'],\n"
    "        ['branching', 'diameter'],\n"
    "        ['branching', 'mean_eccentricity', 'effective_graph_resistance'],\n"
    "    ]:\n"
    "        cfgs.append(('Combos', 'tc+[' + ','.join(combo) + ']',\n"
    "                     dict(count_mode='total_component', graph_metrics=combo)))\n"
    "    for m in ALL_MOL_METRICS:\n"
    "        cfgs.append(('Mol metrics', f'b+mol:{m}',\n"
    "                     dict(count_mode='binary', molecule_metrics=[m])))\n"
    "    cfgs.append(('Combos', 'b+ALL_mol',\n"
    "                 dict(count_mode='binary', molecule_metrics=list(ALL_MOL_METRICS))))\n"
    "    cfgs.append(('Combos', 'tc+ALL_graph+ALL_mol',\n"
    "                 dict(count_mode='total_component',\n"
    "                      graph_metrics=list(ALL_GRAPH_METRICS),\n"
    "                      molecule_metrics=list(ALL_MOL_METRICS))))\n"
    "    return cfgs\n"
    "\n"
    "FP_CONFIGS = build_configs()\n"
    "print(f'{len(FP_CONFIGS)} fingerprint configurations defined')\n",
    "cell-fp-configs",
))

# ── Cell 5 ── compute fingerprints ───────────────────────────────────────────
cells.append(md("## 3. Compute all fingerprints & Tanimoto matrices\n", "md-compute"))

cells.append(code(
    "smiles_list  = list(MOL_DF['smiles'])\n"
    "labels_list  = list(MOL_DF['label'])\n"
    "branching_list = list(MOL_DF['branching'])\n"
    "chain_list   = list(MOL_DF['n_fluoro_c'])\n"
    "\n"
    "# ── Parse once ────────────────────────────────────────────────────────\n"
    "print('Parsing …')\n"
    "parsed = []\n"
    "for smi in smiles_list:\n"
    "    try:    parsed.append(parse_smiles(smi, halogens='F'))\n"
    "    except Exception as e:\n"
    "        print(f'  WARN {smi}: {e}'); parsed.append(None)\n"
    "print(f'{sum(p is not None for p in parsed)}/{len(parsed)} OK')\n"
    "\n"
    "# ── Tanimoto helpers ─────────────────────────────────────────────────\n"
    "def tanimoto(a, b):\n"
    "    a, b = np.asarray(a, float), np.asarray(b, float)\n"
    "    num = np.sum(np.minimum(a, b)); den = np.sum(np.maximum(a, b))\n"
    "    return float(num / den) if den > 0 else 0.0\n"
    "\n"
    "def pairwise_tanimoto(X):\n"
    "    n = len(X); S = np.eye(n)\n"
    "    for i in range(n):\n"
    "        for j in range(i+1, n):\n"
    "            v = tanimoto(X[i], X[j]); S[i,j] = S[j,i] = v\n"
    "    return S\n"
    "\n"
    "def off_diag(S):\n"
    "    mask = ~np.eye(len(S), dtype=bool)\n"
    "    v = S[mask]; return v.mean(), v.min(), v.std()\n"
    "\n"
    "# ── Compute all configs ───────────────────────────────────────────────\n"
    "def compute_config(parsed_list, section, label, kwargs):\n"
    "    vecs = []\n"
    "    for res in parsed_list:\n"
    "        if res is None: vecs.append(None); continue\n"
    "        try:\n"
    "            fp  = res.to_fingerprint(**kwargs)\n"
    "            arr = np.asarray(fp.fingerprints, dtype=float)\n"
    "            vecs.append(arr[0] if arr.ndim == 2 else arr)\n"
    "        except Exception: vecs.append(None)\n"
    "    valid = [v for v in vecs if v is not None]\n"
    "    if not valid: return None\n"
    "    w = max(len(v) for v in valid)\n"
    "    mat = np.zeros((len(vecs), w))\n"
    "    for i, v in enumerate(vecs):\n"
    "        if v is not None: mat[i, :len(v)] = v\n"
    "    S = pairwise_tanimoto(mat)\n"
    "    mean_t, min_t, std_t = off_diag(S)\n"
    "    return dict(section=section, label=label, kwargs=kwargs,\n"
    "                mat=mat, sim=S, mean_t=mean_t, min_t=min_t, std_t=std_t, n_cols=w)\n"
    "\n"
    "ALL_RESULTS = []\n"
    "for i, (sec, lbl, kw) in enumerate(FP_CONFIGS):\n"
    "    r = compute_config(parsed, sec, lbl, kw)\n"
    "    if r: ALL_RESULTS.append(r)\n"
    "    if (i+1) % 10 == 0: print(f'  {i+1}/{len(FP_CONFIGS)}')\n"
    "print(f'Done — {len(ALL_RESULTS)} configs computed')\n",
    "cell-compute",
))

# ── Cell 6 ── overall ranking ────────────────────────────────────────────────
cells.append(md("## 4. Overall discrimination ranking (all 39 configs)\n"
                "Lower mean Tanimoto = more discriminating.\n", "md-rank"))

cells.append(code(
    "_CMAP = LinearSegmentedColormap.from_list('T', ['#f7fbff','#c6dbef','#6baed6','#2171b5','#08306b'])\n"
    "\n"
    "def sec_colours(sections):\n"
    "    uniq = list(dict.fromkeys(sections))\n"
    "    cmap_ = plt.cm.tab10\n"
    "    return {s: cmap_(i/(max(len(uniq)-1,1))) for i,s in enumerate(uniq)}\n"
    "\n"
    "srt = sorted(ALL_RESULTS, key=lambda d: d['mean_t'])\n"
    "labs   = [d['label'] for d in srt]\n"
    "means  = [d['mean_t'] for d in srt]\n"
    "mins   = [d['min_t']  for d in srt]\n"
    "secs   = [d['section'] for d in srt]\n"
    "sc     = sec_colours(secs)\n"
    "cols   = [sc[s] for s in secs]\n"
    "\n"
    "n = len(labs)\n"
    "fig, ax = plt.subplots(figsize=(10, max(6, n*0.28)))\n"
    "y = np.arange(n)\n"
    "ax.barh(y, means, color=cols, alpha=0.85)\n"
    "ax.scatter(mins, y, marker='|', s=70, color='k', zorder=3)\n"
    "ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=7)\n"
    "ax.set_xlabel('Mean Tanimoto (lower = more discriminating)')\n"
    "ax.set_title(f'All {len(ALL_RESULTS)} FP configs — {len(MOL_DF)} molecules')\n"
    "ax.legend(handles=[mpatches.Patch(color=c,label=s) for s,c in sc.items()],\n"
    "          title='Section', loc='lower right', fontsize=8)\n"
    "for k in range(min(5,n)):\n"
    "    ax.annotate(f'#{k+1}', xy=(means[k],k), xytext=(means[k]+0.003,k),\n"
    "                fontsize=7, color='#d73027', fontweight='bold')\n"
    "ax.grid(axis='x', alpha=0.3); ax.invert_yaxis()\n"
    "plt.tight_layout(); plt.show()\n"
    "\n"
    "print('\\nTop-10:')\n"
    "print(pd.DataFrame([dict(rank=i+1, label=d['label'], section=d['section'],\n"
    "                         mean_T=round(d['mean_t'],4), min_T=round(d['min_t'],4),\n"
    "                         n_cols=d['n_cols'])\n"
    "                    for i,d in enumerate(srt[:10])]).to_string(index=False))\n",
    "cell-rank",
))

# ── Cell 7 ── chain-length sensitivity ──────────────────────────────────────
cells.append(md("## 5. Chain-length sensitivity\n"
                "Mean Tanimoto between **adjacent homologues** in the linear PFCA series (Cn ↔ Cn+1).\n"
                "Lower = fingerprint better distinguishes consecutive chain lengths.\n", "md-chain"))

cells.append(code(
    "linear_idx = [\n"
    "    MOL_DF.index[MOL_DF['label'] == f'C{n}-linear'][0] for n in range(2,9)\n"
    "]\n"
    "linear_labels = [MOL_DF.loc[i,'label'] for i in linear_idx]\n"
    "\n"
    "chain_rows = []\n"
    "for d in ALL_RESULTS:\n"
    "    S = d['sim']\n"
    "    adj_ts = [S[linear_idx[k], linear_idx[k+1]] for k in range(len(linear_idx)-1)]\n"
    "    chain_rows.append(dict(label=d['label'], section=d['section'],\n"
    "                           adj_mean_T=np.mean(adj_ts), adj_min_T=np.min(adj_ts)))\n"
    "chain_df = pd.DataFrame(chain_rows).sort_values('adj_mean_T')\n"
    "\n"
    "n=len(chain_df); fig,ax=plt.subplots(figsize=(10,max(5,n*0.28)))\n"
    "y=np.arange(n); sc2=sec_colours(list(chain_df['section'])); cols2=[sc2[s] for s in chain_df['section']]\n"
    "ax.barh(y,chain_df['adj_mean_T'],color=cols2,alpha=0.85)\n"
    "ax.scatter(chain_df['adj_min_T'],y,marker='|',s=70,color='k',zorder=3)\n"
    "ax.set_yticks(y); ax.set_yticklabels(chain_df['label'],fontsize=7)\n"
    "ax.set_xlabel('Mean Tanimoto between Cn ↔ Cn+1 (lower = more chain-length sensitive)')\n"
    "ax.set_title('Chain-length sensitivity — linear PFCA C2–C8')\n"
    "ax.legend(handles=[mpatches.Patch(color=c,label=s) for s,c in sc2.items()],\n"
    "          title='Section',loc='lower right',fontsize=8)\n"
    "ax.grid(axis='x',alpha=0.3); ax.invert_yaxis(); plt.tight_layout(); plt.show()\n"
    "\n"
    "# Show top-5 + binary heatmap for the linear series\n"
    "top5_lbl = [d['label'] for d in srt[:5]]\n"
    "show_lbl  = list(dict.fromkeys(['binary'] + top5_lbl))[:6]\n"
    "show_data = sorted([d for d in ALL_RESULTS if d['label'] in show_lbl],\n"
    "                   key=lambda d: show_lbl.index(d['label']))\n"
    "\n"
    "nc=len(show_data); fig,axes=plt.subplots(1,nc,figsize=(4*nc,4))\n"
    "axes=[axes] if nc==1 else list(axes)\n"
    "for ax,d in zip(axes,show_data):\n"
    "    sub=d['sim'][np.ix_(linear_idx,linear_idx)]\n"
    "    im=ax.imshow(sub,cmap=_CMAP,vmin=0,vmax=1)\n"
    "    ax.set_xticks(range(len(linear_labels))); ax.set_xticklabels(linear_labels,rotation=90,fontsize=8)\n"
    "    ax.set_yticks(range(len(linear_labels))); ax.set_yticklabels(linear_labels,fontsize=8)\n"
    "    cadj=chain_df.set_index('label').loc[d['label'],'adj_mean_T'] if d['label'] in chain_df['label'].values else float('nan')\n"
    "    ax.set_title(f'{d[\"label\"]}\\nadj-T={cadj:.3f}',fontsize=8)\n"
    "    plt.colorbar(im,ax=ax,fraction=0.046)\n"
    "fig.suptitle('Linear PFCA series — binary + top-5 configs',fontsize=10)\n"
    "plt.tight_layout(); plt.show()\n",
    "cell-chain",
))

# ── Cell 8 ── branching discrimination ──────────────────────────────────────
cells.append(md("## 6. Branching discrimination\n"
                "Mean Tanimoto between **linear and branched isomers** (same C-count, n=5–8).\n"
                "Lower = fingerprint better distinguishes branching.\n", "md-branch"))

cells.append(code(
    "BTYPES = ['linear','alpha','beta','gamma','gem','double_ab']\n"
    "\n"
    "def get_idx(n, btype):\n"
    "    rows = MOL_DF.index[(MOL_DF['n_fluoro_c']==n)&(MOL_DF['branching']==btype)]\n"
    "    return rows[0] if len(rows) else None\n"
    "\n"
    "branch_rows = []\n"
    "for d in ALL_RESULTS:\n"
    "    S = d['sim']; pairs = []\n"
    "    for n in range(5,9):\n"
    "        li = get_idx(n,'linear')\n"
    "        if li is None: continue\n"
    "        for bt in BTYPES[1:]:\n"
    "            bi = get_idx(n,bt)\n"
    "            if bi is not None: pairs.append(S[li,bi])\n"
    "    branch_rows.append(dict(label=d['label'],section=d['section'],\n"
    "                            branch_mean_T=np.mean(pairs) if pairs else np.nan,\n"
    "                            branch_min_T=np.min(pairs) if pairs else np.nan))\n"
    "branch_df = pd.DataFrame(branch_rows).dropna().sort_values('branch_mean_T')\n"
    "\n"
    "n=len(branch_df); fig,ax=plt.subplots(figsize=(10,max(5,n*0.28)))\n"
    "y=np.arange(n); sc3=sec_colours(list(branch_df['section'])); cols3=[sc3[s] for s in branch_df['section']]\n"
    "ax.barh(y,branch_df['branch_mean_T'],color=cols3,alpha=0.85)\n"
    "ax.scatter(branch_df['branch_min_T'],y,marker='|',s=70,color='k',zorder=3)\n"
    "ax.set_yticks(y); ax.set_yticklabels(branch_df['label'],fontsize=7)\n"
    "ax.set_xlabel('Mean T(linear, branched isomer) — lower = better branching discrimination')\n"
    "ax.set_title('Branching discrimination — linear vs α/β/γ/gem/α+β, C5–C8')\n"
    "ax.legend(handles=[mpatches.Patch(color=c,label=s) for s,c in sc3.items()],\n"
    "          title='Section',loc='lower right',fontsize=8)\n"
    "ax.grid(axis='x',alpha=0.3); ax.invert_yaxis(); plt.tight_layout(); plt.show()\n"
    "\n"
    "# C6 branching × branching matrix for binary + top-5\n"
    "N_DEMO=6\n"
    "demo_bt=[b for b in BTYPES if get_idx(N_DEMO,b) is not None]\n"
    "demo_idx=[get_idx(N_DEMO,b) for b in demo_bt]\n"
    "demo_lbl=[f'C{N_DEMO}-{b}' for b in demo_bt]\n"
    "\n"
    "nc=len(show_data); fig,axes=plt.subplots(1,nc,figsize=(4*nc,4))\n"
    "axes=[axes] if nc==1 else list(axes)\n"
    "for ax,d in zip(axes,show_data):\n"
    "    sub=d['sim'][np.ix_(demo_idx,demo_idx)]\n"
    "    im=ax.imshow(sub,cmap=_CMAP,vmin=0,vmax=1)\n"
    "    ax.set_xticks(range(len(demo_lbl))); ax.set_xticklabels(demo_lbl,rotation=90,fontsize=8)\n"
    "    ax.set_yticks(range(len(demo_lbl))); ax.set_yticklabels(demo_lbl,fontsize=8)\n"
    "    ax.set_title(f'{d[\"label\"]}',fontsize=8)\n"
    "    plt.colorbar(im,ax=ax,fraction=0.046)\n"
    "    for i in range(len(demo_lbl)):\n"
    "        for j in range(len(demo_lbl)):\n"
    "            ax.text(j,i,f'{sub[i,j]:.2f}',ha='center',va='center',fontsize=6.5,\n"
    "                    color='white' if sub[i,j]>0.6 else 'black')\n"
    "fig.suptitle(f'C{N_DEMO} branching×branching Tanimoto — binary + top-5',fontsize=10)\n"
    "plt.tight_layout(); plt.show()\n"
    "\n"
    "print('\\nTop-10 branching discriminators:')\n"
    "print(branch_df.head(10)[['label','branch_mean_T','branch_min_T']].to_string(index=False))\n",
    "cell-branch",
))

# ── Cell 9 ── MDS ────────────────────────────────────────────────────────────
cells.append(md("## 7. MDS projections\n"
                "Left: coloured by chain length.  Right: coloured by branching type.\n"
                "Shown for binary + top-5 most discriminating configs.\n", "md-mds"))

cells.append(code(
    "def classical_mds(dist, k=2):\n"
    "    n=dist.shape[0]; D2=dist**2\n"
    "    H=np.eye(n)-np.ones((n,n))/n; B=-0.5*H@D2@H\n"
    "    vals,vecs=np.linalg.eigh(B); idx=np.argsort(vals)[::-1]\n"
    "    vals,vecs=vals[idx],vecs[:,idx]\n"
    "    return vecs[:,:k]*np.sqrt(np.maximum(vals[:k],0))\n"
    "\n"
    "CHAIN_CMAP = plt.cm.viridis\n"
    "BRANCH_COL = {'linear':'#2166ac','alpha':'#4dac26','beta':'#e08214',\n"
    "              'gamma':'#762a83','gem':'#d73027','double_ab':'#1b9e77'}\n"
    "cmin,cmax = min(chain_list),max(chain_list)\n"
    "chain_colours  = [CHAIN_CMAP((c-cmin)/(cmax-cmin+1e-9)) for c in chain_list]\n"
    "branch_colours = [BRANCH_COL.get(b,'#888') for b in branching_list]\n"
    "\n"
    "nc=len(show_data); fig,axes=plt.subplots(nc,2,figsize=(10,4.5*nc))\n"
    "if nc==1: axes=[axes]\n"
    "for row,d in enumerate(show_data):\n"
    "    dist=np.clip(1-d['sim'],0,None); np.fill_diagonal(dist,0)\n"
    "    try: coords=classical_mds(dist)\n"
    "    except: [axes[row][k].axis('off') for k in range(2)]; continue\n"
    "    for col,(clrs,title_suf) in enumerate([(chain_colours,'chain-length'),(branch_colours,'branching')]):\n"
    "        ax=axes[row][col]\n"
    "        ax.scatter(coords[:,0],coords[:,1],c=clrs,s=50,zorder=3,\n"
    "                   edgecolors='white',linewidths=0.4)\n"
    "        for j,lbl in enumerate(labels_list):\n"
    "            ax.annotate(lbl,coords[j],xytext=(2,2),textcoords='offset points',\n"
    "                        fontsize=5,zorder=4)\n"
    "        ax.set_title(f'{d[\"label\"]} — {title_suf}',fontsize=7)\n"
    "        ax.axis('equal'); ax.grid(alpha=0.2)\n"
    "\n"
    "sm=plt.cm.ScalarMappable(cmap=CHAIN_CMAP,norm=plt.Normalize(cmin,cmax))\n"
    "sm.set_array([]); fig.colorbar(sm,ax=np.array(axes)[:,0],label='n_fluoro_c',shrink=0.6)\n"
    "fig.legend(handles=[mpatches.Patch(color=c,label=b) for b,c in BRANCH_COL.items()],\n"
    "           title='Branching',loc='lower right',fontsize=8,ncol=3)\n"
    "fig.suptitle('MDS — binary + top-5 configs',fontsize=10)\n"
    "plt.tight_layout(); plt.show()\n",
    "cell-mds",
))

# ── Cell 10 ── per-metric impact ─────────────────────────────────────────────
cells.append(md("## 8. Per-metric impact (vs. binary baseline)\n"
                "Δmean_T = binary_mean_T − config_mean_T.  Positive = improvement.\n"
                "Three panels: overall, chain-length, branching axes.\n", "md-impact"))

cells.append(code(
    "def vec(df, lbl, col):\n"
    "    row=df[df['label']==lbl]; return float(row[col].iloc[0]) if len(row) else float('nan')\n"
    "\n"
    "binary_d = next(d for d in ALL_RESULTS if d['label']=='binary')\n"
    "B_mean   = binary_d['mean_t']\n"
    "B_chain  = vec(chain_df,  'binary','adj_mean_T')\n"
    "B_branch = vec(branch_df, 'binary','branch_mean_T') if 'binary' in branch_df['label'].values else float('nan')\n"
    "\n"
    "def plot_impact(rows, title):\n"
    "    labs_ = [d['label'].replace('b+mol:','').replace('b+','') for d in rows]\n"
    "    Δo = [B_mean  - d['mean_t']                for d in rows]\n"
    "    Δc = [B_chain - vec(chain_df,d['label'],'adj_mean_T') for d in rows]\n"
    "    Δb = [B_branch- vec(branch_df,d['label'],'branch_mean_T') for d in rows]\n"
    "    fig,axes=plt.subplots(1,3,figsize=(15,max(4,len(rows)*0.36)),sharey=True)\n"
    "    for ax,deltas,sub in zip(axes,[Δo,Δc,Δb],\n"
    "            ['Overall\\n(all pairs)','Chain-length\\n(adj homologues)',\n"
    "             'Branching\\n(lin vs branched)']):\n"
    "        y_=np.arange(len(labs_))\n"
    "        clrs=['#2166ac' if v>0 else '#d73027' for v in deltas]\n"
    "        ax.barh(y_,deltas,color=clrs,alpha=0.85)\n"
    "        ax.set_yticks(y_); ax.set_yticklabels(labs_,fontsize=8)\n"
    "        ax.axvline(0,color='k',lw=1)\n"
    "        ax.set_xlabel('Δ mean T vs binary (positive=improvement)')\n"
    "        ax.set_title(sub); ax.grid(axis='x',alpha=0.3); ax.invert_yaxis()\n"
    "    fig.suptitle(title,fontsize=10); plt.tight_layout(); plt.show()\n"
    "\n"
    "graph_rows=[d for d in ALL_RESULTS if d['section']=='Graph metrics']\n"
    "mol_rows  =[d for d in ALL_RESULTS if d['section']=='Mol metrics']\n"
    "plot_impact(sorted(graph_rows,key=lambda d:d['mean_t']),\n"
    "            'Individual graph metrics — impact decomposition')\n"
    "plot_impact(sorted(mol_rows,key=lambda d:d['mean_t']),\n"
    "            'Individual mol metrics — impact decomposition')\n",
    "cell-impact",
))

# ── Cell 11 ── heatmaps ──────────────────────────────────────────────────────
cells.append(md("## 9. Tanimoto heatmaps (binary + top-10)\n\n"
                "Rows/columns ordered by hierarchical clustering.\n", "md-heatmaps"))

cells.append(code(
    "def cluster_order(sim):\n"
    "    dist=np.clip(1-sim,0,None); np.fill_diagonal(dist,0)\n"
    "    cond=squareform(dist,checks=False); cond=np.clip(cond,0,None)\n"
    "    try:    return leaves_list(linkage(cond,method='average'))\n"
    "    except: return np.arange(len(sim))\n"
    "\n"
    "top10=[next(r for r in ALL_RESULTS if r['label']=='binary')]+ srt[:10]\n"
    "ncols=3; nrows=int(np.ceil(len(top10)/ncols))\n"
    "fig,axes=plt.subplots(nrows,ncols,figsize=(5.5*ncols,5*nrows))\n"
    "axes_flat=np.array(axes).flatten()\n"
    "for i,d in enumerate(top10):\n"
    "    ax=axes_flat[i]\n"
    "    order=cluster_order(d['sim'])\n"
    "    s=d['sim'][np.ix_(order,order)]\n"
    "    lbl=[labels_list[k] for k in order]\n"
    "    im=ax.imshow(s,cmap=_CMAP,vmin=0,vmax=1,aspect='auto')\n"
    "    ax.set_xticks(range(len(lbl))); ax.set_xticklabels(lbl,rotation=90,fontsize=4.5)\n"
    "    ax.set_yticks(range(len(lbl))); ax.set_yticklabels(lbl,fontsize=4.5)\n"
    "    ax.set_title(f'{d[\"label\"]}\\nmean T={d[\"mean_t\"]:.3f}',fontsize=7)\n"
    "    plt.colorbar(im,ax=ax,fraction=0.046,pad=0.04)\n"
    "for ax in axes_flat[len(top10):]:\n"
    "    ax.axis('off')\n"
    "fig.suptitle(f'Tanimoto heatmaps — binary + top-10 configs\\n{len(MOL_DF)} molecules',fontsize=10)\n"
    "plt.tight_layout(); plt.show()\n",
    "cell-heatmaps",
))

# ── Cell 12 ── summary ───────────────────────────────────────────────────────
cells.append(md("## 10. Summary\n", "md-summary"))

cells.append(code(
    "summary = pd.DataFrame([dict(\n"
    "    rank=i+1, label=d['label'], section=d['section'],\n"
    "    mean_T=round(d['mean_t'],4), min_T=round(d['min_t'],4), n_cols=d['n_cols'],\n"
    "    chain_adj_T=round(vec(chain_df,d['label'],'adj_mean_T'),4),\n"
    "    branch_T   =round(vec(branch_df,d['label'],'branch_mean_T'),4),\n"
    ") for i,d in enumerate(srt)])\n"
    "\n"
    "display(summary.style\n"
    "    .background_gradient(subset=['mean_T'],      cmap='RdYlGn_r', vmin=0, vmax=1)\n"
    "    .background_gradient(subset=['chain_adj_T'], cmap='RdYlGn_r', vmin=0, vmax=1)\n"
    "    .background_gradient(subset=['branch_T'],    cmap='RdYlGn_r', vmin=0, vmax=1)\n"
    "    .format(precision=4)\n"
    ")\n"
    "\n"
    "out_csv = NOTEBOOK_DIR / 'results' / 'tanimoto' / 'structure_analysis_summary.csv'\n"
    "out_csv.parent.mkdir(parents=True, exist_ok=True)\n"
    "summary.to_csv(out_csv, index=False)\n"
    "print(f'Saved → {out_csv}')\n"
    "\n"
    "# Key findings\n"
    "best_all    = summary.nsmallest(1,'mean_T').iloc[0]['label']\n"
    "best_chain  = summary.nsmallest(1,'chain_adj_T').iloc[0]['label']\n"
    "best_branch = summary.nsmallest(1,'branch_T').iloc[0]['label']\n"
    "print(f'\\nBest overall:   {best_all}')\n"
    "print(f'Best chain-len: {best_chain}')\n"
    "print(f'Best branching: {best_branch}')\n"
    "tradeoff = not (best_all==best_chain==best_branch)\n"
    "print('Trade-off detected.' if tradeoff else 'No trade-off: one config wins all axes.')\n",
    "cell-summary",
))

# ── Write notebook ────────────────────────────────────────────────────────────
nb = {"nbformat": 4, "nbformat_minor": 5, "metadata": {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python"}
}, "cells": cells}

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
print(f"Written: {OUT}  ({OUT.stat().st_size:,} bytes, {len(cells)} cells)")
