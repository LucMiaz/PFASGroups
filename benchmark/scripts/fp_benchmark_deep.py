"""
Fingerprint Deep Analysis – Stratified Pairwise Similarity
===========================================================

Exhaustive benchmark comparing seven PFASGroups fingerprint variants and
(optionally) the CSRML / TxP_PFAS fingerprint on three compound datasets:

Dataset A – F-only PFAS (600 compounds)
    6 FG types × 10 chain sizes × 10 k-variants
    (PFCA, PFSA, PFAL, PFAM, Ketone, PFPE)

Dataset B – Cl-only analogs (75 compounds)
    Same 3 FG types as above but F→Cl (PFCA, PFSA, PFAL)
    3 FG types × 5 n × 5 k

Dataset C – Mixed F+Cl PFAS (90 compounds)
    One F-chain + one Cl-chain with same functional group
    3 FG types × 5 n_F × 3 n_Cl × 2 k

Fingerprint variants
--------------------
1. binary_F_per         : binary, F only, saturation='per'       (116d)
2. count_F_per          : count, F only, saturation='per'        (116d)
3. maxcomp_F_per        : max_component, F only, sat='per'       (116d)
4. binary_F_poly        : binary, F only, saturation='poly'      (116d)
5. binary_F_none        : binary, F only, saturation=None        (116d)
6. binary_FCl_per       : binary, F+Cl stacked, sat='per'        (232d)
7. descriptive          : 18-dim numerical descriptor vector

Pair stratification (Dataset A, 5 categories)
----------------------------------------------
A – Same FG, same n, diff k  (only count mode differs)
B – Same FG, diff n, same k  (only max_component mode differs)
C – Same FG, diff n, diff k  (both count & max-comp differ)
D – Diff FG, same n, same k  (all modes differ)
E – Diff FG, diff n, diff k  (all modes differ, max contrast)

For multi-halogen (Dataset B+C), category:
F – Diff halogen, same FG (F vs Cl compound, same n & k)

Statistical tests
-----------------
Kruskal–Wallis (non-parametric one-way test across all categories)
followed by pairwise Mann–Whitney U with Bonferroni correction.
Effect size: η² = (H − k + 1) / (n − k) for KW.

CSRML comparison (optional)
----------------------------
Run with --txppfas_csv pointing to the exported SI Table S2 CSV from
  Richard et al., Chem. Res. Toxicol. 2023, 36, 508-534
  https://pubs.acs.org/doi/suppl/10.1021/acs.chemrestox.2c00403
  (Sheet: 'Table S2', exported to CSV with DTXSID in index column)
A matching SMILES TSV must also be supplied via --smiles_tsv.

Outputs
-------
fp_deep_split_violin.pdf / .png   – main split plot
fp_deep_heatmap_summary.pdf / .png  – mean Jaccard heatmap
fp_deep_significance.pdf / .png   – KW + MW post-hoc significance heatmap
fp_deep_multihalogen.pdf / .png   – multi-halogen F vs Cl analysis
fp_deep_pca_compare.pdf / .png    – PCA overlay across FP variants
fp_deep_results.csv               – full pairwise result table
fp_deep_stats.csv                 – statistical test summary

Usage
-----
python fp_benchmark_deep.py [--outdir PATH] [--article_imgs PATH] [--quiet]
python fp_benchmark_deep.py --txppfas_csv Richard2023_SI_TableS2.csv --smiles_tsv molecules.tsv
"""
from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from scipy.stats import kruskal, mannwhitneyu
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib.ticker as mticker

# ── Path setup ────────────────────────────────────────────────────────────────
SCRIPT_DIR    = Path(__file__).resolve().parent
REPO_ROOT     = SCRIPT_DIR.parent.parent
BENCHMARK_DIR = SCRIPT_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

_ARTICLE_IMGS = REPO_ROOT.parent / 'overleaf_PFASgroups_article' / 'imgs'
warnings.filterwarnings('ignore')

# ── Compound design ───────────────────────────────────────────────────────────
FG_FULL = [
    ("PFCA",   "Carboxylic acid",  "C(=O)O"),
    ("PFSA",   "Sulfonic acid",    "S(=O)(=O)O"),
    ("PFAL",   "Alcohol",          "O"),
    ("PFAM",   "Sulfonamide",      "S(=O)(=O)N"),
    ("Ketone", "Ketone",           "C(=O)C"),
    ("PFPE",   "Ether",            "OC"),
]
FG_MIXED = [("PFCA", "Carboxylic acid", "C(=O)O"),
            ("PFSA", "Sulfonic acid", "S(=O)(=O)O"),
            ("PFAL", "Alcohol", "O")]

COMPONENT_SIZES_A = list(range(2, 12))   # n_cf2 ∈ {2..11}, 10 sizes
K_VARIANTS_A      = list(range(1, 11))   # k ∈ {1..10}
COMPONENT_SIZES_B = [2, 4, 6, 8, 10]
K_VARIANTS_B      = [1, 2, 4, 6, 8]
N_F_MIXED  = [2, 4, 6, 8, 10]
N_CL_MIXED = [2, 4, 6]
K_MIXED    = [1, 2]

FG_COLOURS = ['#1976D2', '#D32F2F', '#388E3C', '#F57C00', '#7B1FA2', '#0288D1']
FG_ID_TO_COLOUR = {fid: FG_COLOURS[i] for i, (fid, *_) in enumerate(FG_FULL)}
FP_VARIANT_COLOURS = {
    'binary_F_per'  : '#1976D2',
    'count_F_per'   : '#D32F2F',
    'maxcomp_F_per' : '#388E3C',
    'binary_F_poly' : '#F57C00',
    'binary_F_none' : '#7B1FA2',
    'binary_FCl_per': '#0288D1',
    'descriptive'   : '#795548',
}

# ── SMILES builders ───────────────────────────────────────────────────────────

def _pf_chain(n_cf2: int) -> str: return 'FC(F)(F)' + 'C(F)(F)' * n_cf2
def _pc_chain(n_cf2: int) -> str: return 'ClC(Cl)(Cl)' + 'C(Cl)(Cl)' * n_cf2

def build_f(fg_suffix: str, n: int, k: int) -> str:
    return '.'.join([_pf_chain(n) + fg_suffix] * k)

def build_cl(fg_suffix: str, n: int, k: int) -> str:
    return '.'.join([_pc_chain(n) + fg_suffix] * k)

def build_mixed(fg_suffix: str, n_f: int, n_cl: int, k: int) -> str:
    f_frag  = _pf_chain(n_f)  + fg_suffix
    cl_frag = _pc_chain(n_cl) + fg_suffix
    return '.'.join([f_frag] * k + [cl_frag] * k)


def validate_smiles(smi: str) -> bool:
    from rdkit import Chem
    return Chem.MolFromSmiles(smi) is not None


# ── Dataset construction ──────────────────────────────────────────────────────

def build_dataset(verbose: bool = True) -> pd.DataFrame:
    """Return a DataFrame with all 765+ compounds (meta + SMILES)."""
    rows = []

    # Dataset A – F-only
    for fid, fname, fg_suffix in FG_FULL:
        for n in COMPONENT_SIZES_A:
            for k in K_VARIANTS_A:
                smi = build_f(fg_suffix, n, k)
                rows.append({'dataset': 'A_F',
                             'fg_id': fid, 'fg_name': fname,
                             'n': n, 'k': k, 'halogen': 'F',
                             'smiles': smi, 'valid': validate_smiles(smi)})

    # Dataset B – Cl-only
    for fid, fname, fg_suffix in FG_MIXED:
        for n in COMPONENT_SIZES_B:
            for k in K_VARIANTS_B:
                smi = build_cl(fg_suffix, n, k)
                rows.append({'dataset': 'B_Cl',
                             'fg_id': fid, 'fg_name': fname,
                             'n': n, 'k': k, 'halogen': 'Cl',
                             'smiles': smi, 'valid': validate_smiles(smi)})

    # Dataset C – mixed F+Cl
    for fid, fname, fg_suffix in FG_MIXED:
        for n_f in N_F_MIXED:
            for n_cl in N_CL_MIXED:
                for k in K_MIXED:
                    smi = build_mixed(fg_suffix, n_f, n_cl, k)
                    rows.append({'dataset': 'C_mixed',
                                 'fg_id': fid, 'fg_name': fname,
                                 'n': n_f, 'k': k, 'halogen': 'F+Cl',
                                 'n_cl': n_cl, 'smiles': smi,
                                 'valid': validate_smiles(smi)})

    df = pd.DataFrame(rows).fillna(0)
    n_valid = df['valid'].sum()
    if verbose:
        print(f"  Built dataset: {len(df)} total ({n_valid} valid SMILES)")
        for ds in ('A_F', 'B_Cl', 'C_mixed'):
            n = (df['dataset'] == ds).sum()
            print(f"    {ds}: {n} compounds")
    return df[df['valid']].reset_index(drop=True)


# ── Fingerprint computation ───────────────────────────────────────────────────

def _dims_and_names(halogens, saturation='per') -> tuple[int, list]:
    from HalogenGroups import generate_fingerprint
    fp, info = generate_fingerprint(
        'FC(F)(F)C(=O)O', representation='vector', count_mode='binary',
        halogens=halogens, saturation=saturation)
    return len(fp), info['group_names']


def compute_fp_variant(smiles_list: list, halogen,
                       saturation, count_mode, verbose=True) -> np.ndarray:
    from HalogenGroups import generate_fingerprint
    if isinstance(halogen, str):
        halogen_arg = halogen
    else:
        halogen_arg = list(halogen)

    n_dim, _ = _dims_and_names(halogen_arg, saturation)
    n = len(smiles_list)
    X = np.zeros((n, n_dim), dtype=np.float32)
    for i, smi in enumerate(smiles_list):
        try:
            fp, _ = generate_fingerprint(
                smi, representation='vector', count_mode=count_mode,
                halogens=halogen_arg, saturation=saturation)
            X[i] = fp
        except Exception:
            pass
    return X


def compute_descriptive(smiles_list: list, X_bin: np.ndarray,
                         X_cnt: np.ndarray, X_mx: np.ndarray) -> np.ndarray:
    """18-d numerical descriptor derived from RDKit + FP vectors."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Descriptors

    n = len(smiles_list)
    feats = np.zeros((n, 18), dtype=np.float32)

    # Group category boundaries in the 116-bit FP
    # Based on default group ordering: OECD 0-27 (28), Generic 28-72 (45), Telomer 73-115 (43)
    OE, GE, TE = 28, 73, 116

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        n_H = mol.GetNumHeavyAtoms()
        n_F = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
        n_Cl= sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 17)
        n_C = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        mw  = Descriptors.MolWt(mol)
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        n_frag  = smi.count('.') + 1   # dot-separated components

        bin_i = X_bin[i]
        cnt_i = X_cnt[i]
        mx_i  = X_mx[i]
        mx_nz = mx_i[mx_i > 0]

        feats[i] = [
            n_H,                                    # 0 heavy atoms
            n_F,                                    # 1 F count
            n_Cl,                                   # 2 Cl count
            n_C,                                    # 3 C count
            n_F / n_H if n_H > 0 else 0,            # 4 F ratio
            n_F / n_C if n_C > 0 else 0,            # 5 F/C ratio
            mw,                                     # 6 mol weight
            n_rings,                                # 7 ring count
            n_frag,                                 # 8 fragment count
            float(bin_i.sum()),                     # 9 total active bits
            float(cnt_i.sum()),                     # 10 total count
            float(mx_i.max()) if len(mx_i) > 0 else 0,  # 11 max component
            float(mx_nz.mean()) if len(mx_nz) > 0 else 0,  # 12 mean active comp
            float(bin_i[:OE].sum()),                # 13 n OECD bits active
            float(bin_i[OE:GE].sum()),              # 14 n generic bits active
            float(bin_i[GE:TE].sum()),              # 15 n telomer bits active
            float(mx_i[:OE].max()),                 # 16 max OECD component
            float(mx_i[OE:GE].max()),               # 17 max generic component
        ]
    return feats


# ── Similarity / distance ─────────────────────────────────────────────────────

def gen_jaccard(a: np.ndarray, b: np.ndarray) -> float:
    """Generalised Tanimoto–Jaccard: Σmin(aᵢ,bᵢ)/Σmax(aᵢ,bᵢ)."""
    num = float(np.sum(np.minimum(a, b)))
    den = float(np.sum(np.maximum(a, b)))
    return num / den if den > 0 else 0.0


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else 0.0


# ── Pair stratification ────────────────────────────────────────────────────────

CAT_LABELS = {
    'A': 'Same FG\nSame n\nDiff k',
    'B': 'Same FG\nDiff n\nSame k',
    'C': 'Same FG\nDiff n\nDiff k',
    'D': 'Diff FG\nSame n\nSame k',
    'E': 'Diff FG\nDiff n\nDiff k',
}
CAT_COLOURS = {
    'A': '#1976D2', 'B': '#388E3C', 'C': '#F57C00',
    'D': '#D32F2F', 'E': '#7B1FA2',
}

def categorise_pair(r1, r2) -> str | None:
    """Return category string for a pair (Dataset A only)."""
    same_fg = r1['fg_id'] == r2['fg_id']
    same_n  = r1['n']     == r2['n']
    same_k  = r1['k']     == r2['k']
    if same_fg and same_n and same_k:
        return None   # trivially identical
    if same_fg and same_n and not same_k:
        return 'A'
    if same_fg and not same_n and same_k:
        return 'B'
    if same_fg and not same_n and not same_k:
        return 'C'
    if not same_fg and same_n and same_k:
        return 'D'
    if not same_fg and not same_n and not same_k:
        return 'E'
    return None   # mixed (Diff FG, one same attribute) – skip


def sample_pairs(df_A: pd.DataFrame, max_pairs: int = 6000,
                 seed: int = 42) -> list[tuple[int, int, str]]:
    """Return a list of (i, j, category) pairs sampled from Dataset A."""
    rng = np.random.default_rng(seed)
    n   = len(df_A)
    idx = df_A.index.to_numpy()
    records = df_A.to_dict('records')

    # Build category index for fast lookup
    cat_pairs: dict[str, list[tuple[int, int]]] = {c: [] for c in CAT_LABELS}

    # Sample random pairs until we have enough in each category
    target_per_cat = max_pairs // len(CAT_LABELS)
    attempts = 0
    while min(len(v) for v in cat_pairs.values()) < target_per_cat and attempts < max_pairs * 20:
        attempts += 1
        ii, jj = rng.integers(0, n, size=2)
        if ii >= jj:
            continue
        cat = categorise_pair(records[ii], records[jj])
        if cat is not None and len(cat_pairs[cat]) < target_per_cat * 3:
            cat_pairs[cat].append((idx[ii], idx[jj]))

    pairs = []
    for cat, ps in cat_pairs.items():
        rng.shuffle(ps)
        for i, j in ps[:target_per_cat]:
            pairs.append((i, j, cat))
    return pairs


# ── Utility ───────────────────────────────────────────────────────────────────

def _save(fig, outdir: Path, stem: str, verbose: bool = True) -> None:
    for ext in ('.pdf', '.png'):
        fig.savefig(outdir / (stem + ext), dpi=150, bbox_inches='tight')
    if verbose:
        print(f"  Saved {stem}.pdf / .png")
    plt.close(fig)


def _numpy_pca(X: np.ndarray, n_components: int = 2):
    Xc = X - X.mean(axis=0)
    _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    variance  = s ** 2 / max(X.shape[0] - 1, 1)
    total_var = variance.sum()
    ev_ratio  = variance[:n_components] / total_var if total_var > 0 else np.zeros(n_components)
    return Xc @ Vt[:n_components].T, ev_ratio


def _zscore(X: np.ndarray) -> np.ndarray:
    mu, sigma = X.mean(axis=0), X.std(axis=0)
    sigma[sigma < 1e-9] = 1.0
    return (X - mu) / sigma


# ── Statistical tests ─────────────────────────────────────────────────────────

def run_stats(cat_data: dict[str, list[float]]) -> dict:
    """KW test + pairwise MW U with Bonferroni. Returns summary dict."""
    cats  = [c for c in CAT_LABELS if cat_data.get(c)]
    arrs  = [np.array(cat_data[c]) for c in cats]
    result = {'categories': cats, 'n': {c: len(cat_data[c]) for c in cats}}

    # KW
    if len(arrs) >= 2 and all(len(a) > 0 for a in arrs):
        H, p_kw = kruskal(*arrs)
        n_total = sum(len(a) for a in arrs)
        k       = len(arrs)
        eta2    = (H - k + 1) / (n_total - k) if n_total > k else 0.0
        result['kruskal_H'] = float(H)
        result['kruskal_p'] = float(p_kw)
        result['eta2']      = float(eta2)
    else:
        result['kruskal_H'] = float('nan')
        result['kruskal_p'] = float('nan')
        result['eta2']      = float('nan')

    # Pairwise MW U (Bonferroni-corrected)
    n_tests = len(cats) * (len(cats) - 1) // 2
    pairwise = {}
    for ci, cj in combinations(range(len(cats)), 2):
        ai, aj = arrs[ci], arrs[cj]
        if len(ai) > 0 and len(aj) > 0:
            _, p = mannwhitneyu(ai, aj, alternative='two-sided')
            p_bonf = min(p * n_tests, 1.0)
        else:
            p_bonf = float('nan')
        pairwise[f'{cats[ci]}_vs_{cats[cj]}'] = p_bonf
    result['pairwise_bonf'] = pairwise

    # Summary per category
    result['medians'] = {c: float(np.median(cat_data[c])) if cat_data.get(c) else float('nan')
                         for c in CAT_LABELS}
    result['means']   = {c: float(np.mean(cat_data[c])) if cat_data.get(c) else float('nan')
                         for c in CAT_LABELS}
    result['stds']    = {c: float(np.std(cat_data[c])) if cat_data.get(c) else float('nan')
                         for c in CAT_LABELS}
    return result


# ── Plotting ──────────────────────────────────────────────────────────────────

FP_DISPLAY_NAMES = {
    'binary_F_per'  : 'Binary F\n(sat=per)',
    'count_F_per'   : 'Count F\n(sat=per)',
    'maxcomp_F_per' : 'Max-comp F\n(sat=per)',
    'binary_F_poly' : 'Binary F\n(sat=poly)',
    'binary_F_none' : 'Binary F\n(sat=none)',
    'binary_FCl_per': 'Binary F+Cl\n(sat=per)',
    'descriptive'   : 'Descriptive\n(cosine)',
}


def plot_split_violin(
    all_cat_data: dict[str, dict[str, list[float]]],
    outdir: Path, verbose: bool = True,
) -> None:
    """
    Large multi-panel violin plot.
    Rows = FP variants. Columns = pair categories.
    """
    fp_names  = list(FP_DISPLAY_NAMES.keys())
    cat_names = list(CAT_LABELS.keys())
    n_fp, n_cat = len(fp_names), len(cat_names)

    fig, axes = plt.subplots(n_fp, 1, figsize=(14, 2.2 * n_fp), constrained_layout=True)
    if n_fp == 1:
        axes = [axes]

    for ax, fp_name in zip(axes, fp_names):
        data_per_cat = [all_cat_data[fp_name].get(c, []) for c in cat_names]
        non_empty    = [d for d in data_per_cat if len(d) > 0]
        if not non_empty:
            ax.set_visible(False)
            continue

        positions = list(range(n_cat))
        parts = ax.violinplot(
            data_per_cat,
            positions=positions,
            showmedians=True,
            showextrema=False,
            widths=0.7,
        )
        for pc in parts['bodies']:
            pc.set_facecolor(FP_VARIANT_COLOURS[fp_name])
            pc.set_alpha(0.65)
        parts['cmedians'].set_color('black')
        parts['cmedians'].set_linewidth(1.5)

        # Add box-plot quartiles as thin boxes
        for ci, cat in enumerate(cat_names):
            d = all_cat_data[fp_name].get(cat, [])
            if len(d) < 4:
                continue
            q1, med, q3 = np.percentile(d, [25, 50, 75])
            iqr = q3 - q1
            ax.add_patch(plt.Rectangle((ci - 0.12, q1), 0.24, iqr,
                                       facecolor='white', edgecolor='black',
                                       linewidth=0.8, zorder=3))
            ax.axhline(med, xmin=(ci)/n_cat + 0.01,
                       xmax=(ci + 1)/n_cat - 0.01,
                       color='black', linewidth=1.2, zorder=4, linestyle='-')

        # Significance stars between adjacent categories
        y_max = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1.0
        ax.set_ylim(-0.05, 1.1)
        ax.set_xticks(positions)
        ax.set_xticklabels([CAT_LABELS[c] for c in cat_names], fontsize=8)
        ax.set_ylabel('Jaccard\nsimilarity', fontsize=8)
        ax.set_title(FP_DISPLAY_NAMES[fp_name], fontsize=10, fontweight='bold',
                     loc='left', pad=2)
        ax.grid(axis='y', alpha=0.25, linewidth=0.5)

    fig.suptitle(
        'Pairwise Fingerprint Similarity by Structural Category\n'
        '(6 FG types × 10 chain sizes × 10 k-variants; sampled pairs)',
        fontsize=13, fontweight='bold',
    )
    _save(fig, outdir, 'fp_deep_split_violin', verbose=verbose)


def plot_heatmap_summary(
    all_cat_data: dict[str, dict[str, list[float]]],
    stats_dict: dict[str, dict],
    outdir: Path, verbose: bool = True,
) -> None:
    """Mean Jaccard heatmap (FP variants × categories) with η² side bar."""
    fp_names  = list(FP_DISPLAY_NAMES.keys())
    cat_names = list(CAT_LABELS.keys())
    n_fp, n_cat = len(fp_names), len(cat_names)

    mean_mat = np.zeros((n_fp, n_cat), dtype=np.float32)
    std_mat  = np.zeros((n_fp, n_cat), dtype=np.float32)
    for fi, fp in enumerate(fp_names):
        for ci, cat in enumerate(cat_names):
            d = all_cat_data[fp].get(cat, [])
            mean_mat[fi, ci] = np.mean(d) if d else np.nan
            std_mat[fi, ci]  = np.std(d)  if d else np.nan

    eta2_vals = np.array([stats_dict.get(fp, {}).get('eta2', np.nan) for fp in fp_names])

    fig, (ax, ax_eta) = plt.subplots(1, 2, figsize=(13, 0.9 * n_fp + 2),
                                      gridspec_kw={'width_ratios': [6, 1]})
    im = ax.imshow(mean_mat, cmap='Blues', vmin=0, vmax=1, aspect='auto')
    for fi in range(n_fp):
        for ci in range(n_cat):
            v = mean_mat[fi, ci]
            txt = f'{v:.2f}' if not np.isnan(v) else '–'
            ax.text(ci, fi, txt, ha='center', va='center',
                    fontsize=8.5,
                    color='white' if v > 0.6 else 'black')
    ax.set_xticks(range(n_cat))
    ax.set_xticklabels([CAT_LABELS[c] for c in cat_names], fontsize=9)
    ax.set_yticks(range(n_fp))
    ax.set_yticklabels([FP_DISPLAY_NAMES[fp].replace('\n', ' ') for fp in fp_names], fontsize=9)
    ax.set_title('Mean pairwise Jaccard similarity', fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label='Mean Jaccard')

    # η² bar chart
    eta2_vals_clean = np.nan_to_num(eta2_vals, nan=0)
    ax_eta.barh(range(n_fp), eta2_vals_clean,
                color=[FP_VARIANT_COLOURS[fp] for fp in fp_names], alpha=0.8)
    ax_eta.set_yticks(range(n_fp))
    ax_eta.set_yticklabels([])
    ax_eta.set_xlabel('η²\n(KW effect)', fontsize=9)
    ax_eta.set_xlim(0, max(1, float(np.nanmax(eta2_vals_clean)) * 1.1))
    ax_eta.axvline(0.06, linestyle='--', color='grey', alpha=0.5, linewidth=0.8, label='small')
    ax_eta.axvline(0.14, linestyle='--', color='orange', alpha=0.5, linewidth=0.8, label='medium')
    ax_eta.legend(fontsize=7, loc='lower right')
    ax_eta.set_title('Effect\nsize', fontsize=9)

    fig.suptitle('Mean Similarity Heatmap + Kruskal–Wallis Effect Size', fontsize=12, fontweight='bold')
    fig.tight_layout()
    _save(fig, outdir, 'fp_deep_heatmap_summary', verbose=verbose)


def plot_significance(
    stats_dict: dict[str, dict],
    outdir: Path, verbose: bool = True,
) -> None:
    """Post-hoc Bonferroni p-value heatmap for all FP types."""
    fp_names = list(FP_DISPLAY_NAMES.keys())
    cat_names = list(CAT_LABELS.keys())
    cat_pairs = [(ci, cj) for ci in range(len(cat_names)) for cj in range(ci + 1, len(cat_names))]
    pair_labels = [f'{cat_names[ci]}\nvs\n{cat_names[cj]}' for ci, cj in cat_pairs]
    n_fp, n_pairs = len(fp_names), len(cat_pairs)

    pval_mat = np.ones((n_fp, n_pairs))
    for fi, fp in enumerate(fp_names):
        pw = stats_dict.get(fp, {}).get('pairwise_bonf', {})
        for pi, (ci, cj) in enumerate(cat_pairs):
            key = f'{cat_names[ci]}_vs_{cat_names[cj]}'
            pval_mat[fi, pi] = pw.get(key, 1.0)

    # Significance coding: -log10(p)
    logp_mat = -np.log10(np.clip(pval_mat, 1e-15, 1.0))

    fig, ax = plt.subplots(figsize=(max(8, n_pairs * 1.2), n_fp * 0.75 + 2))
    im = ax.imshow(logp_mat, cmap='YlOrRd', vmin=0, vmax=6, aspect='auto')
    for fi in range(n_fp):
        for pi in range(n_pairs):
            p = pval_mat[fi, pi]
            if p < 0.001:
                s = '***'
            elif p < 0.01:
                s = '**'
            elif p < 0.05:
                s = '*'
            else:
                s = 'ns'
            ax.text(pi, fi, s, ha='center', va='center', fontsize=8)
    ax.set_xticks(range(n_pairs))
    ax.set_xticklabels(pair_labels, fontsize=7)
    ax.set_yticks(range(n_fp))
    ax.set_yticklabels([FP_DISPLAY_NAMES[fp].replace('\n', ' ') for fp in fp_names], fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label='−log₁₀(p) Bonferroni')
    ax.set_title('Pairwise Mann–Whitney U – Bonferroni-corrected significance\n'
                 '(*p<0.05, **p<0.01, ***p<0.001; ns = non-significant)',
                 fontsize=11, fontweight='bold')
    fig.tight_layout()
    _save(fig, outdir, 'fp_deep_significance', verbose=verbose)


def plot_multihalogen(
    df: pd.DataFrame,
    fps: dict[str, np.ndarray],
    outdir: Path, verbose: bool = True,
) -> None:
    """
    Analyse multi-halogen fingerprint ability to discriminate F vs Cl compounds.
    Uses Dataset A (F-only) and Dataset B (Cl-only), sub-sampled for speed.
    """
    mask_F  = (df['dataset'] == 'A_F') & (df['fg_id'].isin([f for f, *_ in FG_MIXED]))
    mask_Cl = (df['dataset'] == 'B_Cl')
    idx_F   = np.where(mask_F.values)[0][:100]   # cap at 100 each
    idx_Cl  = np.where(mask_Cl.values)[0][:100]

    categories = {
        'F–F (same FG)' : [],
        'Cl–Cl (same FG)': [],
        'F–Cl (same FG)' : [],
        'F–Cl (diff FG)' : [],
    }
    fg_F  = df.loc[df.index[idx_F], 'fg_id'].values
    fg_Cl = df.loc[df.index[idx_Cl], 'fg_id'].values

    rng  = np.random.default_rng(42)
    pairs_to_eval = {
        'F–F (same FG)': [(i, j) for i, j in zip(rng.choice(len(idx_F), 200, replace=True),
                                                    rng.choice(len(idx_F), 200, replace=True))
                          if i != j and fg_F[i] == fg_F[j]][:200],
        'Cl–Cl (same FG)': [(i, j) for i, j in zip(rng.choice(len(idx_Cl), 200, replace=True),
                                                      rng.choice(len(idx_Cl), 200, replace=True))
                             if i != j and fg_Cl[i] == fg_Cl[j]][:200],
        'F–Cl (same FG)': [(i, j) for i in range(len(idx_F))
                            for j in range(len(idx_Cl))
                            if fg_F[i] == fg_Cl[j]][:200],
        'F–Cl (diff FG)': [(i, j) for i in range(min(20, len(idx_F)))
                            for j in range(min(20, len(idx_Cl)))
                            if fg_F[i] != fg_Cl[j]][:200],
    }

    # For each FP type that has Cl sensitivity
    fp_compare = ['binary_F_per', 'binary_FCl_per']
    result = {fp: {cat: [] for cat in categories} for fp in fp_compare}

    for fp_name in fp_compare:
        X = fps[fp_name]
        X_F  = X[idx_F]
        X_Cl = X[idx_Cl]

        for cat, ps in pairs_to_eval.items():
            for (i, j) in ps:
                if 'F–F' in cat:
                    sim = gen_jaccard(X_F[i].astype(float), X_F[j].astype(float))
                elif 'Cl–Cl' in cat:
                    sim = gen_jaccard(X_Cl[i].astype(float), X_Cl[j].astype(float))
                else:  # F-Cl cross
                    sim = gen_jaccard(X_F[i].astype(float), X_Cl[j].astype(float))
                result[fp_name][cat].append(sim)

    cat_order = ['F–F (same FG)', 'Cl–Cl (same FG)', 'F–Cl (same FG)', 'F–Cl (diff FG)']
    cat_cols  = ['#1976D2', '#D32F2F', '#7B1FA2', '#888888']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, fp_name in zip(axes, fp_compare):
        data_list = [result[fp_name][c] for c in cat_order]
        parts = ax.violinplot(data_list, positions=range(len(cat_order)),
                              showmedians=True, showextrema=False, widths=0.7)
        for pc, c in zip(parts['bodies'], cat_cols):
            pc.set_facecolor(c); pc.set_alpha(0.7)
        parts['cmedians'].set_color('black')
        ax.set_xticks(range(len(cat_order)))
        ax.set_xticklabels(cat_order, rotation=25, ha='right', fontsize=9)
        ax.set_ylabel('Jaccard similarity', fontsize=10)
        ax.set_title(FP_DISPLAY_NAMES[fp_name].replace('\n', ' '), fontsize=11, fontweight='bold')
        ax.set_ylim(-0.05, 1.1)
        ax.grid(axis='y', alpha=0.25)

        # Print medians
        if verbose:
            for cat in cat_order:
                vals = result[fp_name][cat]
                if vals:
                    print(f"  [{fp_name}] {cat}: n={len(vals)} median={np.median(vals):.3f}")

    fig.suptitle(
        'Multi-Halogen Fingerprint: F vs Cl Compound Discrimination\n'
        'Only binary_FCl_per (F+Cl stacked) can separate F- from Cl-compounds sharing a functional group',
        fontsize=12, fontweight='bold',
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    _save(fig, outdir, 'fp_deep_multihalogen', verbose=verbose)


def plot_pca_compare(
    df_A: pd.DataFrame,
    fps: dict[str, np.ndarray],
    idx_A: np.ndarray,
    outdir: Path, verbose: bool = True,
) -> None:
    """PCA of Dataset A compounds for each FP type (2×4 grid)."""
    fp_names = [fp for fp in FP_DISPLAY_NAMES if fps.get(fp) is not None][:8]
    n_panels = len(fp_names)
    n_cols   = 4
    n_rows   = (n_panels + n_cols - 1) // n_cols

    fg_ids    = df_A['fg_id'].values
    unique_fg = list(dict.fromkeys(fg_ids))
    fg_to_col = {fid: FG_COLOURS[i] for i, fid in enumerate(unique_fg)}

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(4.5 * n_cols, 4 * n_rows),
                              constrained_layout=True)
    axes_flat = axes.flatten() if n_rows > 1 else list(axes)

    for ax, fp_name in zip(axes_flat[:n_panels], fp_names):
        X = fps[fp_name][idx_A].astype(float)
        # Normalise descriptive FP; others use raw values
        if fp_name == 'descriptive':
            X = _zscore(X)
        else:
            nrm = X.max(axis=0); nrm[nrm == 0] = 1
            X = X / nrm

        Z, ev = _numpy_pca(X)
        for fid in unique_fg:
            mask = fg_ids == fid
            ax.scatter(Z[mask, 0], Z[mask, 1],
                       c=fg_to_col[fid], s=8, alpha=0.6, linewidths=0)
        ax.set_title(FP_DISPLAY_NAMES[fp_name].replace('\n', ' '), fontsize=9)
        ax.set_xlabel(f'PC1 {ev[0]*100:.0f}%', fontsize=7)
        ax.set_ylabel(f'PC2 {ev[1]*100:.0f}%', fontsize=7)
        ax.tick_params(labelsize=6)

    # Legend
    handles = [mpatches.Patch(color=fg_to_col[fid], label=fid) for fid in unique_fg]
    fig.legend(handles=handles, fontsize=8, loc='lower right',
               ncol=3, title='FG type', title_fontsize=8)

    # Hide unused panels
    for ax in axes_flat[n_panels:]:
        ax.set_visible(False)

    fig.suptitle('PCA comparison across fingerprint variants (Dataset A, F-only)',
                 fontsize=13, fontweight='bold')
    _save(fig, outdir, 'fp_deep_pca_compare', verbose=verbose)


# ── CSRML comparison (optional) ───────────────────────────────────────────────

def run_csrml_comparison(
    txppfas_csv: str,
    smiles_tsv: str | None,
    fps_all: dict,
    outdir: Path,
    verbose: bool = True,
) -> None:
    """
    Compare PFASGroups pair-similarity distributions vs CSRML.

    Strategy
    --------
    1. Load CSRML fingerprints (TxP_PFAS SI Table S2 CSV).
    2. If smiles_tsv is provided, compute PFASGroups FPs for those molecules.
    3. Stratify pairs by primary PFASGroups group assignment (majority OECD bit).
    4. Compare similarity distributions for PFASGroups binary vs CSRML binary.
    """
    try:
        import sys
        sys.path.insert(0, str(SCRIPT_DIR))
        from compare_fingerprints_vs_txppfas import load_txppfas_csv, load_molecules, compute_pfasgroups_fps
    except ImportError as e:
        print(f"  Cannot import comparison utilities: {e}")
        return

    if verbose:
        print("\n--- CSRML comparison ---")
    try:
        X_csrml, chemotypes, dtxsids = load_txppfas_csv(txppfas_csv)
    except ValueError as e:
        print(f"  Skipping CSRML comparison: {e}")
        return
    if verbose:
        print(f"  CSRML matrix: {X_csrml.shape} ({len(chemotypes)} chemotypes)")

    smiles_list = None
    if smiles_tsv:
        df_smi = load_molecules(smiles_tsv)
        # Align by DTXSID if possible
        dtxsid_to_smi = dict(zip(df_smi.get('DTXSID', []), df_smi.get('SMILES', [])))
        smiles_list = [dtxsid_to_smi.get(d, None) for d in dtxsids]
        n_smiles = sum(1 for s in smiles_list if s)
        if verbose:
            print(f"  Matched {n_smiles}/{len(dtxsids)} SMILES to DTXSIDs")

    # Use DTXSIDs to define structural categories via DTXSID-to-group from
    # PFASGroups if smiles are available, else skip the PG stratification
    if smiles_list and sum(1 for s in smiles_list if s) > 50:
        valid_smiles = [(i, s) for i, s in enumerate(smiles_list) if s][:500]
        idx_valid = [i for i, _ in valid_smiles]
        smi_valid  = [s for _, s in valid_smiles]
        X_pg, gnames = compute_pfasgroups_fps(smi_valid, verbose=verbose)
        X_csrml_valid = X_csrml[idx_valid]

        # Stratify: primary PFAS group (column with highest sum in top-20 activated)
        primary_group = X_pg.argmax(axis=1)
        n_cats = min(10, len(np.unique(primary_group)))
        top_groups = np.argsort(X_pg.sum(axis=0))[::-1][:n_cats]

        # Sample pairs: same primary group vs different primary group
        rng = np.random.default_rng(42)
        same_pg, diff_pg = [], []
        n_mol = len(smi_valid)
        for _ in range(3000):
            i, j = rng.integers(0, n_mol, 2)
            if i == j:
                continue
            if primary_group[i] == primary_group[j]:
                same_pg.append((i, j))
            else:
                diff_pg.append((i, j))
        same_pg = same_pg[:500]
        diff_pg = diff_pg[:500]

        cat_data_csrml = {'Same PG': [], 'Diff PG': []}
        cat_data_pg    = {'Same PG': [], 'Diff PG': []}
        for cat, pairs in [('Same PG', same_pg), ('Diff PG', diff_pg)]:
            for i, j in pairs:
                cat_data_csrml[cat].append(gen_jaccard(X_csrml_valid[i].astype(float),
                                                        X_csrml_valid[j].astype(float)))
                cat_data_pg[cat].append(gen_jaccard(X_pg[i].astype(float),
                                                     X_pg[j].astype(float)))
        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        for ax, (name, data) in zip(axes, [('PFASGroups binary', cat_data_pg),
                                             ('CSRML TxP_PFAS', cat_data_csrml)]):
            parts = ax.violinplot(
                [data['Same PG'], data['Diff PG']],
                positions=[0, 1], showmedians=True, showextrema=False, widths=0.6)
            for pc, c in zip(parts['bodies'], ['#1976D2', '#D32F2F']):
                pc.set_facecolor(c); pc.set_alpha(0.7)
            parts['cmedians'].set_color('black')
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['Same primary\nPFAS group', 'Diff primary\nPFAS group'])
            ax.set_ylabel('Jaccard similarity')
            ax.set_title(name, fontsize=11, fontweight='bold')
            ax.set_ylim(-0.05, 1.1)
            ax.grid(axis='y', alpha=0.25)
            if verbose:
                for cat in ('Same PG', 'Diff PG'):
                    print(f"  [{name}] {cat}: n={len(data[cat])} "
                          f"median={np.median(data[cat]):.3f}")

        fig.suptitle('PFASGroups vs CSRML: Similarity within/across primary PFAS group',
                     fontsize=12, fontweight='bold')
        fig.tight_layout()
        _save(fig, outdir, 'fp_deep_csrml_compare', verbose=verbose)
    else:
        if verbose:
            print("  Not enough matched SMILES for detailed CSRML pair comparison.")
            print("  Provide --smiles_tsv with DTXSIDs matching the CSRML CSV.")


# ── Main ──────────────────────────────────────────────────────────────────────

def run(outdir: Path, article_imgs: Path | None,
        txppfas_csv: str | None = None, smiles_tsv: str | None = None,
        verbose: bool = True) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Build dataset ──────────────────────────────────────────────────────
    if verbose:
        print("\n=== Building compound dataset ===")
    df = build_dataset(verbose=verbose)
    df_A  = df[df['dataset'] == 'A_F'].reset_index(drop=True)
    idx_A = np.where(df['dataset'].values == 'A_F')[0]
    n_A   = len(df_A)
    if verbose:
        print(f"  Dataset A (F-only): {n_A} compounds")

    # ── Compute fingerprint variants ───────────────────────────────────────
    if verbose:
        print("\n=== Computing fingerprint variants ===")
    smiles_all = df['smiles'].tolist()

    fps: dict[str, np.ndarray] = {}

    # 1. binary_F_per
    if verbose: print("  [1/7] binary_F_per …")
    fps['binary_F_per'] = compute_fp_variant(smiles_all, 'F', 'per', 'binary', verbose)

    # 2. count_F_per
    if verbose: print("  [2/7] count_F_per …")
    fps['count_F_per'] = compute_fp_variant(smiles_all, 'F', 'per', 'count', verbose)

    # 3. maxcomp_F_per
    if verbose: print("  [3/7] maxcomp_F_per …")
    fps['maxcomp_F_per'] = compute_fp_variant(smiles_all, 'F', 'per', 'max_component', verbose)

    # 4. binary_F_poly
    if verbose: print("  [4/7] binary_F_poly …")
    fps['binary_F_poly'] = compute_fp_variant(smiles_all, 'F', 'poly', 'binary', verbose)

    # 5. binary_F_none (no saturation filter)
    if verbose: print("  [5/7] binary_F_none …")
    fps['binary_F_none'] = compute_fp_variant(smiles_all, 'F', None, 'binary', verbose)

    # 6. binary_FCl_per (F+Cl stacked)
    if verbose: print("  [6/7] binary_FCl_per (F+Cl) …")
    fps['binary_FCl_per'] = compute_fp_variant(smiles_all, ['F', 'Cl'], 'per', 'binary', verbose)

    # 7. descriptive (compute from A-subset for speed)
    if verbose: print("  [7/7] descriptive descriptor vector …")
    X_bin_all = fps['binary_F_per']
    X_cnt_all = fps['count_F_per']
    X_mx_all  = fps['maxcomp_F_per']
    fps['descriptive'] = compute_descriptive(smiles_all, X_bin_all, X_cnt_all, X_mx_all)

    # ── Sample pairs from Dataset A ────────────────────────────────────────
    if verbose:
        print("\n=== Sampling pairs from Dataset A ===")
    pairs = sample_pairs(df_A, max_pairs=8000)
    if verbose:
        from collections import Counter
        cnt = Counter(c for _, _, c in pairs)
        for cat, n in sorted(cnt.items()):
            print(f"  Category {cat}: {n} pairs")

    # ── Compute similarities for all FP types ─────────────────────────────
    if verbose:
        print("\n=== Computing pairwise similarities ===")
    all_cat_data: dict[str, dict[str, list]] = {fp: {c: [] for c in CAT_LABELS} for fp in fps}
    stats_dict:  dict[str, dict]             = {}

    for fp_name, X_all in fps.items():
        X_A = X_all[idx_A]
        for pi, (gi, gj, cat) in enumerate(pairs):
            a = X_A[gi].astype(float)
            b = X_A[gj].astype(float)
            if fp_name == 'descriptive':
                # Use cosine similarity after z-score (already continuous)
                a_z = _zscore(X_A)[gi]
                b_z = _zscore(X_A)[gj]
                sim = (cosine_similarity(a_z, b_z) + 1.0) / 2.0  # map [-1,1]→[0,1]
            else:
                sim = gen_jaccard(a, b)
            all_cat_data[fp_name][cat].append(float(sim))

        # Statistics
        stats_dict[fp_name] = run_stats(all_cat_data[fp_name])

    # ── Print summary ──────────────────────────────────────────────────────
    if verbose:
        print("\n=== Summary ===")
        print(f"{'FP Variant':24s}  {'η²':>5s}  {'KW-p':>10s}  "
              + "  ".join(f'{c:>6s}' for c in CAT_LABELS))
        print("-" * (24 + 5 + 10 + 6 * len(CAT_LABELS) + 20))
        for fp in fps:
            st = stats_dict[fp]
            medians = st.get('medians', {})
            row = f"{fp:24s}  {st.get('eta2', 0):5.3f}  {st.get('kruskal_p', 1):10.2e}  "
            row += "  ".join(f"{medians.get(c, float('nan')):6.3f}" for c in CAT_LABELS)
            print(row)

    # ── Save result CSV ────────────────────────────────────────────────────
    rows = []
    for fp_name in fps:
        for cat in CAT_LABELS:
            for sim in all_cat_data[fp_name][cat]:
                rows.append({'fp_variant': fp_name, 'category': cat, 'similarity': sim})
    result_df = pd.DataFrame(rows)
    result_df.to_csv(outdir / 'fp_deep_results.csv', index=False)

    stats_rows = []
    for fp_name in fps:
        st = stats_dict[fp_name]
        row = {'fp_variant': fp_name,
               'kruskal_H': st.get('kruskal_H'),
               'kruskal_p': st.get('kruskal_p'),
               'eta2': st.get('eta2')}
        for cat in CAT_LABELS:
            row[f'median_{cat}'] = st['medians'].get(cat)
            row[f'mean_{cat}']   = st['means'].get(cat)
            row[f'std_{cat}']    = st['stds'].get(cat)
        stats_rows.append(row)
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(outdir / 'fp_deep_stats.csv', index=False)
    if verbose:
        print(f"\n  Saved fp_deep_results.csv and fp_deep_stats.csv")

    # ── Generate figures ───────────────────────────────────────────────────
    if verbose:
        print("\n=== Generating figures ===")
    plot_split_violin(all_cat_data, outdir, verbose)
    plot_heatmap_summary(all_cat_data, stats_dict, outdir, verbose)
    plot_significance(stats_dict, outdir, verbose)
    plot_multihalogen(df, fps, outdir, verbose)
    plot_pca_compare(df_A, fps, idx_A, outdir, verbose)

    # ── Optional CSRML ────────────────────────────────────────────────────
    if txppfas_csv:
        run_csrml_comparison(txppfas_csv, smiles_tsv, fps, outdir, verbose)

    # ── Copy article figures ───────────────────────────────────────────────
    if article_imgs is not None and article_imgs.is_dir():
        import shutil
        for stem in ('fp_deep_split_violin', 'fp_deep_heatmap_summary',
                     'fp_deep_significance', 'fp_deep_multihalogen',
                     'fp_deep_pca_compare'):
            for ext in ('.pdf', '.png'):
                src = outdir / (stem + ext)
                dst = article_imgs / (stem + ext)
                if src.exists():
                    shutil.copy2(src, dst)
        if verbose:
            print(f"  Copied main figures to {article_imgs}")

    if verbose:
        print("\nDone.")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', type=Path,
                   default=BENCHMARK_DIR / 'data' / 'figures_fp_benchmarks')
    p.add_argument('--article_imgs', type=Path,
                   default=_ARTICLE_IMGS if _ARTICLE_IMGS.is_dir() else None)
    p.add_argument('--txppfas_csv', type=str, default=None,
                   help='Path to TxP_PFAS SI Table S2 CSV for CSRML comparison')
    p.add_argument('--smiles_tsv', type=str, default=None,
                   help='TSV with DTXSID and SMILES columns for CSRML comparison')
    p.add_argument('--quiet', action='store_true')
    args = p.parse_args(argv)
    run(args.outdir, args.article_imgs,
        txppfas_csv=args.txppfas_csv, smiles_tsv=args.smiles_tsv,
        verbose=not args.quiet)


if __name__ == '__main__':
    main()
