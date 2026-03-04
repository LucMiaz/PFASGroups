"""
Fingerprint Benchmark – Experiment 3: Multi-Functional Compounds
=================================================================

Constructs a 10 × 10 × 10 design:
  • 10 functional-group types (same as Experiment 1)
  • 10 component sizes  : CF₂ units n ∈ {2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
  • 10 count variants   : k ∈ {1, 2, …, 10} copies of the FG on the backbone

Each compound is built as k copies of CF₃-(CF₂)n-[FG] connected as a
multi-component SMILES (dot-separated), so PFASGroups processes them as a
single molecular record with multiple matches.

Design rationale
----------------
Binary mode   : Hamming weight driven by FG type, k- and n-invariant
Count mode    : Hamming weight scales linearly with k  (expected r² ≈ high)
Max-component : Hamming weight reflects n (chain-length class)

Outputs
-------
fp_multifunc_clustering.pdf / .png
    Ward-clustered heatmap of binary fingerprints for the 100-compound
    representative subset (k = 1).  Row colour = FG type; col density plot.
fp_multifunc_distance_violin.pdf / .png
    Pairwise Jaccard distance violin plots (three modes × four pair categories).
fp_multifunc_results.csv
    Metadata + fingerprints for all 1 000 compounds.

Usage
-----
python fp_benchmark_exp3_multifunc.py [--outdir PATH] [--article_imgs PATH]
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr

# ── Path setup ────────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
REPO_ROOT    = SCRIPT_DIR.parent.parent
BENCHMARK_DIR = SCRIPT_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

_ARTICLE_IMGS = REPO_ROOT.parent / 'overleaf_PFASgroups_article' / 'imgs'

warnings.filterwarnings('ignore')

# ── Grid definition ───────────────────────────────────────────────────────────
FG_TYPES = [
    ("PFCA",   "Carboxylic acid",  "C(=O)O"),
    ("PFSA",   "Sulfonic acid",    "S(=O)(=O)O"),
    ("PFAL",   "Alcohol",          "O"),
    ("PFAM",   "Sulfonamide",      "S(=O)(=O)N"),
    ("PFPA",   "Phosphonic acid",  "P(=O)(O)O"),
    ("Ketone", "Ketone",           "C(=O)C"),
    ("Ald",    "Aldehyde",         "C=O"),
    ("PFPE",   "Ether",            "OC"),
    ("Amine",  "Amine",            "N"),
    ("Nitro",  "Nitro",            "[N+](=O)[O-]"),
]

COMPONENT_SIZES  = list(range(2, 12))    # n_cf2 ∈ {2,...,11}  (10 sizes)
FG_COUNTS        = list(range(1, 11))    # k ∈ {1,...,10}      (10 variants)

# Colour palette – one per FG type
FG_COLOURS = [
    '#1976D2', '#D32F2F', '#388E3C', '#F57C00', '#7B1FA2',
    '#0288D1', '#795548', '#00796B', '#AFB42B', '#E91E63',
]


# ── SMILES builder ────────────────────────────────────────────────────────────

def build_multifunc(fg_suffix: str, n_cf2: int, k: int) -> str:
    """k copies of CF3-(CF2)n-[fg] as a multi-component SMILES."""
    single = 'FC(F)(F)' + 'C(F)(F)' * n_cf2 + fg_suffix
    return '.'.join([single] * k)


# ── Fingerprint computation ───────────────────────────────────────────────────

def compute_fps_mode(smiles_list: list, count_mode: str) -> tuple:
    """Compute fingerprints for one count_mode; return (X, group_names)."""
    from HalogenGroups import generate_fingerprint

    _fp, _info = generate_fingerprint(
        'FC(F)(F)C(=O)O',
        representation='vector', count_mode='binary', halogens='F', saturation='per',
    )
    n_groups, group_names = len(_fp), _info['group_names']
    n = len(smiles_list)
    X = np.zeros((n, n_groups), dtype=np.float32)

    for i, smi in enumerate(smiles_list):
        try:
            fp, _ = generate_fingerprint(
                smi, representation='vector', count_mode=count_mode,
                halogens='F', saturation='per',
            )
            X[i] = fp
        except Exception:
            pass
    return X, group_names


# ── Pairwise Jaccard ──────────────────────────────────────────────────────────

def generalised_jaccard(a: np.ndarray, b: np.ndarray) -> float:
    num = float(np.sum(np.minimum(a, b)))
    den = float(np.sum(np.maximum(a, b)))
    return num / den if den > 0 else 0.0


def pairwise_jaccard_matrix(X: np.ndarray) -> np.ndarray:
    """n×n pairwise generalised Jaccard similarity matrix."""
    n = X.shape[0]
    J = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i, n):
            v = generalised_jaccard(X[i], X[j])
            J[i, j] = J[j, i] = v
    return J


def pairwise_jaccard_distance(X: np.ndarray) -> np.ndarray:
    """Return condensed distance array (1 - Jaccard) computed row-by-row."""
    n = X.shape[0]
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            dists.append(1.0 - generalised_jaccard(X[i], X[j]))
    return np.array(dists, dtype=np.float32)


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
    Z = Xc @ Vt[:n_components].T
    return Z, ev_ratio


# ── Silhouette ────────────────────────────────────────────────────────────────

def silhouette_score_manual(X: np.ndarray, labels: np.ndarray) -> tuple:
    """Compute mean ± std silhouette from pairwise cosine distance."""
    from scipy.spatial.distance import cdist
    D = cdist(X, X, metric='cosine')
    n = len(labels)
    s = np.zeros(n)
    unique_labels = np.unique(labels)
    for i in range(n):
        same   = D[i, labels == labels[i]]
        others = [D[i, labels == lbl].mean() for lbl in unique_labels if lbl != labels[i]]
        a = same[same > 0].mean() if (same > 0).any() else 0.0
        b = min(others) if others else 0.0
        denom = max(a, b)
        s[i] = (b - a) / denom if denom > 0 else 0.0
    return float(s.mean()), float(s.std())


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_clustering_heatmap(
    X_binary: np.ndarray,
    fg_ids: list,
    ns: list,
    outdir: Path,
    verbose: bool = True,
) -> None:
    """Ward-clustered heatmap for the representative 100-compound subset (k=1)."""
    n_comps, p = X_binary.shape

    # Cosine distance → Ward linkage
    row_dists = pdist(X_binary.astype(float) + 1e-8, metric='cosine')
    row_link  = linkage(row_dists, method='ward')
    row_order = leaves_list(row_link)

    # Column ordering: by activation density
    col_order = np.argsort(X_binary.sum(axis=0))[::-1]

    X_ord = X_binary[np.ix_(row_order, col_order)]
    fg_ord = [fg_ids[i] for i in row_order]
    n_ord  = [ns[i]     for i in row_order]

    # Row-colour strip = FG type
    fg_id_list = [fid for fid, *_ in FG_TYPES]
    row_colours = [FG_COLOURS[fg_id_list.index(f)] if f in fg_id_list else '#888888'
                   for f in fg_ord]

    fig, (ax, ax_bar) = plt.subplots(
        1, 2, figsize=(14, 8),
        gridspec_kw={'width_ratios': [1, 0.08]},
    )

    cmap = ListedColormap(['#FFFFFF', '#1565C0'])
    ax.imshow(X_ord, aspect='auto', cmap=cmap, vmin=0, vmax=1, interpolation='none')

    # Colour strip on the right (component size)
    n_unique = sorted(set(ns))
    n_cmap   = plt.cm.get_cmap('YlOrRd', len(n_unique))
    ax_bar.set_xlim(0, 1)
    ax_bar.set_ylim(0, n_comps)
    ax_bar.axis('off')
    for row_i, n_val in enumerate(n_ord):
        c_idx = n_unique.index(n_val) / max(len(n_unique) - 1, 1)
        ax_bar.add_patch(mpatches.Rectangle(
            (0, n_comps - row_i - 1), 1, 1,
            color=n_cmap(c_idx), linewidth=0,
        ))
    ax_bar.set_title('n', fontsize=9)

    # Left row-colour strip
    ax_strip = ax.inset_axes([-0.025, 0, 0.018, 1])
    ax_strip.set_xlim(0, 1); ax_strip.set_ylim(0, n_comps)
    ax_strip.axis('off')
    for row_i, c in enumerate(row_colours):
        ax_strip.add_patch(mpatches.Rectangle(
            (0, n_comps - row_i - 1), 1, 1, color=c, linewidth=0
        ))

    ax.set_xlabel('Fingerprint bit (sorted by prevalence)', fontsize=11)
    ax.set_ylabel('Compound (Ward-ordered)', fontsize=11)
    ax.set_yticks([])
    ax.set_title(
        'Experiment 3: Ward-Clustered Binary Fingerprints\n'
        '(100 representative compounds; k=1; row colour = FG type)',
        fontsize=11, fontweight='bold',
    )

    # FG-type legend
    handles = [mpatches.Patch(color=FG_COLOURS[i], label=f'{fid} – {fname}')
               for i, (fid, fname, _) in enumerate(FG_TYPES)]
    ax.legend(handles=handles, fontsize=7.5, bbox_to_anchor=(1.15, 1),
              loc='upper left', framealpha=0.9)

    fig.tight_layout()
    _save(fig, outdir, 'fp_multifunc_clustering', verbose=verbose)


def plot_distance_violin(
    D_dict: dict,
    fg_ids_all: list,
    ns_all: list,
    ks_all: list,
    outdir: Path,
    verbose: bool = True,
) -> None:
    """Pairwise Jaccard distance violin plots stratified by pair category."""
    n = len(fg_ids_all)
    # For 1000 compounds, 1000×999/2 = 499500 pairs: subsample to keep it fast
    # Use a deterministic random subsample of up to 20000 pairs
    rng   = np.random.default_rng(42)
    all_i = np.arange(n, dtype=np.int32)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, min(i + 50, n))]
    # Full sampling is too slow for 1000; use 5000 random pairs
    i_idx    = rng.integers(0, n, size=5000)
    j_idx    = rng.integers(0, n, size=5000)
    keep     = i_idx != j_idx
    i_idx, j_idx = i_idx[keep], j_idx[keep]
    # Ensure i < j
    swap = i_idx > j_idx
    i_idx[swap], j_idx[swap] = j_idx[swap], i_idx[swap]
    pair_idx = list(set(zip(i_idx.tolist(), j_idx.tolist())))

    if verbose:
        print(f"  Violin plot: evaluating {len(pair_idx)} random pairs …")

    categories = {
        'Same FG\nSame n':   [],
        'Same FG\nDiff n':   [],
        'Diff FG\nSame n':   [],
        'Diff FG\nDiff n':   [],
    }
    cat_data = {mode: {cat: [] for cat in categories} for mode in D_dict}

    for i, j in pair_idx:
        same_fg = (fg_ids_all[i] == fg_ids_all[j])
        same_n  = (ns_all[i]     == ns_all[j])
        if same_fg and same_n:
            cat = 'Same FG\nSame n'
        elif same_fg:
            cat = 'Same FG\nDiff n'
        elif same_n:
            cat = 'Diff FG\nSame n'
        else:
            cat = 'Diff FG\nDiff n'
        for mode, X in D_dict.items():
            d = 1.0 - generalised_jaccard(X[i], X[j])
            cat_data[mode][cat].append(d)

    # Print medians for reporting
    if verbose:
        for mode in D_dict:
            for cat in list(categories.keys())[:2]:   # focus on same-FG categories
                vals = cat_data[mode][cat]
                if vals:
                    print(f"  [{mode}] {cat:25s}: median={np.median(vals):.3f}, "
                          f"n={len(vals)}")

    # Plot
    cat_labels = list(categories.keys())
    n_cats = len(cat_labels)
    mode_labels = {'binary': 'Binary', 'count': 'Count', 'max_component': 'Max-component'}
    mode_colours = {'binary': '#1976D2', 'count': '#D32F2F', 'max_component': '#388E3C'}

    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
    for ax, mode in zip(axes, D_dict):
        data_per_cat = [cat_data[mode][c] for c in cat_labels]
        parts = ax.violinplot(
            data_per_cat,
            positions=range(n_cats),
            showmedians=True, showextrema=False,
        )
        for pc in parts['bodies']:
            pc.set_facecolor(mode_colours[mode])
            pc.set_alpha(0.7)
        parts['cmedians'].set_color('black')
        ax.set_xticks(range(n_cats))
        ax.set_xticklabels(cat_labels, fontsize=9)
        ax.set_title(f'{mode_labels[mode]} mode', fontsize=12, fontweight='bold')
        ax.set_xlabel('Pair category', fontsize=10)
        ax.grid(axis='y', alpha=0.3)

    axes[0].set_ylabel('Jaccard distance (1 − J)', fontsize=11)
    fig.suptitle(
        'Experiment 3: Pairwise Jaccard Distance by Pair Category\n'
        '(sampled from 1 000-compound multi-functional PFAS dataset)',
        fontsize=12, fontweight='bold',
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    _save(fig, outdir, 'fp_multifunc_distance_violin', verbose=verbose)


# ── Main ──────────────────────────────────────────────────────────────────────

def run(outdir: Path, article_imgs: Path | None, verbose: bool = True) -> None:
    from rdkit import Chem

    outdir.mkdir(parents=True, exist_ok=True)

    # ── Build 1000-compound dataset ───────────────────────────────────────
    if verbose:
        print(f"\nBuilding {len(FG_TYPES)} × {len(COMPONENT_SIZES)} × {len(FG_COUNTS)}"
              " = 1000 multi-functional PFAS …")
    records = []
    for fid, fname, fg_suffix in FG_TYPES:
        for n in COMPONENT_SIZES:
            for k in FG_COUNTS:
                smi = build_multifunc(fg_suffix, n, k)
                mol = Chem.MolFromSmiles(smi)
                if mol is None and verbose:
                    print(f"  Invalid: {fid}, n={n}, k={k}")
                records.append({
                    'fg_id':  fid,
                    'fg_name': fname,
                    'n_cf2':  n,
                    'k':      k,
                    'smiles': smi,
                    'valid':  mol is not None,
                })

    n_total = len(records)
    n_valid = sum(r['valid'] for r in records)
    if verbose:
        print(f"  Valid SMILES: {n_valid}/{n_total}")

    smiles_all = [r['smiles'] for r in records]
    fg_ids_all = [r['fg_id'] for r in records]
    ns_all     = [r['n_cf2'] for r in records]
    ks_all     = [r['k']     for r in records]

    # ── Compute fingerprints ───────────────────────────────────────────────
    if verbose:
        print("Computing fingerprints (binary, count, max_component) …")
    fps = {}
    group_names = None
    for mode in ('binary', 'count', 'max_component'):
        if verbose:
            print(f"  mode={mode} …")
        X, gnames = compute_fps_mode(smiles_all, mode)
        fps[mode] = X
        if group_names is None:
            group_names = gnames

    X_bin = fps['binary']

    # ── Summary statistics ─────────────────────────────────────────────────
    if verbose:
        print("\nSummary:")
        # 1. r² for count-mode Hamming weight vs. k
        hw_count = fps['count'].sum(axis=1)
        r, pval  = pearsonr(np.array(ks_all), hw_count)
        print(f"  Count-mode Hamming weight vs. k: r={r:.3f}, r²={r**2:.3f} (p={pval:.2e})")

        # 2. Binary silhouette by FG type
        fg_label_arr = np.array([list(dict.fromkeys(fg_ids_all)).index(f) for f in fg_ids_all])
        fg_id_order  = list(dict.fromkeys(fg_ids_all))   # unique FG types in insertion order
        g_arr        = np.array([fg_id_order.index(f) for f in fg_ids_all])
        sil_mean, sil_std = silhouette_score_manual(X_bin.astype(float), g_arr)
        print(f"  Binary fingerprint silhouette (by FG type): {sil_mean:.3f} ± {sil_std:.3f}")

    # ── Save CSV ───────────────────────────────────────────────────────────
    meta_df = pd.DataFrame(records)
    fp_df   = pd.DataFrame(X_bin, columns=group_names)
    out_df  = pd.concat([meta_df.reset_index(drop=True), fp_df], axis=1)
    csv_path = outdir / 'fp_multifunc_results.csv'
    out_df.to_csv(csv_path, index=False)
    if verbose:
        print(f"  Saved {csv_path.name}")

    # ── Clustering heatmap (k=1 subset: 100 compounds) ────────────────────
    if verbose:
        print("Generating clustering heatmap (k=1 subset) …")
    k1_mask   = np.array([r['k'] == 1 for r in records])
    X_k1      = X_bin[k1_mask]
    fg_ids_k1 = [fg_ids_all[i] for i in range(n_total) if k1_mask[i]]
    ns_k1     = [ns_all[i]     for i in range(n_total) if k1_mask[i]]
    plot_clustering_heatmap(X_k1, fg_ids_k1, ns_k1, outdir, verbose=verbose)

    # ── Violin plots (all 1000 compounds, sampled pairs) ──────────────────
    if verbose:
        print("Generating violin plots …")
    plot_distance_violin(fps, fg_ids_all, ns_all, ks_all, outdir, verbose=verbose)

    # ── Copy PDFs to article imgs/ ─────────────────────────────────────────
    if article_imgs is not None and article_imgs.is_dir():
        import shutil
        for stem in ('fp_multifunc_clustering', 'fp_multifunc_distance_violin'):
            for ext in ('.pdf', '.png'):
                src = outdir / (stem + ext)
                dst = article_imgs / (stem + ext)
                if src.exists():
                    shutil.copy2(src, dst)
        if verbose:
            print(f"  Copied figures to article imgs: {article_imgs}")

    if verbose:
        print("\n=== Experiment 3 summary ===")
        print(f"  Total compounds: {n_total} ({n_valid} valid)")
        print(f"  Grid: {len(FG_TYPES)} FG types × {len(COMPONENT_SIZES)} sizes "
              f"× {len(FG_COUNTS)} k-variants")
        print("Done.")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        '--outdir', type=Path,
        default=BENCHMARK_DIR / 'data' / 'figures_fp_benchmarks',
    )
    p.add_argument(
        '--article_imgs', type=Path,
        default=_ARTICLE_IMGS if _ARTICLE_IMGS.is_dir() else None,
    )
    p.add_argument('--quiet', action='store_true')
    args = p.parse_args(argv)
    run(args.outdir, args.article_imgs, verbose=not args.quiet)


if __name__ == '__main__':
    main()
