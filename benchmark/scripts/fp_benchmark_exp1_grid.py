"""
Fingerprint Benchmark – Experiment 1: Functional-Group × Component-Size Grid
==============================================================================

Constructs a 10 × 10 grid of synthetic PFAS (10 functional-group types × 10
perfluoroalkyl chain lengths, CF₂ units n ∈ {2,3,4,5,6,7,8,9,10,12}) plus 10
non-fluorinated references (same functional groups, n = 6 alkyl carbons, F→H),
giving 110 structures in total.

For each compound a 116-bit PFASGroups binary fingerprint is computed with
  halogens='F', saturation='per', count_mode='binary'.

Outputs
-------
fp_grid_heatmap.pdf / .png
    Ward-linked binary fingerprint heatmap (rows = compounds, cols = 116 bits).
fp_grid_pca.pdf / .png
    PCA scatter coloured by functional-group type and sized by chain length.
fp_grid_results.csv
    Full fingerprint matrix with metadata columns.

Usage
-----
python fp_benchmark_exp1_grid.py [--outdir PATH] [--article_imgs PATH]
"""
import argparse
import os
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

# ── Path setup ────────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
REPO_ROOT    = SCRIPT_DIR.parent.parent           # .../PFASGroups/
BENCHMARK_DIR = SCRIPT_DIR.parent                  # .../PFASGroups/benchmark/

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Article imgs directory (sibling repo, auto-detected if present)
_ARTICLE_IMGS = REPO_ROOT.parent / 'overleaf_PFASgroups_article' / 'imgs'

warnings.filterwarnings('ignore')

# ── Grid definition ───────────────────────────────────────────────────────────
#  Each entry: (short_id, display_name, fg_smiles_suffix)
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

CHAIN_LENGTHS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]   # CF₂ units in PF chain
REF_N_C      = 7                                     # alkyl carbons in non-F reference

# Colour palette (one per FG type)
FG_COLOURS = [
    '#1976D2', '#D32F2F', '#388E3C', '#F57C00', '#7B1FA2',
    '#0288D1', '#795548', '#00796B', '#AFB42B', '#E91E63',
]


# ── SMILES builders ───────────────────────────────────────────────────────────

def build_pfas(fg_suffix: str, n_cf2: int) -> str:
    """Linear perfluoroalkyl chain: CF3-(CF2)n_cf2-[fg].

    n_cf2 is the number of –CF2– repeating units (excluding the CF3 terminal).
    Total fluorinated carbons = n_cf2 + 1.
    """
    return 'FC(F)(F)' + 'C(F)(F)' * n_cf2 + fg_suffix


def build_nonfluor(fg_suffix: str, n_c: int = REF_N_C) -> str:
    """Non-fluorinated reference: straight alkyl chain -(CH2)n-[fg]."""
    return 'C' * n_c + fg_suffix


def validate_smiles(smiles: str) -> bool:
    """Return True if RDKit can parse the SMILES."""
    from rdkit import Chem
    return Chem.MolFromSmiles(smiles) is not None


# ── Fingerprint computation ───────────────────────────────────────────────────

def compute_fingerprints(
    smiles_list: list,
    count_mode: str = 'binary',
    verbose: bool = True,
) -> tuple:
    """Return (X: ndarray [n, 116], group_names: list[str])."""
    from HalogenGroups import generate_fingerprint

    # Probe one molecule to get dimension and group names
    _fp, _info = generate_fingerprint(
        'FC(F)(F)C(=O)O',
        representation='vector',
        count_mode='binary',
        halogens='F',
        saturation='per',
    )
    n_groups    = len(_fp)
    group_names = _info['group_names']
    n           = len(smiles_list)
    X           = np.zeros((n, n_groups), dtype=np.uint8 if count_mode == 'binary' else np.int32)

    errors = 0
    for i, smi in enumerate(smiles_list):
        try:
            fp, _ = generate_fingerprint(
                smi,
                representation='vector',
                count_mode=count_mode,
                halogens='F',
                saturation='per',
            )
            X[i] = fp
        except Exception as exc:
            errors += 1
            if verbose:
                print(f"  Warning [{i}] {smi!r}: {exc}")
    if verbose and errors:
        print(f"  Total errors: {errors}/{n}")
    return X, group_names


# ── Utility ───────────────────────────────────────────────────────────────────

def _save(fig, outdir: Path, stem: str, verbose: bool = True) -> None:
    for ext in ('.pdf', '.png'):
        path = outdir / (stem + ext)
        fig.savefig(path, dpi=150, bbox_inches='tight')
    if verbose:
        print(f"  Saved {stem}.pdf / .png in {outdir}")
    plt.close(fig)


def _numpy_pca(X: np.ndarray, n_components: int = 2):
    """Return (Z [n, n_components], ev_ratio [n_components]) via NumPy SVD."""
    Xc = X - X.mean(axis=0)
    _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    variance  = s ** 2 / max(X.shape[0] - 1, 1)
    total_var = variance.sum()
    ev_ratio  = variance[:n_components] / total_var if total_var > 0 else np.zeros(n_components)
    Z = Xc @ Vt[:n_components].T
    return Z, ev_ratio


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_heatmap(X: np.ndarray, labels: list, is_ref: list, outdir: Path) -> None:
    """Ward-clustered binary fingerprint heatmap."""
    n, p = X.shape

    # Ward linkage on rows (compounds)
    X_float = X.astype(float)
    # add tiny jitter to avoid zero-distance duplicates in Ward
    rng = np.random.default_rng(42)
    X_jitter = X_float + rng.normal(0, 1e-6, X_float.shape)
    row_linkage  = linkage(X_jitter, method='ward', metric='euclidean')
    row_order    = leaves_list(row_linkage)

    # Column ordering: group by activation density
    col_order = np.argsort(X.sum(axis=0))[::-1]

    X_ordered = X[np.ix_(row_order, col_order)]
    labels_ordered  = [labels[i] for i in row_order]
    is_ref_ordered  = [is_ref[i] for i in row_order]

    # Assign row colours based on FG type (first substring of label)
    row_colours = []
    for lbl in labels_ordered:
        fg = lbl.split('_')[0]
        idx = next((i for i, (fid, *_) in enumerate(FG_TYPES) if fid == fg), 10)
        colour = FG_COLOURS[idx % len(FG_COLOURS)] if idx < 10 else '#888888'
        row_colours.append(colour)

    fig, ax = plt.subplots(figsize=(14, max(6, n * 0.18)))
    cmap = ListedColormap(['#FFFFFF', '#1565C0'])
    ax.imshow(X_ordered, aspect='auto', cmap=cmap, vmin=0, vmax=1,
              interpolation='none')

    # Colour strip on the left (FG identity)
    ax_strip = ax.inset_axes([-0.03, 0, 0.02, 1])
    ax_strip.set_xlim(0, 1)
    ax_strip.set_ylim(0, n)
    ax_strip.axis('off')
    for row_i, c in enumerate(row_colours):
        ax_strip.add_patch(mpatches.Rectangle(
            (0, n - row_i - 1), 1, 1, color=c, linewidth=0
        ))

    ax.set_yticks(range(n))
    ax.set_yticklabels(labels_ordered, fontsize=5)
    ax.set_xlabel('Fingerprint bit (sorted by prevalence)', fontsize=11)
    ax.set_ylabel('Compound', fontsize=11)
    ax.set_title(
        'Experiment 1: Binary Fingerprint Heatmap\n'
        '(Ward-clustered rows; 100 PFAS + 10 non-fluorinated references)',
        fontsize=12, fontweight='bold',
    )

    # Legend for FG types
    handles = [
        mpatches.Patch(color=FG_COLOURS[i], label=f'{fid} – {fname}')
        for i, (fid, fname, _) in enumerate(FG_TYPES)
    ]
    handles.append(mpatches.Patch(color='#888888', label='Non-fluorinated ref.'))
    ax.legend(handles=handles, fontsize=7, bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0, framealpha=0.9)

    fig.tight_layout()
    _save(fig, outdir, 'fp_grid_heatmap')


def plot_pca(X: np.ndarray, labels: list, is_ref: list,
             chain_lengths: list, fg_ids: list, outdir: Path) -> float:
    """PCA scatter: colour = FG type, size = chain length. Returns PC1 variance."""
    Z, ev = _numpy_pca(X)
    pc1_var = float(ev[0]) * 100

    fig, ax = plt.subplots(figsize=(9, 7))
    size_map = {n: 20 + (n / 12) * 120 for n in CHAIN_LENGTHS}

    # Non-fluorinated references
    ref_mask = np.array(is_ref, dtype=bool)
    ax.scatter(Z[ref_mask, 0], Z[ref_mask, 1],
               c='#888888', marker='x', s=60, linewidths=1.5,
               zorder=3, label='Non-fluorinated ref.')

    # PFAS grid points
    for i, (fid, fname, _) in enumerate(FG_TYPES):
        mask = np.array([fg_ids[j] == fid and not is_ref[j]
                         for j in range(len(labels))], dtype=bool)
        sizes = np.array([size_map.get(chain_lengths[j], 50)
                          for j in range(len(labels)) if mask[j]])
        ax.scatter(Z[mask, 0], Z[mask, 1],
                   c=FG_COLOURS[i], s=sizes, alpha=0.85, zorder=4,
                   label=f'{fid} – {fname}', edgecolors='white', linewidths=0.3)

    # Chain-length size legend (manual)
    for n_val in [2, 6, 12]:
        ax.scatter([], [], c='grey', s=size_map[n_val],
                   label=f'n = {n_val} CF₂ units', edgecolors='white', linewidths=0.3)

    ax.set_xlabel(f'PC1 ({ev[0]*100:.1f}% variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({ev[1]*100:.1f}% variance)', fontsize=12)
    ax.set_title(
        'Experiment 1: PCA of Binary Fingerprints\n'
        '(colour = functional group type; size = chain length)',
        fontsize=12, fontweight='bold',
    )
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.02, 1), loc='upper left',
              borderaxespad=0, framealpha=0.9)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    _save(fig, outdir, 'fp_grid_pca')
    print(f"  PC1 explains {pc1_var:.1f}% of variance; PC2: {ev[1]*100:.1f}%")
    return pc1_var


# ── Main ──────────────────────────────────────────────────────────────────────

def run(outdir: Path, article_imgs: Path | None, verbose: bool = True) -> None:
    from rdkit import Chem

    outdir.mkdir(parents=True, exist_ok=True)

    # ── Build molecule list ────────────────────────────────────────────────
    smiles_list  = []
    labels       = []
    is_ref       = []
    chain_lens   = []
    fg_ids_list  = []

    invalid = []
    for fid, fname, fg_suffix in FG_TYPES:
        for n in CHAIN_LENGTHS:
            smi = build_pfas(fg_suffix, n)
            if not validate_smiles(smi):
                print(f"  Invalid SMILES [{fid}, n={n}]: {smi}")
                invalid.append(smi)
                continue
            smiles_list.append(smi)
            labels.append(f'{fid}_n{n}')
            is_ref.append(False)
            chain_lens.append(n)
            fg_ids_list.append(fid)

    # Non-fluorinated references (one per FG type, chain length = REF_N_C)
    for fid, fname, fg_suffix in FG_TYPES:
        smi = build_nonfluor(fg_suffix, REF_N_C)
        if not validate_smiles(smi):
            print(f"  Invalid non-F SMILES [{fid}]: {smi}")
            invalid.append(smi)
            continue
        smiles_list.append(smi)
        labels.append(f'{fid}_ref')
        is_ref.append(True)
        chain_lens.append(REF_N_C)
        fg_ids_list.append(fid)

    n_pfas = sum(not r for r in is_ref)
    n_ref  = sum(is_ref)
    if verbose:
        print(f"\nBuilt {n_pfas} PFAS + {n_ref} non-fluorinated references"
              f" ({len(invalid)} invalid SMILES skipped)")

    # ── Compute fingerprints ───────────────────────────────────────────────
    if verbose:
        print("Computing PFASGroups binary fingerprints …")
    X, group_names = compute_fingerprints(smiles_list, count_mode='binary', verbose=verbose)

    # Statistics
    ref_mask  = np.array(is_ref, dtype=bool)
    pfas_mask = ~ref_mask
    n_active  = (X[pfas_mask].sum(axis=1) > 0).sum()
    n_zero    = (X[pfas_mask].sum(axis=1) == 0).sum()
    n_ref_zero = (X[ref_mask].sum(axis=1) == 0).sum()
    if verbose:
        print(f"  PFAS with ≥1 active bit: {n_active}/{pfas_mask.sum()}")
        print(f"  PFAS with all-zero fingerprint: {n_zero}")
        print(f"  Non-fluorinated refs with all-zero fingerprint: {n_ref_zero}/{ref_mask.sum()}")

    # ── Save CSV ───────────────────────────────────────────────────────────
    df = pd.DataFrame(X, columns=group_names)
    df.insert(0, 'label',      labels)
    df.insert(1, 'fg_type',    fg_ids_list)
    df.insert(2, 'chain_n',    chain_lens)
    df.insert(3, 'is_ref',     is_ref)
    df.insert(4, 'smiles',     smiles_list)
    csv_path = outdir / 'fp_grid_results.csv'
    df.to_csv(csv_path, index=False)
    if verbose:
        print(f"  Saved {csv_path.name}")

    # ── Generate figures ───────────────────────────────────────────────────
    if verbose:
        print("Generating heatmap …")
    plot_heatmap(X, labels, is_ref, outdir)

    if verbose:
        print("Generating PCA scatter …")
    pc1_var = plot_pca(X, labels, is_ref, chain_lens, fg_ids_list, outdir)

    # ── Copy PDFs to article imgs/ ─────────────────────────────────────────
    if article_imgs is not None and article_imgs.is_dir():
        import shutil
        for stem in ('fp_grid_heatmap', 'fp_grid_pca'):
            for ext in ('.pdf', '.png'):
                src = outdir / (stem + ext)
                dst = article_imgs / (stem + ext)
                if src.exists():
                    shutil.copy2(src, dst)
        if verbose:
            print(f"  Copied figures to article imgs: {article_imgs}")

    if verbose:
        print("\n=== Experiment 1 summary ===")
        print(f"  Compounds: {n_pfas} PFAS + {n_ref} non-F refs")
        print(f"  PFAS all active: {n_active}/{pfas_mask.sum()}")
        print(f"  Non-F refs silent: {n_ref_zero}/{ref_mask.sum()}")
        print(f"  PC1 variance explained: {pc1_var:.1f}%")
        print("Done.")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        '--outdir', type=Path,
        default=BENCHMARK_DIR / 'data' / 'figures_fp_benchmarks',
        help='Output directory for figures and CSV.',
    )
    p.add_argument(
        '--article_imgs', type=Path,
        default=_ARTICLE_IMGS if _ARTICLE_IMGS.is_dir() else None,
        help='Article imgs/ directory to copy PDFs to.',
    )
    p.add_argument('--quiet', action='store_true', help='Suppress progress output.')
    args = p.parse_args(argv)
    run(args.outdir, args.article_imgs, verbose=not args.quiet)


if __name__ == '__main__':
    main()
