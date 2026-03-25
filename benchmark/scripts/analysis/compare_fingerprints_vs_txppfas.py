#!/usr/bin/env python3
"""
Compare PFASGroups fingerprints (PFASGroups) against TxP_PFAS (CSRML) chemotype fingerprints.

Reference:
  Richard et al., "A New CSRML Structure-Based Fingerprint Method for Profiling and
  Categorizing Per- and Polyfluoroalkyl Substances (PFAS)"
  Chem. Res. Toxicol. 2023, 36, 508-534. DOI: 10.1021/acs.chemrestox.2c00403

TxP_PFAS (129 chemotypes) consists of:
  - ~56 modified ToxPrint bond-type functional groups (enforced proximity to CF/F)
  - ~73 PFAS-specific chain features: linear perfluoro chains (capped/uncapped, C1-C11+),
    branching, fluorotelomers (n1/n2/n3), alternate halogenation (Cl/Br/I), cyclics/aromatics

PFASGroups / PFASGroups: algorithmic detection of perhalogenated/polyhalogenated
  structural components (88 groups); binary fingerprint optionally per halogen.

What this script computes
-------------------------
  1. PFASGroups binary fingerprints via generate_fingerprint()
  2. (Optional) TxP_PFAS fingerprints loaded from the SI Table S2 CSV
     (exported from ChemoTyper for PFASSTRUCTV5; 14 735 rows x 129 columns)
  3. Coverage statistics: % molecules hitting each feature
  4. Per-feature information content (Shannon entropy of the bit distribution)
  5. Pairwise Jaccard / Tanimoto similarity distributions
  6. Cross-fingerprint pairwise similarity correlation (Pearson r, Spearman rho)
  7. PCA scatter plots of each fingerprint space
  8. Coverage concentration curve (Lorenz-style: cumulative coverage vs. features used)
  9. Summary CSV

Usage
-----
  # PFASGroups-only analysis (no TxP_PFAS CSV needed):
  python compare_pfasgroups_vs_txppfas.py

  # Full comparison (requires SI Table S2 CSV):
  python compare_pfasgroups_vs_txppfas.py --txppfas_csv /path/to/txp_pfas_fingerprints.csv

  # Custom molecule list:
  python compare_pfasgroups_vs_txppfas.py --smiles /path/to/molecules.tsv [--txppfas_csv ...]

  # Quick test with a small subset:
  python compare_pfasgroups_vs_txppfas.py --max_mols 200

How to get the TxP_PFAS CSV
----------------------------
  1. Download the SI zip:
       https://pubs.acs.org/doi/suppl/10.1021/acs.chemrestox.2c00403/suppl_file/tx2c00403_si_001.zip
  2. Extract SuppInfo_TablesS1-S5_*.xlsx
  3. Save the "Table S2" sheet (14 735 rows, 129 TxP_PFAS columns, indexed by DTXSID) as CSV.
     Example with pandas:
       import pandas as pd
       df = pd.read_excel('SuppInfo_TablesS1-S5_23Jan2023.xlsx', sheet_name='Table S2', index_col=0)
       df.to_csv('txp_pfas_fingerprints_si_s2.csv')
  4. Re-run this script with --txppfas_csv txp_pfas_fingerprints_si_s2.csv
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path

import pytest
import numpy as np
import pandas as pd

# ── matplotlib backend must be set before any pyplot import ──────────────────
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr

# ── Path setup ────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent          # .../PFASGroups/
BENCHMARK_DIR = SCRIPT_DIR.parents[1]              # .../PFASGroups/benchmark/
DATA_DIR  = BENCHMARK_DIR / 'test_data'        # .../PFASGroups/benchmark/test_data

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

warnings.filterwarnings('ignore')

try:
    import tqdm
    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False

# ══════════════════════════════════════════════════════════════════════════════
# Utility functions
# ══════════════════════════════════════════════════════════════════════════════

def jaccard_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """Binary Tanimoto / Jaccard similarity between two 1-D binary arrays."""
    and_ = int(np.dot(a, b))
    or_  = int(np.sum(a | b))
    return and_ / or_ if or_ > 0 else 0.0


def pairwise_jaccard(X: np.ndarray) -> np.ndarray:
    """Return n×n pairwise Jaccard similarity matrix for binary matrix X (n_mols × n_feat)."""
    n = X.shape[0]
    sim = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i, n):
            s = jaccard_similarity(X[i].astype(bool), X[j].astype(bool))
            sim[i, j] = sim[j, i] = s
    return sim


def information_content(freq: np.ndarray) -> np.ndarray:
    """Shannon IC per feature bit: H(p) = -p*log2(p) - (1-p)*log2(1-p)."""
    eps = 1e-12
    p = np.clip(freq, eps, 1 - eps)
    return -(p * np.log2(p)) - ((1 - p) * np.log2(1 - p))


# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════

def load_molecules(path: str) -> pd.DataFrame:
    """
    Load molecules from a tab-separated file.
    Required column:  SMILES
    Optional columns: DTXSID, name, Included
    Returns a DataFrame (rows with missing SMILES dropped).
    """
    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = [c.strip() for c in df.columns]

    if 'SMILES' not in df.columns:
        raise ValueError(
            f"No 'SMILES' column found in {path}. Available: {df.columns.tolist()}"
        )

    before = len(df)
    df = df.dropna(subset=['SMILES'])
    df = df[df['SMILES'].str.strip().astype(bool)]
    after  = len(df)
    if before != after:
        print(f"  (Dropped {before - after} rows with missing SMILES)")
    df = df.reset_index(drop=True)
    return df


def compute_pfasgroups_fps(
    smiles_list: list[str],
    halogens='F',
    verbose: bool = True,
) -> tuple[np.ndarray, list[str]]:
    """
    Generate binary PFASGroups fingerprints for a list of SMILES.

    generate_fingerprint is decorated with @load_PFASGroups() which automatically
    injects the pfas_groups (list of HalogenGroup objects); do NOT pass pfas_groups
    explicitly to avoid overriding with the wrong type.

    Parameters
    ----------
    smiles_list : list of SMILES strings
    halogens    : str or list-of-str passed to generate_fingerprint (default 'F')
    verbose     : show progress bar

    Returns
    -------
    X           : np.ndarray of shape (n, n_groups) – binary (uint8)
    group_names : list of str, length n_groups
    """
    from PFASGroups import generate_fingerprint

    # Probe one molecule to get n_groups and group_names (use a simple valid SMILES)
    _probe_fp, _probe_info = generate_fingerprint(
        'FC(F)(F)C(=O)O',
        representation='vector',
        count_mode='binary',
        halogens=halogens,
        saturation=None,
    )
    n_groups    = len(_probe_fp)
    group_names = _probe_info['group_names']
    n           = len(smiles_list)
    X           = np.zeros((n, n_groups), dtype=np.uint8)
    errors      = 0

    iterable = (
        enumerate(tqdm.tqdm(smiles_list, desc='PFASGroups FP', leave=False))
        if _HAS_TQDM and verbose
        else enumerate(smiles_list)
    )

    for i, smi in iterable:
        try:
            fp, _ = generate_fingerprint(
                smi,
                representation='vector',
                count_mode='binary',
                halogens=halogens,
                saturation=None,
            )
            X[i] = fp
        except Exception:
            errors += 1

    if verbose and errors:
        print(f"  Warning: {errors} molecules failed (invalid SMILES / RDKit error)")

    return X, group_names


def _find_dtxsid_header_row(path: str) -> int:
    """
    Scan a CSV file for the row index where the first cell is 'DTXSID'.
    Returns the 0-based row index, or -1 if not found.
    """
    with open(path, 'r', encoding='utf-8-sig') as fh:
        for i, line in enumerate(fh):
            first_cell = line.split(',')[0].strip().strip('"').upper()
            if first_cell == 'DTXSID':
                return i
    return -1


def load_txppfas_csv(path: str) -> tuple[np.ndarray, list[str], list[str]]:
    """
    Load the TxP_PFAS fingerprint CSV (SI Table S2).

    Automatically skips metadata / description rows that appear before
    the actual data matrix (e.g. when the CSV was exported directly from
    the Excel workbook rather than from the data sheet alone).

    Expects a row whose first cell is 'DTXSID' followed by chemotype-name
    columns, then numeric (0/1) data rows.
    Returns: (X, chemotype_names, dtxsids)

    Raises
    ------
    ValueError
        If no 'DTXSID' header row is found.  This typically means the file
        is the metadata placeholder rather than the actual fingerprint matrix.
        Extract the correct sheet from the SI Excel file first::

            import pandas as pd
            df = pd.read_excel('SuppInfo_TablesS1-S5_23Jan2023.xlsx',
                               sheet_name='Table S2', index_col=0)
            df.to_csv('Richard2023_SI_TableS2.csv')
    """
    header_row = _find_dtxsid_header_row(path)
    if header_row < 0:
        raise ValueError(
            f"No 'DTXSID' header row found in {path}.\n"
            "The TxP_PFAS CSV must be the actual fingerprint matrix from SI Table S2,"
            " not the metadata/description file.\n"
            "Export it from the SI Excel file:\n"
            "  import pandas as pd\n"
            "  df = pd.read_excel('SuppInfo_TablesS1-S5_23Jan2023.xlsx',\n"
            "                     sheet_name='Table S2', index_col=0)\n"
            "  df.to_csv('Richard2023_SI_TableS2.csv')"
        )

    df = pd.read_csv(path, index_col=0, skiprows=header_row, encoding='utf-8-sig')
    df = df.dropna(how='all')        # drop blank spacer rows Excel sometimes adds
    dtxsids    = list(df.index)
    chemotypes = list(df.columns)
    X          = df.values.astype(np.uint8)
    return X, chemotypes, dtxsids


# ══════════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════════

_BLUE = '#1976D2'   # PFASGroups colour
_RED  = '#D32F2F'   # TxP_PFAS colour


def _save(fig: plt.Figure, outdir: str, fname: str, verbose: bool = True) -> None:
    path = os.path.join(outdir, fname)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    if verbose:
        print(f"  Saved {fname}")


def plot_coverage(freq_pg, names_pg, outdir, freq_tx=None, names_tx=None):
    """Feature prevalence bar chart (one or two panels)."""
    n_panels = 2 if freq_tx is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(9 * n_panels, 6))
    if n_panels == 1:
        axes = [axes]

    datasets = [(freq_pg, names_pg, 'PFASGroups (PFASGroups)', _BLUE)]
    if freq_tx is not None:
        datasets.append((freq_tx, names_tx, 'TxP_PFAS (CSRML)', _RED))

    for ax, (freq, names, title, color) in zip(axes, datasets):
        order = np.argsort(freq)[::-1]
        ax.bar(range(len(freq)), freq[order] * 100, color=color, alpha=0.85)
        ax.set_xlabel('Feature (sorted by prevalence)', fontsize=11)
        ax.set_ylabel('% molecules', fontsize=11)
        ax.set_title(f'{title}\n({len(freq)} features)', fontsize=11)
        ax.set_xlim(-0.5, len(freq) - 0.5)

    label = 'comparison' if n_panels == 2 else 'pfasgroups'
    fig.suptitle('Feature Coverage Across PFAS Dataset', fontsize=13, fontweight='bold')
    fig.tight_layout()
    _save(fig, outdir, f'{label}_coverage.png')


def plot_info_content(freq_pg, names_pg, outdir, freq_tx=None, names_tx=None):
    """Per-feature Shannon information content."""
    ic_pg = information_content(freq_pg)
    n_panels = 2 if freq_tx is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(9 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    datasets = [(ic_pg, 'PFASGroups', _BLUE)]
    if freq_tx is not None:
        ic_tx = information_content(freq_tx)
        datasets.append((ic_tx, 'TxP_PFAS', _RED))
    else:
        ic_tx = None

    for ax, (ic, label, color) in zip(axes, datasets):
        order = np.argsort(ic)[::-1]
        ax.bar(range(len(ic)), ic[order], color=color, alpha=0.85)
        ax.axhline(1.0, linestyle='--', color='black', alpha=0.45, label='Max IC = 1 bit')
        ax.set_xlabel('Feature (sorted by IC)', fontsize=11)
        ax.set_ylabel('IC (bits)', fontsize=11)
        ax.set_title(f'{label}: per-feature IC\nMean = {ic.mean():.3f} bits, '
                     f'{(ic > 0.5).sum()} bits with IC > 0.5', fontsize=10)
        ax.legend(fontsize=9)

    title_suffix = '(comparison)' if n_panels == 2 else ''
    fig.suptitle(f'Per-Feature Information Content {title_suffix}', fontsize=12, fontweight='bold')
    fig.tight_layout()
    _save(fig, outdir, f'{"comparison" if n_panels == 2 else "pfasgroups"}_info_content.png')

    return ic_pg, ic_tx


def plot_lorenz(freq_pg, outdir, freq_tx=None):
    """
    Coverage concentration curve: cumulative fraction of total coverage
    as a function of fraction of features used (sorted by descending prevalence).
    """
    fig, ax = plt.subplots(figsize=(7, 6))

    for freq, label, color in [
        (freq_pg, 'PFASGroups', _BLUE),
        *( [(freq_tx, 'TxP_PFAS', _RED)] if freq_tx is not None else [] ),
    ]:
        s = np.sort(freq)[::-1]
        total = s.sum()
        cum = np.cumsum(s) / total if total > 0 else np.zeros_like(s)
        cum = np.clip(cum, 0, 1)
        x = np.arange(1, len(s) + 1) / len(s)
        ax.plot(x, cum, label=f'{label} ({len(s)} features)', color=color, lw=2)

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.4, label='Uniform distribution')
    ax.set_xlabel('Fraction of features (sorted by prevalence)', fontsize=12)
    ax.set_ylabel('Cumulative fraction of total coverage', fontsize=12)
    ax.set_title('Coverage Concentration Curve\n'
                 '(fewer features = more concentrated; better for identification)', fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, outdir, 'coverage_lorenz.png')


def _numpy_pca(X: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """
    Minimal PCA using NumPy SVD.
    Returns (transformed coordinates, explained_variance_ratio).
    """
    Xc = X - X.mean(axis=0)
    _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    variance = s ** 2 / (X.shape[0] - 1)
    total_var = variance.sum()
    ev_ratio  = variance[:n_components] / total_var if total_var > 0 else np.zeros(n_components)
    Z = Xc @ Vt[:n_components].T
    return Z, ev_ratio


def plot_pca(X_pg, outdir, X_tx=None):
    """PCA scatter: one or two panels (NumPy-based, no sklearn required)."""
    n_panels = 2 if X_tx is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 6))
    if n_panels == 1:
        axes = [axes]

    datasets = [(X_pg, 'PFASGroups (PFASGroups)', _BLUE)]
    if X_tx is not None:
        datasets.append((X_tx, 'TxP_PFAS (CSRML)', _RED))

    for ax, (X, title, color) in zip(axes, datasets):
        var_mask = X.var(axis=0) > 0
        Xv = X[:, var_mask].astype(float)

        if Xv.shape[1] < 2:
            ax.text(0.5, 0.5, 'Insufficient feature variance\n(too few molecules?)',
                    ha='center', va='center', transform=ax.transAxes, fontsize=11)
            ax.set_title(title, fontsize=11)
            continue

        Z, ev = _numpy_pca(Xv, n_components=2)

        ax.scatter(Z[:, 0], Z[:, 1], c=color, alpha=0.5, s=10, edgecolors='none')
        ax.set_xlabel(f'PC1 ({ev[0]*100:.1f}%)', fontsize=11)
        ax.set_ylabel(f'PC2 ({ev[1]*100:.1f}%)', fontsize=11)
        ax.set_title(f'{title}\n({Xv.shape[1]} non-zero-variance features)', fontsize=11)

    fig.suptitle('PCA of Fingerprint Spaces', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fname = 'comparison_pca.png' if n_panels == 2 else 'pfasgroups_pca.png'
    _save(fig, outdir, fname)


def plot_jaccard_dist(sim_pg, sim_tx, outdir):
    """Histogram of pairwise Jaccard similarities (upper triangle) for both FPs."""
    pg_vals = sim_pg[np.triu_indices_from(sim_pg, k=1)]
    tx_vals = sim_tx[np.triu_indices_from(sim_tx, k=1)]

    fig, ax = plt.subplots(figsize=(9, 5))
    bins = np.linspace(0, 1, 51)
    ax.hist(pg_vals, bins=bins, alpha=0.6, color=_BLUE, density=True,
            label=f'PFASGroups  (median={np.median(pg_vals):.2f})')
    ax.hist(tx_vals, bins=bins, alpha=0.6, color=_RED,  density=True,
            label=f'TxP_PFAS  (median={np.median(tx_vals):.2f})')
    ax.set_xlabel('Pairwise Jaccard similarity', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Distribution of Pairwise Molecular Similarities', fontsize=12)
    ax.legend(fontsize=10)
    fig.tight_layout()
    _save(fig, outdir, 'jaccard_distributions.png')


def plot_cross_similarity(sim_pg, sim_tx, outdir) -> tuple[float, float]:
    """
    Scatter: PFASGroups pairwise Jaccard vs TxP_PFAS pairwise Jaccard.
    Returns (pearson_r, spearman_rho).
    """
    pg_vals = sim_pg[np.triu_indices_from(sim_pg, k=1)]
    tx_vals = sim_tx[np.triu_indices_from(sim_tx, k=1)]

    # Sub-sample if too many pairs
    max_pts = 8000
    n_pairs = len(pg_vals)
    if n_pairs > max_pts:
        rng = np.random.default_rng(0)
        idx     = rng.choice(n_pairs, max_pts, replace=False)
        pg_vals = pg_vals[idx]
        tx_vals = tx_vals[idx]

    r,   _ = pearsonr(pg_vals, tx_vals)
    rho, _ = spearmanr(pg_vals, tx_vals)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(pg_vals, tx_vals, alpha=0.25, s=6, c='#424242', edgecolors='none')
    ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, lw=1, label='x = y')
    ax.set_xlabel('PFASGroups pairwise Jaccard', fontsize=12)
    ax.set_ylabel('TxP_PFAS pairwise Jaccard', fontsize=12)
    ax.set_title(
        f'Cross-Fingerprint Pairwise Similarity\n'
        f'Pearson r = {r:.3f}   Spearman ρ = {rho:.3f}',
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.tight_layout()
    _save(fig, outdir, 'cross_similarity.png')
    return r, rho


# ══════════════════════════════════════════════════════════════════════════════
# Summary CSV
# ══════════════════════════════════════════════════════════════════════════════

def save_feature_summary(
    freq_pg, names_pg, ic_pg,
    outdir,
    freq_tx=None, names_tx=None, ic_tx=None,
):
    """Write per-feature statistics to CSV files."""
    df_pg = pd.DataFrame({
        'group_name': names_pg,
        'prevalence': freq_pg,
        'info_content_bits': ic_pg,
    }).sort_values('prevalence', ascending=False)
    df_pg.to_csv(os.path.join(outdir, 'pfasgroups_feature_stats.csv'), index=False)
    print(f"  Saved pfasgroups_feature_stats.csv ({len(df_pg)} features)")

    if freq_tx is not None:
        df_tx = pd.DataFrame({
            'chemotype': names_tx,
            'prevalence': freq_tx,
            'info_content_bits': ic_tx,
        }).sort_values('prevalence', ascending=False)
        df_tx.to_csv(os.path.join(outdir, 'txppfas_feature_stats.csv'), index=False)
        print(f"  Saved txppfas_feature_stats.csv ({len(df_tx)} features)")


def save_comparison_summary(
    n_mol,
    X_pg, names_pg, freq_pg,
    outdir,
    n_total_pg: int = None,
    n_total_tx: int = None,
    X_tx=None, names_tx=None, freq_tx=None,
    r=None, rho=None,
):
    """Print and save scalar comparison metrics."""
    n_active_pg = (X_pg.any(axis=0)).sum()
    cov_pg      = (X_pg.any(axis=1)).mean() * 100
    mean_fp_pg  = X_pg.sum(axis=1).mean()

    rows = []
    rows.append({'metric': 'n_molecules', 'PFASGroups': n_mol,
                 'TxP_PFAS': n_mol if X_tx is not None else ''})
    rows.append({'metric': 'total_features', 'PFASGroups': n_total_pg or X_pg.shape[1],
                 'TxP_PFAS': (n_total_tx or X_tx.shape[1]) if X_tx is not None else ''})
    rows.append({'metric': 'active_features_(hit_ge1_mol)',
                 'PFASGroups': int(n_active_pg),
                 'TxP_PFAS': int((X_tx.any(axis=0)).sum()) if X_tx is not None else ''})
    rows.append({'metric': 'pct_mols_with_ge1_feature',
                 'PFASGroups': f'{cov_pg:.1f}',
                 'TxP_PFAS': f'{(X_tx.any(axis=1)).mean()*100:.1f}' if X_tx is not None else ''})
    rows.append({'metric': 'mean_prevalence_pct',
                 'PFASGroups': f'{freq_pg.mean()*100:.2f}',
                 'TxP_PFAS': f'{freq_tx.mean()*100:.2f}' if freq_tx is not None else ''})
    rows.append({'metric': 'median_prevalence_pct',
                 'PFASGroups': f'{np.median(freq_pg)*100:.2f}',
                 'TxP_PFAS': f'{np.median(freq_tx)*100:.2f}' if freq_tx is not None else ''})
    rows.append({'metric': 'mean_bits_per_mol',
                 'PFASGroups': f'{mean_fp_pg:.2f}',
                 'TxP_PFAS': f'{X_tx.sum(axis=1).mean():.2f}' if X_tx is not None else ''})
    rows.append({'metric': 'mean_IC_bits',
                 'PFASGroups': f'{information_content(freq_pg).mean():.3f}',
                 'TxP_PFAS': f'{information_content(freq_tx).mean():.3f}' if freq_tx is not None else ''})
    if r is not None:
        rows.append({'metric': 'pairwise_sim_pearson_r', 'PFASGroups': '', 'TxP_PFAS': f'{r:.3f}'})
        rows.append({'metric': 'pairwise_sim_spearman_rho', 'PFASGroups': '', 'TxP_PFAS': f'{rho:.3f}'})

    df = pd.DataFrame(rows)
    csv_path = os.path.join(outdir, 'comparison_summary.csv')
    df.to_csv(csv_path, index=False)

    print('\n── Comparison Summary ──')
    print(df.to_string(index=False))
    print(f'\n  Saved comparison_summary.csv')


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--smiles', default=None,
        help=(
            'TSV file with SMILES column (and optionally DTXSID). '
            f'Default: {DATA_DIR / "test_set_for_PFASSTRUCTv5.tsv"}'
        ),
    )
    parser.add_argument(
        '--txppfas_csv', default=None,
        help=(
            'TxP_PFAS fingerprint CSV (SI Table S2): index=DTXSID, '
            '129 chemotype columns. Without this, only PFASGroups analysis is run.'
            f'Default: {DATA_DIR / "Richard2023_SI_TableS2.csv"}'
        ),
    )
    parser.add_argument(
        '--halogens', default='F',
        help=(
            "Comma-separated halogens for PFASGroups FP (e.g. 'F' or 'F,Cl'). "
            "Use 'F' to match TxP_PFAS scope (default: F)."
        ),
    )
    parser.add_argument(
        '--outdir', default=None,
        help=(
            'Output directory for figures and CSV files. '
            f'Default: {DATA_DIR / "figures_csrml_comparison"}'
        ),
    )
    parser.add_argument(
        '--max_mols', type=int, default=None,
        help='Limit the number of molecules processed (for quick testing).',
    )
    parser.add_argument(
        '--no_pairwise', action='store_true',
        help='Skip pairwise Jaccard computation (can be slow for >200 molecules).',
    )
    args = parser.parse_args()

    # ── Defaults ─────────────────────────────────────────────────────────────
    smiles_path = args.smiles or str(DATA_DIR / 'test_set_for_PFASSTRUCTv5.tsv')
    outdir      = args.outdir or str(DATA_DIR / 'figures_csrml_comparison')
    txppfas_csv = args.txppfas_csv or str(DATA_DIR / "Richard2023_SI_TableS2.csv")
    os.makedirs(outdir, exist_ok=True)

    halogens_raw = [h.strip() for h in args.halogens.split(',')]
    halogens     = halogens_raw[0] if len(halogens_raw) == 1 else halogens_raw

    print(f'\n{"="*70}')
    print('PFASGroups vs TxP_PFAS (CSRML) Fingerprint Comparison')
    print(f'{"="*70}\n')

    # ── Load molecules ────────────────────────────────────────────────────────
    print(f'Loading molecules from: {smiles_path}')
    df = load_molecules(smiles_path)
    if args.max_mols:
        df = df.head(args.max_mols).reset_index(drop=True)
    n_mol = len(df)
    smiles_list = df['SMILES'].tolist()
    print(f'  {n_mol} valid molecules')
    if n_mol < 5:
        print('  WARNING: Very few molecules – statistical results will not be meaningful.')

    # ── Compute PFASGroups fingerprints ───────────────────────────────────────
    print(f'\nComputing PFASGroups fingerprints (halogens={args.halogens}) ...')
    X_pg_full, names_pg_full = compute_pfasgroups_fps(smiles_list, halogens=halogens)

    # Active features (hit at least once in the dataset)
    active_pg = X_pg_full.any(axis=0)
    X_pg      = X_pg_full[:, active_pg]
    names_pg  = [names_pg_full[i] for i, m in enumerate(active_pg) if m]
    freq_pg   = X_pg.mean(axis=0)

    n_total_pg = X_pg_full.shape[1]
    n_active_pg = active_pg.sum()
    cov_pg     = (X_pg.any(axis=1)).mean() * 100
    print(f'  Total groups: {n_total_pg}  |  Active: {n_active_pg}  |  '
          f'Mol coverage: {cov_pg:.1f}%  |  '
          f'Mean bits/mol: {X_pg.sum(axis=1).mean():.2f}')

    # ── Optionally load TxP_PFAS fingerprints ────────────────────────────────
    X_tx     = None
    names_tx = None
    freq_tx  = None

    if txppfas_csv:
        print(f'\nLoading TxP_PFAS fingerprints from: {txppfas_csv}')
        try:
            X_tx_full, names_tx_full, dtxsids_tx = load_txppfas_csv(txppfas_csv)
        except (ValueError, FileNotFoundError) as exc:
            print(f'  WARNING: Cannot load TxP_PFAS CSV – {exc}')
            print('  Proceeding with PFASGroups-only analysis.')
            txppfas_csv = None
            X_tx_full = names_tx_full = dtxsids_tx = None
        else:
            print(f'  TxP_PFAS matrix: {X_tx_full.shape[0]} rows × {X_tx_full.shape[1]} features')

        # Align on DTXSID if available in the SMILES file
        if 'DTXSID' in df.columns and dtxsids_tx is not None:
            shared_ids = [d for d in df['DTXSID'].tolist() if d in dtxsids_tx]
            n_shared   = len(shared_ids)
            print(f'  {n_shared} DTXSIDs overlap between molecule file and TxP_PFAS CSV')
            if n_shared == 0:
                print('  ERROR: No overlapping DTXSIDs. Cannot align fingerprints.')
                print('  Proceeding with PFASGroups-only analysis.')
                txppfas_csv = None
            else:
                if n_shared < 10:
                    print('  WARNING: Fewer than 10 overlapping molecules.')

                # Keep only overlapping molecules (both FPs aligned to same row order)
                keep_mask  = df['DTXSID'].isin(shared_ids)
                df_shared  = df[keep_mask].reset_index(drop=True)
                X_pg_match = X_pg_full[keep_mask.values, :][:, active_pg]

                tx_row_idx = {d: i for i, d in enumerate(dtxsids_tx)}
                tx_rows    = [tx_row_idx[d] for d in df_shared['DTXSID'].tolist()]
                X_tx_match = X_tx_full[tx_rows, :]

                active_tx   = X_tx_match.any(axis=0)
                X_tx        = X_tx_match[:, active_tx]
                names_tx    = [names_tx_full[i] for i, m in enumerate(active_tx) if m]
                freq_tx     = X_tx.mean(axis=0)

                # Align X_pg to the shared subset
                X_pg    = X_pg_match
                freq_pg = X_pg.mean(axis=0)
                n_mol   = len(df_shared)
                smiles_list = df_shared['SMILES'].tolist()

                print(f'  Active TxP_PFAS features: {active_tx.sum()} / {X_tx_match.shape[1]}')
                print(f'  Mol coverage: {(X_tx.any(axis=1)).mean()*100:.1f}%  |  '
                      f'Mean bits/mol: {X_tx.sum(axis=1).mean():.2f}')
                print(f'  Using {n_mol} aligned molecules for comparison')

        else:
            # No DTXSID → assume row order matches (truncate to shorter)
            min_n    = min(n_mol, X_tx_full.shape[0])
            X_tx_sub = X_tx_full[:min_n, :]
            X_pg     = X_pg[:min_n, :]
            freq_pg  = X_pg.mean(axis=0)

            active_tx = X_tx_sub.any(axis=0)
            X_tx      = X_tx_sub[:, active_tx]
            names_tx  = [names_tx_full[i] for i, m in enumerate(active_tx) if m]
            freq_tx   = X_tx.mean(axis=0)

            if min_n < n_mol:
                print(f'  WARNING: Row-order alignment truncated to {min_n} molecules.')
    else:
        print('\nNo TxP_PFAS CSV provided – running PFASGroups-only analysis.')
        print('  To enable full comparison, see the docstring at the top of this file.')

    # ── Generate figures ──────────────────────────────────────────────────────
    print('\n── Generating figures ──')

    plot_coverage(freq_pg, names_pg, outdir,
                  freq_tx=freq_tx, names_tx=names_tx)

    ic_pg, ic_tx = plot_info_content(freq_pg, names_pg, outdir,
                                     freq_tx=freq_tx, names_tx=names_tx)

    plot_lorenz(freq_pg, outdir, freq_tx=freq_tx)

    plot_pca(X_pg, outdir, X_tx=X_tx)

    r = rho = None
    if X_tx is not None and not args.no_pairwise and n_mol <= 600:
        print(f'  Computing pairwise Jaccard ({n_mol}×{n_mol}) ...')
        sim_pg = pairwise_jaccard(X_pg.astype(bool))
        sim_tx = pairwise_jaccard(X_tx.astype(bool))
        plot_jaccard_dist(sim_pg, sim_tx, outdir)
        r, rho = plot_cross_similarity(sim_pg, sim_tx, outdir)
        print(f'  Pairwise sim:  Pearson r = {r:.3f},  Spearman ρ = {rho:.3f}')

    elif X_tx is not None and n_mol > 600 and not args.no_pairwise:
        print(f'  Skipping pairwise sim (n={n_mol} > 600). Use --max_mols 500 or --no_pairwise.')

    # ── Feature-level statistics CSV ─────────────────────────────────────────
    print('\n── Saving feature statistics ──')
    save_feature_summary(freq_pg, names_pg, ic_pg, outdir,
                         freq_tx=freq_tx, names_tx=names_tx, ic_tx=ic_tx)

    # ── Scalar summary ────────────────────────────────────────────────────────
    save_comparison_summary(
        n_mol, X_pg, names_pg, freq_pg, outdir,
        n_total_pg=n_total_pg,
        n_total_tx=(X_tx_full.shape[1] if args.txppfas_csv and X_tx is not None else None),
        X_tx=X_tx, names_tx=names_tx, freq_tx=freq_tx,
        r=r, rho=rho,
    )

    # ── Top-20 features by prevalence ─────────────────────────────────────────
    order_pg = np.argsort(freq_pg)[::-1]
    print('\n── Top 20 PFASGroups features by prevalence ──')
    for rank in range(min(20, len(freq_pg))):
        i = order_pg[rank]
        print(f'  {rank+1:2d}. {freq_pg[i]*100:5.1f}%  {names_pg[i]}')

    if freq_tx is not None:
        order_tx = np.argsort(freq_tx)[::-1]
        print('\n── Top 20 TxP_PFAS chemotypes by prevalence ──')
        for rank in range(min(20, len(freq_tx))):
            i = order_tx[rank]
            print(f'  {rank+1:2d}. {freq_tx[i]*100:5.1f}%  {names_tx[i]}')

    print(f'\nAll outputs saved to: {outdir}\n')
    print('Done.\n')


# ══════════════════════════════════════════════════════════════════════════════
# Pytest entry points
# ══════════════════════════════════════════════════════════════════════════════

def test_jaccard_similarity_basic():
    """Jaccard similarity: known values and degenerate all-zero case."""
    a = np.array([1, 0, 1, 1, 0], dtype=bool)
    b = np.array([1, 1, 1, 0, 0], dtype=bool)
    # intersection=2, union=4
    assert jaccard_similarity(a, b) == pytest.approx(2 / 4, abs=1e-9)
    # both zero vectors → 0
    assert jaccard_similarity(np.zeros(5, dtype=bool), np.zeros(5, dtype=bool)) == 0.0
    # identical → 1
    assert jaccard_similarity(a, a) == pytest.approx(1.0, abs=1e-9)


def test_information_content_known_values():
    """Shannon IC: max at p=0.5 (1 bit), near-zero at extreme probabilities."""
    freq = np.array([0.5, 0.5])
    ic = information_content(freq)
    assert ic == pytest.approx([1.0, 1.0], abs=1e-6)

    freq_low = np.array([1e-9, 1.0 - 1e-9])
    ic_low = information_content(freq_low)
    assert ic_low[0] < 0.01


_KNOWN_PFAS = [
    ('FC(F)(F)C(=O)O',           'perfluoroacetic acid',  True),
    ('FC(F)(F)C(F)(F)C(F)(F)F', 'perfluorobutane',       True),
    ('CCCC',                     'n-butane (not PFAS)',   False),
]


@pytest.mark.parametrize('smiles,name,expect_nonzero', _KNOWN_PFAS)
def test_pfasgroups_fp_known_molecules(smiles, name, expect_nonzero):
    """PFASGroups FP should be non-zero for PFAS and zero for non-PFAS."""
    from PFASGroups import generate_fingerprint
    fp, _ = generate_fingerprint(
        smiles, representation='vector', count_mode='binary',
        halogens='F', saturation=None,
    )
    nonzero = any(v > 0 for v in fp)
    assert nonzero == expect_nonzero, (
        f"Expected {'non-zero' if expect_nonzero else 'zero'} FP "
        f"for {name} ({smiles})"
    )


def test_compute_pfasgroups_fps_shape():
    """compute_pfasgroups_fps returns a matrix with the right shape and dtype."""
    smiles = ['FC(F)(F)C(=O)O', 'FC(F)(F)C(F)(F)F', 'CCCC']
    X, names = compute_pfasgroups_fps(smiles, verbose=False)
    assert X.shape[0] == 3
    assert X.shape[1] == len(names)
    assert X.dtype == np.uint8
    # PFAS molecules should have at least one bit set
    assert X[0].any() or X[1].any()


def test_load_txppfas_csv_skips_metadata(tmp_path):
    """load_txppfas_csv correctly skips metadata rows to find the DTXSID header."""
    content = (
        "Sheet description,some text\n"
        "More metadata,\n"
        "DTXSID,feat_A,feat_B\n"
        "DTXSID001,1,0\n"
        "DTXSID002,0,1\n"
    )
    csv_file = tmp_path / 'test_txp.csv'
    csv_file.write_text(content)
    X, chemotypes, dtxsids = load_txppfas_csv(str(csv_file))
    assert X.shape == (2, 2)
    assert chemotypes == ['feat_A', 'feat_B']
    assert dtxsids == ['DTXSID001', 'DTXSID002']


def test_load_txppfas_csv_no_dtxsid_raises(tmp_path):
    """load_txppfas_csv raises ValueError with helpful message when no DTXSID row exists."""
    csv_file = tmp_path / 'bad.csv'
    csv_file.write_text("Table S1.,description\nMore text,here\n")
    with pytest.raises(ValueError, match='DTXSID'):
        load_txppfas_csv(str(csv_file))


def _txppfas_csv_has_actual_data(path: Path) -> bool:
    """Return True only if path exists and contains a valid DTXSID header row."""
    if not path.exists():
        return False
    return _find_dtxsid_header_row(str(path)) >= 0


@pytest.mark.skipif(
    not _txppfas_csv_has_actual_data(DATA_DIR / 'Richard2023_SI_TableS2.csv'),
    reason='Richard2023_SI_TableS2.csv is absent or is the metadata placeholder',
)
def test_txppfas_csv_loads_full_matrix():
    """If the actual TxP_PFAS CSV is present it should load ≥100 rows × ≥10 cols."""
    X, chemotypes, dtxsids = load_txppfas_csv(str(DATA_DIR / 'Richard2023_SI_TableS2.csv'))
    assert X.shape[0] >= 100, f'Expected ≥100 rows, got {X.shape[0]}'
    assert X.shape[1] >= 10,  f'Expected ≥10 chemotype columns, got {X.shape[1]}'
    assert np.isin(X.flatten(), [0, 1]).all(), 'All fingerprint values should be 0 or 1'


if __name__ == '__main__':
    main()
