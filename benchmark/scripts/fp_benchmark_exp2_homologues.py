"""
Fingerprint Benchmark – Experiment 2: Homologue Series
=======================================================

Generates PFCA, PFSA, and PFAL homologue series (C2–C8, i.e. 2–8 fluorinated
carbons in the perfluoroalkyl chain) using generate_homologues from a C8
seed compound, yielding 21 structures (7 per series).

Fingerprints are computed in three modes:
  - binary    : 1 if group present, 0 otherwise
  - count     : number of matched substructures
  - max_component : largest perfluoroalkyl component size matched

Pairwise Tanimoto–Jaccard similarity is computed for each mode:
  J(a,b) = Σ min(aᵢ,bᵢ) / Σ max(aᵢ,bᵢ)

which collapses to standard set-Jaccard for binary vectors.

Outputs
-------
fp_homologue_jaccard.pdf / .png
    Three-panel pairwise Jaccard heatmap (binary / count / max-component).
fp_homologue_results.csv
    Fingerprint matrix (all three modes, stacked) with metadata.

Usage
-----
python fp_benchmark_exp2_homologues.py [--outdir PATH] [--article_imgs PATH]
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
import matplotlib.colors as mcolors
from scipy.stats import pearsonr

# ── Path setup ────────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
REPO_ROOT    = SCRIPT_DIR.parent.parent
BENCHMARK_DIR = SCRIPT_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

_ARTICLE_IMGS = REPO_ROOT.parent / 'overleaf_PFASgroups_article' / 'imgs'

warnings.filterwarnings('ignore')

# ── Series definition ─────────────────────────────────────────────────────────
# Each C8 seed: n_cf2 = 7  →  CF3-(CF2)7-[FG]  (8 fluorinated carbons)
SEEDS = {
    'PFCA': 'FC(F)(F)' + 'C(F)(F)' * 7 + 'C(=O)O',    # C8 perfluoroalkyl carboxylic acid
    'PFSA': 'FC(F)(F)' + 'C(F)(F)' * 7 + 'S(=O)(=O)O', # C8 perfluoroalkyl sulfonic acid  (PFOS)
    'PFAL': 'FC(F)(F)' + 'C(F)(F)' * 7 + 'O',          # C8 perfluoroalkyl alcohol
}

# Colour per series
SERIES_COLOURS = {'PFCA': '#1976D2', 'PFSA': '#D32F2F', 'PFAL': '#388E3C'}


# ── SMILES builders (fallback, used when generate_homologues cannot be loaded) ──

def build_pfas(fg_suffix: str, n_cf2: int) -> str:
    """Linear perfluoroalkyl chain: CF3-(CF2)n_cf2-[fg]."""
    return 'FC(F)(F)' + 'C(F)(F)' * n_cf2 + fg_suffix


FG_SUFFIX = {'PFCA': 'C(=O)O', 'PFSA': 'S(=O)(=O)O', 'PFAL': 'O'}


# ── Chain-length counting ─────────────────────────────────────────────────────

def count_f_carbons(mol) -> int:
    """Count carbons bearing ≥2 fluorine atoms (CF₂ and CF₃ groups)."""
    return sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6
        and sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 9) >= 2
    )


# ── Homologue generation ──────────────────────────────────────────────────────

def build_homologue_series(series_name: str, verbose: bool = True) -> list[dict]:
    """
    Build the C2–C8 homologue series for *series_name* (PFCA, PFSA, or PFAL).

    Strategy
    --------
    1. Try to use generate_homologues on the C8 seed to obtain the shorter
       homologues; combine with the C8 parent.
    2. Fall back to manual SMILES construction if generate_homologues is
       unavailable or returns an incomplete set.

    Returns
    -------
    list of dicts with keys: series, cn, n_cf2, smiles
    Sorted by ascending chain length (C2 first, C8 last).
    """
    from rdkit import Chem

    seed_smi  = SEEDS[series_name]
    fg_suffix = FG_SUFFIX[series_name]
    records   = []

    # ── Attempt generate_homologues ───────────────────────────────────────
    generated = {}
    try:
        from PFASGroups import generate_homologues as _gen_hom
        homologue_dict = _gen_hom(seed_smi)
        for inchikey, formula_mol in homologue_dict.items():
            for formula, mol in formula_mol.items():
                smi = Chem.MolToSmiles(mol)
                cn  = count_f_carbons(mol)
                if cn >= 2:
                    generated[cn] = smi
        if verbose:
            print(f"  {series_name}: generate_homologues returned {len(generated)} homologues "
                  f"(chain lengths: {sorted(generated.keys())})")
    except Exception as exc:
        if verbose:
            print(f"  {series_name}: generate_homologues unavailable ({exc}); using manual build")

    # Merge: manual build for C2–C8 (n_cf2 1→7), override with generate_homologues where available
    for n_cf2 in range(1, 8):   # n_cf2 = 1..7  →  Cn where n = n_cf2+1 (CF3 + n_cf2 × CF2)
        cn  = n_cf2 + 1         # total fluorinated carbons = 1 CF3 + n_cf2 CF2
        manual_smi = build_pfas(fg_suffix, n_cf2)

        if cn in generated:
            smi = generated[cn]
        else:
            smi = manual_smi

        # Validate
        mol_check = Chem.MolFromSmiles(smi)
        if mol_check is None:
            if verbose:
                print(f"  WARNING: invalid SMILES for {series_name} C{cn}: {smi!r}; using manual")
            smi = manual_smi
            mol_check = Chem.MolFromSmiles(smi)

        records.append({
            'series': series_name,
            'cn'    : cn,
            'n_cf2' : n_cf2,
            'smiles': smi,
        })

    return records    # C2 first, C8 last


# ── Fingerprint computation ───────────────────────────────────────────────────

def compute_fps_all_modes(smiles_list: list) -> dict:
    """Return dict {'binary': X, 'count': X, 'max_component': X} (float32 arrays)."""
    from HalogenGroups import generate_fingerprint

    # Probe dims
    _fp, _info = generate_fingerprint(
        'FC(F)(F)C(=O)O',
        representation='vector', count_mode='binary', halogens='F', saturation='per',
    )
    n_groups, group_names = len(_fp), _info['group_names']
    n = len(smiles_list)

    result = {mode: np.zeros((n, n_groups), dtype=np.float32)
              for mode in ('binary', 'count', 'max_component')}

    for mode in ('binary', 'count', 'max_component'):
        for i, smi in enumerate(smiles_list):
            try:
                fp, _ = generate_fingerprint(
                    smi, representation='vector', count_mode=mode,
                    halogens='F', saturation='per',
                )
                result[mode][i] = fp
            except Exception:
                pass   # leave as zeros

    return result, group_names


# ── Jaccard ───────────────────────────────────────────────────────────────────

def generalised_jaccard(a: np.ndarray, b: np.ndarray) -> float:
    """Tanimoto-Jaccard: Σmin(aᵢ,bᵢ) / Σmax(aᵢ,bᵢ). Reduces to set-Jaccard for binary."""
    num = np.sum(np.minimum(a, b))
    den = np.sum(np.maximum(a, b))
    return float(num / den) if den > 0 else 0.0


def pairwise_jaccard(X: np.ndarray) -> np.ndarray:
    """n×n pairwise generalised Jaccard matrix."""
    n = X.shape[0]
    J = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i, n):
            v = generalised_jaccard(X[i], X[j])
            J[i, j] = J[j, i] = v
    return J


# ── Utility ───────────────────────────────────────────────────────────────────

def _save(fig, outdir: Path, stem: str, verbose: bool = True) -> None:
    for ext in ('.pdf', '.png'):
        fig.savefig(outdir / (stem + ext), dpi=150, bbox_inches='tight')
    if verbose:
        print(f"  Saved {stem}.pdf / .png")
    plt.close(fig)


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_jaccard_heatmaps(
    Js: dict, records: list, outdir: Path, verbose: bool = True
) -> None:
    """Three-panel Jaccard heatmap (binary / count / max-component)."""
    n = len(records)
    labels = [f"{r['series']} C{r['cn']}" for r in records]
    # Series boundary lines
    series_names = [r['series'] for r in records]
    boundaries = [i + 0.5 for i in range(n - 1) if series_names[i] != series_names[i + 1]]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    mode_titles = {
        'binary'       : 'Binary mode\n(set-Jaccard)',
        'count'        : 'Count mode\n(generalised Jaccard)',
        'max_component': 'Max-component mode\n(generalised Jaccard)',
    }

    for ax, (mode, J) in zip(axes, Js.items()):
        im = ax.imshow(J, vmin=0, vmax=1, cmap='Blues', aspect='auto')
        ax.set_xticks(range(n))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_title(mode_titles[mode], fontsize=11, fontweight='bold')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Jaccard similarity')
        # Series boundary lines
        for b in boundaries:
            ax.axhline(b, color='red', lw=1.2, alpha=0.7)
            ax.axvline(b, color='red', lw=1.2, alpha=0.7)

    fig.suptitle(
        'Experiment 2: Pairwise Jaccard Similarity – Homologue Series (C2–C8)\n'
        'Red lines separate PFCA / PFSA / PFAL series',
        fontsize=13, fontweight='bold',
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, outdir, 'fp_homologue_jaccard', verbose=verbose)


# ── Main ──────────────────────────────────────────────────────────────────────

def run(outdir: Path, article_imgs: Path | None, verbose: bool = True) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Build 21-compound dataset ─────────────────────────────────────────
    if verbose:
        print("\nBuilding homologue series …")
    all_records = []
    for series in ('PFCA', 'PFSA', 'PFAL'):
        recs = build_homologue_series(series, verbose=verbose)
        all_records.extend(recs)
        if verbose:
            print(f"  {series}: {len(recs)} compounds (C{recs[0]['cn']}–C{recs[-1]['cn']})")

    n_total = len(all_records)
    if verbose:
        print(f"  Total: {n_total} compounds")

    smiles_list = [r['smiles'] for r in all_records]
    chain_ns    = [r['cn']    for r in all_records]
    series_list = [r['series'] for r in all_records]

    # ── Compute fingerprints in all three modes ───────────────────────────
    if verbose:
        print("Computing fingerprints (binary, count, max_component) …")
    fps_dict, group_names = compute_fps_all_modes(smiles_list)

    # ── Compute pairwise Jaccard for each mode ────────────────────────────
    if verbose:
        print("Computing pairwise Jaccard matrices …")
    Js = {mode: pairwise_jaccard(fps_dict[mode]) for mode in fps_dict}

    # ── Summary statistics ────────────────────────────────────────────────
    if verbose:
        for mode, J in Js.items():
            # Within-series: pairs where series is the same
            within_vals = [
                J[i, j]
                for i in range(n_total) for j in range(i + 1, n_total)
                if series_list[i] == series_list[j]
            ]
            across_vals = [
                J[i, j]
                for i in range(n_total) for j in range(i + 1, n_total)
                if series_list[i] != series_list[j]
            ]
            print(f"  [{mode}]  within-series Jaccard: "
                  f"mean={np.mean(within_vals):.3f} ± {np.std(within_vals):.3f}  "
                  f"| across: mean={np.mean(across_vals):.3f} ± {np.std(across_vals):.3f}")

        # Pearson r: chain-length difference vs Jaccard for max_component within-series
        chain_arr    = np.array(chain_ns)
        maxcomp_J    = Js['max_component']
        delta_n_mc   = []
        sim_vals_mc  = []
        delta_n_cnt  = []
        sim_vals_cnt = []
        for i in range(n_total):
            for j in range(i + 1, n_total):
                if series_list[i] == series_list[j]:
                    dn = abs(chain_arr[i] - chain_arr[j])
                    delta_n_mc.append(dn)
                    sim_vals_mc.append(maxcomp_J[i, j])
                    delta_n_cnt.append(dn)
                    sim_vals_cnt.append(Js['count'][i, j])

        r_mc, pval_mc = (pearsonr(delta_n_mc, sim_vals_mc)
                         if len(set(sim_vals_mc)) > 1 else (float('nan'), float('nan')))
        r_cnt, pval_cnt = (pearsonr(delta_n_cnt, sim_vals_cnt)
                           if len(set(sim_vals_cnt)) > 1 else (float('nan'), float('nan')))
        print(f"\n  Max-component within-series: Pearson r(Δn, Jaccard) = {r_mc:.3f} (p={pval_mc:.3g})")
        print(f"  Count-mode within-series:    Pearson r(Δn, Jaccard) = {r_cnt:.3f}")
        print(f"  (max-component encodes chain length; count mode is chain-length-invariant)")
        print(f"\n  Max-component within-series Jaccard: {np.mean(sim_vals_mc):.3f} ± {np.std(sim_vals_mc):.3f}")
        print(f"  Max-component across-series Jaccard: ", end='')
        across_mc = [maxcomp_J[i, j]
                     for i in range(n_total) for j in range(i+1, n_total)
                     if series_list[i] != series_list[j]]
        print(f"{np.mean(across_mc):.3f} ± {np.std(across_mc):.3f}")

    # ── Save CSV ───────────────────────────────────────────────────────────
    rows = []
    for i, rec in enumerate(all_records):
        row = {**rec}
        for mode in fps_dict:
            for j, gname in enumerate(group_names):
                row[f'{mode}__{gname}'] = fps_dict[mode][i, j]
        rows.append(row)
    df = pd.DataFrame(rows)
    csv_path = outdir / 'fp_homologue_results.csv'
    df.to_csv(csv_path, index=False)
    if verbose:
        print(f"  Saved {csv_path.name}")

    # ── Generate figure ───────────────────────────────────────────────────
    if verbose:
        print("Generating Jaccard heatmap …")
    plot_jaccard_heatmaps(Js, all_records, outdir, verbose=verbose)

    # ── Copy PDFs to article imgs/ ─────────────────────────────────────────
    if article_imgs is not None and article_imgs.is_dir():
        import shutil
        for stem in ('fp_homologue_jaccard',):
            for ext in ('.pdf', '.png'):
                src = outdir / (stem + ext)
                dst = article_imgs / (stem + ext)
                if src.exists():
                    shutil.copy2(src, dst)
        if verbose:
            print(f"  Copied figures to article imgs: {article_imgs}")

    if verbose:
        print("\n=== Experiment 2 summary ===")
        print(f"  Total compounds: {n_total}")
        for series in ('PFCA', 'PFSA', 'PFAL'):
            n_s = sum(1 for r in all_records if r['series'] == series)
            print(f"  {series}: {n_s} compounds")
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
