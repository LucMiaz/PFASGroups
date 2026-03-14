#!/usr/bin/env python3
"""
Timing performance comparison: PFASGroups vs PFAS-Atlas
on the OECD PFAS dataset (CSV) and the clinventory database.

For each dataset the script:
  • Classifies every molecule with PFASGroups (parse_mol) and PFAS-Atlas
  • Records per-molecule wall-clock time (ms, via time.perf_counter)
  • Computes molecular complexity metrics from RDKit + PFASGroups:
      – n_heavy_atoms (size)
      – fluorine_ratio  (n_F / n_heavy)
      – branching_index (avg excess degree of C atoms beyond 2)
      – n_rings, n_rot_bonds
      – component sizes reported by PFASGroups
  • Saves results to data/oecd_clinventory_timing_<TIMESTAMP>.json
  • Saves all figures to imgs/ (PDF + PNG)

Usage (from the benchmark/ directory, chem conda env unless --no-atlas is set):
    python scripts/compare_oecd_clinventory_timing.py [--no-clinventory] [--no-atlas]
                                                       [--limit N] [--oecd-csv PATH]
"""

import argparse
import csv
import getpass
import json
import math
import os
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
SCRIPT_DIR   = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent
REPO_ROOT    = BENCHMARK_DIR.parent
ATLAS_DIR    = REPO_ROOT.parent / "PFAS-atlas"

for _p in (str(REPO_ROOT), str(ATLAS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ---------------------------------------------------------------------------
# Optional dependencies
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    RDKIT = True
except ImportError as _e:
    print(f"ERROR: rdkit not available: {_e}")
    sys.exit(1)

try:
    from PFASGroups import parse_mol
    from PFASGroups.core import rdkit_disable_log
    rdkit_disable_log()
    PFASGROUPS_AVAILABLE = True
except ImportError as _e:
    print(f"ERROR: PFASGroups not available: {_e}")
    sys.exit(1)

try:
    from classification_helper.classify_pfas import classify_pfas_molecule as _cp
except ImportError:
    try:
        sys.path.insert(0, str(ATLAS_DIR / "classification_helper"))
        from classify_pfas import classify_pfas_molecule as _cp
    except ImportError:
        _cp = None  # handled in main

def _atlas_classify(smiles: str):
    """Return (class1, is_pfas) or raise."""
    result = _cp(smiles)  # type: ignore[misc]
    class1 = result[0] if result else "Unknown"
    is_pfas = class1 not in ("Not PFAS", "Unknown", "", None)
    return class1, is_pfas

try:
    import psycopg2
    PSYCOPG2 = True
except ImportError:
    PSYCOPG2 = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.ticker import MaxNLocator
    import seaborn as sns
    plt.style.use("seaborn-v0_8-whitegrid")
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False
    print("WARNING: matplotlib/seaborn not available — plots will be skipped")

try:
    import numpy as np
    NUMPY = True
except ImportError:
    NUMPY = False


# ---------------------------------------------------------------------------
# Colour palette (stdlib only, no pyyaml)
# ---------------------------------------------------------------------------
def _load_palette() -> List[str]:
    import re
    defaults = ["#E15D0B", "#306DBA", "#9D206C", "#51127C"]
    try:
        _p = REPO_ROOT / "PFASGroups" / "data" / "color_scheme.yaml"
        found = re.findall(r'"(#[0-9A-Fa-f]{6})"', _p.read_text())
        if len(found) >= 4:
            return found[:4]
    except Exception:
        pass
    return defaults

_PALETTE = _load_palette()
_C0, _C1, _C2, _C3 = _PALETTE          # orange, blue, magenta, dark-purple

COLORS = {
    "PFASGroups": _C0,
    "PFAS-Atlas": _C1,
    "OECD":       _C2,
    "clinventory": _C3,
}

DATASET_MARKERS = {"OECD": "o", "clinventory": "s"}
DATASET_LS      = {"OECD": "-", "clinventory": "--"}

STYLE = {
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "font.size": 10,
    "figure.dpi": 150,
}

def apply_style():
    if MATPLOTLIB:
        plt.rcParams.update(STYLE)


# ---------------------------------------------------------------------------
# Bracket helpers (reused from compare_clinventory_classifiers.py)
# ---------------------------------------------------------------------------
BRACKETS = [
    "tiny (<10)", "small (10–19)", "medium (20–34)",
    "large (35–59)", "very large (≥60)",
]

def atom_bracket(n: Optional[int]) -> str:
    if n is None: return "unknown"
    if n < 10:    return "tiny (<10)"
    if n < 20:    return "small (10–19)"
    if n < 35:    return "medium (20–34)"
    if n < 60:    return "large (35–59)"
    return "very large (≥60)"

def stats(vals: list) -> dict:
    if not vals:
        return {}
    n = len(vals)
    s = sorted(vals)
    mu = sum(vals) / n
    var = sum((v - mu) ** 2 for v in vals) / n
    return {
        "n": n, "mean": mu, "std": math.sqrt(var),
        "median": s[n // 2], "min": s[0], "max": s[-1],
        "p25":  s[n // 4]        if n >= 4 else s[0],
        "p75":  s[3 * n // 4]    if n >= 4 else s[-1],
        "p95":  s[int(0.95 * n)] if n >= 2 else s[-1],
    }


# ---------------------------------------------------------------------------
# Molecular complexity metrics
# ---------------------------------------------------------------------------
def compute_complexity(mol) -> dict:
    """Compute complexity metrics from an RDKit mol (without PFASGroups)."""
    heavy = mol.GetNumHeavyAtoms()
    # elemental counts
    n_f  = sum(a.GetAtomicNum() == 9  for a in mol.GetAtoms())
    n_cl = sum(a.GetAtomicNum() == 17 for a in mol.GetAtoms())
    n_br = sum(a.GetAtomicNum() == 35 for a in mol.GetAtoms())
    n_i  = sum(a.GetAtomicNum() == 53 for a in mol.GetAtoms())
    n_hal = n_f + n_cl + n_br + n_i
    # topological
    n_rings   = mol.GetRingInfo().NumRings()
    n_rot     = rdMolDescriptors.CalcNumRotatableBonds(mol)
    # branching index: per C-atom excess over 2 (i.e. each branch beyond linear chain)
    c_atoms   = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    n_carbons = len(c_atoms)
    branching = (
        sum(max(0, a.GetDegree() - 2) for a in c_atoms) / n_carbons
        if n_carbons > 0 else 0.0
    )
    return {
        "n_heavy_atoms":   heavy,
        "n_fluorine":      n_f,
        "n_halogens":      n_hal,
        "fluorine_ratio":  n_f / heavy if heavy > 0 else 0.0,
        "n_rings":         n_rings,
        "n_rot_bonds":     n_rot,
        "branching_index": branching,
        "n_carbons":       n_carbons,
    }


def extract_pfasgroups_complexity(pg_result: dict) -> dict:
    """Extract component-size info from a parse_mol result dict."""
    all_comp_sizes: List[int] = []
    all_eccentricities: List[float] = []
    n_groups = 0
    for match in pg_result.get("matches", []):
        if match.get("type") == "HalogenGroup":
            n_groups += 1
            all_comp_sizes.extend(match.get("components_sizes") or [])
            ecc = match.get("mean_eccentricity")
            if ecc is not None:
                all_eccentricities.append(ecc)
    return {
        "n_pfas_groups":       n_groups,
        "n_components_total":  len(all_comp_sizes),
        "max_component_size":  max(all_comp_sizes) if all_comp_sizes else 0,
        "mean_component_size": (sum(all_comp_sizes) / len(all_comp_sizes)
                                if all_comp_sizes else 0.0),
        "mean_eccentricity":   (sum(all_eccentricities) / len(all_eccentricities)
                                if all_eccentricities else 0.0),
    }


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------
def run_pfasgroups(mol) -> Tuple[dict, float, Optional[str]]:
    """Return (result_dict, elapsed_ms, error_str)."""
    try:
        t0 = time.perf_counter()
        result = parse_mol(mol, include_PFAS_definitions=False)
        elapsed = (time.perf_counter() - t0) * 1000.0
        return result, elapsed, None
    except Exception as exc:
        return {}, 0.0, str(exc)


def run_atlas(smiles: str) -> Tuple[str, bool, float, Optional[str]]:
    """Return (class1, is_pfas, elapsed_ms, error_str)."""
    if _cp is None:
        return "unavailable", False, 0.0, "atlas not installed"
    try:
        t0 = time.perf_counter()
        class1, is_pfas = _atlas_classify(smiles)
        elapsed = (time.perf_counter() - t0) * 1000.0
        return class1, is_pfas, elapsed, None
    except Exception as exc:
        return "error", False, 0.0, str(exc)


# ---------------------------------------------------------------------------
# Dataset loaders
# ---------------------------------------------------------------------------
def load_oecd_csv(csv_path: Path, limit: Optional[int] = None) -> List[dict]:
    """Load molecules from the OECD PFAS CSV (SMILES column)."""
    rows = []
    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            smiles = row.get("SMILES", "").strip()
            if smiles in ("", "-"):
                continue
            rows.append({
                "id":    row.get("DTXSID", f"OECD_{len(rows)}"),
                "name":  row.get("PREFERRED_NAME", ""),
                "smiles": smiles,
                "formula": row.get("MOLECULAR_FORMULA", ""),
            })
            if limit and len(rows) >= limit:
                break
    return rows


def fetch_clinventory(db_params: dict, limit: Optional[int],
                      all_mols: bool) -> List[dict]:
    """Fetch molecules from the clinventory PostgreSQL database."""
    if not PSYCOPG2:
        raise RuntimeError("psycopg2 not installed")
    conn = psycopg2.connect(**db_params)
    where = "WHERE smiles IS NOT NULL AND smiles != ''" if all_mols else (
        "WHERE smiles IS NOT NULL AND smiles != '' "
        "AND (smiles LIKE '%F%' OR smiles LIKE '%Cl%' "
        "     OR smiles LIKE '%Br%' OR smiles LIKE '%I%')"
    )
    lim = f"LIMIT {limit}" if limit else ""
    sql = f"SELECT id, smiles, formula FROM molecules {where} ORDER BY id {lim}"
    with conn.cursor() as cur:
        cur.execute(sql)
        rows = [{"id": r[0], "name": "", "smiles": r[1], "formula": r[2] or ""}
                for r in cur.fetchall()]
    conn.close()
    return rows


# ---------------------------------------------------------------------------
# Main classification loop for one dataset
# ---------------------------------------------------------------------------
def classify_dataset(
    molecules: List[dict],
    dataset_name: str,
    atlas_enabled: bool,
) -> List[dict]:
    """Run PFASGroups and optionally PFAS-Atlas on every molecule.

    Returns a list of per-molecule result dicts.
    """
    records = []
    n = len(molecules)
    n_errors_pg = 0
    n_errors_at = 0
    report_every = max(1, n // 20)

    print(f"  Processing {n:,} molecules …")
    for i, mol_info in enumerate(molecules, 1):
        if i % report_every == 0 or i == n:
            print(f"    {i:,}/{n:,} …", end="\r")

        smiles = mol_info["smiles"]
        base = {
            "id":      mol_info["id"],
            "name":    mol_info.get("name", ""),
            "smiles":  smiles,
            "formula": mol_info.get("formula", ""),
            "dataset": dataset_name,
        }

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            base.update({
                "error": "invalid_smiles",
                "pg_time_ms": None, "atlas_time_ms": None,
            })
            records.append(base)
            continue

        # Complexity metrics
        base.update(compute_complexity(mol))

        # PFASGroups
        pg_result, pg_ms, pg_err = run_pfasgroups(mol)
        if pg_err:
            n_errors_pg += 1
        pg_comp = extract_pfasgroups_complexity(pg_result)
        base.update({
            "pg_time_ms":          pg_ms,
            "pg_error":            pg_err,
            "pg_is_pfas":          pg_comp["n_pfas_groups"] > 0,
            "pg_n_groups":         pg_comp["n_pfas_groups"],
            "pg_n_components":     pg_comp["n_components_total"],
            "pg_max_comp_size":    pg_comp["max_component_size"],
            "pg_mean_comp_size":   pg_comp["mean_component_size"],
            "pg_mean_eccentricity": pg_comp["mean_eccentricity"],
        })

        # PFAS-Atlas (optional)
        if atlas_enabled and _cp is not None:
            canon = Chem.MolToSmiles(mol)
            at_class1, at_is_pfas, at_ms, at_err = run_atlas(canon)
            if at_err:
                n_errors_at += 1
            base.update({
                "atlas_time_ms": at_ms,
                "atlas_error":   at_err,
                "atlas_is_pfas": at_is_pfas,
                "atlas_class1":  at_class1,
            })
        else:
            base.update({
                "atlas_time_ms": None,
                "atlas_error":   "disabled",
                "atlas_is_pfas": None,
                "atlas_class1":  None,
            })

        records.append(base)

    print()   # end the \r line
    print(f"    PFASGroups errors:  {n_errors_pg}")
    if atlas_enabled:
        print(f"    PFAS-Atlas errors: {n_errors_at}")
    return records


# ---------------------------------------------------------------------------
# Aggregate helpers
# ---------------------------------------------------------------------------
def timing_by_bracket(records: List[dict], tool_key: str) -> dict:
    """Group timing values by atom-count bracket; return stats-dict per bracket."""
    groups: dict = defaultdict(list)
    for r in records:
        t = r.get(tool_key)
        if t is None or t == 0.0:
            continue
        b = atom_bracket(r.get("n_heavy_atoms"))
        groups[b].append(t)
    return {b: {**stats(v), "median_ms": stats(v).get("median", 0)} for b, v in groups.items()}


# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------
def save_fig(fig, name: str, imgs_dir: Path):
    imgs_dir.mkdir(exist_ok=True, parents=True)
    for ext in ("pdf", "png"):
        fig.savefig(imgs_dir / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: imgs/{name}.pdf / .png")


# ---------------------------------------------------------------------------
# ── TIMING PLOTS ────────────────────────────────────────────────────────────
# ---------------------------------------------------------------------------
def plot_timing_box(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """Box plot of per-molecule timing for PFASGroups and Atlas, per dataset."""
    if not MATPLOTLIB:
        return
    apply_style()
    n_ds = len(datasets)
    fig, axes = plt.subplots(1, n_ds, figsize=(5 * n_ds, 4), sharey=True)
    if n_ds == 1:
        axes = [axes]

    for ax, (ds_name, records) in zip(axes, datasets.items()):
        pg_times    = [r["pg_time_ms"]    for r in records if r.get("pg_time_ms")]
        atlas_times = [r.get("atlas_time_ms") for r in records
                       if r.get("atlas_time_ms") is not None and r.get("atlas_time_ms") != 0]

        data, labels, colors = [], [], []
        if pg_times:
            data.append(pg_times);    labels.append("PFASGroups"); colors.append(_C0)
        if atlas_times:
            data.append(atlas_times); labels.append("PFAS-Atlas");  colors.append(_C1)
        if not data:
            ax.set_title(ds_name); continue

        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True,
                        showfliers=False, widths=0.5)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color); patch.set_alpha(0.75)
        ax.set_ylabel("Time per molecule (ms)" if ax is axes[0] else "")
        ax.set_title(ds_name)
        ax.yaxis.set_major_locator(MaxNLocator(5))

    fig.suptitle("Per-molecule classification timing", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_timing_box", imgs_dir)


def plot_timing_cdf(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """CDF of per-molecule timing on log scale for all tool×dataset combos."""
    if not MATPLOTLIB or not NUMPY:
        return
    apply_style()
    fig, ax = plt.subplots(figsize=(8, 4))

    linestyles = {"OECD": "-", "clinventory": "--"}
    for ds_name, records in datasets.items():
        ls = linestyles.get(ds_name, "-")
        for tool, key, color in [
            ("PFASGroups", "pg_time_ms",    _C0),
            ("PFAS-Atlas", "atlas_time_ms", _C1),
        ]:
            times = [r[key] for r in records
                     if r.get(key) is not None and r.get(key, 0) > 0]
            if not times:
                continue
            t_sorted = np.sort(times)
            cdf = np.arange(1, len(t_sorted) + 1) / len(t_sorted)
            ax.plot(t_sorted, cdf, color=color, ls=ls, lw=1.8,
                    label=f"{tool} / {ds_name}")

    ax.set_xscale("log")
    ax.set_xlabel("Time per molecule (ms, log scale)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("CDF of per-molecule classification time")
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    save_fig(fig, "oecd_clinventory_timing_cdf", imgs_dir)


def plot_timing_by_bracket(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """Grouped bar chart of median timing by size bracket and dataset."""
    if not MATPLOTLIB:
        return
    apply_style()
    n_ds = len(datasets)
    fig, axes = plt.subplots(1, n_ds, figsize=(7 * n_ds, 4), sharey=False)
    if n_ds == 1:
        axes = [axes]

    for ax, (ds_name, records) in zip(axes, datasets.items()):
        pg_bkt    = timing_by_bracket(records, "pg_time_ms")
        atlas_bkt = timing_by_bracket(records, "atlas_time_ms")
        brackets  = [b for b in BRACKETS if b in pg_bkt or b in atlas_bkt]
        x = range(len(brackets))
        pg_med    = [pg_bkt.get(b,    {}).get("median_ms", 0) for b in brackets]
        atlas_med = [atlas_bkt.get(b, {}).get("median_ms", 0) for b in brackets]
        w = 0.35
        ax.bar([i - w/2 for i in x], pg_med,    w, label="PFASGroups",
               color=_C0, alpha=0.85)
        ax.bar([i + w/2 for i in x], atlas_med, w, label="PFAS-Atlas",
               color=_C1, alpha=0.85)
        ax.set_xticks(list(x))
        ax.set_xticklabels(brackets, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel("Median time (ms)" if ax is axes[0] else "")
        ax.set_title(ds_name)
        ax.legend(fontsize=8)

    fig.suptitle("Median timing by molecular size", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_timing_by_bracket", imgs_dir)


def plot_timing_ratio(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """PFASGroups/Atlas median ratio per size bracket for each dataset."""
    if not MATPLOTLIB:
        return
    apply_style()
    n_ds = len(datasets)
    fig, axes = plt.subplots(1, n_ds, figsize=(7 * n_ds, 4), sharey=True)
    if n_ds == 1:
        axes = [axes]

    for ax, (ds_name, records) in zip(axes, datasets.items()):
        pg_bkt    = timing_by_bracket(records, "pg_time_ms")
        atlas_bkt = timing_by_bracket(records, "atlas_time_ms")
        brackets  = [b for b in BRACKETS if b in pg_bkt and b in atlas_bkt]
        ratios = [
            pg_bkt[b]["median_ms"] / atlas_bkt[b]["median_ms"]
            if atlas_bkt[b].get("median_ms", 0) > 0 else float("nan")
            for b in brackets
        ]
        bar_colors = [_C0 if (r >= 1 and not math.isnan(r)) else _C1 for r in ratios]
        ax.bar(range(len(brackets)), ratios, color=bar_colors, alpha=0.85)
        ax.axhline(1.0, color="gray", ls="--", lw=1)
        ax.set_xticks(range(len(brackets)))
        ax.set_xticklabels(brackets, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel("PFASGroups / PFAS-Atlas timing ratio" if ax is axes[0] else "")
        ax.set_title(ds_name)
        for j, (bar, r) in enumerate(zip(ax.patches, ratios)):
            if not math.isnan(r):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() * 1.02, f"{r:.2f}×",
                        ha="center", va="bottom", fontsize=7)
        pg_patch  = mpatches.Patch(color=_C0, label="PFASGroups slower (>1)")
        at_patch  = mpatches.Patch(color=_C1, label="PFAS-Atlas slower (<1)")
        ax.legend(handles=[pg_patch, at_patch], fontsize=7)

    fig.suptitle("Timing ratio: PFASGroups / PFAS-Atlas by size bracket",
                 fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_timing_ratio", imgs_dir)


def plot_timing_scatter(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """Scatter: timing vs n_heavy_atoms for both tools and datasets."""
    if not MATPLOTLIB or not NUMPY:
        return
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    titles = ["PFASGroups", "PFAS-Atlas"]
    keys   = ["pg_time_ms", "atlas_time_ms"]
    colors_ds = {"OECD": _C2, "clinventory": _C3}

    for ax, title, key in zip(axes, titles, keys):
        for ds_name, records in datasets.items():
            xs = [r["n_heavy_atoms"] for r in records
                  if r.get("n_heavy_atoms") and r.get(key)]
            ys = [r[key] for r in records
                  if r.get("n_heavy_atoms") and r.get(key)]
            if not xs:
                continue
            ax.scatter(xs, ys, s=6, alpha=0.35,
                       color=colors_ds.get(ds_name, "#888"),
                       label=ds_name, marker=DATASET_MARKERS.get(ds_name, "o"))
        ax.set_xlabel("Heavy atom count")
        ax.set_ylabel("Time (ms)")
        ax.set_yscale("log")
        ax.set_title(title)
        handles, _ = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)

    fig.suptitle("Classification time vs. molecular size", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_timing_vs_size", imgs_dir)


# ---------------------------------------------------------------------------
# ── COMPLEXITY PLOTS ─────────────────────────────────────────────────────
# ---------------------------------------------------------------------------
def _plot_hist_compare(ax, datasets: Dict[str, List[dict]], key: str,
                       xlabel: str, bins: int = 40, log: bool = False):
    """Overlapping histograms for a complexity metric across datasets."""
    colors_ds = {"OECD": _C2, "clinventory": _C3}
    for ds_name, records in datasets.items():
        vals = [r[key] for r in records if r.get(key) is not None]
        if not vals:
            continue
        ax.hist(vals, bins=bins, density=True, alpha=0.55,
                color=colors_ds.get(ds_name, "#888"), label=ds_name,
                log=log)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)


def plot_complexity_summary(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """2×2 panel: size, fluorine ratio, branching index, component sizes."""
    if not MATPLOTLIB:
        return
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    _plot_hist_compare(axes[0, 0], datasets, "n_heavy_atoms",
                       "Heavy atom count", bins=50)
    axes[0, 0].set_title("Molecular size")

    _plot_hist_compare(axes[0, 1], datasets, "fluorine_ratio",
                       "Fluorine / heavy atoms", bins=30)
    axes[0, 1].set_title("Fluorine ratio")

    _plot_hist_compare(axes[1, 0], datasets, "branching_index",
                       "Branching index (avg excess C degree)", bins=40)
    axes[1, 0].set_title("Branching")

    _plot_hist_compare(axes[1, 1], datasets, "pg_max_comp_size",
                       "Largest component (atoms)", bins=40)
    axes[1, 1].set_title("Largest PFASGroups component size")

    fig.suptitle("Dataset complexity comparison: OECD vs clinventory",
                 fontweight="bold", y=1.01)
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_complexity_summary", imgs_dir)


def plot_complexity_boxplots(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """Box plots comparing OECD vs clinventory for each complexity metric."""
    if not MATPLOTLIB:
        return
    apply_style()
    metrics = [
        ("n_heavy_atoms",   "Heavy atom count"),
        ("fluorine_ratio",  "Fluorine / heavy atoms"),
        ("branching_index", "Branching index"),
        ("n_rings",         "Number of rings"),
        ("n_rot_bonds",     "Rotatable bonds"),
        ("pg_max_comp_size","Largest component size"),
    ]
    n = len(metrics)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    axes_flat = axes.flatten()
    colors_ds = {"OECD": _C2, "clinventory": _C3}

    for ax, (key, label) in zip(axes_flat, metrics):
        data, tick_labels, bcolors = [], [], []
        for ds_name, records in datasets.items():
            vals = [r[key] for r in records if r.get(key) is not None]
            if vals:
                data.append(vals)
                tick_labels.append(ds_name)
                bcolors.append(colors_ds.get(ds_name, "#888"))
        if not data:
            ax.set_visible(False); continue
        bp = ax.boxplot(data, tick_labels=tick_labels, patch_artist=True,
                        showfliers=False, widths=0.5)
        for patch, c in zip(bp["boxes"], bcolors):
            patch.set_facecolor(c); patch.set_alpha(0.75)
        ax.set_ylabel(label, fontsize=8)
        ax.set_title(label, fontsize=9)

    if n < len(axes_flat):
        for ax in axes_flat[n:]:
            ax.set_visible(False)

    fig.suptitle("Molecular complexity: OECD vs clinventory",
                 fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_complexity_boxplots", imgs_dir)


def plot_fluorine_ratio_vs_timing(datasets: Dict[str, List[dict]],
                                  imgs_dir: Path):
    """Scatter: fluorine ratio vs PFASGroups timing for both datasets."""
    if not MATPLOTLIB or not NUMPY:
        return
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    titles = ["PFASGroups", "PFAS-Atlas"]
    keys   = ["pg_time_ms", "atlas_time_ms"]
    colors_ds = {"OECD": _C2, "clinventory": _C3}

    for ax, title, key in zip(axes, titles, keys):
        for ds_name, records in datasets.items():
            xs = [r["fluorine_ratio"] for r in records
                  if r.get("fluorine_ratio") is not None and r.get(key)]
            ys = [r[key] for r in records
                  if r.get("fluorine_ratio") is not None and r.get(key)]
            if not xs:
                continue
            ax.scatter(xs, ys, s=6, alpha=0.35,
                       color=colors_ds.get(ds_name, "#888"),
                       label=ds_name, marker=DATASET_MARKERS.get(ds_name, "o"))
        ax.set_xlabel("Fluorine ratio")
        ax.set_ylabel("Time (ms)")
        ax.set_yscale("log")
        ax.set_title(title)
        ax.legend(fontsize=8)

    fig.suptitle("Classification time vs. fluorine ratio", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_timing_vs_fluorine", imgs_dir)


def plot_component_size_distribution(datasets: Dict[str, List[dict]],
                                     imgs_dir: Path):
    """Distribution of largest component sizes from PFASGroups per dataset."""
    if not MATPLOTLIB:
        return
    apply_style()
    fig, ax = plt.subplots(figsize=(7, 4))
    colors_ds = {"OECD": _C2, "clinventory": _C3}

    for ds_name, records in datasets.items():
        vals = [r["pg_max_comp_size"] for r in records
                if r.get("pg_max_comp_size", 0) > 0]
        if not vals:
            continue
        ax.hist(vals, bins=40, density=True, alpha=0.6,
                color=colors_ds.get(ds_name, "#888"), label=ds_name)

    ax.set_xlabel("Largest PFASGroups component size (atoms)")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of largest detected component size")
    ax.legend()
    save_fig(fig, "oecd_clinventory_component_sizes", imgs_dir)


# ---------------------------------------------------------------------------
# ── AGREEMENT PLOT ───────────────────────────────────────────────────────
# ---------------------------------------------------------------------------
def plot_agreement(datasets: Dict[str, List[dict]], imgs_dir: Path):
    """Bar chart of PFASGroups / Atlas agreement per dataset."""
    if not MATPLOTLIB:
        return
    apply_style()
    n_ds = len(datasets)
    fig, axes = plt.subplots(1, n_ds, figsize=(6 * n_ds, 4))
    if n_ds == 1:
        axes = [axes]

    keys   = ["both_pfas", "pg_only", "atlas_only", "neither"]
    labels = ["Both PFAS", "PFASGroups\nonly", "PFAS-Atlas\nonly", "Neither"]
    bar_colors = ["#2ca02c", _C0, _C1, "#7f7f7f"]

    for ax, (ds_name, records) in zip(axes, datasets.items()):
        eligible = [r for r in records
                    if r.get("pg_is_pfas") is not None
                    and r.get("atlas_is_pfas") is not None]
        mat = {
            "both_pfas":  sum(1 for r in eligible if r["pg_is_pfas"]  and r["atlas_is_pfas"]),
            "pg_only":    sum(1 for r in eligible if r["pg_is_pfas"]  and not r["atlas_is_pfas"]),
            "atlas_only": sum(1 for r in eligible if not r["pg_is_pfas"] and r["atlas_is_pfas"]),
            "neither":    sum(1 for r in eligible if not r["pg_is_pfas"] and not r["atlas_is_pfas"]),
        }
        total = sum(mat.values())
        values = [mat[k] for k in keys]
        bars = ax.bar(range(4), values, color=bar_colors, alpha=0.85)
        ax.set_xticks(range(4))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("Molecules")
        ax.set_title(ds_name)
        for bar, val in zip(bars, values):
            pct = 100 * val / total if total else 0
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() * 1.01,
                    f"{val:,}\n({pct:.1f}%)",
                    ha="center", va="bottom", fontsize=7)

    fig.suptitle("PFASGroups vs PFAS-Atlas classification agreement",
                 fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "oecd_clinventory_agreement", imgs_dir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Compare PFASGroups vs PFAS-Atlas timing on OECD CSV + clinventory DB"
        )
    )
    p.add_argument("--oecd-csv", type=Path, default=None,
                   help="Path to OECD PFAS CSV (auto-detected in test_data/ if omitted)")
    p.add_argument("--no-clinventory", action="store_true",
                   help="Skip clinventory (useful when DB is unavailable)")
    p.add_argument("--no-atlas", action="store_true",
                   help="Skip PFAS-Atlas (only run PFASGroups)")
    p.add_argument("--limit", type=int, default=None,
                   help="Max molecules per dataset (for quick testing)")
    p.add_argument("--all-molecules", action="store_true",
                   help="Include non-halogenated molecules from clinventory")
    p.add_argument("--db-name",     default=os.environ.get("DATABASE_NAME", "clinventory"))
    p.add_argument("--db-user",     default=os.environ.get("DATABASE_USER", "luc"))
    p.add_argument("--db-password", default=os.environ.get("DATABASE_PASSWORD", None))
    p.add_argument("--db-host",     default=os.environ.get("DATABASE_HOST", "localhost"))
    p.add_argument("--db-port",     type=int,
                   default=int(os.environ.get("DATABASE_PORT", "5432")))
    p.add_argument("--output", type=Path, default=None,
                   help="Override output JSON path")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    imgs_dir = BENCHMARK_DIR / "imgs"
    data_dir = BENCHMARK_DIR / "data"
    data_dir.mkdir(exist_ok=True)

    atlas_enabled = not args.no_atlas and _cp is not None
    if not atlas_enabled:
        if args.no_atlas:
            print("PFAS-Atlas: disabled by --no-atlas flag")
        else:
            print("WARNING: PFAS-Atlas not available — only PFASGroups timings will be produced")

    print("=" * 70)
    print("  OECD × clinventory timing comparison")
    print("=" * 70)

    # ── Locate OECD CSV ──────────────────────────────────────────────────
    oecd_csv = args.oecd_csv
    if oecd_csv is None:
        candidates = sorted(
            (BENCHMARK_DIR / "test_data").glob("S25_OECDPFAS_list_*.csv"),
            key=lambda p: p.stat().st_mtime, reverse=True,
        )
        if candidates:
            oecd_csv = candidates[0]
        else:
            print("ERROR: OECD CSV not found.  Pass --oecd-csv PATH.")
            sys.exit(1)
    print(f"\n  OECD CSV  : {oecd_csv}")

    # ── OECD dataset ─────────────────────────────────────────────────────
    print("\n[1/2]  Loading OECD dataset …")
    oecd_mols = load_oecd_csv(oecd_csv, limit=args.limit)
    print(f"  {len(oecd_mols):,} molecules with SMILES found")

    print("  Classifying …")
    oecd_records = classify_dataset(oecd_mols, "OECD", atlas_enabled)

    # ── clinventory dataset ───────────────────────────────────────────────
    clin_records: List[dict] = []
    if args.no_clinventory:
        print("\n[2/2]  clinventory: skipped (--no-clinventory)")
    elif not PSYCOPG2:
        print("\n[2/2]  clinventory: skipped (psycopg2 not installed)")
    else:
        print(f"\n[2/2]  clinventory DB: {args.db_name}@{args.db_host}:{args.db_port}")
        if args.db_password:
            pw = args.db_password
        else:
            try:
                pw = getpass.getpass(f"  Password for '{args.db_user}': ")
            except (EOFError, KeyboardInterrupt):
                print("\n  No password entered — skipping clinventory.")
                pw = None
        if pw is not None:
            try:
                db_params = dict(
                    dbname=args.db_name, user=args.db_user, password=pw,
                    host=args.db_host,   port=args.db_port,
                )
                clin_mols = fetch_clinventory(db_params, args.limit, args.all_molecules)
                print(f"  {len(clin_mols):,} molecules fetched from DB")
                print("  Classifying …")
                clin_records = classify_dataset(clin_mols, "clinventory", atlas_enabled)
            except Exception as exc:
                print(f"  WARNING: DB unavailable ({exc}) — clinventory skipped")

    # ── Assemble datasets dict ────────────────────────────────────────────
    all_records = oecd_records + clin_records
    datasets: Dict[str, List[dict]] = {"OECD": oecd_records}
    if clin_records:
        datasets["clinventory"] = clin_records

    # ── Save JSON ─────────────────────────────────────────────────────────
    out_path = args.output or (data_dir / f"oecd_clinventory_timing_{ts}.json")

    def _clean(obj):
        if isinstance(obj, float) and math.isnan(obj): return None
        if isinstance(obj, dict):  return {k: _clean(v) for k, v in obj.items()}
        if isinstance(obj, list):  return [_clean(v) for v in obj]
        return obj

    summary = {
        "timestamp": ts,
        "datasets": {
            name: {
                "n_total":   len(recs),
                "n_valid":   sum(1 for r in recs if r.get("pg_time_ms")),
                "pg_timing": _clean(stats([r["pg_time_ms"] for r in recs if r.get("pg_time_ms")])),
                **({"atlas_timing": _clean(stats(
                    [r["atlas_time_ms"] for r in recs if r.get("atlas_time_ms")]
                ))} if atlas_enabled else {}),
            }
            for name, recs in datasets.items()
        },
        "molecules": _clean(all_records),
    }
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\n  Results saved to: {out_path}")

    # ── Plots ────────────────────────────────────────────────────────────
    if not MATPLOTLIB:
        print("\nPlots skipped (matplotlib not available).")
        return

    print("\nGenerating plots …")
    plot_timing_box(datasets, imgs_dir)
    plot_timing_cdf(datasets, imgs_dir)
    plot_timing_by_bracket(datasets, imgs_dir)
    if atlas_enabled:
        plot_timing_ratio(datasets, imgs_dir)
        plot_agreement(datasets, imgs_dir)
        plot_fluorine_ratio_vs_timing(datasets, imgs_dir)
    plot_timing_scatter(datasets, imgs_dir)
    plot_complexity_summary(datasets, imgs_dir)
    plot_complexity_boxplots(datasets, imgs_dir)
    plot_component_size_distribution(datasets, imgs_dir)

    # ── Summary table ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("TIMING SUMMARY  (ms per molecule)")
    print(f"  {'Dataset':<15} {'Tool':<15} {'n':>6} {'Median':>8} {'Mean':>8} {'P95':>8}")
    print("  " + "-" * 60)
    for ds_name, recs in datasets.items():
        for tool, key in [("PFASGroups", "pg_time_ms"), ("PFAS-Atlas", "atlas_time_ms")]:
            times = [r[key] for r in recs if r.get(key)]
            if not times:
                continue
            s = stats(times)
            print(f"  {ds_name:<15} {tool:<15} {s['n']:>6,} "
                  f"{s['median']:>8.3f} {s['mean']:>8.3f} {s['p95']:>8.3f}")
    print()

    print("COMPLEXITY SUMMARY")
    for ds_name, recs in datasets.items():
        valid = [r for r in recs if r.get("n_heavy_atoms")]
        if not valid:
            continue
        ha = [r["n_heavy_atoms"]   for r in valid]
        fr = [r["fluorine_ratio"]  for r in valid]
        bi = [r["branching_index"] for r in valid]
        mcs = [r["pg_max_comp_size"] for r in valid if r.get("pg_max_comp_size")]
        def _med(v): return sorted(v)[len(v) // 2] if v else float("nan")
        print(f"  {ds_name}: n={len(valid):,}  "
              f"size_med={_med(ha):.0f}  "
              f"F_ratio_med={_med(fr):.2f}  "
              f"branching_med={_med(bi):.2f}  "
              f"comp_size_med={_med(mcs):.0f}")
    print("=" * 70)


if __name__ == "__main__":
    main()
