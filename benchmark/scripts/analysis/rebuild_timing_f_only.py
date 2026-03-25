#!/usr/bin/env python3
"""
Rebuild the timing comparison JSONs restricted to PFASGroups(halogens='F').

Produces:
  data/oecd_clinventory_timing_f_only_<TS>.json
      – summary stats for OECD, CLInventory (F) and stress datasets
      – used by generate_si_figures.py for timing_box / timing_by_bracket

  data/clinventory_comparison_f_only_<TS>.json
      – timing_overall, timing_by_atom_bracket, agreement_matrix
      – used by generate_si_figures.py for atlas_agreement_bar / atlas_disagreement

Sources:
  • OECD molecules         : existing oecd_clinventory_timing_20260315_062059.json
                             (SMILES + Atlas timing reused; PFASGroups re-run F-only)
  • CLInventory F-only PFG : pfasgroups_clinventory_20260325T105742.json  (today's run)
  • CLInventory Atlas      : pfasatlas_clinventory_20260313T231144.json
  • Stress molecules       : pfas_timing_benchmark_full_20260313_201139.json
                             (SMILES extracted; PFASGroups re-run F-only, 5 iters)

Usage (from benchmark/):
    python scripts/analysis/rebuild_timing_f_only.py [--stress-iters N]
"""

import argparse
import json
import math
import os
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT     = BENCHMARK_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    print(f"ERROR: rdkit not available: {exc}")
    sys.exit(1)

try:
    from PFASGroups import parse_mol
    from PFASGroups.core import rdkit_disable_log
    rdkit_disable_log()
except ImportError as exc:
    print(f"ERROR: PFASGroups not available: {exc}")
    sys.exit(1)

try:
    import numpy as np
    NUMPY = True
except ImportError:
    NUMPY = False
    print("WARNING: numpy not available — using stdlib stats only")


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------
def _stats_np(vals):
    """Return rich stats dict using numpy."""
    arr = np.array(vals)
    return {
        "n":      int(len(arr)),
        "mean":   float(arr.mean()),
        "std":    float(arr.std()),
        "median": float(np.median(arr)),
        "min":    float(arr.min()),
        "max":    float(arr.max()),
        "p25":    float(np.percentile(arr, 25)),
        "p75":    float(np.percentile(arr, 75)),
        "p95":    float(np.percentile(arr, 95)),
    }


def _stats_stdlib(vals):
    """Pure-stdlib stats fallback."""
    if not vals:
        return {}
    s = sorted(vals)
    n = len(s)
    mu = sum(vals) / n
    var = sum((v - mu) ** 2 for v in vals) / n
    return {
        "n":      n,
        "mean":   mu,
        "std":    math.sqrt(var),
        "median": s[n // 2],
        "min":    s[0],
        "max":    s[-1],
        "p25":    s[n // 4]        if n >= 4 else s[0],
        "p75":    s[3 * n // 4]    if n >= 4 else s[-1],
        "p95":    s[int(0.95 * n)] if n >= 2 else s[-1],
    }


def make_stats(vals):
    if not vals:
        return {"n": 0}
    return _stats_np(vals) if NUMPY else _stats_stdlib(vals)


# ---------------------------------------------------------------------------
# PFASGroups F-only runner
# ---------------------------------------------------------------------------
def run_pfasgroups_f(mol) -> tuple:
    """Return (is_pfas, elapsed_ms, error)."""
    try:
        t0 = time.perf_counter()
        result = parse_mol(mol, halogens="F", include_PFAS_definitions=False)
        elapsed = (time.perf_counter() - t0) * 1000.0
        n_groups = sum(
            1 for m in result.get("matches", []) if m.get("type") == "HalogenGroup"
        ) if isinstance(result, dict) else 0
        return n_groups > 0, elapsed, None
    except Exception as exc:
        return False, 0.0, str(exc)


def atom_bracket(n):
    if n is None:  return "unknown"
    if n < 10:     return "tiny (<10)"
    if n < 20:     return "small (10-19)"
    if n < 35:     return "medium (20-34)"
    if n < 60:     return "large (35-59)"
    return "very large (>=60)"


# ---------------------------------------------------------------------------
# Data file paths (latest versions)
# ---------------------------------------------------------------------------
DATA_DIR = BENCHMARK_DIR / "data"


def _latest(glob_pat: str) -> Path:
    """Return newest file matching glob under DATA_DIR."""
    candidates = sorted(DATA_DIR.glob(glob_pat),
                        key=lambda p: p.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No file matching {DATA_DIR / glob_pat}")
    return candidates[0]


# ---------------------------------------------------------------------------
# 1.  Re-run PFASGroups(F) on OECD molecules
# ---------------------------------------------------------------------------
@rdkit_disable_log()
def rebuild_oecd(oecd_timing_path: Path) -> dict:
    print(f"\n[1/3] OECD re-run (halogens='F')")
    print(f"      source: {oecd_timing_path.name}")

    with open(oecd_timing_path, encoding="utf-8") as fh:
        src = json.load(fh)

    oecd_mols = [m for m in src.get("molecules", [])
                 if m.get("dataset") == "OECD"]
    print(f"      {len(oecd_mols):,} OECD molecules found")

    pg_times, at_times = [], []
    n_err = 0
    total = len(oecd_mols)
    report_every = max(1, total // 20)

    for i, rec in enumerate(oecd_mols, 1):
        if i % report_every == 0 or i == total:
            print(f"      {i}/{total} …", end="\r")
        mol = Chem.MolFromSmiles(rec.get("smiles", ""))
        if mol is None:
            n_err += 1
            continue
        _, pg_ms, err = run_pfasgroups_f(mol)
        if err:
            n_err += 1
        else:
            pg_times.append(pg_ms)
        at_ms = rec.get("atlas_time_ms")
        if at_ms:
            at_times.append(at_ms)

    print(f"\n      errors: {n_err}, valid pg: {len(pg_times)}, valid atlas: {len(at_times)}")
    return {
        "pg_timing":    make_stats(pg_times),
        "atlas_timing": make_stats(at_times),
        "n_total":      len(oecd_mols),
        "n_valid":      len(pg_times),
    }


# ---------------------------------------------------------------------------
# 2.  Build CLInventory(F) summary from existing JSONs
# ---------------------------------------------------------------------------
def build_clinventory_f(pfg_path: Path, atlas_path: Path) -> dict:
    print(f"\n[2/3] ClinInventory (F) from pre-computed data")
    print(f"      pfg  : {pfg_path.name}")
    print(f"      atlas: {atlas_path.name}")

    with open(pfg_path,   encoding="utf-8") as fh:
        pfg_data = json.load(fh)
    with open(atlas_path, encoding="utf-8") as fh:
        at_data  = json.load(fh)

    pfg_mols  = {r["id"]: r for r in pfg_data.get("molecules", [])}
    at_mols   = {r["id"]: r for r in at_data.get("molecules",  [])}

    pg_times, at_times = [], []
    for mid, rec in pfg_mols.items():
        t = rec.get("time_ms")
        if t is not None:
            pg_times.append(t)
    for mid, rec in at_mols.items():
        t = rec.get("time_ms")
        if t is not None:
            at_times.append(t)

    n = len(pfg_mols)
    print(f"      n_pfg={n:,}, valid pg times={len(pg_times):,}, valid atlas times={len(at_times):,}")
    return {
        "pg_timing":    make_stats(pg_times),
        "atlas_timing": make_stats(at_times),
        "n_total":      n,
        "n_valid":      len(pg_times),
    }


# ---------------------------------------------------------------------------
# 3.  Re-run stress molecules (halogens='F'), N iterations
# ---------------------------------------------------------------------------
@rdkit_disable_log()
def rebuild_stress(stress_path: Path, n_iters: int = 5) -> tuple:
    """Returns (pg_times_ms, at_times_ms) for ≥35-atom molecules."""
    print(f"\n[3/3] Stress dataset re-run (halogens='F', {n_iters} iters)")
    print(f"      source: {stress_path.name}")

    with open(stress_path, encoding="utf-8") as fh:
        recs = json.load(fh)

    big_mols = [r for r in recs if isinstance(r, dict) and r.get("num_atoms", 0) >= 35]
    print(f"      {len(big_mols):,} molecules with ≥35 heavy atoms")

    pg_times, at_times = [], []
    n_err = 0
    total = len(big_mols)
    report_every = max(1, total // 20)

    for i, rec in enumerate(big_mols, 1):
        if i % report_every == 0 or i == total:
            print(f"      {i}/{total} …", end="\r")
        smiles = rec.get("smiles", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            n_err += 1
            continue
        # F-only timing: multiple iterations
        iter_times = []
        for _ in range(n_iters):
            _, ms, err = run_pfasgroups_f(mol)
            if err:
                n_err += 1
                break
            iter_times.append(ms)
        if iter_times:
            pg_times.append(sum(iter_times) / len(iter_times))
        # Keep original atlas timing
        at_avg = rec.get("atlas_time_avg")
        if at_avg is not None:
            at_times.append(at_avg * 1000.0)

    print(f"\n      errors: {n_err}, valid pg: {len(pg_times)}, valid atlas: {len(at_times)}")
    return pg_times, at_times


# ---------------------------------------------------------------------------
# 4.  Build clinventory_comparison JSON (for agreement_matrix + timing CDF)
# ---------------------------------------------------------------------------
def build_clinventory_comparison(pfg_path: Path, atlas_path: Path) -> dict:
    """Build the clinventory_comparison structure used by generate_si_figures.

    Includes timing_overall, timing_by_atom_bracket, agreement_matrix,
    and group_name_counts (count of each PFASGroups group name detected).
    """
    from collections import Counter
    print(f"\n[4/4] Building clinventory_comparison (F-only)")

    with open(pfg_path,   encoding="utf-8") as fh:
        pfg_data = json.load(fh)
    with open(atlas_path, encoding="utf-8") as fh:
        at_data  = json.load(fh)

    pfg_mols   = {r["id"]: r for r in pfg_data.get("molecules", [])}
    at_mols    = {r["id"]: r for r in at_data.get("molecules",  [])}
    common_ids = sorted(set(pfg_mols) & set(at_mols))

    # ── timing ────────────────────────────────────────────────────────────
    pg_times: list  = []
    at_times: list  = []
    bracket_pg: dict  = defaultdict(list)
    bracket_at: dict  = defaultdict(list)

    for mid in common_ids:
        pr  = pfg_mols[mid]
        ar  = at_mols[mid]
        pg_t = pr.get("time_ms")
        at_t = ar.get("time_ms")
        n_ha = pr.get("n_heavy_atoms") or ar.get("n_heavy_atoms")
        b    = atom_bracket(n_ha)
        if pg_t is not None:
            pg_times.append(pg_t)
            bracket_pg[b].append(pg_t)
        if at_t is not None:
            at_times.append(at_t)
            bracket_at[b].append(at_t)

    timing_overall = {
        "hg":    {**make_stats(pg_times)},
        "atlas": {**make_stats(at_times)},
        "ratio_hg_over_atlas": (
            (make_stats(pg_times).get("mean", 0) /
             make_stats(at_times).get("mean", 1))
            if at_times else None
        ),
    }

    # Bracket ordering expected by generate_si_figures
    bracket_order = [
        "tiny (<10)", "small (10-19)", "medium (20-34)",
        "large (35-59)", "very large (>=60)",
    ]
    timing_by_bracket = {}
    for b in bracket_order:
        pg_b = bracket_pg.get(b, [])
        at_b = bracket_at.get(b, [])
        pg_med = make_stats(pg_b).get("median", 0) if pg_b else 0
        at_med = make_stats(at_b).get("median", 0) if at_b else 0
        timing_by_bracket[b] = {
            "n":                   len(pg_b),
            "hg_median":           pg_med,
            "atlas_median":        at_med,
            "ratio_hg_over_atlas": pg_med / at_med if at_med > 0 else None,
        }

    # ── agreement matrix ─────────────────────────────────────────────────
    both_pfas  = 0
    only_hg    = 0
    only_atlas = 0
    neither    = 0

    for mid in common_ids:
        pr = pfg_mols[mid]
        ar = at_mols[mid]
        pfg_pfas   = pr.get("has_fluorine_group", False)
        atlas_pfas = ar.get("is_pfas", False)
        if pfg_pfas and atlas_pfas:
            both_pfas += 1
        elif pfg_pfas and not atlas_pfas:
            only_hg += 1
        elif not pfg_pfas and atlas_pfas:
            only_atlas += 1
        else:
            neither += 1

    total       = len(common_ids)
    po          = (both_pfas + neither) / total if total else 0.0
    pfg_pos     = (both_pfas + only_hg)    / total if total else 0.0
    at_pos      = (both_pfas + only_atlas) / total if total else 0.0
    pfg_neg     = 1.0 - pfg_pos
    at_neg      = 1.0 - at_pos
    pe          = pfg_pos * at_pos + pfg_neg * at_neg
    kappa       = (po - pe) / (1.0 - pe) if pe < 1.0 else 1.0
    agree_pct   = 100.0 * po

    print(f"      n={total:,}  both={both_pfas:,}  pfg_only={only_hg:,}  "
          f"atlas_only={only_atlas:,}  neither={neither:,}")
    print(f"      agreement={agree_pct:.2f}%  kappa={kappa:.4f}")

    agreement_matrix = {
        "both_pfas":    both_pfas,
        "only_hg":      only_hg,
        "only_atlas":   only_atlas,
        "neither":      neither,
        "total":        total,
        "agreement_pct": agree_pct,
        "cohen_kappa":  kappa,
    }

    return {
        "timing_overall":        timing_overall,
        "timing_by_atom_bracket": timing_by_bracket,
        "agreement_matrix":      agreement_matrix,
        "common_molecules":      total,
        "hg_only_molecules":     only_hg,
        "atlas_only_molecules":  only_atlas,
        "group_name_counts":     _build_group_name_counts(pfg_path),
    }


def _build_group_name_counts(pfg_path: Path) -> dict:
    """Count occurrences of each group name across all classified molecules."""
    try:
        from collections import Counter
        with open(pfg_path, encoding="utf-8") as fh:
            pfg_data = json.load(fh)
        gnc: Counter = Counter()
        for r in pfg_data.get("molecules", []):
            for gname in r.get("group_names", []):
                if gname:
                    gnc[gname] += 1
        return dict(gnc)
    except Exception as exc:
        print(f"  WARNING: could not build group_name_counts: {exc}")
        return {}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Rebuild timing JSONs with PFASGroups halogens='F'")
    p.add_argument("--stress-iters", type=int, default=5,
                   help="Number of timing iterations for stress molecules (default: 5)")
    p.add_argument("--no-stress", action="store_true",
                   help="Skip stress dataset re-run (use existing Atlas times only)")
    return p.parse_args()


def _clean_json(obj):
    """Replace float NaN/Inf with None for JSON serialisation."""
    if isinstance(obj, float):
        return None if (math.isnan(obj) or math.isinf(obj)) else obj
    if isinstance(obj, dict):
        return {k: _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_json(v) for v in obj]
    return obj


def main():
    args = parse_args()
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")

    # ── Locate input files ────────────────────────────────────────────────
    oecd_timing_path = _latest("oecd_clinventory_timing_20260315_062059.json")
    pfg_clin_path    = _latest("pfasgroups_clinventory_20260325T*.json")
    at_clin_path     = _latest("pfasatlas_clinventory_*.json")
    stress_path      = _latest("pfas_timing_benchmark_full_*.json")

    print("=" * 70)
    print("  Rebuild timing JSONs — PFASGroups(halogens='F')")
    print(f"  Timestamp: {ts}")
    print("=" * 70)
    print(f"  OECD timing source  : {oecd_timing_path.name}")
    print(f"  PFG clinventory     : {pfg_clin_path.name}")
    print(f"  Atlas clinventory   : {at_clin_path.name}")
    print(f"  Stress benchmark    : {stress_path.name}")

    # ── Process datasets ─────────────────────────────────────────────────
    oecd_stats   = rebuild_oecd(oecd_timing_path)
    clin_stats   = build_clinventory_f(pfg_clin_path, at_clin_path)

    if args.no_stress:
        print("\n[3/3] Stress dataset: skipped (--no-stress)")
        stress_pg_times, stress_at_times = [], []
        # Fall back to existing Atlas times
        with open(stress_path, encoding="utf-8") as fh:
            stress_raw = json.load(fh)
        stress_at_times = [r["atlas_time_avg"] * 1000
                           for r in stress_raw
                           if isinstance(r, dict) and r.get("num_atoms", 0) >= 35]
    else:
        stress_pg_times, stress_at_times = rebuild_stress(stress_path, args.stress_iters)

    # ── Assemble oecd_clinventory_timing JSON ─────────────────────────────
    stress_stats_pg  = make_stats(stress_pg_times)  if stress_pg_times  else {}
    stress_stats_at  = make_stats(stress_at_times)  if stress_at_times  else {}

    oecd_timing_out = {
        "timestamp": ts,
        "note": "PFASGroups run with halogens='F' on all datasets",
        "datasets": {
            "OECD": {
                "n_total":      oecd_stats["n_total"],
                "n_valid":      oecd_stats["n_valid"],
                "pg_timing":    oecd_stats["pg_timing"],
                "atlas_timing": oecd_stats["atlas_timing"],
            },
            "clinventory (F)": {
                "n_total":      clin_stats["n_total"],
                "n_valid":      clin_stats["n_valid"],
                "pg_timing":    clin_stats["pg_timing"],
                "atlas_timing": clin_stats["atlas_timing"],
            },
            "stress (≥35 atoms)": {
                "n_total":      len(stress_pg_times) if stress_pg_times else len(stress_at_times),
                "n_valid":      len(stress_pg_times),
                "pg_timing":    stress_stats_pg,
                "atlas_timing": stress_stats_at,
            },
        },
    }

    out1 = DATA_DIR / f"oecd_clinventory_timing_f_only_{ts}.json"
    with open(out1, "w", encoding="utf-8") as fh:
        json.dump(_clean_json(oecd_timing_out), fh, indent=2)
    print(f"\n  Saved: {out1.name}")

    # ── Build clinventory_comparison JSON ─────────────────────────────────
    comparison = build_clinventory_comparison(pfg_clin_path, at_clin_path)
    comparison["timestamp"] = ts
    comparison["note"]      = "Rebuilt with PFASGroups(halogens='F') from pfasgroups_clinventory_20260325T105742.json"

    out2 = DATA_DIR / f"clinventory_comparison_f_only_{ts}.json"
    with open(out2, "w", encoding="utf-8") as fh:
        json.dump(_clean_json(comparison), fh, indent=2)
    print(f"  Saved: {out2.name}")

    # ── Print summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("TIMING SUMMARY  (ms per molecule, PFASGroups halogens='F')")
    print(f"  {'Dataset':<22} {'n':>6} {'PG median':>10} {'Atlas median':>13}")
    print("  " + "-" * 55)
    for ds, pg, at in [
        ("OECD", oecd_stats["pg_timing"], oecd_stats["atlas_timing"]),
        ("CLInventory (F)", clin_stats["pg_timing"], clin_stats["atlas_timing"]),
        ("Stress (≥35 atoms)", stress_stats_pg, stress_stats_at),
    ]:
        n   = pg.get("n", 0)
        pgm = pg.get("median", float("nan"))
        atm = at.get("median", float("nan"))
        print(f"  {ds:<22} {n:>6,} {pgm:>10.3f} {atm:>13.3f}")
    print()
    am = comparison["agreement_matrix"]
    print("AGREEMENT (CLinventory F-only PFASGroups vs PFAS-Atlas)")
    print(f"  both PFAS        : {am['both_pfas']:>6,}  ({100*am['both_pfas']/am['total']:.1f}%)")
    print(f"  PFASGroups only  : {am['only_hg']:>6,}  ({100*am['only_hg']/am['total']:.1f}%)")
    print(f"  PFAS-Atlas only  : {am['only_atlas']:>6,}  ({100*am['only_atlas']/am['total']:.1f}%)")
    print(f"  neither          : {am['neither']:>6,}  ({100*am['neither']/am['total']:.1f}%)")
    print(f"  Cohen's κ        : {am['cohen_kappa']:.4f}")
    print()
    print(f"Output files (update generate_si_figures.py to these names):")
    print(f"  oecd_timing  → {out1.name}")
    print(f"  clinv_comp   → {out2.name}")


if __name__ == "__main__":
    main()
