#!/usr/bin/env python3
"""
Classify all molecules in the clinventory database using PFAS-Atlas.

Fetches every molecule from the `molecules` table, runs the PFAS-Atlas
classify_pfas_molecule() on each, records per-molecule timing, and writes
a JSON result file to:

    data/pfasatlas_clinventory_<TIMESTAMP>.json

Must be run inside the `pfasatlas` conda environment:
    mamba activate pfasatlas
    python scripts/classify_pfasatlas_clinventory.py [OPTIONS]

Environment variables:
    DATABASE_NAME      default: clinventory
    DATABASE_USER      default: luc
    DATABASE_PASSWORD  default: (prompted)
    DATABASE_HOST      default: localhost
    DATABASE_PORT      default: 5432
"""

import argparse
import datetime as dt
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Path setup — locate PFAS-atlas relative to this repo
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT = BENCHMARK_DIR.parent
ATLAS_DIR = REPO_ROOT.parent / "PFAS-atlas"

if str(ATLAS_DIR) not in sys.path:
    sys.path.insert(0, str(ATLAS_DIR))

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    print(f"ERROR: rdkit not available: {exc}")
    sys.exit(1)

try:
    from classification_helper.classify_pfas import classify_pfas_molecule
    ATLAS_AVAILABLE = True
except ImportError:
    # Try alternative import path
    try:
        sys.path.insert(0, str(ATLAS_DIR / "classification_helper"))
        from classify_pfas import classify_pfas_molecule
        ATLAS_AVAILABLE = True
    except ImportError as exc:
        print(f"ERROR: PFAS-Atlas not available at {ATLAS_DIR}: {exc}")
        print("  → Make sure you have activated the 'pfasatlas' environment")
        print("  → mamba activate pfasatlas")
        sys.exit(1)

try:
    import psycopg2
except ImportError:
    print("ERROR: psycopg2 not installed.  Run: pip install psycopg2-binary")
    sys.exit(1)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Classify clinventory molecules with PFAS-Atlas"
    )
    p.add_argument("--db-name",     default=os.environ.get("DATABASE_NAME", "clinventory"))
    p.add_argument("--db-user",     default=os.environ.get("DATABASE_USER", "luc"))
    p.add_argument("--db-password", default=os.environ.get("DATABASE_PASSWORD", None),
                   help="DB password (env: DATABASE_PASSWORD; prompted if omitted)")
    p.add_argument("--db-host",     default=os.environ.get("DATABASE_HOST", "localhost"))
    p.add_argument("--db-port",     type=int,
                   default=int(os.environ.get("DATABASE_PORT", "5432")))
    p.add_argument("--limit",       type=int, default=None)
    p.add_argument("--offset",      type=int, default=0)
    p.add_argument("--all-molecules", action="store_true",
                   help="Process ALL molecules (default: only fluorine-containing ones)")
    p.add_argument("--output", type=Path, default=None)
    return p.parse_args()


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def get_password(args: argparse.Namespace) -> str:
    if args.db_password:
        return args.db_password
    import getpass
    return getpass.getpass(f"Password for DB user '{args.db_user}': ")


def connect(args: argparse.Namespace, password: str):
    return psycopg2.connect(
        dbname=args.db_name,
        user=args.db_user,
        password=password,
        host=args.db_host,
        port=args.db_port,
    )


def fetch_molecules(conn, all_molecules: bool, limit: Optional[int],
                    offset: int) -> List[Tuple[int, str, Optional[str]]]:
    if all_molecules:
        where = "WHERE smiles IS NOT NULL AND smiles != ''"
    else:
        where = (
            "WHERE smiles IS NOT NULL AND smiles != '' "
            "AND formula IS NOT NULL "
            "AND formula ~ '(F[0-9]*)'"
        )
    limit_clause = f"LIMIT {limit}" if limit else ""
    offset_clause = f"OFFSET {offset}" if offset else ""
    query = f"SELECT id, smiles, formula FROM molecules {where} ORDER BY id {limit_clause} {offset_clause}"
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

# PFAS-Atlas top-level classes considered as "PFAS"
_PFAS_CLASSES = {
    "PFAAs", "PFAA precursors", "Fluorotelomers", "Other PFASs",
    "PolyFAAs", "Complex structure", "Silicon PFASs", "Side-chain aromatics",
}

def is_pfas(class1: str) -> bool:
    """Return True if PFAS-Atlas classified the molecule as any kind of PFAS."""
    return class1 != "Not PFAS"


def atom_bracket(n: Optional[int]) -> str:
    if n is None:
        return "unknown"
    if n < 10:
        return "tiny (<10)"
    if n < 20:
        return "small (10-19)"
    if n < 35:
        return "medium (20-34)"
    if n < 60:
        return "large (35-59)"
    return "very large (>=60)"


def stats(vals: list) -> dict:
    if not vals:
        return {}
    n = len(vals)
    vals_s = sorted(vals)
    mu = sum(vals) / n
    var = sum((v - mu) ** 2 for v in vals) / n
    return {
        "n": n,
        "mean_ms": mu,
        "std_ms": math.sqrt(var),
        "median_ms": vals_s[n // 2],
        "min_ms": vals_s[0],
        "max_ms": vals_s[-1],
        "p95_ms": vals_s[int(0.95 * n)],
        "total_ms": sum(vals),
    }


def classify_molecule(mol_id: int, smiles: str, formula: Optional[str]) -> Dict:
    """Run PFAS-Atlas classify_pfas_molecule and return a result dict."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "id": mol_id, "smiles": smiles, "formula": formula,
            "time_ms": None, "n_heavy_atoms": None,
            "n_fluorine": None, "n_halogens": None,
            "class1": None, "class2": None,
            "is_pfas": False, "error": "invalid_smiles",
        }

    n_heavy = mol.GetNumHeavyAtoms()
    n_f     = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
    n_hal   = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9, 17, 35, 53))

    # PFAS-Atlas expects canonical RDKit SMILES
    canon_smiles = Chem.MolToSmiles(mol)

    t0 = time.perf_counter()
    try:
        result = classify_pfas_molecule(canon_smiles)
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        c1 = result[0] if result else "error"
        c2 = result[1] if (result and len(result) > 1) else ""
        return {
            "id": mol_id, "smiles": smiles, "formula": formula,
            "time_ms": elapsed_ms, "n_heavy_atoms": n_heavy,
            "n_fluorine": n_f, "n_halogens": n_hal,
            "class1": c1, "class2": c2,
            "is_pfas": is_pfas(c1), "error": None,
        }
    except Exception as exc:
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        return {
            "id": mol_id, "smiles": smiles, "formula": formula,
            "time_ms": elapsed_ms, "n_heavy_atoms": n_heavy,
            "n_fluorine": n_f, "n_halogens": n_hal,
            "class1": "error", "class2": str(exc),
            "is_pfas": False, "error": str(exc),
        }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    password = get_password(args)

    print("=" * 70)
    print("PFAS-Atlas — clinventory classification benchmark")
    print("=" * 70)
    print(f"  DB   : {args.db_user}@{args.db_host}:{args.db_port}/{args.db_name}")
    print(f"  Atlas: {ATLAS_DIR}")

    conn = connect(args, password)
    print("  DB connection: OK")

    print("\nFetching molecules …")
    t0 = time.perf_counter()
    rows = fetch_molecules(conn, args.all_molecules, args.limit, args.offset)
    conn.close()
    print(f"  {len(rows):,} molecules fetched in {time.perf_counter()-t0:.1f}s")

    print("\nClassifying …")
    timestamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    results = []
    errors = 0
    n_total = len(rows)
    t_start = time.perf_counter()

    for i, (mol_id, smiles, formula) in enumerate(rows):
        rec = classify_molecule(mol_id, smiles, formula)
        results.append(rec)
        if rec["error"]:
            errors += 1
        if (i + 1) % 500 == 0 or (i + 1) == n_total:
            elapsed = time.perf_counter() - t_start
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            print(f"  {i+1:>7,}/{n_total:,}  ({rate:,.0f} mol/s)  errors: {errors}", end="\r")

    t_total = time.perf_counter() - t_start
    print(f"\n  Done in {t_total:.1f}s  ({n_total/t_total:.0f} mol/s)  errors: {errors}")

    # Summary
    valid       = [r for r in results if r["time_ms"] is not None]
    pfas_mols   = [r for r in valid if r["is_pfas"]]
    times       = [r["time_ms"] for r in valid]

    # Class distribution
    class_dist: Dict[str, int] = {}
    for r in valid:
        k = r["class1"] or "None"
        class_dist[k] = class_dist.get(k, 0) + 1
    class2_dist: Dict[str, int] = {}
    for r in valid:
        k = r["class2"] or "None"
        class2_dist[k] = class2_dist.get(k, 0) + 1

    # Timing by size bracket
    bracket_times: Dict[str, list] = {}
    for r in valid:
        b = atom_bracket(r["n_heavy_atoms"])
        bracket_times.setdefault(b, []).append(r["time_ms"])
    timing_by_bracket = {b: stats(ts) for b, ts in bracket_times.items()}

    # Timing by halogen type
    halogen_cls_times: Dict[str, list] = {}
    for r in valid:
        nf = r["n_fluorine"] or 0
        nh = (r["n_halogens"] or 0) - nf
        if nf > 0 and nh == 0:   cls = "F only"
        elif nf > 0 and nh > 0:  cls = "Mixed (F + other)"
        elif nf == 0 and nh > 0: cls = "Other halogen (no F)"
        else:                     cls = "No halogen"
        halogen_cls_times.setdefault(cls, []).append(r["time_ms"])
    timing_by_halogen_class = {c: stats(ts) for c, ts in halogen_cls_times.items()}

    # Timing by PFAS-Atlas class1
    class_timing: Dict[str, list] = {}
    for r in valid:
        k = r["class1"] or "None"
        class_timing.setdefault(k, []).append(r["time_ms"])
    timing_by_class = {k: stats(ts) for k, ts in class_timing.items()}

    output = {
        "tool": "PFAS-Atlas",
        "timestamp": timestamp,
        "db": {"name": args.db_name, "host": args.db_host, "port": args.db_port, "user": args.db_user},
        "n_molecules_fetched": n_total,
        "n_valid": len(valid),
        "n_errors": errors,
        "n_pfas": len(pfas_mols),
        "pfas_rate_pct": 100.0 * len(pfas_mols) / len(valid) if valid else 0,
        "total_classification_time_s": t_total,
        "class1_distribution": class_dist,
        "class2_distribution": class2_dist,
        "timing_overall": stats(times),
        "timing_by_atom_bracket": timing_by_bracket,
        "timing_by_halogen_class": timing_by_halogen_class,
        "timing_by_pfas_class": timing_by_class,
        "molecules": results,
    }

    data_dir = BENCHMARK_DIR / "data"
    data_dir.mkdir(exist_ok=True)
    out_path = args.output or (data_dir / f"pfasatlas_clinventory_{timestamp}.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to: {out_path}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Molecules processed : {n_total:,}")
    print(f"  Valid (parsed)      : {len(valid):,}")
    print(f"  Classified as PFAS  : {len(pfas_mols):,}  ({100*len(pfas_mols)/max(len(valid),1):.1f}%)")
    if times:
        s = stats(times)
        print(f"  Timing median       : {s['median_ms']:.3f} ms/mol")
        print(f"  Timing mean ± std   : {s['mean_ms']:.3f} ± {s['std_ms']:.3f} ms/mol")
    print("\n  Class distribution (top-10):")
    for cls, cnt in sorted(class_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"    {cls:<40s}: {cnt:>7,}")
    print(f"\nOutput: {out_path}")


if __name__ == "__main__":
    main()
