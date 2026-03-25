#!/usr/bin/env python3
"""
Classify all molecules in the clinventory database using PFASGroups.

Fetches every molecule that has at least one halogen (F, Cl, Br, I) from the
`molecules` table, runs PFASGroups parse_mol on each, records per-molecule
timing, and writes a JSON result file to:

    data/PFASGroups_clinventory_<TIMESTAMP>.json

Usage (from the benchmark/ directory):
    python scripts/classify_PFASGroups_clinventory.py [OPTIONS]

Environment variables (can also be passed as CLI flags):
    DATABASE_NAME      default: clinventory
    DATABASE_USER      default: luc
    DATABASE_PASSWORD  default: (prompted)
    DATABASE_HOST      default: localhost
    DATABASE_PORT      default: 5432
"""

import argparse
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple, Dict

# ---------------------------------------------------------------------------
# Path setup — prefer local PFASGroups repo
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT = BENCHMARK_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")
except ImportError as exc:
    print(f"ERROR: rdkit not available: {exc}")
    sys.exit(1)

try:
    from PFASGroups import parse_mol, get_HalogenGroups
    from PFASGroups.core import rdkit_disable_log
    rdkit_disable_log()
    PFASGroups_AVAILABLE = True
except ImportError as exc:
    print(f"ERROR: PFASGroups not available: {exc}")
    sys.exit(1)

try:
    import psycopg2
    from psycopg2.extras import RealDictCursor
except ImportError:
    print("ERROR: psycopg2 not installed.  Run: pip install psycopg2-binary")
    sys.exit(1)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Classify clinventory molecules with PFASGroups"
    )
    p.add_argument("--db-name",     default=os.environ.get("DATABASE_NAME", "clinventory"))
    p.add_argument("--db-user",     default=os.environ.get("DATABASE_USER", "luc"))
    p.add_argument("--db-password", default=os.environ.get("DATABASE_PASSWORD", None),
                   help="DB password (env: DATABASE_PASSWORD; prompted if omitted)")
    p.add_argument("--db-host",     default=os.environ.get("DATABASE_HOST", "localhost"))
    p.add_argument("--db-port",     type=int,
                   default=int(os.environ.get("DATABASE_PORT", "5432")))
    p.add_argument("--limit",       type=int, default=None,
                   help="Maximum number of molecules to process")
    p.add_argument("--offset",      type=int, default=0)
    p.add_argument("--batch-size",  type=int, default=500)
    p.add_argument("--all-molecules", action="store_true",
                   help="Process ALL molecules (not just halogenated ones)")
    p.add_argument("--output", type=Path, default=None,
                   help="Override output JSON path")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Database helpers
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


def fetch_molecules(conn, all_molecules: bool, limit: Optional[int], offset: int,
                    batch_size: int) -> List[Tuple[int, str, Optional[str]]]:
    """Return list of (id, smiles, formula) tuples."""
    if all_molecules:
        where = "WHERE smiles IS NOT NULL AND smiles != ''"
    else:
        where = (
            "WHERE smiles IS NOT NULL AND smiles != '' "
            "AND formula IS NOT NULL "
            "AND formula ~ '(F[0-9]*|Cl[0-9]*|Br[0-9]*|I[0-9]*)'"
        )
    limit_clause = f"LIMIT {limit}" if limit else ""
    offset_clause = f"OFFSET {offset}" if offset else ""
    query = f"SELECT id, smiles, formula FROM molecules {where} ORDER BY id {limit_clause} {offset_clause}"

    with conn.cursor() as cur:
        cur.execute(query)
        rows = cur.fetchall()
    return rows  # list of (id, smiles, formula)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify_molecule(mol_id: int, smiles: str, formula: Optional[str]) -> Dict:
    """Run PFASGroups parse_mol and return a result dict."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "id": mol_id,
            "smiles": smiles,
            "formula": formula,
            "time_ms": None,
            "n_heavy_atoms": None,
            "n_fluorine": None,
            "n_halogens": None,
            "classified": False,
            "has_fluorine_group": False,
            "n_groups": 0,
            "group_ids": [],
            "group_names": [],
            "group_categories": [],
            "error": "invalid_smiles",
        }

    # Count atoms
    n_heavy = mol.GetNumHeavyAtoms()
    n_f = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
    n_halogen = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9, 17, 35, 53))

    t0 = time.perf_counter()
    try:
        result = parse_mol(mol, include_PFAS_definitions=False)
        elapsed_ms = (time.perf_counter() - t0) * 1000.0

        if result is None:
            return {
                "id": mol_id,
                "smiles": smiles,
                "formula": formula,
                "time_ms": elapsed_ms,
                "n_heavy_atoms": n_heavy,
                "n_fluorine": n_f,
                "n_halogens": n_halogen,
                "classified": False,
                "has_fluorine_group": False,
                "n_groups": 0,
                "group_ids": [],
                "group_names": [],
                "group_categories": [],
                "error": "no_result",
            }

        # MoleculeResult is a dict-like with a 'matches' list.
        # Each match dict has: type ('HalogenGroup'|'PFASDefinition'),
        #   id, group_name, components, category, ...
        all_matches = result.get("matches", []) if isinstance(result, dict) else []
        groups = [m for m in all_matches if m.get("type") == "HalogenGroup"]

        # Unique groups by ID (a molecule can have multiple components of the same group)
        seen_ids = set()
        unique_groups = []
        for g in groups:
            gid = g.get("id")
            if gid not in seen_ids:
                seen_ids.add(gid)
                unique_groups.append(g)

        group_ids   = [g.get("id")         for g in unique_groups]
        group_names = [g.get("group_name") for g in unique_groups]
        group_cats  = [g.get("category")   for g in unique_groups]

        # "has_fluorine_group": at least one group detected AND molecule contains F
        classified  = len(unique_groups) > 0
        has_f_group = classified and n_f > 0

        return {
            "id": mol_id,
            "smiles": smiles,
            "formula": formula,
            "time_ms": elapsed_ms,
            "n_heavy_atoms": n_heavy,
            "n_fluorine": n_f,
            "n_halogens": n_halogen,
            "classified": classified,
            "has_fluorine_group": has_f_group,
            "n_groups": len(unique_groups),
            "group_ids": group_ids,
            "group_names": group_names,
            "group_categories": group_cats,
            "error": None,
        }

    except Exception as exc:
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        return {
            "id": mol_id,
            "smiles": smiles,
            "formula": formula,
            "time_ms": elapsed_ms,
            "n_heavy_atoms": n_heavy,
            "n_fluorine": n_f,
            "n_halogens": n_halogen,
            "classified": False,
            "has_fluorine_group": False,
            "n_groups": 0,
            "group_ids": [],
            "group_names": [],
            "group_categories": [],
            "error": str(exc),
        }


# ---------------------------------------------------------------------------
# Size bracket helper
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    password = get_password(args)

    print("=" * 70)
    print("PFASGroups — clinventory classification benchmark")
    print("=" * 70)
    print(f"  DB   : {args.db_user}@{args.db_host}:{args.db_port}/{args.db_name}")

    # Connect
    conn = connect(args, password)
    print("  DB connection: OK")

    # Fetch
    print("\nFetching molecules …")
    t_fetch0 = time.perf_counter()
    rows = fetch_molecules(conn, args.all_molecules, args.limit, args.offset, args.batch_size)
    conn.close()
    t_fetch = time.perf_counter() - t_fetch0
    print(f"  {len(rows):,} molecules fetched in {t_fetch:.1f}s")

    # Classify
    print("\nClassifying …")
    timestamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    results = []
    errors = 0
    n_total = len(rows)
    t_class_start = time.perf_counter()

    for i, (mol_id, smiles, formula) in enumerate(rows):
        rec = classify_molecule(mol_id, smiles, formula)
        results.append(rec)
        if rec["error"]:
            errors += 1
        if (i + 1) % 500 == 0 or (i + 1) == n_total:
            elapsed = time.perf_counter() - t_class_start
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            print(f"  {i+1:>7,}/{n_total:,}  ({rate:,.0f} mol/s)  errors so far: {errors}", end="\r")

    t_class_total = time.perf_counter() - t_class_start
    print(f"\n  Done in {t_class_total:.1f}s  ({n_total / t_class_total:.0f} mol/s)  errors: {errors}")

    # Summary stats
    valid = [r for r in results if r["time_ms"] is not None]
    classified = [r for r in valid if r["classified"]]
    f_group = [r for r in valid if r["has_fluorine_group"]]
    times = [r["time_ms"] for r in valid]

    def stats(vals):
        if not vals:
            return {}
        n = len(vals)
        vals_s = sorted(vals)
        mu = sum(vals) / n
        var = sum((v - mu)**2 for v in vals) / n
        import math
        return {
            "n": n, "mean_ms": mu, "std_ms": math.sqrt(var),
            "median_ms": vals_s[n // 2],
            "min_ms": vals_s[0], "max_ms": vals_s[-1],
            "p95_ms": vals_s[int(0.95 * n)], "total_ms": sum(vals),
        }

    # Timing by size bracket
    bracket_times: Dict[str, list] = {}
    for r in valid:
        b = atom_bracket(r["n_heavy_atoms"])
        bracket_times.setdefault(b, []).append(r["time_ms"])
    timing_by_bracket = {b: stats(ts) for b, ts in bracket_times.items()}

    # Timing by halogen class (F only / Cl, Br, I / mixed / none)
    halogen_cls_times: Dict[str, list] = {}
    for r in valid:
        nf = r["n_fluorine"] or 0
        nh = (r["n_halogens"] or 0) - nf
        if nf > 0 and nh == 0:
            cls = "F only"
        elif nf > 0 and nh > 0:
            cls = "Mixed (F + other)"
        elif nf == 0 and nh > 0:
            cls = "Other halogen (no F)"
        else:
            cls = "No halogen"
        halogen_cls_times.setdefault(cls, []).append(r["time_ms"])
    timing_by_halogen_class = {c: stats(ts) for c, ts in halogen_cls_times.items()}

    summary = {
        "tool": "PFASGroups",
        "timestamp": timestamp,
        "db": {"name": args.db_name, "host": args.db_host, "port": args.db_port, "user": args.db_user},
        "n_molecules_fetched": n_total,
        "n_valid": len(valid),
        "n_errors": errors,
        "n_classified": len(classified),
        "n_with_fluorine_group": len(f_group),
        "classification_rate_pct": 100.0 * len(classified) / len(valid) if valid else 0,
        "fluorine_group_rate_pct": 100.0 * len(f_group) / len(valid) if valid else 0,
        "total_classification_time_s": t_class_total,
        "timing_overall": stats(times),
        "timing_by_atom_bracket": timing_by_bracket,
        "timing_by_halogen_class": timing_by_halogen_class,
        "molecules": results,
    }

    # Save
    data_dir = BENCHMARK_DIR / "data"
    data_dir.mkdir(exist_ok=True)
    out_path = args.output or (data_dir / f"PFASGroups_clinventory_{timestamp}.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to: {out_path}")

    # Quick human-readable summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Molecules processed   : {n_total:,}")
    print(f"  Valid (parsed)        : {len(valid):,}")
    print(f"  Classified (any group): {len(classified):,}  ({100*len(classified)/max(len(valid),1):.1f}%)")
    print(f"  Has F-containing group: {len(f_group):,}  ({100*len(f_group)/max(len(valid),1):.1f}%)")
    if times:
        s = stats(times)
        print(f"  Timing median         : {s['median_ms']:.3f} ms/mol")
        print(f"  Timing mean ± std     : {s['mean_ms']:.3f} ± {s['std_ms']:.3f} ms/mol")
        print(f"  Total wall time       : {t_class_total:.1f} s")
    print(f"\nOutput: {out_path}")


if __name__ == "__main__":
    main()
