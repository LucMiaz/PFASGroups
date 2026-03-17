"""assess_valence_fragmentation.py
----------------------------------
Count how many molecules in the clinventory database trigger
``fragment_until_valence_is_correct`` during sanitisation.

Usage
-----
    python benchmark/assess_valence_fragmentation.py

Environment variables (all optional, fall back to localhost defaults)
---------------------------------------------------------------------
    CLINVENTORY_DB_URL      full SQLAlchemy connection string, e.g.
                            ``postgresql://user:pwd@host:5432/clinventory``
    CLINVENTORY_DB_NAME     database name          (default: clinventory)
    CLINVENTORY_DB_HOST     host                   (default: localhost)
    CLINVENTORY_DB_PORT     port                   (default: 5432)
    CLINVENTORY_DB_USER     username               (default: current OS user)
    CLINVENTORY_DB_PASSWORD password               (default: empty)
    CLINVENTORY_SMILES_TABLE  table name           (default: substances)
    CLINVENTORY_SMILES_COL  column with SMILES     (default: smiles)
    CLINVENTORY_BATCH_SIZE  rows fetched at a time (default: 5000)
    CLINVENTORY_LIMIT       max rows to process, 0 = all  (default: 0)
"""

from __future__ import annotations

import os
import sys
import time
import textwrap
from collections import Counter

from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog("rdApp.warning")
rdBase.DisableLog("rdApp.error")

# ── local import ──────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from PFASGroups.core import fragment_until_valence_is_correct  # noqa: E402


# ─────────────────────────────────────────────────────────────────────────────
# DB helpers
# ─────────────────────────────────────────────────────────────────────────────

def _build_engine():
    """Return a SQLAlchemy engine for the clinventory database."""
    try:
        from sqlalchemy import create_engine
    except ImportError as exc:
        raise SystemExit("SQLAlchemy is required: pip install sqlalchemy psycopg2-binary") from exc

    url = os.getenv("CLINVENTORY_DB_URL")
    if url:
        return create_engine(url)

    db   = os.getenv("CLINVENTORY_DB_NAME",     "clinventory")
    host = os.getenv("CLINVENTORY_DB_HOST",     "localhost")
    port = os.getenv("CLINVENTORY_DB_PORT",     "5432")
    user = os.getenv("CLINVENTORY_DB_USER",     os.getenv("USERNAME", os.getenv("USER", "postgres")))
    pwd  = os.getenv("CLINVENTORY_DB_PASSWORD", "")

    # prefer Unix socket when host is localhost and no password is set
    if host == "localhost" and not pwd:
        try:
            return create_engine(f"postgresql+psycopg2://{user}@/{db}")
        except Exception:
            pass

    url = f"postgresql+psycopg2://{user}:{pwd}@{host}:{port}/{db}"
    return create_engine(url)


def _iter_smiles(engine, table: str, col: str, batch: int, limit: int):
    """Yield individual SMILES strings from the database in batches."""
    from sqlalchemy import text
    limit_clause = f"LIMIT {limit}" if limit > 0 else ""
    query = text(f'SELECT "{col}" FROM "{table}" WHERE "{col}" IS NOT NULL {limit_clause}')
    with engine.connect() as conn:
        result = conn.execute(query)
        while True:
            rows = result.fetchmany(batch)
            if not rows:
                break
            for row in rows:
                yield str(row[0])


# ─────────────────────────────────────────────────────────────────────────────
# Analysis
# ─────────────────────────────────────────────────────────────────────────────

def assess(table: str = "substances",
           col: str = "smiles",
           batch: int = 5_000,
           limit: int = 0) -> dict:
    """Scan the clinventory and report fragmentation statistics.

    Parameters
    ----------
    table : str
        Table to query.
    col : str
        Column containing SMILES.
    batch : int
        Rows fetched per round-trip.
    limit : int
        Maximum molecules to process (0 = all).

    Returns
    -------
    dict with keys:

    ``n_total``
        Molecules attempted.
    ``n_invalid``
        Molecules RDKit could not parse at all.
    ``n_fragmented``
        Molecules that triggered ``fragment_until_valence_is_correct``.
    ``n_clean``
        Molecules that sanitised without error.
    ``events``
        List of verbose event dicts for each fragmentation step.
    ``atom_counter``
        :class:`Counter` of problematic atom indices (element symbol).
    """
    engine = _build_engine()
    print(f"Connected to:  {engine.url!s}")

    table    = os.getenv("CLINVENTORY_SMILES_TABLE", table)
    col      = os.getenv("CLINVENTORY_SMILES_COL",   col)
    batch    = int(os.getenv("CLINVENTORY_BATCH_SIZE", batch))
    limit    = int(os.getenv("CLINVENTORY_LIMIT",      limit))

    n_total      = 0
    n_invalid    = 0
    n_fragmented = 0
    n_clean      = 0
    all_events: list[dict] = []
    atom_counter: Counter  = Counter()

    t0 = time.perf_counter()
    print(f"Table: {table!r}  |  SMILES column: {col!r}  "
          f"|  batch={batch:,}  |  limit={'all' if limit == 0 else limit:,}\n")

    for smi in _iter_smiles(engine, table, col, batch, limit):
        n_total += 1
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            n_invalid += 1
            continue

        mol = Chem.AddHs(mol)
        try:
            Chem.SanitizeMol(mol)
            n_clean += 1
        except Chem.AtomValenceException:
            # Use verbose=True to collect event details
            _frags, events = fragment_until_valence_is_correct(mol, [], verbose=True)
            if events:
                n_fragmented += 1
                for ev in events:
                    all_events.append(ev)
                    # Resolve atom index → element symbol from the original mol
                    try:
                        sym = mol.GetAtomWithIdx(ev["atom_idx"]).GetSymbol()
                    except Exception:
                        sym = "unknown"
                    atom_counter[sym] += 1
            else:
                # fragment_until_valence returned empty events (raised internally)
                n_fragmented += 1

        if n_total % 10_000 == 0:
            elapsed = time.perf_counter() - t0
            rate    = n_total / elapsed
            print(f"  {n_total:>8,}  processed  |  fragmented so far: "
                  f"{n_fragmented:,} ({100*n_fragmented/n_total:.2f}%)  "
                  f"|  {rate:,.0f} mol/s")

    elapsed = time.perf_counter() - t0
    print()
    return dict(
        n_total=n_total,
        n_invalid=n_invalid,
        n_fragmented=n_fragmented,
        n_clean=n_clean,
        events=all_events,
        atom_counter=atom_counter,
        elapsed_s=elapsed,
    )


def print_report(stats: dict) -> None:
    nt = stats["n_total"]
    nf = stats["n_fragmented"]
    ni = stats["n_invalid"]
    nc = stats["n_clean"]

    pct = lambda n: f"{100 * n / nt:.2f}%" if nt else "n/a"

    print("=" * 60)
    print("  fragment_until_valence_is_correct — clinventory audit")
    print("=" * 60)
    print(f"  Total molecules processed : {nt:>10,}")
    print(f"  Invalid SMILES (skipped)  : {ni:>10,}  ({pct(ni)})")
    print(f"  Clean sanitisation        : {nc:>10,}  ({pct(nc)})")
    print(f"  Needed valence fix        : {nf:>10,}  ({pct(nf)})")
    print(f"  Total fragmentation steps : {len(stats['events']):>10,}")
    print(f"  Elapsed                   : {stats['elapsed_s']:>10.1f} s")

    if stats["atom_counter"]:
        print()
        print("  Problematic atoms (element → fragmentation events):")
        for sym, cnt in stats["atom_counter"].most_common():
            print(f"    {sym:>4s}  {cnt:,}")

    if stats["events"]:
        print()
        print("  Sample events (first 5):")
        for ev in stats["events"][:5]:
            smi_preview = (ev["smiles"] or "")[:60]
            print(textwrap.indent(
                f"atom_idx={ev['atom_idx']}  n_frags={ev['n_fragments']}\n"
                f"  error  : {ev['error']}\n"
                f"  SMILES : {smi_preview}",
                "    "
            ))
    print("=" * 60)


# ─────────────────────────────────────────────────────────────────────────────
# Entry-point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Assess how many clinventory molecules need valence fragmentation."
    )
    parser.add_argument("--table",  default="substances",
                        help="DB table name (default: substances)")
    parser.add_argument("--col",    default="smiles",
                        help="SMILES column (default: smiles)")
    parser.add_argument("--batch",  type=int, default=5_000,
                        help="Rows per DB fetch (default: 5000)")
    parser.add_argument("--limit",  type=int, default=0,
                        help="Max molecules to process, 0 = all (default: 0)")
    args = parser.parse_args()

    stats = assess(table=args.table, col=args.col, batch=args.batch, limit=args.limit)
    print_report(stats)
