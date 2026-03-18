"""assess_valence_fragmentation.py
----------------------------------
Count how many molecules trigger ``fragment_until_valence_is_correct``
during sanitisation.  By default reads from the ECHA CLP Excel file;
pass ``--db`` to switch to the legacy clinventory PostgreSQL database.

Usage
-----
    # Default: read from Excel file
    python benchmark/scripts/assess_valence_fragmentation.py

    # Explicit file path / column
    python benchmark/scripts/assess_valence_fragmentation.py \\
        --file /path/to/substances.xlsx --col smiles

    # Legacy: use the PostgreSQL database
    python benchmark/scripts/assess_valence_fragmentation.py --db

File source (default)
---------------------
    --file      path to an Excel (.xlsx/.xls) or CSV file
                (default: <repo>/zeropm_db/database/lists/data_lists/echa_clp_with_smiles.xlsx)
    --col       SMILES column name in the file   (default: smiles)
    --limit     max rows to process, 0 = all     (default: 0)
    --batch     rows processed per parse call    (default: 5000)

DB source (only when --db is passed)
-------------------------------------
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
from pathlib import Path
from collections import Counter

# Default Excel file shipped with the zeropm_db repo
_REPO_ROOT = Path(__file__).resolve().parents[3]  # …/PFASGroups/../
_DEFAULT_FILE = (
    _REPO_ROOT / "zeropm_db" / "database" / "lists" / "data_lists"
    / "echa_clp_with_smiles.xlsx"
)

from rdkit import rdBase
rdBase.DisableLog("rdApp.warning")
rdBase.DisableLog("rdApp.error")
from getpass import getpass


# ── local import ──────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from PFASGroups.parser import parse_smiles  # noqa: E402


# ─────────────────────────────────────────────────────────────────────────────
# File helpers
# ─────────────────────────────────────────────────────────────────────────────

def _iter_smiles_from_file(filepath: str | Path, col: str, limit: int):
    """Yield individual SMILES strings from an Excel or CSV file."""
    try:
        import pandas as pd
    except ImportError as exc:
        raise SystemExit("pandas is required: pip install pandas openpyxl") from exc

    filepath = Path(filepath)
    suffix = filepath.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        df = pd.read_excel(filepath, usecols=[col])
    elif suffix == ".csv":
        df = pd.read_csv(filepath, usecols=[col])
    else:
        raise SystemExit(f"Unsupported file format: {suffix!r}  (use .xlsx, .xls, or .csv)")

    series = df[col].dropna()
    if limit > 0:
        series = series.iloc[:limit]
    yield from series.astype(str)


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
    pwd  = getpass('DB password (leave empty for none): ') if os.getenv("CLINVENTORY_DB_PASSWORD") is None else os.getenv("CLINVENTORY_DB_PASSWORD")

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

def assess(file: str | Path | None = None,
           col: str = "smiles",
           batch: int = 5_000,
           limit: int = 0,
           use_db: bool = False,
           table: str = "substances") -> dict:
    """Scan a substance list and report fragmentation statistics.

    Parameters
    ----------
    file : path-like or None
        Path to an Excel (.xlsx/.xls) or CSV file containing SMILES.
        Ignored when *use_db* is True.  Defaults to the ECHA CLP file
        shipped with the zeropm_db repository.
    col : str
        Column containing SMILES (applies to both file and DB sources).
    batch : int
        Rows processed per :func:`parse_smiles` call.
    limit : int
        Maximum molecules to process (0 = all).
    use_db : bool
        When True, read from the clinventory PostgreSQL database instead
        of a file (legacy behaviour).
    table : str
        DB table name — only used when *use_db* is True.

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
    if use_db:
        engine = _build_engine()
        print(f"Connected to:  {engine.url!s}")
        table  = os.getenv("CLINVENTORY_SMILES_TABLE", table)
        col    = os.getenv("CLINVENTORY_SMILES_COL",   col)
        batch  = int(os.getenv("CLINVENTORY_BATCH_SIZE", batch))
        limit  = int(os.getenv("CLINVENTORY_LIMIT",      limit))
        smiles_iter = _iter_smiles(engine, table, col, batch, limit)
        source_label = f"DB table {table!r}  |  SMILES column: {col!r}"
    else:
        if file is None:
            file = _DEFAULT_FILE
        file = Path(file)
        if not file.exists():
            raise SystemExit(f"File not found: {file}")
        smiles_iter = _iter_smiles_from_file(file, col, limit)
        source_label = f"File: {file}  |  SMILES column: {col!r}"

    n_total      = 0
    n_invalid    = 0
    n_fragmented = 0
    n_clean      = 0
    all_events: list[dict] = []
    atom_counter: Counter  = Counter()

    t0 = time.perf_counter()
    print(f"{source_label}  |  batch={batch:,}  |  limit={'all' if limit == 0 else limit}\n")

    batch_smiles: list[str] = []

    def _process_batch(smis: list[str]) -> None:
        """Call parse_smiles(verbose=True) on a batch and accumulate stats."""
        nonlocal n_invalid, n_fragmented, n_clean
        import warnings as _warnings
        with _warnings.catch_warnings(record=True):
            _warnings.simplefilter("always")
            _, verbose_info = parse_smiles(
                smis,
                verbose=True,
                halogens=None,              # assess all compounds, not just PFAS
                compute_component_metrics=False,  # skip graph metrics for speed
            )
        n_frag_batch = len(verbose_info['fragmented'])
        n_inv_batch  = verbose_info['n_invalid']
        n_invalid    += n_inv_batch
        n_fragmented += n_frag_batch
        n_clean      += len(smis) - n_inv_batch - n_frag_batch
        for fe in verbose_info['fragmented']:
            for ev in fe['events']:
                all_events.append(ev)
                atom_counter[ev.get('atom_symbol', 'unknown')] += 1

    for smi in smiles_iter:
        batch_smiles.append(smi)
        n_total += 1

        if len(batch_smiles) >= batch:
            _process_batch(batch_smiles)
            batch_smiles = []

        if n_total % 10_000 == 0:
            elapsed = time.perf_counter() - t0
            rate    = n_total / elapsed
            print(f"  {n_total:>8,}  processed  |  fragmented so far: "
                  f"{n_fragmented:,} ({100*n_fragmented/n_total:.2f}%)  "
                  f"|  {rate:,.0f} mol/s")

    # flush the last (possibly partial) batch
    if batch_smiles:
        _process_batch(batch_smiles)

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
        description=(
            "Assess how many molecules need valence fragmentation.  "
            "By default reads from the ECHA CLP Excel file; use --db to "
            "switch to the clinventory PostgreSQL database."
        )
    )
    # ── source selection ───────────────────────────────────────────────
    src = parser.add_mutually_exclusive_group()
    src.add_argument(
        "--file", default=None, metavar="PATH",
        help=(
            f"Excel (.xlsx/.xls) or CSV file with SMILES "
            f"(default: {_DEFAULT_FILE})"
        ),
    )
    src.add_argument(
        "--db", action="store_true",
        help="Read from the clinventory PostgreSQL database instead of a file.",
    )
    # ── shared options ─────────────────────────────────────────────────
    parser.add_argument("--col",   default="smiles",
                        help="SMILES column name (default: smiles)")
    parser.add_argument("--batch", type=int, default=5_000,
                        help="Rows per processing batch (default: 5000)")
    parser.add_argument("--limit", type=int, default=0,
                        help="Max molecules to process, 0 = all (default: 0)")
    # ── DB-only options ────────────────────────────────────────────────
    parser.add_argument("--table", default="substances",
                        help="DB table name, only used with --db (default: substances)")
    args = parser.parse_args()

    stats = assess(
        file=args.file,
        col=args.col,
        batch=args.batch,
        limit=args.limit,
        use_db=args.db,
        table=args.table,
    )
    print_report(stats)
