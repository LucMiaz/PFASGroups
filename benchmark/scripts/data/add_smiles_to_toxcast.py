"""
Populate the `smiles`, `inchi`, `inchikey` and `formula` columns of the
`chemical` table in a local ToxCast MariaDB database.

Strategy (run in order):
1. If ``--source-db`` is given, copy structures from that existing database
   for all chemicals matched by ``dsstox_substance_id``.  This avoids
   unnecessary PubChem API calls for chemicals that are already known.
2. For any chemical still lacking a SMILES string after step 1, fall back to
   the PubChem REST API (pubchem.dtxsid_to_mols).

Usage
-----
    # Fill from previous DB first, then PubChem for the rest:
    python add_smiles_to_toxcast.py --target invitrodb_v4_3 --source-db toxcast

    # PubChem only (original behaviour):
    python add_smiles_to_toxcast.py --target invitrodb_v4_3
"""

import argparse
from getpass import getpass

from pubchem import dtxsid_to_mols
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from sqlalchemy import create_engine, text
from tqdm import tqdm

_STRUCTURE_COLS = ("smiles", "inchi", "inchikey", "formula")


def _make_engine(user: str, password: str, db: str):
    return create_engine(
        f"mariadb+mariadbconnector://{user}:{password}@localhost/{db}",
        pool_pre_ping=True,
    )


def _ensure_structure_columns(engine) -> None:
    """Add smiles/inchi/inchikey/formula columns to `chemical` if they don't exist yet."""
    col_defs = {
        "smiles":    "TEXT",
        "inchi":     "TEXT",
        "inchikey":  "VARCHAR(27)",
        "formula":   "VARCHAR(100)",
    }
    with engine.connect() as conn:
        existing = {
            row[0]
            for row in conn.execute(
                text("SHOW COLUMNS FROM chemical")
            ).fetchall()
        }
        for col, col_type in col_defs.items():
            if col not in existing:
                conn.execute(text(f"ALTER TABLE chemical ADD COLUMN {col} {col_type}"))
                print(f"[schema] Added column '{col}' to chemical table")
        conn.commit()


def _write_structures(conn, rows: list[dict]) -> None:
    """UPDATE structure columns for a list of dicts with keys matching _STRUCTURE_COLS + dtxsid."""
    for row in rows:
        conn.execute(
            text(
                "UPDATE chemical"
                " SET smiles=:smiles, inchi=:inchi, inchikey=:inchikey, formula=:formula"
                " WHERE dsstox_substance_id=:dtxsid"
            ),
            row,
        )


# ---------------------------------------------------------------------------
# Step 1 – copy from an existing database
# ---------------------------------------------------------------------------

def fill_from_source_db(target_engine, source_engine) -> list[str]:
    """Copy structures from *source_engine* into *target_engine*.

    Returns the list of DTXSIDs that were successfully populated so the caller
    knows which chemicals still need PubChem resolution.
    """
    # Fetch all structure data from the source DB (indexed by DTXSID).
    with source_engine.connect() as conn:
        src_rows = conn.execute(
            text(
                "SELECT dsstox_substance_id, smiles, inchi, inchikey, formula"
                " FROM chemical"
                " WHERE smiles IS NOT NULL AND smiles != ''"
            )
        ).fetchall()

    src_map = {
        row[0]: {
            "smiles":    row[1],
            "inchi":     row[2],
            "inchikey":  row[3],
            "formula":   row[4],
        }
        for row in src_rows
    }
    print(f"[source DB] {len(src_map):,} chemicals with SMILES available")

    # Find chemicals in the target DB that are still missing a SMILES.
    with target_engine.connect() as conn:
        missing = conn.execute(
            text(
                "SELECT dsstox_substance_id FROM chemical"
                " WHERE smiles IS NULL OR smiles = ''"
            )
        ).fetchall()

    missing_dtxsids = [row[0] for row in missing if row[0]]
    print(f"[target DB] {len(missing_dtxsids):,} chemicals need a SMILES")

    filled: list[str] = []
    batch: list[dict] = []

    for dtxsid in tqdm(missing_dtxsids, desc="Copying from source DB"):
        if dtxsid not in src_map:
            continue
        batch.append({"dtxsid": dtxsid, **src_map[dtxsid]})
        if len(batch) >= 500:
            with target_engine.connect() as conn:
                _write_structures(conn, batch)
                conn.commit()
            filled.extend(b["dtxsid"] for b in batch)
            batch.clear()

    if batch:
        with target_engine.connect() as conn:
            _write_structures(conn, batch)
            conn.commit()
        filled.extend(b["dtxsid"] for b in batch)

    print(f"[source DB] populated {len(filled):,} chemicals from source DB")
    return filled


# ---------------------------------------------------------------------------
# Step 2 – fill remaining chemicals via PubChem
# ---------------------------------------------------------------------------

def fill_from_pubchem(target_engine, batch_size: int = 10) -> None:
    """Fetch structure data from PubChem for chemicals still missing a SMILES."""
    with target_engine.connect() as conn:
        rows = conn.execute(
            text(
                "SELECT chid, dsstox_substance_id FROM chemical"
                " WHERE smiles IS NULL OR smiles = ''"
            )
        ).fetchall()

    if not rows:
        print("[PubChem] Nothing left to fetch – all chemicals have a SMILES.")
        return

    print(f"[PubChem] Fetching structures for {len(rows):,} remaining chemicals …")

    for i in tqdm(range(0, len(rows), batch_size), desc="PubChem batches"):
        batch = rows[i : i + batch_size]
        dtxsid_to_chid = {row[1]: row[0] for row in batch if row[1]}
        dtxsids = list(dtxsid_to_chid.keys())

        mols_dict, failed = dtxsid_to_mols(dtxsids)
        if failed:
            tqdm.write(f"  Failed DTXSIDs: {failed}")

        write_batch: list[dict] = []
        for dtxsid, mol_list in mols_dict.items():
            if not mol_list:
                continue
            mol = mol_list[0]
            if mol is None:
                continue
            write_batch.append(
                {
                    "dtxsid":   dtxsid,
                    "smiles":   Chem.MolToSmiles(mol),
                    "inchi":    Chem.MolToInchi(mol),
                    "inchikey": Chem.MolToInchiKey(mol),
                    "formula":  CalcMolFormula(mol),
                }
            )

        if write_batch:
            with target_engine.connect() as conn:
                _write_structures(conn, write_batch)
                conn.commit()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Populate structure columns in a ToxCast-schema MariaDB database."
    )
    parser.add_argument(
        "--target", "-t",
        default=None,
        help="Target database name (default: 'toxcast').",
    )
    parser.add_argument(
        "--source-db", "-s",
        default=None,
        metavar="DB_NAME",
        help=(
            "Name of an existing database to copy structures from before "
            "falling back to PubChem.  Uses the same host/credentials."
        ),
    )
    parser.add_argument(
        "--batch-size", "-b",
        type=int,
        default=10,
        help="Number of DTXSIDs per PubChem API call (default: 10).",
    )
    parser.add_argument("--user",     "-u", default=None, help="DB username.")
    parser.add_argument("--password", "-p", default=None, help="DB password.")
    args = parser.parse_args()

    user     = args.user     or input("Enter username for database connection: ")
    password = args.password or getpass(f"Enter password for database connection for {user}: ")
    target   = args.target   or "toxcast"

    target_engine = _make_engine(user, password, target)
    _ensure_structure_columns(target_engine)

    if args.source_db:
        source_engine = _make_engine(user, password, args.source_db)
        fill_from_source_db(target_engine, source_engine)

    fill_from_pubchem(target_engine, batch_size=args.batch_size)
    print("Done.")


if __name__ == "__main__":
    main()
