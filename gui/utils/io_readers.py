"""
I/O readers — load compound datasets from various file formats.

Supported sources
-----------------
  CSV          : any delimiter auto-detected by pandas
  Excel        : .xlsx / .xls (requires openpyxl)
  SQLite       : local .db / .sqlite file (stdlib sqlite3)
  MariaDB/MySQL: via SQLAlchemy + PyMySQL
  PostgreSQL   : via SQLAlchemy + psycopg2

All readers return a ``pandas.DataFrame``.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_csv(path: str | Path) -> pd.DataFrame:
    """Read a CSV file, auto-detecting delimiter."""
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def read_excel(path: str | Path, sheet: str | int = 0) -> pd.DataFrame:
    """Read an Excel file (.xlsx / .xls)."""
    try:
        import openpyxl  # noqa: F401
    except ImportError as exc:
        raise ImportError("openpyxl is required to read Excel files. "
                          "Install it with:  pip install openpyxl") from exc
    return pd.read_excel(path, sheet_name=sheet, dtype=str)


def read_sqlite(path: str | Path, query: str = "", table: str = "") -> pd.DataFrame:
    """Read from a SQLite database file.

    Supply either ``query`` (arbitrary SQL) or ``table`` (table name).
    """
    import sqlite3

    conn = sqlite3.connect(str(path))
    try:
        sql = query or f'SELECT * FROM "{table}"'
        df = pd.read_sql_query(sql, conn)
    finally:
        conn.close()
    return df.astype(str)


def read_database(
    dialect: str,        # "mariadb", "postgresql", "sqlite"
    host: str = "",
    port: int = 0,
    database: str = "",
    username: str = "",
    password: str = "",
    query: str = "",
    table: str = "",
    sqlite_path: str = "",
) -> pd.DataFrame:
    """Read from a remote or local database via SQLAlchemy.

    For SQLite supply ``sqlite_path`` instead of host/port/credentials.
    """
    try:
        from sqlalchemy import create_engine, text
    except ImportError as exc:
        raise ImportError(
            "sqlalchemy is required for database connections. "
            "Install it with:  pip install sqlalchemy"
        ) from exc

    if dialect == "sqlite":
        url = f"sqlite:///{sqlite_path}"
    elif dialect in ("mariadb", "mysql"):
        try:
            import pymysql  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "PyMySQL is required for MariaDB/MySQL connections. "
                "Install it with:  pip install PyMySQL"
            ) from exc
        url = (f"mysql+pymysql://{username}:{password}@{host}:{port or 3306}"
               f"/{database}")
    elif dialect == "postgresql":
        try:
            import psycopg2  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "psycopg2 is required for PostgreSQL connections. "
                "Install it with:  pip install psycopg2-binary"
            ) from exc
        url = (f"postgresql+psycopg2://{username}:{password}@{host}:{port or 5432}"
               f"/{database}")
    else:
        raise ValueError(f"Unsupported dialect: {dialect!r}")

    engine = create_engine(url)
    sql = query or f'SELECT * FROM "{table}"'
    with engine.connect() as conn:
        df = pd.read_sql_query(text(sql), conn)
    return df.astype(str)


# ---------------------------------------------------------------------------
# Auto-detect and dispatch
# ---------------------------------------------------------------------------

def get_excel_sheets(path: str | Path) -> list[str]:
    """Return the list of sheet names in an Excel workbook without reading data."""
    try:
        import openpyxl
        wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
        sheets = wb.sheetnames
        wb.close()
        return sheets
    except Exception:
        return []


def read_smiles_file(path: str | Path) -> pd.DataFrame:
    """Read a plain SMILES file (.smi / .smiles / .txt).

    Each non-empty line is expected to be either::

        SMILES
        SMILES name
        name SMILES       (if second token looks like a SMILES string)

    Returns a DataFrame with columns ``smiles`` and optionally ``name``.
    """
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Split on whitespace (tab, space)
            parts = line.split(None, 1)
            if len(parts) == 1:
                rows.append({"smiles": parts[0]})
            else:
                # Heuristic: the longer / more bracket-containing token is the SMILES
                a, b = parts[0], parts[1]
                if len(a) >= len(b) or any(c in a for c in "()[]@#"):
                    rows.append({"smiles": a, "name": b})
                else:
                    rows.append({"smiles": b, "name": a})
    return pd.DataFrame(rows, columns=[c for c in ["name", "smiles"] if c in (rows[0] if rows else {})])


def read_file(path: str | Path, sheet: str | int | None = None) -> pd.DataFrame:
    """Dispatch to the appropriate reader based on file extension."""
    p = Path(path)
    ext = p.suffix.lower()
    if ext == ".csv":
        return read_csv(p)
    if ext in (".xlsx", ".xls", ".xlsm"):
        return read_excel(p, sheet=sheet)
    if ext in (".db", ".sqlite", ".sqlite3"):
        return read_sqlite(p)
    if ext in (".smi", ".smiles", ".txt"):
        return read_smiles_file(p)
    # Fallback: try CSV
    return read_csv(p)


def get_smiles_column(df: pd.DataFrame) -> str | None:
    """Heuristically identify the SMILES column in a DataFrame."""
    candidates = ["smiles", "SMILES", "Smiles", "canonical_smiles",
                  "CanonicalSMILES", "smi", "SMI", "structure"]
    for c in candidates:
        if c in df.columns:
            return c
    # Try any column that contains mostly valid-looking SMILES characters
    for col in df.columns:
        sample = df[col].dropna().head(5).astype(str)
        if sample.str.match(r"^[A-Za-z0-9@\[\]()=#\-\\\/\.\+%]+$").mean() > 0.6:
            return col
    return None


def get_name_column(df: pd.DataFrame) -> str | None:
    """Heuristically identify a compound name column."""
    candidates = ["name", "Name", "compound_name", "CompoundName",
                  "substance_name", "preferred_name", "dtxsid", "DTXSID",
                  "id", "ID", "CAS"]
    for c in candidates:
        if c in df.columns:
            return c
    return None
