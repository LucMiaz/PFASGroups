#!/usr/bin/env python3
"""
Build an ML-ready dataset from the ToxCast MariaDB database.

Steps
-----
1. Pull chemicals with SMILES from the toxcast.chemical table.
2. Generate PFASGroups binary fingerprints (115 columns) for every chemical.
3. Build a binary label matrix from mc5 hit-calls for a curated set of
   high-coverage, biologically meaningful endpoints.
4. Save features (X) and labels (y) as parquet files next to this script
   under ../data/toxcast_*.parquet, plus a CSV summary.

Usage
-----
    conda activate chem
    cd benchmark
    python scripts/build_toxcast_dataset.py
"""

from __future__ import annotations

import sys
from getpass import getpass
from pathlib import Path

import numpy as np
import pandas as pd
from sqlalchemy import create_engine, text
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
DATA_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Curated endpoint list
# aeid → human-readable label (name kept short for column headers)
# Criteria: burst_assay=0, cell_viability_assay=0, n_chems ≥ 7000,
#            hit-rate 2–50 %, biologically interpretable target.
# ---------------------------------------------------------------------------
ENDPOINTS = {
    762:  "AR_antagonist",          # androgen receptor antagonism
    761:  "AR_agonist",             # androgen receptor agonism
    786:  "ERa_antagonist",         # estrogen receptor α antagonism
    788:  "ERa_agonist",            # estrogen receptor α agonism
    806:  "AhR_agonist",            # aryl hydrocarbon receptor
    767:  "Aromatase_antagonist",   # CYP19 aromatase inhibition
    804:  "TR_antagonist",          # thyroid receptor antagonism
    1134: "DT40_genotoxicity",      # DNA damage (DT40 cell line)
    1854: "MMP_ratio",              # mitochondrial membrane potential
    3185: "CYP2D6_antagonist",      # CYP2D6 inhibition
    3186: "CYP2C19_antagonist",     # CYP2C19 inhibition
    3187: "CYP2C9_antagonist",      # CYP2C9 inhibition
    2544: "CYP3A4_antagonist",      # CYP3A4 inhibition
    1116: "p53_ratio",              # p53 activation (DNA stress)
    2366: "Caspase3_HEPG2",        # apoptosis (caspase-3, HepG2)
}

HITC_THRESHOLD = 0.9  # ToxCast convention: hitc ≥ 0.9 → active


def get_engine(user: str, password: str) -> object:
    return create_engine(
        f"mariadb+mariadbconnector://{user}:{password}@127.0.0.1/toxcast",
        pool_pre_ping=True,
    )


# ---------------------------------------------------------------------------
# Step 1 – chemicals
# ---------------------------------------------------------------------------

def load_chemicals(engine) -> pd.DataFrame:
    """Return all chemicals that have a SMILES string."""
    query = text(
        "SELECT chid, casn, chnm, dsstox_substance_id, smiles "
        "FROM chemical "
        "WHERE smiles IS NOT NULL AND smiles != ''"
    )
    with engine.connect() as conn:
        df = pd.read_sql(query, conn)
    print(f"[chemicals] {len(df):,} chemicals with SMILES")
    return df


# ---------------------------------------------------------------------------
# Step 2 – PFASGroups fingerprints
# ---------------------------------------------------------------------------

def build_fingerprint_matrix(
    smiles_list: list[str],
    halogens: str = "F",
) -> tuple[np.ndarray, list[str]]:
    """
    Generate binary PFASGroups fingerprints.

    Returns
    -------
    X          : np.ndarray, shape (n, 115), dtype uint8
                 115 = all compute=True groups (116 after excluding the
                 aggregate-telomers group) minus the fluoride ionic halide
                 group excluded when scanning F only.
    group_names: list[str] of length 115
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from PFASGroups import generate_fingerprint  # noqa: PLC0415

    # Probe to get dimensions
    _probe, _info = generate_fingerprint(
        "FC(F)(F)C(=O)O",
        representation="vector",
        count_mode="binary",
        halogens=halogens,
        saturation=None,
    )
    n_groups    = len(_probe)
    group_names = _info["group_names"]

    X      = np.zeros((len(smiles_list), n_groups), dtype=np.uint8)
    errors = 0
    for i, smi in enumerate(tqdm(smiles_list, desc="PFASGroups fingerprints")):
        try:
            fp, _ = generate_fingerprint(
                smi,
                representation="vector",
                count_mode="binary",
                halogens=halogens,
                saturation=None,
            )
            X[i] = fp
        except Exception:  # noqa: BLE001
            errors += 1

    if errors:
        print(f"  Warning: {errors} molecules failed fingerprint generation")

    pfas_count = int((X.sum(axis=1) > 0).sum())
    print(
        f"[fingerprints] shape={X.shape}  "
        f"PFAS (≥1 group): {pfas_count} ({100*pfas_count/len(smiles_list):.1f}%)"
    )
    return X, group_names


# ---------------------------------------------------------------------------
# Step 3 – label matrix
# ---------------------------------------------------------------------------

def build_label_matrix(engine, chids: list[int]) -> pd.DataFrame:
    """
    Build a (n_chems × n_endpoints) binary hit-call matrix.

    Only the most recent mc5 result per chemical×endpoint is kept
    (MAX(m5id) tie-break).
    """
    aeid_list = ",".join(str(a) for a in ENDPOINTS)
    chid_list = ",".join(str(c) for c in chids)

    query = text(f"""
        SELECT s.chid,
               m5.aeid,
               MAX(CASE WHEN m5.hitc >= {HITC_THRESHOLD} THEN 1 ELSE 0 END) AS active
        FROM mc5 m5
        JOIN mc4  m4 USING (m4id)
        JOIN sample s USING (spid)
        WHERE m5.aeid IN ({aeid_list})
          AND s.chid   IN ({chid_list})
        GROUP BY s.chid, m5.aeid
    """)

    with engine.connect() as conn:
        df = pd.read_sql(query, conn)

    # Pivot to wide format: rows = chid, columns = endpoint label
    label_map = {aeid: label for aeid, label in ENDPOINTS.items()}
    df["endpoint"] = df["aeid"].map(label_map)
    y = df.pivot_table(index="chid", columns="endpoint", values="active", aggfunc="max")
    y = y.reindex(columns=list(ENDPOINTS.values()))  # consistent column order

    # Coverage report
    print("[labels] coverage per endpoint:")
    for col in y.columns:
        n_tested = y[col].notna().sum()
        n_active = y[col].eq(1).sum()
        print(f"  {col:<30s}  tested={n_tested:5d}  active={n_active:5d}  ({100*n_active/max(n_tested,1):.1f}%)")

    return y


# ---------------------------------------------------------------------------
# Step 4 – assemble & save
# ---------------------------------------------------------------------------

def build_and_save(user: str, password: str) -> None:
    engine = get_engine(user, password)

    # 1 – chemicals
    chem_df = load_chemicals(engine)

    # 2 – fingerprints
    X, group_names = build_fingerprint_matrix(chem_df["smiles"].tolist())

    fp_df = pd.DataFrame(X, columns=group_names, dtype=np.uint8)
    fp_df.insert(0, "chid",               chem_df["chid"].values)
    fp_df.insert(1, "casn",               chem_df["casn"].values)
    fp_df.insert(2, "chnm",               chem_df["chnm"].values)
    fp_df.insert(3, "dsstox_substance_id", chem_df["dsstox_substance_id"].values)
    fp_df.insert(4, "smiles",              chem_df["smiles"].values)

    fp_path = DATA_DIR / "toxcast_fingerprints.parquet"
    fp_df.to_parquet(fp_path, index=False)
    print(f"[saved] features → {fp_path}")

    # 3 – labels
    y = build_label_matrix(engine, chem_df["chid"].tolist())

    # Merge fingerprints + labels on chid, keep only rows that have ≥1 label
    merged = fp_df.merge(
        y.reset_index(),
        on="chid",
        how="left",
    )
    # Report overall matrix density
    label_cols = list(ENDPOINTS.values())
    y_sub = merged[label_cols]
    tested = y_sub.notna().sum().sum()
    active = y_sub.eq(1).sum().sum()
    total  = y_sub.size
    print(
        f"[matrix] {len(merged):,} chemicals × {len(label_cols)} endpoints  "
        f"fill={100*tested/total:.1f}%  "
        f"overall hit-rate={100*active/max(tested,1):.1f}%"
    )

    dataset_path = DATA_DIR / "toxcast_dataset.parquet"
    merged.to_parquet(dataset_path, index=False)
    print(f"[saved] full dataset → {dataset_path}")

    # Compact CSV summary (metadata + labels only, no fingerprint columns)
    summary = merged[["chid", "casn", "chnm", "dsstox_substance_id", "smiles"] + label_cols]
    summary_path = DATA_DIR / "toxcast_labels_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[saved] label summary → {summary_path}")

    print("\nDone. Suggested next step:")
    print("  python scripts/train_toxcast_models.py")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    _user     = input("DB username: ")
    _password = getpass(f"DB password for {_user}: ")
    build_and_save(_user, _password)
