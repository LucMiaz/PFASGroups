#!/usr/bin/env python3
"""
Compare ToxCast v4.2 (cached parquet) with v4.3 (invitro_v4_3 database).

Checks
------
1. Chemical coverage  – same set of DTXSID/CASN/SMILES?
2. Label matrix       – same hit-calls for shared chemicals?
3. C-F chemicals      – focused diff for fluorinated compounds.

Usage
-----
    conda activate chem
    cd benchmark
    python scripts/data/compare_toxcast_versions.py
"""

from __future__ import annotations

import sys
from getpass import getpass
from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine, text

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR   = Path(__file__).resolve().parent
DATA_DIR     = SCRIPT_DIR.parents[1] / "data"
PARQUET_PATH = DATA_DIR / "toxcast_dataset.parquet"

# ---------------------------------------------------------------------------
# Same endpoints as build_toxcast_dataset.py
# ---------------------------------------------------------------------------
ENDPOINTS = {
    762:  "AR_antagonist",
    761:  "AR_agonist",
    786:  "ERa_antagonist",
    788:  "ERa_agonist",
    806:  "AhR_agonist",
    767:  "Aromatase_antagonist",
    804:  "TR_antagonist",
    1134: "DT40_genotoxicity",
    1854: "MMP_ratio",
    3185: "CYP2D6_antagonist",
    3186: "CYP2C19_antagonist",
    3187: "CYP2C9_antagonist",
    2544: "CYP3A4_antagonist",
    1116: "p53_ratio",
    2366: "Caspase3_HEPG2",
}
LABEL_COLS   = list(ENDPOINTS.values())
HITC_THRESHOLD = 0.9

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_engine(user: str, password: str, db: str):
    return create_engine(
        f"mariadb+mariadbconnector://{user}:{password}@127.0.0.1/{db}",
        pool_pre_ping=True,
    )


def load_chemicals_from_db(engine) -> pd.DataFrame:
    """Pull chemicals with SMILES from a ToxCast-schema database."""
    q = text(
        "SELECT chid, casn, chnm, dsstox_substance_id, smiles "
        "FROM chemical "
        "WHERE smiles IS NOT NULL AND smiles != ''"
    )
    with engine.connect() as conn:
        return pd.read_sql(q, conn)


def load_labels_from_db(engine, chids: list[int]) -> pd.DataFrame:
    """Return a (chid x endpoint) binary hit-call DataFrame."""
    aeid_list = ",".join(str(a) for a in ENDPOINTS)
    chid_list = ",".join(str(c) for c in chids)
    q = text(f"""
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
        df = pd.read_sql(q, conn)

    label_map = {aeid: label for aeid, label in ENDPOINTS.items()}
    df["endpoint"] = df["aeid"].map(label_map)
    y = df.pivot_table(index="chid", columns="endpoint", values="active", aggfunc="max")
    y = y.reindex(columns=LABEL_COLS)
    return y


def is_cf_smiles(smiles: str) -> bool:
    """Heuristic: SMILES contains at least one C-F bond (capital F, attached to carbon)."""
    try:
        # Simple string heuristic – works for standard SMILES
        return "F" in smiles
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def compare(user: str, password: str) -> None:
    # ── Load v4.2 cache ──────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Loading v4.2 cache from {PARQUET_PATH.name} …")
    old = pd.read_parquet(PARQUET_PATH)
    print(f"  {len(old):,} rows, columns: {list(old.columns[:8])} …")

    meta_cols = ["chid", "casn", "chnm", "dsstox_substance_id", "smiles"]
    old_meta  = old[meta_cols].copy()
    old_labels = old.set_index("chid")[LABEL_COLS]

    # ── Load v4.3 from DB ────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Connecting to invitrodb_v4_3 …")
    engine43 = get_engine(user, password, "invitrodb_v4_3")
    new_chem  = load_chemicals_from_db(engine43)
    print(f"  {len(new_chem):,} chemicals with SMILES in v4.3")

    new_labels = load_labels_from_db(engine43, new_chem["chid"].tolist())
    new_meta   = new_chem[meta_cols].copy()

    # ── 1. Chemical coverage ─────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("1. CHEMICAL COVERAGE")
    old_dtxsids = set(old_meta["dsstox_substance_id"].dropna())
    new_dtxsids = set(new_meta["dsstox_substance_id"].dropna())
    only_old    = old_dtxsids - new_dtxsids
    only_new    = new_dtxsids - old_dtxsids
    shared      = old_dtxsids & new_dtxsids
    print(f"  v4.2 total  : {len(old_dtxsids):,}")
    print(f"  v4.3 total  : {len(new_dtxsids):,}")
    print(f"  Shared      : {len(shared):,}")
    print(f"  Only in v4.2: {len(only_old):,}")
    print(f"  Only in v4.3: {len(only_new):,}")

    # ── 2. SMILES changes for shared chemicals ───────────────────────────────
    print(f"\n{'='*60}")
    print("2. SMILES CHANGES FOR SHARED CHEMICALS")
    merged_meta = old_meta.merge(
        new_meta,
        on="dsstox_substance_id",
        suffixes=("_v42", "_v43"),
    )
    smiles_changed = merged_meta[
        merged_meta["smiles_v42"].fillna("") != merged_meta["smiles_v43"].fillna("")
    ]
    print(f"  Shared chemicals with different SMILES: {len(smiles_changed):,}")
    if not smiles_changed.empty:
        print(smiles_changed[["dsstox_substance_id", "chnm_v42", "smiles_v42", "smiles_v43"]].head(10).to_string())

    # ── 3. Label matrix diff (shared chemicals, all endpoints) ───────────────
    print(f"\n{'='*60}")
    print("3. LABEL MATRIX DIFF (shared chemicals only)")

    # Align by chid – use dsstox -> chid mapping
    dtxsid_to_chid_old = dict(zip(old_meta["dsstox_substance_id"], old_meta["chid"]))
    dtxsid_to_chid_new = dict(zip(new_meta["dsstox_substance_id"], new_meta["chid"]))

    shared_list = list(shared)
    old_chids   = [dtxsid_to_chid_old[d] for d in shared_list if d in dtxsid_to_chid_old]
    new_chids   = [dtxsid_to_chid_new[d] for d in shared_list if d in dtxsid_to_chid_new]

    ol = old_labels.loc[old_labels.index.isin(old_chids)].copy()
    nl = new_labels.loc[new_labels.index.isin(new_chids)].copy()

    # Re-index both on DTXSID for a fair comparison
    old_dtxsid_index = old_meta.set_index("chid")["dsstox_substance_id"]
    new_dtxsid_index = new_meta.set_index("chid")["dsstox_substance_id"]

    ol.index = ol.index.map(old_dtxsid_index)
    nl.index = nl.index.map(new_dtxsid_index)
    ol = ol[ol.index.isin(shared)]
    nl = nl[nl.index.isin(shared)]
    ol.sort_index(inplace=True)
    nl.sort_index(inplace=True)

    # Align on common index
    common_idx = ol.index.intersection(nl.index)
    ol_aligned = ol.loc[common_idx]
    nl_aligned = nl.loc[common_idx]

    print(f"  Chemicals with label data in both versions: {len(common_idx):,}")

    for col in LABEL_COLS:
        if col not in ol_aligned.columns or col not in nl_aligned.columns:
            print(f"  {col}: MISSING IN ONE VERSION")
            continue
        a = ol_aligned[col].fillna(-1)
        b = nl_aligned[col].fillna(-1)
        n_diff = (a != b).sum()
        n_old_tested = (ol_aligned[col] != -1).sum() if hasattr(ol_aligned[col], "__len__") else 0
        n_old_tested = ol_aligned[col].notna().sum()
        n_new_tested = nl_aligned[col].notna().sum()
        print(f"  {col:<30s}  tested_v42={n_old_tested:5d}  tested_v43={n_new_tested:5d}  changed={n_diff:4d}")

    # ── 4. C-F chemicals focused diff ────────────────────────────────────────
    print(f"\n{'='*60}")
    print("4. C-F (FLUORINATED) CHEMICALS DIFF")

    # Flag CF chemicals in both versions
    old_meta["is_cf"] = old_meta["smiles"].fillna("").apply(is_cf_smiles)
    new_meta["is_cf"] = new_meta["smiles"].fillna("").apply(is_cf_smiles)

    n_cf_old = old_meta["is_cf"].sum()
    n_cf_new = new_meta["is_cf"].sum()
    print(f"  C-F chemicals in v4.2 cache : {n_cf_old:,}")
    print(f"  C-F chemicals in v4.3 DB    : {n_cf_new:,}")

    cf_dtxsids_old = set(old_meta.loc[old_meta["is_cf"], "dsstox_substance_id"].dropna())
    cf_dtxsids_new = set(new_meta.loc[new_meta["is_cf"], "dsstox_substance_id"].dropna())
    cf_only_old    = cf_dtxsids_old - cf_dtxsids_new
    cf_only_new    = cf_dtxsids_new - cf_dtxsids_old
    cf_shared      = cf_dtxsids_old & cf_dtxsids_new
    print(f"  CF shared between versions  : {len(cf_shared):,}")
    print(f"  CF only in v4.2             : {len(cf_only_old):,}")
    print(f"  CF only in v4.3             : {len(cf_only_new):,}")

    if cf_only_old:
        lost = old_meta.loc[old_meta["dsstox_substance_id"].isin(cf_only_old),
                             ["dsstox_substance_id", "casn", "chnm", "smiles"]]
        print(f"\n  CF chemicals LOST from v4.2→v4.3 (first 20):")
        print(lost.head(20).to_string(index=False))

    if cf_only_new:
        gained = new_meta.loc[new_meta["dsstox_substance_id"].isin(cf_only_new),
                               ["dsstox_substance_id", "casn", "chnm", "smiles"]]
        print(f"\n  CF chemicals GAINED in v4.3 (first 20):")
        print(gained.head(20).to_string(index=False))

    # Label diff for CF chemicals
    cf_common_idx = ol.index.intersection(nl.index).intersection(list(cf_shared))
    print(f"\n  Label diff for shared CF chemicals ({len(cf_common_idx):,} with data):")
    ol_cf = ol.loc[cf_common_idx]
    nl_cf = nl.loc[cf_common_idx]
    for col in LABEL_COLS:
        if col not in ol_cf.columns or col not in nl_cf.columns:
            continue
        a = ol_cf[col].fillna(-1)
        b = nl_cf[col].fillna(-1)
        n_diff = (a != b).sum()
        n_old_tested = ol_cf[col].notna().sum()
        n_new_tested = nl_cf[col].notna().sum()
        print(f"    {col:<30s}  tested_v42={n_old_tested:4d}  tested_v43={n_new_tested:4d}  changed={n_diff:3d}")

        # Show individual changed rows
        if n_diff > 0:
            mask   = a != b
            idx    = cf_common_idx[mask]
            chnms  = new_meta.set_index("dsstox_substance_id").loc[idx, "chnm"].values
            a_vals = a[mask].values
            b_vals = b[mask].values
            for dtxsid, chnm, va, vb in zip(idx, chnms, a_vals, b_vals):
                tag_v = lambda v: "active" if v == 1 else ("inactive" if v == 0 else "untested")
                print(f"      {dtxsid}  {chnm[:40]:<40s}  {tag_v(va)} → {tag_v(vb)}")

    # ── 5. Summary ───────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"  Chemicals: {len(old_dtxsids):,} (v4.2) → {len(new_dtxsids):,} (v4.3)  [{len(only_new):+,} new, {-len(only_old):+,} removed]")
    print(f"  CF chemicals: {n_cf_old:,} (v4.2) → {n_cf_new:,} (v4.3)  [{len(cf_only_new):+,} new, {-len(cf_only_old):+,} removed]")
    all_old_tested = ol_aligned.notna().sum().sum()
    all_new_tested = nl_aligned.notna().sum().sum()
    print(f"  Label coverage: {all_old_tested:,} (v4.2) → {all_new_tested:,} (v4.3)")
    print(f"  SMILES changes in shared chemicals: {len(smiles_changed):,}")
    print("="*60)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compare ToxCast v4.2 cache with v4.3 database")
    parser.add_argument("--user", "-u", default=None, help="DB username")
    parser.add_argument("--password", "-p", default=None, help="DB password")
    args = parser.parse_args()
    _user     = args.user     or input("DB username: ")
    _password = args.password or getpass(f"DB password for {_user}: ")
    compare(_user, _password)
