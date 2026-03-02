#!/usr/bin/env python3
"""
extract_and_fix_unparsed_sdf.py
================================
Finds PFASSTRUCTV5 SDF entries that produced empty SMILES in the unmatched CSV
(i.e. RDKit could not parse them), extracts the raw SDF blocks, and attempts to
recover SMILES via PubChem REST API using CASRN or preferred name.

These 8 entries have no embedded structure (SMILES = N/A, MOL_FILE = MOL_FILE)
in the original SDF, so RDKit cannot parse them. We fall back to PubChem.

Output
------
  tests/results/pfasstructv5_unparsed.sdf     — raw SDF blocks for these entries
  tests/results/pfasstructv5_repair_report.csv — per-DTXSID recovery result
"""

from __future__ import annotations

import csv
import re
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
import json
from pathlib import Path

# Add PFASGroups to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from rdkit import Chem
from PFASGroups import PFASDefinition
from PFASGroups.core import PFAS_DEFINITIONS_FILE

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE        = Path(__file__).parent
CSV_PATH    = HERE / "results" / "pfasstructv5_unmatched.csv"
SDF_PATH    = HERE / "test_data" / "PFASSTRUCTv5_2022_from_DSSTox.sdf"
OUT_SDF     = HERE / "results" / "pfasstructv5_unparsed.sdf"
OUT_REPORT  = HERE / "results" / "pfasstructv5_repair_report.csv"

PUBCHEM_BASE  = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_DELAY = 0.3   # seconds between PubChem requests (rate limit courtesy)
PFAS_DEFINITION_ID = 5

# Load PFAS definition once
with open(PFAS_DEFINITIONS_FILE) as _fh:
    _defs = json.load(_fh)
_entry = next(d for d in _defs if d["id"] == PFAS_DEFINITION_ID)
PFAS_DEF = PFASDefinition(**_entry)

# ---------------------------------------------------------------------------
# Load DTXSIDs that had empty SMILES
# ---------------------------------------------------------------------------
empty_dtxsids: set[str] = set()
with open(CSV_PATH, newline="") as fh:
    for row in csv.DictReader(fh):
        if not row["SMILES"].strip():
            empty_dtxsids.add(row["DTXSID"].strip())

print(f"Target DTXSIDs with empty SMILES: {sorted(empty_dtxsids)}")

# ---------------------------------------------------------------------------
# Iterate SDF entries and extract metadata + raw block
# ---------------------------------------------------------------------------
def iter_entries(sdf_path: Path):
    """Yield metadata dicts (+ raw block text) for every SDF entry."""
    current_lines: list[str] = []
    with open(sdf_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() == "$$$$":
                current_lines.append(line)
                block = "".join(current_lines)

                def prop(tag: str, _block: str = block) -> str:
                    m = re.search(rf">\s*<{re.escape(tag)}>\s*\n(.+)", _block)
                    return m.group(1).strip() if m else ""

                yield {
                    "dtxsid":     prop("DTXSID") or current_lines[0].strip(),
                    "casrn":      prop("CASRN"),
                    "name":       prop("PREFERRED_NAME"),
                    "smiles_sdf": prop("SMILES"),
                    "inchikey":   prop("INCHIKEY"),
                    "formula":    prop("MOLECULAR_FORMULA"),
                    "raw_block":  block,
                }
                current_lines = []
            else:
                current_lines.append(line)

# ---------------------------------------------------------------------------
# PubChem helpers
# ---------------------------------------------------------------------------
def _pubchem_get(url: str) -> dict | None:
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    except Exception:
        return None


def pubchem_smiles_from_cas(casrn: str) -> str:
    if not casrn or casrn in ("N/A", ""):
        return ""
    url = (f"{PUBCHEM_BASE}/compound/name/"
           f"{urllib.parse.quote(casrn)}/property/IsomericSMILES/JSON")
    data = _pubchem_get(url)
    if data:
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("SMILES", "") or props[0].get("IsomericSMILES", "")
    return ""


def pubchem_smiles_from_inchikey(inchikey: str) -> str:
    if not inchikey or inchikey in ("N/A", ""):
        return ""
    url = (f"{PUBCHEM_BASE}/compound/inchikey/"
           f"{urllib.parse.quote(inchikey)}/property/IsomericSMILES/JSON")
    data = _pubchem_get(url)
    if data:
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("SMILES", "") or props[0].get("IsomericSMILES", "")
    return ""


def pubchem_smiles_from_name(name: str) -> str:
    if not name or name in ("N/A", ""):
        return ""
    url = (f"{PUBCHEM_BASE}/compound/name/"
           f"{urllib.parse.quote(name)}/property/IsomericSMILES/JSON")
    data = _pubchem_get(url)
    if data:
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("SMILES", "") or props[0].get("IsomericSMILES", "")
    return ""


# ---------------------------------------------------------------------------
# Main pass
# ---------------------------------------------------------------------------
results: list[dict] = []
raw_blocks: list[str] = []

print(f"\nScanning {SDF_PATH.name} for target DTXSIDs ...")
for entry in iter_entries(SDF_PATH):
    dtxsid = entry["dtxsid"]
    if dtxsid not in empty_dtxsids:
        continue

    print(f"\n  {dtxsid}  ({entry['name']})")
    print(f"    CASRN={entry['casrn']}  InChIKey={entry['inchikey']}")
    print(f"    SDF SMILES field: {entry['smiles_sdf']}")

    raw_blocks.append(entry["raw_block"])

    strategy, smiles = "", ""
    for label, fn, val in [
        ("pubchem_cas",      pubchem_smiles_from_cas,     entry["casrn"]),
        ("pubchem_inchikey", pubchem_smiles_from_inchikey, entry["inchikey"]),
        ("pubchem_name",     pubchem_smiles_from_name,     entry["name"]),
    ]:
        time.sleep(REQUEST_DELAY)
        smi = fn(val)
        if smi:
            strategy, smiles = label, smi
            break

    fixed = bool(smiles)

    # Compute F ratio and PFASGroups match for recovered SMILES
    f_ratio_with_h = f_ratio_without_h = matched = ""
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol_h = Chem.AddHs(mol)
            nF = sum(1 for a in mol_h.GetAtoms() if a.GetAtomicNum() == 9)
            f_ratio_with_h    = round(nF / mol_h.GetNumAtoms(), 6) if mol_h.GetNumAtoms() else 0.0
            heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
            nF_h  = sum(1 for a in heavy if a.GetAtomicNum() == 9)
            f_ratio_without_h = round(nF_h / len(heavy), 6) if heavy else 0.0
            matched = PFAS_DEF.applies_to_molecule(mol)

    print(f"    -> strategy={strategy!r}  smiles={smiles[:80]}  fixed={fixed}")
    if smiles:
        print(f"       F/heavy={f_ratio_without_h}  matched_by_pfasgroups={matched}")
    results.append({
        "DTXSID":                dtxsid,
        "CASRN":                 entry["casrn"],
        "preferred_name":        entry["name"],
        "formula":               entry["formula"],
        "strategy":              strategy or "no_structure_in_sdf",
        "recovered_smiles":      smiles,
        "fixed":                 fixed,
        "F_ratio_with_H":        f_ratio_with_h,
        "F_ratio_without_H":     f_ratio_without_h,
        "matched_by_pfasgroups": matched,
    })

missing = empty_dtxsids - {r["DTXSID"] for r in results}
if missing:
    print(f"\n  WARNING: DTXSIDs not found in SDF: {missing}")

# ---------------------------------------------------------------------------
# Write raw blocks to output SDF
# ---------------------------------------------------------------------------
OUT_SDF.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_SDF, "w", encoding="utf-8") as fh:
    fh.writelines(raw_blocks)
print(f"\nRaw SDF blocks written -> {OUT_SDF}  ({len(raw_blocks)} entries)")

# ---------------------------------------------------------------------------
# Write repair report
# ---------------------------------------------------------------------------
fieldnames = ["DTXSID", "CASRN", "preferred_name", "formula",
              "strategy", "recovered_smiles", "fixed",
              "F_ratio_with_H", "F_ratio_without_H", "matched_by_pfasgroups"]
with open(OUT_REPORT, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)
print(f"Repair report written -> {OUT_REPORT}")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
n_fixed  = sum(1 for r in results if r["fixed"])
n_failed = len(results) - n_fixed
print(f"\nSummary: {n_fixed}/{len(results)} recovered via PubChem, {n_failed} still unresolvable.")
