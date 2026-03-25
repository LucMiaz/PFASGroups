#!/usr/bin/env python3
"""
Re-classify molecules from an existing pfasgroups_clinventory JSON using the
current version of PFASGroups (useful after updating group definitions without
needing the database).

Usage (from benchmark/):
    conda run -n chem python scripts/classify/reclassify_from_json.py \
        --input data/pfasgroups_clinventory_20260313T231416.json
"""

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT = BENCHMARK_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

from PFASGroups import parse_mol
from PFASGroups.core import rdkit_disable_log
rdkit_disable_log()


def classify_molecule(mol_id, smiles, formula):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "id": mol_id, "smiles": smiles, "formula": formula,
            "time_ms": None, "n_heavy_atoms": None, "n_fluorine": None,
            "n_halogens": None, "classified": False, "has_fluorine_group": False,
            "n_groups": 0, "group_ids": [], "group_names": [],
            "group_categories": [], "error": "invalid_smiles",
        }

    n_heavy   = mol.GetNumHeavyAtoms()
    n_f       = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
    n_halogen = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9, 17, 35, 53))

    t0 = time.perf_counter()
    try:
        result = parse_mol(mol, include_PFAS_definitions=False)
        elapsed_ms = (time.perf_counter() - t0) * 1000.0

        if result is None:
            return {
                "id": mol_id, "smiles": smiles, "formula": formula,
                "time_ms": elapsed_ms, "n_heavy_atoms": n_heavy,
                "n_fluorine": n_f, "n_halogens": n_halogen,
                "classified": False, "has_fluorine_group": False,
                "n_groups": 0, "group_ids": [], "group_names": [],
                "group_categories": [], "error": "no_result",
            }

        all_matches = result.get("matches", []) if isinstance(result, dict) else []
        groups = [m for m in all_matches if m.get("type") == "HalogenGroup"]
        seen_ids, unique_groups = set(), []
        for g in groups:
            gid = g.get("id")
            if gid not in seen_ids:
                seen_ids.add(gid)
                unique_groups.append(g)

        group_ids   = [g.get("id")         for g in unique_groups]
        group_names = [g.get("group_name") for g in unique_groups]
        group_cats  = [g.get("category")   for g in unique_groups]
        classified  = len(unique_groups) > 0
        has_f_group = classified and n_f > 0

        return {
            "id": mol_id, "smiles": smiles, "formula": formula,
            "time_ms": elapsed_ms, "n_heavy_atoms": n_heavy,
            "n_fluorine": n_f, "n_halogens": n_halogen,
            "classified": classified, "has_fluorine_group": has_f_group,
            "n_groups": len(unique_groups), "group_ids": group_ids,
            "group_names": group_names, "group_categories": group_cats,
            "error": None,
        }
    except Exception as exc:
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        return {
            "id": mol_id, "smiles": smiles, "formula": formula,
            "time_ms": elapsed_ms, "n_heavy_atoms": n_heavy,
            "n_fluorine": n_f, "n_halogens": n_halogen,
            "classified": False, "has_fluorine_group": False,
            "n_groups": 0, "group_ids": [], "group_names": [],
            "group_categories": [], "error": str(exc),
        }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    src = json.loads(args.input.read_text())
    molecules_in = src["molecules"]
    print(f"Re-classifying {len(molecules_in)} molecules ...")

    results = []
    for i, m in enumerate(molecules_in, 1):
        rec = classify_molecule(m["id"], m["smiles"], m.get("formula"))
        results.append(rec)
        if i % 1000 == 0:
            print(f"  {i}/{len(molecules_in)}")

    n_classified = sum(1 for r in results if r["classified"])
    n_has_f      = sum(1 for r in results if r["has_fluorine_group"])

    timestamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    out = args.output or (BENCHMARK_DIR / "data" / f"pfasgroups_clinventory_{timestamp}.json")
    out.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "tool": "PFASGroups (reclassified)",
        "source_file": str(args.input),
        "timestamp": timestamp,
        "n_molecules_fetched": len(molecules_in),
        "n_valid": len(molecules_in),
        "n_classified": n_classified,
        "n_with_fluorine_group": n_has_f,
        "molecules": results,
    }
    out.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {out}")
    print(f"Classified: {n_classified}/{len(molecules_in)}  has_fluorine_group: {n_has_f}")


if __name__ == "__main__":
    main()
