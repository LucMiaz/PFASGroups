#!/usr/bin/env python3
"""
Enrich timing JSON with atlas_class2.

The timing JSON produced by compare_oecd_clinventory_timing.py should contain
both atlas_class1 and atlas_class2, but older runs may only have atlas_class1.
This script adds atlas_class2 (and fixes atlas_class1 if needed) in-place.

Usage (from benchmark/ directory, pfasatlas conda env):
    conda run -n pfasatlas python scripts/enrich_atlas_class2.py
    conda run -n pfasatlas python scripts/enrich_atlas_class2.py --input data/oecd_clinventory_timing_20260314_224739.json
    conda run -n pfasatlas python scripts/enrich_atlas_class2.py --no-cache   # force re-classify
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT     = BENCHMARK_DIR.parent

ATLAS_DIR = Path("/home/luc/git/PFAS-atlas")
for _p in (str(ATLAS_DIR), str(ATLAS_DIR / "classification_helper")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError as e:
    sys.exit(f"ERROR: rdkit not available: {e}")

try:
    from classification_helper.classify_pfas import classify_pfas_molecule as _cp
except ImportError:
    try:
        from classify_pfas import classify_pfas_molecule as _cp  # type: ignore
    except ImportError as e:
        sys.exit(
            f"ERROR: PFAS-Atlas not available: {e}\n"
            f"Run this script in the 'pfasatlas' conda environment:\n"
            f"  conda run -n pfasatlas python scripts/enrich_atlas_class2.py"
        )


def _atlas_classify(smiles: str):
    """Return (class1, class2) for a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "error", "error"
        canon = Chem.MolToSmiles(mol)
        result = _cp(canon)
        class1 = result[0] if result else "Unknown"
        class2 = result[1] if result and len(result) > 1 else class1
        return class1, class2
    except Exception as exc:
        return "error", str(exc)[:120]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Add atlas_class2 to timing JSON (in-place)"
    )
    p.add_argument("--input", type=Path, default=None,
                   help="Timing JSON to enrich (auto-detects latest if omitted)")
    p.add_argument("--no-cache", action="store_true",
                   help="Re-classify even if atlas_class2 is already present")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ── Locate timing JSON ────────────────────────────────────────────────
    input_path = args.input
    if input_path is None:
        candidates = sorted(
            (BENCHMARK_DIR / "data").glob("oecd_clinventory_timing_*.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if not candidates:
            sys.exit("ERROR: no oecd_clinventory_timing_*.json found in data/")
        input_path = candidates[0]
    print(f"Input: {input_path}")

    # ── Load data ─────────────────────────────────────────────────────────
    with input_path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    molecules = data.get("molecules", [])
    print(f"Loaded {len(molecules):,} molecules")

    # ── Check if enrichment is needed ─────────────────────────────────────
    already_done = sum(1 for m in molecules if m.get("atlas_class2") is not None)
    if already_done == len(molecules) and not args.no_cache:
        print(f"All {len(molecules):,} molecules already have atlas_class2 — nothing to do.")
        return

    n_todo = len(molecules) - already_done if not args.no_cache else len(molecules)
    print(f"Need to add atlas_class2 for {n_todo:,} molecules …")

    # ── Classify ──────────────────────────────────────────────────────────
    t0 = time.perf_counter()
    report_every = max(1, len(molecules) // 20)
    n_errors = 0

    for i, mol in enumerate(molecules, 1):
        if i % report_every == 0 or i == len(molecules):
            elapsed = time.perf_counter() - t0
            rate = i / elapsed if elapsed > 0 else 0
            eta = (len(molecules) - i) / rate if rate > 0 else 0
            print(f"    {i:,}/{len(molecules):,}  ({rate:.0f} mol/s, ETA {eta:.0f}s) …", end="\r")

        if not args.no_cache and mol.get("atlas_class2") is not None:
            continue  # already enriched

        smiles = mol.get("smiles", "")
        if not smiles or mol.get("atlas_class1") in ("error", "disabled", None):
            mol["atlas_class2"] = mol.get("atlas_class1")
            continue

        c1, c2 = _atlas_classify(smiles)
        if c1 == "error":
            n_errors += 1
        mol["atlas_class1"] = c1
        mol["atlas_class2"] = c2

    print()
    elapsed_total = time.perf_counter() - t0
    print(f"Done in {elapsed_total:.1f}s  ({len(molecules) / elapsed_total:.0f} mol/s)")
    if n_errors:
        print(f"  Errors: {n_errors:,}")

    # ── Save in-place ─────────────────────────────────────────────────────
    with input_path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    print(f"Saved: {input_path}")


if __name__ == "__main__":
    main()
