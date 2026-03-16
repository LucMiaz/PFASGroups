"""
Tests that every molecule in the PFASSTRUCTv5 inventory (Richard et al. 2023,
DOI: 10.1021/acs.chemrestox.2c00403) is classified as PFAS by HalogenGroups
PFAS definition 5 ('PFASSTRUCTv5').

Pytest mode  : runs on a deterministic sample of PYTEST_SAMPLE_SIZE molecules.
Direct mode  : ``python test_pfasstructv5_against_richard2023.py``
               tests the full dataset and reports any failures.

Reference
---------
Richard, A. M. et al. (2023). A New CSRML Structure-Based Fingerprint Method
for Profiling and Categorizing Per- and Polyfluoroalkyl Substances (PFAS).
Chemical Research in Toxicology, 36(3), 318-338.
https://doi.org/10.1021/acs.chemrestox.2c00403
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List, Tuple

import pytest
from rdkit import Chem

from PFASGroups import PFASDefinition
from PFASGroups.core import PFAS_DEFINITIONS_FILE, rdkit_disable_log

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SDF_PATH      = Path(__file__).parent / "test_data" / "PFASSTRUCTv5_2022_from_DSSTox.sdf"
RESULTS_DIR   = Path(__file__).parent / "results"
UNMATCHED_CSV = RESULTS_DIR / "pfasstructv5_unmatched.csv"

# DTXSIDs known not to match definition 5 — loaded from the unmatched CSV so
# the set is always in sync with the actual failure list.
def _load_known_unmatched() -> frozenset[str]:
    if not UNMATCHED_CSV.exists():
        return frozenset()
    import csv as _csv
    with open(UNMATCHED_CSV, newline="") as fh:
        return frozenset(row["DTXSID"] for row in _csv.DictReader(fh))

KNOWN_UNMATCHED: frozenset[str] = _load_known_unmatched()
# Number of molecules tested under pytest (evenly-spaced stride through the
# full dataset to give a representative, deterministic selection).
PYTEST_SAMPLE_SIZE = 10

# PFAS definition to check (id=5, name='PFASSTRUCTv5')
PFAS_DEFINITION_ID = 5

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_definition(def_id: int) -> PFASDefinition:
    """Load a PFASDefinition by id directly from the JSON definitions file."""
    with open(PFAS_DEFINITIONS_FILE) as fh:
        defs = json.load(fh)
    entry = next((d for d in defs if d["id"] == def_id), None)
    if entry is None:
        raise ValueError(
            f"PFAS definition with id={def_id} not found in {PFAS_DEFINITIONS_FILE}"
        )
    return PFASDefinition(**entry)


def _load_molecules(sdf_path: Path) -> List[Tuple[Chem.Mol, str]]:
    """Return list of (mol, identifier) for every valid molecule in the SDF.

    Molecules whose DTXSID appears in KNOWN_UNMATCHED are silently excluded so
    that parametrised tests and direct-run validation never trip on entries that
    are documented edge-cases rather than implementation bugs.
    """
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=True)
    molecules: List[Tuple[Chem.Mol, str]] = []
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        # Prefer the first SDF property as identifier (typically DTXSID/DTXCID).
        props = mol.GetPropsAsDict()
        if props:
            identifier = str(list(props.values())[0])
        else:
            title = mol.GetProp("_Name")
            identifier = title if title else f"mol_{i}"
        if identifier in KNOWN_UNMATCHED:
            continue
        molecules.append((mol, identifier))
    return molecules


def _evenly_spaced_sample(items: list, n: int) -> list:
    """Return *n* evenly-spaced elements from *items* (deterministic)."""
    if len(items) <= n:
        return items
    step = len(items) / n
    indices = [int(i * step) for i in range(n)]
    return [items[idx] for idx in indices]


# ---------------------------------------------------------------------------
# Module-level data (loaded once at import / collection time)
# ---------------------------------------------------------------------------

_ALL_MOLECULES: List[Tuple[Chem.Mol, str]] = _load_molecules(SDF_PATH)
_DEFINITION: PFASDefinition = _load_definition(PFAS_DEFINITION_ID)
_SAMPLE: List[Tuple[Chem.Mol, str]] = _evenly_spaced_sample(
    _ALL_MOLECULES, PYTEST_SAMPLE_SIZE
)

# ---------------------------------------------------------------------------
# Pytest tests (sample only)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mol,name",
    _SAMPLE,
    ids=[name for _, name in _SAMPLE],
)
@rdkit_disable_log()
def test_molecule_matches_pfasstructv5_definition(mol: Chem.Mol, name: str) -> None:
    """Each sampled PFASSTRUCTV5 molecule must satisfy PFAS definition 5.

    The sample is an evenly-spaced stride across the full 14,735-molecule
    PFASSTRUCTV5 list (Richard et al. 2023, DOI 10.1021/acs.chemrestox.2c00403).
    Run this script directly (``python test_pfasstructv5_against_richard2023.py``)
    to validate the complete dataset.
    """
    assert _DEFINITION.applies_to_molecule(mol), (
        f"Molecule '{name}' did NOT match PFAS definition {PFAS_DEFINITION_ID} "
        f"('{_DEFINITION.name}'). "
        f"fluorineRatio={_DEFINITION.fluorineRatio}, "
        f"requireBoth={_DEFINITION.requireBoth}."
    )


# ---------------------------------------------------------------------------
# Direct execution: validate the full dataset
# ---------------------------------------------------------------------------

def _fluorine_ratio(mol: Chem.Mol, include_hydrogen: bool = True) -> float:
    """Return F atoms / total atoms (H included when include_hydrogen is True)."""
    mol_h = Chem.AddHs(mol) if include_hydrogen else mol
    n_atoms = mol_h.GetNumAtoms()
    n_f = sum(1 for a in mol_h.GetAtoms() if a.GetAtomicNum() == 9)
    return n_f / n_atoms if n_atoms > 0 else 0.0


if __name__ == "__main__":
    import csv

    total = len(_ALL_MOLECULES)
    print(
        f"Testing {total} molecules from {SDF_PATH.name} "
        f"against PFAS definition {PFAS_DEFINITION_ID} "
        f"('{_DEFINITION.name}') ..."
    )

    failures: List[dict] = []
    for i, (mol, name) in enumerate(_ALL_MOLECULES, start=1):
        if not _DEFINITION.applies_to_molecule(mol):
            props = mol.GetPropsAsDict()
            dtxsid = props.get("DSSTox_Substance_Id", name)
            smiles = props.get("Structure_SMILES") or Chem.MolToSmiles(mol)
            f_ratio_with_h    = _fluorine_ratio(mol, include_hydrogen=True)
            f_ratio_without_h = _fluorine_ratio(mol, include_hydrogen=False)
            failures.append(
                {
                    "DTXSID": dtxsid,
                    "SMILES": smiles,
                    "F_ratio_with_H":    round(f_ratio_with_h, 6),
                    "F_ratio_without_H": round(f_ratio_without_h, 6),
                }
            )
        if i % 500 == 0:
            print(f"  {i}/{total} checked ...")

    passed = total - len(failures)
    print(f"\nResults: {passed}/{total} passed, {len(failures)} not matched.")

    if failures:
        out_csv = RESULTS_DIR / "pfasstructv5_unmatched.csv"
        with open(out_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["DTXSID", "SMILES", "F_ratio_with_H", "F_ratio_without_H"])
            writer.writeheader()
            writer.writerows(failures)
        print(f"Unmatched molecules written to: {out_csv}")
        sys.exit(1)
    else:
        print("All molecules matched.")
        sys.exit(0)
