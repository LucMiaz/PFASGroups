"""
gui/utils/fingerprints.py
─────────────────────────
Compute fingerprint matrices for ML modelling benchmarks.

Fingerprint sets:
  • PFASGroups presets (best / best_2 / … / binary / count / max_component)
  • Morgan (RDKit, radius 2, 512 bits)
  • ToxPrint (pyCSRML ToxPrintFingerprinter, 729 bits) — optional
  • TxP_PFAS (pyCSRML PFASFingerprinter, 129 bits) — optional
  • Custom TSV (user-uploaded, rows = molecules, cols = bits)
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import pandas as pd


def get_pfasgroups_fingerprints(
    embedding_set, preset: str
) -> np.ndarray:
    """Return (n_mols, n_features) array for the given FINGERPRINT_PRESETS key."""
    arr = embedding_set.to_array(preset=preset, progress=False)
    X = np.asarray(arr, dtype=float)
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


def get_morgan_fingerprints(smiles_list: list[str], radius: int = 2,
                             n_bits: int = 512) -> np.ndarray:
    """Compute Morgan fingerprints via RDKit."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi or "")
        if mol is None:
            rows.append(np.zeros(n_bits, dtype=int))
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            rows.append(np.array(fp))
    return np.vstack(rows).astype(float)


def get_toxprint_fingerprints(smiles_list: list[str]) -> Optional[np.ndarray]:
    """Compute ToxPrint 729-bit fingerprints (requires pyCSRML)."""
    try:
        from pyCSRML import ToxPrintFingerprinter
    except ImportError:
        return None
    fp = ToxPrintFingerprinter()
    rows = [fp.fingerprint(smi) for smi in smiles_list]
    return np.vstack(rows).astype(float)


def get_txp_pfas_fingerprints(smiles_list: list[str]) -> Optional[np.ndarray]:
    """Compute TxP_PFAS 129-bit fingerprints (requires pyCSRML)."""
    try:
        from pyCSRML import PFASFingerprinter
    except ImportError:
        return None
    fp = PFASFingerprinter()
    rows = [fp.fingerprint(smi) for smi in smiles_list]
    return np.vstack(rows).astype(float)


def load_custom_fingerprints(path: str,
                              smiles_list: list[str]) -> Optional[np.ndarray]:
    """Load a TSV/CSV where rows correspond to molecules (same order as smiles_list)."""
    try:
        sep = "\t" if path.endswith(".tsv") else ","
        df = pd.read_csv(path, sep=sep, header=None)
        X = df.values.astype(float)
        if X.shape[0] != len(smiles_list):
            raise ValueError(
                f"Custom fingerprint file has {X.shape[0]} rows "
                f"but {len(smiles_list)} molecules."
            )
        return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    except Exception:
        return None


def is_pycsrml_available() -> bool:
    try:
        import pyCSRML  # noqa: F401
        return True
    except ImportError:
        return False
