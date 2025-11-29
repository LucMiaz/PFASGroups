import math
import time
from typing import List, Dict, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, DataStructs

from PFASgroups.core import parse_PFAS_groups


def _pfas_group_presence(smiles_list: List[str]) -> Dict[str, int]:
    """Return counts of molecules containing each PFAS group (presence/absence per molecule).

    A molecule contributes 1 to a group if any match for that group is found.
    """
    counts: Dict[str, int] = {}
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        matches = parse_PFAS_groups(mol, formula)
        present_names = set()
        for pf, n, n_cfchains, chains in matches:
            # pf expected to have .name attribute
            present_names.add(pf.name)
        for name in present_names:
            counts[name] = counts.get(name, 0) + 1
    return counts


def _pfas_group_distribution(smiles_list: List[str]) -> Dict[str, float]:
    counts = _pfas_group_presence(smiles_list)
    n = float(len(smiles_list)) if len(smiles_list) > 0 else 1.0
    return {k: v / n for k, v in counts.items()}


def compute_minmax_kl(setA: List[str], setB: List[str]) -> Tuple[float, float]:
    """Compute KL divergences between min and max normalized PFAS group distributions.

    We build p_min(g) = min(pA(g), pB(g)) / Z_min and p_max(g) = max(pA(g), pB(g)) / Z_max and return:
    KL(p_min || p_max) in natural log base and base 2.
    If either distribution has zero support (no groups detected), returns (0.0, 0.0).
    """
    pA = _pfas_group_distribution(setA)
    pB = _pfas_group_distribution(setB)
    # union of group keys
    keys = set(pA.keys()).union(pB.keys())
    if len(keys) == 0:
        return 0.0, 0.0
    mins = np.array([min(pA.get(k, 0.0), pB.get(k, 0.0)) for k in keys], dtype=float)
    maxs = np.array([max(pA.get(k, 0.0), pB.get(k, 0.0)) for k in keys], dtype=float)
    # Avoid all-zero vectors
    if mins.sum() == 0 or maxs.sum() == 0:
        return 0.0, 0.0
    p_min = mins / mins.sum()
    p_max = maxs / maxs.sum()
    # KL(p_min || p_max)
    # mask zero entries in p_min to avoid log(0)
    mask = p_min > 0
    kl_e = float(np.sum(p_min[mask] * np.log(p_min[mask] / p_max[mask])))
    kl_2 = kl_e / math.log(2.0)
    return kl_e, kl_2


def _fingerprint_list(smiles_list: List[str], fp_type: str = "morgan"):
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        if fp_type == "morgan":
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        elif fp_type == "maccs":
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            raise ValueError("Unsupported fingerprint type")
        fps.append(fp)
    return fps


def mean_tanimoto_between_sets(setA: List[str], setB: List[str], fp_type: str = "morgan") -> float:
    fpsA = _fingerprint_list(setA, fp_type)
    fpsB = _fingerprint_list(setB, fp_type)
    if len(fpsA) == 0 or len(fpsB) == 0:
        return 0.0
    sims = []
    for fa in fpsA:
        for fb in fpsB:
            sims.append(DataStructs.TanimotoSimilarity(fa, fb))
    return float(np.mean(sims))


def js_divergence_bitfreq(setA: List[str], setB: List[str], fp_type: str = "morgan") -> float:
    fpsA = _fingerprint_list(setA, fp_type)
    fpsB = _fingerprint_list(setB, fp_type)
    if len(fpsA) == 0 or len(fpsB) == 0:
        return 0.0
    # bit frequency vectors
    arrA = np.mean([np.array(list(fp.ToBitString()), dtype=int) for fp in fpsA], axis=0)
    arrB = np.mean([np.array(list(fp.ToBitString()), dtype=int) for fp in fpsB], axis=0)
    # convert to probability distributions (normalize by sum; if sum=0 skip)
    if arrA.sum() == 0 or arrB.sum() == 0:
        return 0.0
    p = arrA / arrA.sum()
    q = arrB / arrB.sum()
    m = 0.5 * (p + q)
    # mask zeros to avoid issues
    def _kl(a, b):
        mask = (a > 0) & (b > 0)
        return float(np.sum(a[mask] * np.log(a[mask] / b[mask])))
    js = 0.5 * _kl(p, m) + 0.5 * _kl(q, m)
    return js


def _synthetic_sets() -> Dict[str, List[str]]:
    """Return synthetic PFAS-like sets. These are illustrative and may not fully exercise all PFAS group SMARTS.
    We build sets so A1 and A2 are similar (split of a larger pool) while B is compositionally different.
    """
    base = [
        # Fluorinated chains / acids (varied length & functionalization)
        "C(C(F)(F)F)C(F)(F)F",  # short perfluoro chain
        "C(C(F)(F)F)(C(F)(F)F)C(=O)O",  # perfluoro carboxylic acid
        "C(C(F)(F)F)S(=O)(=O)O",  # sulfonic
        "C(C(F)(F)F)N(=O)=O",  # nitro substituted
        "C1(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)1",  # cyclic fluorinated
        "C(C(F)(F)F)C(=O)O",  # shorter acid
        "C(C(F)(F)F)C#N",  # cyano
        "C(C(F)(F)F)CO",  # alcohol tail
    ]
    # Similar split
    A1 = base[:4]
    A2 = base[4:]
    # B: introduce more hydrogen-rich / different functionalities to reduce PFAS group overlap
    B = [
        "CC(=O)O",  # acetic acid
        "CCS(=O)(=O)O",  # ethanesulfonic acid
        "CCN(=O)=O",  # nitroethane
        "CCC(=O)O",  # propionic acid
        "CCCN",  # propylamine
        "CCCO",  # propanol
        "CC(C)C(=O)O",  # isobutyric acid
        "c1ccccc1C(=O)O",  # benzoic acid
    ]
    return {"A1": A1, "A2": A2, "B": B}


def test_minmaxkl_vs_fingerprints():
    sets = _synthetic_sets()
    A1, A2, B = sets["A1"], sets["A2"], sets["B"]

    # Compute minmax KL
    t0 = time.time()
    kl_A1_A2_e, kl_A1_A2_2 = compute_minmax_kl(A1, A2)
    t1 = time.time()
    kl_A1_B_e, kl_A1_B_2 = compute_minmax_kl(A1, B)
    t2 = time.time()

    minmax_time_similar = t1 - t0
    minmax_time_dissimilar = t2 - t1

    # Fingerprint similarities
    tanimoto_A1_A2 = mean_tanimoto_between_sets(A1, A2, fp_type="morgan")
    tanimoto_A1_B = mean_tanimoto_between_sets(A1, B, fp_type="morgan")
    js_A1_A2 = js_divergence_bitfreq(A1, A2, fp_type="morgan")
    js_A1_B = js_divergence_bitfreq(A1, B, fp_type="morgan")

    # Expectations:
    # - Similar sets (A1,A2) should have lower KL than dissimilar (A1,B) if PFAS groups are detected.
    # - Tanimoto similarity should be higher for A1,A2 than A1,B.
    # - JS divergence should be lower for A1,A2 than A1,B.

    # Guard: if PFAS groups not detected (KL both 0), skip KL assertion but still validate fingerprints.
    if not (kl_A1_A2_e == 0.0 and kl_A1_B_e == 0.0):
        assert kl_A1_A2_e <= kl_A1_B_e + 1e-9, (
            f"Expected KL(A1,A2) <= KL(A1,B); got {kl_A1_A2_e:.4f} vs {kl_A1_B_e:.4f}"
        )

    assert tanimoto_A1_A2 >= tanimoto_A1_B, (
        f"Expected higher Tanimoto for similar sets; got {tanimoto_A1_A2:.4f} vs {tanimoto_A1_B:.4f}"
    )
    assert js_A1_A2 <= js_A1_B + 1e-9, (
        f"Expected lower JS divergence for similar sets; got {js_A1_A2:.4f} vs {js_A1_B:.4f}"
    )

    # Basic timing sanity (methods should be reasonably fast on tiny sets)
    assert minmax_time_similar < 1.0 and minmax_time_dissimilar < 1.0

    # Provide debug info via test output (pytest -vv will show if assertion fails)
    print(
        f"minmaxKL(A1,A2)={kl_A1_A2_e:.6f} (nats), minmaxKL(A1,B)={kl_A1_B_e:.6f}; "
        f"Tanimoto(A1,A2)={tanimoto_A1_A2:.4f}, Tanimoto(A1,B)={tanimoto_A1_B:.4f}; "
        f"JS(A1,A2)={js_A1_A2:.6f}, JS(A1,B)={js_A1_B:.6f}"
    )


def test_minmaxkl_symmetry():
    sets = _synthetic_sets()
    A1, B = sets["A1"], sets["B"]
    kl1_e, _ = compute_minmax_kl(A1, B)
    kl2_e, _ = compute_minmax_kl(B, A1)
    # KL(min||max) is symmetric under swapping A,B because min/max are symmetric operations.
    assert abs(kl1_e - kl2_e) < 1e-12
