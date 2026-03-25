#!/usr/bin/env python3
"""
Two-phase Tanimoto benchmark for PFASGroups fingerprint variants.

Phase 1 – FP Selection
  Tests all 39 fingerprint configurations on a curated structural series
  (~30 PFAS) spanning chain-length homologues, branching position/size, and
  functional-group variants.  Selects the top-5 most discriminating configs.

Phase 2 – Discrimination benchmark vs TxP-PFAS
  Runs the top-5 PFASGroups configs AND the TxP-PFAS (129-bit CSRML)
  fingerprint from Richard et al. 2023 on the same 18-compound set.
  Compares discrimination power: lower mean off-diagonal Tanimoto = more
  discriminating.  TxP-PFAS is a competitor, not a ground-truth reference.

Usage
-----
    python fingerprint_tanimoto_benchmark.py
    python fingerprint_tanimoto_benchmark.py --outdir /path/to/output --dpi 130

Outputs (all in <outdir>/)
--------------------------
    tanimoto_p1_discrimination.png   Bar chart – all FP configs sorted by mean T
    tanimoto_p1_metric_impact.png    Per-metric impact relative to binary baseline
    tanimoto_p1_heatmaps.png         Heatmap grid (selected configs)
    tanimoto_p1_mds.png              2-D MDS grid (selected configs)
    tanimoto_p1_ranking.png          Ranked summary heatmap
    tanimoto_p2_heatmaps.png         Phase-2 heatmaps: top-5 PFASGroups configs + TxP-PFAS
    tanimoto_p2_discrimination.png   Bar chart comparing discrimination (mean T) for all
    tanimoto_p2_mds.png              MDS overlay (Phase 2)
    tanimoto_summary.csv             Per-FP statistics (Phase 1)
    tanimoto_report.html             Self-contained HTML report
"""

from __future__ import annotations

import argparse
import base64
import io
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.spatial import ConvexHull

try:
    from rdkit import Chem as _Chem
    from rdkit.Chem import Draw as _RDDraw
    _HAS_RDKIT_DRAW = True
except ImportError:
    _HAS_RDKIT_DRAW = False


def _draw_group_hulls(
    ax: plt.Axes,
    coords: np.ndarray,
    series: list[str],
    colour_map: dict[str, str],
    alpha_fill: float = 0.13,
    alpha_edge: float = 0.60,
    pad: float = 0.10,
) -> None:
    """Draw a filled convex-hull region around each group of MDS points.

    Groups with <3 points (or collinear points) get a circle instead of a hull.
    ``pad`` expands the hull outward by this fraction of the centroid-vertex distance.
    """
    coord_range = float(max(np.ptp(coords[:, 0]), np.ptp(coords[:, 1]), 1e-9))
    base_r = 0.05 * coord_range  # fallback radius for degenerate groups

    groups: dict[str, list[int]] = {}
    for i, s in enumerate(series):
        groups.setdefault(s, []).append(i)

    for grp, idx in groups.items():
        col = colour_map.get(grp, '#888888')
        pts = coords[np.array(idx)]

        centroid = pts.mean(axis=0)
        drawn = False
        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hull_pts = pts[hull.vertices]
                padded = centroid + (1.0 + pad) * (hull_pts - centroid)
                poly = plt.Polygon(
                    padded, closed=True,
                    facecolor=col, edgecolor=col,
                    alpha=alpha_fill, zorder=0, linewidth=0)
                ax.add_patch(poly)
                closed_padded = np.vstack([padded, padded[:1]])
                ax.plot(closed_padded[:, 0], closed_padded[:, 1],
                        color=col, alpha=alpha_edge, linewidth=1.2, zorder=1)
                drawn = True
            except Exception:
                pass  # fall through to circle
        if not drawn:
            r = float(np.sqrt(((pts - centroid) ** 2).sum(axis=1)).max())
            r = max(r * (1.0 + pad), base_r)
            circle = plt.Circle(centroid, r, color=col, fill=True,
                                 alpha=alpha_fill, zorder=0)
            ax.add_patch(circle)
            circle_edge = plt.Circle(centroid, r, color=col, fill=False,
                                     alpha=alpha_edge, linewidth=1.2, zorder=1)
            ax.add_patch(circle_edge)


def _group_separation_score(coords: np.ndarray, series: list[str]) -> float:
    """Group-separation index = mean inter-centroid distance / mean within-group spread.

    Higher values mean the fingerprint places chemical groups further apart
    relative to their internal scatter.  Shown in each MDS panel title.
    """
    groups: dict[str, list[int]] = {}
    for i, s in enumerate(series):
        groups.setdefault(s, []).append(i)

    centroids: dict[str, np.ndarray] = {}
    spreads: list[float] = []
    for g, idx in groups.items():
        pts = coords[np.array(idx)]
        c = pts.mean(axis=0)
        centroids[g] = c
        if len(pts) > 1:
            spreads.append(float(np.sqrt(((pts - c) ** 2).sum(axis=1)).mean()))

    names = list(centroids.keys())
    if len(names) < 2:
        return 0.0

    inter_dists = [
        float(np.sqrt(((centroids[names[i]] - centroids[names[j]]) ** 2).sum()))
        for i in range(len(names))
        for j in range(i + 1, len(names))
    ]
    mean_inter = float(np.mean(inter_dists))
    mean_spread = float(np.mean(spreads)) if spreads else 1e-9
    return mean_inter / (mean_spread + 1e-9)


def _per_group_gsi(coords: np.ndarray, series: list[str]) -> dict[str, float]:
    """Per-group separability contribution.

    For each group g:
        separability(g) = mean distance from centroid(g) to all other centroids
                          / (within-group spread of g + epsilon)

    This lets you see which groups drive the overall GSI and which are poorly
    separated.  The 'FT group drives GSI' hypothesis can be directly checked
    by comparing each group's score across fingerprint configs.
    """
    groups: dict[str, list[int]] = {}
    for i, s in enumerate(series):
        groups.setdefault(s, []).append(i)

    centroids: dict[str, np.ndarray] = {}
    spreads: dict[str, float] = {}
    for g, idx in groups.items():
        pts = coords[np.array(idx)]
        c = pts.mean(axis=0)
        centroids[g] = c
        if len(pts) > 1:
            spreads[g] = float(np.sqrt(((pts - c) ** 2).sum(axis=1)).mean())
        else:
            spreads[g] = 0.0

    result: dict[str, float] = {}
    names = list(centroids.keys())
    for g in names:
        other_dists = [
            float(np.sqrt(((centroids[g] - centroids[h]) ** 2).sum()))
            for h in names if h != g
        ]
        mean_dist = float(np.mean(other_dists)) if other_dists else 0.0
        result[g] = mean_dist / (spreads[g] + 1e-9)
    return result


def _smiles_grid_b64(
    smiles: list[str],
    labels: list[str],
    n_cols: int = 6,
    sub_img_size: tuple[int, int] = (220, 160),
) -> str | None:
    """Return base64-encoded PNG of a 2-D structure grid, or None if RDKit Draw unavailable."""
    if not _HAS_RDKIT_DRAW:
        return None
    mols = [_Chem.MolFromSmiles(s) for s in smiles]
    try:
        img_bytes = _RDDraw.MolsToGridImage(
            mols,
            molsPerRow=n_cols,
            subImgSize=sub_img_size,
            legends=labels,
            returnPNG=True,
        )
        return base64.b64encode(img_bytes).decode()
    except Exception:
        return None


def _classical_mds(dist: np.ndarray, n_components: int = 2) -> np.ndarray:
    """Classical (metric) MDS via double-centring and eigendecomposition."""
    n = dist.shape[0]
    D2 = dist ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D2 @ H
    vals, vecs = np.linalg.eigh(B)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    vals_pos = np.maximum(vals[:n_components], 0.0)
    return vecs[:, :n_components] * np.sqrt(vals_pos)

# ── repo paths ────────────────────────────────────────────────────────────────
SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT     = BENCHMARK_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════════════
# Phase-1 compound library – curated structural series
# ══════════════════════════════════════════════════════════════════════════════
# Groups: PFAA, PFSA, iso-PFCA, iso-PFSA, HFPO, FTOH, FTS
P1_COMPOUNDS = [
    # ── PFAA (perfluoroalkyl carboxylic acid) chain-length series C2–C11 ────
    ("OC(=O)C(F)(F)F",                                                                                 "TFA (C2)",    "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)F",                                                                          "PFPrA (C3)",  "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)F",                                                                   "PFBA (C4)",   "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                                            "PFPeA (C5)",  "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                                     "PFHxA (C6)",  "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                             "PFHpA (C7)",  "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                      "PFOA (C8)",   "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                               "PFNA (C9)",   "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                        "PFDA (C10)",  "PFAA"),
    ("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                 "PFUnDA (C11)","PFAA"),

    # ── PFSA (perfluoroalkyl sulfonic acid) series C4 / C6 / C8 / C10 ───────
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",                                                        "PFBS (C4)",   "PFSA"),
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",                                          "PFHxS (C6)",  "PFSA"),
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",                            "PFOS (C8)",   "PFSA"),
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",              "PFDS (C10)",  "PFSA"),

    # ── iso-PFCA  (CF3 branch at alpha-C, varying chain length) ─────────────
    ("OC(=O)C(F)(C(F)(F)F)C(F)(F)F",                                                                   "iso-PFBA",    "iso-PFCA"),
    ("OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)F",                                                             "iso-PFPeA",   "iso-PFCA"),
    ("OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F",                                                      "iso-PFHxA",   "iso-PFCA"),
    ("OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                              "iso-PFHpA",   "iso-PFCA"),
    # CF3 branch at beta-C (different branch position, same total chain)
    ("OC(=O)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F",                                                             "iso-PFBA-b",  "iso-PFCA"),
    ("OC(=O)C(F)(F)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F",                                              "iso-PFHxA-b", "iso-PFCA"),

    # ── iso-PFSA  (CF3 branch, varying chain) ────────────────────────────────
    ("OS(=O)(=O)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F",                                                        "iso-PFBS",    "iso-PFSA"),
    ("OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F",                                          "iso-PFHxS",   "iso-PFSA"),

    # ── HFPO-DA (GenX) – cyclic ether PFCA ───────────────────────────────────
    ("OC(=O)C(F)(C(F)(F)F)OC(F)(F)F",                                                                  "HFPO-DA",     "HFPO"),

    # ── FTOH (fluorotelomer alcohols) series 4:2 / 6:2 / 8:2 / 10:2 ─────────
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                                               "4:2 FTOH",    "FTOH"),
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                                 "6:2 FTOH",    "FTOH"),
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                   "8:2 FTOH",    "FTOH"),
    ("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                     "10:2 FTOH",   "FTOH"),

    # ── FTS (fluorotelomer sulfonic acids) 6:2 / 8:2 ─────────────────────────
    ("OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",                                        "6:2 FTS",     "FTS"),
    ("OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",                           "8:2 FTS",     "FTS"),

    # ── FTCA (fluorotelomer carboxylic acids) – same FT scaffold as FTOH, ─────
    #    COOH instead of OH  →  tests FG discrimination (size-matched with FTOH)
    ("OC(=O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                                            "4:2 FTCA",    "FTCA"),
    ("OC(=O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                              "6:2 FTCA",    "FTCA"),
    ("OC(=O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                               "8:2 FTCA",    "FTCA"),

    # ── PFPA (perfluoroalkyl phosphonic acids) – same PF chain as PFCA/PFSA, ──
    #    phosphonate head group  →  FG variation at constant component size
    ("OP(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                                         "PFBA-PA (C4)",   "PFPA"),
    ("OP(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                                           "PFHxA-PA (C6)",  "PFPA"),

    # ── iso-heavy (complex branching: gem-CF3, large branch, double branch) ───
    #    explores branching *degree* beyond the mono-CF3 in iso-PFCA
    ("OC(=O)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F",                                                            "gem-CF3 C4",     "iso-heavy"),
    ("OC(=O)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F",                                                     "gem-CF3 C5",     "iso-heavy"),
    ("OC(=O)C(F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F",                                              "C2F5-br alpha",  "iso-heavy"),
    ("OC(=O)C(F)(C(F)(F)F)C(F)(C(F)(F)F)C(F)(F)F",                                                     "double-CF3 br",  "iso-heavy"),

    # ── FG-eccentric – same PF component, COOH at different topological ────────
    #    distance from molecule centroid (central star vs. peripheral end)
    ("OC(=O)C(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)F",                                                     "C5-star (ctr)",  "FG-eccentric"),
    ("OC(=O)C(C(F)(F)C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)C(F)(F)F)F",                                       "C7-star (ctr)",  "FG-eccentric"),
]

P1_SMILES  = [c[0] for c in P1_COMPOUNDS]
P1_LABELS  = [c[1] for c in P1_COMPOUNDS]
P1_SERIES  = [c[2] for c in P1_COMPOUNDS]

# ══════════════════════════════════════════════════════════════════════════════
# Phase-2 compound library
# 18 compounds present in both test_set_for_PFASSTRUCTv5.tsv AND
# Richard2023_SI_TableS2.csv (verified by _check_overlap.py)
# ══════════════════════════════════════════════════════════════════════════════
P2_COMPOUNDS = [
    ("DTXSID6062599",   "PFPeA (C5 PFCA)",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID3031862",   "PFHxA (C6 PFCA)",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID8031865",   "PFOA  (C8 PFCA)",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID8031863",   "PFNA  (C9 PFCA)",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID3031860",   "PFDA  (C10 PFCA)",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID5030030",   "PFBS  (C4 PFSA)",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID7040150",   "PFHxS (C6 PFSA)",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID3031864",   "PFOS  (C8 PFSA)",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID30880251",  "PFBS-amide (C4 SA)", "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID3038939",   "PFOSA (C8 SA)",       "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("DTXSID90896595",  "iso-PFOA (C8 br)",    "OC(=O)C(F)(F)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
]

P2_DTXSIDS = [c[0] for c in P2_COMPOUNDS]
P2_LABELS  = [c[1] for c in P2_COMPOUNDS]
P2_SMILES  = [c[2] for c in P2_COMPOUNDS]
# Chemical-family groupings for Phase-2 MDS hull visualisation
# Series 1 (PFCA): PFPeA, PFHxA, PFOA, PFNA, PFDA  → chain-length axis
# Series 2 (PFSA): PFBS, PFHxS, PFOS                → same chain_size axis but SO3H
# Series 3 (PFSA-amide): PFBS-amide, PFOSA          → SO2NH2 head group
# Series 4 (Branched): iso-PFOA                     → branching axis
P2_SERIES = [
    "PFCA",             # PFPeA  (C5)
    "PFCA",             # PFHxA  (C6)
    "PFCA",             # PFOA   (C8)
    "PFCA",             # PFNA   (C9)
    "PFCA",             # PFDA   (C10)
    "PFSA",             # PFBS   (C4)
    "PFSA",             # PFHxS  (C6)
    "PFSA",             # PFOS   (C8)
    "PFSA-amide",       # PFBS-amide (C4)
    "PFSA-amide",       # PFOSA  (C8)
    "Branched",         # iso-PFOA (C8 branched)
]

# ══════════════════════════════════════════════════════════════════════════════
# Series colours
# ══════════════════════════════════════════════════════════════════════════════
_SERIES_COLOURS = {
    "PFAA":        "#2166ac",
    "PFSA":        "#d73027",
    "iso-PFCA":    "#4dac26",
    "iso-PFSA":    "#7fbf7b",
    "HFPO":        "#762a83",
    "FTOH":        "#e08214",
    "FTS":         "#f4a582",
    # new series
    "FTCA":        "#a65628",
    "PFPA":        "#984ea3",
    "iso-heavy":   "#1b9e77",
    "FG-eccentric":"#e41a1c",
}
_P2_SERIES_COLOURS = {
    "PFCA":       "#2166ac",   # blue  — perfluorocarboxylic acids
    "PFSA":       "#d73027",   # red   — perfluorosulfonic acids
    "PFSA-amide": "#984ea3",   # purple — perfluorosulfonamides
    "Branched":   "#4dac26",   # green — branched PFCA
}

# ══════════════════════════════════════════════════════════════════════════════
# Extended compound set (Phase 2B) – PFASGroups-only (no TxP-PFAS CSV needed)
# Groups span classic ionic PFAS, telomer family, sulfonamide PFAS, and
# branched variants to test discrimination scaling with more chemical families.
# ══════════════════════════════════════════════════════════════════════════════
P2B_COMPOUNDS = [
    # ── PFAA (carboxylic acids) C4 / C6 / C8 ─────────────────────────────────
    ("PFBA (C4 CA)",    "PFAA",        "OC(=O)C(F)(F)C(F)(F)C(F)(F)F"),
    ("PFHxA (C6 CA)",   "PFAA",        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("PFOA (C8 CA)",    "PFAA",        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    # ── PFSA (sulfonic acids) C4 / C6 / C8 ───────────────────────────────────
    ("PFBS (C4 SA)",    "PFSA",        "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("PFHxS (C6 SA)",   "PFSA",        "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("PFOS (C8 SA)",    "PFSA",        "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    # ── Fluorotelomer family (FT-chain, various head groups) ──────────────────
    ("4:2 FTOH",        "Telomer",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("6:2 FTOH",        "Telomer",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("6:2 FTCA",        "Telomer",     "OC(=O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("6:2 FTS",         "Telomer",     "OCCS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    # ── Sulfonamide-PFAS (SO2-N headgroup – structurally distinct) ────────────
    ("FOSA (C8)",       "Sulfonamide", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)N"),
    ("N-MeFOSAA",       "Sulfonamide", "OC(=O)CN(C)S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    # ── Branched PFAS (same chain count as linear counterparts) ───────────────
    ("iso-PFOA (α-CF3)","Branched",    "OC(=O)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("gem-PFCA (C5)",   "Branched",    "OC(=O)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F"),
    # ── Perfluoroalkyl phosphonic acid (PFPA) ─────────────────────────────────
    ("C6 PFPA",         "PFPA",        "OP(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    ("C8 PFPA",         "PFPA",        "OP(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    # ── HFPO-type (ether-chain PFAS, e.g. GenX extended) ─────────────────────
    ("HFPO-DA (GenX)",  "Ether-PFAS",  "OC(=O)C(F)(C(F)(F)F)OC(F)(F)F"),
    ("HFPO-TA (trimer)","Ether-PFAS",  "OC(=O)C(F)(OC(F)(F)C(F)(OC(F)(F)F)(F)F)F"),
]

P2B_LABELS  = [c[0] for c in P2B_COMPOUNDS]
P2B_SERIES  = [c[1] for c in P2B_COMPOUNDS]
P2B_SMILES  = [c[2] for c in P2B_COMPOUNDS]

_P2B_SERIES_COLOURS = {
    "PFAA":        "#2166ac",
    "PFSA":        "#d73027",
    "Telomer":     "#e08214",
    "Sulfonamide": "#984ea3",
    "Branched":    "#4dac26",
    "PFPA":        "#a65628",
    "Ether-PFAS":  "#762a83",
}

# ══════════════════════════════════════════════════════════════════════════════
# All supported metrics (must match fingerprints.py constants)
# ══════════════════════════════════════════════════════════════════════════════
ALL_GRAPH_METRICS = [
    'branching', 'mean_eccentricity', 'median_eccentricity', 'diameter',
    'radius', 'effective_graph_resistance', 'component_fraction',
    'min_dist_to_center', 'max_dist_to_periphery', 'min_dist_to_barycenter',
    'min_resistance_dist_to_barycenter', 'min_resistance_dist_to_center',
    'max_resistance_dist_to_periphery', 'size',
]

ALL_MOL_METRICS = [
    'n_components', 'total_size', 'mean_size', 'max_size',
    'mean_branching', 'max_branching', 'mean_eccentricity',
    'max_diameter', 'mean_component_fraction', 'max_component_fraction',
]

# ══════════════════════════════════════════════════════════════════════════════
# Fingerprint configurations (39 variants)
# ══════════════════════════════════════════════════════════════════════════════
def _build_configs() -> list[tuple[str, str, dict]]:
    """Return list of (section, label, kwargs)."""
    cfgs: list[tuple[str, str, dict]] = []
    for mode in ('binary', 'count', 'max_component', 'total_component'): #
        cfgs.append(('Count modes', mode, dict(component_metrics=[mode])))
    for m in ALL_GRAPH_METRICS:
        cfgs.append(('Individual graph metrics', f'binary+{m}',
                     dict(component_metrics=['binary', m])))
    cfgs.append(('Metric combinations', f'binary+ALL_graph ({len(ALL_GRAPH_METRICS)})',
                 dict(component_metrics=['binary'] + list(ALL_GRAPH_METRICS))))
    for combo in [
        ['branching', 'diameter'],
        ['branching', 'effective_graph_resistance'],
        ['branching', 'mean_eccentricity', 'diameter'],
        ['branching', 'mean_eccentricity', 'effective_graph_resistance'],
        ['branching', 'mean_eccentricity', 'diameter', 'radius', 'effective_graph_resistance'],
    ]:
        cfgs.append(('Metric combinations', 'binary+[' + ','.join(combo) + ']',
                     dict(component_metrics=['binary'] + combo)))
    for combo in [['branching'], ['branching', 'diameter'],
                  ['branching', 'mean_eccentricity', 'effective_graph_resistance']]:
        cfgs.append(('Count + graph combos', 'total+[' + ','.join(combo) + ']',
                     dict(component_metrics=['total_component'] + combo)))
    for m in ALL_MOL_METRICS:
        cfgs.append(('Individual mol metrics', f'binary+mol:{m}',
                     dict(component_metrics=['binary'], molecule_metrics=[m])))
    cfgs.append(('Metric combinations', f'binary+ALL_mol ({len(ALL_MOL_METRICS)})',
                 dict(component_metrics=['binary'], molecule_metrics=list(ALL_MOL_METRICS))))
    cfgs.append(('Metric combinations', 'FULL (total+ALL_graph+ALL_mol)',
                 dict(component_metrics=['total_component'] + list(ALL_GRAPH_METRICS),
                      molecule_metrics=list(ALL_MOL_METRICS))))
    return cfgs


FP_CONFIGS = _build_configs()

# ══════════════════════════════════════════════════════════════════════════════
# Tanimoto / similarity helpers
# ══════════════════════════════════════════════════════════════════════════════
def tanimoto_continuous(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a.astype(float), b.astype(float)
    num = float(np.sum(np.minimum(a, b)))
    den = float(np.sum(np.maximum(a, b)))
    return num / den if den > 0 else 0.0


def pairwise_tanimoto(X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    sim = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        sim[i, i] = 1.0
        for j in range(i + 1, n):
            s = tanimoto_continuous(X[i], X[j])
            sim[i, j] = sim[j, i] = s
    return sim


def jaccard_binary(a: np.ndarray, b: np.ndarray) -> float:
    a_b, b_b = a.astype(bool), b.astype(bool)
    inter = int(np.sum(a_b & b_b))
    union = int(np.sum(a_b | b_b))
    return inter / union if union > 0 else 0.0


def pairwise_jaccard(X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    sim = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        sim[i, i] = 1.0
        for j in range(i + 1, n):
            s = jaccard_binary(X[i], X[j])
            sim[i, j] = sim[j, i] = s
    return sim


def cluster_order(sim: np.ndarray) -> np.ndarray:
    dist = np.clip(1.0 - sim, 0.0, None)
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    condensed = np.clip(condensed, 0.0, None)
    try:
        order = leaves_list(linkage(condensed, method='average'))
    except Exception:
        order = np.arange(len(sim))
    return order


def off_diag_stats(sim: np.ndarray) -> tuple[float, float, float]:
    """Return (mean, min, std) of off-diagonal entries."""
    n = sim.shape[0]
    mask = ~np.eye(n, dtype=bool)
    vals = sim[mask]
    return float(vals.mean()), float(vals.min()), float(vals.std())


# ══════════════════════════════════════════════════════════════════════════════
# Colour / style helpers
# ══════════════════════════════════════════════════════════════════════════════
_CMAP = LinearSegmentedColormap.from_list(
    "tanimoto", ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"])


def _fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('ascii')


def _section_colours(sections: list[str]) -> dict[str, tuple]:
    unique = list(dict.fromkeys(sections))
    cmap_ = plt.cm.tab10
    return {s: cmap_(i / max(len(unique) - 1, 1)) for i, s in enumerate(unique)}


# ══════════════════════════════════════════════════════════════════════════════
# Fingerprint computation
# ══════════════════════════════════════════════════════════════════════════════
def compute_all_fingerprints(smiles: list[str],
                             fp_configs: list[tuple[str, str, dict]]) -> list[dict]:
    """Parse all SMILES once, then compute each FP config."""
    from PFASGroups import parse_smiles

    parsed = []
    for smi in smiles:
        try:
            parsed.append(parse_smiles(smi, halogens='F'))
        except Exception as e:
            print(f"  WARN: could not parse '{smi[:30]}': {e}", file=sys.stderr)
            parsed.append(None)

    results = []
    for section, label, kwargs in fp_configs:
        vecs = []
        for i, res in enumerate(parsed):
            if res is None:
                vecs.append(None)
                continue
            try:
                fp_obj = res.to_fingerprint(**kwargs)
                # to_fingerprint() returns a ResultsFingerprint; extract the 1-D row vector
                arr = np.array(fp_obj.fingerprints, dtype=np.float64)
                vecs.append(arr[0] if arr.ndim == 2 else arr)
            except Exception as e:
                print(f"  WARN: FP '{label}' mol {i}: {e}", file=sys.stderr)
                vecs.append(None)
        valid = [v for v in vecs if v is not None]
        if not valid:
            continue
        width = max(len(v) for v in valid)
        mat = np.zeros((len(smiles), width), dtype=np.float64)
        for i, v in enumerate(vecs):
            if v is not None:
                mat[i, :len(v)] = np.array(v, dtype=np.float64)
        sim = pairwise_tanimoto(mat)
        mean_t, min_t, std_t = off_diag_stats(sim)
        results.append(dict(
            section=section, label=label, kwargs=kwargs,
            mat=mat, sim=sim,
            mean_t=mean_t, min_t=min_t, std_t=std_t,
            n_cols=width,
        ))
    return results


# ══════════════════════════════════════════════════════════════════════════════
# TxP-PFAS CSV loader
# ══════════════════════════════════════════════════════════════════════════════
def _find_dtxsid_header_row(csvpath: Path) -> int:
    with open(csvpath, encoding='utf-8', errors='replace') as fh:
        for i, line in enumerate(fh):
            if 'DTXSID' in line and 'pfas' in line.lower():
                return i
    return 0


def load_txppfas_csv(csvpath: Path) -> tuple[np.ndarray, list[str], list[str]]:
    skip = _find_dtxsid_header_row(csvpath)
    df = pd.read_csv(csvpath, skiprows=skip, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    dtx_col = next(c for c in df.columns if 'DTXSID' in c.upper())
    feature_cols = [c for c in df.columns if c != dtx_col]
    dtxsids = df[dtx_col].tolist()
    X = df[feature_cols].fillna('0').replace('', '0').astype(float).values.astype(np.float32)
    return X, feature_cols, dtxsids


def load_smiles_lookup(csvpath: Path) -> dict[str, str]:
    """Load a DTXSID -> SMILES mapping from a CSV that has DTXSID and SMILES columns.

    Supports the Richard2023_SI_TableS5.csv format (and similar CompTox exports).
    Returns a dict keyed by DTXSID with canonical SMILES as values.
    """
    skip = _find_dtxsid_header_row(csvpath)
    df = pd.read_csv(csvpath, skiprows=skip, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    dtx_col = next(c for c in df.columns if 'DTXSID' in c.upper())
    smi_col  = next(
        (c for c in df.columns if c.upper() == 'SMILES'),
        next((c for c in df.columns if 'SMILES' in c.upper()), None)
    )
    if smi_col is None:
        raise ValueError(f"No SMILES column found in {csvpath.name}")
    return {
        dtx: smi
        for dtx, smi in zip(df[dtx_col], df[smi_col])
        if pd.notna(dtx) and pd.notna(smi) and str(smi).strip()
    }


# ══════════════════════════════════════════════════════════════════════════════
# Phase-1 figures
# ══════════════════════════════════════════════════════════════════════════════
def plot_p1_discrimination(all_data: list[dict], dpi: int = 130) -> plt.Figure:
    """Bar chart: all FP configs sorted by mean off-diagonal Tanimoto (lower = better)."""
    srt = sorted(all_data, key=lambda d: d['mean_t'])
    labels   = [d['label']   for d in srt]
    means    = [d['mean_t']  for d in srt]
    mins     = [d['min_t']   for d in srt]
    sections = [d['section'] for d in srt]
    sec_col  = _section_colours(sections)
    colours  = [sec_col[s] for s in sections]
    n = len(labels)
    h = max(5, n * 0.28)
    fig, ax = plt.subplots(figsize=(11, h), dpi=dpi)
    y = np.arange(n)
    ax.barh(y, means, color=colours, alpha=0.85)
    ax.scatter(mins, y, marker='|', s=80, color='k', zorder=3, label='min T')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Tanimoto similarity (lower = more discriminating)')
    ax.set_title(f'Phase 1 – Fingerprint discrimination (all {len(FP_CONFIGS)} configs)\n'
                 f'n={len(P1_SMILES)} compounds spanning chain length, branching, and functional group')
    patches = [mpatches.Patch(color=col, label=sec) for sec, col in sec_col.items()]
    ax.legend(handles=patches, title='Section', loc='lower right', fontsize=7)
    for rank in range(5):
        ax.annotate(f'#{rank+1}', xy=(means[rank], rank),
                    xytext=(means[rank] + 0.005, rank),
                    fontsize=7, color='#d73027', fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    ax.invert_yaxis()
    fig.tight_layout()
    return fig


def plot_p1_metric_impact(all_data: list[dict], dpi: int = 130) -> plt.Figure:
    """Compare each individual graph / mol metric against binary baseline."""
    base      = next(d for d in all_data if d['label'] == 'binary')
    base_mean = base['mean_t']
    graph_rows = sorted([d for d in all_data if d['section'] == 'Individual graph metrics'],
                        key=lambda d: d['mean_t'])
    mol_rows   = sorted([d for d in all_data if d['section'] == 'Individual mol metrics'],
                        key=lambda d: d['mean_t'])
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=dpi)
    for ax, rows, title in [
        (axes[0], graph_rows, 'Individual graph metrics\n(appended to binary FP)'),
        (axes[1], mol_rows,   'Individual mol metrics\n(appended to binary FP)'),
    ]:
        labels = [d['label'].split('binary+')[-1].replace('mol:', '') for d in rows]
        delta  = [base_mean - d['mean_t'] for d in rows]
        cols   = ['#2166ac' if v > 0 else '#d73027' for v in delta]
        y = np.arange(len(rows))
        ax.barh(y, delta, color=cols, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color='k', lw=1)
        ax.set_xlabel('Delta mean Tanimoto vs binary baseline\n'
                      '(positive = more discriminating; negative = less)')
        ax.set_title(title)
        ax.grid(axis='x', alpha=0.3)
        ax.invert_yaxis()
    fig.suptitle('Phase 1 – Per-metric impact relative to binary baseline', fontsize=11)
    fig.tight_layout()
    return fig


def _single_heatmap(ax: plt.Axes, sim: np.ndarray, labels: list[str], title: str) -> None:
    order = cluster_order(sim)
    s = sim[np.ix_(order, order)]
    lbl = [labels[i] for i in order]
    n = len(lbl)
    im = ax.imshow(s, cmap=_CMAP, vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(n))
    ax.set_xticklabels(lbl, rotation=90, fontsize=5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(lbl, fontsize=5)
    ax.set_title(title, fontsize=7)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_p1_heatmaps(all_data: list[dict], top_n: int = 8, dpi: int = 100) -> plt.Figure:
    """Heatmap grid for binary baseline + top-N configs."""
    srt  = sorted(all_data, key=lambda d: d['mean_t'])
    base = next(d for d in all_data if d['label'] == 'binary')
    selected = [base] + [d for d in srt if d['label'] != 'binary'][:top_n]
    ncols = 3
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), dpi=dpi)
    axes_flat = np.array(axes).flatten()
    for i, data in enumerate(selected):
        _single_heatmap(axes_flat[i], data['sim'], P1_LABELS,
                        f"{data['label']}\nmean T={data['mean_t']:.3f}")
    for ax in axes_flat[len(selected):]:
        ax.axis('off')
    fig.suptitle('Phase 1 – Tanimoto heatmaps (binary baseline + top 8 configs)', fontsize=11)
    fig.tight_layout()
    return fig


def plot_p1_mds(all_data: list[dict], top_n: int = 6, dpi: int = 100) -> plt.Figure:
    """2-D MDS projections for binary + top-N configs, with per-series convex-hull shading."""
    srt  = sorted(all_data, key=lambda d: d['mean_t'])
    base = next(d for d in all_data if d['label'] == 'binary')
    selected = [base] + [d for d in srt if d['label'] != 'binary'][:top_n]
    ncols = min(4, len(selected))
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 4.2 * nrows), dpi=dpi)
    axes_flat = np.array(axes).flatten()
    point_colours = [_SERIES_COLOURS.get(s, '#888888') for s in P1_SERIES]
    for i, data in enumerate(selected):
        dist = np.clip(1.0 - data['sim'], 0, None)
        np.fill_diagonal(dist, 0.0)
        try:
            coords = _classical_mds(dist)
        except Exception:
            axes_flat[i].axis('off')
            continue
        ax = axes_flat[i]
        # Convex-hull shading per series
        _draw_group_hulls(ax, coords, P1_SERIES, _SERIES_COLOURS)
        ax.scatter(coords[:, 0], coords[:, 1], c=point_colours,
                   s=55, zorder=3, edgecolors='white', linewidths=0.4)
        for j, lbl in enumerate(P1_LABELS):
            ax.annotate(lbl, coords[j], textcoords='offset points',
                        xytext=(3, 3), fontsize=5, zorder=4)
        gsi = _group_separation_score(coords, P1_SERIES)
        ax.set_title(
            f"{data['label']}\nmean T={data['mean_t']:.3f}  |  GSI={gsi:.2f}",
            fontsize=7)
        ax.axis('equal')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(alpha=0.25)
    for ax in axes_flat[len(selected):]:
        ax.axis('off')
    # Shared legend on the first visible panel
    handles = [mpatches.Patch(color=col, label=nm) for nm, col in _SERIES_COLOURS.items()
               if nm in set(P1_SERIES)]
    axes_flat[0].legend(handles=handles, fontsize=6, loc='best',
                        framealpha=0.85, edgecolor='#cccccc')
    fig.suptitle(
        'Phase 1 – MDS projections (binary baseline + top 6 configs)\n'
        'Shaded regions = series hulls.  GSI = group-separation index (higher = better).',
        fontsize=10)
    fig.tight_layout()
    return fig


def plot_p1_ranking(all_data: list[dict], dpi: int = 130) -> plt.Figure:
    """Compound-pair Tanimoto heatmap for top-20 configs."""
    top20 = sorted(all_data, key=lambda d: d['mean_t'])[:20]
    labels = [d['label'] for d in top20]
    n_mols = len(P1_SMILES)
    pairs = [(i, j) for i in range(n_mols) for j in range(i + 1, n_mols)]
    pair_labels = [f"{P1_LABELS[i][:9]}|{P1_LABELS[j][:9]}" for i, j in pairs]
    mat = np.zeros((len(top20), len(pairs)))
    for ci, data in enumerate(top20):
        for pi, (i, j) in enumerate(pairs):
            mat[ci, pi] = data['sim'][i, j]
    pair_order = np.argsort(mat.mean(axis=0))
    fig, ax = plt.subplots(figsize=(max(13, len(pairs) * 0.15), 6), dpi=dpi)
    im = ax.imshow(mat[:, pair_order], cmap=_CMAP, vmin=0, vmax=1, aspect='auto')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    step = max(1, len(pairs) // 40)
    shown = list(range(0, len(pairs), step))
    tick_pos = [int(np.where(pair_order == p)[0][0]) for p in shown if p in pair_order]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([pair_labels[p] for p in shown[:len(tick_pos)]],
                       rotation=90, fontsize=5)
    plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    ax.set_title('Phase 1 – Per-pair Tanimoto for top-20 configs (sorted by pair mean similarity)',
                 fontsize=9)
    fig.tight_layout()
    return fig


def plot_p1_gsi_breakdown(all_data: list[dict], top_n: int = 6,
                           dpi: int = 110) -> plt.Figure:
    """Per-group separability breakdown for the top-N configs + binary baseline.

    Shows which chemical series each config separates well or poorly.
    Answers: 'is one config's high global GSI driven by just one series?'
    """
    srt  = sorted(all_data, key=lambda d: d['mean_t'])
    base = next(d for d in all_data if d['label'] == 'binary')
    selected = [base] + [d for d in srt if d['label'] != 'binary'][:top_n]

    # Compute per-group GSI for each selected config
    all_group_names: list[str] = []
    records: list[dict] = []     # list of {label, group, score}
    for data in selected:
        dist = np.clip(1.0 - data['sim'], 0, None)
        np.fill_diagonal(dist, 0.0)
        try:
            coords = _classical_mds(dist)
        except Exception:
            continue
        pg = _per_group_gsi(coords, P1_SERIES)
        for grp, score in pg.items():
            records.append({'label': data['label'], 'group': grp, 'score': score})
            if grp not in all_group_names:
                all_group_names.append(grp)

    if not records:
        fig, ax = plt.subplots(figsize=(6, 3), dpi=dpi)
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return fig

    import pandas as _pd
    df = _pd.DataFrame(records)
    groups   = sorted(df['group'].unique())
    configs  = [d['label'] for d in selected]

    n_groups  = len(groups)
    n_configs = len(configs)
    x = np.arange(n_groups)
    bar_w = 0.8 / n_configs
    cmap_ = plt.cm.tab10

    fig, ax = plt.subplots(figsize=(max(10, n_groups * 1.2), 5.5), dpi=dpi)
    for ci, cfg_label in enumerate(configs):
        sub = df[df['label'] == cfg_label].set_index('group')
        scores = [float(sub.loc[g, 'score']) if g in sub.index else 0.0 for g in groups]
        offset = (ci - n_configs / 2 + 0.5) * bar_w
        bars = ax.bar(x + offset, scores, bar_w * 0.9,
                      color=cmap_(ci / max(n_configs - 1, 1)),
                      alpha=0.85, label=cfg_label)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=25, ha='right', fontsize=8)
    ax.set_ylabel('Per-group separability\n(mean dist to other centroids / within-group spread)')
    ax.set_title(
        f'Phase 1 – Per-group GSI breakdown (binary + top {top_n} configs)\n'
        'Higher = that series is well separated by this fingerprint config.',
        fontsize=9)
    ax.legend(fontsize=7, ncol=2, loc='upper right', framealpha=0.85)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Phase-2 figures
# ══════════════════════════════════════════════════════════════════════════════
def plot_p2_heatmaps(p2_results: list[dict], txp_sim: np.ndarray,
                     mol_labels: list[str], dpi: int = 100) -> plt.Figure:
    txp_mean_t = off_diag_stats(txp_sim)[0]
    all_items = p2_results + [dict(label='TxP-PFAS (129-bit)', sim=txp_sim,
                                   mean_t=txp_mean_t)]
    ncols = 3
    nrows = int(np.ceil(len(all_items) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 5 * nrows), dpi=dpi)
    axes_flat = np.array(axes).flatten()
    for i, item in enumerate(all_items):
        _single_heatmap(axes_flat[i], item['sim'], mol_labels,
                        f"{item['label']}\nmean T={item['mean_t']:.3f}")
    for ax in axes_flat[len(all_items):]:
        ax.axis('off')
    fig.suptitle(
        f'Phase 2 – Pairwise Tanimoto heatmaps\n'
        f'top-5 PFASGroups configs vs TxP-PFAS  |  n={len(mol_labels)} compounds',
        fontsize=11)
    fig.tight_layout()
    return fig


def plot_p2_discrimination(p2_results: list[dict], txp_sim: np.ndarray,
                           dpi: int = 130) -> tuple[plt.Figure, dict[str, tuple[float, float]]]:
    """Bar chart comparing discrimination (mean & min off-diagonal Tanimoto).

    Both PFASGroups configs and TxP-PFAS are compared on the same compound set.
    Lower mean T = more discriminating.  TxP-PFAS is shown as a competitor.
    Returns (fig, disc_dict) where disc_dict maps label -> (mean_t, min_t).
    """
    txp_mean, txp_min, _ = off_diag_stats(txp_sim)
    all_items = p2_results + [dict(label='TxP-PFAS (129-bit)',
                                   mean_t=txp_mean, min_t=txp_min)]
    srt = sorted(all_items, key=lambda d: d['mean_t'])
    labels = [d['label'] for d in srt]
    means  = [d['mean_t'] for d in srt]
    mins   = [d['min_t']  for d in srt]
    colours = ['#d62728' if 'TxP-PFAS' in lb else '#2166ac' for lb in labels]

    disc_dict: dict[str, tuple[float, float]] = {
        d['label']: (d['mean_t'], d['min_t']) for d in all_items
    }

    n = len(labels)
    fig, ax = plt.subplots(figsize=(11, max(4, n * 0.55)), dpi=dpi)
    y = np.arange(n)
    ax.barh(y, means, color=colours, alpha=0.85)
    ax.scatter(mins, y, marker='|', s=90, color='k', zorder=3, label='min T')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('Mean off-diagonal Tanimoto (lower = more discriminating)')
    ax.set_title(
        f'Phase 2 – Discrimination comparison: top-5 PFASGroups configs vs TxP-PFAS\n'
        f'n={p2_results[0]["sim"].shape[0] if p2_results else 0} compounds  |  '
        f'lower bar = more discriminating'
    )
    patches = [
        mpatches.Patch(color='#2166ac', label='PFASGroups config'),
        mpatches.Patch(color='#d62728', label='TxP-PFAS (129-bit)'),
    ]
    ax.legend(handles=patches + [mpatches.Patch(color='k', label='min T marker')],
              fontsize=8, loc='lower right')
    ax.grid(axis='x', alpha=0.3)
    ax.invert_yaxis()
    fig.tight_layout()
    return fig, disc_dict


def plot_p2_mds(p2_results: list[dict], txp_sim: np.ndarray,
                mol_labels: list[str], dpi: int = 100) -> plt.Figure:
    """Phase-2 MDS with per-family convex-hull shading and group-separation index."""
    all_items = p2_results + [dict(label='TxP-PFAS (129-bit)', sim=txp_sim,
                                   mean_t=off_diag_stats(txp_sim)[0])]
    # Build series list aligned to mol_labels (may be a matched subset of P2_COMPOUNDS)
    label_to_series: dict[str, str] = dict(zip(P2_LABELS, P2_SERIES))
    series_for_labels = [label_to_series.get(lbl, 'Other') for lbl in mol_labels]
    point_colours = [_P2_SERIES_COLOURS.get(s, '#888888') for s in series_for_labels]

    ncols = 3
    nrows = int(np.ceil(len(all_items) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.8 * nrows), dpi=dpi)
    axes_flat = np.array(axes).flatten()
    for i, data in enumerate(all_items):
        dist = np.clip(1.0 - data['sim'], 0, None)
        np.fill_diagonal(dist, 0.0)
        try:
            coords = _classical_mds(dist)
        except Exception:
            axes_flat[i].axis('off')
            continue
        ax = axes_flat[i]
        if 'TxP-PFAS' in data['label']:
            ax.set_facecolor('#fff5f5')
        # Convex-hull shading per family
        _draw_group_hulls(ax, coords, series_for_labels, _P2_SERIES_COLOURS)
        ax.scatter(coords[:, 0], coords[:, 1], c=point_colours,
                   s=55, zorder=3, edgecolors='white', linewidths=0.4)
        for j, lbl in enumerate(mol_labels):
            ax.annotate(lbl[:22], coords[j], textcoords='offset points',
                        xytext=(3, 3), fontsize=4.5, zorder=4)
        gsi = _group_separation_score(coords, series_for_labels)
        ax.set_title(
            f"{data['label']}\nmean T={data['mean_t']:.3f}  |  GSI={gsi:.2f}",
            fontsize=7)
        ax.axis('equal')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(alpha=0.25)
    for ax in axes_flat[len(all_items):]:
        ax.axis('off')
    # Shared legend on first panel
    hull_handles = [
        mpatches.Patch(facecolor=col, edgecolor=col, alpha=0.6, label=nm)
        for nm, col in _P2_SERIES_COLOURS.items()
        if nm in set(series_for_labels)
    ]
    axes_flat[0].legend(handles=hull_handles, fontsize=6, loc='best',
                        framealpha=0.85, edgecolor='#cccccc')
    fig.suptitle(
        'Phase 2 – MDS projections: top-5 PFASGroups configs vs TxP-PFAS\n'
        'Shaded regions = chemical-family hulls.  GSI = group-separation index (higher = better).',
        fontsize=10)
    fig.tight_layout()
    return fig


def plot_p2_gsi_breakdown(p2_results: list[dict], txp_sim: np.ndarray,
                           series_for_labels: list[str],
                           dpi: int = 110) -> plt.Figure:
    """Per-group separability breakdown for Phase 2 (TxP-PFAS vs top-5 PFASGroups).

    Lets you verify whether TxP-PFAS's higher global GSI is driven by
    a single well-separated group (e.g. FT surfactant) or is broadly superior.
    """
    txp_mean_t = off_diag_stats(txp_sim)[0]
    all_items = p2_results + [dict(label='TxP-PFAS (129-bit)', sim=txp_sim,
                                   mean_t=txp_mean_t)]

    import pandas as _pd
    records: list[dict] = []
    group_names: list[str] = []
    for item in all_items:
        dist = np.clip(1.0 - item['sim'], 0, None)
        np.fill_diagonal(dist, 0.0)
        try:
            coords = _classical_mds(dist)
        except Exception:
            continue
        pg = _per_group_gsi(coords, series_for_labels)
        for grp, score in pg.items():
            records.append({'label': item['label'], 'group': grp, 'score': score})
            if grp not in group_names:
                group_names.append(grp)

    if not records:
        fig, ax = plt.subplots(figsize=(6, 3), dpi=dpi)
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return fig

    df = _pd.DataFrame(records)
    groups  = sorted(df['group'].unique())
    configs = [d['label'] for d in all_items]
    n_groups  = len(groups)
    n_configs = len(configs)
    x = np.arange(n_groups)
    bar_w = 0.8 / n_configs

    cmap_ = plt.cm.tab10
    txp_color = '#d62728'

    fig, ax = plt.subplots(figsize=(max(9, n_groups * 1.3), 5.5), dpi=dpi)
    for ci, cfg_label in enumerate(configs):
        sub = df[df['label'] == cfg_label].set_index('group')
        scores = [float(sub.loc[g, 'score']) if g in sub.index else 0.0 for g in groups]
        offset = (ci - n_configs / 2 + 0.5) * bar_w
        color = txp_color if 'TxP-PFAS' in cfg_label else cmap_(ci / max(n_configs - 1, 1))
        ax.bar(x + offset, scores, bar_w * 0.9,
               color=color, alpha=0.85, label=cfg_label,
               edgecolor='k' if 'TxP-PFAS' in cfg_label else 'none',
               linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=20, ha='right', fontsize=8)
    ax.set_ylabel('Per-group separability\n(mean dist to other centroids / within-group spread)')
    ax.set_title(
        'Phase 2 – Per-group GSI breakdown: top-5 PFASGroups configs vs TxP-PFAS\n'
        'If TxP-PFAS (red, black outline) is only taller for one group, its overall\n'
        'higher GSI is driven by that group alone.',
        fontsize=9)
    ax.legend(fontsize=7, ncol=2, loc='upper right', framealpha=0.85)
    ax.grid(axis='y', alpha=0.3)

    # Clip y-axis so outlier groups don't squash the rest; annotate clipped bars
    all_scores = df['score'].values
    sorted_scores = np.sort(all_scores)
    clip_top = sorted_scores[int(len(sorted_scores) * 0.85)] * 1.35
    clip_top = max(clip_top, 0.1)
    ax.set_ylim(0, clip_top)
    # Annotate bars whose actual value exceeds the clip
    for bar_patch in ax.patches:
        actual = bar_patch.get_height()
        if actual > clip_top:
            ax.text(bar_patch.get_x() + bar_patch.get_width() / 2, clip_top * 0.97,
                    f'{actual:.2f}', ha='center', va='top', fontsize=6,
                    color='black', rotation=90,
                    bbox=dict(boxstyle='round,pad=0.2', fc='lightyellow', ec='grey', alpha=0.85))
            bar_patch.set_height(clip_top)
            bar_patch.set_edgecolor('red')
            bar_patch.set_linewidth(1.2)

    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Phase 2B – Extended compound set (PFASGroups-only, no TxP-PFAS CSV needed)
# ══════════════════════════════════════════════════════════════════════════════
def plot_p2b_mds(p2b_results: list[dict], dpi: int = 100) -> plt.Figure:
    """MDS on extended P2B compound set, PFASGroups configs only."""
    point_colours = [_P2B_SERIES_COLOURS.get(s, '#888888') for s in P2B_SERIES]
    ncols = min(3, len(p2b_results))
    nrows = int(np.ceil(len(p2b_results) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(5.2 * ncols, 4.8 * nrows), dpi=dpi)
    axes_flat = np.array(axes).flatten() if nrows * ncols > 1 else [axes]
    for i, data in enumerate(p2b_results):
        dist = np.clip(1.0 - data['sim'], 0, None)
        np.fill_diagonal(dist, 0.0)
        try:
            coords = _classical_mds(dist)
        except Exception:
            axes_flat[i].axis('off')
            continue
        ax = axes_flat[i]
        _draw_group_hulls(ax, coords, P2B_SERIES, _P2B_SERIES_COLOURS)
        ax.scatter(coords[:, 0], coords[:, 1], c=point_colours,
                   s=55, zorder=3, edgecolors='white', linewidths=0.4)
        for j, lbl in enumerate(P2B_LABELS):
            ax.annotate(lbl[:20], coords[j], textcoords='offset points',
                        xytext=(3, 3), fontsize=4.5, zorder=4)
        gsi = _group_separation_score(coords, P2B_SERIES)
        pg  = _per_group_gsi(coords, P2B_SERIES)
        worst = min(pg, key=pg.get) if pg else '?'
        best  = max(pg, key=pg.get) if pg else '?'
        ax.set_title(
            f"{data['label']}\nmean T={data['mean_t']:.3f}  GSI={gsi:.2f}\n"
            f"worst-sep: {worst}  best-sep: {best}",
            fontsize=6.5)
        ax.axis('equal')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(alpha=0.25)
    for ax in axes_flat[len(p2b_results):]:
        ax.axis('off')
    hull_handles = [
        mpatches.Patch(facecolor=col, edgecolor=col, alpha=0.6, label=nm)
        for nm, col in _P2B_SERIES_COLOURS.items()
        if nm in set(P2B_SERIES)
    ]
    if hull_handles:
        axes_flat[0].legend(handles=hull_handles, fontsize=6, loc='best',
                            framealpha=0.85, edgecolor='#cccccc')
    fig.suptitle(
        f'Phase 2B – Extended compound set ({len(P2B_SMILES)} cpds, '
        f'{len(set(P2B_SERIES))} families) – PFASGroups configs only\n'
        'Tests discrimination scaling: PFAA, PFSA, Telomer, Sulfonamide, Branched, PFPA, Ether-PFAS.',
        fontsize=9)
    fig.tight_layout()
    return fig


def plot_p2b_gsi_breakdown(p2b_results: list[dict], dpi: int = 110) -> plt.Figure:
    """Per-group GSI breakdown on the P2B extended compound set."""
    import pandas as _pd
    records: list[dict] = []
    for data in p2b_results:
        dist = np.clip(1.0 - data['sim'], 0, None)
        np.fill_diagonal(dist, 0.0)
        try:
            coords = _classical_mds(dist)
        except Exception:
            continue
        pg = _per_group_gsi(coords, P2B_SERIES)
        for grp, score in pg.items():
            records.append({'label': data['label'], 'group': grp, 'score': score})

    if not records:
        fig, ax = plt.subplots(figsize=(6, 3), dpi=dpi)
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        return fig

    df      = _pd.DataFrame(records)
    groups  = sorted(df['group'].unique())
    configs = df['label'].unique().tolist()
    x       = np.arange(len(groups))
    bar_w   = 0.8 / len(configs)
    cmap_   = plt.cm.tab10

    fig, ax = plt.subplots(figsize=(max(9, len(groups) * 1.3), 5.5), dpi=dpi)
    for ci, cfg_label in enumerate(configs):
        sub = df[df['label'] == cfg_label].set_index('group')
        scores = [float(sub.loc[g, 'score']) if g in sub.index else 0.0 for g in groups]
        offset = (ci - len(configs) / 2 + 0.5) * bar_w
        ax.bar(x + offset, scores, bar_w * 0.9,
               color=cmap_(ci / max(len(configs) - 1, 1)),
               alpha=0.85, label=cfg_label)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=20, ha='right', fontsize=9)
    ax.set_ylabel('Per-group separability score')
    ax.set_title(
        f'Phase 2B – Per-group GSI: {len(set(P2B_SERIES))} chemical families × top-5 PFASGroups configs\n'
        'Shows which families are consistently well/poorly separated across all configs.',
        fontsize=9)
    ax.legend(fontsize=7, ncol=2, loc='upper right', framealpha=0.85)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    return fig


def plot_p2b_discrimination(p2b_results: list[dict], dpi: int = 110) -> plt.Figure:
    """Bar chart – mean off-diagonal Tanimoto for top-5 configs on P2B compound set."""
    srt    = sorted(p2b_results, key=lambda d: d['mean_t'])
    labels = [d['label'] for d in srt]
    means  = [d['mean_t'] for d in srt]
    mins   = [d['min_t']  for d in srt]

    fig, ax = plt.subplots(figsize=(10, max(3, len(srt) * 0.6)), dpi=dpi)
    y = np.arange(len(srt))
    ax.barh(y, means, color='#2166ac', alpha=0.85)
    ax.scatter(mins, y, marker='|', s=90, color='k', zorder=3, label='min T')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Mean off-diagonal Tanimoto (lower = more discriminating)')
    ax.set_title(
        f'Phase 2B – Discrimination on extended compound set\n'
        f'n={len(P2B_SMILES)} cpds across {len(set(P2B_SERIES))} chemical families',
        fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(axis='x', alpha=0.3)
    ax.invert_yaxis()
    fig.tight_layout()
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# HTML report
# ══════════════════════════════════════════════════════════════════════════════
def build_html_report(
    p1_figs: dict[str, plt.Figure],
    p2_figs: dict[str, plt.Figure],
    p1_top5: list[dict],
    summary_df: pd.DataFrame,
    p2_disc: dict[str, tuple[float, float]],
    txp_mean_t: float,
    n_p2: int,
    p2_smiles: list[str] | None = None,
    p2_labels: list[str] | None = None,
    p2b_figs: dict[str, plt.Figure] | None = None,
) -> str:
    def embed(fig: plt.Figure) -> str:
        return (f'<img src="data:image/png;base64,{_fig_to_b64(fig)}" '
                f'style="max-width:100%;margin:8px 0">')

    def embed_struct(b64: str | None, caption: str = '') -> str:
        if b64 is None:
            return '<p class="note">Structure images require RDKit Draw (not available).</p>'
        out = f'<img src="data:image/png;base64,{b64}" style="max-width:100%;margin:8px 0">'
        if caption:
            out += f'<p style="font-size:.8em;color:#555;margin-top:2px">{caption}</p>'
        return out

    top5_rows = "".join(
        f"<tr><td>{r+1}</td><td><code>{d['label']}</code></td>"
        f"<td>{d['section']}</td><td>{d['mean_t']:.4f}</td>"
        f"<td>{d['min_t']:.4f}</td><td>{d['std_t']:.4f}</td><td>{d['n_cols']}</td></tr>"
        for r, d in enumerate(p1_top5)
    )
    disc_rows = "".join(
        f"<tr><td>{r+1}</td><td><code>{label}</code></td>"
        f"<td>{mean_t:.4f}</td><td>{min_t:.4f}</td></tr>"
        for r, (label, (mean_t, min_t)) in enumerate(
            sorted(p2_disc.items(), key=lambda kv: kv[1][0]))
    )
    summary_html = (
        summary_df.sort_values('mean_T')
        .to_html(index=False, float_format='%.4f', classes='summary-table')
    )

    # Pre-render structure grids
    p1_struct_b64 = _smiles_grid_b64(P1_SMILES, P1_LABELS, n_cols=6, sub_img_size=(220, 160))
    p2_struct_b64 = (_smiles_grid_b64(p2_smiles, p2_labels, n_cols=6, sub_img_size=(220, 160))
                     if p2_smiles and p2_labels else None)

    _mds_note = """
<div class="note">
<strong>How to read MDS plots:</strong>
Multi-Dimensional Scaling (MDS) converts the pairwise Tanimoto <em>distance</em> matrix
(1&nbsp;&minus;&nbsp;similarity) into a 2-D scatter plot while preserving inter-point distances
as faithfully as possible.  Two compounds placed close together means the fingerprint
cannot distinguish them; two compounds far apart means the fingerprint sees them as
structurally different.
<ul style="margin:.4em 0 0 1em">
<li><strong>Series separation</strong> &ndash; each chemical series (PFAA, PFSA, FTOH&hellip;) should form a
  compact, well-separated cluster; overlap between series = lost discrimination.</li>
<li><strong>Chain-length ordering</strong> &ndash; within a homologous series the points should be roughly
  collinear with C-number increasing monotonically along the axis.</li>
<li><strong>Branching resolution</strong> &ndash; iso-PFCA / iso-PFSA variants should sit away from their
  linear homologues, not merged with them.</li>
<li><strong>Fingerprint comparison (Phase 2)</strong> &ndash; a panel where all 18 compounds are
  spread far apart is more discriminating than one where everything clusters in the centre;
  this correlates directly with a low mean off-diagonal Tanimoto.</li>
</ul>
</div>"""

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>PFASGroups Tanimoto Benchmark</title>
<style>
  body{{font-family:system-ui,sans-serif;max-width:1400px;margin:0 auto;padding:20px;background:#fafbfc}}
  h1{{color:#1a237e}}h2{{color:#283593;border-bottom:2px solid #c5cae9;padding-bottom:4px}}
  h3{{color:#3949ab}}
  .meta{{color:#555;font-size:.9em;margin-bottom:20px}}
  table.summary-table{{border-collapse:collapse;font-size:.78em;width:100%}}
  table.summary-table th,table.summary-table td{{border:1px solid #ddd;padding:4px 8px}}
  table.summary-table th{{background:#3949ab;color:#fff}}
  table.summary-table tr:nth-child(even){{background:#f0f4ff}}
  table{{border-collapse:collapse;font-size:.85em;margin:12px 0}}
  th,td{{border:1px solid #ccc;padding:5px 10px}}
  th{{background:#e8eaf6}}
  .phase{{background:#fff;border:1px solid #c5cae9;border-radius:6px;padding:16px;margin:20px 0}}
  code{{background:#e8eaf6;padding:2px 5px;border-radius:3px;font-size:.88em}}
  .note{{background:#fff8e1;border-left:4px solid #ffc107;padding:8px 14px;margin:10px 0}}
  details summary{{cursor:pointer;font-weight:600;color:#3949ab;padding:4px 0}}
  details{{margin:10px 0}}
</style></head><body>
<h1>PFASGroups Fingerprint Tanimoto Benchmark</h1>
<p class="meta">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp;|&nbsp;
  <strong>Phase-1 compounds:</strong> {len(P1_SMILES)} &nbsp;|&nbsp;
  <strong>Phase-2 compounds:</strong> {n_p2} &nbsp;|&nbsp;
  <strong>Phase-2B compounds:</strong> {len(P2B_SMILES)} &nbsp;|&nbsp;
  <strong>FP configs tested:</strong> {len(FP_CONFIGS)}
  ({len(ALL_GRAPH_METRICS)} graph metrics, {len(ALL_MOL_METRICS)} mol metrics)</p>
<div class="note"><strong>Phase-1 compound series (11 chemical families):</strong>
  PFAA chain C2–C11 (10), PFSA C4/C6/C8/C10 (4), iso-PFCA alpha/beta-branch (6),
  iso-PFSA (2), HFPO-DA / GenX (1), FTOH 4:2–10:2 (4), FTS 6:2/8:2 (2),
  <strong>FTCA</strong> 4:2–8:2 (3) &ndash; FTOH scaffold with COOH head,
  <strong>PFPA</strong> C4/C6 phosphonic acids (2),
  <strong>iso-heavy</strong> gem/double/large-branch (4) &ndash; branching degree,
  <strong>FG-eccentric</strong> star-topology C5/C7 (2) &ndash; FG topological position.
  Lower mean off-diagonal Tanimoto = better discrimination.<br>
  <strong>Phase-2 compound set:</strong> {n_p2} structurally diverse PFAS present in both the
  Richard&nbsp;2023 validation TSV and the TxP-PFAS SI Table S2 CSV.
  TxP-PFAS (129-bit CSRML) is benchmarked <em>as a competitor</em>.
  TxP-PFAS mean T on this set: <strong>{txp_mean_t:.4f}</strong><br>
  <strong>Phase-2B extended set:</strong> {len(P2B_SMILES)} compounds across
  {len(set(P2B_SERIES))} families (PFAA, PFSA, Telomer, Sulfonamide, Branched, PFPA,
  Ether-PFAS) &ndash; PFASGroups only, no TxP-PFAS CSV required.</div>

<div class="phase">
<h2>Phase 1 – Fingerprint Selection</h2>
<h3>1.0 Phase-1 compound structures</h3>
<p>The {len(P1_SMILES)}&nbsp;compounds span 11 structural series chosen to challenge
each fingerprint on: (a)&nbsp;chain-length homologues (PFAA C2–C11, PFSA C4–C10),
(b)&nbsp;branching position (iso-PFCA α/β) and degree (iso-heavy: gem-CF₃, double-CF₃,
C₂F₅-branch), (c)&nbsp;functional-group variation at constant component size (FTCA vs FTOH
vs FTS; PFPA vs PFCA vs PFSA), (d)&nbsp;topological FG position / eccentricity
(FG-eccentric: star vs. linear).</p>
{embed_struct(p1_struct_b64, f'{len(P1_SMILES)} Phase-1 compounds &ndash; 6 per row')}
<h3>1.1 All configs ranked by discrimination</h3>{embed(p1_figs['discrimination'])}
<h3>1.2 Per-metric impact relative to binary baseline</h3>{embed(p1_figs['metric_impact'])}
<h3>1.3 Tanimoto heatmaps (binary + top 8)</h3>{embed(p1_figs['heatmaps'])}
<h3>1.4 MDS projections (binary + top 6)</h3>
{_mds_note}
{embed(p1_figs['mds'])}
<h3>1.5 Per-pair similarity for top-20 configs</h3>{embed(p1_figs['ranking'])}
<h3>1.6 Per-group GSI breakdown</h3>
<div class="note">This chart decomposes the overall Group Separation Index into
per-series contributions. A config that scores high globally but only for one series
may be less generically useful than one that scores more evenly across all series.</div>
{embed(p1_figs['gsi_breakdown'])}
<h3>Top-5 most discriminating configs</h3>
<table><tr><th>#</th><th>Label</th><th>Section</th>
<th>Mean T &darr;</th><th>Min T</th><th>Std T</th><th>n cols</th></tr>
{top5_rows}</table>
<details><summary>Full Phase-1 statistics table (all {len(summary_df)} configs)</summary>
{summary_html}
</details>
</div>

<div class="phase">
<h2>Phase 2 – Discrimination benchmark vs TxP-PFAS</h2>
<div class="note">TxP-PFAS (129-bit CSRML, Richard 2023) is included as a competitor.
  Lower mean T = more discriminating fingerprint.</div>
<h3>2.0 Phase-2 compound structures</h3>
<p>The {n_p2}&nbsp;compounds were selected as the overlap between the Richard&nbsp;2023
validation set and the TxP-PFAS SI Table S2 (14&thinsp;735 compounds).</p>
{embed_struct(p2_struct_b64, f'{n_p2} Phase-2 compounds &ndash; 6 per row')}
<h3>2.1 Pairwise Tanimoto heatmaps (top-5 PFASGroups configs + TxP-PFAS)</h3>{embed(p2_figs['heatmaps'])}
<h3>2.2 Discrimination comparison (mean T, lower = better)</h3>{embed(p2_figs['discrimination'])}
<table><tr><th>#</th><th>Fingerprint</th><th>Mean T &darr;</th><th>Min T</th></tr>
{disc_rows}</table>
<h3>2.3 MDS projections</h3>
{_mds_note}
{embed(p2_figs['mds'])}
<h3>2.4 Per-group GSI breakdown (PFASGroups configs vs TxP-PFAS)</h3>
<div class="note">Answers: <em>"Is TxP-PFAS's higher overall GSI driven by just one
well-separated group (e.g. FT surfactant), or does it broadly outperform?"</em>
If only one group's bar is taller for TxP-PFAS (red, black outline) while others
are comparable or lower, the global GSI advantage is group-specific, not general.</div>
{embed(p2_figs['gsi_breakdown'])}
</div>

{f'''<div class="phase">
<h2>Phase 2B – Extended compound set ({len(P2B_SMILES)} compounds, PFASGroups only)</h2>
<div class="note">Phase&nbsp;2B tests discrimination on a larger, more diverse set of
{len(P2B_SMILES)}&nbsp;compounds spanning {len(set(P2B_SERIES))}&nbsp;chemical families
(PFAA, PFSA, Telomer, Sulfonamide, Branched, PFPA, Ether-PFAS). No TxP-PFAS CSV
is needed. This section assesses how fingerprint discrimination scales with more compound
families, and identifies which families remain hard to separate.</div>
<h3>2B.1 Discrimination (mean T)</h3>{embed(p2b_figs["discrimination"])}
<h3>2B.2 MDS projections</h3>{embed(p2b_figs["mds"])}
<h3>2B.3 Per-group GSI breakdown</h3>{embed(p2b_figs["gsi_breakdown"])}
</div>''' if p2b_figs else ''}
</body></html>"""


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description='Two-phase Tanimoto benchmark')
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--dpi', type=int, default=110)
    ap.add_argument('--txppfas_csv', default=None)
    ap.add_argument('--smiles_csv', default=None,
                    help='CSV with DTXSID+SMILES columns used as a SMILES lookup table '
                         '(e.g. Richard2023_SI_TableS5.csv). '
                         'Resolves SMILES for P2 compounds by DTXSID; '
                         'falls back to the hardcoded SMILES when a DTXSID is absent.')
    args = ap.parse_args(argv)

    outdir = Path(args.outdir) if args.outdir else BENCHMARK_DIR / 'results' / 'tanimoto'
    outdir.mkdir(parents=True, exist_ok=True)

    txp_csv = (Path(args.txppfas_csv) if args.txppfas_csv
               else BENCHMARK_DIR / 'test_data' / 'Richard2023_SI_TableS2.csv')

    # ── SMILES lookup (Richard2023_SI_TableS5 or user-supplied) ──────────────
    _default_smiles_csv = BENCHMARK_DIR / 'test_data' / 'Richard2023_SI_TableS5.csv'
    smiles_csv = Path(args.smiles_csv) if args.smiles_csv else _default_smiles_csv
    smiles_lookup: dict[str, str] = {}
    if smiles_csv.exists():
        smiles_lookup = load_smiles_lookup(smiles_csv)
        print(f"Loaded SMILES lookup: {len(smiles_lookup):,} entries from {smiles_csv.name}")
    else:
        print(f"SMILES lookup CSV not found ({smiles_csv.name}); using hardcoded SMILES only")

    dpi = args.dpi

    # ── Phase 1 ──────────────────────────────────────────────────────────────
    print(f"Phase 1: computing {len(FP_CONFIGS)} FP configs on {len(P1_SMILES)} compounds …")
    p1_data   = compute_all_fingerprints(P1_SMILES, FP_CONFIGS)
    p1_sorted = sorted(p1_data, key=lambda d: d['mean_t'])
    top5      = p1_sorted[:5]

    summary_df = pd.DataFrame([
        dict(section=d['section'], label=d['label'],
             mean_T=d['mean_t'], min_T=d['min_t'], std_T=d['std_t'], n_cols=d['n_cols'])
        for d in p1_sorted
    ])
    summary_df.to_csv(outdir / 'tanimoto_summary.csv', index=False)
    print(f"  Top-5: {[d['label'] for d in top5]}")

    print("Phase 1: generating figures …")
    p1_figs = {
        'discrimination': plot_p1_discrimination(p1_data, dpi=dpi),
        'metric_impact':  plot_p1_metric_impact(p1_data, dpi=dpi),
        'heatmaps':       plot_p1_heatmaps(p1_data, top_n=8, dpi=dpi),
        'mds':            plot_p1_mds(p1_data, top_n=6, dpi=dpi),
        'ranking':        plot_p1_ranking(p1_data, dpi=dpi),
    }
    for name, fig in p1_figs.items():
        path = outdir / f'tanimoto_p1_{name}.png'
        fig.savefig(path, bbox_inches='tight')
        print(f"  Saved {path.name}")

    # ── Phase 2 ──────────────────────────────────────────────────────────────
    print(f"\nPhase 2: loading TxP-PFAS CSV from {txp_csv.name} …")
    txp_X, _txp_features, txp_dtxsids = load_txppfas_csv(txp_csv)
    dtx_index = {d: i for i, d in enumerate(txp_dtxsids)}

    # Match P2 compounds to CSV rows, resolving SMILES from lookup when available
    valid_rows: list[int] = []
    valid_smiles: list[str] = []
    valid_labels: list[str] = []
    for dtx, lbl, smi_fallback in zip(P2_DTXSIDS, P2_LABELS, P2_SMILES):
        if dtx in dtx_index:
            valid_rows.append(dtx_index[dtx])
            resolved_smi = smiles_lookup.get(dtx, smi_fallback)
            if resolved_smi != smi_fallback:
                print(f"  SMILES for {dtx} ({lbl}): using lookup ({resolved_smi[:40]}…)")
            valid_smiles.append(resolved_smi)
            valid_labels.append(lbl)
        else:
            print(f"  WARN: {dtx} ({lbl}) not found in CSV – skipped")

    print(f"  {len(valid_rows)}/{len(P2_DTXSIDS)} compounds matched in TxP-PFAS CSV")

    txp_rows = txp_X[valid_rows]
    txp_sim  = pairwise_jaccard(txp_rows)
    txp_mean_t = off_diag_stats(txp_sim)[0]
    print(f"  TxP-PFAS mean T on matched set: {txp_mean_t:.4f}")

    print("Phase 2: computing top-5 PFASGroups FPs …")
    top5_cfgs  = [(d['section'], d['label'], d['kwargs']) for d in top5]
    p2_results = compute_all_fingerprints(valid_smiles, top5_cfgs)

    print("Phase 2: generating figures …")
    # Series aligned to the matched label list
    label_to_series: dict[str, str] = dict(zip(P2_LABELS, P2_SERIES))
    valid_series = [label_to_series.get(lbl, 'Other') for lbl in valid_labels]

    disc_fig, p2_disc = plot_p2_discrimination(p2_results, txp_sim, dpi=dpi)
    p2_figs = {
        'heatmaps':       plot_p2_heatmaps(p2_results, txp_sim, valid_labels, dpi=dpi),
        'discrimination': disc_fig,
        'mds':            plot_p2_mds(p2_results, txp_sim, valid_labels, dpi=dpi),
        'gsi_breakdown':  plot_p2_gsi_breakdown(p2_results, txp_sim, valid_series, dpi=dpi),
    }
    for name, fig in p2_figs.items():
        path = outdir / f'tanimoto_p2_{name}.png'
        fig.savefig(path, bbox_inches='tight')
        print(f"  Saved {path.name}")

    # ── Phase 2B – extended compound set (PFASGroups only) ────────────────────
    print(f"\nPhase 2B: computing top-5 FPs on extended set ({len(P2B_SMILES)} compounds) …")
    p2b_results = compute_all_fingerprints(P2B_SMILES, top5_cfgs)

    print("Phase 2B: generating figures …")
    p2b_figs = {
        'mds':           plot_p2b_mds(p2b_results, dpi=dpi),
        'gsi_breakdown': plot_p2b_gsi_breakdown(p2b_results, dpi=dpi),
        'discrimination':plot_p2b_discrimination(p2b_results, dpi=dpi),
    }
    for name, fig in p2b_figs.items():
        path = outdir / f'tanimoto_p2b_{name}.png'
        fig.savefig(path, bbox_inches='tight')
        print(f"  Saved {path.name}")

    # ── HTML report ───────────────────────────────────────────────────────────
    print("\nBuilding HTML report …")
    # Phase 1 GSI breakdown figure (add after p2b so it can be included)
    p1_figs['gsi_breakdown'] = plot_p1_gsi_breakdown(p1_data, dpi=dpi)
    p1_figs['gsi_breakdown'].savefig(outdir / 'tanimoto_p1_gsi_breakdown.png', bbox_inches='tight')
    print("  Saved tanimoto_p1_gsi_breakdown.png")

    html = build_html_report(
        p1_figs, p2_figs, top5, summary_df, p2_disc,
        txp_mean_t, len(valid_rows),
        p2_smiles=valid_smiles,
        p2_labels=valid_labels,
        p2b_figs=p2b_figs,
    )
    report_path = outdir / 'tanimoto_report.html'
    report_path.write_text(html, encoding='utf-8')

    print("\n" + "=" * 70)
    print(f"Phase-1 top-5 configs (on {len(P1_SMILES)} compounds, 11 series):")
    for r, d in enumerate(top5):
        print(f"  #{r+1}  {d['label']:45s}  mean T={d['mean_t']:.4f}  min T={d['min_t']:.4f}")
    print(f"\nPhase-2 discrimination (mean T, lower = more discriminating):")
    for lbl, (mean_t, min_t) in sorted(p2_disc.items(), key=lambda kv: kv[1][0]):
        marker = '  <-- TxP-PFAS' if 'TxP-PFAS' in lbl else ''
        print(f"  {lbl:45s}  mean T={mean_t:.4f}  min T={min_t:.4f}{marker}")
    print(f"\nPhase-2B discrimination ({len(P2B_SMILES)} compounds, {len(set(P2B_SERIES))} families):")
    for d in sorted(p2b_results, key=lambda x: x['mean_t']):
        print(f"  {d['label']:45s}  mean T={d['mean_t']:.4f}")
    print(f"\nAll outputs -> {outdir}")
    print(f"HTML report -> {report_path}")
    print("=" * 70)
    plt.close('all')


if __name__ == '__main__':
    main()
