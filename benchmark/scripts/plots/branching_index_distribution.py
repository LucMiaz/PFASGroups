#!/usr/bin/env python3
"""
Distribution of branching index for the 300 PFAS molecules in
branching_test_results_2025.json, with 3 structural examples.

Layout
──────
  Top row   : distribution histogram + KDE spanning full width.
  Bottom row : 3 molecule structure panels (low / median / high branching).
  Dotted lines connect each vertical marker in the distribution to its
  corresponding molecule panel.

Saves
─────
  benchmark/imgs/branching_index_distribution.svg           — combined figure (generation BI)
  benchmark/imgs/branching_index_distribution_pfasgroups.svg — variant using PFASGroups mean_branching
  benchmark/imgs/branching_mol_{low,median,high}.svg         — individual structure SVGs

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/plots/branching_index_distribution.py
"""

from __future__ import annotations

import ast
import io
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
IMGS_DIR   = SCRIPT_DIR.parents[1] / "imgs"
IMGS_DIR.mkdir(exist_ok=True)

# ── PFASGroups import (soft dependency for variant plot) ──────────────────────
_PFASGROUPS_ROOT = str(SCRIPT_DIR.parents[2])  # plots -> scripts -> benchmark -> PFASGroups/ root
if _PFASGROUPS_ROOT not in sys.path:
    sys.path.insert(0, _PFASGROUPS_ROOT)
try:
    from PFASGroups.parser import parse_mol as _parse_mol  # type: ignore
    _PFASGROUPS_AVAILABLE = True
except ImportError:
    _PFASGROUPS_AVAILABLE = False
    print("Warning: PFASGroups not importable — variant plot will be skipped")

# ── Color scheme (from PFASGroups/data/color_scheme.yaml) ────────────────────
TINTS = {
    "orange_red":  ["#E15D0B", "#E67834", "#EB935C", "#F0AE85"],
    "blue":        ["#306DBA", "#5385C6", "#759ED1", "#97B6DD"],
    "magenta":     ["#9D206C", "#AD4585", "#BE6A9D", "#CE8FB5"],
    "dark_purple": ["#51127C", "#6E3A92", "#8B61A8", "#A888BD"],
}

# main color, light fill  —  one per example (low, median, high)
EX_PALETTE = [
    (TINTS["blue"][0],       TINTS["blue"][3]),
    (TINTS["orange_red"][0], TINTS["orange_red"][3]),
    (TINTS["magenta"][0],    TINTS["magenta"][3]),
]

# ── Seaborn / matplotlib style ────────────────────────────────────────────────
sns.set_theme(style="whitegrid", font="Ubuntu")
plt.rcParams.update(
    {
        "font.family":       "Ubuntu",
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "figure.facecolor":  "white",
        "axes.facecolor":    "white",
    }
)

# ── Load data ─────────────────────────────────────────────────────────────────
with open(DATA_DIR / "branching_test_results_2025.json") as fh:
    raw = json.load(fh)

mols = [r for r in raw if "branching_index" in r]
bi   = np.array([m["branching_index"] for m in mols])
print(f"Loaded {len(mols)} molecules  |  BI range [{bi.min():.4f}, {bi.max():.4f}]")

# ── Select 3 representative examples (10th / 50th / 90th percentile) ─────────
pcts   = [10, 50, 90]
labels = ["High branching", "Median branching", "Low branching"]
keys   = ["high", "median", "low"]

def _mol_atom_count(mol_data: dict) -> int:
    """Number of heavy atoms, used to prefer smaller molecules in example selection."""
    mol = Chem.MolFromSmiles(mol_data["smiles"])
    return mol.GetNumHeavyAtoms() if mol else 9999


def _pick_in_window(values: np.ndarray, mol_list: list, pct: int,
                    window: float, used: set[int]) -> int:
    """Return the index of the molecule closest to *pct* within ±window/2 percentile
    points, skipping indices already in *used*.  Among molecules equally close to
    the target value, the one with fewer heavy atoms is preferred."""
    lo = np.percentile(values, max(0,  pct - window / 2))
    hi = np.percentile(values, min(100, pct + window / 2))
    target = np.percentile(values, pct)
    candidates = [
        i for i in range(len(values))
        if lo <= values[i] <= hi and i not in used
    ]
    if not candidates:  # fall back to nearest if window is empty after exclusions
        candidates = [
            i for i in np.argsort(np.abs(values - target))
            if i not in used
        ]
    # Primary: closest to the percentile target value; secondary: fewest heavy atoms
    return min(candidates, key=lambda i: (
        abs(values[i] - target),
        _mol_atom_count(mol_list[i]),
    ))


used_idx: set[int] = set()
examples: list[dict] = []
for p, label, key, (main_c, fill_c) in zip(pcts, labels, keys, EX_PALETTE):
    idx = _pick_in_window(bi, mols, p, window=10, used=used_idx)
    used_idx.add(idx)
    examples.append(
        {
            "key":      key,
            "label":    label,
            "bi":       float(bi[idx]),
            "mol_data": mols[idx],
            "color":    main_c,
            "fill":     fill_c,
        }
    )

# ── Molecule rendering helpers ────────────────────────────────────────────────
MOL_W, MOL_H = 300, 240  # pixels for rdkit rendering


def _prepare_mol(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES: {smiles!r}")
    rdDepictor.Compute2DCoords(mol)
    return mol


def mol_to_svg(smiles: str) -> str:
    mol = _prepare_mol(smiles)
    d = rdMolDraw2D.MolDraw2DSVG(MOL_W, MOL_H)
    opts = d.drawOptions()
    opts.padding = 0.12
    opts.bondLineWidth = 1.5
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def mol_to_pil(smiles: str):
    mol = _prepare_mol(smiles)
    return Draw.MolToImage(mol, size=(MOL_W, MOL_H))


# ── Save individual molecule SVGs ─────────────────────────────────────────────
print("\nSaving individual molecule SVGs …")
for ex in examples:
    svg = mol_to_svg(ex["mol_data"]["smiles"])
    out = IMGS_DIR / f"branching_mol_{ex['key']}.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"  ✓  {out.name}")

# ── Build combined figure ─────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 7.5), facecolor="white")
gs = gridspec.GridSpec(
    2, 3,
    figure=fig,
    height_ratios=[1.4, 1.0],
    hspace=0.10,
    wspace=0.10,
    left=0.07, right=0.97,
    top=0.93,  bottom=0.04,
)

# ── Distribution (top, full width) ───────────────────────────────────────────
ax_dist = fig.add_subplot(gs[0, :])

sns.histplot(
    bi,
    ax=ax_dist,
    bins=26,
    stat="density",
    color=TINTS["blue"][3],
    edgecolor="white",
    linewidth=0.3,
    alpha=0.75,
)
sns.kdeplot(
    bi,
    ax=ax_dist,
    color=TINTS["blue"][0],
    linewidth=2.4,
)

for ex in examples:
    ax_dist.axvline(
        ex["bi"],
        color=ex["color"],
        linewidth=1.8,
        linestyle="--",
        alpha=0.85,
        zorder=5,
    )

ax_dist.set_xlabel("Branching Index", labelpad=5, fontsize=12)
ax_dist.set_ylabel("Density",         labelpad=5, fontsize=12)
ax_dist.set_title(
    f"Distribution of Branching Index — {len(mols)} PFAS molecules",
    fontsize=13, fontweight="bold", pad=8,
)
sns.despine(ax=ax_dist, top=True, right=True)

# ── Molecule panels (bottom row) ──────────────────────────────────────────────
ax_mols: list = []
for col, ex in enumerate(examples):
    ax = fig.add_subplot(gs[1, col])

    try:
        img_arr = np.array(mol_to_pil(ex["mol_data"]["smiles"]))
        ax.imshow(img_arr)
    except Exception as err:  # noqa: BLE001
        ax.text(0.5, 0.5, "Render error", transform=ax.transAxes,
                ha="center", va="center", color="red")
        print(f"  Warning [{ex['key']}]: {err}")

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(ex["color"])
        spine.set_linewidth(2.5)

    pfas_group = ex["mol_data"].get("PFASGroup", "")
    ax.set_title(
        f"{ex['label']}\n{pfas_group}  —  BI = {ex['bi']:.3f}",
        fontsize=10, color=ex["color"], fontweight="bold", pad=5,
    )
    ax_mols.append(ax)

# ── Connecting dotted lines ───────────────────────────────────────────────────
fig.canvas.draw()  # finalise all positions before reading transforms

for ex, ax_mol in zip(examples, ax_mols):
    # Source: bottom of the vertical dashed marker in the distribution axes
    y_bottom = ax_dist.get_ylim()[0]
    disp_src = ax_dist.transData.transform((ex["bi"], y_bottom))
    fig_src  = fig.transFigure.inverted().transform(disp_src)

    # Target: top-centre of the molecule subplot (figure fractions)
    pos     = ax_mol.get_position()
    fig_tgt = np.array([(pos.x0 + pos.x1) / 2.0, pos.y1])

    fig.add_artist(
        plt.Line2D(
            [fig_src[0], fig_tgt[0]],
            [fig_src[1], fig_tgt[1]],
            transform=fig.transFigure,
            color=ex["color"],
            linewidth=1.2,
            linestyle=":",
            alpha=0.65,
            zorder=10,
            clip_on=False,
        )
    )

# ── Save ──────────────────────────────────────────────────────────────────────
out = IMGS_DIR / "branching_index_distribution.svg"
fig.savefig(out, format="svg", bbox_inches="tight")
plt.close(fig)
print(f"\n✅  Saved combined figure: {out}")

# ── PFASGroups branching variant ──────────────────────────────────────────────

def _get_pfasgroups_branching(mol_data: dict) -> "float | None":
    """Return mean_branching from the largest-component group match.

    Tries stored PFASGroups_result first (valid when benchmark was run with
    the full _SUMMARY_KEYS serialisation); falls back to a live parse_mol call.
    """
    raw = mol_data.get("PFASGroups_result")
    if raw:
        try:
            result = ast.literal_eval(raw)
            group_matches = [
                m for m in result.get("matches", [])
                if m.get("type") == "group" and "mean_branching" in m
            ]
            if group_matches:
                best = max(group_matches, key=lambda m: max(m.get("components_sizes") or [0]))
                return float(best["mean_branching"])
        except Exception:
            pass
    # Fall back: live computation
    if not _PFASGROUPS_AVAILABLE:
        return None
    smiles = mol_data.get("smiles")
    if not smiles:
        return None
    try:
        result = _parse_mol(smiles)
        group_matches = [
            m for m in result.get("matches", [])
            if m.get("type") == "group" and m.get("mean_branching") is not None
        ]
        if not group_matches:
            return None
        best = max(group_matches, key=lambda m: max(m.get("components_sizes") or [0]))
        return float(best["mean_branching"])
    except Exception:
        return None


print("\nComputing PFASGroups branching (may call parse_mol for stale records) …")
_pfas_pairs: list[tuple[dict, float]] = []
for _i, _m in enumerate(mols, 1):
    _val = _get_pfasgroups_branching(_m)
    if _val is not None:
        _pfas_pairs.append((_m, _val))
    if _i % 50 == 0:
        print(f"  {_i}/{len(mols)} done")

pfas_mols  = [p[0] for p in _pfas_pairs]
pfas_bi    = np.array([p[1] for p in _pfas_pairs])
print(f"PFASGroups branching computed for {len(pfas_mols)} molecules  "
      f"|  range [{pfas_bi.min():.4f}, {pfas_bi.max():.4f}]")

if len(pfas_mols) >= 3:
    # After fixing calculate_branching to count only C–C bonds, mean_branching
    # has the same direction as branching_index: 1.0 = linear, 0.0 = highly branched.
    # So 10th percentile = most branched = "High branching", same as the first figure.
    _pfas_pcts = [10, 50, 90]
    _pfas_used: set[int] = set()
    pfas_examples: list[dict] = []
    for _p, _label, _key, (_main_c, _fill_c) in zip(_pfas_pcts, labels, keys, EX_PALETTE):
        _idx = _pick_in_window(pfas_bi, pfas_mols, _p, window=10, used=_pfas_used)
        _pfas_used.add(_idx)
        pfas_examples.append(
            {
                "key":      _key,
                "label":    _label,
                "bi":       float(pfas_bi[_idx]),
                "mol_data": pfas_mols[_idx],
                "color":    _main_c,
                "fill":     _fill_c,
            }
        )

    fig2 = plt.figure(figsize=(12, 7.5), facecolor="white")
    gs2  = gridspec.GridSpec(
        2, 3,
        figure=fig2,
        height_ratios=[1.4, 1.0],
        hspace=0.10, wspace=0.10,
        left=0.07, right=0.97, top=0.93, bottom=0.04,
    )

    ax_dist2 = fig2.add_subplot(gs2[0, :])
    sns.histplot(
        pfas_bi,
        ax=ax_dist2,
        bins=26,
        stat="density",
        color=TINTS["dark_purple"][3],
        edgecolor="white",
        linewidth=0.3,
        alpha=0.75,
    )
    sns.kdeplot(pfas_bi, ax=ax_dist2, color=TINTS["dark_purple"][0], linewidth=2.4)
    for _ex in pfas_examples:
        ax_dist2.axvline(
            _ex["bi"], color=_ex["color"], linewidth=1.8, linestyle="--", alpha=0.85, zorder=5
        )
    ax_dist2.set_xlabel("PFASGroups Mean Branching (largest component, 1.0 = linear)", labelpad=5, fontsize=12)
    ax_dist2.set_ylabel("Density", labelpad=5, fontsize=12)
    ax_dist2.set_title(
        f"Distribution of PFASGroups Branching — {len(pfas_mols)} PFAS molecules",
        fontsize=13, fontweight="bold", pad=8,
    )
    sns.despine(ax=ax_dist2, top=True, right=True)

    ax_mols2: list = []
    for _col, _ex in enumerate(pfas_examples):
        _ax = fig2.add_subplot(gs2[1, _col])
        try:
            _img_arr = np.array(mol_to_pil(_ex["mol_data"]["smiles"]))
            _ax.imshow(_img_arr)
        except Exception as _err:
            _ax.text(0.5, 0.5, "Render error", transform=_ax.transAxes,
                     ha="center", va="center", color="red")
            print(f"  Warning [{_ex['key']}]: {_err}")
        _ax.set_xticks([])
        _ax.set_yticks([])
        for _spine in _ax.spines.values():
            _spine.set_visible(True)
            _spine.set_color(_ex["color"])
            _spine.set_linewidth(2.5)
        _pfas_group = _ex["mol_data"].get("PFASGroup", "")
        _ax.set_title(
            f"{_ex['label']}\n{_pfas_group}  —  linearity = {_ex['bi']:.3f}",
            fontsize=10, color=_ex["color"], fontweight="bold", pad=5,
        )
        ax_mols2.append(_ax)

    fig2.canvas.draw()
    for _ex, _ax_mol in zip(pfas_examples, ax_mols2):
        _y_bot   = ax_dist2.get_ylim()[0]
        _d_src   = ax_dist2.transData.transform((_ex["bi"], _y_bot))
        _f_src   = fig2.transFigure.inverted().transform(_d_src)
        _pos     = _ax_mol.get_position()
        _f_tgt   = np.array([(_pos.x0 + _pos.x1) / 2.0, _pos.y1])
        fig2.add_artist(
            plt.Line2D(
                [_f_src[0], _f_tgt[0]], [_f_src[1], _f_tgt[1]],
                transform=fig2.transFigure,
                color=_ex["color"], linewidth=1.2, linestyle=":",
                alpha=0.65, zorder=10, clip_on=False,
            )
        )

    out2 = IMGS_DIR / "branching_index_distribution_pfasgroups.svg"
    fig2.savefig(out2, format="svg", bbox_inches="tight")
    plt.close(fig2)
    print(f"\n✅  Saved PFASGroups variant: {out2}")
else:
    print("Not enough molecules with PFASGroups branching — skipping variant plot")
