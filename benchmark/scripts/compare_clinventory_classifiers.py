#!/usr/bin/env python3
"""
Compare PFASGroups vs PFAS-Atlas clinventory classification results.

Loads the two JSON files produced by classify_PFASGroups_clinventory.py and
classify_pfasatlas_clinventory.py, then:
  • Generates timing comparison plots (overall, by atom-count bracket, by halogen class)
  • Generates PFAS-detection agreement plots
  • Saves a unified comparison JSON to data/
  • Saves all figures to imgs/

Usage (from the benchmark/ directory):
    python scripts/compare_clinventory_classifiers.py \\
        [--hg-file  data/PFASGroups_clinventory_*.json] \\
        [--atlas-file data/pfasatlas_clinventory_*.json]
"""

import argparse
import glob
import json
import math
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.ticker import MaxNLocator
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False
    print("WARNING: matplotlib not available — skipping plots")

try:
    import numpy as np
    NUMPY = True
except ImportError:
    NUMPY = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BRACKETS = [
    "tiny (<10)",
    "small (10-19)",
    "medium (20-34)",
    "large (35-59)",
    "very large (>=60)",
    "unknown",
]

HALOGEN_CLASSES = [
    "F only",
    "Mixed (F + other)",
    "Other halogen (no F)",
    "No halogen",
]

COLORS = {
    "PFASGroups": "#1f77b4",   # blue
    "PFAS-Atlas":    "#d62728",   # red
    "PFASGroups":    "#ff7f0e",   # orange  (F-only PFASGroups)
}

def latest_file(pattern: str) -> Optional[Path]:
    files = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)
    return Path(files[0]) if files else None


def stats(vals: list) -> dict:
    if not vals:
        return {}
    n = len(vals)
    s = sorted(vals)
    mu = sum(vals) / n
    var = sum((v - mu) ** 2 for v in vals) / n
    return {
        "n": n, "mean": mu, "std": math.sqrt(var),
        "median": s[n // 2], "min": s[0], "max": s[-1],
        "p25": s[n // 4] if n >= 4 else s[0],
        "p75": s[3 * n // 4] if n >= 4 else s[-1],
        "p95": s[int(0.95 * n)],
    }


def atom_bracket(n_heavy: Optional[int]) -> str:
    if n_heavy is None:
        return "unknown"
    if n_heavy < 10:   return "tiny (<10)"
    if n_heavy < 20:   return "small (10-19)"
    if n_heavy < 35:   return "medium (20-34)"
    if n_heavy < 60:   return "large (35-59)"
    return "very large (>=60)"


def halogen_class(n_f: Optional[int], n_hal: Optional[int]) -> str:
    nf = n_f or 0
    nh = (n_hal or 0) - nf
    if nf > 0 and nh == 0:   return "F only"
    if nf > 0 and nh > 0:    return "Mixed (F + other)"
    if nf == 0 and nh > 0:   return "Other halogen (no F)"
    return "No halogen"


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

STYLE = {
    "figure.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 10,
}

def apply_style():
    if MATPLOTLIB:
        plt.rcParams.update(STYLE)


def save_fig(fig, name: str, imgs_dir: Path):
    imgs_dir.mkdir(exist_ok=True, parents=True)
    for ext in ("pdf", "png"):
        fig.savefig(imgs_dir / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: imgs/{name}.pdf / .png")


# ---------------------------------------------------------------------------
# Individual plot functions
# ---------------------------------------------------------------------------

def plot_timing_box(hg_times: list, atlas_times: list, imgs_dir: Path):
    """Box plot of timing distributions for both tools."""
    if not MATPLOTLIB:
        return
    apply_style()
    fig, ax = plt.subplots(figsize=(6, 4))
    data = []
    labels = []
    colors = []
    if hg_times:
        data.append(hg_times)
        labels.append("PFASGroups")
        colors.append(COLORS["PFASGroups"])
    if atlas_times:
        data.append(atlas_times)
        labels.append("PFAS-Atlas")
        colors.append(COLORS["PFAS-Atlas"])
    bp = ax.boxplot(data, labels=labels, patch_artist=True,
                    showfliers=False, widths=0.5)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel("Time per molecule (ms)")
    ax.set_title("Classification timing — overall distribution")
    ax.yaxis.set_major_locator(MaxNLocator(5))
    save_fig(fig, "clinventory_timing_box", imgs_dir)


def plot_timing_by_bracket(hg_data: dict, atlas_data: dict, imgs_dir: Path):
    """Grouped bar chart of median timing by atom-count bracket."""
    if not MATPLOTLIB:
        return
    apply_style()
    brackets = [b for b in BRACKETS if b in hg_data or b in atlas_data]
    x = range(len(brackets))
    hg_med   = [hg_data.get(b, {}).get("median_ms", 0) for b in brackets]
    atlas_med= [atlas_data.get(b, {}).get("median_ms", 0) for b in brackets]
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 4))
    bars1 = ax.bar([i - width/2 for i in x], hg_med,   width, label="PFASGroups",
                   color=COLORS["PFASGroups"], alpha=0.8)
    bars2 = ax.bar([i + width/2 for i in x], atlas_med, width, label="PFAS-Atlas",
                   color=COLORS["PFAS-Atlas"],    alpha=0.8)

    ax.set_xticks(list(x))
    ax.set_xticklabels(brackets, rotation=25, ha="right")
    ax.set_ylabel("Median time per molecule (ms)")
    ax.set_title("Timing by molecular size (heavy atom count)")
    ax.legend()
    save_fig(fig, "clinventory_timing_by_bracket", imgs_dir)


def plot_timing_by_halogen(hg_data: dict, atlas_data: dict, imgs_dir: Path):
    """Grouped bar chart of median timing by halogen class."""
    if not MATPLOTLIB:
        return
    apply_style()
    classes = [c for c in HALOGEN_CLASSES if c in hg_data or c in atlas_data]
    x = range(len(classes))
    hg_med    = [hg_data.get(c, {}).get("median_ms", 0) for c in classes]
    atlas_med = [atlas_data.get(c, {}).get("median_ms", 0) for c in classes]
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar([i - width/2 for i in x], hg_med,    width, label="PFASGroups",
           color=COLORS["PFASGroups"], alpha=0.8)
    ax.bar([i + width/2 for i in x], atlas_med, width, label="PFAS-Atlas",
           color=COLORS["PFAS-Atlas"],    alpha=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(classes, rotation=20, ha="right")
    ax.set_ylabel("Median time per molecule (ms)")
    ax.set_title("Timing by halogen composition")
    ax.legend()
    save_fig(fig, "clinventory_timing_by_halogen", imgs_dir)


def plot_timing_ratio_by_bracket(hg_data: dict, atlas_data: dict, imgs_dir: Path):
    """Bar chart of PFASGroups/PFAS-Atlas timing ratio per bracket."""
    if not MATPLOTLIB:
        return
    apply_style()
    brackets = [b for b in BRACKETS if b in hg_data and b in atlas_data]
    ratios = []
    for b in brackets:
        hg_m = hg_data[b].get("median_ms", 0)
        at_m = atlas_data[b].get("median_ms", 0)
        ratios.append(hg_m / at_m if at_m > 0 else float("nan"))

    fig, ax = plt.subplots(figsize=(9, 4))
    colors_r = [COLORS["PFASGroups"] if r >= 1 else COLORS["PFAS-Atlas"] for r in ratios]
    ax.bar(range(len(brackets)), ratios, color=colors_r, alpha=0.8)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(brackets)))
    ax.set_xticklabels(brackets, rotation=25, ha="right")
    ax.set_ylabel("Timing ratio  (PFASGroups / PFAS-Atlas)")
    ax.set_title("Relative timing by molecular size (ratio > 1 → PFASGroups slower)")
    hg_patch  = mpatches.Patch(color=COLORS["PFASGroups"], label="PFASGroups faster")
    at_patch  = mpatches.Patch(color=COLORS["PFAS-Atlas"], label="PFAS-Atlas faster")
    ax.legend(handles=[hg_patch, at_patch])
    save_fig(fig, "clinventory_timing_ratio_by_bracket", imgs_dir)


def plot_classification_agreement(matrix: dict, imgs_dir: Path):
    """2×2 agreement bar chart: both PFAS / only HG / only Atlas / neither."""
    if not MATPLOTLIB:
        return
    apply_style()
    keys   = ["both_pfas", "only_hg", "only_atlas", "neither"]
    labels = ["Both classify\nas PFAS",
              "PFASGroups\nonly",
              "PFAS-Atlas\nonly",
              "Neither"]
    colors = ["#2ca02c", COLORS["PFASGroups"], COLORS["PFAS-Atlas"], "#7f7f7f"]
    values = [matrix.get(k, 0) for k in keys]
    total  = sum(values)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Counts
    ax = axes[0]
    bars = ax.bar(range(4), values, color=colors, alpha=0.8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Number of molecules")
    ax.set_title("Classification agreement — counts")
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.01,
                f"{val:,}", ha="center", va="bottom", fontsize=8)

    # Percentages
    ax = axes[1]
    pcts = [100 * v / total if total > 0 else 0 for v in values]
    bars = ax.bar(range(4), pcts, color=colors, alpha=0.8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Percentage of molecules (%)")
    ax.set_title("Classification agreement — percentages")
    for bar, pct in zip(bars, pcts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.01,
                f"{pct:.1f}%", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    save_fig(fig, "clinventory_classification_agreement", imgs_dir)


def plot_atlas_class_distribution(class_dist: dict, imgs_dir: Path):
    """Horizontal bar chart of PFAS-Atlas class distribution."""
    if not MATPLOTLIB:
        return
    apply_style()
    items = sorted(class_dist.items(), key=lambda x: -x[1])[:20]
    labels = [k for k, _ in items]
    values = [v for _, v in items]

    fig, ax = plt.subplots(figsize=(8, max(4, len(labels) * 0.45)))
    colors_c = [COLORS["PFAS-Atlas"] if k != "Not PFAS" else "#7f7f7f" for k, _ in items]
    bars = ax.barh(range(len(labels)), values, color=colors_c, alpha=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Number of molecules")
    ax.set_title("PFAS-Atlas classification distribution")
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() * 1.005, bar.get_y() + bar.get_height() / 2,
                f"{val:,}", va="center", fontsize=7)
    plt.tight_layout()
    save_fig(fig, "clinventory_atlas_class_distribution", imgs_dir)


def plot_hg_groups_distribution(group_name_counts: dict, imgs_dir: Path):
    """Horizontal bar chart of top-20 PFASGroups group matches."""
    if not MATPLOTLIB:
        return
    apply_style()
    items = sorted(group_name_counts.items(), key=lambda x: -x[1])[:20]
    labels = [k for k, _ in items]
    values = [v for _, v in items]

    fig, ax = plt.subplots(figsize=(8, max(4, len(labels) * 0.45)))
    ax.barh(range(len(labels)), values, color=COLORS["PFASGroups"], alpha=0.8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Number of molecules matched")
    ax.set_title("PFASGroups — top-20 most frequent groups")
    plt.tight_layout()
    save_fig(fig, "clinventory_hg_top_groups", imgs_dir)


# ---------------------------------------------------------------------------
# Three-way plots (PFASGroups / PFASGroups / PFAS-Atlas)
# ---------------------------------------------------------------------------

def plot_timing_box_three(hg_times: list, pfg_times: list, atlas_times: list,
                          imgs_dir: Path):
    """Box plot of timing distributions for all three tools."""
    if not MATPLOTLIB:
        return
    apply_style()
    data, labels, colors = [], [], []
    for times, label, color in [
        (hg_times,   "PFASGroups",    COLORS["PFASGroups"]),
        (pfg_times,  "PFASGroups (F)",   COLORS["PFASGroups"]),
        (atlas_times, "PFAS-Atlas",      COLORS["PFAS-Atlas"]),
    ]:
        if times:
            data.append(times)
            labels.append(label)
            colors.append(color)
    if len(data) < 2:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    bp = ax.boxplot(data, labels=labels, patch_artist=True,
                    showfliers=False, widths=0.5)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel("Time per molecule (ms)")
    ax.set_title("Classification timing — three-tool comparison")
    save_fig(fig, "clinventory_timing_box_three", imgs_dir)


def plot_timing_cdf_three(hg_times: list, pfg_times: list, atlas_times: list,
                          imgs_dir: Path):
    """CDF of per-molecule timing for all three tools (log scale)."""
    if not MATPLOTLIB or not NUMPY:
        return
    apply_style()
    fig, ax = plt.subplots(figsize=(7, 4))
    for times, label, color in [
        (hg_times,   "PFASGroups",   COLORS["PFASGroups"]),
        (pfg_times,  "PFASGroups (F)",  COLORS["PFASGroups"]),
        (atlas_times, "PFAS-Atlas",     COLORS["PFAS-Atlas"]),
    ]:
        if not times:
            continue
        t_sorted = np.sort(times)
        cdf = np.arange(1, len(t_sorted) + 1) / len(t_sorted)
        ax.plot(t_sorted, cdf, label=label, color=color, lw=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("Time per molecule (ms, log scale)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("CDF of per-molecule classification time — three tools")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    save_fig(fig, "clinventory_timing_cdf_three", imgs_dir)


def plot_timing_by_bracket_three(hg_data: dict, pfg_data: dict, atlas_data: dict,
                                  imgs_dir: Path):
    """Grouped bar chart of median timing by atom-count bracket for three tools."""
    if not MATPLOTLIB:
        return
    apply_style()
    brackets = [b for b in BRACKETS
                if b in hg_data or b in pfg_data or b in atlas_data]
    x = range(len(brackets))
    hg_med   = [hg_data.get(b,   {}).get("median_ms", 0) for b in brackets]
    pfg_med  = [pfg_data.get(b,  {}).get("median_ms", 0) for b in brackets]
    atlas_med= [atlas_data.get(b, {}).get("median_ms", 0) for b in brackets]
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar([i - width for i in x],     hg_med,    width, label="PFASGroups",
           color=COLORS["PFASGroups"], alpha=0.8)
    ax.bar([i           for i in x],   pfg_med,   width, label="PFASGroups (F)",
           color=COLORS["PFASGroups"],   alpha=0.8)
    ax.bar([i + width for i in x],     atlas_med, width, label="PFAS-Atlas",
           color=COLORS["PFAS-Atlas"],   alpha=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(brackets, rotation=25, ha="right")
    ax.set_ylabel("Median time per molecule (ms)")
    ax.set_title("Timing by molecular size — three tools")
    ax.legend()
    save_fig(fig, "clinventory_timing_by_bracket_three", imgs_dir)


def plot_f_only_speedup(hg_data: dict, pfg_data: dict, imgs_dir: Path):
    """Bar chart showing per-bracket speed gain when restricting to halogens='F'.

    Ratio = PFASGroups median / PFASGroups median; values > 1 mean
    PFASGroups (F-only) is faster than full PFASGroups.
    """
    if not MATPLOTLIB:
        return
    apply_style()
    brackets = [b for b in BRACKETS if b in hg_data and b in pfg_data]
    ratios = []
    for b in brackets:
        hg_m  = hg_data[b].get("median_ms", 0)
        pfg_m = pfg_data[b].get("median_ms", 0)
        ratios.append(hg_m / pfg_m if pfg_m > 0 else float("nan"))

    fig, ax = plt.subplots(figsize=(9, 4))
    bar_colors = ["#2ca02c" if r >= 1 else "#d62728" for r in ratios]
    bars = ax.bar(range(len(brackets)), ratios, color=bar_colors, alpha=0.85)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(brackets)))
    ax.set_xticklabels(brackets, rotation=25, ha="right")
    ax.set_ylabel("Speed-up factor  (PFASGroups / PFASGroups)")
    ax.set_title(
        "Speed gain from restricting to F-only (halogens=\'F\')\n"
        "Values > 1 mean PFASGroups is faster"
    )
    for bar, ratio in zip(bars, ratios):
        if not math.isnan(ratio):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() * 1.01,
                f"{ratio:.2f}×",
                ha="center", va="bottom", fontsize=8,
            )
    import matplotlib.patches as _mp
    ax.legend(
        handles=[
            _mp.Patch(color="#2ca02c", label="PFASGroups faster (ratio>1)"),
            _mp.Patch(color="#d62728", label="PFASGroups faster (ratio<1)"),
        ]
    )
    save_fig(fig, "clinventory_f_only_speedup", imgs_dir)


def plot_three_way_agreement(hg_map: dict, pfg_map: dict, atlas_map: dict,
                              imgs_dir: Path):
    """Side-by-side agreement matrices: HG vs Atlas and PFASGroups vs Atlas."""
    if not MATPLOTLIB:
        return
    apply_style()

    def _matrix(tool_map, tool_key):
        common = set(tool_map) & set(atlas_map)
        tp = fp = fn = tn = 0
        for mid in common:
            t_f = tool_map[mid].get("has_fluorine_group", False)
            a_f = atlas_map[mid].get("is_pfas", False)
            if t_f and a_f:     tp += 1
            elif t_f and not a_f: fp += 1
            elif not t_f and a_f: fn += 1
            else:               tn += 1
        n = tp + fp + fn + tn
        kappa = _cohen_kappa(tp, fp, fn, tn)
        return {
            "both_pfas": tp, "only_tool": fp,
            "only_atlas": fn, "neither": tn,
            "total": n, "kappa": kappa,
            "agreement_pct": 100.0 * (tp + tn) / n if n else 0,
            "label": tool_key,
        }

    matrices = [("PFASGroups", _matrix(hg_map,  "PFASGroups")),
                ("PFASGroups (F)", _matrix(pfg_map, "PFASGroups"))]

    keys   = ["both_pfas", "only_tool", "only_atlas", "neither"]
    xlabel = ["Both PFAS", "Tool only", "Atlas only", "Neither"]
    tool_colors = {"PFASGroups": COLORS["PFASGroups"],
                   "PFASGroups (F)": COLORS["PFASGroups"]}

    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    for ax, (name, mat) in zip(axes, matrices):
        values = [mat[k] for k in keys]
        total  = mat["total"]
        colors_b = ["#2ca02c", tool_colors[name], COLORS["PFAS-Atlas"], "#7f7f7f"]
        bars = ax.bar(range(4), values, color=colors_b, alpha=0.8)
        ax.set_xticks(range(4))
        ax.set_xticklabels(xlabel, fontsize=9)
        ax.set_ylabel("Number of molecules")
        ax.set_title(
            f"{name} vs PFAS-Atlas agreement\n"
            f"κ = {mat['kappa']:.3f}  |  agree {mat['agreement_pct']:.1f}%"
        )
        for bar, val in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() * 1.01,
                f"{val:,}\n({100*val/max(total,1):.1f}%)",
                ha="center", va="bottom", fontsize=7,
            )
    plt.tight_layout()
    save_fig(fig, "clinventory_three_way_agreement", imgs_dir)


def plot_timing_cdf(hg_times: list, atlas_times: list, imgs_dir: Path):
    """CDF of per-molecule timing for both tools on log scale."""
    if not MATPLOTLIB or not NUMPY:
        return
    apply_style()
    fig, ax = plt.subplots(figsize=(7, 4))
    for times, label, color in [
        (hg_times, "PFASGroups", COLORS["PFASGroups"]),
        (atlas_times, "PFAS-Atlas", COLORS["PFAS-Atlas"]),
    ]:
        if not times:
            continue
        t_sorted = np.sort(times)
        cdf = np.arange(1, len(t_sorted) + 1) / len(t_sorted)
        ax.plot(t_sorted, cdf, label=label, color=color, lw=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("Time per molecule (ms, log scale)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("CDF of per-molecule classification time")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    save_fig(fig, "clinventory_timing_cdf", imgs_dir)


# ---------------------------------------------------------------------------
# Comparison analysis
# ---------------------------------------------------------------------------

def build_comparison(hg: dict, atlas: dict,
                      pfasgroups: Optional[dict] = None) -> dict:
    """Build merged per-molecule comparison table indexed by mol ID."""
    hg_map    = {r["id"]: r for r in hg.get("molecules", [])    if r["time_ms"] is not None}
    atlas_map = {r["id"]: r for r in atlas.get("molecules", []) if r["time_ms"] is not None}
    pfg_map   = (
        {r["id"]: r for r in pfasgroups.get("molecules", []) if r["time_ms"] is not None}
        if pfasgroups else {}
    )
    common_ids = set(hg_map) & set(atlas_map)

    rows = []
    for mid in common_ids:
        hg_r    = hg_map[mid]
        atlas_r = atlas_map[mid]
        rows.append({
            "id":             mid,
            "smiles":         hg_r["smiles"],
            "formula":        hg_r["formula"],
            "n_heavy_atoms":  hg_r["n_heavy_atoms"],
            "n_fluorine":     hg_r["n_fluorine"],
            "n_halogens":     hg_r["n_halogens"],
            "hg_time_ms":     hg_r["time_ms"],
            "atlas_time_ms":  atlas_r["time_ms"],
            "hg_classified":  hg_r["classified"],
            "hg_f_group":     hg_r["has_fluorine_group"],
            "hg_n_groups":    hg_r["n_groups"],
            "atlas_is_pfas":  atlas_r["is_pfas"],
            "atlas_class1":   atlas_r["class1"],
            "atlas_class2":   atlas_r["class2"],
        })

    # Agreement matrix
    both_pfas   = sum(1 for r in rows if r["hg_f_group"] and r["atlas_is_pfas"])
    only_hg     = sum(1 for r in rows if r["hg_f_group"] and not r["atlas_is_pfas"])
    only_atlas  = sum(1 for r in rows if not r["hg_f_group"] and r["atlas_is_pfas"])
    neither     = sum(1 for r in rows if not r["hg_f_group"] and not r["atlas_is_pfas"])
    total       = len(rows)

    agreement_matrix = {
        "both_pfas":  both_pfas,
        "only_hg":    only_hg,
        "only_atlas": only_atlas,
        "neither":    neither,
        "total":      total,
        "agreement_pct": 100.0 * (both_pfas + neither) / total if total else 0,
        "cohen_kappa": _cohen_kappa(both_pfas, only_hg, only_atlas, neither),
    }

    # Timing comparison per bracket
    bracket_timing: dict = defaultdict(lambda: {"hg": [], "atlas": []})
    for r in rows:
        b = atom_bracket(r["n_heavy_atoms"])
        bracket_timing[b]["hg"].append(r["hg_time_ms"])
        bracket_timing[b]["atlas"].append(r["atlas_time_ms"])
    timing_by_bracket = {}
    for b, d in bracket_timing.items():
        hg_s    = stats(d["hg"])
        atlas_s = stats(d["atlas"])
        ratio   = (hg_s.get("median", 0) / atlas_s.get("median", 1)
                   if atlas_s.get("median", 0) > 0 else float("nan"))
        timing_by_bracket[b] = {
            "n": len(d["hg"]),
            "hg": hg_s, "atlas": atlas_s, "ratio_hg_over_atlas": ratio
        }

    # Timing comparison per halogen class
    hcls_timing: dict = defaultdict(lambda: {"hg": [], "atlas": []})
    for r in rows:
        c = halogen_class(r["n_fluorine"], r["n_halogens"])
        hcls_timing[c]["hg"].append(r["hg_time_ms"])
        hcls_timing[c]["atlas"].append(r["atlas_time_ms"])
    timing_by_halogen = {}
    for c, d in hcls_timing.items():
        hg_s    = stats(d["hg"])
        atlas_s = stats(d["atlas"])
        ratio   = (hg_s.get("median", 0) / atlas_s.get("median", 1)
                   if atlas_s.get("median", 0) > 0 else float("nan"))
        timing_by_halogen[c] = {
            "n": len(d["hg"]), "hg": hg_s, "atlas": atlas_s,
            "ratio_hg_over_atlas": ratio
        }

    # Overall timing
    hg_times    = [r["hg_time_ms"]    for r in rows]
    atlas_times = [r["atlas_time_ms"] for r in rows]
    timing_overall = {
        "hg":    stats(hg_times),
        "atlas": stats(atlas_times),
        "ratio_hg_over_atlas": (
            stats(hg_times).get("median", 0) / stats(atlas_times).get("median", 1)
            if stats(atlas_times).get("median", 0) > 0 else float("nan")
        ),
    }

    # PFASGroups group name frequency
    group_name_counts: dict = defaultdict(int)
    for hg_r in hg_map.values():
        for name in hg_r.get("group_names", []):
            if name:
                group_name_counts[str(name)] += 1

    # ── PFASGroups (F-only) timing comparison on molecules present in both ──
    pfg_common = set(pfg_map) & set(hg_map)  # F-molecules seen by both HG and PFG
    pfg_timing_overall: dict = {}
    pfg_timing_by_bracket: dict = {}
    pfg_speedup_overall: float = float("nan")
    if pfg_map:
        hg_paired   = [hg_map[mid]["time_ms"]  for mid in pfg_common]
        pfg_paired  = [pfg_map[mid]["time_ms"] for mid in pfg_common]
        hg_pfg_s    = stats(hg_paired)
        pfg_s       = stats(pfg_paired)
        pfg_timing_overall = {
            "hg": hg_pfg_s, "pfasgroups": pfg_s,
            "ratio_hg_over_pfasgroups": (
                hg_pfg_s.get("median", 0) / pfg_s.get("median", 1)
                if pfg_s.get("median", 0) > 0 else float("nan")
            ),
        }
        pfg_speedup_overall = pfg_timing_overall["ratio_hg_over_pfasgroups"]

        br_hg: dict  = defaultdict(list)
        br_pfg: dict = defaultdict(list)
        for mid in pfg_common:
            b = atom_bracket(hg_map[mid].get("n_heavy_atoms"))
            br_hg[b].append(hg_map[mid]["time_ms"])
            br_pfg[b].append(pfg_map[mid]["time_ms"])
        for b in set(br_hg) | set(br_pfg):
            hb = stats(br_hg.get(b, []))
            pb = stats(br_pfg.get(b, []))
            ratio = (hb.get("median", 0) / pb.get("median", 1)
                     if pb.get("median", 0) > 0 else float("nan"))
            pfg_timing_by_bracket[b] = {
                "n": len(br_hg.get(b, [])),
                "hg": hb, "pfasgroups": pb,
                "ratio_hg_over_pfasgroups": ratio,
            }

    # PFASGroups vs Atlas agreement
    pfg_atlas_common = set(pfg_map) & set(atlas_map)
    pfg_both   = sum(1 for m in pfg_atlas_common
                     if pfg_map[m]["has_fluorine_group"] and atlas_map[m]["is_pfas"])
    pfg_only   = sum(1 for m in pfg_atlas_common
                     if pfg_map[m]["has_fluorine_group"] and not atlas_map[m]["is_pfas"])
    atlas_only_pfg = sum(1 for m in pfg_atlas_common
                         if not pfg_map[m]["has_fluorine_group"] and atlas_map[m]["is_pfas"])
    pfg_neither = sum(1 for m in pfg_atlas_common
                      if not pfg_map[m]["has_fluorine_group"] and not atlas_map[m]["is_pfas"])
    pfg_total   = len(pfg_atlas_common)
    pfg_agreement_matrix = {
        "both_pfas":  pfg_both, "only_pfasgroups": pfg_only,
        "only_atlas": atlas_only_pfg, "neither": pfg_neither,
        "total": pfg_total,
        "agreement_pct": 100.0 * (pfg_both + pfg_neither) / pfg_total if pfg_total else 0,
        "cohen_kappa": _cohen_kappa(pfg_both, pfg_only, atlas_only_pfg, pfg_neither),
    }

    return {
        "common_molecules": total,
        "hg_only_molecules": len(hg_map) - total,
        "atlas_only_molecules": len(atlas_map) - total,
        "agreement_matrix": agreement_matrix,
        "timing_overall": timing_overall,
        "timing_by_atom_bracket": timing_by_bracket,
        "timing_by_halogen_class": timing_by_halogen,
        "group_name_counts": dict(group_name_counts),
        "atlas_class1_distribution": hg.get("class1_distribution", {}),
        "pfasgroups_timing_overall": pfg_timing_overall,
        "pfasgroups_timing_by_bracket": pfg_timing_by_bracket,
        "pfasgroups_speedup_overall": pfg_speedup_overall,
        "pfasgroups_agreement_matrix": pfg_agreement_matrix if pfg_map else {},
        # raw timing lists for plots
        "_hg_times":  hg_times,
        "_atlas_times": atlas_times,
        "_pfg_map": pfg_map,
        "_hg_map":  hg_map,
        "_atlas_map": atlas_map,
    }


def _cohen_kappa(tp, fp, fn, tn) -> float:
    """Cohen's kappa between HG (F group) and Atlas (PFAS) binary labels."""
    n = tp + fp + fn + tn
    if n == 0:
        return float("nan")
    p_o = (tp + tn) / n
    p_e = ((tp + fp) * (tp + fn) + (fp + tn) * (fn + tn)) / (n * n)
    return (p_o - p_e) / (1 - p_e) if (1 - p_e) > 0 else float("nan")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare PFASGroups / PFASGroups / PFAS-Atlas clinventory results"
    )
    p.add_argument("--hg-file",        type=Path, default=None,
                   help="Path to PFASGroups JSON.  Auto-detected if omitted.")
    p.add_argument("--pfasgroups-file", type=Path, default=None,
                   help="Path to PFASGroups (F-only) JSON.  Auto-detected if omitted.")
    p.add_argument("--atlas-file",     type=Path, default=None,
                   help="Path to PFAS-Atlas JSON.  Auto-detected if omitted.")
    p.add_argument("--output",         type=Path, default=None,
                   help="Output comparison JSON path")
    return p.parse_args()


def main():
    args = parse_args()
    data_dir = BENCHMARK_DIR / "data"
    imgs_dir = BENCHMARK_DIR / "imgs"

    # Load files
    hg_path  = args.hg_file   or latest_file(str(data_dir / "PFASGroups_clinventory_*.json"))
    pfg_path = args.pfasgroups_file or latest_file(str(data_dir / "pfasgroups_clinventory_*.json"))
    atlas_path = args.atlas_file or latest_file(str(data_dir / "pfasatlas_clinventory_*.json"))

    if not hg_path or not Path(hg_path).exists():
        print(f"ERROR: PFASGroups result file not found (searched: {data_dir / 'PFASGroups_clinventory_*.json'})")
        sys.exit(1)
    if not atlas_path or not Path(atlas_path).exists():
        print(f"ERROR: PFAS-Atlas result file not found (searched: {data_dir / 'pfasatlas_clinventory_*.json'})")
        sys.exit(1)

    print("=" * 70)
    print("Clinventory classifier comparison")
    print("=" * 70)
    print(f"  PFASGroups : {hg_path}")
    if pfg_path and Path(pfg_path).exists():
        print(f"  PFASGroups(F) : {pfg_path}")
    else:
        print("  PFASGroups(F) : not found — three-way comparison skipped")
        pfg_path = None
    print(f"  PFAS-Atlas    : {atlas_path}")

    with open(hg_path)    as f: hg    = json.load(f)
    with open(atlas_path) as f: atlas = json.load(f)
    pfasgroups = None
    if pfg_path:
        with open(pfg_path) as f:
            pfasgroups = json.load(f)

    print(f"\n  PFASGroups : {hg['n_molecules_fetched']:,} molecules")
    if pfasgroups:
        print(f"  PFASGroups(F) : {pfasgroups['n_molecules_fetched']:,} molecules")
    print(f"  PFAS-Atlas    : {atlas['n_molecules_fetched']:,} molecules")

    # -----------------------------------------------------------------------
    print("\nBuilding comparison …")
    cmp = build_comparison(hg, atlas, pfasgroups)
    print(f"  Common molecules: {cmp['common_molecules']:,}")

    # Pull out timing lists / maps (not serializable)
    hg_times    = cmp.pop("_hg_times")
    atlas_times = cmp.pop("_atlas_times")
    pfg_map_raw = cmp.pop("_pfg_map", {})
    hg_map_raw  = cmp.pop("_hg_map",  {})
    atlas_map_raw = cmp.pop("_atlas_map", {})
    pfg_times = [pfg_map_raw[mid]["time_ms"]
                 for mid in sorted(pfg_map_raw)
                 if pfg_map_raw[mid]["time_ms"] is not None]

    # -----------------------------------------------------------------------
    print("\nGenerating plots …")

    # Two-way plots (PFASGroups vs PFAS-Atlas) — always produced
    plot_timing_box(hg_times, atlas_times, imgs_dir)
    plot_timing_cdf(hg_times, atlas_times, imgs_dir)
    plot_timing_by_bracket(
        hg["timing_by_atom_bracket"],
        atlas["timing_by_atom_bracket"],
        imgs_dir,
    )
    plot_timing_by_halogen(
        hg["timing_by_halogen_class"],
        atlas["timing_by_halogen_class"],
        imgs_dir,
    )
    plot_timing_ratio_by_bracket(
        hg["timing_by_atom_bracket"],
        atlas["timing_by_atom_bracket"],
        imgs_dir,
    )
    plot_classification_agreement(cmp["agreement_matrix"], imgs_dir)
    plot_atlas_class_distribution(atlas.get("class1_distribution", {}), imgs_dir)
    plot_hg_groups_distribution(cmp["group_name_counts"], imgs_dir)

    # Three-way plots — only when PFASGroups data is available
    if pfasgroups:
        print("  Three-way comparison plots …")
        pfg_bracket = pfasgroups.get("timing_by_atom_bracket", {})
        plot_timing_box_three(hg_times, pfg_times, atlas_times, imgs_dir)
        plot_timing_cdf_three(hg_times, pfg_times, atlas_times, imgs_dir)
        plot_timing_by_bracket_three(
            hg["timing_by_atom_bracket"], pfg_bracket,
            atlas["timing_by_atom_bracket"], imgs_dir,
        )
        plot_f_only_speedup(hg["timing_by_atom_bracket"], pfg_bracket, imgs_dir)
        plot_three_way_agreement(hg_map_raw, pfg_map_raw, atlas_map_raw, imgs_dir)

    # -----------------------------------------------------------------------
    # Save comparison JSON (remove non-serializable nan floats)
    import math as _math
    def _clean(obj):
        if isinstance(obj, float) and _math.isnan(obj):
            return None
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_clean(v) for v in obj]
        return obj

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    data_dir.mkdir(exist_ok=True)
    out_path = args.output or (data_dir / f"clinventory_comparison_{timestamp}.json")
    with open(out_path, "w") as f:
        json.dump(_clean(cmp), f, indent=2)
    print(f"\nComparison data saved to: {out_path}")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)
    am = cmp["agreement_matrix"]
    to = cmp["timing_overall"]
    print(f"\n  Timing (per molecule, ms):")
    print(f"    {'Tool':<20s}  {'Median':>8s}  {'Mean':>8s}  {'Std':>8s}  {'P95':>8s}")
    print(f"    {'-'*56}")
    for tool, key in [("PFASGroups", "hg"), ("PFAS-Atlas", "atlas")]:
        s = to[key]
        if s:
            print(f"    {tool:<20s}  {s['median']:>8.3f}  {s['mean']:>8.3f}"
                  f"  {s['std']:>8.3f}  {s['p95']:>8.3f}")
    try:
        if not math.isnan(to.get("ratio_hg_over_atlas", float("nan"))):
            print(f"\n  Timing ratio (HG/Atlas median): {to['ratio_hg_over_atlas']:.2f}×")
    except Exception:
        pass

    if pfasgroups and cmp.get("pfasgroups_timing_overall"):
        pto = cmp["pfasgroups_timing_overall"]
        hg_med_pfg = pto.get("hg", {}).get("median", float("nan"))
        pfg_med    = pto.get("pfasgroups", {}).get("median", float("nan"))
        speedup    = cmp.get("pfasgroups_speedup_overall", float("nan"))
        print(f"\n  PFASGroups (F-only) on {pfasgroups['n_molecules_fetched']:,} F-molecules:")
        print(f"    PFASGroups median  : {hg_med_pfg:.3f} ms/mol")
        print(f"    PFASGroups    median  : {pfg_med:.3f} ms/mol")
        try:
            if not math.isnan(speedup):
                print(f"    Speed gain (HG/PFG)   : {speedup:.2f}×")
        except Exception:
            pass
        pam = cmp.get("pfasgroups_agreement_matrix", {})
        if pam.get("total"):
            print(f"\n  PFASGroups vs PFAS-Atlas ({pam['total']:,} paired molecules):")
            print(f"    Both classify as PFAS  : {pam['both_pfas']:>7,}")
            print(f"    PFASGroups only        : {pam['only_pfasgroups']:>7,}")
            print(f"    PFAS-Atlas only        : {pam['only_atlas']:>7,}")
            print(f"    Neither                : {pam['neither']:>7,}")
            print(f"    Agreement              : {pam['agreement_pct']:.1f}%")
            print(f"    Cohen's kappa          : {pam['cohen_kappa']:.3f}")

    print(f"\n  Classification agreement ({am['total']:,} paired molecules):")
    print(f"    Both classify as F/PFAS  : {am['both_pfas']:>7,}  ({100*am['both_pfas']/max(am['total'],1):.1f}%)")
    print(f"    PFASGroups only (F)   : {am['only_hg']:>7,}  ({100*am['only_hg']/max(am['total'],1):.1f}%)")
    print(f"    PFAS-Atlas only          : {am['only_atlas']:>7,}  ({100*am['only_atlas']/max(am['total'],1):.1f}%)")
    print(f"    Neither classifies       : {am['neither']:>7,}  ({100*am['neither']/max(am['total'],1):.1f}%)")
    print(f"    Overall agreement        : {am['agreement_pct']:.1f}%")
    print(f"    Cohen's kappa            : {am['cohen_kappa']:.3f}")

    print(f"\nOutput: {out_path}")
    print(f"Plots : imgs/clinventory_*.pdf / .png")


if __name__ == "__main__":
    main()
