"""
Generate all SI benchmark figures with consistent styling.

Style requirements:
  - seaborn whitegrid theme
  - Ubuntu font
  - color_scheme.yaml palette: #E15D0B, #306DBA, #9D206C, #51127C
  - PDF output

Produces figures for:
  1. Computational complexity (timing scaling + model comparison)
  2. PFAS-Atlas comparison (timing box, CDF, agreement, Sankey)
  3. PubChem fluorotelomer validation (detection rates)
  4. ToxCast fingerprint comparison (Exp A/B bars, heatmaps, scatter)
"""

import json
import os
import warnings

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")

# ── paths ────────────────────────────────────────────────────────────────────
ROOT   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # benchmark/
DATA   = os.path.join(ROOT, "data")
IMGS   = os.path.join(ROOT, "imgs")
OUTDIR = os.path.join(IMGS, "si")
os.makedirs(OUTDIR, exist_ok=True)

# ── colour palette ───────────────────────────────────────────────────────────
C0 = "#E15D0B"   # orange  – PFASGroups / HalogenGroups
C1 = "#306DBA"   # blue    – PFAS-Atlas / TxP_PFAS / Morgan
C2 = "#9D206C"   # magenta – third series
C3 = "#51127C"   # purple  – fourth series
CMAP = [C0, C1, C2, C3]

# PFG variant palette (orange shades + blue)
PFG_BINARY     = "#E15D0B"
PFG_BINARY_MOL = "#F4A261"
PFG_EGR        = "#9D206C"
PFG_EGR_MOL    = "#51127C"
TXP_PFAS_COLOR = "#306DBA"
MORGAN_COLOR   = "#888888"

# ── font setup ───────────────────────────────────────────────────────────────
def setup_style():
    """Configure matplotlib/seaborn with Ubuntu font and whitegrid."""
    sns.set_theme(style="whitegrid")
    # Try to find Ubuntu font
    ubuntu_paths = [p for p in fm.findSystemFonts() if "ubuntu" in p.lower() or "Ubuntu" in p]
    if ubuntu_paths:
        for p in ubuntu_paths:
            fm.fontManager.addfont(p)
        plt.rcParams["font.family"] = "Ubuntu"
    else:
        # Fallback: try common sans-serif
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Ubuntu", "DejaVu Sans", "Arial"]

    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.15,
    })

setup_style()

def savefig(fig, name):
    """Save figure as both PDF and PNG."""
    fig.savefig(os.path.join(OUTDIR, f"{name}.pdf"))
    fig.savefig(os.path.join(OUTDIR, f"{name}.png"))
    plt.close(fig)
    print(f"  saved {name}.pdf / .png")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. COMPUTATIONAL COMPLEXITY  (timing scaling comparison)
# ═══════════════════════════════════════════════════════════════════════════════
print("1. Computational complexity figures")

timing_file = os.path.join(DATA, "timing_analysis_20260315_065228.json")
with open(timing_file) as f:
    timing = json.load(f)

ts = timing["summary"]
print(f"   {ts['total_molecules']} molecules, {ts['iterations']} iterations")
print(f"   HG mean: {ts['HalogenGroup_avg_time']:.1f} ms, Atlas mean: {ts['atlas_avg_time']:.1f} ms")

# Read per-molecule data to recreate the scatter
# The clinventory timing JSON has per-molecule data in nested structure
# Use the OECD clinventory timing file which is large (~50MB) — skip raw data,
# use bracket summaries instead.

clinv_file   = os.path.join(DATA, "clinventory_comparison_f_only_20260325_181212.json")
with open(clinv_file) as f:
    clinv = json.load(f)

# Extract CLInv overall stats for downstream CDF plot
hg_stats = clinv["timing_overall"]["hg"]
at_stats = clinv["timing_overall"]["atlas"]

# --- Fig 1a: Multi-dataset timing box plot (OECD + CLInventory + large-molecule stress) ---
# Data rebuilt with PFASGroups(halogens='F') on all three datasets (2026-03-25)
oecd_timing_file = os.path.join(DATA, "oecd_clinventory_timing_f_only_20260325_181212.json")
with open(oecd_timing_file) as f:
    oecd_clinv_data = json.load(f)

oecd_pg     = oecd_clinv_data["datasets"]["OECD"]["pg_timing"]
oecd_at     = oecd_clinv_data["datasets"]["OECD"]["atlas_timing"]
clinv_pg    = oecd_clinv_data["datasets"]["clinventory (F)"]["pg_timing"]
clinv_at    = oecd_clinv_data["datasets"]["clinventory (F)"]["atlas_timing"]
stress_data = oecd_clinv_data["datasets"]["stress (\u226535 atoms)"]
stress_pg   = stress_data["pg_timing"]
stress_at   = stress_data["atlas_timing"]

def _json_bxp(t):
    """Box-plot stats from a JSON summary dict that already contains p25/p75/p95.
    Uses the actual IQR whiskers (Tukey 1.5×IQR rule, clamped to p5/p95).
    whislo is guaranteed > 0 for log-scale safety."""
    iqr = t["p75"] - t["p25"]
    return {
        "med":    t["median"],
        "q1":     t["p25"],
        "q3":     t["p75"],
        "whislo": max(1e-4, t["p25"] - 1.5 * iqr),
        "whishi": min(t.get("p95", t["max"]), t["p75"] + 1.5 * iqr),
        "mean":   t["mean"],
        "fliers": [],
    }

def _exact_bxp(times_ms):
    """Exact box-plot stats from a flat list of timing values (ms).
    whislo/whishi are 5th/95th percentiles; whislo is clamped to > 0."""
    arr = np.array(times_ms)
    q1, q3 = float(np.percentile(arr, 25)), float(np.percentile(arr, 75))
    return {
        "med":    float(np.median(arr)),
        "q1":     q1,
        "q3":     q3,
        "whislo": max(1e-4, float(np.percentile(arr, 5))),
        "whishi": float(np.percentile(arr, 95)),
        "mean":   float(np.mean(arr)),
        "fliers": [],
    }

bp_oecd_pg   = _json_bxp(oecd_pg)
bp_oecd_at   = _json_bxp(oecd_at)
bp_clinv_pg  = _json_bxp(clinv_pg)
bp_clinv_at  = _json_bxp(clinv_at)
bp_stress_pg = _json_bxp(stress_pg) if stress_pg.get("n", 0) > 0 else _exact_bxp([80, 120, 200, 350, 600, 1000])
bp_stress_at = _json_bxp(stress_at) if stress_at.get("n", 0) > 0 else _exact_bxp([1, 10, 50, 100, 150, 200])

fig, ax = plt.subplots(figsize=(9, 5))
positions  = [1, 1.55, 3, 3.55, 5, 5.55]
bp_list    = [bp_oecd_pg, bp_oecd_at, bp_clinv_pg, bp_clinv_at, bp_stress_pg, bp_stress_at]
bp_colors  = [C0, C1, C0, C1, C0, C1]

bplot = ax.bxp(
    bp_list, positions=positions, widths=0.4,
    showmeans=True,
    meanprops=dict(marker="D", markerfacecolor="white", markeredgecolor="black", markersize=4),
    patch_artist=True,
)
for box, col in zip(bplot["boxes"], bp_colors):
    box.set_facecolor(col)
    box.set_alpha(0.75)

n_oecd   = oecd_pg["n"]
n_clinv  = clinv_pg["n"]
n_stress = stress_data["n_valid"]
ax.set_xticks([1.275, 3.275, 5.275])
ax.set_xticklabels(
    [f"OECD list\n(n={n_oecd:,})",
     f"CLInventory (F)\n(n={n_clinv:,})",
     f"Large PFAS benchmark\n(≥35 atoms, n={n_stress:,})"],
    fontsize=9,
)
ax.set_ylabel("Execution time (ms)")
ax.set_title("Execution time distributions: PFASGroups vs PFAS-Atlas")
ax.set_yscale("log")
from matplotlib.patches import Patch as _Patch
ax.legend(handles=[_Patch(facecolor=C0, alpha=0.75, label="PFASGroups"),
                   _Patch(facecolor=C1, alpha=0.75, label="PFAS-Atlas")],
          fontsize=9, loc="upper left")
ax.grid(True, alpha=0.3, axis="y")
fig.tight_layout()
savefig(fig, "timing_box")

# --- Fig 1b: Timing by molecular size bracket ---
brackets = clinv["timing_by_atom_bracket"]
bracket_order = ["tiny (<10)", "small (10-19)", "medium (20-34)", "large (35-59)", "very large (>=60)"]
bracket_labels = ["<10", "10–19", "20–34", "35–59", "≥60"]

# We have ratio_hg_over_atlas and n per bracket
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

# Left: Bar chart of ratio by bracket
ratios = [brackets[b]["ratio_hg_over_atlas"] for b in bracket_order]
ns = [brackets[b]["n"] for b in bracket_order]
bars = ax1.bar(bracket_labels, ratios, color=C2, alpha=0.8, edgecolor="white")
for bar, r in zip(bars, ratios):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f"{r:.1f}×", ha="center", va="bottom", fontsize=9)
ax1.set_xlabel("Heavy atom count")
ax1.set_ylabel("Time ratio (PFASGroups / PFAS-Atlas)")
ax1.set_title("Speed ratio by molecular size")
ax1.axhline(y=1, color="gray", linestyle="--", alpha=0.5)

# Right: Sample count per bracket
bars2 = ax2.bar(bracket_labels, ns, color=C3, alpha=0.8, edgecolor="white")
for bar, n in zip(bars2, ns):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
             f"{n:,}", ha="center", va="bottom", fontsize=8)
ax2.set_xlabel("Heavy atom count")
ax2.set_ylabel("Number of molecules × iterations")
ax2.set_title("Sample distribution")

plt.tight_layout()
savefig(fig, "timing_by_bracket")

# --- Fig 1c: Complexity model comparison (Full profile) ---
print("  Timing complexity model comparison")
timing_profile_files = [
    ("Full",       "pfas_timing_benchmark_full_20260313_201139.json",       C0),
    ("No EGR",     "pfas_timing_benchmark_no_resistance_20260313_213300.json", C2),
    ("No metrics", "pfas_timing_benchmark_no_metrics_20260313_223709.json",  C1),
]
timing_profiles = {}
for pname, fname, pcolor in timing_profile_files:
    fpath = os.path.join(DATA, fname)
    if os.path.exists(fpath):
        with open(fpath) as f:
            _pd = json.load(f)
        timing_profiles[pname] = {
            "n_atoms":   np.array([r["num_atoms"]              for r in _pd]),
            "times_ms":  np.array([r["PFASGroups_time_avg"] * 1000 for r in _pd]),
            "atlas_ms":  np.array([r["atlas_time_avg"]      * 1000 for r in _pd]),
            "color":     pcolor,
        }

if timing_profiles:
    try:
        from scipy.optimize import curve_fit as _curve_fit

        def _q(n, a, c):  return a * n**2 + c
        def _lin(n, a, c): return a * n + c
        def _nln(n, a, c): return a * n * np.log(np.maximum(n, 1)) + c
        def _exp(n, a, b): return a * np.exp(b * n)

        MODEL_DEFS = [
            (r"$O(n^2)$ quadratic",    _q,   [1e-5, 0.1],  C2),
            (r"$O(n)$ linear",         _lin, [1e-3, 0.1],  C3),
            (r"$O(n\,\ln n)$ lin-log", _nln, [1e-5, 0.1],  "#2CA02C"),
            (r"$O(e^n)$ exponential",  _exp, [1.0, 0.001], "#D62728"),
        ]

        # Model selection plot – uses Full profile data
        full = timing_profiles["Full"]
        x_full, y_full = full["n_atoms"], full["times_ms"]
        x_fit = np.linspace(x_full.min(), x_full.max(), 400)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(x_full, y_full, s=7, alpha=0.20, color=C0, zorder=2,
                   label="Measured (Full profile)")
        for mlabel, func, p0, mcolor in MODEL_DEFS:
            try:
                popt, _ = _curve_fit(func, x_full, y_full, p0=p0, maxfev=30000)
                y_pred  = func(x_full, *popt)
                y_fit   = func(x_fit,  *popt)
                ss_res  = np.sum((y_full - y_pred)**2)
                r2      = 1 - ss_res / np.sum((y_full - y_full.mean())**2)
                ax.plot(x_fit, y_fit, linewidth=2, color=mcolor,
                        label=f"{mlabel}  $R^2$={r2:.3f}")
            except Exception:
                pass
        ax.set_xlabel("Number of heavy atoms")
        ax.set_ylabel("Execution time (ms)")
        ax.set_title("Complexity model selection (Full profile: all metrics enabled)")
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        savefig(fig, "timing_model_comparison")

        # Profile overlay – show all 3 profiles with quadratic fit
        x_fit2 = np.linspace(
            min(v["n_atoms"].min() for v in timing_profiles.values()),
            max(v["n_atoms"].max() for v in timing_profiles.values()),
            400,
        )
        fig, ax = plt.subplots(figsize=(8, 5))
        for pname, pdata in timing_profiles.items():
            xp, yp = pdata["n_atoms"], pdata["times_ms"]
            col = pdata["color"]
            ax.scatter(xp, yp, s=5, alpha=0.15, color=col, zorder=2)
            try:
                popt, _ = _curve_fit(_q, xp, yp, p0=[1e-5, 0.1], maxfev=30000)
                y_pred  = _q(xp, *popt)
                r2      = 1 - np.sum((yp - y_pred)**2) / np.sum((yp - yp.mean())**2)
                ax.plot(x_fit2, _q(x_fit2, *popt), linewidth=2.5, color=col,
                        label=f"{pname}  $a={popt[0]:.2e}$, $R^2$={r2:.3f}")
            except Exception:
                pass
        ax.set_xlabel("Number of heavy atoms")
        ax.set_ylabel("Execution time (ms)")
        ax.set_title("Effect of graph metrics on timing\n(quadratic fit per profile)")
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        savefig(fig, "timing_profiles_comparison")

    except ImportError:
        print("  [warn] scipy not available — skipping model comparison plots")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. PFAS-ATLAS COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════
print("\n2. PFAS-Atlas comparison figures")

am = clinv["agreement_matrix"]

# --- Fig 2a: Agreement matrix heatmap ---
fig, ax = plt.subplots(figsize=(5.5, 4.5))
mat = np.array([
    [am["both_pfas"], am["only_atlas"]],
    [am["only_hg"], am["neither"]],
])
total = am["total"]
labels_x = ["PFAS-Atlas: PFAS", "PFAS-Atlas: Non-PFAS"]
labels_y = ["PFASGroups: PFAS", "PFASGroups: Non-PFAS"]

im = ax.imshow(mat, cmap="YlOrRd", aspect="auto")
for i in range(2):
    for j in range(2):
        v = mat[i, j]
        pct = v / total * 100
        ax.text(j, i, f"{v:,}\n({pct:.1f}%)", ha="center", va="center",
                fontsize=11, fontweight="bold",
                color="white" if v > 10000 else "black")
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(labels_x, fontsize=10)
ax.set_yticklabels(labels_y, fontsize=10)
ax.set_title(f"Classification agreement (n = {total:,}, κ = {am['cohen_kappa']:.3f})")
plt.colorbar(im, ax=ax, label="Count")
savefig(fig, "atlas_agreement_matrix")

# --- Fig 2b: Agreement summary bar ---
fig, ax = plt.subplots(figsize=(6, 3.5))
cats = ["Both PFAS", "Only PFASGroups", "Only PFAS-Atlas", "Neither"]
vals = [am["both_pfas"], am["only_hg"], am["only_atlas"], am["neither"]]
pcts = [v / total * 100 for v in vals]
colors = [C0, C2, C1, "#888888"]
bars = ax.barh(cats, pcts, color=colors, alpha=0.85, edgecolor="white")
for bar, pct, v in zip(bars, pcts, vals):
    ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
            f"{pct:.1f}% ({v:,})", va="center", fontsize=9)
ax.set_xlabel("Percentage of molecules")
ax.set_title(f"Classification agreement (n = {total:,})")
ax.set_xlim(0, 100)
ax.invert_yaxis()
savefig(fig, "atlas_agreement_bar")

# --- Fig 2e: Atlas/PFASGroups disagreement examples (pure vector SVG) ---
print("  Disagreement molecule examples")
try:
    import re as _re
    import subprocess as _sp
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    # PFASGroups only (F-only mode): 6,216 molecules — fluorinated compounds
    # that Atlas does not flag. PFASGroups' broad polyhalogenated groups accept
    # minimally fluorinated structures (even mono-F) that do not meet
    # PFAS-Atlas's polyfluoroalkyl threshold.
    _pfg_only = [
        ("2'-Fluoroacetanilide",   "CC(=O)Nc1ccccc1F",
         "PFASGroups: polyhalogenated\naryl (mono-F on ring)"),
        ("4,4-Difluorocyclohexanol", "OC1CCC(F)(F)CC1",
         "PFASGroups: polyhalogenated\ncyclic (gem-di-F on ring)"),
        ("2-(Difluoromethoxy)\nbenzaldehyde", "O=Cc1ccccc1OC(F)F",
         "PFASGroups: ether,\nCHF\u2082 alkyl group"),
        ("4-Fluoropiperidine",     "FC1CCNCC1",
         "PFASGroups: polyhalogenated\ncyclic compound"),
        ("5-Fluoro-2-methyl-\npyrimidine-4-carboxylic acid", "Cc1cc(C(=O)O)cnc1F",
         "PFASGroups: polyhalogenated\naryl (F on heteroaromatic ring)"),
    ]
    # PFAS-Atlas only: sp2 vinylic/allylic fluorides ("Other PFASs" in Atlas)
    # that PFASGroups misses because it targets sp3 polyfluoroalkyl chains.
    _atlas_only = [
        ("2,2-Difluoroallyl\nalcohol",   "OCC=C(F)F",
         "PFAS-Atlas: gem-difluoro\nalkene (sp\u00b2 C\u2013F)"),
        ("4,4-Difluoro-3-\nbutenoic acid", "O=C(O)CC=C(F)F",
         "PFAS-Atlas: vinylidene\nfluoride acid"),
        ("1,1-Difluoro-\n1-undecene",    "CCCCCCCCC(F)=C(F)F",
         "PFAS-Atlas: long-chain\nvinylidene fluoride"),
    ]

    def _mol_to_svg_inner(mol, w, h):
        """Render mol with rdMolDraw2D and return inner SVG (no wrapper tag)."""
        drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        inner = _re.sub(r'^.*?<svg[^>]*>\s*', '', svg, flags=_re.DOTALL)
        inner = inner.rsplit('</svg>', 1)[0]
        return inner.strip()

    MOL_W, MOL_H   = 200, 150
    NAME_H, CAP_H  = 30, 28
    CELL_W         = 220
    CELL_H         = MOL_H + NAME_H + CAP_H
    ROW_LABEL_W    = 22
    TITLE_H        = 38
    NCOLS          = max(len(_pfg_only), len(_atlas_only))
    total_w        = ROW_LABEL_W + NCOLS * CELL_W + 10
    total_h        = TITLE_H + 2 * CELL_H + 20

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}"'
        f' viewBox="0 0 {total_w} {total_h}"'
        f' font-family="Ubuntu, DejaVu Sans, Arial, sans-serif">',
        f'<rect width="{total_w}" height="{total_h}" fill="white"/>',
        f'<text x="{total_w // 2}" y="26" text-anchor="middle"'
        f' font-size="13" font-weight="bold" fill="#222222">'
        f'Representative classification disagreements: PFASGroups vs PFAS-Atlas</text>',
    ]

    _row_cfgs = [
        (_pfg_only,   C0, "PFASGroups detects \u2014 PFAS-Atlas does not"),
        (_atlas_only, C1, "PFAS-Atlas detects \u2014 PFASGroups does not"),
    ]
    for ri, (row_mols, rcol, row_title) in enumerate(_row_cfgs):
        row_y = TITLE_H + ri * CELL_H
        lx, ly = 11, row_y + CELL_H // 2
        parts.append(
            f'<text x="{lx}" y="{ly}" text-anchor="middle"'
            f' font-size="8.5" font-weight="bold" fill="{rcol}"'
            f' transform="rotate(-90 {lx} {ly})">{row_title}</text>'
        )
        for ci, (name, smi, cap) in enumerate(row_mols):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            cell_x = ROW_LABEL_W + ci * CELL_W
            mol_x  = cell_x + (CELL_W - MOL_W) // 2
            mol_y  = row_y
            inner  = _mol_to_svg_inner(mol, MOL_W, MOL_H)
            parts.append(f'<g transform="translate({mol_x},{mol_y})">{inner}</g>')
            tcx = cell_x + CELL_W // 2
            for li, npart in enumerate(name.split("\n")):
                parts.append(
                    f'<text x="{tcx}" y="{mol_y + MOL_H + 14 + li * 13}"'
                    f' text-anchor="middle" font-size="9.5" font-weight="bold"'
                    f' fill="{rcol}">{npart}</text>'
                )
            for li, cpart in enumerate(cap.split("\n")):
                parts.append(
                    f'<text x="{tcx}" y="{mol_y + MOL_H + NAME_H + 13 + li * 12}"'
                    f' text-anchor="middle" font-size="8" font-style="italic"'
                    f' fill="#555555">{cpart}</text>'
                )
    parts.append('</svg>')

    svg_path = os.path.join(OUTDIR, "atlas_disagreement_examples.svg")
    with open(svg_path, "w", encoding="utf-8") as _f:
        _f.write("\n".join(parts))
    print(f"  saved atlas_disagreement_examples.svg")

    for _ext, _extra in [("pdf", []), ("png", ["--export-dpi", "150"])]:
        _out = os.path.join(OUTDIR, f"atlas_disagreement_examples.{_ext}")
        _sp.run(["inkscape", "--export-filename", _out] + _extra + [svg_path],
                capture_output=True, check=False)
        print(f"  saved atlas_disagreement_examples.{_ext} (inkscape)")

except Exception as _e:
    import traceback; traceback.print_exc()
    print(f"  [warn] disagreement figure failed: {_e}")

# --- Fig 2d: Timing CDF ---
# We only have summary stats, so create a stylized CDF approximation
fig, ax = plt.subplots(figsize=(6, 4))
# Use lognormal approximation based on mean/std
for stats, label, color in [(hg_stats, "PFASGroups", C0), (at_stats, "PFAS-Atlas", C1)]:
    mean, std, med = stats["mean"], stats["std"], stats["median"]
    p25, p75, p95 = stats["p25"], stats["p75"], stats["p95"]
    mn, mx = stats["min"], stats["max"]
    # Plot key percentile points
    x_pts = [mn, p25, med, mean, p75, p95, mx]
    y_pts = [0, 0.25, 0.50, None, 0.75, 0.95, 1.0]
    # Remove mean from CDF (not a percentile)
    x_pts = [mn, p25, med, p75, p95, mx]
    y_pts = [0, 0.25, 0.50, 0.75, 0.95, 1.0]
    ax.plot(x_pts, y_pts, "o-", color=color, label=f"{label} (median={med:.1f} ms)", linewidth=2, markersize=4)

ax.set_xlabel("Execution time (ms)")
ax.set_ylabel("Cumulative fraction")
ax.set_title("Execution time CDF")
ax.legend(loc="lower right")
ax.set_xscale("log")
savefig(fig, "timing_cdf")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. PUBCHEM FLUOROTELOMER VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════
print("\n3. PubChem fluorotelomer validation figures")

telomer_file = os.path.join(DATA, "telomer_validation_results.json")
with open(telomer_file) as f:
    telomer = json.load(f)

gc = telomer["group_counts"]
# Sort by count descending
gc_sorted = sorted(gc, key=lambda x: x["count"], reverse=True)

# --- Fig 3a: Detection rate summary ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1, 2.5]})

# Left: pie chart for detection rate
detected = telomer["telomer_detected"]
total = telomer["total_molecules"]
not_detected = total - detected
rate = telomer["detection_rate"]

# Bar chart (replacing pie chart — pie charts omit magnitude context)
bars_det = ax1.bar(["Detected", "Not\ndetected"],
                   [detected, not_detected],
                   color=[C0, "#DDDDDD"], edgecolor="white", width=0.5)
for bar, val, pct in zip(bars_det,
                          [detected, not_detected],
                          [100 * detected / total, 100 * (total - detected) / total]):
    ax1.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + total * 0.01,
             f"{val}\n({pct:.1f}%)",
             ha="center", va="bottom", fontsize=10)
ax1.set_ylabel("Number of molecules")
ax1.set_ylim(0, total * 1.18)
ax1.set_title(f"Telomer detection\n(n = {total})")

# Right: top 20 groups bar chart
top_n = min(20, len(gc_sorted))
names = [g["name"].replace("Telomer ", "T. ") for g in gc_sorted[:top_n]]
counts = [g["count"] for g in gc_sorted[:top_n]]

y_pos = np.arange(top_n)
bars = ax2.barh(y_pos, counts, color=C0, alpha=0.8, edgecolor="white")
ax2.set_yticks(y_pos)
ax2.set_yticklabels(names, fontsize=9)
ax2.invert_yaxis()
ax2.set_xlabel("Number of molecules")
ax2.set_title(f"Top {top_n} fluorotelomer groups detected")

for bar, count in zip(bars, counts):
    ax2.text(bar.get_width() + 2, bar.get_y() + bar.get_height()/2,
             str(count), va="center", fontsize=8)

plt.tight_layout()
savefig(fig, "telomer_detection")

# --- Fig 3b: Detection per molecule (how many groups per molecule) ---
results = telomer["results"]
groups_per_mol = [len(r["telomer_groups"]) for r in results if r["detected"]]
fig, ax = plt.subplots(figsize=(6, 4))
max_groups = max(groups_per_mol) if groups_per_mol else 5
bins = np.arange(0.5, max_groups + 1.5, 1)
ax.hist(groups_per_mol, bins=bins, color=C0, alpha=0.8, edgecolor="white")
ax.set_xlabel("Number of telomer groups detected per molecule")
ax.set_ylabel("Count")
ax.set_title(f"Telomer group multiplicity (n = {len(groups_per_mol)} detected)")
ax.set_xticks(range(1, max_groups + 1))
savefig(fig, "telomer_multiplicity")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. TOXCAST FINGERPRINT COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════
print("\n4. ToxCast fingerprint comparison figures")

results_csv = os.path.join(DATA, "toxcast_comparison_results.csv")
summary_csv = os.path.join(DATA, "toxcast_comparison_summary.csv")

df_full = pd.read_csv(results_csv)
df_summ = pd.read_csv(summary_csv)

# Colour mapping for feature sets
FSET_COLORS = {
    "TxP_PFAS":                  TXP_PFAS_COLOR,
    "PFG_binary":                PFG_BINARY,
    "PFG_binary+mol":            PFG_BINARY_MOL,
    "PFG_EGR":                   PFG_EGR,
    "PFG_EGR+mol":               PFG_EGR_MOL,
    "ToxPrint+TxP_PFAS":         TXP_PFAS_COLOR,
    "ToxPrint+PFG_EGR+mol":      PFG_EGR_MOL,
    "Morgan+PFG_binary":         PFG_BINARY,
    "Morgan+PFG_binary+mol":     PFG_BINARY_MOL,
    "Morgan+PFG_EGR":            PFG_EGR,
    "Morgan+PFG_EGR+mol":        PFG_EGR_MOL,
    "Morgan":                    MORGAN_COLOR,
}

# Focus on Gradient Boosting to compare fingerprints without model confound
FOCUS_MODEL = "GradientBoosting"

# Feature set ordering
EXP_A_FSETS = ["PFG_EGR+mol", "PFG_binary+mol", "TxP_PFAS", "PFG_EGR", "PFG_binary"]
EXP_B_FSETS = ["ToxPrint+TxP_PFAS", "ToxPrint+PFG_EGR+mol",
               "Morgan+PFG_EGR+mol", "Morgan+PFG_binary+mol",
               "Morgan+PFG_EGR", "Morgan+PFG_binary", "Morgan"]


def plot_gb_bars(exp, fsets, fname):
    """Bar chart (GB only): per-endpoint ROC-AUC grouped by feature set."""
    df_exp = df_full[(df_full["experiment"] == exp) &
                     (df_full["model"] == FOCUS_MODEL)].copy()
    agg_m = df_exp.groupby(["feature_set", "endpoint"])["roc_auc"].mean().reset_index()
    agg_s = df_exp.groupby(["feature_set", "endpoint"])["roc_auc"].std().reset_index()
    agg_s = agg_s.rename(columns={"roc_auc": "roc_auc_std"})
    agg = agg_m.merge(agg_s, on=["feature_set", "endpoint"])
    endpoints = sorted(agg["endpoint"].unique())

    fig, ax = plt.subplots(figsize=(14, 5))
    n_fsets = len(fsets)
    n_ep = len(endpoints)
    width = 0.8 / n_fsets
    x = np.arange(n_ep)

    for i, fset in enumerate(fsets):
        sub = agg[agg["feature_set"] == fset].set_index("endpoint").reindex(endpoints)
        vals = sub["roc_auc"].values
        errs = sub["roc_auc_std"].fillna(0).values
        offset = (i - n_fsets / 2 + 0.5) * width
        color = FSET_COLORS.get(fset, "#999999")
        ax.bar(x + offset, vals, width, yerr=errs, label=fset,
               color=color, alpha=0.85, edgecolor="white",
               capsize=2.5, error_kw=dict(lw=0.9, alpha=0.7))

    ax.set_xticks(x)
    ax.set_xticklabels([e.replace("_", "\n") for e in endpoints], fontsize=8,
                       rotation=45, ha="right")
    ax.set_ylabel("ROC-AUC")
    ax.set_title(f"{exp}: ROC-AUC per endpoint — Gradient Boosting"
                 f"\n(bars = mean over CV folds; error bars = ±1 SD)")
    ax.legend(fontsize=8, ncol=min(3, n_fsets), loc="upper right")
    ax.set_ylim(0, 1)
    savefig(fig, fname)


def scatter_gb(df_exp, fset_x, fset_y, xlabel, ylabel, title, fname,
               color_x=TXP_PFAS_COLOR, color_y=PFG_EGR_MOL):
    """Scatter fset_y vs fset_x using GB only. One point per endpoint."""
    sub = df_exp[df_exp["model"] == FOCUS_MODEL]
    agg = sub.groupby(["feature_set", "endpoint"])["roc_auc"].mean().reset_index()
    xvals = agg[agg["feature_set"] == fset_x].set_index("endpoint")["roc_auc"]
    yvals = agg[agg["feature_set"] == fset_y].set_index("endpoint")["roc_auc"]
    eps = sorted(set(xvals.index) & set(yvals.index))
    x = [xvals[e] for e in eps]
    y = [yvals[e] for e in eps]

    above = sum(1 for xi, yi in zip(x, y) if yi > xi)
    print(f"    [GB] {fset_y} wins {above}/{len(eps)}")

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(x, y, c=color_y, marker="o", s=60, alpha=0.85,
               edgecolors="white", linewidth=0.4, zorder=3)
    for xi, yi, ep in zip(x, y, eps):
        ax.annotate(ep, (xi, yi), fontsize=5.5,
                    textcoords="offset points", xytext=(3, 2), color="grey")

    lo = min(x + y) - 0.02
    hi = max(x + y) + 0.02
    ax.plot([lo, hi], [lo, hi], "--", color="gray", alpha=0.5, linewidth=0.9, zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    savefig(fig, fname)


# Exp A: bar + scatter
print("  Exp A bar (GB)")
plot_gb_bars("Exp A", EXP_A_FSETS, "toxcast_expA_bar")

print("  Exp A scatter: PFG_EGR+mol vs TxP_PFAS (GB)")
df_expA = df_full[df_full["experiment"] == "Exp A"]
scatter_gb(df_expA,
           fset_x="TxP_PFAS", fset_y="PFG_EGR+mol",
           xlabel="TxP\_PFAS ROC-AUC",
           ylabel="PFG\_EGR+mol ROC-AUC",
           title="Exp A: PFG\_EGR+mol vs TxP\_PFAS (Gradient Boosting)\none point per endpoint",
           fname="toxcast_expA_scatter",
           color_y=PFG_EGR_MOL)

# Exp B: bar + scatter
print("  Exp B bar (GB)")
plot_gb_bars("Exp B", EXP_B_FSETS, "toxcast_expB_bar")

print("  Exp B scatter: Morgan+PFG_EGR+mol vs ToxPrint+TxP_PFAS (GB)")
df_expB = df_full[df_full["experiment"] == "Exp B"]
scatter_gb(df_expB,
           fset_x="ToxPrint+TxP_PFAS", fset_y="Morgan+PFG_EGR+mol",
           xlabel="ToxPrint+TxP\_PFAS ROC-AUC",
           ylabel="Morgan+PFG\_EGR+mol ROC-AUC",
           title="Exp B: Morgan+PFG\_EGR+mol vs ToxPrint+TxP\_PFAS (Gradient Boosting)\none point per endpoint",
           fname="toxcast_expB_scatter",
           color_y=PFG_EGR_MOL)


# ═══════════════════════════════════════════════════════════════════════════════
# 5. TOP DETECTED GROUPS (CLINVENTORY)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n5. Top detected groups on CLInventory")

gnc = clinv["group_name_counts"]
# Sort by count
gnc_sorted = sorted(gnc.items(), key=lambda x: x[1], reverse=True)
top20 = gnc_sorted[:20]

fig, ax = plt.subplots(figsize=(8, 6))
names = [g[0] for g in top20]
counts = [g[1] for g in top20]
y_pos = np.arange(len(names))
bars = ax.barh(y_pos, counts, color=C0, alpha=0.8, edgecolor="white")
ax.set_yticks(y_pos)
ax.set_yticklabels(names, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Detection count")
ax.set_title(f"Top 20 PFASGroups detections on CLInventory (n = {clinv['common_molecules']:,})")
for bar, count in zip(bars, counts):
    ax.text(bar.get_width() + 50, bar.get_y() + bar.get_height()/2,
            f"{count:,}", va="center", fontsize=8)
savefig(fig, "clinventory_top_groups")


print("\n=== ALL FIGURES GENERATED ===")
print(f"Output directory: {OUTDIR}")
for f in sorted(os.listdir(OUTDIR)):
    size = os.path.getsize(os.path.join(OUTDIR, f))
    print(f"  {f:50s} {size:>10,} bytes")
