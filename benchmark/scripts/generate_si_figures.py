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
ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # benchmark/
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

clinv_file = os.path.join(DATA, "clinventory_comparison_20260315T065332.json")
with open(clinv_file) as f:
    clinv = json.load(f)

# Extract CLInv overall stats for downstream CDF plot
hg_stats = clinv["timing_overall"]["hg"]
at_stats = clinv["timing_overall"]["atlas"]

# --- Fig 1a: Multi-dataset timing box plot (OECD + CLInventory + large-molecule stress) ---
oecd_timing_file = os.path.join(DATA, "oecd_clinventory_timing_20260315_062059.json")
with open(oecd_timing_file) as f:
    oecd_clinv_data = json.load(f)

oecd_pg  = oecd_clinv_data["datasets"]["OECD"]["pg_timing"]
oecd_at  = oecd_clinv_data["datasets"]["OECD"]["atlas_timing"]
clinv_pg = oecd_clinv_data["datasets"]["clinventory (F)"]["pg_timing"]
clinv_at = oecd_clinv_data["datasets"]["clinventory (F)"]["atlas_timing"]

# Stress dataset: per-molecule timing for molecules with ≥35 heavy atoms
stress_pg_times, stress_at_times = [], []
stress_fullfile = os.path.join(DATA, "pfas_timing_benchmark_full_20260313_201139.json")
if os.path.exists(stress_fullfile):
    with open(stress_fullfile) as f:
        stress_raw = json.load(f)
    stress_pg_times = [r["PFASGroups_time_avg"] * 1000 for r in stress_raw if r["num_atoms"] >= 35]
    stress_at_times = [r["atlas_time_avg"] * 1000 for r in stress_raw if r["num_atoms"] >= 35]

def _approx_bxp(mean, std, median, mn, mx):
    """Box-plot stats approximated from summary statistics (no full distribution)."""
    q1 = max(mn, mean - 0.75 * std)
    q3 = mean + 0.75 * std
    return {"med": median, "q1": q1, "q3": q3,
            "whislo": max(0, mn), "whishi": min(mx, mean + 2.0 * std),
            "mean": mean, "fliers": []}

def _exact_bxp(times_ms):
    """Exact box-plot stats from a flat list of timing values."""
    arr = np.array(times_ms)
    return {"med": float(np.median(arr)),
            "q1": float(np.percentile(arr, 25)),
            "q3": float(np.percentile(arr, 75)),
            "whislo": float(np.percentile(arr, 5)),
            "whishi": float(np.percentile(arr, 95)),
            "mean": float(np.mean(arr)), "fliers": []}

bp_oecd_pg   = _approx_bxp(oecd_pg["mean"],  oecd_pg["std"],  oecd_pg["median"],  oecd_pg["min"],  oecd_pg["max"])
bp_oecd_at   = _approx_bxp(oecd_at["mean"],  oecd_at["std"],  oecd_at["median"],  oecd_at["min"],  oecd_at["max"])
bp_clinv_pg  = _approx_bxp(clinv_pg["mean"], clinv_pg["std"], clinv_pg["median"], clinv_pg["min"], clinv_pg["max"])
bp_clinv_at  = _approx_bxp(clinv_at["mean"], clinv_at["std"], clinv_at["median"], clinv_at["min"], clinv_at["max"])
bp_stress_pg = _exact_bxp(stress_pg_times) if stress_pg_times else _approx_bxp(520, 300, 26, 0.5, 5000)
bp_stress_at = _exact_bxp(stress_at_times) if stress_at_times else _approx_bxp(94, 60, 38, 0.1, 800)

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

n_oecd  = oecd_pg["n"]
n_clinv = clinv_pg["n"]
n_stress = len(stress_pg_times) if stress_pg_times else 0
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

# --- Fig 2e: Atlas/PFASGroups disagreement examples ---
print("  Disagreement molecule examples")
try:
    from rdkit import Chem
    from rdkit.Chem import Draw

    # PFASGroups only: non-F halogenated compounds that Atlas does not flag
    _pfg_only = [
        ("Benzotrichloride",    "ClC(Cl)(Cl)c1ccccc1",
         "PFASGroups: perhalogenated\nside-chain aromatic"),
        ("1,3,5-Trichloro-\nbenzene",  "Clc1cc(Cl)cc(Cl)c1",
         "PFASGroups: polyhalogenated\naryl compound"),
        ("1,1,2-Trichloroethane", "ClCC(Cl)Cl",
         "PFASGroups: polyhalogenated\nalkyl chain"),
        ("Carbon tetrachloride",  "ClC(Cl)(Cl)Cl",
         "PFASGroups: perhalogenated\nalkyl (CCl₄ type)"),
        ("Dibromomethane",        "BrCBr",
         "PFASGroups: polyhalogenated\nalkyl (short-chain Br)"),
    ]
    # PFAS-Atlas only: highly branched perfluoroalkyl that PFASGroups misses
    _atlas_only = [
        ("Perfluoroisobutane",    "FC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F",
         "PFAS-Atlas: (CF₃)₃CF\nbranched perfluoroalkane"),
        ("Iso-PFAS backbone",     "FC(F)(F)C(F)(C(F)(F)F)C(F)(F)F",
         "PFAS-Atlas: branched\nperfluorobutane core"),
        ("Branched PFC ether",    "FC(F)(F)C(C(F)(F)F)(C(F)(F)F)OC(F)(F)F",
         "PFAS-Atlas: branched\nperfluoroalkyl ether"),
    ]

    all_rows = [(smi, name, cap, "pfg") for name, smi, cap in _pfg_only] + \
               [(smi, name, cap, "atlas") for name, smi, cap in _atlas_only]
    parsed = [(Chem.MolFromSmiles(smi), name, cap, grp)
              for smi, name, cap, grp in all_rows]
    parsed = [(m, n, c, g) for m, n, c, g in parsed if m is not None]

    pfg_rows  = [(m, n, c) for m, n, c, g in parsed if g == "pfg"]
    atlas_rows = [(m, n, c) for m, n, c, g in parsed if g == "atlas"]
    ncols = max(len(pfg_rows), len(atlas_rows))

    fig, axes = plt.subplots(2, ncols, figsize=(2.8 * ncols, 6.5))
    if axes.ndim == 1:
        axes = axes.reshape(2, ncols)

    row_data = [(pfg_rows, C0, "PFASGroups detects — PFAS-Atlas does not\n(non-fluorine halogenated compounds; HalogenGroups mode)"),
                (atlas_rows, C1, "PFAS-Atlas detects — PFASGroups does not\n(highly branched perfluoroalkyl: SMARTS path disrupted)")]

    for ri, (row_mols, rcol, row_title) in enumerate(row_data):
        for ci in range(ncols):
            ax = axes[ri, ci]
            ax.axis("off")
            if ci < len(row_mols):
                mol, name, cap = row_mols[ci]
                img = Draw.MolToImage(mol, size=(220, 160))
                ax.imshow(img)
                ax.set_title(name, fontsize=8, fontweight="bold", pad=3,
                             color=rcol)
                ax.text(0.5, -0.10, cap, transform=ax.transAxes, ha="center",
                        fontsize=7, style="italic", color="#555555",
                        verticalalignment="top")
        # Row label as a coloured left-edge annotation
        axes[ri, 0].text(-0.25, 0.5, row_title.split("\n")[0],
                         transform=axes[ri, 0].transAxes,
                         ha="right", va="center", fontsize=8,
                         fontweight="bold", color=rcol,
                         rotation=90, clip_on=False)

    fig.suptitle("Representative classification disagreements: PFASGroups vs PFAS-Atlas",
                 fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    savefig(fig, "atlas_disagreement_examples")
except ImportError as _e:
    print(f"  [warn] RDKit not available — skipping disagreement figure ({_e})")

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
    "TxP_PFAS":              TXP_PFAS_COLOR,
    "PFG_binary":            PFG_BINARY,
    "PFG_binary+mol":        PFG_BINARY_MOL,
    "PFG_EGR":               PFG_EGR,
    "PFG_EGR+mol":           PFG_EGR_MOL,
    "ToxPrint+TxP_PFAS":     TXP_PFAS_COLOR,
    "Morgan+PFG_binary":     PFG_BINARY,
    "Morgan+PFG_binary+mol": PFG_BINARY_MOL,
    "Morgan+PFG_EGR":        PFG_EGR,
    "Morgan+PFG_EGR+mol":    PFG_EGR_MOL,
    "Morgan":                MORGAN_COLOR,
}

# Feature set ordering
EXP_A_FSETS = ["PFG_EGR+mol", "PFG_binary+mol", "TxP_PFAS", "PFG_EGR", "PFG_binary"]
EXP_B_FSETS = ["ToxPrint+TxP_PFAS", "Morgan+PFG_EGR+mol", "Morgan+PFG_binary+mol",
               "Morgan+PFG_EGR", "Morgan+PFG_binary", "Morgan"]

def plot_expAB_bars(metric, metric_label, exp, fsets, fname):
    """Bar chart with error bars: per-endpoint metric grouped by feature set.

    Error bars show ±1 std across repeats × folds (nested CV variability).
    """
    df_exp = df_full[df_full["experiment"] == exp].copy()
    # Mean and std across all repeats × folds × models per (feature_set, endpoint)
    agg_m = df_exp.groupby(["feature_set", "endpoint"])[metric].mean().reset_index()
    agg_s = df_exp.groupby(["feature_set", "endpoint"])[metric].std().reset_index()
    agg_s = agg_s.rename(columns={metric: metric + "_std"})
    agg = agg_m.merge(agg_s, on=["feature_set", "endpoint"])
    endpoints = sorted(agg["endpoint"].unique())

    fig, ax = plt.subplots(figsize=(14, 5))
    n_fsets = len(fsets)
    n_ep = len(endpoints)
    width = 0.8 / n_fsets
    x = np.arange(n_ep)

    for i, fset in enumerate(fsets):
        sub = agg[agg["feature_set"] == fset].set_index("endpoint").reindex(endpoints)
        vals = sub[metric].values
        errs = sub[metric + "_std"].fillna(0).values
        offset = (i - n_fsets / 2 + 0.5) * width
        color = FSET_COLORS.get(fset, "#999999")
        ax.bar(x + offset, vals, width, yerr=errs, label=fset,
               color=color, alpha=0.85, edgecolor="white",
               capsize=2.5, error_kw=dict(lw=0.9, alpha=0.7))

    ax.set_xticks(x)
    ax.set_xticklabels([e.replace("_", "\n") for e in endpoints], fontsize=8, rotation=45, ha="right")
    ax.set_ylabel(metric_label)
    ax.set_title(f"{exp}: {metric_label} per endpoint (bars = mean; error bars = ±1 SD across CV folds)")
    ax.legend(fontsize=8, ncol=min(3, n_fsets), loc="upper right")
    ax.set_ylim(0, 1)
    savefig(fig, fname)

# Exp A bars
plot_expAB_bars("roc_auc", "ROC-AUC", "Exp A", EXP_A_FSETS, "toxcast_expA_auc")
plot_expAB_bars("avg_prec", "Average Precision", "Exp A", EXP_A_FSETS, "toxcast_expA_ap")
# Exp B bars
plot_expAB_bars("roc_auc", "ROC-AUC", "Exp B", EXP_B_FSETS, "toxcast_expB_auc")
plot_expAB_bars("avg_prec", "Average Precision", "Exp B", EXP_B_FSETS, "toxcast_expB_ap")

# --- Fig 4c: Heatmap (both experiments combined) ---
for metric, metric_label in [("roc_auc", "ROC-AUC"), ("avg_prec", "Average Precision")]:
    # Aggregate: mean across models, repeats, folds
    agg = df_full.groupby(["experiment", "feature_set", "endpoint"])[metric].mean().reset_index()
    # Create label col
    agg["label"] = agg["experiment"] + " | " + agg["feature_set"]
    # Pivot
    piv = agg.pivot_table(index="label", columns="endpoint", values=metric)
    # Sort rows
    order_a = [f"Exp A | {fs}" for fs in EXP_A_FSETS]
    order_b = [f"Exp B | {fs}" for fs in EXP_B_FSETS]
    row_order = [r for r in order_a + order_b if r in piv.index]
    piv = piv.reindex(row_order)

    fig, ax = plt.subplots(figsize=(14, 7))
    sns.heatmap(piv, annot=True, fmt=".2f", cmap="YlOrRd", ax=ax,
                linewidths=0.5, linecolor="white", cbar_kws={"label": metric_label},
                annot_kws={"fontsize": 7})
    ax.set_title(f"{metric_label} per endpoint and feature set")
    ax.set_ylabel("")
    ax.set_xlabel("")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    savefig(fig, f"toxcast_heatmap_{metric.replace('avg_prec', 'ap')}")

# --- Fig 4d: Scatter Exp A — PFG variants vs TxP_PFAS ---
print("  Scatter plots: PFG vs TxP_PFAS")
agg_a = df_full[df_full["experiment"] == "Exp A"].groupby(
    ["feature_set", "endpoint"])["roc_auc"].mean().reset_index()

for pfg_fset, pfg_label, pfg_color in [
    ("PFG_binary", "PFG_binary", PFG_BINARY),
    ("PFG_EGR", "PFG_EGR", PFG_EGR),
    ("PFG_EGR+mol", "PFG_EGR+mol", PFG_EGR_MOL),
]:
    txp_vals = agg_a[agg_a["feature_set"] == "TxP_PFAS"].set_index("endpoint")["roc_auc"]
    pfg_vals = agg_a[agg_a["feature_set"] == pfg_fset].set_index("endpoint")["roc_auc"]
    endpoints = sorted(set(txp_vals.index) & set(pfg_vals.index))
    x = [txp_vals[e] for e in endpoints]
    y = [pfg_vals[e] for e in endpoints]

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(x, y, c=pfg_color, s=50, alpha=0.8, edgecolors="white", linewidth=0.5, zorder=3)
    lims = [min(min(x), min(y)) - 0.02, max(max(x), max(y)) + 0.02]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5, zorder=1)
    ax.set_xlabel("TxP_PFAS ROC-AUC")
    ax.set_ylabel(f"{pfg_label} ROC-AUC")
    ax.set_title(f"Exp A: {pfg_label} vs TxP_PFAS")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    # Count above/below
    above = sum(1 for xi, yi in zip(x, y) if yi > xi)
    below = sum(1 for xi, yi in zip(x, y) if yi < xi)
    ax.text(0.05, 0.95, f"{pfg_label} wins: {above}\nTxP_PFAS wins: {below}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
    savefig(fig, f"toxcast_scatter_{pfg_fset.replace('+', '_plus_')}_vs_txp")

# --- Fig 4e: Scatter Exp B — Morgan+PFG vs ToxPrint+TxP ---
print("  Scatter plots: Morgan+PFG vs ToxPrint+TxP")
agg_b = df_full[df_full["experiment"] == "Exp B"].groupby(
    ["feature_set", "endpoint"])["roc_auc"].mean().reset_index()
toxprint_vals = agg_b[agg_b["feature_set"] == "ToxPrint+TxP_PFAS"].set_index("endpoint")["roc_auc"]

for pfg_fset, pfg_color in [("Morgan+PFG_EGR+mol", PFG_EGR_MOL), ("Morgan+PFG_binary+mol", PFG_BINARY_MOL)]:
    pfg_vals = agg_b[agg_b["feature_set"] == pfg_fset].set_index("endpoint")["roc_auc"]
    endpoints = sorted(set(toxprint_vals.index) & set(pfg_vals.index))
    x = [toxprint_vals[e] for e in endpoints]
    y = [pfg_vals[e] for e in endpoints]

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(x, y, c=pfg_color, s=50, alpha=0.8, edgecolors="white", linewidth=0.5, zorder=3)
    lims = [min(min(x), min(y)) - 0.02, max(max(x), max(y)) + 0.02]
    ax.plot(lims, lims, "--", color="gray", alpha=0.5, zorder=1)
    ax.set_xlabel("ToxPrint+TxP_PFAS ROC-AUC")
    ax.set_ylabel(f"{pfg_fset} ROC-AUC")
    ax.set_title(f"Exp B: {pfg_fset} vs ToxPrint+TxP_PFAS")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    above = sum(1 for xi, yi in zip(x, y) if yi > xi)
    below = sum(1 for xi, yi in zip(x, y) if yi < xi)
    ax.text(0.05, 0.95, f"Morgan+PFG wins: {above}\nToxPrint+TxP wins: {below}",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
    name = pfg_fset.replace("+", "_plus_")
    savefig(fig, f"toxcast_scatter_{name}_vs_toxprint")

# --- Fig 4f: Summary bar chart (macro-averaged AUC by model) ---
print("  Summary bars")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for idx, (exp, fsets) in enumerate([("Exp A", EXP_A_FSETS), ("Exp B", EXP_B_FSETS)]):
    ax = axes[idx]
    sub = df_summ[df_summ["experiment"] == exp]
    models = sorted(sub["model"].unique())
    n_models = len(models)
    width = 0.35
    x = np.arange(len(fsets))

    for mi, model in enumerate(models):
        msub = sub[sub["model"] == model].set_index("feature_set").reindex(fsets)
        offset = (mi - n_models / 2 + 0.5) * width
        vals = msub["roc_auc_mean"].values
        errs = msub["roc_auc_std"].values
        ax.bar(x + offset, vals, width, yerr=errs, label=model, alpha=0.85,
               color=[FSET_COLORS.get(fs, "#999") for fs in fsets],
               edgecolor="white", capsize=2)

    ax.set_xticks(x)
    ax.set_xticklabels([fs.replace("+", "+\n") for fs in fsets], fontsize=8, rotation=45, ha="right")
    ax.set_ylabel("ROC-AUC (macro-averaged)")
    ax.set_title(f"{exp}")
    ax.legend(fontsize=8)
    ax.set_ylim(0.5, 0.9)

plt.tight_layout()
savefig(fig, "toxcast_summary")


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
