"""
update_toxcast_si.py
====================
Post-processing step after compare_fingerprints_toxcast.py completes.

Actions
-------
1. Read new toxcast_comparison_summary.csv + toxcast_comparison_results.csv
2. Regenerate all SI publication figures (generate_si_figures.py)
3. Copy new toxcast PDFs/PNGs from imgs/si/ → classification_article/imgs/
4. Rewrite the two result tables and narrative paragraphs in supplementary_file_1.tex

Usage
-----
    python benchmark/scripts/update_toxcast_si.py
    # or from the benchmark dir:
    python scripts/update_toxcast_si.py
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
BENCH_DIR    = SCRIPT_DIR.parent
DATA_DIR     = BENCH_DIR / "data"
SI_IMGS_DIR  = BENCH_DIR / "imgs" / "si"
ARTICLE_IMGS = Path("/home/luc/git/classification_article/imgs")
TEX_FILE     = Path("/home/luc/git/classification_article/supplementary_files/supplementary_file_1.tex")

SUMMARY_CSV  = DATA_DIR / "toxcast_comparison_summary.csv"
RESULTS_CSV  = DATA_DIR / "toxcast_comparison_results.csv"
SI_SCRIPT    = SCRIPT_DIR / "generate_si_figures.py"

# Toxcast figure stems to copy (benchmark/imgs/si/<stem>.{pdf,png} → article/imgs/)
TOXCAST_FIGS = [
    "toxcast_expA_auc",
    "toxcast_expA_ap",
    "toxcast_expB_auc",
    "toxcast_expB_ap",
    "toxcast_heatmap_roc_auc",
    "toxcast_heatmap_ap",
    "toxcast_scatter_Morgan_plus_PFG_EGR_plus_mol_vs_toxprint",
    "toxcast_scatter_Morgan_plus_PFG_binary_plus_mol_vs_toxprint",
    "toxcast_scatter_PFG_binary_vs_txp",
    "toxcast_scatter_PFG_EGR_vs_txp",
    "toxcast_scatter_PFG_EGR_plus_mol_vs_txp",
    "toxcast_summary",
]

# ── Exp A / B feature set ordered rows ───────────────────────────────────────
EXP_A_FSETS = ["PFG_EGR+mol", "PFG_binary+mol", "TxP_PFAS", "PFG_EGR", "PFG_binary"]

# For Exp B we show: best PFG+mol variants (RF), ToxPrint both models, Morgan baseline
EXP_B_SELECTED = [
    ("Morgan+PFG_EGR+mol",    "RandomForest"),
    ("Morgan+PFG_binary+mol", "RandomForest"),
    ("ToxPrint+TxP_PFAS",     "RandomForest"),
    ("ToxPrint+TxP_PFAS",     "GradientBoosting"),
    ("Morgan",                "RandomForest"),
]

MODEL_ABBREV = {"RandomForest": "RF", "GradientBoosting": "GB"}
METRICS = ["roc_auc_mean", "avg_prec_mean", "mcc_mean", "bal_acc_mean"]
METRIC_LABELS = ["AUC", "AP", "MCC", "Bal Acc"]


def fmt(v):
    """Format a float to 3 decimal places."""
    if pd.isna(v):
        return "---"
    return f"{v:.3f}"


def bold_max(values, idx):
    """Return values with the maximum entry wrapped in \\textbf{}."""
    max_v = max(v for v in values if not pd.isna(v))
    return [f"\\textbf{{{fmt(v)}}}" if (not pd.isna(v) and abs(v - max_v) < 1e-9)
            else fmt(v) for v in values]


# ── 1. Load CSVs ─────────────────────────────────────────────────────────────
def load_data():
    print(f"[1] Reading {SUMMARY_CSV.name} and {RESULTS_CSV.name}")
    summ = pd.read_csv(SUMMARY_CSV)
    res  = pd.read_csv(RESULTS_CSV)
    print(f"    summary: {len(summ)} rows; results: {len(res)} rows")
    return summ, res


# ── 2. Build Exp A table TeX  ─────────────────────────────────────────────────
def build_expA_table(summ: pd.DataFrame) -> str:
    """Best model per Exp-A feature set, ranked by AUC, 4 metrics bolded."""
    rows_data = []
    for fset in EXP_A_FSETS:
        sub = summ[(summ["experiment"] == "Exp A") & (summ["feature_set"] == fset)]
        if sub.empty:
            print(f"    [WARN] Exp A / {fset} not found in summary")
            continue
        best = sub.loc[sub["roc_auc_mean"].idxmax()]
        rows_data.append({
            "feature_set": fset,
            "model": MODEL_ABBREV.get(best["model"], best["model"]),
            **{m: best[m] for m in METRICS},
        })

    if not rows_data:
        return ""

    df = pd.DataFrame(rows_data)
    # Bold per-column max
    bold_cols = {}
    for m in METRICS:
        bold_cols[m] = bold_max(df[m].tolist(), None)

    lines = []
    for i, row in df.iterrows():
        vals = [bold_cols[m][i] for m in METRICS]
        fs = row["feature_set"].replace("+", r"+").replace("_", r"\_")
        lines.append(
            f"{fs:<24} & {row['model']:<2} & "
            + " & ".join(vals)
            + r" \\"
        )
    return "\n".join(lines)


# ── 3. Build Exp B table TeX ──────────────────────────────────────────────────
def build_expB_table(summ: pd.DataFrame) -> str:
    """Selected Exp-B rows with 4 metrics bolded."""
    rows_data = []
    for fset, model_full in EXP_B_SELECTED:
        sub = summ[
            (summ["experiment"] == "Exp B")
            & (summ["feature_set"] == fset)
            & (summ["model"] == model_full)
        ]
        if sub.empty:
            print(f"    [WARN] Exp B / {fset} / {model_full} not found")
            continue
        row = sub.iloc[0]
        rows_data.append({
            "feature_set": fset,
            "model": MODEL_ABBREV.get(model_full, model_full),
            **{m: row[m] for m in METRICS},
        })

    if not rows_data:
        return ""

    df = pd.DataFrame(rows_data)
    bold_cols = {}
    for m in METRICS:
        bold_cols[m] = bold_max(df[m].tolist(), None)

    lines = []
    for i, row in df.iterrows():
        vals = [bold_cols[m][i] for m in METRICS]
        fs = row["feature_set"].replace("+", r"+").replace("_", r"\_")
        lines.append(
            f"{fs:<28} & {row['model']:<2} & "
            + " & ".join(vals)
            + r" \\"
        )
    return "\n".join(lines)


# ── 4. Build narrative paragraphs ─────────────────────────────────────────────
def build_expA_narrative(summ: pd.DataFrame) -> str:
    """One paragraph describing Exp A results."""
    s = summ[summ["experiment"] == "Exp A"]

    def best(fset):
        sub = s[s["feature_set"] == fset]
        return sub.loc[sub["roc_auc_mean"].idxmax()] if not sub.empty else None

    pfg_mol_egr   = best("PFG_EGR+mol")
    pfg_mol_bin   = best("PFG_binary+mol")
    txp           = best("TxP_PFAS")
    pfg_egr       = best("PFG_EGR")
    pfg_bin       = best("PFG_binary")

    if any(x is None for x in [pfg_mol_egr, pfg_mol_bin, txp, pfg_egr, pfg_bin]):
        return ""  # fallback: don't overwrite if data missing

    best_pfg_auc  = max(pfg_mol_egr["roc_auc_mean"], pfg_mol_bin["roc_auc_mean"])
    txp_auc_range = f"{txp['roc_auc_mean']:.3f}"
    # Check if both GB and RF for TxP_PFAS exist
    txp_both = s[s["feature_set"] == "TxP_PFAS"]
    if len(txp_both) > 1:
        lo = txp_both["roc_auc_mean"].min()
        hi = txp_both["roc_auc_mean"].max()
        txp_auc_range = f"{lo:.3f}--{hi:.3f}"

    gain_vs_txp   = best_pfg_auc - txp_both["roc_auc_mean"].max()
    mol_boost     = pfg_mol_bin["roc_auc_mean"] - pfg_bin["roc_auc_mean"]
    egr_boost     = pfg_mol_egr["roc_auc_mean"] - pfg_egr["roc_auc_mean"]
    egr_vs_bin    = pfg_egr["roc_auc_mean"] - pfg_bin["roc_auc_mean"]

    return (
        f"PFASGroups embeddings with molecule-wide metrics "
        f"(PFG\\_binary+mol, PFG\\_EGR+mol) outperform TxP\\_PFAS by "
        f"+{gain_vs_txp:.3f} AUC. "
        f"The molecule-wide metrics provide the dominant uplift: adding them to PFG\\_binary "
        f"boosts AUC from {pfg_bin['roc_auc_mean']:.3f} to {pfg_mol_bin['roc_auc_mean']:.3f} "
        f"(+{mol_boost:.3f}), whereas adding EGR to PFG\\_binary yields only "
        f"+{egr_vs_bin:.3f}. "
        f"Without molecule-wide metrics, PFG\\_binary ({pfg_bin['roc_auc_mean']:.3f}) "
        f"and PFG\\_EGR ({pfg_egr['roc_auc_mean']:.3f}) perform comparably to "
        f"TxP\\_PFAS ({txp_auc_range}). "
        f"See Figures~\\ref{{fig:toxcast_expA_auc}} and~\\ref{{fig:toxcast_expA_ap}}."
    )


def build_expB_narrative(summ: pd.DataFrame, res: pd.DataFrame) -> str:
    """One paragraph describing Exp B results (RF/GB split design)."""
    s = summ[summ["experiment"] == "Exp B"]

    def get(fset, model):
        sub = s[(s["feature_set"] == fset) & (s["model"] == model)]
        return sub.iloc[0] if not sub.empty else None

    rf_pfg = get("Morgan+PFG_EGR+mol", "RandomForest")
    rf_txp = get("ToxPrint+TxP_PFAS",  "RandomForest")
    gb_pfg = get("Morgan+PFG_EGR+mol", "GradientBoosting")
    gb_txp = get("ToxPrint+TxP_PFAS",  "GradientBoosting")

    if any(x is None for x in [rf_pfg, rf_txp, gb_pfg, gb_txp]):
        return ""

    r = res[res["experiment"] == "Exp B"]

    def ep_aucs(fset, model):
        return (r[(r["feature_set"] == fset) & (r["model"] == model)]
                  .groupby("endpoint")["roc_auc"].mean())

    ep_rf_pfg = ep_aucs("Morgan+PFG_EGR+mol", "RandomForest")
    ep_rf_txp = ep_aucs("ToxPrint+TxP_PFAS",  "RandomForest")
    ep_gb_pfg = ep_aucs("Morgan+PFG_EGR+mol", "GradientBoosting")
    ep_gb_txp = ep_aucs("ToxPrint+TxP_PFAS",  "GradientBoosting")

    rf_ep = sorted(ep_rf_pfg.index.intersection(ep_rf_txp.index))
    gb_ep = sorted(ep_gb_pfg.index.intersection(ep_gb_txp.index))
    n_ep = len(rf_ep)

    n_rf_wins = (ep_rf_pfg[rf_ep] > ep_rf_txp[rf_ep]).sum()
    n_gb_txp_wins = (ep_gb_txp[gb_ep] > ep_gb_pfg[gb_ep]).sum()
    mean_rf_delta = (ep_rf_pfg[rf_ep] - ep_rf_txp[rf_ep]).mean()

    return (
        f"On the full library the two classifiers tell opposite stories. "
        f"With Random Forest, Morgan+PFG\\_EGR+mol ({rf_pfg['roc_auc_mean']:.3f} AUC) "
        f"and ToxPrint+TxP\\_PFAS ({rf_txp['roc_auc_mean']:.3f} AUC) are essentially identical, "
        f"with Morgan+PFG improving {n_rf_wins}/{n_ep} endpoints "
        f"(mean $\\Delta$AUC $\\approx${mean_rf_delta:+.3f}). "
        f"With Gradient Boosting, ToxPrint+TxP\\_PFAS has a clear advantage "
        f"({gb_txp['roc_auc_mean']:.3f} vs.\\ {gb_pfg['roc_auc_mean']:.3f} AUC, "
        f"improving {n_gb_txp_wins}/{len(gb_ep)} endpoints), "
        f"because ToxPrint's 729 general-purpose structural bits provide a richer encoding "
        f"of non-PFAS chemistry than ECFP4 alone; PFASGroups bits add little signal when "
        f"GB is already exploiting ToxPrint's breadth. "
        f"The per-endpoint scatter in Figure~\\ref{{fig:toxcast_scatter_B}} makes this split "
        f"explicit: RF points (circles) cluster near the diagonal, while GB points (squares) "
        f"systematically fall below it. "
        f"See also Figures~\\ref{{fig:toxcast_expB_auc}} and~\\ref{{fig:toxcast_expB_ap}}."
    )


# ── 5. Patch supplementary_file_1.tex ─────────────────────────────────────────
# We replace two anchored blocks, identified by unique labels.

TABLE_A_PATTERN = re.compile(
    r"(\\label\{tab:toxcast_expA_summary\}.*?\\small\s*\\begin\{tabular\}\{llcccc\}\s*"
    r"\\toprule\s*Feature Set.*?\\midrule\s*)"
    r"(.*?)"
    r"(\\bottomrule\s*\\end\{tabular\})",
    re.DOTALL,
)

TABLE_B_PATTERN = re.compile(
    r"(\\label\{tab:toxcast_expB_summary\}.*?\\small\s*\\begin\{tabular\}\{llcccc\}\s*"
    r"\\toprule\s*Feature Set.*?\\midrule\s*)"
    r"(.*?)"
    r"(\\bottomrule\s*\\end\{tabular\})",
    re.DOTALL,
)

# Narrative: replace the paragraph just after expA table \end{table}
# (ends at \begin{figure} for expA_auc)
NARRATIVE_A_PATTERN = re.compile(
    r"(\\end\{table\}\s*\n)"
    r"(PFASGroups embeddings with molecule-wide metrics.*?See Figures~\\ref\{fig:toxcast_expA_auc\} and~\\ref\{fig:toxcast_expA_ap\}\.)"
    r"(\s*\n\\begin\{figure\})",
    re.DOTALL,
)

NARRATIVE_B_PATTERN = re.compile(
    r"(\\end\{table\}\s*\n)"
    r"(On the full library the two classifiers tell opposite stories\..*?"
    r"See also Figures~\\ref\{fig:toxcast_expB_auc\} and~\\ref\{fig:toxcast_expB_ap\}\.)"
    r"(\s*\n\\begin\{figure\})",
    re.DOTALL,
)


def patch_tex(tex_path: Path, table_a: str, table_b: str, narr_a: str, narr_b: str):
    print(f"[4] Patching {tex_path.name}")
    content = tex_path.read_text(encoding="utf-8")
    original = content

    if table_a:
        content, n = TABLE_A_PATTERN.subn(
            lambda m: m.group(1) + table_a + "\n" + m.group(3), content
        )
        print(f"    Exp A table: {n} replacement(s)")

    if table_b:
        content, n = TABLE_B_PATTERN.subn(
            lambda m: m.group(1) + table_b + "\n" + m.group(3), content
        )
        print(f"    Exp B table: {n} replacement(s)")

    if narr_a:
        content, n = NARRATIVE_A_PATTERN.subn(
            lambda m: m.group(1) + narr_a + m.group(3), content
        )
        print(f"    Exp A narrative: {n} replacement(s)")

    if narr_b:
        content, n = NARRATIVE_B_PATTERN.subn(
            lambda m: m.group(1) + narr_b + m.group(3), content
        )
        print(f"    Exp B narrative: {n} replacement(s)")

    # Always save (idempotent writes are fine); warn only if nothing matched at all
    total_matches = sum(
        1 for pattern, text in [
            (TABLE_A_PATTERN, table_a), (TABLE_B_PATTERN, table_b),
            (NARRATIVE_A_PATTERN, narr_a), (NARRATIVE_B_PATTERN, narr_b),
        ]
        if text and pattern.search(content)
    )
    if total_matches == 0:
        print("    [WARN] No patterns matched — check anchors")
    else:
        tex_path.write_text(content, encoding="utf-8")
        print(f"    Saved {tex_path}")


# ── 6. Regenerate SI figures ──────────────────────────────────────────────────
def regen_figures():
    print(f"\n[2] Running {SI_SCRIPT.name}")
    # Use the same Python that's running this script
    result = subprocess.run(
        [sys.executable, str(SI_SCRIPT)],
        cwd=str(BENCH_DIR),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print("[ERROR] generate_si_figures.py failed:")
        print(result.stderr[-3000:])
        raise RuntimeError("SI figure generation failed")
    print(result.stdout[-1000:])


# ── 7. Copy figures to article ────────────────────────────────────────────────
def copy_figures():
    print(f"\n[3] Copying toxcast figures → {ARTICLE_IMGS}")
    ARTICLE_IMGS.mkdir(parents=True, exist_ok=True)
    for stem in TOXCAST_FIGS:
        for ext in (".pdf", ".png"):
            src = SI_IMGS_DIR / (stem + ext)
            if src.exists():
                dst = ARTICLE_IMGS / src.name
                shutil.copy2(src, dst)
                print(f"    {src.name} → {dst}")
            else:
                print(f"    [WARN] not found: {src}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("update_toxcast_si.py")
    print("=" * 60)

    if not SUMMARY_CSV.exists():
        sys.exit(f"[ERROR] {SUMMARY_CSV} not found — run compare_fingerprints_toxcast.py first")
    if not RESULTS_CSV.exists():
        sys.exit(f"[ERROR] {RESULTS_CSV} not found")

    summ, res = load_data()

    # Build new TeX content
    table_a  = build_expA_table(summ)
    table_b  = build_expB_table(summ)
    narr_a   = build_expA_narrative(summ)
    narr_b   = build_expB_narrative(summ, res)

    print("\n--- Exp A table rows ---")
    print(table_a)
    print("\n--- Exp B table rows ---")
    print(table_b)
    print("\n--- Exp A narrative ---")
    print(narr_a)
    print("\n--- Exp B narrative ---")
    print(narr_b)

    # Regenerate SI figures
    regen_figures()

    # Copy to article
    copy_figures()

    # Patch TeX
    patch_tex(TEX_FILE, table_a, table_b, narr_a, narr_b)

    print("\n[DONE] All steps complete.")


if __name__ == "__main__":
    main()
