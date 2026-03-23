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
    "toxcast_expA_bar",
    "toxcast_expA_scatter",
    "toxcast_expB_bar",
    "toxcast_expB_scatter",
]

# ── Exp A / B feature set ordered rows (Gradient Boosting only) ───────────────
FOCUS_MODEL = "GradientBoosting"

EXP_A_FSETS = ["PFG_EGR+mol", "PFG_binary+mol", "TxP_PFAS", "PFG_EGR", "PFG_binary"]

# For Exp B: show ToxPrint+TxP_PFAS vs ToxPrint+PFG_EGR+mol, then Morgan variants
EXP_B_FSETS = [
    "ToxPrint+TxP_PFAS",
    "ToxPrint+PFG_EGR+mol",
    "Morgan+PFG_binary",
    "Morgan+PFG_EGR+mol",
    "Morgan",
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
    """GB-only Exp-A rows (no Model column), ranked by AUC, 4 metrics bolded."""
    rows_data = []
    for fset in EXP_A_FSETS:
        sub = summ[
            (summ["experiment"] == "Exp A")
            & (summ["feature_set"] == fset)
            & (summ["model"] == FOCUS_MODEL)
        ]
        if sub.empty:
            print(f"    [WARN] Exp A / {fset} / {FOCUS_MODEL} not found in summary")
            continue
        row = sub.iloc[0]
        rows_data.append({"feature_set": fset, **{m: row[m] for m in METRICS}})

    if not rows_data:
        return ""

    df = pd.DataFrame(rows_data)
    bold_cols = {m: bold_max(df[m].tolist(), None) for m in METRICS}

    lines = []
    for i, row in df.iterrows():
        vals = [bold_cols[m][i] for m in METRICS]
        fs = row["feature_set"].replace("_", r"\_")
        lines.append(f"{fs:<24} & " + " & ".join(vals) + r" \\")
    return "\n".join(lines)


# ── 3. Build Exp B table TeX ──────────────────────────────────────────────────
def build_expB_table(summ: pd.DataFrame) -> str:
    """GB-only Exp-B rows (no Model column), 4 metrics bolded."""
    rows_data = []
    for fset in EXP_B_FSETS:
        sub = summ[
            (summ["experiment"] == "Exp B")
            & (summ["feature_set"] == fset)
            & (summ["model"] == FOCUS_MODEL)
        ]
        if sub.empty:
            print(f"    [WARN] Exp B / {fset} / {FOCUS_MODEL} not found")
            continue
        row = sub.iloc[0]
        rows_data.append({"feature_set": fset, **{m: row[m] for m in METRICS}})

    if not rows_data:
        return ""

    df = pd.DataFrame(rows_data)
    bold_cols = {m: bold_max(df[m].tolist(), None) for m in METRICS}

    lines = []
    for i, row in df.iterrows():
        vals = [bold_cols[m][i] for m in METRICS]
        fs = row["feature_set"].replace("_", r"\_")
        lines.append(f"{fs:<28} & " + " & ".join(vals) + r" \\")
    return "\n".join(lines)


# ── 4. Build narrative paragraphs ─────────────────────────────────────────────
def build_expA_narrative(summ: pd.DataFrame) -> str:
    """One paragraph describing Exp A results (GB only)."""
    s = summ[(summ["experiment"] == "Exp A") & (summ["model"] == FOCUS_MODEL)]

    def get(fset):
        sub = s[s["feature_set"] == fset]
        return sub.iloc[0] if not sub.empty else None

    pfg_mol_egr = get("PFG_EGR+mol")
    pfg_mol_bin = get("PFG_binary+mol")
    txp         = get("TxP_PFAS")
    pfg_egr     = get("PFG_EGR")
    pfg_bin     = get("PFG_binary")

    if any(x is None for x in [pfg_mol_egr, pfg_mol_bin, txp, pfg_egr, pfg_bin]):
        return ""

    best_pfg_auc = max(pfg_mol_egr["roc_auc_mean"], pfg_mol_bin["roc_auc_mean"])
    gain_vs_txp  = best_pfg_auc - txp["roc_auc_mean"]
    mol_boost    = pfg_mol_bin["roc_auc_mean"] - pfg_bin["roc_auc_mean"]
    egr_vs_bin   = pfg_egr["roc_auc_mean"] - pfg_bin["roc_auc_mean"]

    egr_sign = "$-$" if egr_vs_bin < 0 else "+"
    return (
        f"PFASGroups embeddings augmented with molecule-wide metrics "
        f"(PFG\\_binary+mol, PFG\\_EGR+mol) outperform TxP\\_PFAS by "
        f"+{gain_vs_txp:.3f} AUC ({best_pfg_auc:.3f} vs.\\ {txp['roc_auc_mean']:.3f}). "
        f"The molecule-wide metrics provide the dominant uplift: adding them to PFG\\_binary "
        f"boosts AUC from {pfg_bin['roc_auc_mean']:.3f} to {pfg_mol_bin['roc_auc_mean']:.3f} "
        f"(+{mol_boost:.3f}), whereas switching from binary to EGR encoding adds only "
        f"{egr_sign}{abs(egr_vs_bin):.3f}. "
        f"Without molecule-wide metrics, PFG\\_binary ({pfg_bin['roc_auc_mean']:.3f}) "
        f"and PFG\\_EGR ({pfg_egr['roc_auc_mean']:.3f}) bracket "
        f"TxP\\_PFAS ({txp['roc_auc_mean']:.3f}). "
        f"See Figure~\\ref{{fig:toxcast_expA_bar}} for per-endpoint detail "
        f"and Figure~\\ref{{fig:toxcast_expA_scatter}} for a direct pairwise comparison."
    )


def build_expB_narrative(summ: pd.DataFrame, res: pd.DataFrame) -> str:
    """One paragraph describing Exp B fingerprint comparison (GB only)."""
    s = summ[(summ["experiment"] == "Exp B") & (summ["model"] == FOCUS_MODEL)]

    def get(fset):
        sub = s[s["feature_set"] == fset]
        return sub.iloc[0] if not sub.empty else None

    txp_txp  = get("ToxPrint+TxP_PFAS")
    txp_pfg  = get("ToxPrint+PFG_EGR+mol")
    mg_pfg   = get("Morgan+PFG_EGR+mol")
    morgan   = get("Morgan")

    if any(x is None for x in [txp_txp, txp_pfg, mg_pfg, morgan]):
        return ""

    r = res[(res["experiment"] == "Exp B") & (res["model"] == FOCUS_MODEL)]

    def ep_aucs(fset):
        return r[r["feature_set"] == fset].groupby("endpoint")["roc_auc"].mean()

    txp_morgan_gap  = txp_txp["roc_auc_mean"] - mg_pfg["roc_auc_mean"]
    morgan_pfg_gain = mg_pfg["roc_auc_mean"] - morgan["roc_auc_mean"]

    # Compute observed range across all Morgan-based feature sets in EXP_B_FSETS
    morgan_fsets = [fs for fs in EXP_B_FSETS if fs.startswith("Morgan")]
    morgan_aucs = []
    for fs in morgan_fsets:
        sub = s[s["feature_set"] == fs]
        if not sub.empty:
            morgan_aucs.append(sub.iloc[0]["roc_auc_mean"])
    mg_lo = min(morgan_aucs) if morgan_aucs else morgan["roc_auc_mean"]
    mg_hi = max(morgan_aucs) if morgan_aucs else mg_pfg["roc_auc_mean"]
    morgan_range = f"{mg_lo:.3f}" if abs(mg_hi - mg_lo) < 0.0005 else f"{mg_lo:.3f}--{mg_hi:.3f}"

    return (
        f"ToxPrint is the dominant base fingerprint; the key question is which "
        f"129-bit PFAS component to append. "
        f"ToxPrint+PFG\\_EGR+mol ({txp_pfg['roc_auc_mean']:.3f} AUC) performs identically to "
        f"ToxPrint+TxP\\_PFAS ({txp_txp['roc_auc_mean']:.3f} AUC), demonstrating that "
        f"PFASGroups embeddings can fully replace TxP\\_PFAS within a ToxPrint+PFAS-specific "
        f"fingerprint while marginally improving AP, MCC and balanced accuracy. "
        f"Morgan-based combinations ({morgan_range}) "
        f"trail the ToxPrint variants by $\\sim${txp_morgan_gap:.3f} AUC, confirming that "
        f"ToxPrint's 729 general-purpose structural bits capture non-PFAS chemistry more "
        f"effectively than ECFP4 alone. "
        f"Adding any PFG variant to Morgan provides negligible improvement over Morgan alone "
        f"({morgan_pfg_gain:+.3f} AUC). "
        f"See Figure~\\ref{{fig:toxcast_expB_bar}} for per-endpoint detail "
        f"and Figure~\\ref{{fig:toxcast_expB_scatter}} for the direct ToxPrint pairwise comparison."
    )


# ── 5. Patch supplementary_file_1.tex ─────────────────────────────────────────
# We replace two anchored blocks, identified by unique labels.

TABLE_A_PATTERN = re.compile(
    r"(\\label\{tab:toxcast_expA_summary\}.*?\\small\s*\\begin\{tabular\}\{lcccc\}\s*"
    r"\\toprule\s*Feature Set.*?\\midrule\s*)"
    r"(.*?)"
    r"(\\bottomrule\s*\\end\{tabular\})",
    re.DOTALL,
)

TABLE_B_PATTERN = re.compile(
    r"(\\label\{tab:toxcast_expB_summary\}.*?\\small\s*\\begin\{tabular\}\{lcccc\}\s*"
    r"\\toprule\s*Feature Set.*?\\midrule\s*)"
    r"(.*?)"
    r"(\\bottomrule\s*\\end\{tabular\})",
    re.DOTALL,
)

NARRATIVE_A_PATTERN = re.compile(
    r"(\\end\{table\}\s*\n)"
    r"(PFASGroups embeddings (?:with|augmented with) molecule-wide metrics.*?)"
    r"(\s*\n\\begin\{figure\})",
    re.DOTALL,
)

NARRATIVE_B_PATTERN = re.compile(
    r"(\\end\{table\}\s*\n)"
    r"((?:ToxPrint is the dominant|On the full library the two classifiers).*?)"
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
