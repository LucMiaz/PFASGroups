#!/usr/bin/env python3
"""
Generate LaTeX tables and narrative text for the clinventory comparison.

Reads the three JSON files produced by:
  • classify_halogengroups_clinventory.py
  • classify_pfasatlas_clinventory.py
  • compare_clinventory_classifiers.py

and writes publication-ready .tex fragments to:
  reports/clinventory_comparison_tables.tex
  reports/clinventory_comparison_text.tex

Usage (from the benchmark/ directory):
    python scripts/generate_clinventory_latex.py \\
        [--hg-file  data/halogengroups_clinventory_*.json] \\
        [--atlas-file data/pfasatlas_clinventory_*.json]   \\
        [--cmp-file  data/clinventory_comparison_*.json]
"""

import argparse
import glob
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

SCRIPT_DIR   = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BRACKETS = [
    "tiny (<10)",
    "small (10-19)",
    "medium (20-34)",
    "large (35-59)",
    "very large (>=60)",
]

HALOGEN_CLASSES = [
    "F only",
    "Mixed (F + other)",
    "Other halogen (no F)",
    "No halogen",
]


def latest_file(pattern: str) -> Optional[Path]:
    files = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)
    return Path(files[0]) if files else None


def esc(text: str) -> str:
    """Escape LaTeX special characters."""
    for old, new in [
        ("\\", r"\textbackslash{}"),
        ("&", r"\&"), ("%", r"\%"), ("$", r"\$"),
        ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"),
        ("~", r"\textasciitilde{}"), ("^", r"\^{}"),
    ]:
        text = text.replace(old, new)
    return text


def fmt(val, decimals: int = 3) -> str:
    """Format a float for LaTeX, handling None/nan."""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return "--"
    return f"{val:.{decimals}f}"


def fmt_int(val) -> str:
    if val is None:
        return "--"
    return f"{int(val):,}"


def pct(num, den) -> str:
    if not den:
        return "--"
    return f"{100.0 * num / den:.1f}\\%"


# ---------------------------------------------------------------------------
# Table generators
# ---------------------------------------------------------------------------

def table_overall_timing(hg: dict, atlas: dict) -> str:
    """Overall timing statistics table."""
    hg_s    = hg.get("timing_overall", {})
    atlas_s = atlas.get("timing_overall", {})

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Overall per-molecule classification timing (ms) for",
        r"  HalogenGroups and PFAS-Atlas on the C\&L inventory.}",
        r"\label{tab:clinventory_timing_overall}",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"\textbf{Method} & \textbf{$n$} & \textbf{Median} & \textbf{Mean}"
        r" & \textbf{Std} & \textbf{P95} & \textbf{Total (s)} \\",
        r"\midrule",
    ]

    for label, s, data in [
        ("HalogenGroups", hg_s, hg),
        ("PFAS-Atlas",    atlas_s, atlas),
    ]:
        n     = fmt_int(s.get("n"))
        med   = fmt(s.get("median_ms"), 3)
        mean  = fmt(s.get("mean_ms"), 3)
        std   = fmt(s.get("std_ms"), 3)
        p95   = fmt(s.get("p95_ms"), 3)
        total = fmt(data.get("total_classification_time_s"), 1)
        lines.append(f"  {esc(label)} & {n} & {med} & {mean} & {std} & {p95} & {total} \\\\")

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def table_timing_by_bracket(hg: dict, atlas: dict) -> str:
    """Timing by heavy-atom count bracket."""
    hg_br    = hg.get("timing_by_atom_bracket",    {})
    atlas_br = atlas.get("timing_by_atom_bracket", {})

    brackets = [b for b in BRACKETS if b in hg_br or b in atlas_br]

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Median per-molecule timing (ms) by molecular size (heavy-atom count)",
        r"  for both classifiers.}",
        r"\label{tab:clinventory_timing_bracket}",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r" & \multicolumn{2}{c}{\textbf{HalogenGroups}} &"
        r" \multicolumn{2}{c}{\textbf{PFAS-Atlas}} &"
        r" \multicolumn{2}{c}{\textbf{Ratio (HG/Atlas)}} \\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
        r"\textbf{Bracket} & \textbf{$n$} & \textbf{Median (ms)}"
        r" & \textbf{$n$} & \textbf{Median (ms)}"
        r" & \textbf{Median} & \\",
        r"\midrule",
    ]

    for b in brackets:
        hg_s    = hg_br.get(b, {})
        atlas_s = atlas_br.get(b, {})
        n_hg    = fmt_int(hg_s.get("n"))
        med_hg  = fmt(hg_s.get("median_ms"), 3)
        n_at    = fmt_int(atlas_s.get("n"))
        med_at  = fmt(atlas_s.get("median_ms"), 3)
        # ratio
        hm = hg_s.get("median_ms") or 0
        am = atlas_s.get("median_ms") or 0
        ratio = fmt(hm / am if am > 0 else float("nan"), 2)
        lines.append(f"  {esc(b)} & {n_hg} & {med_hg} & {n_at} & {med_at} & {ratio} \\\\")

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def table_timing_by_halogen(hg: dict, atlas: dict) -> str:
    """Timing by halogen composition class."""
    hg_h    = hg.get("timing_by_halogen_class",    {})
    atlas_h = atlas.get("timing_by_halogen_class", {})

    classes = [c for c in HALOGEN_CLASSES if c in hg_h or c in atlas_h]

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Median per-molecule timing (ms) by halogen composition"
        r" for both classifiers.}",
        r"\label{tab:clinventory_timing_halogen}",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r" & \multicolumn{2}{c}{\textbf{HalogenGroups}} &"
        r" \multicolumn{2}{c}{\textbf{PFAS-Atlas}} &"
        r" \multicolumn{2}{c}{\textbf{Ratio (HG/Atlas)}} \\",
        r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
        r"\textbf{Halogen class} & \textbf{$n$} & \textbf{Median (ms)}"
        r" & \textbf{$n$} & \textbf{Median (ms)}"
        r" & \textbf{Median} & \\",
        r"\midrule",
    ]

    for c in classes:
        hg_s    = hg_h.get(c, {})
        atlas_s = atlas_h.get(c, {})
        n_hg    = fmt_int(hg_s.get("n"))
        med_hg  = fmt(hg_s.get("median_ms"), 3)
        n_at    = fmt_int(atlas_s.get("n"))
        med_at  = fmt(atlas_s.get("median_ms"), 3)
        hm = hg_s.get("median_ms") or 0
        am = atlas_s.get("median_ms") or 0
        ratio = fmt(hm / am if am > 0 else float("nan"), 2)
        lines.append(f"  {esc(c)} & {n_hg} & {med_hg} & {n_at} & {med_at} & {ratio} \\\\")

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def table_agreement(cmp: dict) -> str:
    """2×2 agreement contingency table."""
    am    = cmp.get("agreement_matrix", {})
    total = am.get("total", 0)

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Classification agreement between HalogenGroups (presence of",
        r"  F-containing group) and PFAS-Atlas (any PFAS classification) across",
        r"  the common set of molecules.  Cohen's $\kappa$ measures agreement",
        r"  beyond chance.}",
        r"\label{tab:clinventory_agreement}",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"\textbf{Category} & \textbf{Molecules} & \textbf{\%} \\",
        r"\midrule",
        f"  Both classify as F/PFAS & {fmt_int(am.get('both_pfas'))} & {pct(am.get('both_pfas',0), total)} \\\\",
        f"  HalogenGroups only (F detected) & {fmt_int(am.get('only_hg'))} & {pct(am.get('only_hg',0), total)} \\\\",
        f"  PFAS-Atlas only & {fmt_int(am.get('only_atlas'))} & {pct(am.get('only_atlas',0), total)} \\\\",
        f"  Neither classifies & {fmt_int(am.get('neither'))} & {pct(am.get('neither',0), total)} \\\\",
        r"\midrule",
        f"  \\textbf{{Total}} & \\textbf{{{fmt_int(total)}}} & 100.0\\% \\\\",
        r"\midrule",
        f"  Overall agreement & \\multicolumn{{2}}{{r}}{{{fmt(am.get('agreement_pct'), 1)}\\%}} \\\\",
        f"  Cohen's $\\kappa$ & \\multicolumn{{2}}{{r}}{{{fmt(am.get('cohen_kappa'), 3)}}} \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ]
    return "\n".join(lines)


def table_atlas_classes(atlas: dict) -> str:
    """PFAS-Atlas class-1 distribution table."""
    dist  = atlas.get("class1_distribution", {})
    total = sum(dist.values())
    items = sorted(dist.items(), key=lambda x: -x[1])

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Distribution of PFAS-Atlas primary classifications on"
        r"  the fluorine-containing molecules of the C\&L inventory.}",
        r"\label{tab:clinventory_atlas_classes}",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"\textbf{PFAS-Atlas class} & \textbf{Molecules} & \textbf{\%} \\",
        r"\midrule",
    ]
    for cls, cnt in items:
        lines.append(f"  {esc(cls)} & {fmt_int(cnt)} & {pct(cnt, total)} \\\\")
    lines += [
        r"\midrule",
        f"  \\textbf{{Total}} & \\textbf{{{fmt_int(total)}}} & 100.0\\% \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ]
    return "\n".join(lines)


def table_hg_top_groups(cmp: dict, n: int = 20) -> str:
    """Top-N most frequently matched HalogenGroups groups."""
    counts = cmp.get("group_name_counts", {})
    items  = sorted(counts.items(), key=lambda x: -x[1])[:n]
    total  = sum(counts.values())

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption{{Top {n} most frequently matched HalogenGroups functional groups",
        r"  in the C\&L inventory (by number of molecules matched).}",
        r"\label{tab:clinventory_hg_top_groups}",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"\textbf{Group name} & \textbf{Molecules} & \textbf{\%\ of matches} \\",
        r"\midrule",
    ]
    for name, cnt in items:
        lines.append(f"  {esc(str(name))} & {fmt_int(cnt)} & {pct(cnt, total)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Narrative text
# ---------------------------------------------------------------------------

def table_speedup_pfasgroups(hg: dict, pfasgroups: dict, cmp: dict) -> str:
    """Speed-gain table: HalogenGroups vs PFASGroups (F-only) per bracket."""
    pfg_br = cmp.get("pfasgroups_timing_by_bracket", {})
    if not pfg_br:
        # Reconstruct from the individual JSON files
        pfg_br = {}
        for b in BRACKETS:
            hb = hg.get("timing_by_atom_bracket", {}).get(b, {})
            pb = pfasgroups.get("timing_by_atom_bracket", {}).get(b, {})
            if hb or pb:
                pfg_br[b] = {"hg": hb, "pfasgroups": pb}

    brackets = [b for b in BRACKETS if b in pfg_br]
    if not brackets:
        return ""

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Speed gain from restricting HalogenGroups to fluorine-only",
        r"  patterns (\texttt{halogens=`F'}, PFASGroups mode) vs.\ full",
        r"  HalogenGroups on the F-containing molecules of the C\&L inventory.",
        r"  Speed-up factor $=$ HalogenGroups median / PFASGroups median;",
        r"  values $>1$ indicate PFASGroups is faster.}",
        r"\label{tab:clinventory_pfasgroups_speedup}",
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
        r"\textbf{Bracket} & \textbf{$n$}"
        r" & \textbf{HG median (ms)}"
        r" & \textbf{PFASGroups median (ms)}"
        r" & \textbf{Speed-up} \\",
        r"\midrule",
    ]
    for b in brackets:
        d    = pfg_br[b]
        n_b  = fmt_int(d.get("n") or (d.get("hg") or {}).get("n"))
        hg_d = d.get("hg") or {}
        pfg_d = d.get("pfasgroups") or {}
        hm   = hg_d.get("median_ms") or hg_d.get("median")
        pm   = pfg_d.get("median_ms") or pfg_d.get("median")
        r_v  = d.get("ratio_hg_over_pfasgroups")
        if not hm or not pm:
            continue
        if r_v is None:
            r_v = hm / pm if pm > 0 else float("nan")
        lines.append(
            f"  {esc(b)} & {n_b} & {fmt(hm, 3)} & {fmt(pm, 3)} & {fmt(r_v, 2)} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def table_pfasgroups_agreement(cmp: dict) -> str:
    """2x2 agreement table: PFASGroups (F-only) vs PFAS-Atlas."""
    am    = cmp.get("pfasgroups_agreement_matrix", {})
    if not am:
        return ""
    total = am.get("total", 0)
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Classification agreement between PFASGroups",
        r"  (\texttt{halogens=`F'}, F-group detection) and PFAS-Atlas",
        r"  (any PFAS classification) on fluorine-containing molecules.}",
        r"\label{tab:clinventory_pfasgroups_agreement}",
        r"\begin{tabular}{lrr}",
        r"\toprule",
        r"\textbf{Category} & \textbf{Molecules} & \textbf{\%} \\",
        r"\midrule",
        f"  Both classify as F/PFAS & {fmt_int(am.get('both_pfas'))} & {pct(am.get('both_pfas',0), total)} \\\\",
        f"  PFASGroups only (F detected) & {fmt_int(am.get('only_pfasgroups'))} & {pct(am.get('only_pfasgroups',0), total)} \\\\",
        f"  PFAS-Atlas only & {fmt_int(am.get('only_atlas'))} & {pct(am.get('only_atlas',0), total)} \\\\",
        f"  Neither classifies & {fmt_int(am.get('neither'))} & {pct(am.get('neither',0), total)} \\\\",
        r"\midrule",
        f"  \\textbf{{Total}} & \\textbf{{{fmt_int(total)}}} & 100.0\\% \\\\",
        r"\midrule",
        f"  Overall agreement & \\multicolumn{{2}}{{r}}{{{fmt(am.get('agreement_pct'), 1)}\\%}} \\\\",
        f"  Cohen's $\\kappa$ & \\multicolumn{{2}}{{r}}{{{fmt(am.get('cohen_kappa'), 3)}}} \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ]
    return "\n".join(lines)


def narrative_text(hg: dict, atlas: dict, cmp: dict,
                   pfasgroups: Optional[dict] = None) -> str:
    """Plain LaTeX narrative paragraph for the comparison section."""
    am    = cmp.get("agreement_matrix", {})
    to    = cmp.get("timing_overall", {})
    hg_t  = to.get("hg", {})
    at_t  = to.get("atlas", {})
    kappa = am.get("cohen_kappa", float("nan"))
    ratio = to.get("ratio_hg_over_atlas", float("nan"))

    n_hg_db   = hg.get("n_molecules_fetched", 0)
    n_at_db   = atlas.get("n_molecules_fetched", 0)
    n_common  = am.get("total", 0)
    n_both    = am.get("both_pfas", 0)
    n_only_hg = am.get("only_hg", 0)
    n_only_at = am.get("only_atlas", 0)
    n_neither = am.get("neither", 0)
    agree_pct = am.get("agreement_pct", 0)

    hg_f_pct  = hg.get("fluorine_group_rate_pct", 0)
    at_p_pct  = atlas.get("pfas_rate_pct", 0)

    kappa_str = fmt(kappa, 3)
    ratio_str = fmt(ratio, 2)
    hg_med    = fmt(hg_t.get("median"), 3) if hg_t else "--"
    at_med    = fmt(at_t.get("median"), 3) if at_t else "--"
    hg_p95    = fmt(hg_t.get("p95"),    3) if hg_t else "--"
    at_p95    = fmt(at_t.get("p95"),    3) if at_t else "--"

    # Optional PFASGroups paragraph
    pfg_para = ""
    if pfasgroups:
        pfg_n      = pfasgroups.get("n_molecules_fetched", 0)
        pfg_rate   = pfasgroups.get("fluorine_group_rate_pct", 0)
        pfg_pto    = cmp.get("pfasgroups_timing_overall", {})
        pfg_speedup = cmp.get("pfasgroups_speedup_overall")
        pfg_am     = cmp.get("pfasgroups_agreement_matrix", {})
        pfg_kappa  = pfg_am.get("cohen_kappa", float("nan"))
        pfg_agree  = pfg_am.get("agreement_pct", 0)
        pfg_paired = pfg_am.get("total", 0)
        pfg_s      = (pfg_pto.get("pfasgroups") or {}).get("median")
        hg_pfg_s   = (pfg_pto.get("hg") or {}).get("median")
        pfg_para = (
            f"\n\\paragraph{{PFASGroups (F-only) mode.}}\n"
            f"To quantify the speed gain from restricting group-matching to fluorine-only\n"
            f"patterns, we ran HalogenGroups with \\texttt{{halogens=`F'}} (PFASGroups mode)\n"
            f"on the \\num{{{pfg_n:,}}} fluorine-containing molecules.\n"
            f"PFASGroups identified a fluorine-containing functional group in\n"
            f"\\SI{{{pfg_rate:.1f}}}{{\\percent}} of these molecules.\n"
            f"The median per-molecule time dropped from {fmt(hg_pfg_s, 3)}\\,ms\n"
            f"(full HalogenGroups on the same set) to {fmt(pfg_s, 3)}\\,ms,\n"
            f"a speed-up of {fmt(pfg_speedup, 2)}$\\times$.\n"
            f"Agreement with PFAS-Atlas across \\num{{{pfg_paired:,}}} paired molecules\n"
            f"was \\SI{{{pfg_agree:.1f}}}{{\\percent}} "
            f"(Cohen's $\\kappa = {fmt(pfg_kappa, 3)}$;\n"
            f"see Table~\\ref{{tab:clinventory_pfasgroups_agreement}} and\n"
            f"Figure~\\ref{{fig:clinventory_three_way_agreement}})."
        )

    return f"""\
% ----  Auto-generated by generate_clinventory_latex.py  ----

\\subsection{{Comparison of HalogenGroups and PFAS-Atlas on the C\\&L Inventory}}
\\label{{subsec:clinventory_comparison}}

We applied both HalogenGroups and PFAS-Atlas to the C\\&L inventory
(clinventory database, \\num{{{n_hg_db:,}}} molecules for HalogenGroups,
\\num{{{n_at_db:,}}} for PFAS-Atlas).  The two classifiers were evaluated on a
common set of \\num{{{n_common:,}}} molecules for which both tools successfully
produced a result.

\\paragraph{{PFAS detection rates.}}
HalogenGroups identified a fluorine-containing functional group in
\\SI{{{hg_f_pct:.1f}}}{{\\percent}} of the processed molecules, while PFAS-Atlas
returned a PFAS classification in \\SI{{{at_p_pct:.1f}}}{{\\percent}}.
The two classifiers agreed on \\SI{{{agree_pct:.1f}}}{{\\percent}} of the common
molecules (Cohen's $\\kappa = {kappa_str}$).
Of the {n_common:,} paired molecules, {n_both:,} were classified as
fluorinated/PFAS by both tools ({pct(n_both, n_common)}), {n_only_hg:,} were
detected exclusively by HalogenGroups ({pct(n_only_hg, n_common)}), {n_only_at:,}
exclusively by PFAS-Atlas ({pct(n_only_at, n_common)}), and {n_neither:,} were
not classified as fluorinated by either tool ({pct(n_neither, n_common)}).

\\paragraph{{Timing performance.}}
HalogenGroups achieved a median per-molecule classification time of
{hg_med}\\,ms (95th percentile: {hg_p95}\\,ms), while PFAS-Atlas
required {at_med}\\,ms at the median ({at_p95}\\,ms at the 95th percentile).
The overall timing ratio (HalogenGroups / PFAS-Atlas) was {ratio_str}$\\times$,
indicating that
{"HalogenGroups was faster" if (not math.isnan(float(ratio_str.replace("--", "nan"))) and float(ratio_str.replace("--","nan")) < 1) else "PFAS-Atlas was faster"}.
Detailed breakdowns by molecular size and halogen composition are provided in
Tables~\\ref{{tab:clinventory_timing_bracket}} and
\\ref{{tab:clinventory_timing_halogen}}.{{pfg_para if pfasgroups else ''}}
% ---- end auto-generated text ----
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate LaTeX for clinventory comparison")
    p.add_argument("--hg-file",        type=Path, default=None)
    p.add_argument("--pfasgroups-file", type=Path, default=None,
                   help="Path to PFASGroups (F-only) JSON. Auto-detected if omitted.")
    p.add_argument("--atlas-file",     type=Path, default=None)
    p.add_argument("--cmp-file",       type=Path, default=None)
    return p.parse_args()


def main():
    args = parse_args()
    data_dir    = BENCHMARK_DIR / "data"
    reports_dir = BENCHMARK_DIR / "reports"
    reports_dir.mkdir(exist_ok=True, parents=True)

    hg_path    = args.hg_file    or latest_file(str(data_dir / "halogengroups_clinventory_*.json"))
    atlas_path = args.atlas_file or latest_file(str(data_dir / "pfasatlas_clinventory_*.json"))
    cmp_path   = args.cmp_file   or latest_file(str(data_dir / "clinventory_comparison_*.json"))
    pfg_path   = getattr(args, "pfasgroups_file", None) or latest_file(
        str(data_dir / "pfasgroups_clinventory_*.json")
    )

    missing = False
    for label, path in [("HalogenGroups", hg_path), ("PFAS-Atlas", atlas_path),
                         ("Comparison",   cmp_path)]:
        if not path or not Path(path).exists():
            print(f"ERROR: {label} file not found")
            missing = True
    if missing:
        sys.exit(1)

    print("=" * 70)
    print("Generating LaTeX for clinventory comparison")
    print("=" * 70)
    print(f"  HalogenGroups  : {hg_path}")
    if pfg_path and Path(pfg_path).exists():
        print(f"  PFASGroups (F) : {pfg_path}")
    else:
        print("  PFASGroups (F) : not found — skipping PFASGroups tables/text")
        pfg_path = None
    print(f"  PFAS-Atlas     : {atlas_path}")
    print(f"  Comparison     : {cmp_path}")

    with open(hg_path)    as f: hg    = json.load(f)
    with open(atlas_path) as f: atlas = json.load(f)
    with open(cmp_path)   as f: cmp   = json.load(f)
    pfasgroups = None
    if pfg_path:
        with open(pfg_path) as f:
            pfasgroups = json.load(f)

    # -----------------------------------------------------------------------
    # Tables file
    tables_parts = [
        f"% Clinventory comparison tables \u2014 generated {datetime.now().isoformat()}\n",
        "% Requires: booktabs, siunitx\n\n",
        table_overall_timing(hg, atlas),
        table_timing_by_bracket(hg, atlas),
        table_timing_by_halogen(hg, atlas),
        table_agreement(cmp),
        table_atlas_classes(atlas),
        table_hg_top_groups(cmp),
    ]
    if pfasgroups:
        spd = table_speedup_pfasgroups(hg, pfasgroups, cmp)
        if spd:
            tables_parts.append(spd)
        agr = table_pfasgroups_agreement(cmp)
        if agr:
            tables_parts.append(agr)
    tables_out = reports_dir / "clinventory_comparison_tables.tex"
    with open(tables_out, "w") as f:
        f.write("\n".join(tables_parts))
    print(f"\n  Tables  → {tables_out}")

    # -----------------------------------------------------------------------
    # Narrative text file
    text_out = reports_dir / "clinventory_comparison_text.tex"
    with open(text_out, "w") as f:
        f.write(narrative_text(hg, atlas, cmp, pfasgroups))
    print(f"  Narrative → {text_out}")

    # -----------------------------------------------------------------------
    # Figure-reference snippets
    pfg_figs = ""
    if pfasgroups:
        pfg_figs = (
            "\n\\begin{figure}[htbp]\n"
            "\\centering\n"
            "\\includegraphics[width=0.75\\textwidth]{imgs/clinventory_timing_box_three.pdf}\n"
            "\\caption{Box plots of per-molecule timing for HalogenGroups, PFASGroups\n"
            "  (F-only, \\texttt{halogens=`F'}), and PFAS-Atlas. Outliers omitted.}\n"
            "\\label{fig:clinventory_timing_box_three}\n"
            "\\end{figure}\n"
            "\n"
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            "\\includegraphics[width=0.9\\textwidth]{imgs/clinventory_timing_by_bracket_three.pdf}\n"
            "\\caption{Median classification time by heavy-atom-count bracket for\n"
            "  all three classifiers.}\n"
            "\\label{fig:clinventory_timing_bracket_three}\n"
            "\\end{figure}\n"
            "\n"
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            "\\includegraphics[width=0.9\\textwidth]{imgs/clinventory_f_only_speedup.pdf}\n"
            "\\caption{Speed-up factor from restricting HalogenGroups to F-only SMARTS\n"
            "  patterns (\\texttt{halogens=`F'}, PFASGroups mode) vs.\\ full HalogenGroups,\n"
            "  by heavy-atom-count bracket. Values above~1 indicate PFASGroups is faster.}\n"
            "\\label{fig:clinventory_f_only_speedup}\n"
            "\\end{figure}\n"
            "\n"
            "\\begin{figure}[htbp]\n"
            "\\centering\n"
            "\\includegraphics[width=0.95\\textwidth]{imgs/clinventory_three_way_agreement.pdf}\n"
            "\\caption{PFAS-detection agreement with PFAS-Atlas for HalogenGroups (left)\n"
            "  and PFASGroups F-only mode (right), showing Cohen's $\\kappa$ and agreement\n"
            "  percentage.}\n"
            "\\label{fig:clinventory_three_way_agreement}\n"
            "\\end{figure}\n"
        )

    fig_snippets = (
        "% ---------- Figure snippets (paste into your article) ----------\n"
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        "\\includegraphics[width=0.7\\textwidth]{imgs/clinventory_timing_box.pdf}\n"
        "\\caption{Box plots of per-molecule classification time (ms) for HalogenGroups\n"
        "  and PFAS-Atlas on the C\\&L inventory. Outliers omitted.}\n"
        "\\label{fig:clinventory_timing_box}\n"
        "\\end{figure}\n"
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        "\\includegraphics[width=0.85\\textwidth]{imgs/clinventory_timing_cdf.pdf}\n"
        "\\caption{Empirical CDF of per-molecule classification time (log scale).}\n"
        "\\label{fig:clinventory_timing_cdf}\n"
        "\\end{figure}\n"
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        "\\includegraphics[width=0.9\\textwidth]{imgs/clinventory_timing_by_bracket.pdf}\n"
        "\\caption{Median classification time by heavy-atom-count bracket.}\n"
        "\\label{fig:clinventory_timing_bracket}\n"
        "\\end{figure}\n"
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        "\\includegraphics[width=0.9\\textwidth]{imgs/clinventory_classification_agreement.pdf}\n"
        "\\caption{Classification agreement between HalogenGroups and PFAS-Atlas.}\n"
        "\\label{fig:clinventory_agreement}\n"
        "\\end{figure}\n"
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        "\\includegraphics[width=0.75\\textwidth]{imgs/clinventory_atlas_class_distribution.pdf}\n"
        "\\caption{Distribution of PFAS-Atlas primary classifications.}\n"
        "\\label{fig:clinventory_atlas_dist}\n"
        "\\end{figure}\n"
        "\n"
        "\\begin{figure}[htbp]\n"
        "\\centering\n"
        "\\includegraphics[width=0.75\\textwidth]{imgs/clinventory_hg_top_groups.pdf}\n"
        "\\caption{Top-20 most frequently detected HalogenGroups functional groups.}\n"
        "\\label{fig:clinventory_hg_groups}\n"
        "\\end{figure}\n"
        + pfg_figs
    )
    figs_out = reports_dir / "clinventory_comparison_figures.tex"
    with open(figs_out, "w") as f:
        f.write(fig_snippets)
    print(f"  Figures → {figs_out}")

    print("\nDone.")
    print(f"\nOutput directory: {reports_dir}")
    print("Files generated:")
    print(f"  {tables_out.name}")
    print(f"  {text_out.name}")
    print(f"  {figs_out.name}")


if __name__ == "__main__":
    main()
