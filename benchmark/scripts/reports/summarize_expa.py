#!/usr/bin/env python3
"""
Summarize Experiment A results from toxcast_comparison_results.csv.

Two tables are produced:

1. GLOBAL summary  (expa_summary.csv / .md)
   For each (feature_set, model): mean +/- SD aggregated over all folds,
   repeats, AND endpoints.

2. PER-ENDPOINT summary  (expa_per_endpoint_summary.csv / .md)
   For each (endpoint, model): one row per feature_set showing
   mean +/- SD across folds & repeats, plus a "best_fs" column
   highlighting the winner per metric.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
METRICS  = ["roc_auc", "avg_prec", "mcc", "bal_acc"]
METRIC_LABELS = {
    "roc_auc":  "ROC-AUC",
    "avg_prec": "Avg. Precision",
    "mcc":      "Matthews Correlation Coefficient",
    "bal_acc":  "Balanced Accuracy",
}

# Feature-set display order: TxP baseline first, then PFG variants by complexity
FS_ORDER = [
    "TxP_PFAS",
    "TxP_PFAS+g51g52_total",
    "PFG_binary",
    "PFG_binary+mol",
    "PFG_EGR",
    "PFG_EGR+mol",
    "PFG_EGR+branch+mol",
    "PFG_EGR+spacer+mol",
    "PFG_EGR+ring+mol",
    "PFG_EGR+spacer+ring+mol",
    "PFG_EGR_4X",
    "PFG_EGR_4X+mol",
    "PFG_total_component",
    "PFG_total_component+mol",
    "PFG_min_dist_to_barycenter",
    "PFG_min_dist_to_barycenter+mol",
]


def _sort_fs(series: pd.Index) -> list[str]:
    present = set(series)
    ordered = [fs for fs in FS_ORDER if fs in present]
    ordered += sorted(present - set(ordered))
    return ordered


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Mean +/- SD across all folds*repeats*endpoints per (feature_set, model)."""
    rows = []
    for (fs, model), grp in df.groupby(["feature_set", "model"]):
        row: dict = {"feature_set": fs, "model": model}
        for m in METRICS:
            vals = grp[m].dropna()
            row[f"{m}_mean"] = vals.mean()
            row[f"{m}_std"]  = vals.std()
        rows.append(row)
    return pd.DataFrame(rows)


def print_table(summary: pd.DataFrame, model: str) -> None:
    sub = summary[summary["model"] == model].copy()
    fs_order = _sort_fs(sub["feature_set"])
    sub = sub.set_index("feature_set").loc[fs_order]
    best = {m: sub[f"{m}_mean"].idxmax() for m in METRICS}

    col_w = 22
    header = f"{'Feature set':<28}" + "".join(
        f"{METRIC_LABELS[m]:>{col_w}}" for m in METRICS
    )
    sep = "-" * len(header)
    print(f"\n{'='*len(header)}")
    print(f"  Experiment A  |  {model}")
    print(f"{'='*len(header)}")
    print(header)
    print(sep)
    for fs, row in sub.iterrows():
        cells = []
        for m in METRICS:
            val = f"{row[f'{m}_mean']:.4f} +/- {row[f'{m}_std']:.4f}"
            marker = " *" if best[m] == fs else "  "
            cells.append(f"{val + marker:>{col_w}}")
        print(f"{fs:<28}" + "".join(cells))
    print(sep)
    print("  * best mean for that metric")


def _md_table(summary: pd.DataFrame, model: str) -> str:
    """GitHub-flavoured Markdown table for one model (global summary)."""
    sub = summary[summary["model"] == model].copy()
    fs_order = _sort_fs(sub["feature_set"])
    sub = sub.set_index("feature_set").loc[fs_order]
    best = {m: sub[f"{m}_mean"].idxmax() for m in METRICS}

    header = "| Feature set | " + " | ".join(METRIC_LABELS[m] for m in METRICS) + " |"
    sep    = "|---|" + "|".join(["---"] * len(METRICS)) + "|"
    lines  = [f"### {model}\n", header, sep]
    for fs, row in sub.iterrows():
        cells = []
        for m in METRICS:
            val = f"{row[f'{m}_mean']:.4f} +/- {row[f'{m}_std']:.4f}"
            if best[m] == fs:
                val = f"**{val}**"
            cells.append(val)
        lines.append("| " + fs + " | " + " | ".join(cells) + " |")
    lines.append("\n_Bold = best mean for that metric._")
    return "\n".join(lines)


def build_per_endpoint(df: pd.DataFrame) -> pd.DataFrame:
    """Mean +/- SD across folds*repeats per (endpoint, feature_set, model)."""
    rows = []
    for (ep, fs, model), grp in df.groupby(["endpoint", "feature_set", "model"]):
        row: dict = {"endpoint": ep, "feature_set": fs, "model": model}
        for m in METRICS:
            vals = grp[m].dropna()
            row[f"{m}_mean"] = vals.mean()
            row[f"{m}_std"]  = vals.std()
        rows.append(row)
    return pd.DataFrame(rows)


def _md_per_endpoint(per_ep: pd.DataFrame, model: str) -> str:
    """
    One Markdown sub-table per endpoint.
    Rows = feature_sets, columns = all metrics (mean +/- SD).
    Best feature set per metric is bolded.
    """
    sub = per_ep[per_ep["model"] == model].copy()
    endpoints = sorted(sub["endpoint"].unique())
    fs_order  = _sort_fs(sub["feature_set"])

    sections = [f"### {model}\n"]
    for ep in endpoints:
        ep_sub = sub[sub["endpoint"] == ep].set_index("feature_set").reindex(fs_order)
        best = {m: ep_sub[f"{m}_mean"].idxmax() for m in METRICS}

        sections.append(f"#### {ep}\n")
        sections.append("| Feature set | " + " | ".join(METRIC_LABELS[m] for m in METRICS) + " |")
        sections.append("|---|" + "|".join(["---:"] * len(METRICS)) + "|")
        for fs in fs_order:
            row = ep_sub.loc[fs]
            cells = []
            for m in METRICS:
                mean_v = row[f"{m}_mean"]
                std_v  = row[f"{m}_std"]
                if pd.isna(mean_v):
                    cell = "-"
                else:
                    cell = f"{mean_v:.4f} +/- {std_v:.4f}"
                    if best[m] == fs:
                        cell = f"**{cell}**"
                cells.append(cell)
            sections.append("| " + fs + " | " + " | ".join(cells) + " |")
        sections.append("\n_Bold = best mean for that metric._\n")

    return "\n".join(sections)


def main() -> None:
    raw = pd.read_csv(DATA_DIR / "toxcast_comparison_results.csv")
    df  = raw[raw["experiment"] == "Exp A"].copy()

    print(f"Exp A rows: {len(df)}")
    print(f"Feature sets ({len(df['feature_set'].unique())}): {sorted(df['feature_set'].unique())}")
    print(f"Endpoints   ({len(df['endpoint'].unique())}): {sorted(df['endpoint'].unique())}")

    # ------------------------------------------------------------------ #
    # 1. Global summary (aggregated over endpoints)                        #
    # ------------------------------------------------------------------ #
    summary = build_summary(df)

    for model in sorted(summary["model"].unique()):
        print_table(summary, model)

    avg_rows = []
    for fs, grp in summary.groupby("feature_set"):
        row: dict = {"feature_set": fs, "model": "both"}
        for m in METRICS:
            row[f"{m}_mean"] = grp[f"{m}_mean"].mean()
            row[f"{m}_std"]  = grp[f"{m}_std"].mean()
        avg_rows.append(row)
    avg = pd.DataFrame(avg_rows)
    summary_all = pd.concat([summary, avg], ignore_index=True)
    print_table(summary_all, "both")

    out_csv = DATA_DIR / "expa_summary.csv"
    summary_all.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"\n[saved] {out_csv.relative_to(DATA_DIR.parent)}")

    md_sections = [
        "# Experiment A — Feature Set Summary (global)\n",
        "> Mean +/- SD across all CV folds x repeats x endpoints.\n",
    ]
    for model in sorted(summary["model"].unique()) + ["both"]:
        md_sections.append(_md_table(summary_all, model))
        md_sections.append("")
    out_md = DATA_DIR / "expa_summary.md"
    out_md.write_text("\n".join(md_sections), encoding="utf-8")
    print(f"[saved] {out_md.relative_to(DATA_DIR.parent)}")

    # ------------------------------------------------------------------ #
    # 2. Per-endpoint summary                                              #
    # ------------------------------------------------------------------ #
    per_ep = build_per_endpoint(df)

    # Add "both" model (average of GB and RF means)
    avg_ep_rows = []
    for (ep, fs), grp in per_ep.groupby(["endpoint", "feature_set"]):
        row2: dict = {"endpoint": ep, "feature_set": fs, "model": "both"}
        for m in METRICS:
            row2[f"{m}_mean"] = grp[f"{m}_mean"].mean()
            row2[f"{m}_std"]  = grp[f"{m}_std"].mean()
        avg_ep_rows.append(row2)
    per_ep_all = pd.concat([per_ep, pd.DataFrame(avg_ep_rows)], ignore_index=True)

    out_ep_csv = DATA_DIR / "expa_per_endpoint_summary.csv"
    per_ep_all.to_csv(out_ep_csv, index=False, float_format="%.6f")
    print(f"[saved] {out_ep_csv.relative_to(DATA_DIR.parent)}")

    md_ep_sections = [
        "# Experiment A — Performance per Endpoint\n",
        "> Each sub-table = one endpoint. Rows = feature sets.",
        "> Cells = mean +/- SD across CV folds & repeats. Bold = best mean for that metric.\n",
    ]
    for model in sorted(per_ep["model"].unique()) + ["both"]:
        md_ep_sections.append(_md_per_endpoint(per_ep_all, model))
        md_ep_sections.append("")
    out_ep_md = DATA_DIR / "expa_per_endpoint_summary.md"
    out_ep_md.write_text("\n".join(md_ep_sections), encoding="utf-8")
    print(f"[saved] {out_ep_md.relative_to(DATA_DIR.parent)}")


if __name__ == "__main__":
    main()
