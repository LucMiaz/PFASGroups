#!/usr/bin/env python3
"""
Experiment A  --  filtered analysis and bar-plot figures.

Steps
-----
1. Identify "reasonably predictable" endpoints:
   max ROC-AUC across all feature_sets (both models avg) >= ROC_THRESHOLD.

2. For those endpoints report the top-3 feature_sets (by mean ROC-AUC, both
   models averaged) with all performance statistics incl. SD.
   Saved to:
       data/expa_top3_per_endpoint.csv
       data/expa_top3_per_endpoint.md

3. Collect the union of feature_sets that appear in the top-3 for at least
   one endpoint.  Draw bar plots using only those feature_sets and only
   the predictable endpoints.
   - One figure per model (GradientBoosting / RandomForest)
   - One subplot per endpoint (3-col grid)
   - Horizontal bars = feature sets, length = ROC-AUC mean, error bar = SD
   - Gold star next to 1st-best feature set label, silver next to 2nd-best
   Saved to:
       imgs/expa_barplot_GradientBoosting.{png,pdf}
       imgs/expa_barplot_RandomForest.{png,pdf}

Usage
-----
    conda activate chem
    cd benchmark
    python scripts/plot_expa_filtered.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
ROC_THRESHOLD = 0.72
N_PLOTS      = 4   # number of endpoint subplots (top by best ROC-AUC)
N_FS         = 6   # number of feature sets shown (top by mean ROC-AUC; TxP_PFAS always included)
TOP_N_REPORT = 3   # feature sets per endpoint in the summary CSV/MD
PRIMARY      = "roc_auc"

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
IMGS_DIR   = SCRIPT_DIR.parents[1] / "imgs"
IMGS_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Font
# ---------------------------------------------------------------------------
_ubuntu_ok = any("Ubuntu" in f.name for f in fm.fontManager.ttflist)
_FONT = "Ubuntu" if _ubuntu_ok else "sans-serif"
if not _ubuntu_ok:
    print("  [warning] Ubuntu font not found – falling back to sans-serif")
plt.rcParams.update({"font.family": _FONT, "font.size": 9})

METRICS = ["roc_auc", "avg_prec", "mcc", "bal_acc"]
METRIC_LABELS = {
    "roc_auc":  "ROC-AUC",
    "avg_prec": "Avg. Precision",
    "mcc":      "MCC",
    "bal_acc":  "Bal. Accuracy",
}
# Color scheme from PFASGroups/data/color_scheme.yaml
METRIC_COLORS = {
    "roc_auc":  "#EB935C",   # orange-red
    "avg_prec": "#759ED1",   # blue
    "mcc":      "#BE6A9D",   # magenta
    "bal_acc":  "#8B61A8",   # dark purple
}

# Feature-set display order
FS_ORDER = [
    "TxP_PFAS",
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

EXCLUDED_FS = {"TxP_PFAS+g51g52_total"}


def _sort_fs(fsets):
    present = set(fsets) - EXCLUDED_FS
    ordered = [fs for fs in FS_ORDER if fs in present]
    ordered += sorted(present - set(ordered))
    return ordered


def _md_top2(top2_df: pd.DataFrame) -> str:
    lines = [
        f"# Experiment A -- Top-{TOP_N_REPORT} Feature Sets per Endpoint (summary)\n",
        "> Ranked by mean ROC-AUC (both models averaged). All metrics shown as mean +/- SD.\n",
    ]
    for ep in sorted(top2_df["endpoint"].unique()):
        sub = top2_df[top2_df["endpoint"] == ep].reset_index(drop=True)
        lines.append(f"## {ep}\n")
        lines.append("| Rank | Feature set | ROC-AUC | Avg. Precision | MCC | Bal. Accuracy |")
        lines.append("|:---:|---|---:|---:|---:|---:|")
        for _, row in sub.iterrows():
            rank_sym = {1: "**1st**", 2: "2nd", 3: "3rd"}.get(row["rank"], str(row["rank"]))
            cells = [f"{row[f'{m}_mean']:.4f} +/- {row[f'{m}_std']:.4f}" for m in METRICS]
            lines.append("| " + rank_sym + " | " + row["feature_set"] + " | " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines)


def _draw_barplot_figure(per_ep_model: pd.DataFrame, fs_filtered: list,
                         good_eps: list, model: str,
                         ep_pos_rate: dict | None = None) -> plt.Figure:
    """
    One subplot per endpoint.  For each feature set, four grouped horizontal
    bars (one per metric: ROC-AUC, Avg. Precision, MCC, Bal. Accuracy).
    Feature-set order is FIXED (same in every subplot, from FS_ORDER).
    Gold star = best ROC-AUC label; silver star = 2nd-best ROC-AUC label.
    """
    n_ep    = len(good_eps)
    n_cols  = 3
    n_rows  = math.ceil(n_ep / n_cols)
    n_fs    = len(fs_filtered)
    n_m     = len(METRICS)

    # Bar geometry (data-coordinate units)
    bar_h   = 0.15                    # height of one metric bar
    group_h = n_m * bar_h + 0.14     # total height per feature-set group

    # Display order: first in fs_filtered -> top of chart -> highest y index
    fs_display = list(reversed(fs_filtered))   # index 0 at bottom
    y_centers  = [i * group_h for i in range(n_fs)]

    # Figure height: scale with number of feature sets
    sp_h  = max(4.0, n_fs * group_h * 0.58)  # inches per subplot row
    fig_h = n_rows * sp_h + 1.8              # + title + legend

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 6.5, fig_h),
        constrained_layout=True,
    )
    axes_flat = axes.ravel()

    fig.suptitle(
        f"Experiment A  --  {model}\n"
        "ROC-AUC / Avg. Precision / MCC / Bal. Accuracy"
        "  (predictable endpoints, top feature sets)",
        fontsize=12, fontweight="bold",
    )

    for ax_idx, ep in enumerate(sorted(good_eps)):
        ax = axes_flat[ax_idx]
        ep_df = per_ep_model[per_ep_model["endpoint"] == ep].copy()

        # Determine gold / silver by ROC-AUC
        roc_map = {
            row["feature_set"]: row["roc_auc_mean"]
            for _, row in ep_df.iterrows()
            if row["feature_set"] in set(fs_filtered)
        }
        ranked    = sorted(roc_map, key=lambda f: roc_map[f], reverse=True)
        gold_fs   = ranked[0] if len(ranked) >= 1 else None
        silver_fs = ranked[1] if len(ranked) >= 2 else None

        all_vals = []
        for i, fs in enumerate(fs_display):
            y_c = y_centers[i]
            fd  = ep_df[ep_df["feature_set"] == fs]
            if fd.empty:
                continue
            r = fd.iloc[0]
            for j, metric in enumerate(METRICS):
                y_pos  = y_c + (j - (n_m - 1) / 2.0) * bar_h
                mean_v = r[f"{metric}_mean"]
                std_v  = r[f"{metric}_std"]
                all_vals.append(mean_v)
                ax.barh(
                    y_pos, mean_v, xerr=std_v,
                    height=bar_h,
                    color=METRIC_COLORS[metric],
                    edgecolor="white", linewidth=0.3,
                    error_kw=dict(elinewidth=0.7, capsize=1.5,
                                  capthick=0.7, ecolor="#333"),
                )

        # y-tick labels with stars
        ylabels = []
        for fs in fs_display:
            if fs == gold_fs:
                ylabels.append(f"\u2605 {fs}")
            elif fs == silver_fs:
                ylabels.append(f"\u2605 {fs}")
            else:
                ylabels.append(f"  {fs}")

        ax.set_yticks(y_centers)
        ax.set_yticklabels(ylabels, fontsize=7.0)

        for text_obj, fs in zip(ax.get_yticklabels(), fs_display):
            if fs == gold_fs:
                text_obj.set_color("#B8860B")   # dark gold
            elif fs == silver_fs:
                text_obj.set_color("#808080")   # silver

        # x-axis range: always start at 0
        xmax = min(1.02,  (max(all_vals) if all_vals else 1.0) + 0.08)
        ax.axvline(0.5, color="#aaa", linewidth=0.7, linestyle="--", zorder=0)
        ax.set_xlabel("Metric value", fontsize=8)
        pos_str = ""
        if ep_pos_rate and ep in ep_pos_rate:
            pos_str = f"  [pos: {ep_pos_rate[ep]*100:.1f}%]"
        ax.set_title(ep + pos_str, fontsize=9, fontweight="bold", pad=4)
        ax.set_xlim(0.0, xmax)
        ax.set_ylim(-group_h * 0.7, (n_fs - 1) * group_h + group_h * 0.7)
        ax.grid(axis="x", linewidth=0.3, alpha=0.5)
        ax.spines[["top", "right"]].set_visible(False)

    # Hide unused subplots
    for ax_idx in range(n_ep, n_rows * n_cols):
        axes_flat[ax_idx].set_visible(False)

    # Legend: metrics
    legend_handles = [
        Patch(color=METRIC_COLORS[m], label=METRIC_LABELS[m])
        for m in METRICS
    ]
    fig.legend(
        handles=legend_handles,
        title="Metric",
        loc="lower center",
        ncol=len(METRICS),
        fontsize=9,
        title_fontsize=10,
        bbox_to_anchor=(0.5, -0.03),
        frameon=True,
    )
    return fig


def main(dataset: str = "toxcast") -> None:
    suffix = f"_{dataset}" if dataset != "toxcast" else ""

    # ------------------------------------------------------------------ #
    # Load comparison results                                              #
    # ------------------------------------------------------------------ #
    raw = pd.read_csv(DATA_DIR / f"{dataset}_comparison_results.csv")
    df  = raw[raw["experiment"] == "Exp A"].copy()

    # ------------------------------------------------------------------ #
    # Step 1: identify predictable endpoints                               #
    # Derive ep_max directly from comparison results (avg over folds &     #
    # models, then max over feature sets).                                 #
    # ------------------------------------------------------------------ #
    per_ep_summary = DATA_DIR / f"{dataset}_per_endpoint_summary.csv"
    if not per_ep_summary.exists():
        per_ep_summary = DATA_DIR / "expa_per_endpoint_summary.csv"
        if per_ep_summary.exists():
            print(f"  [info] Using fallback per-endpoint summary: {per_ep_summary.name}")

    if per_ep_summary.exists():
        per_ep = pd.read_csv(per_ep_summary)
        both   = per_ep[per_ep["model"] == "both"].copy()
        ep_max = both.groupby("endpoint")[f"{PRIMARY}_mean"].max()
    else:
        # Derive ep_max inline from comparison results
        ep_max = (
            df.groupby(["endpoint", "feature_set", "model"])[PRIMARY]
            .mean()
            .groupby(["endpoint", "feature_set"])
            .mean()
            .reset_index()
            .groupby("endpoint")[PRIMARY]
            .max()
        )

    good_eps = sorted(ep_max.sort_values(ascending=False).head(N_PLOTS).index.tolist())

    print(f"Top {N_PLOTS} endpoints by best ROC-AUC:")
    for ep in sorted(good_eps, key=lambda e: ep_max[e], reverse=True):
        print(f"  {ep:35s}  best ROC-AUC = {ep_max[ep]:.4f}")
    print(f"  ({len(good_eps)} / {len(ep_max)} endpoints)")

    # ------------------------------------------------------------------ #
    # Step 2: top-2 feature sets per predictable endpoint                  #
    # ------------------------------------------------------------------ #

    agg = (
        df.groupby(["endpoint", "feature_set", "model"])[METRICS]
        .agg(["mean", "std"])
        .reset_index()
    )
    agg.columns = (
        ["endpoint", "feature_set", "model"]
        + [f"{m}_{s}" for m in METRICS for s in ["mean", "std"]]
    )

    avg_rows = []
    for (ep, fs), grp in agg.groupby(["endpoint", "feature_set"]):
        row = {"endpoint": ep, "feature_set": fs, "model": "both"}
        for m in METRICS:
            row[f"{m}_mean"] = grp[f"{m}_mean"].mean()
            row[f"{m}_std"]  = grp[f"{m}_std"].mean()
        avg_rows.append(row)
    agg_both = pd.DataFrame(avg_rows)

    report_rows = []
    for ep in good_eps:
        ep_both = agg_both[agg_both["endpoint"] == ep].copy()
        ep_both = ep_both[~ep_both["feature_set"].isin(EXCLUDED_FS)]
        ep_both = ep_both.sort_values(f"{PRIMARY}_mean", ascending=False).reset_index(drop=True)
        top_fsets = ep_both.head(TOP_N_REPORT)["feature_set"].tolist()
        for rank, fs in enumerate(top_fsets, start=1):
            for model in ["both", "GradientBoosting", "RandomForest"]:
                src = agg_both if model == "both" else agg
                row_data = src[(src["endpoint"] == ep) & (src["feature_set"] == fs)]
                if row_data.empty:
                    continue
                r = row_data.iloc[0].to_dict()
                r["rank"] = rank
                report_rows.append(r)

    report_df = pd.DataFrame(report_rows)
    col_order = (["endpoint", "rank", "feature_set", "model"]
                 + [f"{m}_{s}" for m in METRICS for s in ["mean", "std"]])
    report_df = report_df[col_order].sort_values(["endpoint", "rank", "model"])

    out_csv = DATA_DIR / f"expa_top_per_endpoint{suffix}.csv"
    report_df.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"\n[saved] {out_csv.relative_to(DATA_DIR.parent)}")

    out_md = DATA_DIR / f"expa_top_per_endpoint{suffix}.md"
    out_md.write_text(_md_top2(report_df[report_df["model"] == "both"]), encoding="utf-8")
    print(f"[saved] {out_md.relative_to(DATA_DIR.parent)}")

    # ------------------------------------------------------------------ #
    # Step 3: top N_FS feature sets and bar-plot figures                   #
    # ------------------------------------------------------------------ #
    fs_roc = (
        agg_both[agg_both["endpoint"].isin(good_eps)]
        .groupby("feature_set")[f"{PRIMARY}_mean"]
        .mean()
    )
    fs_roc = fs_roc[~fs_roc.index.isin(EXCLUDED_FS)].sort_values(ascending=False)
    top_fsets = fs_roc.head(N_FS).index.tolist()
    # Always include TxP_PFAS (swap out lowest-ranked if needed)
    if "TxP_PFAS" not in top_fsets:
        top_fsets[-1] = "TxP_PFAS"
    fs_filtered = _sort_fs(top_fsets)

    print(f"\nTop {N_FS} feature sets by mean ROC-AUC across selected endpoints (TxP_PFAS always included):")
    for fs in fs_filtered:
        val = fs_roc.get(fs, float("nan"))
        print(f"  {fs:35s}  mean ROC-AUC = {val:.4f}")

    # pre-compute per-endpoint aggregates for plotting (mean + std per fold*repeat)
    agg_plot = agg[agg["endpoint"].isin(good_eps) & agg["feature_set"].isin(fs_filtered)].copy()

    # class imbalance: mean positive rate per endpoint
    ep_pos_rate = (
        df[df["endpoint"].isin(good_eps)]
        .groupby("endpoint")["pos_rate"]
        .mean()
        .to_dict()
    )

    for model in ["GradientBoosting", "RandomForest"]:
        per_ep_model = agg_plot[agg_plot["model"] == model].copy()
        fig = _draw_barplot_figure(per_ep_model, fs_filtered, good_eps, model,
                                   ep_pos_rate=ep_pos_rate)

        stem = f"expa_barplot_{model}{suffix}"
        for ext in ("png", "pdf"):
            out = IMGS_DIR / f"{stem}.{ext}"
            fig.savefig(out, dpi=150, bbox_inches="tight")
            print(f"[saved] {out.relative_to(IMGS_DIR.parent)}")
        plt.close(fig)

    print("\nDone.")


if __name__ == "__main__":
    import argparse
    _parser = argparse.ArgumentParser(
        description="Experiment A filtered bar-plot figures."
    )
    _parser.add_argument(
        "--dataset", "-d", default="toxcast",
        help="Dataset prefix used in input CSV and output filenames (default: 'toxcast').",
    )
    _args = _parser.parse_args()
    main(_args.dataset)

