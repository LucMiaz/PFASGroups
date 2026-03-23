#!/usr/bin/env python3
"""
Plot three metrics (balanced accuracy, average precision, ROC AUC) for each
feature set across all ToxCast endpoints, separately per experiment and model.
Uses grouped bar charts (endpoints are unordered categories).

Two sets of figures are produced per experiment × model:
  - Full:       all feature sets present in that experiment
  - Restricted: a curated subset (see RESTRICTED_FSETS_A / _B below)

Layout (one PNG per experiment × model × variant, 8 files total):
  3 subplots stacked vertically: Balanced Accuracy / Average Precision / ROC AUC.
  X-axis  – ToxCast endpoints (alphabetical).
  Bars    – one group per endpoint, one bar per feature_set.
            Error bars = ±1 std from nested-CV folds/repeats.
  Colours – 4-colour palette from color_scheme.yaml, cycling;
            hatch patterns distinguish feature sets sharing a colour.

Output files (benchmark/imgs/)
───────────────────────────────
  expA_GradientBoosting_metrics.png
  expA_RandomForest_metrics.png
  expB_GradientBoosting_metrics.png
  expB_RandomForest_metrics.png
  expA_GradientBoosting_metrics_restricted.png
  expA_RandomForest_metrics_restricted.png
  expB_GradientBoosting_metrics_restricted.png
  expB_RandomForest_metrics_restricted.png

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/plot_performance_by_endpoint.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parent / "data"
IMGS_DIR   = SCRIPT_DIR.parent / "imgs"
COLOR_YAML = SCRIPT_DIR.resolve().parents[1] / "PFASGroups" / "data" / "color_scheme.yaml"

IMGS_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
with open(COLOR_YAML) as fh:
    _cs = yaml.safe_load(fh)

PALETTE = _cs["colorscheme"]  # ["#E15D0B", "#306DBA", "#9D206C", "#51127C"]
# Hatch patterns to visually distinguish bars sharing the same colour
HATCHES = ["", "//", "xx", ".."]

_ubuntu_available = any("Ubuntu" in f.name for f in fm.fontManager.ttflist)
FONT_FAMILY = "Ubuntu" if _ubuntu_available else "sans-serif"
if not _ubuntu_available:
    print("  [warning] Ubuntu font not found in font cache – falling back to sans-serif")

sns.set_theme(style="whitegrid", font=FONT_FAMILY)
plt.rcParams.update({
    "font.family":       FONT_FAMILY,
    "axes.titlesize":    11,
    "axes.labelsize":    10,
    "xtick.labelsize":   8,
    "ytick.labelsize":   9,
    "legend.fontsize":   8,
    "legend.framealpha": 0.9,
    "figure.dpi":        150,
})


def _bar_style(idx: int) -> dict:
    """Return color / hatch for bar at rank *idx* (cycles over palette + hatches)."""
    n = len(PALETTE)
    return dict(
        color=PALETTE[idx % n],
        hatch=HATCHES[(idx // n) % len(HATCHES)],
        edgecolor="#333333",
        linewidth=0.4,
    )


# ---------------------------------------------------------------------------
# Curated feature-set subsets
# ---------------------------------------------------------------------------
RESTRICTED_FSETS_A = [
    "TxP_PFAS",
    "TxP_PFAS+g51g52_total",
    "PFG_binary+mol",
    "PFG_EGR+mol",
    "PFG_EGR+branch+mol",
    "PFG_EGR+spacer+mol",
    "PFG_EGR+ring+mol",
    "PFG_total_component+mol",
    "PFG_min_dist_to_barycenter+mol",
]

RESTRICTED_FSETS_B = [
    "Morgan",
    "ToxPrint+TxP_PFAS",
    "Morgan+PFG_EGR_4X+mol",
    "ToxPrint+PFG_EGR_4X+mol",
]

# ---------------------------------------------------------------------------
# Metrics to plot
# ---------------------------------------------------------------------------
METRICS = [
    ("bal_acc",  "Balanced Accuracy"),
    ("avg_prec", "Average Precision"),
    ("roc_auc",  "ROC AUC"),
]

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
raw = pd.read_csv(DATA_DIR / "toxcast_comparison_results.csv")

# Endpoint order: alphabetical
ep_order = sorted(raw["endpoint"].unique())

# Feature-set rank by mean balanced accuracy (ascending → best last in legend)
fs_rank = (
    raw.groupby("feature_set")["bal_acc"]
    .mean()
    .sort_values(ascending=True)
    .index.tolist()
)

EXPERIMENTS  = [("Exp A", "expA"), ("Exp B", "expB")]
MODELS       = ["GradientBoosting", "RandomForest"]
MODEL_LABELS = {
    "GradientBoosting": "Gradient Boosting",
    "RandomForest":     "Random Forest",
}


# ---------------------------------------------------------------------------
# Core figure builder
# ---------------------------------------------------------------------------
def make_bar_figure(
    df_model: pd.DataFrame,
    fsets: list[str],
    title: str,
) -> plt.Figure:
    """
    Return a figure with 3 stacked bar-chart subplots (one per metric).
    X-axis = endpoint, grouped bars = feature_set, error bars = ±1 std.
    """
    agg = (
        df_model.groupby(["feature_set", "endpoint"])[
            [m for m, _ in METRICS]
        ]
        .agg(["mean", "std"])
        .reset_index()
    )
    agg.columns = pd.Index(
        ["feature_set", "endpoint"]
        + [f"{m}_{stat}" for m, _ in METRICS for stat in ("mean", "std")]
    )

    n_ep  = len(ep_order)
    n_fs  = len(fsets)
    width = 0.82 / max(n_fs, 1)
    x     = np.arange(n_ep)

    fig_w = max(14, n_ep * max(0.7, 0.12 * n_fs + 0.3) + 4)
    fig, axes = plt.subplots(
        len(METRICS), 1,
        figsize=(fig_w, 3.8 * len(METRICS)),
        sharex=True,
        constrained_layout=True,
    )
    fig.suptitle(title, fontsize=13, fontweight="bold")

    for ax, (metric, metric_label) in zip(axes, METRICS):
        ax.axhline(0.5, color="grey", linewidth=0.8, linestyle=":", alpha=0.6, zorder=0)

        for idx, fs in enumerate(fsets):
            sub = (
                agg[agg["feature_set"] == fs]
                .set_index("endpoint")
                .reindex(ep_order)
            )
            means = sub[f"{metric}_mean"].values.astype(float)
            stds  = sub[f"{metric}_std"].values.astype(float)

            offset = (idx - (n_fs - 1) / 2) * width
            style  = _bar_style(idx)
            ax.bar(
                x + offset, means, width * 0.92,
                label=fs,
                yerr=stds,
                capsize=2,
                error_kw={"elinewidth": 0.9, "ecolor": "#111111", "alpha": 0.75},
                alpha=0.85,
                **style,
            )

        ax.set_ylabel(metric_label)
        ax.set_ylim(-0.02, 1.1)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(ep_order, rotation=45, ha="right", fontsize=8)
    axes[-1].legend(
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
        frameon=True,
        ncol=1,
        title="Feature set",
        title_fontsize=8,
    )
    return fig


# ---------------------------------------------------------------------------
# Generate all figures
# ---------------------------------------------------------------------------
for exp_label, exp_slug in EXPERIMENTS:
    df_exp = raw[raw["experiment"] == exp_label].copy()
    exp_fsets_all = [fs for fs in fs_rank if fs in df_exp["feature_set"].unique()]

    restricted_pool = RESTRICTED_FSETS_A if exp_label == "Exp A" else RESTRICTED_FSETS_B
    exp_fsets_restricted = [fs for fs in restricted_pool if fs in df_exp["feature_set"].unique()]

    for model in MODELS:
        df_m = df_exp[df_exp["model"] == model].copy()

        # ── Full set ──────────────────────────────────────────────────
        fig = make_bar_figure(
            df_m, exp_fsets_all,
            f"{exp_label}  —  {MODEL_LABELS[model]}",
        )
        out = IMGS_DIR / f"{exp_slug}_{model}_metrics.png"
        fig.savefig(out, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print(f"  Saved: {out.relative_to(IMGS_DIR.parent)}")

        # ── Restricted set ────────────────────────────────────────────
        if exp_fsets_restricted:
            fig = make_bar_figure(
                df_m, exp_fsets_restricted,
                f"{exp_label}  —  {MODEL_LABELS[model]}  (selected feature sets)",
            )
            out = IMGS_DIR / f"{exp_slug}_{model}_metrics_restricted.png"
            fig.savefig(out, bbox_inches="tight", dpi=150)
            plt.close(fig)
            print(f"  Saved: {out.relative_to(IMGS_DIR.parent)}")

print("Done.")
