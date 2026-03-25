#!/usr/bin/env python3
"""
Plot accuracy (balanced accuracy) and average precision for each feature_set
in Experiment A, broken down by individual endpoint, for each model.

Layout: 2 rows (GradientBoosting / RandomForest) × 2 columns (Balanced
Accuracy / Average Precision), each panel a heatmap of
feature_sets (rows) × endpoints (columns).

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/plot_expa_accuracy_precision.py
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
IMGS_DIR   = SCRIPT_DIR.parents[1] / "imgs"
IMGS_DIR.mkdir(exist_ok=True)

MODELS = ["GradientBoosting", "RandomForest"]
MODEL_LABELS = {"GradientBoosting": "Gradient Boosting", "RandomForest": "Random Forest"}
METRICS = [
    ("bal_acc",  "Balanced Accuracy",       "YlOrRd",  0.45, 0.80),
    ("avg_prec", "Average Precision (AP)",  "YlGnBu",  0.15, 0.65),
]


def main(dataset: str = "toxcast") -> None:
    suffix = f"_{dataset}" if dataset != "toxcast" else ""

    # -----------------------------------------------------------------------
    # Load raw CV results and filter for Exp A
    # -----------------------------------------------------------------------
    raw = pd.read_csv(DATA_DIR / f"{dataset}_comparison_results.csv")
    df  = raw[raw["experiment"] == "Exp A"].copy()

    # Aggregate across repeats and folds → mean per (feature_set, model, endpoint)
    per = (
        df.groupby(["feature_set", "model", "endpoint"])[["bal_acc", "avg_prec"]]
        .mean()
        .reset_index()
    )

    # Sort feature sets by overall mean balanced accuracy (ascending → best at top)
    fs_order = (
        per.groupby("feature_set")["bal_acc"]
        .mean()
        .sort_values(ascending=True)
        .index.tolist()
    )

    # Sort endpoints alphabetically
    ep_order = sorted(per["endpoint"].unique().tolist())

    # -----------------------------------------------------------------------
    # Build figure: 2 rows × 2 cols of heatmaps
    # -----------------------------------------------------------------------
    n_fs = len(fs_order)
    n_ep = len(ep_order)

    fig, axes = plt.subplots(
        2, 2,
        figsize=(max(14, n_ep * 0.85 + 3), max(8, n_fs * 0.55 + 3)),
        constrained_layout=True,
    )
    fig.suptitle("Experiment A — Per-endpoint Performance by Feature Set",
                 fontsize=13, fontweight="bold")

    for row_idx, model in enumerate(MODELS):
        sub = per[per["model"] == model]
        for col_idx, (metric, metric_label, cmap, vmin, vmax) in enumerate(METRICS):
            ax = axes[row_idx][col_idx]

            # Build matrix: rows=feature_sets (bottom to top), cols=endpoints
            mat = np.full((n_fs, n_ep), np.nan)
            for i, fs in enumerate(fs_order):
                for j, ep in enumerate(ep_order):
                    val = sub.loc[(sub["feature_set"] == fs) & (sub["endpoint"] == ep), metric]
                    if not val.empty:
                        mat[i, j] = val.values[0]

            im = ax.imshow(
                mat,
                aspect="auto",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )

            # Annotate cells
            for i in range(n_fs):
                for j in range(n_ep):
                    v = mat[i, j]
                    if not np.isnan(v):
                        norm_v = (v - vmin) / (vmax - vmin)
                        txt_color = "white" if norm_v > 0.65 else "black"
                        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                                fontsize=6.5, color=txt_color)

            # Axes labels and ticks
            ax.set_xticks(range(n_ep))
            ax.set_xticklabels(ep_order, rotation=40, ha="right", fontsize=8)
            ax.set_yticks(range(n_fs))
            ax.set_yticklabels(fs_order, fontsize=8)

            ax.set_title(f"{MODEL_LABELS[model]}  —  {metric_label}", fontsize=10, pad=6)

            cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
            cb.ax.tick_params(labelsize=7)
            cb.set_label(metric_label, fontsize=8)

    out_path = IMGS_DIR / f"expa_accuracy_precision{suffix}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved → {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    _parser = argparse.ArgumentParser(
        description="Heatmap of Exp A accuracy/precision by feature set."
    )
    _parser.add_argument(
        "--dataset", "-d", default="toxcast",
        help="Dataset prefix used in input CSV and output filenames (default: 'toxcast').",
    )
    _args = _parser.parse_args()
    main(_args.dataset)
