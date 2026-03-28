#!/usr/bin/env python3
"""
Bayesian correlated t-test comparison of PFASGroups embeddings vs TxP_PFAS.

For each endpoint in Experiment A, and for each metric, we use the Bayesian
correlated t-test (Benavoli et al., 2017, JMLR) to estimate the posterior
probability that the best PFASGroups feature set outperforms TxP_PFAS.

The test implementation is taken directly from the original authors' code:
https://github.com/BayesianTestsML/tutorial/blob/master/Python/bayesiantests.py

The variance uses the Nadeau-Bengio correction:
  var = s^2 * (1/n + 1/(nfolds - 1))
where n = total CV splits and nfolds = n / runs (runs = number of repetitions).

Reference
---------
Benavoli A, Corani G, Demsar J, Zaffalon M (2017). Time for a change:
a tutorial for comparing multiple classifiers through Bayesian analysis.
J. Machine Learning Research, 18(77), 1-36.

Outputs
-------
data/bayesian_comparison.csv   -- per-endpoint posterior probabilities
imgs/bayesian_comparison.pdf   -- figure: P(PFG > TxP_PFAS) per endpoint
imgs/bayesian_comparison.png
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch

# Use the original Benavoli et al. implementation bundled in this directory
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from bayesiantests import correlated_ttest, correlated_ttest_MC  # noqa: E402
try:
    from bayesiantests import hierarchical_MC_endpoints as _hierarchical_MC_endpoints
    _HIERARCHICAL_OK = True
except ImportError:
    _HIERARCHICAL_OK = False
from scipy.stats import gaussian_kde

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SCRIPT_DIR  = Path(__file__).resolve().parent
DATA_DIR    = SCRIPT_DIR.parents[1] / "data"
IMGS_DIR    = SCRIPT_DIR.parents[1] / "imgs"
IMGS_DIR.mkdir(exist_ok=True)

# Experiment A: 5 repeats × 3 folds = 15 splits total
N_REPEATS = 5

# Region of practical equivalence (ROPE): |delta| < 0.01 considered equivalent
ROPE = 0.01
METRICS = ["roc_auc", "avg_prec", "mcc", "bal_acc"]
METRIC_LABELS = {
    "roc_auc":  "ROC-AUC",
    "avg_prec": "Avg. Precision",
    "mcc":      "MCC",
    "bal_acc":  "Bal. Accuracy",
}
# Project colour scheme
METRIC_COLORS = {
    "roc_auc":  "#EB935C",
    "avg_prec": "#759ED1",
    "mcc":      "#BE6A9D",
    "bal_acc":  "#8B61A8",
}
BASELINE = "TxP_PFAS"

# ---------------------------------------------------------------------------
# Font
# ---------------------------------------------------------------------------
_ubuntu_ok = any("Ubuntu" in f.name for f in fm.fontManager.ttflist)
_FONT = "Ubuntu" if _ubuntu_ok else "sans-serif"
plt.rcParams.update({"font.family": _FONT, "font.size": 9})


# ---------------------------------------------------------------------------
# Posterior distribution figure
# ---------------------------------------------------------------------------

# Significance colours: gold / silver / bronze
_SIG_COLORS = {
    "***": "#B8860B",   # dark gold
    "**":  "#777777",   # silver-grey
    "*":   "#8B4513",   # bronze/brown
}


def _draw_posterior_figure(best_diffs: dict, ep_order: list,
                           best_rows: "pd.DataFrame" = None,
                           suffix: str = "") -> None:
    """
    For each (endpoint, metric) pair, plot the posterior KDE of mean differences
    (best PFASGroups − TxP_PFAS) with ROPE regions shaded.
    Layout: endpoints as rows, metrics as columns.
    If best_rows is provided, annotate each panel with the mean difference
    and significance stars (gold ***, silver **, bronze *).
    """
    n_eps = len(ep_order)
    n_met = len(METRICS)

    fig, axes = plt.subplots(
        n_eps, n_met,
        figsize=(n_met * 3.2, n_eps * 1.35),
        constrained_layout=True,
    )

    fig.suptitle(
        r"Posterior distribution of mean difference (best PFASGroups $-$ TxP\_PFAS)" + "\n"
        r"Orange lines: ROPE boundaries ($\pm$0.01).  "
        "Shading: red = TxP better, grey = equivalent, colour = PFASGroups better.",
        fontsize=14, fontweight="bold",
    )

    for col, metric in enumerate(METRICS):
        axes[0, col].set_title(METRIC_LABELS[metric], fontsize=14, fontweight="bold")

    for row, ep in enumerate(ep_order):
        axes[row, 0].set_ylabel(
            ep.replace("_", " "), fontsize=12, rotation=0, ha="right", labelpad=4,
        )
        for col, metric in enumerate(METRICS):
            ax  = axes[row, col]
            key = (ep, metric)

            if key not in best_diffs:
                ax.axis("off")
                continue

            diffs   = best_diffs[key]
            samples = correlated_ttest_MC(diffs, ROPE, runs=N_REPEATS, nsamples=40000)

            p1, p99 = np.percentile(samples, [0.5, 99.5])
            xs = np.linspace(min(p1, -5 * ROPE), max(p99, 5 * ROPE), 400)
            ys = gaussian_kde(samples)(xs)

            # Shaded regions
            ax.fill_between(xs, ys, where=(xs < -ROPE),
                            color="#E07070", alpha=1, linewidth=0)
            ax.fill_between(xs, ys, where=((xs >= -ROPE) & (xs <= ROPE)),
                            color="#C0C0C0", alpha=1, linewidth=0)
            ax.fill_between(xs, ys, where=(xs > ROPE),
                            color=METRIC_COLORS[metric], alpha=1, linewidth=0)

            ax.plot(xs, ys, color="#222", linewidth=0.9)
            ax.axvline(-ROPE, color="orange", linewidth=0.8)
            ax.axvline( ROPE, color="orange", linewidth=0.8)
            ax.axvline(0,     color="#888",   linewidth=0.5, linestyle="--")

            if row % 2 == 0:
                ax.axhspan(0, max(ys)+0.5, color="#f0f0f0", zorder=0, linewidth=0)

            # ── mean difference annotation ─────────────────────────────
            if best_rows is not None:
                sub = best_rows[
                    (best_rows["endpoint"] == ep) & (best_rows["metric"] == metric)
                ]
                if not sub.empty:
                    mean_d = float(sub.iloc[0]["mean_diff"])
                    p_pfg  = float(sub.iloc[0]["p_pfg_better"])
                    if p_pfg >= 0.999:
                        stars, sig_color = "***", _SIG_COLORS["***"]
                    elif p_pfg >= 0.99:
                        stars, sig_color = "**",  _SIG_COLORS["**"]
                    elif p_pfg >= 0.95:
                        stars, sig_color = "*",   _SIG_COLORS["*"]
                    else:
                        stars, sig_color = "",    "#444"
                    sign  = "+" if mean_d >= 0 else ""
                    label = f"{sign}{mean_d:.3f}{stars}\n ({p_pfg:.0%})"
                    fw    = "bold" if stars else "normal"
                    # Pin label to top-right corner of the axes
                    ax.text(
                        0.98, 0.97, label,
                        ha="right", va="top", fontsize=12,
                        fontweight=fw, color=sig_color,
                        transform=ax.transAxes
                    )
                    # Vertical line at the mean
                    ax.axvline(mean_d, color=sig_color, linewidth=0.8,
                               linestyle="-", alpha=1, zorder=4)

            ax.set_yticks([])
            ax.tick_params(axis="x", labelsize=12)
            ax.spines[["top", "right", "left"]].set_visible(False)

    for ext in ("png", "pdf"):
        out = IMGS_DIR / f"bayesian_posteriors{suffix}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"[saved] {out.relative_to(IMGS_DIR.parent)}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Hierarchical posterior figure
# ---------------------------------------------------------------------------

def _draw_hierarchical_posterior_figure(diffs_store: dict, ep_order: list,
                                        best_rows: "pd.DataFrame" = None,
                                        suffix: str = "") -> None:
    """
    Like _draw_posterior_figure, but uses the hierarchical Bayesian model
    (Benavoli et al. 2017) run jointly over all endpoints per metric.
    Each endpoint's posterior is shrunk toward the global hyperprior mean,
    which better accounts for overlapping training sets across CV folds.
    An extra bottom row shows the global hyperprior mean (delta0) distribution.

    Requires pystan v2.  Skipped gracefully if pystan is unavailable.
    """
    if not _HIERARCHICAL_OK:
        print("[skip] hierarchical figure: pystan / hierarchical_MC_endpoints not available")
        return

    n_eps = len(ep_order)
    n_met = len(METRICS)
    rho = 1.0 / 3  # 3-fold CV   ->  fold correlation = 1 / nfolds

    # n_eps rows for endpoints + 1 row for the global hyperprior mean
    fig, axes = plt.subplots(
        n_eps + 1, n_met,
        figsize=(n_met * 3.2, (n_eps + 1) * 1.35),
        constrained_layout=True,
    )

    fig.suptitle(
        r"Hierarchical Bayesian posterior: best PFASGroups $-$ TxP\_PFAS" + "\n"
        r"Orange lines: ROPE ($\pm$0.01).  "
        "Bottom row = global hyperprior mean $\\delta_0$ (shrinkage regularisation).",
        fontsize=14, fontweight="bold",
    )

    for col, metric in enumerate(METRICS):
        axes[0, col].set_title(METRIC_LABELS[metric], fontsize=14, fontweight="bold")

    for col, metric in enumerate(METRICS):
        # Build diff matrix: one row per endpoint (best feature set for that endpoint)
        diff_rows, ep_rows = [], []
        for ep in ep_order:
            best_fs = None
            if best_rows is not None:
                sub = best_rows[
                    (best_rows["endpoint"] == ep) & (best_rows["metric"] == metric)
                ]
                if not sub.empty:
                    best_fs = sub.iloc[0]["feature_set"]
            if best_fs is None:
                for key in diffs_store:
                    if key[0] == ep and key[1] == metric:
                        best_fs = key[2]
                        break
            if best_fs is None:
                continue
            key = (ep, metric, best_fs)
            if key in diffs_store:
                diff_rows.append(diffs_store[key])
                ep_rows.append(ep)

        if not diff_rows:
            for row in range(n_eps + 1):
                axes[row, col].axis("off")
            continue

        n_folds = min(len(d) for d in diff_rows)
        diff_matrix = np.array([d[:n_folds] for d in diff_rows], dtype=float)  # (q, n_folds)

        try:
            result = _hierarchical_MC_endpoints(diff_matrix, ROPE, rho)
        except Exception as exc:
            print(f"[hierarchical] {metric}: {exc}")
            for row in range(n_eps + 1):
                axes[row, col].axis("off")
            continue

        delta_ep = result["delta"]   # shape (n_mcmc, q)
        delta0   = result["delta0"]  # shape (n_mcmc,)

        # ── Per-endpoint rows ──────────────────────────────────────────
        for row, ep in enumerate(ep_order):
            ax = axes[row, col]
            if ep not in ep_rows:
                ax.axis("off")
                continue
            ep_idx = ep_rows.index(ep)
            samples = delta_ep[:, ep_idx]

            p1, p99 = np.percentile(samples, [0.5, 99.5])
            xs = np.linspace(min(p1, -5 * ROPE), max(p99, 5 * ROPE), 400)
            ys = gaussian_kde(samples)(xs)

            ax.fill_between(xs, ys, where=(xs < -ROPE),
                            color="#E07070", alpha=1, linewidth=0)
            ax.fill_between(xs, ys, where=((xs >= -ROPE) & (xs <= ROPE)),
                            color="#C0C0C0", alpha=1, linewidth=0)
            ax.fill_between(xs, ys, where=(xs > ROPE),
                            color=METRIC_COLORS[metric], alpha=1, linewidth=0)
            ax.plot(xs, ys, color="#222", linewidth=0.9)
            ax.axvline(-ROPE, color="orange", linewidth=0.8)
            ax.axvline( ROPE, color="orange", linewidth=0.8)
            ax.axvline(0,     color="#888",   linewidth=0.5, linestyle="--")

            if row % 2 == 0:
                ax.axhspan(0, max(ys) + 0.5, color="#f0f0f0", zorder=0, linewidth=0)

            mean_d   = float(np.mean(samples))
            p_better = float(np.mean(samples > ROPE))
            if p_better >= 0.999:
                stars, sig_color = "***", _SIG_COLORS["***"]
            elif p_better >= 0.99:
                stars, sig_color = "**",  _SIG_COLORS["**"]
            elif p_better >= 0.95:
                stars, sig_color = "*",   _SIG_COLORS["*"]
            else:
                stars, sig_color = "",    "#444"
            sign  = "+" if mean_d >= 0 else ""
            label = f"{sign}{mean_d:.3f}{stars}\n ({p_better:.0%})"
            fw    = "bold" if stars else "normal"
            ax.text(0.98, 0.97, label, ha="right", va="top", fontsize=12,
                    fontweight=fw, color=sig_color, transform=ax.transAxes)
            ax.axvline(mean_d, color=sig_color, linewidth=0.8, linestyle="-",
                       alpha=1, zorder=4)

            ax.set_yticks([])
            ax.tick_params(axis="x", labelsize=12)
            ax.spines[["top", "right", "left"]].set_visible(False)

        # ── Global delta0 row ──────────────────────────────────────────
        ax_g = axes[n_eps, col]
        p1, p99 = np.percentile(delta0, [0.5, 99.5])
        xs = np.linspace(min(p1, -5 * ROPE), max(p99, 5 * ROPE), 400)
        ys = gaussian_kde(delta0)(xs)
        ax_g.fill_between(xs, ys, where=(xs < -ROPE),
                          color="#E07070", alpha=1, linewidth=0)
        ax_g.fill_between(xs, ys, where=((xs >= -ROPE) & (xs <= ROPE)),
                          color="#C0C0C0", alpha=1, linewidth=0)
        ax_g.fill_between(xs, ys, where=(xs > ROPE),
                          color=METRIC_COLORS[metric], alpha=1, linewidth=0)
        ax_g.plot(xs, ys, color="#222", linewidth=0.9)
        ax_g.axvline(-ROPE, color="orange", linewidth=0.8)
        ax_g.axvline( ROPE, color="orange", linewidth=0.8)
        ax_g.axvline(0,     color="#888",   linewidth=0.5, linestyle="--")
        mean_d   = float(np.mean(delta0))
        p_better = float(np.mean(delta0 > ROPE))
        sign     = "+" if mean_d >= 0 else ""
        ax_g.text(0.98, 0.97, f"{sign}{mean_d:.3f}\n ({p_better:.0%})",
                  ha="right", va="top", fontsize=12, color="#333",
                  transform=ax_g.transAxes)
        ax_g.axvline(mean_d, color="#333", linewidth=0.8, linestyle="-",
                     alpha=1, zorder=4)
        ax_g.axhspan(0, max(ys) + 0.5, color="#FFFACD", zorder=0, linewidth=0)
        ax_g.set_yticks([])
        ax_g.tick_params(axis="x", labelsize=12)
        ax_g.spines[["top", "right", "left"]].set_visible(False)

    # Row labels
    for row, ep in enumerate(ep_order):
        axes[row, 0].set_ylabel(
            ep.replace("_", " "), fontsize=12, rotation=0, ha="right", labelpad=4)
    axes[n_eps, 0].set_ylabel(
        r"Global $\delta_0$", fontsize=12, rotation=0, ha="right", labelpad=4,
        fontweight="bold")

    for ext in ("png", "pdf"):
        out = IMGS_DIR / f"bayesian_posteriors_hierarchical{suffix}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"[saved] {out.relative_to(IMGS_DIR.parent)}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(dataset: str = "toxcast") -> None:
    suffix = f"_{dataset}" if dataset != "toxcast" else ""

    raw = pd.read_csv(DATA_DIR / f"{dataset}_comparison_results.csv")
    df  = raw[(raw["experiment"] == "Exp A") &
              (raw["model"] == "GradientBoosting")].copy()

    endpoints  = sorted(df["endpoint"].unique())
    # All PFASGroups feature sets (exclude baseline)
    pfg_sets   = sorted([fs for fs in df["feature_set"].unique()
                         if fs != BASELINE and "TxP_PFAS" not in fs])

    rows       = []
    diffs_store = {}   # (endpoint, metric, feature_set) -> np.ndarray of per-fold diffs
    for ep in endpoints:
        ep_df = df[df["endpoint"] == ep]
        for metric in METRICS:
            # Baseline scores per fold
            base_df = ep_df[ep_df["feature_set"] == BASELINE][["repeat", "fold", metric]]
            if base_df.empty:
                continue

            for fs in pfg_sets:
                cmp_df = ep_df[ep_df["feature_set"] == fs][["repeat", "fold", metric]]
                if cmp_df.empty:
                    continue

                merged = base_df.merge(cmp_df, on=["repeat", "fold"],
                                       suffixes=("_base", "_cmp"))
                if merged.empty:
                    continue

                diffs = merged[f"{metric}_cmp"].values - merged[f"{metric}_base"].values
                diffs_store[(ep, metric, fs)] = diffs
                # correlated_ttest(x, rope, runs) returns (p_left, p_rope, p_right)
                # where x = diffs = PFG - TxP, so:
                #   p_left  = P(TxP > PFG)  (diff < -rope)
                #   p_right = P(PFG > TxP)  (diff > rope)
                p_txp, p_eq, p_pfg = correlated_ttest(diffs, ROPE, runs=N_REPEATS)
                rows.append({
                    "endpoint":     ep,
                    "metric":       metric,
                    "feature_set":  fs,
                    "mean_diff":    round(float(np.mean(diffs)), 4),
                    "p_pfg_better": round(float(p_pfg), 4),
                    "p_rope":       round(float(p_eq),  4),
                    "p_txp_better": round(float(p_txp), 4),
                    "n_pairs":      len(diffs),
                })

    result_df = pd.DataFrame(rows)
    out_csv = DATA_DIR / f"bayesian_comparison{suffix}.csv"
    result_df.to_csv(out_csv, index=False)
    print(f"[saved] {out_csv.relative_to(DATA_DIR.parent)}")

    # ------------------------------------------------------------------
    # For each endpoint and metric: best PFG set (highest mean_diff)
    # ------------------------------------------------------------------
    best_rows = (
        result_df
        .loc[result_df.groupby(["endpoint", "metric"])["mean_diff"].idxmax()]
        .copy()
    )

    # Build dict of best-feature-set diffs for the posterior figure
    best_diffs = {
        (row["endpoint"], row["metric"]): diffs_store[
            (row["endpoint"], row["metric"], row["feature_set"])
        ]
        for _, row in best_rows.iterrows()
        if (row["endpoint"], row["metric"], row["feature_set"]) in diffs_store
    }

    # ------------------------------------------------------------------
    # Figure 1: P(best PFG > TxP_PFAS) per endpoint and metric
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        1, len(METRICS),
        figsize=(len(METRICS) * 3.8, max(4.0, len(endpoints) * 0.45 + 1.5)),
        constrained_layout=True,
    )

    fig.suptitle(
        "Bayesian Correlated $t$-test: P(best PFASGroups $>$ TxP\\_PFAS)\n"
        "(Experiment~A, Gradient Boosting, ROPE = $\\pm$0.01)",
        fontsize=14, fontweight="bold",
    )

    ep_order = sorted(endpoints)
    n_ep = len(ep_order)

    for col_idx, (ax, metric) in enumerate(zip(axes, METRICS)):
        sub = best_rows[best_rows["metric"] == metric].set_index("endpoint")
        probs       = [sub.loc[ep, "p_pfg_better"] if ep in sub.index else float("nan")
                       for ep in ep_order]
        mean_diffs  = [sub.loc[ep, "mean_diff"]    if ep in sub.index else float("nan")
                       for ep in ep_order]

        y = list(range(n_ep))

        # Light-grey band on even rows (y = 0, 2, 4, …)
        for yy in y:
            if yy % 2 == 0:
                ax.axhspan(yy - 0.5, yy + 0.5, color="#f0f0f0", zorder=0, linewidth=0)

        colors = [METRIC_COLORS[metric]] * n_ep
        ax.barh(y, probs, color=colors, edgecolor="white", linewidth=0.4, zorder=2)

        # Annotate mean difference
        for yy, p, d in zip(y, probs, mean_diffs):
            if not math.isnan(p) and not math.isnan(d):
                sign = "+" if d >= 0 else ""
                fontweight = 'normal'
                fontcolor = "#8C8C8C"
                if p >= 0.999:
                    significance = "***"
                    fontweight = 'bold'
                    fontcolor = "#C0A300" # gold
                elif p >= 0.99:
                    significance = "**"
                    fontweight = 'bold'
                    fontcolor = "#5E5D5D" # silver
                elif p >= 0.95:
                    significance = "*"
                    fontweight = 'bold'
                    fontcolor = "#CD8A46" # bronze
                else:
                    significance = ""
                ax.text(min(p + 0.02, 0.98), yy, f"{sign}{d:.3f}{significance}",
                        va="center", ha="left", fontsize=11, color=fontcolor, fontweight=fontweight, zorder=3)

        ax.axvline(0.5,  color="#aaa", linewidth=0.8, linestyle="--", zorder=1)
        ax.axvline(0.95, color="#888", linewidth=0.6, linestyle=":",  zorder=1)
        ax.set_xlim(0, 1.12)
        ax.set_ylim(-0.5, n_ep - 0.5)
        ax.set_yticks(y)
        ax.set_xlabel("P(PFASGroups > TxP\\_PFAS)", fontsize=12)
        ax.set_title(METRIC_LABELS[metric], fontsize=12, fontweight="bold")
        ax.grid(axis="x", linewidth=0.3, alpha=0.5, zorder=1)
        ax.spines[["top", "right"]].set_visible(False)

        if col_idx == 0:
            # Show y-tick labels only on the leftmost panel
            ax.set_yticklabels(ep_order, fontsize=12)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
            ax.spines["left"].set_visible(False)

    fig.text(
        0.5, -0.02,
        "Bar annotations show the mean difference (best PFASGroups\u2009\u2212\u2009TxP\u2011PFAS). "
        "Stars indicate posterior probability: "
        r"$* \, P\!\geq\!0.95$; "
        r"$** \, P\!\geq\!0.99$;"
        r"$*** \, P\!\geq\!0.999$.",
        ha="center", va="top", fontsize=12, style="normal", color="#000000",
        transform=fig.transFigure,
        wrap=True,
    )

    for ext in ("png", "pdf"):
        out = IMGS_DIR / f"bayesian_comparison{suffix}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"[saved] {out.relative_to(IMGS_DIR.parent)}")
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 2: Posterior distributions with ROPE
    # ------------------------------------------------------------------
    _draw_posterior_figure(best_diffs, ep_order, best_rows=best_rows, suffix=suffix)

    # ------------------------------------------------------------------
    # Figure 3: Hierarchical Bayesian posterior distributions
    # ------------------------------------------------------------------
    _draw_hierarchical_posterior_figure(diffs_store, ep_order,
                                        best_rows=best_rows, suffix=suffix)

    # ------------------------------------------------------------------
    # Print LaTeX table for best PFG per endpoint x metric
    # ------------------------------------------------------------------
    print("\n% === LaTeX table (Bayesian comparison, best PFG vs TxP_PFAS) ===")
    print(r"\begin{table}[h]")
    print(r"\centering")
    print(r"\small")
    print(r"\caption{Bayesian correlated $t$-test: posterior probability that the best "
          r"PFASGroups feature set outperforms TxP\_PFAS (Experiment~A, Gradient Boosting, "
          r"ROPE $=\pm0.01$). Values $\geq 0.95$ are in bold.}")
    print(r"\label{tab:bayesian_comparison}")
    print(r"\begin{tabular}{lcccc}")
    print(r"\toprule")
    print(r"Endpoint & P(PFG>TxP) ROC-AUC & P(PFG>TxP) AP & P(PFG>TxP) MCC & P(PFG>TxP) Bal.Acc \\")
    print(r"\midrule")
    for ep in ep_order:
        vals = []
        best_fsets = []
        for metric in METRICS:
            sub = best_rows[(best_rows["endpoint"] == ep) & (best_rows["metric"] == metric)]
            if sub.empty:
                vals.append("--")
                best_fsets.append("")
            else:
                p  = sub.iloc[0]["p_pfg_better"]
                fs = sub.iloc[0]["feature_set"]
                best_fsets.append(fs)
                s  = f"{p:.2f}"
                vals.append(r"\textbf{" + s + r"}" if p >= 0.95 else s)
        print(f"  {ep.replace('_', r'\_')} & " + " & ".join(vals) + r" \\")
    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")

    print("\nDone.")


if __name__ == "__main__":
    import argparse
    _parser = argparse.ArgumentParser(
        description="Bayesian correlated t-test comparison of PFASGroups vs TxP_PFAS."
    )
    _parser.add_argument(
        "--dataset", "-d", default="toxcast",
        help="Dataset prefix used in input CSV and output filenames (default: 'toxcast').",
    )
    _args = _parser.parse_args()
    main(_args.dataset)
