#!/usr/bin/env python3
"""
Train and evaluate ML models on the ToxCast / PFASGroups dataset.

For each endpoint:
  - Drop chemicals not tested in that assay (NaN label)
  - 5-fold stratified cross-validation
  - Three classifiers: RandomForest, GradientBoosting, LogisticRegression
  - Metrics: ROC-AUC, Average Precision (AUPRC), MCC, Balanced-Accuracy

Outputs (all written to ../data/):
  toxcast_cv_results.csv      – per-endpoint × per-model metrics
  toxcast_auc_heatmap.png     – heatmap of ROC-AUC
  toxcast_auprc_heatmap.png   – heatmap of Average Precision
  models/                     – final models trained on full data (joblib)
  toxcast_feature_importance.csv – mean RF feature importances per endpoint

Usage
-----
    conda activate chem
    cd benchmark
    python scripts/train_toxcast_models.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
MODEL_DIR  = DATA_DIR / "models"
MODEL_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DATASET_PATH   = DATA_DIR / "toxcast_dataset.parquet"
MIN_POSITIVES  = 20      # skip endpoint if fewer than this many positives
N_FOLDS        = 5
RANDOM_STATE   = 42

META_COLS  = ["chid", "casn", "chnm", "dsstox_substance_id", "smiles"]
LABEL_COLS = [
    "AR_antagonist", "AR_agonist",
    "ERa_antagonist", "ERa_agonist",
    "AhR_agonist", "Aromatase_antagonist",
    "TR_antagonist", "DT40_genotoxicity",
    "MMP_ratio",
    "CYP2D6_antagonist", "CYP2C19_antagonist", "CYP2C9_antagonist",
    "CYP3A4_antagonist", "p53_ratio", "Caspase3_HEPG2",
]

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

def make_models() -> dict[str, object]:
    return {
        "RandomForest": RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced",
            max_features="sqrt",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        ),
        "GradientBoosting": GradientBoostingClassifier(
            n_estimators=150,
            learning_rate=0.05,
            max_depth=3,
            subsample=0.8,
            random_state=RANDOM_STATE,
        ),
        "LogisticRegression": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                class_weight="balanced",
                max_iter=2000,
                C=0.1,
                solver="lbfgs",
                random_state=RANDOM_STATE,
            )),
        ]),
    }


# ---------------------------------------------------------------------------
# Training & evaluation
# ---------------------------------------------------------------------------

def evaluate_endpoint(
    X: np.ndarray,
    y: np.ndarray,
    endpoint: str,
    fp_col_names: list[str],
) -> tuple[dict, dict, dict]:
    """
    Run 5-fold stratified CV for all models on one endpoint.

    Returns
    -------
    cv_metrics      : {model_name: {metric: value}}
    trained_models  : {model_name: fitted model}
    importances     : {model_name: array of length n_features}  (RF only, else None)
    """
    skf     = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    models  = make_models()
    cv_metrics     : dict[str, dict] = {}
    trained_models : dict[str, object] = {}
    importances    : dict[str, np.ndarray | None] = {}

    n_pos = int(y.sum())
    n_neg = int((y == 0).sum())

    for name, model in models.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # OOB probability predictions from CV
            y_prob = cross_val_predict(
                model, X, y,
                cv=skf,
                method="predict_proba",
                n_jobs=1,
            )[:, 1]
            y_pred = (y_prob >= 0.5).astype(int)

        auc = roc_auc_score(y, y_prob)
        ap  = average_precision_score(y, y_prob)
        mcc = matthews_corrcoef(y, y_pred)
        bal = balanced_accuracy_score(y, y_pred)

        cv_metrics[name] = {
            "endpoint":   endpoint,
            "model":      name,
            "n_tested":   len(y),
            "n_pos":      n_pos,
            "n_neg":      n_neg,
            "pos_rate":   round(n_pos / len(y), 4),
            "roc_auc":    round(auc, 4),
            "avg_prec":   round(ap, 4),
            "mcc":        round(mcc, 4),
            "bal_acc":    round(bal, 4),
        }

        # Train on full data for final model
        final = make_models()[name]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            final.fit(X, y)
        trained_models[name] = final

        # Feature importance (RandomForest only)
        if name == "RandomForest":
            importances[name] = final.feature_importances_
        else:
            importances[name] = None

    return cv_metrics, trained_models, importances


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_heatmap(
    pivot: pd.DataFrame,
    title: str,
    path: Path,
    vmin: float = 0.5,
    vmax: float = 1.0,
    cmap: str = "YlOrRd",
    annot_fmt: str = ".2f",
) -> None:
    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns) * 1.1), max(4, len(pivot) * 0.6)))
    mask = pivot.isna()
    sns.heatmap(
        pivot,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        annot=True,
        fmt=annot_fmt,
        linewidths=0.5,
        mask=mask,
        cbar_kws={"label": title.split("–")[-1].strip()},
    )
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"[plot] saved → {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Loading dataset from {DATASET_PATH} …")
    df = pd.read_parquet(DATASET_PATH)

    # Fingerprint columns = everything except meta and label cols
    fp_cols = [c for c in df.columns if c not in META_COLS and c not in LABEL_COLS]
    print(f"  {len(df):,} chemicals  |  {len(fp_cols)} fingerprint features  |  {len(LABEL_COLS)} endpoints\n")

    X_all = df[fp_cols].values.astype(np.float32)

    all_metrics   : list[dict] = []
    all_models    : dict[str, dict[str, object]] = {}   # endpoint → {model → fitted}
    all_importance: dict[str, np.ndarray] = {}          # endpoint → RF importances

    for ep in LABEL_COLS:
        if ep not in df.columns:
            print(f"  [skip] {ep} – not in dataset")
            continue

        mask = df[ep].notna()
        y    = df.loc[mask, ep].values.astype(int)
        X    = X_all[mask.values]

        n_pos = int(y.sum())
        if n_pos < MIN_POSITIVES:
            print(f"  [skip] {ep:<30s}  n_pos={n_pos} < {MIN_POSITIVES}")
            continue

        print(f"  [train] {ep:<30s}  n={len(y):5d}  pos={n_pos:4d}  ({100*n_pos/len(y):.1f}%)")

        cv_metrics, trained_models, importances = evaluate_endpoint(X, y, ep, fp_cols)
        all_metrics.extend(cv_metrics.values())
        all_models[ep]     = trained_models
        if importances["RandomForest"] is not None:
            all_importance[ep] = importances["RandomForest"]

    if not all_metrics:
        print("\nNo endpoints had enough positives — nothing to report.")
        return

    # ------------------------------------------------------------------
    # Save CV results
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(all_metrics)
    results_path = DATA_DIR / "toxcast_cv_results.csv"
    results_df.to_csv(results_path, index=False)
    print(f"\n[saved] CV results → {results_path}")

    # ------------------------------------------------------------------
    # Print summary table
    # ------------------------------------------------------------------
    pivot_auc  = results_df.pivot(index="model", columns="endpoint", values="roc_auc")
    pivot_ap   = results_df.pivot(index="model", columns="endpoint", values="avg_prec")
    pivot_mcc  = results_df.pivot(index="model", columns="endpoint", values="mcc")

    print("\n=== ROC-AUC (5-fold CV) ===")
    print(pivot_auc.to_string(float_format="%.3f"))

    print("\n=== Average Precision (5-fold CV) ===")
    print(pivot_ap.to_string(float_format="%.3f"))

    print("\n=== MCC (5-fold CV) ===")
    print(pivot_mcc.to_string(float_format="%.3f"))

    # ------------------------------------------------------------------
    # Heatmaps
    # ------------------------------------------------------------------
    plot_heatmap(pivot_auc, "ToxCast / PFASGroups – ROC-AUC (5-fold CV)",
                 DATA_DIR / "toxcast_auc_heatmap.png", vmin=0.5, vmax=1.0)
    plot_heatmap(pivot_ap,  "ToxCast / PFASGroups – Average Precision (5-fold CV)",
                 DATA_DIR / "toxcast_auprc_heatmap.png", vmin=0.0, vmax=1.0)
    plot_heatmap(pivot_mcc, "ToxCast / PFASGroups – MCC (5-fold CV)",
                 DATA_DIR / "toxcast_mcc_heatmap.png", vmin=-0.2, vmax=0.8, cmap="RdYlGn")

    # ------------------------------------------------------------------
    # Feature importance (aggregate RF importances across endpoints)
    # ------------------------------------------------------------------
    if all_importance:
        imp_df = pd.DataFrame(all_importance, index=fp_cols)
        imp_df["mean"] = imp_df.mean(axis=1)
        imp_df = imp_df.sort_values("mean", ascending=False)
        imp_path = DATA_DIR / "toxcast_feature_importance.csv"
        imp_df.to_csv(imp_path)
        print(f"[saved] feature importances → {imp_path}")

        # Bar plot of top 20 features (mean across endpoints)
        top20 = imp_df["mean"].head(20)
        fig, ax = plt.subplots(figsize=(10, 5))
        top20.plot.bar(ax=ax, color="steelblue")
        ax.set_title("Top 20 PFASGroups features (mean RF importance across all endpoints)")
        ax.set_ylabel("Mean importance")
        ax.tick_params(axis="x", labelsize=8)
        plt.tight_layout()
        fig.savefig(DATA_DIR / "toxcast_top_features.png", dpi=150)
        plt.close(fig)
        print(f"[plot] saved → {DATA_DIR / 'toxcast_top_features.png'}")

    # ------------------------------------------------------------------
    # Save trained models
    # ------------------------------------------------------------------
    for ep, model_dict in all_models.items():
        for model_name, fitted_model in model_dict.items():
            safe_ep   = ep.replace(" ", "_").replace("/", "_")
            save_path = MODEL_DIR / f"{safe_ep}__{model_name}.joblib"
            joblib.dump(fitted_model, save_path)
    print(f"[saved] {sum(len(v) for v in all_models.values())} model files → {MODEL_DIR}/")

    print("\nDone.")
    print("  Outputs:")
    print(f"    {results_path}")
    print(f"    {DATA_DIR / 'toxcast_auc_heatmap.png'}")
    print(f"    {DATA_DIR / 'toxcast_auprc_heatmap.png'}")
    print(f"    {DATA_DIR / 'toxcast_mcc_heatmap.png'}")
    print(f"    {DATA_DIR / 'toxcast_feature_importance.csv'}")
    print(f"    {DATA_DIR / 'toxcast_top_features.png'}")


if __name__ == "__main__":
    main()
