"""
gui/utils/modelling.py
──────────────────────
ML benchmark: RepeatedStratifiedKFold + HistGradientBoostingClassifier
with ROC-AUC / MCC / Balanced Accuracy / Average Precision metrics.
Bayesian correlated t-test between all fingerprint set pairs.

Returns a dict ready for display in ModellingTab.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd


def run_benchmark(
    feature_sets: Dict[str, np.ndarray],
    y: np.ndarray,
    cv_splits: int = 3,
    cv_repeats: int = 5,
    progress_cb: Optional[Callable[[int], None]] = None,
) -> Dict:
    """Run HistGBM benchmark across multiple feature sets.

    Parameters
    ----------
    feature_sets : mapping of {name: (n_mols, n_features) array}
    y : binary (0/1) target vector, length n_mols
    cv_splits / cv_repeats : cross-validation settings
    progress_cb : int callback (0–100)

    Returns
    -------
    dict with keys:
      "scores" : DataFrame (rows = feature sets, cols = metrics, values = mean±std)
      "raw"    : DataFrame (long form, all fold scores)
      "bayes"  : list of dicts with pairwise Bayesian t-test results
    """
    from sklearn.model_selection import RepeatedStratifiedKFold, cross_validate
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import make_scorer, matthews_corrcoef, balanced_accuracy_score

    y = np.asarray(y, dtype=int)
    n_sets = len(feature_sets)
    total_steps = n_sets + n_sets * (n_sets - 1) // 2 + 1
    step = 0

    cv = RepeatedStratifiedKFold(
        n_splits=cv_splits, n_repeats=cv_repeats, random_state=42
    )

    scoring = {
        "roc_auc": "roc_auc",
        "mcc": make_scorer(matthews_corrcoef),
        "bal_acc": make_scorer(balanced_accuracy_score),
        "avg_prec": "average_precision",
    }

    all_raw: List[dict] = []
    set_scores: Dict[str, Dict[str, float]] = {}

    for name, X in feature_sets.items():
        clfs = HistGradientBoostingClassifier(random_state=42)
        cv_result = cross_validate(clfs, X, y, cv=cv, scoring=scoring,
                                   return_train_score=False, n_jobs=1)
        for fold_i in range(len(cv_result["test_roc_auc"])):
            all_raw.append({
                "set": name,
                "roc_auc": cv_result["test_roc_auc"][fold_i],
                "mcc": cv_result["test_mcc"][fold_i],
                "bal_acc": cv_result["test_bal_acc"][fold_i],
                "avg_prec": cv_result["test_avg_prec"][fold_i],
            })
        set_scores[name] = {
            metric: float(np.mean(cv_result[f"test_{metric}"]))
            for metric in ("roc_auc", "mcc", "bal_acc", "avg_prec")
        }
        step += 1
        if progress_cb:
            progress_cb(int(step / total_steps * 90))

    raw_df = pd.DataFrame(all_raw)

    # ── Summary table ─────────────────────────────────────────────────────
    summary_rows = []
    for name, scores in set_scores.items():
        summary_rows.append({"Fingerprint set": name, **scores})
    summary_df = pd.DataFrame(summary_rows).set_index("Fingerprint set")

    # ── Bayesian correlated t-test (Benavoli et al. 2017) ────────────────
    names = list(feature_sets)
    bayes_results = []
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            try:
                bt = _bayesian_correlated_t_test(
                    raw_df[raw_df["set"] == a]["roc_auc"].values,
                    raw_df[raw_df["set"] == b]["roc_auc"].values,
                    cv_splits=cv_splits,
                    ROPE=0.01,
                )
                bayes_results.append({"A": a, "B": b, **bt})
            except Exception:
                bayes_results.append({"A": a, "B": b, "error": "failed"})
            step += 1
            if progress_cb:
                progress_cb(int(step / total_steps * 90))

    if progress_cb:
        progress_cb(100)

    return {
        "scores": summary_df,
        "raw": raw_df,
        "bayes": bayes_results,
    }


# ── Bayesian correlated t-test ─────────────────────────────────────────────

def _bayesian_correlated_t_test(
    scores_a: np.ndarray,
    scores_b: np.ndarray,
    cv_splits: int = 3,
    ROPE: float = 0.01,
) -> dict:
    """Bayesian correlated t-test (Benavoli et al. 2017).

    Returns probabilities: P(A > B), P(ROPE), P(B > A).
    """
    from scipy import stats

    n = len(scores_a)
    diff = scores_a - scores_b
    rho = 1.0 / cv_splits  # block correlation

    mean_diff = float(np.mean(diff))
    std_diff = float(np.std(diff, ddof=1))

    if std_diff == 0.0:
        p_a = 1.0 if mean_diff > ROPE else 0.0
        p_rope = 1.0 if abs(mean_diff) <= ROPE else 0.0
        p_b = 1.0 if mean_diff < -ROPE else 0.0
        return {"p_A_wins": p_a, "p_rope": p_rope, "p_B_wins": p_b,
                "mean_diff": mean_diff, "std_diff": std_diff}

    # Correlated variance correction (Nadeau & Bengio 2003)
    se = std_diff * np.sqrt(1.0 / n + rho / (1.0 - rho))
    df = n - 1
    t_plus_rope = (mean_diff + ROPE) / se
    t_minus_rope = (mean_diff - ROPE) / se

    p_b_wins = float(stats.t.cdf(-abs(mean_diff) / se * np.sign(mean_diff)
                                 if mean_diff <= 0 else
                                 stats.t.cdf(t_minus_rope, df=df), df=df)) \
               if mean_diff < -ROPE else float(stats.t.cdf(-abs(t_minus_rope), df=df))

    # Simpler: integrate t-distribution over 3 regions
    cdf_lo = float(stats.t.cdf((-ROPE - mean_diff) / se, df=df))
    cdf_hi = float(stats.t.cdf((+ROPE - mean_diff) / se, df=df))

    p_b = cdf_lo
    p_rope = max(0.0, cdf_hi - cdf_lo)
    p_a = max(0.0, 1.0 - cdf_hi)

    return {
        "p_A_wins": round(p_a, 4),
        "p_rope": round(p_rope, 4),
        "p_B_wins": round(p_b, 4),
        "mean_diff": round(mean_diff, 5),
        "std_diff": round(std_diff, 5),
    }
