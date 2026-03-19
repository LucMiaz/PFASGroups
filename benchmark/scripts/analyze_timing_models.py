"""
Analyze timing benchmark data and compare multiple complexity models:
  n          linear:         t = a*n + c
  log(n)     logarithmic:    t = a*log(n) + c
  n*log(n)   linearithmic:   t = a*n*log(n) + c
  n^2        quadratic:      t = a*n^2 + c
  n^3        cubic:          t = a*n^3 + c
  exp(n)     exponential:    t = a*exp(b*n)

All models are fitted to the measured time vs. atom-count data.
Results are ranked by R^2 and written to JSON + LaTeX table.
"""

import json
import re
import numpy as np
import pandas as pd
try:
    from scipy.optimize import curve_fit
    from scipy.stats import pearsonr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available, will use numpy polyfit fallback")

try:
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['font.family'] = 'Ubuntu'
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available, skipping plots")

import glob
from pathlib import Path


# ── palette ──────────────────────────────────────────────────────────────────
def _load_palette():
    """Load hex colours from color_scheme.yaml (stdlib only, no pyyaml needed)."""
    _defaults = ["#E15D0B", "#306DBA", "#9D206C", "#51127C",
                 "#2CA02C", "#D62728", "#9467BD"]
    try:
        _p = Path(__file__).parent.parent.parent / "PFASGroups" / "data" / "color_scheme.yaml"
        _colors = re.findall(r'"(#[0-9A-Fa-f]{6})"', _p.read_text())
        if len(_colors) >= 4:
            return _colors[:4] + _defaults[4:]
    except Exception:
        pass
    return _defaults


_PALETTE = _load_palette()
_C_DATA = _PALETTE[0]   # orange  – raw data points

# one colour per model (cycling if needed)
_MODEL_COLORS = _PALETTE[1:]


# ── model definitions ─────────────────────────────────────────────────────────
# Each entry: (key, label, latex_label, func, p0, param_names)
#   func(n, *params) -> predicted time (seconds)
#   p0               -> initial parameter guess for curve_fit
#   param_names      -> list of parameter names for reporting

def _linear(n, a, c):
    return a * n + c

def _log(n, a, c):
    return a * np.log(n) + c

def _nlogn(n, a, c):
    return a * n * np.log(n) + c

def _quadratic(n, a, c):
    return a * n**2 + c

def _cubic(n, a, c):
    return a * n**3 + c

def _exponential(n, a, b):
    return a * np.exp(b * n)

MODELS = [
    ("linear",      "O(n)",         r"$t = a n + c$",            _linear,      [1e-4,  1e-4], ["a", "c"]),
    ("log",         "O(log n)",     r"$t = a \log n + c$",       _log,         [1e-4,  1e-4], ["a", "c"]),
    ("nlogn",       "O(n log n)",   r"$t = a n \log n + c$",     _nlogn,       [1e-6,  1e-4], ["a", "c"]),
    ("quadratic",   "O(n^2)",       r"$t = a n^2 + c$",          _quadratic,   [1e-8,  1e-4], ["a", "c"]),
    ("cubic",       "O(n^3)",       r"$t = a n^3 + c$",          _cubic,       [1e-11, 1e-4], ["a", "c"]),
    ("exponential", "O(exp(b*n))",  r"$t = a e^{b n}$",          _exponential, [1e-3,  0.01], ["a", "b"]),
]


# ── data loading ──────────────────────────────────────────────────────────────
timing_files = sorted(
    glob.glob('data/pfas_timing_benchmark_*.json'),
    key=lambda x: Path(x).stat().st_mtime, reverse=True
)
if not timing_files:
    print("Error: No timing benchmark files found in data/")
    exit(1)

timing_file = timing_files[0]
print(f"Loading latest timing benchmark: {timing_file}")

with open(timing_file, 'r') as f:
    timing_data = json.load(f)

df = pd.DataFrame([{
    'molecule_id':        d['molecule_id'],
    'chain_length':       d['chain_length'],
    'num_atoms':          d['num_atoms'],
    'PFASGroups_time_avg': d['PFASGroups_time_avg'],
    'PFASGroups_time_std': d['PFASGroups_time_std'],
} for d in timing_data])

x_data = df['num_atoms'].values
y_data = df['PFASGroups_time_avg'].values


# ── fitting ───────────────────────────────────────────────────────────────────
def fit_model(key, func, p0, param_names):
    """Fit a single model; return a result dict."""
    try:
        if SCIPY_AVAILABLE:
            popt, pcov = curve_fit(func, x_data, y_data, p0=p0, maxfev=50000)
            perr = np.sqrt(np.diag(pcov))
        else:
            # Fallback: numpy least-squares via 1-D polyfit on the feature vector
            # Build feature matrix for a 2-param affine model on the relevant basis
            raise RuntimeError("scipy required for non-linear fitting")
    except Exception as exc:
        print(f"  [WARN] {key}: fitting failed ({exc}), skipping")
        return None

    y_pred = func(x_data, *popt)
    residuals = y_data - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    # AIC / BIC (Gaussian likelihood)
    n = len(y_data)
    k = len(popt)
    sigma2 = ss_res / n
    log_lik = -n / 2 * np.log(2 * np.pi * sigma2) - ss_res / (2 * sigma2)
    aic = 2 * k - 2 * log_lik
    bic = k * np.log(n) - 2 * log_lik

    params = {name: {'value': float(v), 'error': float(e)}
              for name, v, e in zip(param_names, popt, perr)}

    return {
        'key':        key,
        'params':     params,
        'r_squared':  float(r_squared),
        'rmse_ms':    float(np.sqrt(ss_res / n) * 1000),
        'aic':        float(aic),
        'bic':        float(bic),
        'y_pred':     y_pred,
        'residuals':  residuals,
    }


print("=" * 80)
print("TIMING ANALYSIS – MODEL COMPARISON")
print("=" * 80)
print(f"\nDataset: {len(df)} molecules, {x_data.min()}–{x_data.max()} atoms")
print(f"Mean time: {y_data.mean()*1000:.2f} ms  |  Median: {np.median(y_data)*1000:.2f} ms\n")
print(f"{'Model':<18} {'R²':>8} {'RMSE(ms)':>10} {'AIC':>12} {'BIC':>12}  Parameters")
print("-" * 80)

fit_results = []
for key, label, latex, func, p0, pnames in MODELS:
    res = fit_model(key, func, p0, pnames)
    if res is None:
        continue
    res['label'] = label
    res['latex'] = latex
    fit_results.append(res)

    param_str = "  ".join(
        f"{n}={v['value']:.3e}±{v['error']:.1e}"
        for n, v in res['params'].items()
    )
    print(f"{label:<18} {res['r_squared']:>8.4f} {res['rmse_ms']:>10.3f} "
          f"{res['aic']:>12.1f} {res['bic']:>12.1f}  {param_str}")

# Sort by R² descending
fit_results.sort(key=lambda r: r['r_squared'], reverse=True)
best = fit_results[0]
print("-" * 80)
print(f"\n>>> Best model by R²: {best['label']}  (R²={best['r_squared']:.4f})")
print(f">>> Best model by AIC: "
      f"{min(fit_results, key=lambda r: r['aic'])['label']}")
print(f">>> Best model by BIC: "
      f"{min(fit_results, key=lambda r: r['bic'])['label']}")


# ── plots ─────────────────────────────────────────────────────────────────────
if MATPLOTLIB_AVAILABLE:
    n_models = len(fit_results)
    x_fit = np.linspace(x_data.min(), x_data.max(), 300)

    # 1) Comparison: data + all fit lines
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.scatter(x_data, y_data * 1000, s=12, alpha=0.35, color=_C_DATA,
               zorder=2, label='Measured data')
    for i, res in enumerate(fit_results):
        func = next(f for k, _, _, f, _, _ in MODELS if k == res['key'])
        pvals = [v['value'] for v in res['params'].values()]
        y_fit = func(x_fit, *pvals) * 1000
        color = _MODEL_COLORS[i % len(_MODEL_COLORS)]
        ax.plot(x_fit, y_fit, linewidth=2, color=color,
                label=f"{res['label']}  R²={res['r_squared']:.3f}")
    ax.set_xlabel('Number of atoms', fontsize=12)
    ax.set_ylabel('Execution time (ms)', fontsize=12)
    ax.set_title('PFASGroups timing: complexity model comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('reports/timing_analysis_all_metrics.png', dpi=300, bbox_inches='tight')
    fig.savefig('reports/timing_analysis_all_metrics.pdf', dpi=300, bbox_inches='tight')
    print(f"\nComparison plot saved: reports/timing_analysis_all_metrics.png / .pdf")

    # 2) Residuals grid (one subplot per model)
    ncols = 3
    nrows = (n_models + ncols - 1) // ncols
    fig2, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                               sharex=True)
    axes_flat = list(axes.flat) if hasattr(axes, 'flat') else [axes]
    for i, res in enumerate(fit_results):
        ax = axes_flat[i]
        color = _MODEL_COLORS[i % len(_MODEL_COLORS)]
        ax.scatter(x_data, res['residuals'] * 1000, s=8, alpha=0.35, color=color)
        ax.axhline(0, color='#666666', linewidth=1, linestyle='--')
        ax.set_title(f"{res['label']}  (R²={res['r_squared']:.3f})", fontsize=10)
        ax.set_xlabel('Atoms', fontsize=9)
        ax.set_ylabel('Residual (ms)', fontsize=9)
        ax.grid(True, alpha=0.3)
    # hide unused subplots
    for j in range(n_models, nrows * ncols):
        axes_flat[j].set_visible(False)
    fig2.suptitle('Residuals per complexity model', fontsize=13, fontweight='bold')
    fig2.tight_layout()
    fig2.savefig('reports/timing_analysis_residuals_grid.png', dpi=300, bbox_inches='tight')
    fig2.savefig('reports/timing_analysis_residuals_grid.pdf', dpi=300, bbox_inches='tight')
    print(f"Residuals grid saved:   reports/timing_analysis_residuals_grid.png / .pdf")

    plt.close('all')
else:
    print("\nSkipping plot generation (matplotlib not available)")


# ── JSON output ───────────────────────────────────────────────────────────────
output = {
    'dataset_statistics': {
        'n_molecules':       int(len(df)),
        'atom_count_range':  [int(x_data.min()), int(x_data.max())],
        'chain_length_range':[int(df['chain_length'].min()), int(df['chain_length'].max())],
        'mean_time_ms':      float(y_data.mean() * 1000),
        'median_time_ms':    float(np.median(y_data) * 1000),
        'std_time_ms':       float(y_data.std() * 1000),
    },
    'models': []
}
for res in fit_results:
    output['models'].append({
        'key':       res['key'],
        'label':     res['label'],
        'formula':   res['latex'],
        'r_squared': res['r_squared'],
        'rmse_ms':   res['rmse_ms'],
        'aic':       res['aic'],
        'bic':       res['bic'],
        'parameters': res['params'],
    })

with open('reports/timing_model_all_metrics.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved: reports/timing_model_all_metrics.json")


# ── LaTeX table ───────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("LaTeX table (sorted by R²):")
print("=" * 80)
print(r"\begin{table}[ht]")
print(r"\centering")
print(r"\begin{tabular}{llrrrr}")
print(r"\toprule")
print(r"Model & Formula & $R^2$ & RMSE (ms) & AIC & BIC \\")
print(r"\midrule")
for res in fit_results:
    print(f"{res['label']} & {res['latex']} & "
          f"{res['r_squared']:.4f} & {res['rmse_ms']:.3f} & "
          f"{res['aic']:.1f} & {res['bic']:.1f} \\\\")
print(r"\bottomrule")
print(r"\end{tabular}")
print(r"\caption{Complexity model comparison for PFASGroups timing vs.\ number of atoms.}")
print(r"\end{table}")
print("=" * 80)
