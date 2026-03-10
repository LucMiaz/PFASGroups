"""
Compare timing benchmark profiles and generate exponential fit plots + LaTeX summary.

Profiles expected:
- full
- no_resistance
- no_metrics
"""

import json
import glob
from pathlib import Path

import numpy as np

try:
    from scipy.optimize import curve_fit
    from scipy.stats import pearsonr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available, using log-linear fit")

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available, skipping plots")

PROFILE_ORDER = ["full", "no_resistance", "no_metrics"]
PROFILE_LABELS = {
    "full": "All metrics",
    "no_resistance": "No resistance",
    "no_metrics": "No metrics"
}
PROFILE_COLORS = {
    "full": "#1f77b4",
    "no_resistance": "#ff7f0e",
    "no_metrics": "#2ca02c"
}


def exponential_model(n, a, b):
    return a * np.exp(b * n)


def load_latest_profile_file(profile):
    pattern = f"data/pfas_timing_benchmark_{profile}_*.json"
    files = sorted(glob.glob(pattern), key=lambda x: Path(x).stat().st_mtime, reverse=True)
    if not files:
        return None
    return files[0]


def fit_profile(profile, timing_data):
    x_data = np.array([d["num_atoms"] for d in timing_data])
    y_data = np.array([d["PFASGroups_time_avg"] for d in timing_data])

    if len(x_data) < 2:
        raise ValueError(f"Not enough data points for profile {profile}")

    if SCIPY_AVAILABLE:
        p0 = [0.001, 0.01]
        popt, pcov = curve_fit(exponential_model, x_data, y_data, p0=p0, maxfev=10000)
        a_fit, b_fit = popt
        y_pred = exponential_model(x_data, a_fit, b_fit)
        residuals = y_data - y_pred
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        corr_coef, p_value = pearsonr(x_data, y_data)
        perr = np.sqrt(np.diag(pcov))
        a_err, b_err = perr
    else:
        log_y = np.log(y_data)
        coeffs = np.polyfit(x_data, log_y, 1)
        b_fit = coeffs[0]
        a_fit = np.exp(coeffs[1])
        a_err = 0.0
        b_err = 0.0
        y_pred = exponential_model(x_data, a_fit, b_fit)
        residuals = y_data - y_pred
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
        corr_coef = float(np.corrcoef(x_data, y_data)[0, 1]) if len(x_data) > 1 else 0.0
        p_value = 0.0

    return {
        "profile": profile,
        "n_molecules": len(x_data),
        "atom_min": int(np.min(x_data)),
        "atom_max": int(np.max(x_data)),
        "mean_time_ms": float(np.mean(y_data) * 1000.0),
        "median_time_ms": float(np.median(y_data) * 1000.0),
        "a": float(a_fit),
        "b": float(b_fit),
        "a_err": float(a_err),
        "b_err": float(b_err),
        "r_squared": float(r_squared),
        "pearson_r": float(corr_coef),
        "p_value": float(p_value),
        "x_data": x_data,
        "y_data": y_data,
        "y_pred": y_pred,
        "residuals": residuals
    }


def plot_comparison(results):
    if not MATPLOTLIB_AVAILABLE:
        return None, None

    # Plot comparison with fit lines
    fig, ax = plt.subplots(figsize=(10, 6))
    for res in results:
        profile = res["profile"]
        color = PROFILE_COLORS.get(profile, "#333333")
        label = PROFILE_LABELS.get(profile, profile)
        ax.scatter(res["x_data"], res["y_data"] * 1000.0, s=10, alpha=0.4, color=color, label=label)
        x_fit = np.linspace(res["atom_min"], res["atom_max"], 200)
        y_fit = exponential_model(x_fit, res["a"], res["b"]) * 1000.0
        ax.plot(x_fit, y_fit, color=color, linewidth=2)

    ax.set_xlabel("Number of atoms")
    ax.set_ylabel("PFASGroups time (ms)")
    ax.set_title("Timing comparison with exponential fits")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)

    comparison_path = "reports/timing_profiles_comparison.png"
    fig.tight_layout()
    fig.savefig(comparison_path, dpi=300, bbox_inches="tight")

    # Residuals plot
    fig_res, ax_res = plt.subplots(figsize=(10, 6))
    for res in results:
        profile = res["profile"]
        color = PROFILE_COLORS.get(profile, "#333333")
        label = PROFILE_LABELS.get(profile, profile)
        ax_res.scatter(res["x_data"], res["residuals"] * 1000.0, s=10, alpha=0.4, color=color, label=label)

    ax_res.axhline(y=0, color="#444444", linestyle="--", linewidth=1)
    ax_res.set_xlabel("Number of atoms")
    ax_res.set_ylabel("Residuals (ms)")
    ax_res.set_title("Residuals for exponential fits")
    ax_res.grid(True, alpha=0.3)
    ax_res.legend(fontsize=10)

    residuals_path = "reports/timing_profiles_residuals.png"
    fig_res.tight_layout()
    fig_res.savefig(residuals_path, dpi=300, bbox_inches="tight")

    return comparison_path, residuals_path


def write_latex_summary(results, comparison_path=None, residuals_path=None):
    lines = []
    lines.append("%% Timing profile comparison summary")
    lines.append("\\section*{Timing Profile Comparison}")
    lines.append("Exponential model: $t = a \\cdot e^{b n}$, where $t$ is time (s) and $n$ is atom count.")
    lines.append("\\begin{table}[ht]")
    lines.append("\\centering")
    lines.append("\\begin{tabular}{lrrrrr}")
    lines.append("\\toprule")
    lines.append("Profile & $a$ (s) & $b$ (atoms$^{-1}$) & $R^2$ & Mean (ms) & Median (ms) \\")
    lines.append("\\midrule")

    for res in results:
        label = PROFILE_LABELS.get(res["profile"], res["profile"])
        lines.append(
            f"{label} & {res['a']:.2e} & {res['b']:.4f} & {res['r_squared']:.4f} & "
            f"{res['mean_time_ms']:.1f} & {res['median_time_ms']:.1f} \\")

    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("\\end{table}")

    if comparison_path:
        lines.append("\\begin{figure}[ht]")
        lines.append("\\centering")
        lines.append("\\includegraphics[width=0.95\\linewidth]{" + comparison_path + "}")
        lines.append("\\caption{Timing comparison across profiles with exponential fits.}")
        lines.append("\\end{figure}")

    if residuals_path:
        lines.append("\\begin{figure}[ht]")
        lines.append("\\centering")
        lines.append("\\includegraphics[width=0.95\\linewidth]{" + residuals_path + "}")
        lines.append("\\caption{Residuals for exponential fits across profiles.}")
        lines.append("\\end{figure}")

    latex_path = "reports/timing_profiles_summary.tex"
    with open(latex_path, "w") as f:
        f.write("\n".join(lines))

    return latex_path


def main():
    profile_results = []
    missing_profiles = []

    for profile in PROFILE_ORDER:
        timing_file = load_latest_profile_file(profile)
        if not timing_file:
            missing_profiles.append(profile)
            continue
        print(f"Loading {profile} timing data: {timing_file}")
        with open(timing_file, "r") as f:
            timing_data = json.load(f)
        profile_results.append(fit_profile(profile, timing_data))

    if missing_profiles:
        missing_str = ", ".join(missing_profiles)
        print(f"Error: Missing timing files for profiles: {missing_str}")
        print("Run timing benchmarks for these profiles first.")
        return 1

    comparison_path, residuals_path = plot_comparison(profile_results)
    latex_path = write_latex_summary(profile_results, comparison_path, residuals_path)

    summary = {}
    for res in profile_results:
        summary[res["profile"]] = {
            "profile": res["profile"],
            "n_molecules": res["n_molecules"],
            "atom_min": res["atom_min"],
            "atom_max": res["atom_max"],
            "mean_time_ms": res["mean_time_ms"],
            "median_time_ms": res["median_time_ms"],
            "a": res["a"],
            "b": res["b"],
            "a_err": res["a_err"],
            "b_err": res["b_err"],
            "r_squared": res["r_squared"],
            "pearson_r": res["pearson_r"],
            "p_value": res["p_value"]
        }
    with open("reports/timing_profiles_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\nTiming profile comparison complete:")
    if comparison_path:
        print(f"  - Plot: {comparison_path}")
    if residuals_path:
        print(f"  - Residuals plot: {residuals_path}")
    print(f"  - LaTeX summary: {latex_path}")
    print("  - JSON summary: reports/timing_profiles_summary.json")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
