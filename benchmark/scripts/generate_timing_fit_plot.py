import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib

# Use Type 1 fonts for publication
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Data extracted from HTML timing scatter plot
atom_times = [
    (13, [13.82, 13.87, 15.03, 13.67, 12.97, 14.28, 13.65, 14.94]),
    (16, [17.78, 17.93, 18.55, 17.64, 17.83, 18.20, 17.62, 17.11]),
    (19, [22.26, 22.05, 21.59, 22.77, 22.28, 23.15, 22.88, 21.34]),
    (22, [30.93, 28.39, 27.05, 27.46, 27.75, 29.23, 27.51, 27.89]),
    (34, [54.49, 55.33, 56.19, 53.20, 52.15, 53.67, 53.49, 55.28, 58.22, 58.09, 56.69]),
    (70, [219.68, 227.60, 249.91, 218.59]),
    (73, [238.95, 230.51, 231.27, 234.15]),
    (76, [352.06, 267.30, 257.18, 235.88]),
    (79, [313.77, 284.71, 270.64, 271.29, 278.34, 280.55, 266.94]),
    (127, [374.41, 380.80, 383.19, 379.81]),
    (154, [663.66, 576.62, 626.78]),
]

# Create arrays for regression
atom_counts = []
exec_times = []

for atoms, times in atom_times:
    for time in times:
        atom_counts.append(atoms)
        exec_times.append(time)

atom_counts = np.array(atom_counts)
exec_times = np.array(exec_times)

print(f"Data points: {len(atom_counts)}")
print(f"Atom range: {atom_counts.min()}-{atom_counts.max()}")
print(f"Time range: {exec_times.min():.1f}-{exec_times.max():.1f} ms")

# Define exponential model
def exponential_model(x, a, b):
    """time = a * exp(b * x)"""
    return a * np.exp(b * x)

# Fit exponential model
popt, pcov = curve_fit(exponential_model, atom_counts, exec_times, 
                       p0=[52, 0.016], maxfev=10000)

# Calculate predictions and residuals
pred = exponential_model(atom_counts, *popt)
residuals = exec_times - pred
ss_res = np.sum(residuals**2)
ss_tot = np.sum((exec_times - np.mean(exec_times))**2)
r2 = 1 - (ss_res / ss_tot)
rmse = np.sqrt(np.mean(residuals**2))

# Calculate standard error and confidence intervals
n = len(atom_counts)
p = len(popt)  # number of parameters
dof = n - p  # degrees of freedom
mse = ss_res / dof  # mean squared error

# t-value for 95% confidence interval
t_val = stats.t.ppf(0.975, dof)

print(f"\nExponential Model: t = {popt[0]:.3f} * exp({popt[1]:.5f} * n)")
print(f"R² = {r2:.4f}, RMSE = {rmse:.1f} ms")
print(f"Standard errors: a = {np.sqrt(pcov[0,0]):.3f}, b = {np.sqrt(pcov[1,1]):.6f}")

# Create prediction grid
x_grid = np.linspace(atom_counts.min(), atom_counts.max(), 200)
y_pred = exponential_model(x_grid, *popt)

# Calculate confidence interval using error propagation
# For y = a * exp(b * x), the variance is:
# var(y) ≈ (∂y/∂a)² var(a) + (∂y/∂b)² var(b) + 2(∂y/∂a)(∂y/∂b) cov(a,b)
# ∂y/∂a = exp(b*x)
# ∂y/∂b = a*x*exp(b*x)

dy_da = np.exp(popt[1] * x_grid)
dy_db = popt[0] * x_grid * np.exp(popt[1] * x_grid)

var_y = (dy_da**2 * pcov[0,0] + 
         dy_db**2 * pcov[1,1] + 
         2 * dy_da * dy_db * pcov[0,1])

se_y = np.sqrt(var_y)
ci_lower = y_pred - t_val * se_y
ci_upper = y_pred + t_val * se_y

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot data points with transparency
ax.scatter(atom_counts, exec_times, alpha=0.5, s=50, color='#667eea', 
           label='Measured times', edgecolors='darkblue', linewidth=0.5)

# Plot fit line
ax.plot(x_grid, y_pred, 'r-', linewidth=2, 
        label=f'Exponential fit: $t = {popt[0]:.1f} \\times \\exp({popt[1]:.4f} n)$')

# Plot confidence interval
ax.fill_between(x_grid, ci_lower, ci_upper, alpha=0.2, color='red',
                label='95% confidence interval')

# Labels and formatting
ax.set_xlabel('Atom count ($n$)', fontsize=14, fontweight='bold')
ax.set_ylabel('Execution time (ms)', fontsize=14, fontweight='bold')
ax.set_title('HalogenGroups Timing Dataset: Exponential Scaling with Molecular Size', 
             fontsize=15, fontweight='bold', pad=20)
ax.legend(fontsize=12, loc='upper left', framealpha=0.95)
ax.grid(True, alpha=0.3, linestyle='--')

# Add statistics text box
textstr = f'$R^2 = {r2:.4f}$\nRMSE = {rmse:.1f} ms\n$n$ = {n} measurements\n$\\alpha = {popt[1]:.4f}$ atoms$^{{-1}}$'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.98, 0.35, textstr, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment='right', bbox=props)

# Set y-axis to start at 0 for better visualization
ax.set_ylim(bottom=0)

plt.tight_layout()

# Save figure
output_path = 'timing_exponential_fit.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_path}")

# Also save as PNG for backup
output_path_png = 'timing_exponential_fit.png'
plt.savefig(output_path_png, dpi=300, bbox_inches='tight')
print(f"PNG version saved to: {output_path_png}")

plt.show()

# Print confidence interval statistics
print(f"\n95% Confidence Interval Statistics:")
print(f"At n=20: {exponential_model(20, *popt):.1f} ms ± {t_val * np.sqrt((np.exp(popt[1]*20))**2 * pcov[0,0] + (popt[0]*20*np.exp(popt[1]*20))**2 * pcov[1,1] + 2*np.exp(popt[1]*20)*(popt[0]*20*np.exp(popt[1]*20))*pcov[0,1]):.1f} ms")
print(f"At n=50: {exponential_model(50, *popt):.1f} ms ± {t_val * np.sqrt((np.exp(popt[1]*50))**2 * pcov[0,0] + (popt[0]*50*np.exp(popt[1]*50))**2 * pcov[1,1] + 2*np.exp(popt[1]*50)*(popt[0]*50*np.exp(popt[1]*50))*pcov[0,1]):.1f} ms")
print(f"At n=100: {exponential_model(100, *popt):.1f} ms ± {t_val * np.sqrt((np.exp(popt[1]*100))**2 * pcov[0,0] + (popt[0]*100*np.exp(popt[1]*100))**2 * pcov[1,1] + 2*np.exp(popt[1]*100)*(popt[0]*100*np.exp(popt[1]*100))*pcov[0,1]):.1f} ms")
