"""
Generate exponential fit plot with 95% confidence intervals
Using complete dataset (200 measurements)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats

# Configure matplotlib for publication-quality output
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts
plt.rcParams['ps.fonttype'] = 42   # TrueType fonts

# Load the complete dataset
data = np.load('timing_full_dataset.npz')
atom_counts = data['atom_counts']
exec_times = data['exec_times']

print(f"Data points: {len(atom_counts)}")
print(f"Atom range: {atom_counts.min()}-{atom_counts.max()}")
print(f"Time range: {exec_times.min():.1f}-{exec_times.max():.1f} ms")

# Define exponential model
def exponential_model(x, a, b):
    """Exponential model: t = a * exp(b * x)"""
    return a * np.exp(b * x)

# Fit the exponential model
# Initial guess: a=50, b=0.015 (based on expected exponential growth)
popt, pcov = curve_fit(exponential_model, atom_counts, exec_times, 
                       p0=[50, 0.015], maxfev=10000)

a_fit, b_fit = popt
print(f"\nExponential Model: t = {a_fit:.3f} * exp({b_fit:.5f} * n)")

# Calculate R² and RMSE
y_pred_data = exponential_model(atom_counts, *popt)
ss_res = np.sum((exec_times - y_pred_data) ** 2)
ss_tot = np.sum((exec_times - np.mean(exec_times)) ** 2)
r_squared = 1 - (ss_res / ss_tot)
rmse = np.sqrt(np.mean((exec_times - y_pred_data) ** 2))

print(f"R² = {r_squared:.4f}, RMSE = {rmse:.1f} ms")

# Calculate standard errors from covariance matrix
perr = np.sqrt(np.diag(pcov))
print(f"Standard errors: a = {perr[0]:.3f}, b = {perr[1]:.6f}")

# Generate prediction grid for smooth curve
x_grid = np.linspace(atom_counts.min(), atom_counts.max(), 300)
y_pred = exponential_model(x_grid, *popt)

# Calculate 95% confidence intervals using error propagation
# For y = a * exp(b * x), we have:
# dy/da = exp(b * x)
# dy/db = a * x * exp(b * x)
# 
# Variance propagation:
# var(y) = (dy/da)^2 * var(a) + (dy/db)^2 * var(b) + 2*(dy/da)*(dy/db)*cov(a,b)

dy_da = np.exp(b_fit * x_grid)
dy_db = a_fit * x_grid * np.exp(b_fit * x_grid)

var_y = (dy_da**2 * pcov[0, 0] + 
         dy_db**2 * pcov[1, 1] + 
         2 * dy_da * dy_db * pcov[0, 1])

se_y = np.sqrt(var_y)

# Use t-distribution for 95% confidence interval
dof = len(atom_counts) - 2  # degrees of freedom
t_val = stats.t.ppf(0.975, dof)  # 97.5th percentile for two-tailed test

ci_lower = y_pred - t_val * se_y
ci_upper = y_pred + t_val * se_y

print(f"\n95% Confidence Interval Statistics:")
for test_n in [20, 50, 100]:
    idx = np.argmin(np.abs(x_grid - test_n))
    print(f"At n={test_n}: {y_pred[idx]:.1f} ms ± {t_val * se_y[idx]:.1f} ms")

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot data points
ax.scatter(atom_counts, exec_times, alpha=0.6, s=30, color='#667eea', 
           label='Measured times', zorder=3)

# Plot fitted curve
ax.plot(x_grid, y_pred, 'r-', linewidth=2, label='Exponential fit', zorder=2)

# Plot confidence interval
ax.fill_between(x_grid, ci_lower, ci_upper, alpha=0.2, color='red',
                label='95% confidence interval', zorder=1)

# Add grid
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

# Labels and title
ax.set_xlabel('Number of Atoms', fontsize=12, fontweight='bold')
ax.set_ylabel('Execution Time (ms)', fontsize=12, fontweight='bold')
ax.set_title('HalogenGroups Execution Time: Exponential Scaling with Molecular Size', 
             fontsize=14, fontweight='bold', pad=20)

# Add statistics text box
textstr = f'Model: $t = {a_fit:.1f} \\times \\exp({b_fit:.4f} \\times n)$\n'
textstr += f'$R^2 = {r_squared:.4f}$\n'
textstr += f'RMSE = {rmse:.1f} ms\n'
textstr += f'$n = {len(atom_counts)}$ measurements\n'
textstr += f'$\\alpha = {b_fit:.4f}$ atoms$^{{-1}}$'

props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

# Legend
ax.legend(loc='lower right', framealpha=0.9, fontsize=10)

# Set y-axis to start at 0 for better visualization
ax.set_ylim(bottom=0)

# Tight layout
plt.tight_layout()

# Save figure
output_pdf = 'timing_exponential_fit_v2.pdf'
output_png = 'timing_exponential_fit_v2.png'
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
plt.savefig(output_png, dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {output_pdf}")
print(f"PNG version saved to: {output_png}")

plt.close()
