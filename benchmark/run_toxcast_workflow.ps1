<#
.SYNOPSIS
    Run the full ToxCast fingerprint comparison workflow.

.DESCRIPTION
    Steps:
      0. Ask which database / dataset to use (default: invitrodb_v4_3).
      1. (Optional) Rebuild the ToxCast dataset parquet from MariaDB.
      2. Run the fingerprint comparison: nested CV, Random Forests,
         Gradient Boosting, PFASGroups vs Richard2023.
      3. Print a summary of all output files.

    The dataset rebuild step requires a running MariaDB ToxCast instance and
    is skipped when --skip-build is passed or when the parquet already exists.

.PARAMETER SkipBuild
    Skip the dataset build step (use an existing <db>_dataset.parquet).

.PARAMETER Env
    Name of the conda environment to use (default: chem).

.PARAMETER DB
    Database / dataset name to use (default: invitrodb_v4_3).
    Overrides the interactive prompt when provided.

.EXAMPLE
    # Full run (prompts for DB, builds dataset then runs comparison):
    .\run_toxcast_workflow.ps1

    # Skip dataset rebuild, use a specific DB non-interactively:
    .\run_toxcast_workflow.ps1 -SkipBuild -DB invitrodb_v4_3

    # Use a different conda env:
    .\run_toxcast_workflow.ps1 -Env myenv
#>
param(
    [switch]$SkipBuild,
    [string]$Env = "chem",
    [string]$DB  = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
function Write-Step([string]$msg) {
    Write-Host ""
    Write-Host "=== $msg ===" -ForegroundColor Cyan
}

function Invoke-CondaPython([string]$script, [string[]]$PyArgs, [string]$label) {
    Write-Step $label
    $t0 = Get-Date
    conda run -n $Env python $script @PyArgs
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ERROR] '$script' exited with code $LASTEXITCODE" -ForegroundColor Red
        exit $LASTEXITCODE
    }
    $elapsed = (Get-Date) - $t0
    Write-Host "[OK] finished in $($elapsed.ToString('mm\:ss'))" -ForegroundColor Green
}

# ---------------------------------------------------------------------------
# Locate benchmark directory
# ---------------------------------------------------------------------------
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$DataDir   = Join-Path $ScriptDir "data"

Write-Host "Benchmark dir : $ScriptDir" -ForegroundColor DarkGray
Write-Host "Conda env     : $Env"       -ForegroundColor DarkGray

# ---------------------------------------------------------------------------
# Step 0 – choose database / dataset
# ---------------------------------------------------------------------------
if ($DB -eq "") {
    $DB = Read-Host "Database / dataset name [default: invitrodb_v4_3]"
    if ($DB -eq "") { $DB = "invitrodb_v4_3" }
}
Write-Host "Using database/dataset : $DB" -ForegroundColor DarkGray

$ParquetPath = Join-Path $DataDir "${DB}_dataset.parquet"

# ---------------------------------------------------------------------------
# Step 1 – dataset build
# ---------------------------------------------------------------------------
if ($SkipBuild) {
    Write-Step "Skipping dataset build (-SkipBuild)"
    if (-not (Test-Path $ParquetPath)) {
        Write-Host "[WARN] $ParquetPath not found. The comparison step will fail." -ForegroundColor Yellow
    }
} elseif (Test-Path $ParquetPath) {
    Write-Host ""
    Write-Host "Dataset already exists at $ParquetPath" -ForegroundColor DarkGray
    $answer = Read-Host "Rebuild it? [y/N]"
    if ($answer -match '^[Yy]') {
        Invoke-CondaPython (Join-Path $ScriptDir "scripts\data\build_toxcast_dataset.py") @("--db", $DB) "Building ToxCast dataset ($DB)"
    } else {
        Write-Step "Using existing dataset"
    }
} else {
    Invoke-CondaPython (Join-Path $ScriptDir "scripts\data\build_toxcast_dataset.py") @("--db", $DB) "Building ToxCast dataset ($DB)"
}

# ---------------------------------------------------------------------------
# Step 2 – fingerprint comparison
# ---------------------------------------------------------------------------
Invoke-CondaPython (Join-Path $ScriptDir "scripts\analysis\compare_fingerprints_toxcast.py") @("--db", $DB, "--both") "Running fingerprint comparison (nested CV)"

# ---------------------------------------------------------------------------
# Step 3 – summary
# ---------------------------------------------------------------------------
Write-Step "Output files"
if (Test-Path $DataDir) {
    Get-ChildItem $DataDir -Filter "${DB}_comparison_*" |
        Sort-Object Name |
        ForEach-Object { Write-Host "  $($_.Name)  ($([math]::Round($_.Length/1KB, 1)) KB)" }
} else {
    Write-Host "[WARN] data directory not found: $DataDir" -ForegroundColor Yellow
}

# ---------------------------------------------------------------------------
# Step 4 – performance plots
# ---------------------------------------------------------------------------
Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\plot_performance_by_endpoint.py") @("--dataset", $DB) "Generating performance plots"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\plot_expa_accuracy_precision.py") @("--dataset", $DB) "Generating accuracy and precision plots"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\plot_expa_filtered.py") @("--dataset", $DB) "Generating filtered Exp A bar plots"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\bayesian_comparison.py") @("--dataset", $DB) "Generating Bayesian comparison plots"

# ---------------------------------------------------------------------------
# Step 5 – Unsupervised analysis reports
# ---------------------------------------------------------------------------
Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @() "Generating UMAP plots  --top 20"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "umap", "--fp", "full", "--top", "20") "Generating UMAP plots for fp EGR+cycle+nspacer+branch+mol"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "tsne", "--fp", "umap", "--top", "20") "Generating UMAP plots for fp total_component+spacer+ring+mol"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "tsne", "--fp", "binary", "--top", "20") "Generating t-SNE plots"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "tsne", "--fp", "full", "--top", "20") "Generating t-SNE plots for fp EGR+cycle+nspacer+branch+mol"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "tsne", "--fp", "total", "--top", "20") "Generating t-SNE plots for fp total_component+spacer+ring+mol"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "mds", "--fp", "binary", "--top", "20") "Generating MDS plots"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "mds", "--fp", "full", "--top", "20") "Generating MDS plots for fp EGR+cycle+nspacer+branch+mol"

Invoke-CondaPython (Join-Path $ScriptDir "scripts\plots\umap_oecd_pfas.py") @("--method", "mds", "--fp", "total", "--top", "20") "Generating MDS plots for fp total_component+spacer+ring+mol"

Write-Host ""
Write-Host "Workflow complete." -ForegroundColor Green
