# Run all analysis and report-generation scripts
# Regenerates all figures, HTML reports, and LaTeX tables from existing benchmark data.
# Does NOT re-run data collection (classify_* / benchmark_* scripts).
#
# Usage (from the benchmark/ directory):
#   .\run_analysis_reports.ps1
#
# Requires the following conda environments:
#   pfasgroups  – for most analysis scripts
#   pfasatlas   – for compare_clinventory_classifiers.py and generate_clinventory_latex.py
#
# Both environments must be accessible via `conda run`.

param(
    [string]$Env       = "chem",
    [string]$AtlasEnv  = "chem",
    [switch]$SkipAtlas                   # pass -SkipAtlas to skip scripts that need pfasatlas
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Continue"     # keep going even if one script fails

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
function Step {
    param([string]$Label)
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "  $Label" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
}

function Run-Script {
    param(
        [string]$Env,
        [string]$Script,
        [string[]]$ScriptArgs = @()
    )
    $cmd = @("python", "scripts\$Script") + $ScriptArgs
    Write-Host "> conda run -n $Env $($cmd -join ' ')" -ForegroundColor DarkGray
    conda run -n $Env @cmd
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  [FAILED  exit $LASTEXITCODE]" -ForegroundColor Red
    } else {
        Write-Host "  [OK]" -ForegroundColor Green
    }
}

# ---------------------------------------------------------------------------
# Guard: must be run from benchmark/
# ---------------------------------------------------------------------------
if (-not (Test-Path "scripts\analysis\analyze_timing.py")) {
    Write-Host "ERROR: Run this script from the benchmark/ directory." -ForegroundColor Red
    exit 1
}

# ---------------------------------------------------------------------------
# Resolve the latest timing benchmark file (required by analyze_timing.py)
# ---------------------------------------------------------------------------
$timingFile = Get-ChildItem "data\pfas_timing_benchmark_full_*.json","data\pfas_timing_benchmark_*.json" `
    -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1

if (-not $timingFile) {
    Write-Host "WARNING: No timing benchmark file found in data/ — analyze_timing.py will be skipped." -ForegroundColor Yellow
}

$defsBenchFile = Get-ChildItem "data\pfas_definitions_benchmark_*.json" `
    -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1

# ---------------------------------------------------------------------------
# 1. Timing analysis  (needs explicit file argument)
# ---------------------------------------------------------------------------
Step "1 / 11  analyze_timing.py"
if ($timingFile) {
    Run-Script -Env $Env -Script "analyze_timing.py" -ScriptArgs @($timingFile.FullName)
} else {
    Write-Host "  SKIPPED (no data file)" -ForegroundColor Yellow
}

# ---------------------------------------------------------------------------
# 2. Timing scaling model
# ---------------------------------------------------------------------------
Step "2 / 11  analyze_timing_models.py"
Run-Script -Env $Env -Script "analyze_timing_models.py"

# ---------------------------------------------------------------------------
# 3. Timing profiles comparison
# ---------------------------------------------------------------------------
Step "3 / 11  compare_timing_profiles.py"
Run-Script -Env $Env -Script "compare_timing_profiles.py"

# ---------------------------------------------------------------------------
# 4. Definitions benchmark analysis  (auto-detects latest file)
# ---------------------------------------------------------------------------
Step "4 / 11  analyze_definitions_benchmark.py"
if ($defsBenchFile) {
    Run-Script -Env $Env -Script "analyze_definitions_benchmark.py" -ScriptArgs @($defsBenchFile.FullName)
} else {
    Run-Script -Env $Env -Script "analyze_definitions_benchmark.py"
}

# ---------------------------------------------------------------------------
# 5. Enhanced analysis  (auto-detects latest enhanced + OECD files)
# ---------------------------------------------------------------------------
Step "5 / 11  enhanced_analysis.py"
Run-Script -Env $Env -Script "enhanced_analysis.py" -ScriptArgs @("--combined")

# ---------------------------------------------------------------------------
# 6. Unified HTML report  (auto-detects all benchmark files)
# ---------------------------------------------------------------------------
Step "6 / 11  generate_unified_report.py"
Run-Script -Env $Env -Script "generate_unified_report.py"

# ---------------------------------------------------------------------------
# 7. Telomer HTML report
# ---------------------------------------------------------------------------
Step "7 / 11  generate_telomer_report.py"
Run-Script -Env $Env -Script "generate_telomer_report.py"

# ---------------------------------------------------------------------------
# 8. LaTeX tables  (needs pfasgroups; uses data/ benchmark files)
# ---------------------------------------------------------------------------
Step "8 / 11  generate_latex_tables.py"
Run-Script -Env $Env -Script "generate_latex_tables.py"

# ---------------------------------------------------------------------------
# 9. Telomer LaTeX tables
# ---------------------------------------------------------------------------
Step "9 / 11  generate_telomer_latex.py"
Run-Script -Env $Env -Script "generate_telomer_latex.py"

# ---------------------------------------------------------------------------
# 10 & 11. Clinventory: compare classifiers + LaTeX tables
#          These scripts read clinventory JSON produced by classify_*.py
# ---------------------------------------------------------------------------
if ($SkipAtlas) {
    Write-Host ""
    Write-Host "Skipping clinventory scripts (-SkipAtlas was set)." -ForegroundColor Yellow
} else {
    Step "10 / 11  compare_clinventory_classifiers.py  [env: $AtlasEnv]"
    Run-Script -Env $AtlasEnv -Script "compare_clinventory_classifiers.py"

    Step "11 / 11  generate_clinventory_latex.py  [env: $AtlasEnv]"
    Run-Script -Env $AtlasEnv -Script "generate_clinventory_latex.py"
}

# ---------------------------------------------------------------------------
# 12. OECD × clinventory timing + complexity comparison
#     Classifies fresh on the OECD CSV; connects to clinventory if available.
# ---------------------------------------------------------------------------
Step "12 / 12  compare_oecd_clinventory_timing.py  [env: $AtlasEnv]"
Run-Script -Env $AtlasEnv -Script "compare_oecd_clinventory_timing.py"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
Write-Host ""
Write-Host ("=" * 70) -ForegroundColor Cyan
Write-Host "  Done.  Reports written to benchmark/reports/ and benchmark/imgs/" -ForegroundColor Cyan
Write-Host ("=" * 70) -ForegroundColor Cyan
