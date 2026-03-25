# PFAS Definitions Benchmark Runner
# Comprehensive script to run definitions benchmark and generate reports

Write-Host "`n================================================================================" -ForegroundColor Cyan
Write-Host "   PFAS DEFINITIONS COMPREHENSIVE BENCHMARK SUITE" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan

Write-Host "`nThis script will:" -ForegroundColor Yellow
Write-Host "  1. Run comprehensive PFAS definitions benchmark (~10-15 min)" -ForegroundColor White
Write-Host "  2. Generate analysis and visualizations (~2-3 min)" -ForegroundColor White
Write-Host "  3. Create HTML report with all results" -ForegroundColor White
Write-Host "  4. Open report in browser" -ForegroundColor White

# Check if we're in the right directory
if (-not (Test-Path "scripts\classify\benchmark_pfas_definitions.py")) {
    Write-Host "`n❌ Error: scripts\classify\benchmark_pfas_definitions.py not found" -ForegroundColor Red
    Write-Host "Please run this script from the benchmark directory:" -ForegroundColor Yellow
    Write-Host "  cd c:\Users\luc\git\PFASGroups\benchmark" -ForegroundColor White
    exit 1
}

# Check Python environment
Write-Host "`n🔍 Checking Python environment..." -ForegroundColor Cyan
python --version
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Python not found. Please activate conda environment:" -ForegroundColor Red
    Write-Host "  conda activate chem" -ForegroundColor White
    exit 1
}

# Check required packages
Write-Host "`n🔍 Checking required packages..." -ForegroundColor Cyan
$packages = @("rdkit", "pandas", "numpy", "matplotlib", "plotly")
$missing = @()

foreach ($package in $packages) {
    python -c "import $package" 2>$null
    if ($LASTEXITCODE -ne 0) {
        $missing += $package
    }
}

if ($missing.Count -gt 0) {
    Write-Host "⚠️  Missing packages: $($missing -join ', ')" -ForegroundColor Yellow
    $install = Read-Host "Install missing packages? (y/n)"
    if ($install -eq 'y') {
        pip install $missing
    } else {
        Write-Host "❌ Cannot proceed without required packages" -ForegroundColor Red
        exit 1
    }
}

Write-Host "✅ All dependencies satisfied" -ForegroundColor Green

# Prompt to continue
Write-Host "`n" -NoNewline
$continue = Read-Host "Press Enter to start benchmark or 'q' to quit"
if ($continue -eq 'q') {
    Write-Host "Benchmark cancelled" -ForegroundColor Yellow
    exit 0
}

# Create data directory if it doesn't exist
if (-not (Test-Path "data")) {
    New-Item -ItemType Directory -Path "data" | Out-Null
    Write-Host "✅ Created data directory" -ForegroundColor Green
}

# Step 1: Run benchmark
Write-Host "`n================================================================================" -ForegroundColor Cyan
Write-Host "STEP 1: Running PFAS Definitions Benchmark" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host "Testing:" -ForegroundColor Yellow
Write-Host "  • True Positives (~30 known PFAS compounds)" -ForegroundColor White
Write-Host "  • True Negatives (~15 non-PFAS compounds)" -ForegroundColor White
Write-Host "  • Edge Cases (~15 borderline compounds)" -ForegroundColor White
Write-Host "  • Definition-Specific Tests (~15 targeted tests)" -ForegroundColor White
Write-Host "  • Concordance Analysis (~10 reference compounds)" -ForegroundColor White
Write-Host "  • Performance Testing (~100 random molecules)" -ForegroundColor White
Write-Host "`nEstimated time: 10-15 minutes`n" -ForegroundColor Yellow

$benchmark_start = Get-Date

python scripts\classify\benchmark_pfas_definitions.py

if ($LASTEXITCODE -ne 0) {
    Write-Host "`n❌ Benchmark failed with error code $LASTEXITCODE" -ForegroundColor Red
    exit 1
}

$benchmark_duration = (Get-Date) - $benchmark_start
Write-Host "`n✅ Benchmark completed in $([math]::Round($benchmark_duration.TotalMinutes, 1)) minutes" -ForegroundColor Green

# Find the latest benchmark file
$latest_benchmark = Get-ChildItem -Path "data" -Filter "pfas_definitions_benchmark_*.json" | 
                    Sort-Object LastWriteTime -Descending | 
                    Select-Object -First 1

if (-not $latest_benchmark) {
    Write-Host "`n❌ No benchmark results file found" -ForegroundColor Red
    exit 1
}

Write-Host "📁 Results saved to: $($latest_benchmark.Name)" -ForegroundColor White

# Step 2: Generate analysis
Write-Host "`n================================================================================" -ForegroundColor Cyan
Write-Host "STEP 2: Generating Analysis and Visualizations" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host "Creating:" -ForegroundColor Yellow
Write-Host "  • Summary statistics table" -ForegroundColor White
Write-Host "  • Sensitivity vs Specificity plot" -ForegroundColor White
Write-Host "  • Detection heatmap by category" -ForegroundColor White
Write-Host "  • Inter-definition concordance matrix" -ForegroundColor White
Write-Host "  • Edge case disagreement analysis" -ForegroundColor White
Write-Host "  • Performance analysis" -ForegroundColor White
Write-Host "  • Comprehensive HTML report" -ForegroundColor White
Write-Host "`nEstimated time: 2-3 minutes`n" -ForegroundColor Yellow

$analysis_start = Get-Date

python scripts\analysis\analyze_definitions_benchmark.py $latest_benchmark.FullName

if ($LASTEXITCODE -ne 0) {
    Write-Host "`n❌ Analysis failed with error code $LASTEXITCODE" -ForegroundColor Red
    exit 1
}

$analysis_duration = (Get-Date) - $analysis_start
Write-Host "`n✅ Analysis completed in $([math]::Round($analysis_duration.TotalSeconds, 0)) seconds" -ForegroundColor Green

# Find the HTML report
$html_report = Get-ChildItem -Path "data" -Filter "definitions_benchmark_report_*.html" | 
               Sort-Object LastWriteTime -Descending | 
               Select-Object -First 1

if (-not $html_report) {
    Write-Host "`n⚠️  HTML report not found" -ForegroundColor Yellow
} else {
    Write-Host "📄 Report saved to: $($html_report.Name)" -ForegroundColor White
}

# Summary
$total_duration = (Get-Date) - $benchmark_start

Write-Host "`n================================================================================" -ForegroundColor Cyan
Write-Host "   BENCHMARK COMPLETE!" -ForegroundColor Green
Write-Host "================================================================================" -ForegroundColor Cyan

Write-Host "`n📊 Summary:" -ForegroundColor Yellow
Write-Host "  • Total time: $([math]::Round($total_duration.TotalMinutes, 1)) minutes" -ForegroundColor White
Write-Host "  • Benchmark data: data\$($latest_benchmark.Name)" -ForegroundColor White
if ($html_report) {
    Write-Host "  • HTML report: data\$($html_report.Name)" -ForegroundColor White
}
Write-Host "  • Figures: data\figures_definitions\" -ForegroundColor White

Write-Host "`n📂 Output Files:" -ForegroundColor Yellow
Write-Host "  • Raw JSON results for programmatic access" -ForegroundColor White
Write-Host "  • Interactive HTML visualizations" -ForegroundColor White
Write-Host "  • High-resolution PNG images" -ForegroundColor White

# Offer to open report
if ($html_report) {
    Write-Host "`n" -NoNewline
    $open_report = Read-Host "Open HTML report in browser? (y/n)"
    if ($open_report -eq 'y') {
        Start-Process $html_report.FullName
        Write-Host "✅ Report opened in default browser" -ForegroundColor Green
    }
}

# Offer to open figures directory
Write-Host "`n" -NoNewline
$open_figures = Read-Host "Open figures directory? (y/n)"
if ($open_figures -eq 'y') {
    Start-Process "data\figures_definitions"
    Write-Host "✅ Figures directory opened" -ForegroundColor Green
}

Write-Host "`n================================================================================" -ForegroundColor Cyan
Write-Host "Next Steps:" -ForegroundColor Yellow
Write-Host "  1. Review HTML report for detailed results" -ForegroundColor White
Write-Host "  2. Check sensitivity/specificity for each definition" -ForegroundColor White
Write-Host "  3. Examine edge cases where definitions disagree" -ForegroundColor White
Write-Host "  4. Compare performance metrics" -ForegroundColor White
Write-Host "  5. Document any unexpected findings" -ForegroundColor White
Write-Host "`nSee PFAS_DEFINITIONS_BENCHMARKING_GUIDE.md for interpretation help" -ForegroundColor Cyan
Write-Host "================================================================================" -ForegroundColor Cyan

Write-Host "`n✨ Benchmark suite completed successfully!`n" -ForegroundColor Green
