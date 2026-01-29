#!/usr/bin/env pwsh
# Update all script references to match new directory structure
# New structure:
# - scripts/     : All Python scripts
# - reports/     : All HTML reports  
# - data/        : All JSON, CSV, NPZ data files
# - imgs/        : All PNG, PDF, SVG images

Write-Host "🔧 Updating script path references..." -ForegroundColor Cyan
Write-Host ""

$ErrorActionPreference = "Stop"

# =============================================================================
# Update run_all_benchmarks.sh
# =============================================================================
Write-Host "📝 Updating run_all_benchmarks.sh..." -ForegroundColor Yellow

$bashContent = Get-Content "run_all_benchmarks.sh" -Raw

# Update check for script existence
$bashContent = $bashContent -replace 'if \[ ! -f "enhanced_pfas_benchmark\.py" \]', 'if [ ! -f "scripts/enhanced_pfas_benchmark.py" ]'
$bashContent = $bashContent -replace 'echo "❌ Error: enhanced_pfas_benchmark\.py not found!"', 'echo "❌ Error: scripts/enhanced_pfas_benchmark.py not found!"'

# Update directory creation
$bashContent = $bashContent -replace 'mkdir -p imgs html data', 'mkdir -p scripts reports data imgs'

# Update all python script calls
$bashContent = $bashContent -replace 'python enhanced_pfas_benchmark\.py', 'python scripts/enhanced_pfas_benchmark.py'
$bashContent = $bashContent -replace 'python test_highly_branched\.py', 'python scripts/test_highly_branched.py'
$bashContent = $bashContent -replace 'python generate_unified_report\.py', 'python scripts/generate_unified_report.py'
$bashContent = $bashContent -replace 'python generate_timing_report\.py', 'python scripts/generate_timing_report.py'

# Update file organization - move to reports instead of html
$bashContent = $bashContent -replace 'mv \*\.html html/', 'mv *.html reports/'

Set-Content "run_all_benchmarks.sh" -Value $bashContent
Write-Host "  ✓ run_all_benchmarks.sh updated" -ForegroundColor Green

# =============================================================================
# Update run_all_benchmarks.ps1
# =============================================================================
Write-Host "📝 Updating run_all_benchmarks.ps1..." -ForegroundColor Yellow

$ps1Content = Get-Content "run_all_benchmarks.ps1" -Raw

# Update check for script existence
$ps1Content = $ps1Content -replace 'if \(-not \(Test-Path "enhanced_pfas_benchmark\.py"\)\)', 'if (-not (Test-Path "scripts/enhanced_pfas_benchmark.py"))'
$ps1Content = $ps1Content -replace 'Write-Host "❌ Error: enhanced_pfas_benchmark\.py not found!"', 'Write-Host "❌ Error: scripts/enhanced_pfas_benchmark.py not found!"'

# Update directory creation - change html to reports
$ps1Content = $ps1Content -replace 'New-Item -ItemType Directory -Force -Path imgs', 'New-Item -ItemType Directory -Force -Path scripts, reports, data, imgs'
$ps1Content = $ps1Content -replace 'New-Item -ItemType Directory -Force -Path html.*\r?\n', ''
$ps1Content = $ps1Content -replace 'New-Item -ItemType Directory -Force -Path data.*\r?\n', ''

# Update python script calls
$ps1Content = $ps1Content -replace '\$input \| python enhanced_pfas_benchmark\.py', '$input | python scripts/enhanced_pfas_benchmark.py'
$ps1Content = $ps1Content -replace 'python test_highly_branched\.py', 'python scripts/test_highly_branched.py'
$ps1Content = $ps1Content -replace 'python generate_unified_report\.py', 'python scripts/generate_unified_report.py'
$ps1Content = $ps1Content -replace 'python generate_timing_report\.py', 'python scripts/generate_timing_report.py'

# Update file organization - reports instead of html
$ps1Content = $ps1Content -replace 'Move-Item \*\.html html\\', 'Move-Item *.html reports\'

Set-Content "run_all_benchmarks.ps1" -Value $ps1Content
Write-Host "  ✓ run_all_benchmarks.ps1 updated" -ForegroundColor Green

# =============================================================================
# Update quick-start.ps1
# =============================================================================
Write-Host "📝 Updating quick-start.ps1..." -ForegroundColor Yellow

$quickContent = Get-Content "quick-start.ps1" -Raw

# Update python script references
$quickContent = $quickContent -replace 'python enhanced_pfas_benchmark\.py', 'python scripts/enhanced_pfas_benchmark.py'
$quickContent = $quickContent -replace 'python generate_timing_report\.py benchmark_timing_report\.html', 'python scripts/generate_timing_report.py reports/benchmark_timing_report.html'

Set-Content "quick-start.ps1" -Value $quickContent
Write-Host "  ✓ quick-start.ps1 updated" -ForegroundColor Green

# =============================================================================
# Update generate_unified_report.py
# =============================================================================
Write-Host "📝 Updating scripts/generate_unified_report.py..." -ForegroundColor Yellow

$genReportContent = Get-Content "scripts/generate_unified_report.py" -Raw

# Update file search patterns to look in data/ directory
$genReportContent = $genReportContent -replace "glob\.glob\('pfas_", "glob.glob('data/pfas_"
$genReportContent = $genReportContent -replace "glob\.glob\(pattern\)", "glob.glob('data/' + pattern.split('/')[-1])"

# Update output path to reports/
$genReportContent = $genReportContent -replace 'shutil\.move\(html_filename, f"html/\{html_filename\}"\)', 'shutil.move(html_filename, f"reports/{html_filename}")'
$genReportContent = $genReportContent -replace 'html_filename = f"html/\{html_filename\}"', 'html_filename = f"reports/{html_filename}"'

# Update image save paths to use imgs/
$genReportContent = $genReportContent -replace "plt\.savefig\('", "plt.savefig('imgs/"
$genReportContent = $genReportContent -replace "\.write_image\('", ".write_image('imgs/"

Set-Content "scripts/generate_unified_report.py" -Value $genReportContent
Write-Host "  ✓ scripts/generate_unified_report.py updated" -ForegroundColor Green

# =============================================================================
# Update generate_timing_report.py
# =============================================================================
Write-Host "📝 Updating scripts/generate_timing_report.py..." -ForegroundColor Yellow

$timingReportContent = Get-Content "scripts/generate_timing_report.py" -Raw

# Update default output path
$timingReportContent = $timingReportContent -replace "output_file='benchmark_timing_report\.html'", "output_file='reports/benchmark_timing_report.html'"
$timingReportContent = $timingReportContent -replace "output_file = sys\.argv\[1\] if len\(sys\.argv\) > 1 else 'benchmark_timing_report\.html'", "output_file = sys.argv[1] if len(sys.argv) > 1 else 'reports/benchmark_timing_report.html'"

# Update data file paths
$timingReportContent = $timingReportContent -replace "'timing_analysis\.json'", "'data/timing_analysis.json'"
$timingReportContent = $timingReportContent -replace "'pfas_timing_benchmark_", "'data/pfas_timing_benchmark_"

Set-Content "scripts/generate_timing_report.py" -Value $timingReportContent
Write-Host "  ✓ scripts/generate_timing_report.py updated" -ForegroundColor Green

# =============================================================================
# Update enhanced_pfas_benchmark.py
# =============================================================================
Write-Host "📝 Updating scripts/enhanced_pfas_benchmark.py..." -ForegroundColor Yellow

$enhancedContent = Get-Content "scripts/enhanced_pfas_benchmark.py" -Raw

# Update output file paths to use data/ directory
$enhancedContent = $enhancedContent -replace "output_file = f'pfas_", "output_file = f'data/pfas_"
$enhancedContent = $enhancedContent -replace 'with open\(f"pfas_', 'with open(f"data/pfas_'

# Update HTML report paths to use reports/ directory  
$enhancedContent = $enhancedContent -replace "f'enhanced_pfas_analysis_", "f'reports/enhanced_pfas_analysis_"
$enhancedContent = $enhancedContent -replace "f'complex_benchmark_report_", "f'reports/complex_benchmark_report_"

# Update image paths to use imgs/ directory
$enhancedContent = $enhancedContent -replace "plt\.savefig\(f'enhanced_timing_report_", "plt.savefig(f'imgs/enhanced_timing_report_"

Set-Content "scripts/enhanced_pfas_benchmark.py" -Value $enhancedContent
Write-Host "  ✓ scripts/enhanced_pfas_benchmark.py updated" -ForegroundColor Green

# =============================================================================
# Update other key Python scripts
# =============================================================================
Write-Host "📝 Updating other Python scripts..." -ForegroundColor Yellow

$scriptFiles = @(
    "scripts/generate_enhanced_timing_report.py",
    "scripts/export_timing_figures.py",
    "scripts/analyze_timing.py"
)

foreach ($scriptFile in $scriptFiles) {
    if (Test-Path $scriptFile) {
        $content = Get-Content $scriptFile -Raw
        
        # Update data paths
        $content = $content -replace "'pfas_timing_benchmark_", "'data/pfas_timing_benchmark_"
        $content = $content -replace "'timing_analysis\.json'", "'data/timing_analysis.json'"
        $content = $content -replace '"pfas_timing_benchmark_', '"data/pfas_timing_benchmark_'
        $content = $content -replace '"timing_analysis\.json"', '"data/timing_analysis.json"'
        
        # Update image paths
        $content = $content -replace "plt\.savefig\('enhanced_timing_", "plt.savefig('imgs/enhanced_timing_"
        $content = $content -replace "\.write_image\('enhanced_timing_", ".write_image('imgs/enhanced_timing_"
        
        # Update HTML report paths
        $content = $content -replace "'enhanced_timing_report\.html'", "'reports/enhanced_timing_report.html'"
        $content = $content -replace '"enhanced_timing_report\.html"', '"reports/enhanced_timing_report.html"'
        
        Set-Content $scriptFile -Value $content
        Write-Host "  ✓ $scriptFile updated" -ForegroundColor Green
    }
}

# =============================================================================
# Update review-app import script reference
# =============================================================================
Write-Host "📝 Updating review-app import script..." -ForegroundColor Yellow

if (Test-Path "review-app/scripts/import-benchmark-data.js") {
    $importContent = Get-Content "review-app/scripts/import-benchmark-data.js" -Raw
    
    # Update data directory path
    $importContent = $importContent -replace "path\.join\(__dirname, '\.\./\.\./data'", "path.join(__dirname, '../../data'"
    $importContent = $importContent -replace "glob\.sync\('pfas_", "glob.sync('../data/pfas_"
    
    Set-Content "review-app/scripts/import-benchmark-data.js" -Value $importContent
    Write-Host "  ✓ review-app/scripts/import-benchmark-data.js updated" -ForegroundColor Green
}

Write-Host ""
Write-Host "✅ All path references updated!" -ForegroundColor Green
Write-Host ""
Write-Host "📊 Updated structure:" -ForegroundColor Cyan
Write-Host "  scripts/  - Python scripts (*.py)" -ForegroundColor White
Write-Host "  reports/  - HTML reports (*.html)" -ForegroundColor White
Write-Host "  data/     - Data files (*.json, *.csv, *.npz)" -ForegroundColor White
Write-Host "  imgs/     - Images (*.png, *.pdf, *.svg)" -ForegroundColor White
Write-Host ""
Write-Host "🎯 Next steps:" -ForegroundColor Cyan
Write-Host "  1. Test the reorganization:" -ForegroundColor White
Write-Host "     cd review-app && node scripts/import-benchmark-data.js" -ForegroundColor Gray
Write-Host "  2. Run benchmarks to verify:" -ForegroundColor White  
Write-Host "     .\run_all_benchmarks.ps1" -ForegroundColor Gray
Write-Host ""
