# Generate Analysis Reports for PFAS Benchmarks
# Run this after completing all benchmarks

Write-Host "📊 PFAS Analysis Report Generator" -ForegroundColor Cyan
Write-Host ("=" * 50) -ForegroundColor Cyan
Write-Host ""

# Check if we're in the right directory
if (-not (Test-Path "analyze_timing.py")) {
    Write-Host "❌ Error: Analysis scripts not found!" -ForegroundColor Red
    Write-Host "   Please run this script from the benchmark directory" -ForegroundColor Red
    exit 1
}

# Find the latest benchmark files
$timingFile = Get-ChildItem "data\pfas_timing_benchmark_*.json" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
$complexFile = Get-ChildItem "data\pfas_complex_branched_benchmark_*.json" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
$enhancedFile = Get-ChildItem "data\pfas_enhanced_benchmark_*.json" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
$oecdFile = Get-ChildItem "data\pfas_oecd_benchmark_*.json" | Sort-Object LastWriteTime -Descending | Select-Object -First 1

Write-Host "📁 Found benchmark files:" -ForegroundColor Yellow
if ($timingFile) { Write-Host "  ✅ Timing: $($timingFile.Name)" -ForegroundColor Green } else { Write-Host "  ❌ Timing: Not found" -ForegroundColor Red }
if ($complexFile) { Write-Host "  ✅ Complex: $($complexFile.Name)" -ForegroundColor Green } else { Write-Host "  ❌ Complex: Not found" -ForegroundColor Red }
if ($enhancedFile) { Write-Host "  ✅ Enhanced: $($enhancedFile.Name)" -ForegroundColor Green } else { Write-Host "  ❌ Enhanced: Not found" -ForegroundColor Red }
if ($oecdFile) { Write-Host "  ✅ OECD: $($oecdFile.Name)" -ForegroundColor Green } else { Write-Host "  ⚠️  OECD: Not found (run option 2 in enhanced_pfas_benchmark.py)" -ForegroundColor Yellow }
Write-Host ""

# Create analysis_reports directory
New-Item -ItemType Directory -Force -Path review-app\analysis_reports\figures | Out-Null

# Run timing analysis
if ($timingFile) {
    Write-Host "⏱️  Generating timing analysis..." -ForegroundColor Yellow
    python analyze_timing.py $timingFile.FullName
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Timing analysis completed" -ForegroundColor Green
    } else {
        Write-Host "❌ Timing analysis failed" -ForegroundColor Red
    }
    Write-Host ""
}

# Run complex analysis
if ($complexFile) {
    Write-Host "🧬 Generating complex branched analysis..." -ForegroundColor Yellow
    python analyze_complex.py $complexFile.FullName
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Complex analysis completed" -ForegroundColor Green
    } else {
        Write-Host "❌ Complex analysis failed" -ForegroundColor Red
    }
    Write-Host ""
}

# Run enhanced analysis
if ($enhancedFile -and $oecdFile) {
    Write-Host "🔬 Generating enhanced + OECD analysis..." -ForegroundColor Yellow
    python enhanced_analysis.py $enhancedFile.FullName $oecdFile.FullName
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Enhanced analysis completed" -ForegroundColor Green
    } else {
        Write-Host "❌ Enhanced analysis failed" -ForegroundColor Red
    }
    Write-Host ""
} elseif ($enhancedFile) {
    Write-Host "⚠️  Skipping enhanced analysis - OECD benchmark not found" -ForegroundColor Yellow
    Write-Host "   Run: echo 2 | python enhanced_pfas_benchmark.py" -ForegroundColor Gray
    Write-Host ""
}

# Move analysis files to review-app
Write-Host "📁 Organizing analysis results..." -ForegroundColor Yellow
Move-Item *_analysis.json review-app\analysis_reports\ -Force -ErrorAction SilentlyContinue
Move-Item timing_*.png review-app\analysis_reports\figures\ -Force -ErrorAction SilentlyContinue
Move-Item timing_*.svg review-app\analysis_reports\figures\ -Force -ErrorAction SilentlyContinue
Move-Item complex_*.png review-app\analysis_reports\figures\ -Force -ErrorAction SilentlyContinue
Move-Item complex_*.svg review-app\analysis_reports\figures\ -Force -ErrorAction SilentlyContinue
Move-Item enhanced_*.png review-app\analysis_reports\figures\ -Force -ErrorAction SilentlyContinue
Move-Item enhanced_*.svg review-app\analysis_reports\figures\ -Force -ErrorAction SilentlyContinue

Write-Host ""
Write-Host "✨ Analysis generation complete!" -ForegroundColor Green
Write-Host ""
Write-Host "📊 Analysis reports saved to: review-app\analysis_reports\" -ForegroundColor Cyan
Write-Host "   View them in the Review App at: http://localhost:5000/analysis" -ForegroundColor Cyan
Write-Host ""

# Check if server is running
$serverProcess = Get-Process -Name node -ErrorAction SilentlyContinue | Where-Object { $_.Path -like "*node*" }
if ($serverProcess) {
    Write-Host "✅ Server appears to be running" -ForegroundColor Green
    Write-Host "   Refresh your browser to see the analysis reports" -ForegroundColor Gray
} else {
    Write-Host "⚠️  Server not running" -ForegroundColor Yellow
    Write-Host "   Start with: cd review-app ; node server.js" -ForegroundColor Gray
}
