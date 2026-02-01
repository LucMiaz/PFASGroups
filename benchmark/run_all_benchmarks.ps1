# Unified PFAS Benchmark Runner (Windows PowerShell)
# Runs all benchmark types and generates unified HTML report
# Now includes database integration for the review app

Write-Host "🚀 PFAS BENCHMARK SUITE - UNIFIED RUNNER" -ForegroundColor Cyan
Write-Host ("=" * 50) -ForegroundColor Cyan
Write-Host ""

# Check if we're in the right directory
if (-not (Test-Path "scripts/enhanced_pfas_benchmark.py")) {
    Write-Host "❌ Error: scripts/enhanced_pfas_benchmark.py not found!" -ForegroundColor Red
    Write-Host "   Please run this script from the benchmark directory" -ForegroundColor Red
    exit 1
}

# Create output directories
New-Item -ItemType Directory -Force -Path scripts, reports, data, imgs, html | Out-Null

Write-Host "📋 Running all PFAS benchmarks..." -ForegroundColor Yellow
Write-Host ""

# Function to run benchmark
function Run-Benchmark {
    param($number, $name)
    Write-Host "$number Running $name..." -ForegroundColor Yellow
    $input = "$($number.Substring(0,1))"
    $input | python scripts/enhanced_pfas_benchmark.py
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ $name failed" -ForegroundColor Red
        exit 1
    }
    Write-Host "✅ $name completed" -ForegroundColor Green
    Write-Host ""
}

# Run all benchmarks
Run-Benchmark "1️⃣" "Functional Groups Benchmark"
Run-Benchmark "2️⃣" "OECD Validation Benchmark"
Run-Benchmark "3️⃣" "Timing Performance Benchmark"
Run-Benchmark "4️⃣" "Non-Fluorinated Exclusion Benchmark"
Run-Benchmark "5️⃣" "Complex Branched Structures Benchmark"

# Highly Branched Compounds Test
Write-Host "6️⃣ Running Highly Branched Compounds Test..." -ForegroundColor Yellow
python scripts/test_highly_branched.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Highly Branched Compounds Test failed" -ForegroundColor Red
    exit 1
}
Write-Host "✅ Highly Branched Compounds Test completed" -ForegroundColor Green
Write-Host ""

# Telomer Validation
Write-Host "7️⃣ Running Telomer Detection Validation..." -ForegroundColor Yellow
python scripts/validate_telomers.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Telomer Validation failed" -ForegroundColor Red
    exit 1
}
Write-Host "✅ Telomer Validation completed" -ForegroundColor Green
Write-Host ""

# Generate Unified Report
Write-Host "📊 Generating Unified HTML Report..." -ForegroundColor Yellow
python scripts/generate_unified_report.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Unified Report generation failed" -ForegroundColor Red
    exit 1
}
Write-Host "✅ Unified Report generated" -ForegroundColor Green
Write-Host ""

# Run analysis scripts BEFORE organizing files
Write-Host "📊 Running Analysis Scripts..." -ForegroundColor Cyan
Write-Host ""

# Timing Analysis
Write-Host "⏱️  Analyzing timing performance..." -ForegroundColor Yellow
$timingFile = Get-ChildItem "data\pfas_timing_benchmark_*.json" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($timingFile) {
    python scripts\analyze_timing.py $timingFile.FullName
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Timing analysis completed" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Timing analysis failed" -ForegroundColor Yellow
    }
} else {
    Write-Host "⚠️  No timing benchmark file found" -ForegroundColor Yellow
}
Write-Host ""

# Complex Branched Analysis
Write-Host "🧬 Analyzing complex branched structures..." -ForegroundColor Yellow
$complexFile = Get-ChildItem "data\pfas_complex_branched_benchmark_*.json" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($complexFile) {
    python scripts\analyze_complex.py $complexFile.FullName
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Complex branched analysis completed" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Complex branched analysis failed" -ForegroundColor Yellow
    }
} else {
    Write-Host "⚠️  No complex branched benchmark file found" -ForegroundColor Yellow
}
Write-Host ""

# Enhanced Analysis
Write-Host "🔬 Analyzing enhanced functional groups and OECD..." -ForegroundColor Yellow
$enhancedFile = Get-ChildItem "data\pfas_enhanced_benchmark_*.json" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
$oecdFile = Get-ChildItem "data\pfas_oecd_benchmark_*.json" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if ($enhancedFile -and $oecdFile) {
    python scripts\enhanced_analysis.py $enhancedFile.FullName $oecdFile.FullName
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Enhanced analysis completed" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Enhanced analysis failed" -ForegroundColor Yellow
    }
} else {
    Write-Host "⚠️  Missing enhanced or OECD benchmark files" -ForegroundColor Yellow
    if (-not $enhancedFile) { Write-Host "    Missing: enhanced benchmark" -ForegroundColor Gray }
    if (-not $oecdFile) { Write-Host "    Missing: OECD benchmark" -ForegroundColor Gray }
}
Write-Host ""

# Organize output files
Write-Host "📁 Organizing output files..." -ForegroundColor Yellow
# Keep data files in data/ directory (they're already there)
# Move HTML files
Move-Item *.html reports\ -Force -ErrorAction SilentlyContinue
Move-Item unified_pfas_benchmark_report_*.html reports\ -Force -ErrorAction SilentlyContinue
# Move image files
Move-Item *.png imgs\ -Force -ErrorAction SilentlyContinue
Move-Item *.svg imgs\ -Force -ErrorAction SilentlyContinue
# Move any analysis JSON files to data/
Move-Item *_analysis.json data\ -Force -ErrorAction SilentlyContinue
Move-Item *_results.json data\ -Force -ErrorAction SilentlyContinue
Move-Item *_benchmark_*.json data\ -Force -ErrorAction SilentlyContinue
Write-Host ""
Write-Host "🎉 ALL BENCHMARKS COMPLETED SUCCESSFULLY!" -ForegroundColor Green
Write-Host ""
Write-Host "📋 Results Summary:" -ForegroundColor Cyan
Write-Host "   • Functional Groups: Tests basic PFAS group detection"
Write-Host "   • OECD Validation: Validates against OECD database"
Write-Host "   • Timing Performance: Measures execution speed scaling"
Write-Host "   • Non-Fluorinated: Ensures proper exclusion of non-PFAS"
Write-Host "   • Complex Branched: Tests complex molecular structures"
Write-Host "   • Highly Branched: Tests functional groups on perfluorinated components"
Write-Host "   • Telomer Validation: Tests detection of fluorotelomers on PubChem dataset"
Write-Host "   • Comprehensive Statistics: LaTeX tables and benchmark summary (reports/benchmark_summary.json)"
Write-Host ""

# Telomer Validation Report
Write-Host "🧪 Generating telomer validation report..." -ForegroundColor Yellow
if (Test-Path "data\telomer_validation_results.json") {
    python scripts\generate_telomer_report.py
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Telomer validation report completed" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Telomer validation report failed" -ForegroundColor Yellow
    }
} else {
    Write-Host "⚠️  No telomer validation data found" -ForegroundColor Yellow
}
Write-Host ""

# Comprehensive Benchmark Analysis
Write-Host "📊 Generating comprehensive benchmark statistics..." -ForegroundColor Yellow
python analyze_benchmarks_simple.py
if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Comprehensive benchmark analysis completed" -ForegroundColor Green
} else {
    Write-Host "⚠️  Comprehensive benchmark analysis failed" -ForegroundColor Yellow
}
Write-Host ""

# Organize analysis results
Write-Host "📁 Organizing analysis results..." -ForegroundColor Yellow

Write-Host "✅ Analysis results organized" -ForegroundColor Green
Write-Host ""

# Import data to database
Write-Host "💾 Importing data to Review App database..." -ForegroundColor Cyan
Push-Location review-app

# Check if database exists and backup if needed
if (Test-Path "database\pfas_benchmark.db") {
    Write-Host "📦 Backing up existing database..." -ForegroundColor Yellow
    $backupName = "database\pfas_benchmark.db.backup_$(Get-Date -Format 'yyyyMMdd_HHmmss')"
    Copy-Item "database\pfas_benchmark.db" $backupName
    Write-Host "✅ Backup created: $backupName" -ForegroundColor Green
    
    # Clear existing data
    Write-Host "🗑️  Clearing old data from database..." -ForegroundColor Yellow
    node -e "const db = require('./database/database'); const d = new db(); d.waitForReady().then(() => { return Promise.all([d.run('DELETE FROM manual_reviews'), d.run('DELETE FROM pfasgroups_results'), d.run('DELETE FROM pfasgroups_results_bycomponent'), d.run('DELETE FROM atlas_results'), d.run('DELETE FROM molecules')]); }).then(() => { console.log('✅ Database cleared'); process.exit(0); });"
}

# Import benchmark data
Write-Host "📥 Importing benchmark data..." -ForegroundColor Yellow
node scripts\import-benchmark-data.js
if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Data imported successfully" -ForegroundColor Green
    
    # Calculate molecular formulas
    Write-Host "🧪 Calculating molecular formulas..." -ForegroundColor Yellow
    python scripts\calculate-formulas.py
    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Molecular formulas calculated" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Formula calculation failed (non-critical)" -ForegroundColor Yellow
    }
} else {
    Write-Host "❌ Data import failed" -ForegroundColor Red
    Pop-Location
    exit 1
}

Pop-Location
Write-Host ""

Write-Host "📊 Database Update Complete!" -ForegroundColor Green
Write-Host "   • Database: review-app\database\pfas_benchmark.db"
Write-Host "   • Analysis reports: reports"
Write-Host ""

$reportFile = Get-ChildItem "html\unified_pfas_benchmark_report_*.html" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending | Select-Object -First 1
}
Write-Host "🌐 Open the HTML file in your browser to view detailed results"
Write-Host "📁 Files organized in: data\ (JSON), html\ (reports), imgs\ (plots)"
Write-Host ""
Write-Host "🔬 Review App:" -ForegroundColor Cyan
Write-Host "   • Navigate to: cd review-app" -ForegroundColor Gray
Write-Host "   • Start server: node server.js" -ForegroundColor Gray
Write-Host "   • Open browser: http://localhost:5000" -ForegroundColor Gray
Write-Host "   • View Analysis Reports tab for timing, complex, and enhanced analysis" -ForegroundColor Gray
Write-Host ""
Write-Host "✨ Benchmark Suite Complete!" -ForegroundColor Green

