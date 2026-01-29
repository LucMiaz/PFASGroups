#!/usr/bin/env pwsh

Write-Host "🚀 Quick Start: PFAS Benchmark Reviewer" -ForegroundColor Cyan
Write-Host "======================================" -ForegroundColor Cyan

# Check if we're in the right directory
if (-not (Test-Path "setup-review-app.ps1")) {
    Write-Host "❌ Error: Please run this from the benchmark directory" -ForegroundColor Red
    Write-Host "   cd C:\Users\luc\git\PFASGroups\benchmark" -ForegroundColor Yellow
    exit 1
}

Write-Host ""
Write-Host "📋 Setup checklist:"
Write-Host "  □ Install Node.js dependencies"
Write-Host "  □ Setup SQLite database"  
Write-Host "  □ Import existing benchmark data"
Write-Host "  □ Start the review application"
Write-Host ""

# Check if benchmark data exists
Write-Host "🔍 Checking for benchmark data..." -ForegroundColor Yellow
$dataFiles = @(Get-ChildItem -Path "data" -Filter "pfas*.json" -ErrorAction SilentlyContinue)
if ($dataFiles.Count -eq 0) {
    Write-Host ""
    Write-Host "⚠️  No benchmark data files found!" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Would you like to generate benchmark data now? (Recommended)" -ForegroundColor Cyan
    Write-Host "  This will run: python enhanced_pfas_benchmark.py" -ForegroundColor White
    Write-Host "  Estimated time: 5-15 minutes" -ForegroundColor White
    Write-Host ""
    $response = Read-Host "Generate data now? (Y/n)"
    
    if ($response -eq 'n' -or $response -eq 'N') {
        Write-Host ""
        Write-Host "⚠️  Continuing without data - the review app will be empty." -ForegroundColor Yellow
        Write-Host "   To generate data later, run:" -ForegroundColor White
        Write-Host "   python enhanced_pfas_benchmark.py" -ForegroundColor Gray
        Write-Host ""
    } else {
        Write-Host ""
        Write-Host "📊 Generating OECD benchmark data..." -ForegroundColor Green
        Write-Host "   This may take 5-15 minutes. Please wait..." -ForegroundColor Yellow
        Write-Host ""
        
        # Activate conda environment and run benchmark
        $pythonCmd = "python enhanced_pfas_benchmark.py"
        
        # Check if conda is available
        $condaAvailable = Get-Command conda -ErrorAction SilentlyContinue
        if ($condaAvailable) {
            Write-Host "Using conda environment 'chem'..." -ForegroundColor Cyan
            & conda run -n chem python enhanced_pfas_benchmark.py
        } else {
            Write-Host "Running with default Python..." -ForegroundColor Cyan
            & python enhanced_pfas_benchmark.py
        }
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host ""
            Write-Host "✅ Benchmark data generation completed!" -ForegroundColor Green
            
            # Re-check for data files
            $dataFiles = @(Get-ChildItem -Path "data" -Filter "pfas*.json" -ErrorAction SilentlyContinue)
            Write-Host "✓ Generated $($dataFiles.Count) benchmark data file(s)" -ForegroundColor Green
        } else {
            Write-Host ""
            Write-Host "❌ Data generation failed. Continuing with setup anyway..." -ForegroundColor Red
            Write-Host "   You can generate data manually later." -ForegroundColor Yellow
        }
    }
} else {
    Write-Host "✓ Found $($dataFiles.Count) benchmark data file(s)" -ForegroundColor Green
}

# Run setup
Write-Host ""
Write-Host "1️⃣ Running setup script..." -ForegroundColor Green
& .\setup-review-app.ps1

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "✅ Setup completed successfully!" -ForegroundColor Green
    
    # Import benchmark data if any files exist
    $dataFiles = @(Get-ChildItem -Path "data" -Filter "pfas*.json" -ErrorAction SilentlyContinue)
    if ($dataFiles.Count -gt 0) {
        Write-Host ""
        Write-Host "2️⃣ Importing benchmark data..." -ForegroundColor Green
        Push-Location review-app
        node scripts\import-benchmark-data.js
        $importExitCode = $LASTEXITCODE
        Pop-Location
        
        if ($importExitCode -eq 0) {
            Write-Host ""
            Write-Host "✅ Data import completed successfully!" -ForegroundColor Green
            
            # Generate timing report
            Write-Host ""
            Write-Host "3️⃣ Generating timing analysis report..." -ForegroundColor Green
            
            # Check if conda is available
            $condaAvailable = Get-Command conda -ErrorAction SilentlyContinue
            if ($condaAvailable) {
                & conda run -n chem python generate_timing_report.py benchmark_timing_report.html
            } else {
                & python generate_timing_report.py benchmark_timing_report.html
            }
            
            if ($LASTEXITCODE -eq 0) {
                Write-Host "✅ Timing report generated: benchmark_timing_report.html" -ForegroundColor Green
                Write-Host "   Open this file in your browser to view performance metrics" -ForegroundColor Cyan
            } else {
                Write-Host "⚠️  Could not generate timing report (non-critical)" -ForegroundColor Yellow
            }
        } else {
            Write-Host ""
            Write-Host "⚠️  Data import had issues. Check the output above." -ForegroundColor Yellow
        }
    } else {
        Write-Host ""
        Write-Host "⚠️  No data files to import (data directory is empty)." -ForegroundColor Yellow
        Write-Host "   The review app will be empty without data." -ForegroundColor Yellow
    }
    
    Write-Host "🎯 What's next?" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Option 1: Development mode (recommended for reviewing)" -ForegroundColor Yellow
    Write-Host "  cd review-app"
    Write-Host "  .\start-dev.ps1"
    Write-Host "  Open http://localhost:3000"
    Write-Host ""
    Write-Host "Option 2: Production mode" -ForegroundColor Yellow
    Write-Host "  cd review-app"  
    Write-Host "  .\start-prod.ps1"
    Write-Host "  Open http://localhost:5000"
    Write-Host ""
    Write-Host "📊 Features available:" -ForegroundColor Cyan
    Write-Host "  ✓ Prioritized display of misclassified molecules"
    Write-Host "  ✓ Manual review with click buttons (✅❌🤷)"
    Write-Host "  ✓ Enhanced performance metrics"
    Write-Host "  ✓ Real-time accuracy computation"
    Write-Host "  ✓ Export capabilities (JSON/CSV)"
    Write-Host "  ✓ RDKit.js molecular visualization"
    Write-Host ""
    Write-Host "🔬 Pro tip: Enable 'Prioritize Misclassified' filter to review" -ForegroundColor Magenta
    Write-Host "   the most important molecules first!"
    Write-Host ""
} else {
    Write-Host ""
    Write-Host "❌ Setup failed. Please check the error messages above." -ForegroundColor Red
    Write-Host "   You may need to install Node.js or resolve dependency issues." -ForegroundColor Yellow
    exit 1
}
