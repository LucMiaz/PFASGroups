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

# Run setup
Write-Host "1️⃣ Running setup script..." -ForegroundColor Green
& .\setup-review-app.ps1

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "✅ Setup completed successfully!" -ForegroundColor Green
    Write-Host ""
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
