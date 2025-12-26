#!/usr/bin/env pwsh
Write-Host "🚀 Building and starting PFAS Benchmark Reviewer..." -ForegroundColor Cyan
Set-Location client
npm run build
if ($LASTEXITCODE -eq 0) {
    Write-Host "✓ React app built" -ForegroundColor Green
}
Set-Location ..
Write-Host "🌐 Starting server at http://localhost:5000" -ForegroundColor Green
npm start
