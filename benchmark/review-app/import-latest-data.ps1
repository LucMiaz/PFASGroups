#!/usr/bin/env pwsh
Write-Host "📊 Importing latest benchmark data..." -ForegroundColor Cyan
node scripts/import-benchmark-data.js
Write-Host "✅ Data import completed" -ForegroundColor Green
