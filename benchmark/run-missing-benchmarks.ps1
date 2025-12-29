# Run OECD, Non-Fluorinated, and Complex Branched benchmarks
# Navigate to benchmark directory
Set-Location c:\Users\luc\git\PFASGroups\benchmark

Write-Host "🚀 Running PFAS Benchmarks" -ForegroundColor Cyan
Write-Host ("=" * 50)

# Run OECD Benchmark (option 2)
Write-Host "`n2️⃣  Running OECD Benchmark..." -ForegroundColor Yellow
"2" | python enhanced_pfas_benchmark.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ OECD Benchmark failed" -ForegroundColor Red
    exit 1
}
Write-Host "✅ OECD Benchmark completed" -ForegroundColor Green

# Run Non-Fluorinated Benchmark (option 4)
Write-Host "`n4️⃣  Running Non-Fluorinated Benchmark (50+ molecules)..." -ForegroundColor Yellow
"4" | python enhanced_pfas_benchmark.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Non-Fluorinated Benchmark failed" -ForegroundColor Red
    exit 1
}
Write-Host "✅ Non-Fluorinated Benchmark completed" -ForegroundColor Green

# Run Complex Branched Benchmark (option 5)
Write-Host "`n5️⃣  Running Complex Branched Benchmark (50+ molecules)..." -ForegroundColor Yellow
"5" | python enhanced_pfas_benchmark.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Complex Branched Benchmark failed" -ForegroundColor Red
    exit 1
}
Write-Host "✅ Complex Branched Benchmark completed" -ForegroundColor Green

# Move files to data directory
Write-Host "`n📁 Moving files to data directory..." -ForegroundColor Yellow
Move-Item pfas_*_benchmark_*.json data\ -Force -ErrorAction SilentlyContinue

Write-Host "`n🎉 All benchmarks completed!" -ForegroundColor Green
Write-Host "`nNext step: Run import script to update the database" -ForegroundColor Cyan
Write-Host "  cd review-app" -ForegroundColor Gray
Write-Host "  node scripts/import-benchmark-data.js" -ForegroundColor Gray
