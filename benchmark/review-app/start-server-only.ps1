# Start only the backend server (useful for testing API or if React app is already built)
Write-Host "🚀 Starting PFAS Benchmark Reviewer server..." -ForegroundColor Green
Write-Host "Server: http://localhost:5000" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

Set-Location $PSScriptRoot
node server.js
