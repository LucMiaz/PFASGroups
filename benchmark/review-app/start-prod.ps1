# Start PFAS Benchmark Reviewer in Production Mode
Write-Host "🚀 Building and starting PFAS Benchmark Reviewer..." -ForegroundColor Green

# Build the React app
Write-Host "Building React app..." -ForegroundColor Cyan
Set-Location "$PSScriptRoot\client"
npm run build

if ($LASTEXITCODE -eq 0) {
    Write-Host "✓ React app built successfully" -ForegroundColor Green
    Write-Host ""
    Write-Host "🌐 Starting server at http://localhost:5000" -ForegroundColor Cyan
    Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
    Write-Host ""
    
    Set-Location "$PSScriptRoot"
    node server.js
} else {
    Write-Host "❌ Build failed" -ForegroundColor Red
    exit 1
}
