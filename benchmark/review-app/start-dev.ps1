# Start PFAS Benchmark Reviewer in Development Mode
Write-Host "🚀 Starting PFAS Benchmark Reviewer in development mode..." -ForegroundColor Green
Write-Host "Server: http://localhost:5000" -ForegroundColor Cyan
Write-Host "React Dev Server: http://localhost:3000" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop both servers" -ForegroundColor Yellow
Write-Host ""

# Start the backend server in a new PowerShell window
$serverProcess = Start-Process pwsh -ArgumentList "-NoExit", "-Command", "Set-Location '$PSScriptRoot'; node server.js" -PassThru

# Give the server a moment to start
Start-Sleep -Seconds 2

# Start the React dev server
Set-Location "$PSScriptRoot\client"
npm start

# When React dev server stops, kill the backend server too
Stop-Process -Id $serverProcess.Id -ErrorAction SilentlyContinue
