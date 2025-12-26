#!/usr/bin/env pwsh
$ErrorActionPreference = "Stop"

Write-Host "🚀 Setting up PFAS Benchmark Reviewer Application" -ForegroundColor Cyan
Write-Host "================================================" -ForegroundColor Cyan

# Change to review-app directory
$scriptPath = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location (Join-Path $scriptPath "review-app")

Write-Host "📦 Installing Node.js dependencies..." -ForegroundColor Yellow
npm install
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Failed to install Node.js dependencies" -ForegroundColor Red
    exit 1
}

Write-Host "📦 Installing React app dependencies..." -ForegroundColor Yellow
# Check if the automated React app creation completed
if (Test-Path "client") {
    Write-Host "Using auto-generated React app..." -ForegroundColor Green
    Set-Location client
    
    # Install dependencies from package.json
    npm install
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ Failed to install React dependencies" -ForegroundColor Red
        exit 1
    }
    
    # Copy our custom React components
    if (Test-Path "../client-src/src") {
        Write-Host "📝 Copying custom React components..." -ForegroundColor Yellow
        Copy-Item -Path "../client-src/src/*" -Destination "src/" -Recurse -Force
        Copy-Item -Path "../client-src/public/index.html" -Destination "public/" -Force
        Copy-Item -Path "../client-src/package.json" -Destination "package-temp.json" -Force
        
        # Merge package.json dependencies
        node -e @"
        const fs = require('fs');
        const existing = JSON.parse(fs.readFileSync('package.json', 'utf8'));
        const custom = JSON.parse(fs.readFileSync('package-temp.json', 'utf8'));
        
        existing.dependencies = { ...existing.dependencies, ...custom.dependencies };
        existing.proxy = custom.proxy;
        
        fs.writeFileSync('package.json', JSON.stringify(existing, null, 2));
        fs.unlinkSync('package-temp.json');
        console.log('✓ Package.json updated');
"@
        
        # Install updated dependencies
        npm install
    }
    
    Set-Location ..
} else {
    Write-Host "Using pre-built React components..." -ForegroundColor Green
    # If auto-generation didn't complete, use our pre-built structure
    Move-Item -Path "client-src" -Destination "client" -Force
    Set-Location client
    npm install
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ Failed to install client dependencies" -ForegroundColor Red
        exit 1
    }
    Set-Location ..
}

Write-Host "🗄️  Initializing database..." -ForegroundColor Yellow
node -e @"
const Database = require('./database/database.js');
const db = new Database();
console.log('✓ Database initialized');
"@

Write-Host "📊 Importing existing benchmark data..." -ForegroundColor Yellow
node scripts/import-benchmark-data.js
if ($LASTEXITCODE -ne 0) {
    Write-Host "⚠️  Data import failed - will continue without existing data" -ForegroundColor Yellow
}

Write-Host "🔧 Creating startup scripts..." -ForegroundColor Yellow

# Create development startup script
$devScript = @'
#!/usr/bin/env pwsh
Write-Host "🚀 Starting PFAS Benchmark Reviewer in development mode..." -ForegroundColor Cyan
Write-Host "Server: http://localhost:5000" -ForegroundColor Green
Write-Host "React Dev Server: http://localhost:3000 (if available)" -ForegroundColor Green
Write-Host ""
npm run dev
'@
$devScript | Out-File -FilePath "start-dev.ps1" -Encoding utf8 -Force

# Create production startup script
$prodScript = @'
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
'@
$prodScript | Out-File -FilePath "start-prod.ps1" -Encoding utf8 -Force

# Create data import script
$importScript = @'
#!/usr/bin/env pwsh
Write-Host "📊 Importing latest benchmark data..." -ForegroundColor Cyan
node scripts/import-benchmark-data.js
Write-Host "✅ Data import completed" -ForegroundColor Green
'@
$importScript | Out-File -FilePath "import-latest-data.ps1" -Encoding utf8 -Force

Write-Host ""
Write-Host "✅ Setup completed successfully!" -ForegroundColor Green
Write-Host ""
Write-Host "🎯 Quick Start:" -ForegroundColor Cyan
Write-Host "  Development mode: .\start-dev.ps1" -ForegroundColor Yellow
Write-Host "  Production mode:  .\start-prod.ps1" -ForegroundColor Yellow
Write-Host "  Import data:      .\import-latest-data.ps1" -ForegroundColor Yellow
Write-Host ""
Write-Host "📝 Manual Steps:" -ForegroundColor Cyan
Write-Host "  1. Run .\import-latest-data.ps1 to load existing benchmark data"
Write-Host "  2. Run .\start-dev.ps1 to start the development server"
Write-Host "  3. Open http://localhost:3000 (dev) or http://localhost:5000 (prod)"
Write-Host ""
Write-Host "🔧 Development:" -ForegroundColor Cyan
Write-Host "  - Server runs on port 5000"
Write-Host "  - React dev server on port 3000 (proxies API calls to 5000)"
Write-Host "  - Database: SQLite file in database/pfas_benchmark.db"
Write-Host ""
Write-Host "📊 Features Available:" -ForegroundColor Cyan
Write-Host "  ✓ Molecule visualization with RDKit.js"
Write-Host "  ✓ Pagination and filtering"
Write-Host "  ✓ Manual review interface"
Write-Host "  ✓ Accuracy computation"
Write-Host "  ✓ Data export (JSON/CSV)"
Write-Host "  ✓ Dashboard with statistics"
