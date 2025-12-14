@echo off
REM Quick Unified Report Generator
REM =============================
REM Generates unified HTML report from existing benchmark results
REM Usage: generate_report.bat

echo 🚀 QUICK UNIFIED REPORT GENERATOR
echo ==================================

REM Check if we're in the right directory
if not exist "generate_unified_report.py" (
    echo ERROR: generate_unified_report.py not found!
    echo Please run this script from the streamlined_benchmark directory
    pause
    exit /b 1
)

echo 📊 Generating unified HTML report from existing benchmark data...
python generate_unified_report.py

if errorlevel 1 (
    echo ❌ ERROR: Report generation failed
    echo 📋 Make sure you have benchmark result files (*.json) in this directory
    echo 📋 Run the full benchmark pipeline first if needed: run_enhanced_benchmark.bat
    pause
    exit /b 1
) else (
    REM Find the most recent unified report
    for /f "delims=" %%i in ('dir /b /o-d unified_pfas_benchmark_report_*.html 2^>nul') do (
        set "unified_report=%%i"
        goto :found_report
    )
    
    echo ❌ ERROR: Report file not found after generation
    pause
    exit /b 1
    
    :found_report
    echo.
    echo 🎉 SUCCESS: Unified report generated!
    echo 📄 Report file: %unified_report%
    echo.
    echo 🌐 Open this file in your web browser to view the complete analysis dashboard
    echo 📊 The report includes all available benchmark results in a single view
)

pause