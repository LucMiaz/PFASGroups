@echo off
REM Unified PFAS Benchmark Runner (Windows)
REM Runs all benchmark types and generates unified HTML report

echo 🚀 PFAS BENCHMARK SUITE - UNIFIED RUNNER
echo =========================================
echo.

REM Check if we're in the right directory
if not exist "enhanced_pfas_benchmark.py" (
    echo ❌ Error: enhanced_pfas_benchmark.py not found!
    echo    Please run this script from the benchmark/streamlined_benchmark directory
    pause
    exit /b 1
)

echo 📋 Running all PFAS benchmarks...
echo.

REM Functional Groups Benchmark
echo 1️⃣ Running Functional Groups Benchmark...
echo 1 | python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ❌ Functional Groups Benchmark failed
    pause
    exit /b 1
)
echo ✅ Functional Groups Benchmark completed
echo.

REM OECD Validation Benchmark
echo 2️⃣ Running OECD Validation Benchmark...
echo 2 | python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ❌ OECD Validation Benchmark failed
    pause
    exit /b 1
)
echo ✅ OECD Validation Benchmark completed
echo.

REM Timing Performance Benchmark
echo 3️⃣ Running Timing Performance Benchmark...
echo 3 | python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ❌ Timing Performance Benchmark failed
    pause
    exit /b 1
)
echo ✅ Timing Performance Benchmark completed
echo.

REM Non-Fluorinated Exclusion Benchmark
echo 4️⃣ Running Non-Fluorinated Exclusion Benchmark...
echo 4 | python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ❌ Non-Fluorinated Exclusion Benchmark failed
    pause
    exit /b 1
)
echo ✅ Non-Fluorinated Exclusion Benchmark completed
echo.

REM Complex Branched Structures Benchmark
echo 5️⃣ Running Complex Branched Structures Benchmark...
echo 5 | python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ❌ Complex Branched Structures Benchmark failed
    pause
    exit /b 1
)
echo ✅ Complex Branched Structures Benchmark completed
echo.

REM Generate Unified Report
echo 📊 Generating Unified HTML Report...
python generate_unified_report.py
if errorlevel 1 (
    echo ❌ Unified Report generation failed
    pause
    exit /b 1
)

echo.
echo 🎉 ALL BENCHMARKS COMPLETED SUCCESSFULLY!
echo.
echo 📋 Results Summary:
echo    • Functional Groups: Tests basic PFAS group detection
echo    • OECD Validation: Validates against OECD database
echo    • Timing Performance: Measures execution speed scaling
echo    • Non-Fluorinated: Ensures proper exclusion of non-PFAS
echo    • Complex Branched: Tests complex molecular structures
echo.

REM Find the latest unified report
for /f "delims=" %%i in ('dir /b /od "unified_pfas_benchmark_report_*.html" 2^>nul ^| findstr /r ".*" ^| tail -1') do set "latest_report=%%i"
if defined latest_report (
    echo 📄 Unified Report: %latest_report%
    echo 🌐 Open the HTML file in your browser to view detailed results
) else (
    echo 📄 Unified Report: Check for unified_pfas_benchmark_report_*.html files
)

echo.
echo ✨ Benchmark Suite Complete!
echo.
pause