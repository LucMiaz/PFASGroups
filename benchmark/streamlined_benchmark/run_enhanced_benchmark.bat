@echo off
REM Enhanced PFAS Benchmark - Single Command Execution
REM ===================================================
REM This batch file runs the complete enhanced PFAS benchmarking pipeline
REM Usage: run_enhanced_benchmark.bat

echo.
echo ============================================
echo   ENHANCED PFAS BENCHMARK PIPELINE
echo ============================================
echo.

REM Check if we're in the right directory
if not exist "enhanced_pfas_benchmark.py" (
    echo ERROR: enhanced_pfas_benchmark.py not found!
    echo Please run this script from the streamlined_benchmark directory
    pause
    exit /b 1
)

REM Step 1: Run the enhanced benchmark
echo [1/3] Running Enhanced PFAS Benchmark...
echo ----------------------------------------
python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ERROR: Enhanced benchmark failed!
    pause
    exit /b 1
)

REM Get the most recent benchmark file
for /f "delims=" %%i in ('dir /b /o-d pfas_enhanced_benchmark_*.json 2^>nul') do (
    set "LATEST_BENCHMARK=%%i"
    goto :found_benchmark
)

echo ERROR: No benchmark results file found!
pause
exit /b 1

:found_benchmark
echo Latest benchmark file: %LATEST_BENCHMARK%
echo.

REM Step 2: Run the enhanced analysis
echo [2/3] Running Enhanced Analysis...
echo ----------------------------------
python enhanced_analysis.py "%LATEST_BENCHMARK%"
if errorlevel 1 (
    echo ERROR: Enhanced analysis failed!
    pause
    exit /b 1
)

REM Step 3: Show summary
echo [3/3] Generating Summary...
echo ----------------------------
python benchmark_summary.py "%LATEST_BENCHMARK%"
if errorlevel 1 (
    echo ERROR: Summary generation failed!
    pause
    exit /b 1
)

REM Get the most recent analysis file
for /f "delims=" %%i in ('dir /b /o-d enhanced_pfas_analysis_*.html 2^>nul') do (
    set "LATEST_ANALYSIS=%%i"
    goto :found_analysis
)

:found_analysis
echo.
echo ============================================
echo   ENHANCED BENCHMARK COMPLETE!
echo ============================================
echo.
echo Generated Files:
echo   - Benchmark Data: %LATEST_BENCHMARK%
echo   - Analysis Report: %LATEST_ANALYSIS%
echo   - Performance Heatmaps (PNG/SVG)
echo   - Sankey Diagrams (PNG/SVG)
echo.
echo Opening analysis report...

REM Try to open the HTML report in default browser
start "" "%LATEST_ANALYSIS%" 2>nul
if errorlevel 1 (
    echo Note: Could not auto-open report. Please open %LATEST_ANALYSIS% manually.
)

echo.
echo Press any key to exit...
pause >nul