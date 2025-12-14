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

REM Step 1: Run the enhanced timing benchmark
echo [1/5] Running Enhanced Timing Benchmark...
echo --------------------------------------------
echo This will run 200 molecules with 10 iterations each for statistical analysis...
echo 3| python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ERROR: Enhanced timing benchmark failed!
    pause
    exit /b 1
)

REM Get the most recent timing benchmark file
for /f "delims=" %%i in ('dir /b /o-d pfas_timing_benchmark_*.json 2^>nul') do (
    set "LATEST_TIMING=%%i"
    goto :found_timing
)

echo ERROR: No timing benchmark results file found!
pause
exit /b 1

:found_timing
echo Latest timing benchmark file: %LATEST_TIMING%
echo.

REM Step 2: Run timing analysis
echo [2/5] Running Timing Analysis...
echo ----------------------------------
python analyze_timing.py "%LATEST_TIMING%"
if errorlevel 1 (
    echo ERROR: Timing analysis failed!
    pause
    exit /b 1
)

REM Step 3: Run the comprehensive benchmark (all molecules)
echo [3/5] Running Comprehensive PFAS Benchmark...
echo -----------------------------------------------
echo This will run the full enhanced benchmark with all molecules...
echo 1| python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ERROR: Enhanced benchmark failed!
    pause
    exit /b 1
)

REM Get the most recent comprehensive benchmark file
for /f "delims=" %%i in ('dir /b /o-d pfas_enhanced_benchmark_*.json 2^>nul') do (
    set "LATEST_BENCHMARK=%%i"
    goto :found_benchmark
)

echo ERROR: No benchmark results file found!
pause
exit /b 1

:found_benchmark
echo Latest comprehensive benchmark file: %LATEST_BENCHMARK%
echo.

REM Step 4: Run the enhanced analysis
echo [4/5] Running Enhanced Analysis...
echo ----------------------------------
python enhanced_analysis.py "%LATEST_BENCHMARK%"
if errorlevel 1 (
    echo ERROR: Enhanced analysis failed!
    pause
    exit /b 1
)

REM Step 5: Generate summaries
echo [5/5] Generating Summaries...
echo ------------------------------
python benchmark_summary.py "%LATEST_BENCHMARK%"
if errorlevel 1 (
    echo ERROR: Summary generation failed!
    pause
    exit /b 1
)

python enhanced_summary.py "%LATEST_BENCHMARK%"
if errorlevel 1 (
    echo ERROR: Enhanced summary generation failed!
    pause
    exit /b 1
)

REM Get the most recent analysis files
for /f "delims=" %%i in ('dir /b /o-d enhanced_pfas_analysis_*.html 2^>nul') do (
    set "LATEST_ANALYSIS=%%i"
    goto :found_analysis
)

:found_analysis
for /f "delims=" %%i in ('dir /b /o-d timing_analysis_*.html 2^>nul') do (
    set "LATEST_TIMING_ANALYSIS=%%i"
    goto :found_timing_analysis
)

:found_timing_analysis
echo.
echo ============================================
echo   ENHANCED BENCHMARK COMPLETE!
echo ============================================
echo.
echo Generated Files:
echo   📊 Timing Benchmark Data: %LATEST_TIMING%
echo   📈 Timing Analysis Report: %LATEST_TIMING_ANALYSIS%
echo   🧪 Comprehensive Benchmark Data: %LATEST_BENCHMARK%
echo   📋 Enhanced Analysis Report: %LATEST_ANALYSIS%
echo   🎯 Performance Heatmaps (PNG/SVG)
echo   🔗 Sankey Diagrams (PNG/SVG)
echo   ⏱️  Timing Visualizations (PNG/SVG)
echo   📊 Statistical Error Bar Charts
echo.
echo Main Analysis Reports:
echo   - Timing Analysis: %LATEST_TIMING_ANALYSIS%
echo   - Enhanced Analysis: %LATEST_ANALYSIS%
echo.
echo Opening analysis reports...
echo.

REM Try to open both HTML reports in default browser
start "" "%LATEST_TIMING_ANALYSIS%" 2>nul
if errorlevel 1 (
    echo Note: Could not auto-open timing analysis report.
)

start "" "%LATEST_ANALYSIS%" 2>nul
if errorlevel 1 (
    echo Note: Could not auto-open enhanced analysis report.
)

echo 🎉 Complete pipeline executed successfully!
echo    ✅ Statistical timing analysis with 10x iterations
echo    ✅ Comprehensive molecule testing
echo    ✅ System specifications included
echo    ✅ Enhanced visualizations generated
echo.
echo Press any key to exit...
pause >nul