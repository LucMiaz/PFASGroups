@echo off
REM Enhanced PFAS Benchmark - Complete Pipeline
REM ==========================================
REM Comprehensive PFAS benchmarking suite with all test types
REM Usage: run_enhanced_benchmark.bat

setlocal enabledelayedexpansion

echo.
echo ============================================
echo   ENHANCED PFAS BENCHMARK PIPELINE v2.0
echo ============================================
echo.

REM Check if we're in the right directory
if not exist "enhanced_pfas_benchmark.py" (
    echo ERROR: enhanced_pfas_benchmark.py not found!
    echo Please run this script from the streamlined_benchmark directory
    pause
    exit /b 1
)

set total_errors=0

REM Step 1: Run enhanced functional groups benchmark
echo [1/7] Running Enhanced Functional Groups Benchmark...
echo -------------------------------------------------------
echo Testing comprehensive functional group detection...
echo 1| python enhanced_pfas_benchmark.py
set enhanced_status=!errorlevel!

REM Step 2: Run OECD benchmark
echo [2/7] Running OECD PFAS List Benchmark...
echo -------------------------------------------
echo Validating against official OECD PFAS definitions...
echo 2| python enhanced_pfas_benchmark.py
set oecd_status=!errorlevel!

REM Step 3: Run timing performance benchmark
echo [3/7] Running Timing Performance Benchmark...
echo -----------------------------------------------
echo Testing performance across molecular sizes (200 molecules, 10 iterations)...
echo 3| python enhanced_pfas_benchmark.py
set timing_status=!errorlevel!

REM Step 4: Run non-fluorinated exclusion benchmark
echo [4/7] Running Non-Fluorinated Exclusion Test...
echo -------------------------------------------------
echo Testing specificity by excluding non-PFAS molecules...
echo 4| python enhanced_pfas_benchmark.py
set nonfluor_status=!errorlevel!

REM Step 5: Run complex branched PFAS benchmark
echo [5/7] Running Complex Branched PFAS Benchmark...
echo -------------------------------------------------
echo Testing detection of highly complex branched structures...
echo 5| python enhanced_pfas_benchmark.py
set complex_status=!errorlevel!

REM Step 6: Generate analysis reports
echo [6/7] Generating Comprehensive Analysis Reports...
echo ---------------------------------------------------

REM Find the most recent files
set enhanced_file=
set oecd_file=
set timing_file=
set nonfluor_file=
set complex_file=

for /f "delims=" %%i in ('dir /b /o-d pfas_enhanced_benchmark_*.json 2^>nul') do (
    if not defined enhanced_file set "enhanced_file=%%i"
)

for /f "delims=" %%i in ('dir /b /o-d pfas_oecd_benchmark_*.json 2^>nul') do (
    if not defined oecd_file set "oecd_file=%%i"
)

for /f "delims=" %%i in ('dir /b /o-d pfas_timing_benchmark_*.json 2^>nul') do (
    if not defined timing_file set "timing_file=%%i"
)

for /f "delims=" %%i in ('dir /b /o-d pfas_non_fluorinated_benchmark_*.json 2^>nul') do (
    if not defined nonfluor_file set "nonfluor_file=%%i"
)

for /f "delims=" %%i in ('dir /b /o-d pfas_complex_branched_benchmark_*.json 2^>nul') do (
    if not defined complex_file set "complex_file=%%i"
)

echo Found benchmark files:
if defined enhanced_file echo   ✅ Enhanced: !enhanced_file!
if defined oecd_file echo   ✅ OECD: !oecd_file!
if defined timing_file echo   ✅ Timing: !timing_file!
if defined nonfluor_file echo   ✅ Non-fluorinated: !nonfluor_file!
if defined complex_file echo   ✅ Complex: !complex_file!

set analysis_status=0

REM Run enhanced analysis if both files exist
if defined enhanced_file if defined oecd_file (
    echo   📈 Running enhanced analysis...
    python enhanced_analysis.py "!enhanced_file!" "!oecd_file!"
    if errorlevel 1 set analysis_status=1
)

REM Run timing analysis if file exists
if defined timing_file (
    echo   ⏱️  Running timing analysis...
    python analyze_timing.py "!timing_file!"
    if errorlevel 1 set analysis_status=1
)

REM Run non-fluorinated analysis if file exists
if defined nonfluor_file (
    echo   ❌ Running non-fluorinated analysis...
    python analyze_nonfluorinated.py "!nonfluor_file!"
    if errorlevel 1 set analysis_status=1
)

REM Run complex analysis if file exists
if defined complex_file (
    echo   🌳 Running complex molecules analysis...
    python analyze_complex.py "!complex_file!"
    if errorlevel 1 set analysis_status=1
)

REM Generate unified report
echo   📋 Generating unified HTML report...
python generate_unified_report.py
if errorlevel 1 set analysis_status=1

REM Step 7: Final summary
echo [7/7] Pipeline Summary Report...
echo ==================================

echo 📊 Benchmark Results:
if !enhanced_status! equ 0 (
    echo   ✅ Enhanced functional groups: SUCCESS
) else (
    echo   ❌ Enhanced functional groups: FAILED
    set /a total_errors+=1
)

if !oecd_status! equ 0 (
    echo   ✅ OECD validation: SUCCESS
) else (
    echo   ❌ OECD validation: FAILED
    set /a total_errors+=1
)

if !timing_status! equ 0 (
    echo   ✅ Timing performance: SUCCESS
) else (
    echo   ❌ Timing performance: FAILED
    set /a total_errors+=1
)

if !nonfluor_status! equ 0 (
    echo   ✅ Non-fluorinated exclusion: SUCCESS
) else (
    echo   ❌ Non-fluorinated exclusion: FAILED
    set /a total_errors+=1
)

if !complex_status! equ 0 (
    echo   ✅ Complex branched structures: SUCCESS
) else (
    echo   ❌ Complex branched structures: FAILED
    set /a total_errors+=1
)

if !analysis_status! equ 0 (
    echo   ✅ Analysis reports: SUCCESS
) else (
    echo   ❌ Analysis reports: FAILED
    set /a total_errors+=1
)

echo.
echo 📁 Generated Benchmark Data Files:
if defined enhanced_file echo   • !enhanced_file!
if defined oecd_file echo   • !oecd_file!
if defined timing_file echo   • !timing_file!
if defined nonfluor_file echo   • !nonfluor_file!
if defined complex_file echo   • !complex_file!

echo.
echo 📄 Generated Analysis Reports:
for /f "delims=" %%i in ('dir /b /o-d unified_pfas_benchmark_report_*.html 2^>nul') do (
    echo   🎯 MAIN: %%i (Complete Dashboard)
    goto :found_unified
)
:found_unified
for %%f in (*benchmark*.html *analysis*.html) do (
    echo %%f | findstr /v "unified_pfas_benchmark_report" >nul && echo   • %%f
)

echo.
echo 📊 Generated Visualizations:
for %%f in (*benchmark*.png *analysis*.png) do (
    echo   • %%f
)

echo.
echo ============================================
if !total_errors! equ 0 (
    echo 🎉 SUCCESS: Complete benchmark pipeline executed successfully!
    echo.
    echo ✅ Functional group detection validated
    echo ✅ OECD PFAS list compatibility confirmed
    echo ✅ Performance scaling analyzed
    echo ✅ Specificity verified (non-PFAS exclusion)
    echo ✅ Complex structure handling validated
    echo ✅ Statistical analysis with error bars
    echo ✅ Interactive HTML reports generated
    echo.
    echo 📋 Review the HTML reports for detailed findings!
) else (
    echo ⚠️  WARNING: !total_errors! component(s) failed
    echo.
    echo 📋 Check the error messages above for details
    echo 📋 Partial results may still be available in generated files
)

echo.
pause
if errorlevel 1 (
    echo ERROR: Timing analysis failed!
    pause
    exit /b 1
)

REM Step 3: Run non-fluorinated benchmark
echo [3/6] Running Non-Fluorinated Exclusion Benchmark...
echo -----------------------------------------------------
echo Testing that systems correctly exclude non-PFAS molecules...
echo 4| python enhanced_pfas_benchmark.py
if errorlevel 1 (
    echo ERROR: Non-fluorinated benchmark failed!
    pause
    exit /b 1
)

REM Step 4: Run the comprehensive benchmark (all molecules)
echo [4/6] Running Comprehensive PFAS Benchmark...
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

REM Step 5: Run the enhanced analysis
echo [5/6] Running Enhanced Analysis...
echo ----------------------------------
python enhanced_analysis.py "%LATEST_BENCHMARK%"
if errorlevel 1 (
    echo ERROR: Enhanced analysis failed!
    pause
    exit /b 1
)

REM Step 6: Generate summaries
echo [6/6] Generating Summaries...
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
for /f "delims=" %%i in ('dir /b /o-d pfas_non_fluorinated_benchmark_*.json 2^>nul') do (
    set "LATEST_NONFLUOR=%%i"
    goto :found_nonfluor
)

:found_nonfluor
echo.
echo ============================================
echo   ENHANCED BENCHMARK COMPLETE!
echo ============================================
echo.
echo Generated Files:
echo   📊 Timing Benchmark Data: %LATEST_TIMING%
echo   📈 Timing Analysis Report: %LATEST_TIMING_ANALYSIS%
echo   🚫 Non-Fluorinated Test: %LATEST_NONFLUOR%
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
echo   - Non-Fluorinated Results: %LATEST_NONFLUOR%
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
echo    ✅ Functional group validation included
echo    ✅ Non-fluorinated exclusion testing
echo    ✅ Comprehensive molecule testing
echo    ✅ System specifications included
echo    ✅ Enhanced visualizations generated
echo.
echo Press any key to exit...
pause >nul