#!/bin/bash
# Enhanced PFAS Benchmark - Complete Pipeline
# ==========================================
# Comprehensive PFAS benchmarking suite with all test types
# Usage: ./run_enhanced_benchmark.sh

echo ""
echo "============================================"
echo "  ENHANCED PFAS BENCHMARK PIPELINE v2.0"
echo "============================================"
echo ""

# Check if we're in the right directory
if [ ! -f "enhanced_pfas_benchmark.py" ]; then
    echo "ERROR: enhanced_pfas_benchmark.py not found!"
    echo "Please run this script from the streamlined_benchmark directory"
    exit 1
fi

# Step 1: Run enhanced functional groups benchmark
echo "[1/7] Running Enhanced Functional Groups Benchmark..."
echo "-------------------------------------------------------"
echo "Testing comprehensive functional group detection..."
echo "1" | python enhanced_pfas_benchmark.py
enhanced_status=$?

# Step 2: Run OECD benchmark
echo "[2/7] Running OECD PFAS List Benchmark..."
echo "-------------------------------------------"
echo "Validating against official OECD PFAS definitions..."
echo "2" | python enhanced_pfas_benchmark.py
oecd_status=$?

# Step 3: Run timing performance benchmark
echo "[3/7] Running Timing Performance Benchmark..."
echo "-----------------------------------------------"
echo "Testing performance across molecular sizes (200 molecules, 10 iterations)..."
echo "3" | python enhanced_pfas_benchmark.py
timing_status=$?

# Step 4: Run non-fluorinated exclusion benchmark
echo "[4/7] Running Non-Fluorinated Exclusion Test..."
echo "-------------------------------------------------"
echo "Testing specificity by excluding non-PFAS molecules..."
echo "4" | python enhanced_pfas_benchmark.py
nonfluor_status=$?

# Step 5: Run complex branched PFAS benchmark
echo "[5/7] Running Complex Branched PFAS Benchmark..."
echo "-------------------------------------------------"
echo "Testing detection of highly complex branched structures..."
echo "5" | python enhanced_pfas_benchmark.py
complex_status=$?

# Step 6: Generate analysis reports
echo "[6/7] Generating Individual Analysis Reports..."
echo "---------------------------------------------------"

# Find the most recent files
enhanced_file=$(ls -t pfas_enhanced_benchmark_*.json 2>/dev/null | head -n1)
oecd_file=$(ls -t pfas_oecd_benchmark_*.json 2>/dev/null | head -n1)
timing_file=$(ls -t pfas_timing_benchmark_*.json 2>/dev/null | head -n1)
nonfluor_file=$(ls -t pfas_non_fluorinated_benchmark_*.json 2>/dev/null | head -n1)
complex_file=$(ls -t pfas_complex_branched_benchmark_*.json 2>/dev/null | head -n1)

analysis_status=0

echo "Found benchmark files:"
if [[ -n "$enhanced_file" ]]; then echo "  ✅ Enhanced: $enhanced_file"; fi
if [[ -n "$oecd_file" ]]; then echo "  ✅ OECD: $oecd_file"; fi
if [[ -n "$timing_file" ]]; then echo "  ✅ Timing: $timing_file"; fi
if [[ -n "$nonfluor_file" ]]; then echo "  ✅ Non-fluorinated: $nonfluor_file"; fi
if [[ -n "$complex_file" ]]; then echo "  ✅ Complex: $complex_file"; fi

# Run enhanced analysis if both files exist
if [[ -n "$enhanced_file" && -n "$oecd_file" ]]; then
    echo "  📈 Running enhanced analysis..."
    python enhanced_analysis.py "$enhanced_file" "$oecd_file"
    if [ $? -ne 0 ]; then analysis_status=1; fi
fi

# Run timing analysis if file exists
if [[ -n "$timing_file" ]]; then
    echo "  ⏱️  Running timing analysis..."
    python analyze_timing.py "$timing_file"
    if [ $? -ne 0 ]; then analysis_status=1; fi
fi

# Run non-fluorinated analysis if file exists
if [[ -n "$nonfluor_file" ]]; then
    echo "  ❌ Running non-fluorinated analysis..."
    python analyze_nonfluorinated.py "$nonfluor_file"
    if [ $? -ne 0 ]; then analysis_status=1; fi
fi

# Run complex analysis if file exists
if [[ -n "$complex_file" ]]; then
    echo "  🌳 Running complex molecules analysis..."
    python analyze_complex.py "$complex_file"
    if [ $? -ne 0 ]; then analysis_status=1; fi
fi

# Generate unified report
echo "  📋 Generating unified HTML report..."
python generate_unified_report.py
if [ $? -eq 0 ]; then
    unified_report=$(ls -t unified_pfas_benchmark_report_*.html 2>/dev/null | head -n1)
    if [[ -n "$unified_report" ]]; then
        echo "  ✅ Unified report: $unified_report"
    fi
else
    analysis_status=1
fi

# Step 7: Final summary
echo "[7/7] Pipeline Summary Report..."
echo "=================================="

total_errors=0

echo "📊 Benchmark Results:"
if [ $enhanced_status -eq 0 ]; then
    echo "  ✅ Enhanced functional groups: SUCCESS"
else
    echo "  ❌ Enhanced functional groups: FAILED"
    ((total_errors++))
fi

if [ $oecd_status -eq 0 ]; then
    echo "  ✅ OECD validation: SUCCESS"
else
    echo "  ❌ OECD validation: FAILED"
    ((total_errors++))
fi

if [ $timing_status -eq 0 ]; then
    echo "  ✅ Timing performance: SUCCESS"
else
    echo "  ❌ Timing performance: FAILED"
    ((total_errors++))
fi

if [ $nonfluor_status -eq 0 ]; then
    echo "  ✅ Non-fluorinated exclusion: SUCCESS"
else
    echo "  ❌ Non-fluorinated exclusion: FAILED"
    ((total_errors++))
fi

if [ $complex_status -eq 0 ]; then
    echo "  ✅ Complex branched structures: SUCCESS"
else
    echo "  ❌ Complex branched structures: FAILED"
    ((total_errors++))
fi

if [ $analysis_status -eq 0 ]; then
    echo "  ✅ Analysis reports: SUCCESS"
else
    echo "  ❌ Analysis reports: FAILED"
    ((total_errors++))
fi

echo ""
echo "📁 Generated Benchmark Data Files:"
if [[ -n "$enhanced_file" ]]; then echo "  • $enhanced_file"; fi
if [[ -n "$oecd_file" ]]; then echo "  • $oecd_file"; fi
if [[ -n "$timing_file" ]]; then echo "  • $timing_file"; fi
if [[ -n "$nonfluor_file" ]]; then echo "  • $nonfluor_file"; fi
if [[ -n "$complex_file" ]]; then echo "  • $complex_file"; fi

# List HTML reports  
unified_report=$(ls -t unified_pfas_benchmark_report_*.html 2>/dev/null | head -n1)
other_reports=$(ls -t *benchmark*.html *analysis*.html 2>/dev/null | grep -v "unified_pfas_benchmark_report" | head -5)

echo ""
echo "📄 Generated Analysis Reports:"
if [[ -n "$unified_report" ]]; then
    echo "  🎯 MAIN: $unified_report (Complete Dashboard)"
fi
if [[ -n "$other_reports" ]]; then
    echo "$other_reports" | while read report; do
        echo "  • $report"
    done
fi

# List plot files
plot_files=$(ls -t *benchmark*.png *analysis*.png 2>/dev/null | head -5)
if [[ -n "$plot_files" ]]; then
    echo ""
    echo "📊 Generated Visualizations:"
    echo "$plot_files" | while read plot; do
        echo "  • $plot"
    done
fi

echo ""
echo "============================================"
if [ $total_errors -eq 0 ]; then
    echo "🎉 SUCCESS: Complete benchmark pipeline executed successfully!"
    echo ""
    echo "✅ Functional group detection validated"
    echo "✅ OECD PFAS list compatibility confirmed"
    echo "✅ Performance scaling analyzed" 
    echo "✅ Specificity verified (non-PFAS exclusion)"
    echo "✅ Complex structure handling validated"
    echo "✅ Statistical analysis with error bars"
    echo "✅ Interactive HTML reports generated"
    echo ""
    echo "📋 Review the HTML reports for detailed findings!"
    exit 0
else
    echo "⚠️  WARNING: $total_errors component(s) failed"
    echo ""
    echo "📋 Check the error messages above for details"
    echo "📋 Partial results may still be available in generated files"
    exit 1
fi