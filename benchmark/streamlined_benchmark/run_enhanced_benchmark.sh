#!/bin/bash
# Enhanced PFAS Benchmark - Single Command Execution
# ===================================================
# This shell script runs the complete enhanced PFAS benchmarking pipeline
# Usage: ./run_enhanced_benchmark.sh

echo ""
echo "============================================"
echo "  ENHANCED PFAS BENCHMARK PIPELINE"
echo "============================================"
echo ""

# Check if we're in the right directory
if [ ! -f "enhanced_pfas_benchmark.py" ]; then
    echo "ERROR: enhanced_pfas_benchmark.py not found!"
    echo "Please run this script from the streamlined_benchmark directory"
    exit 1
fi

# Step 1: Run the enhanced timing benchmark
echo "[1/5] Running Enhanced Timing Benchmark..."
echo "--------------------------------------------"
echo "This will run 200 molecules with 10 iterations each for statistical analysis..."
echo "3" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "ERROR: Enhanced timing benchmark failed!"
    exit 1
fi

# Get the most recent timing benchmark file
LATEST_TIMING=$(ls -t pfas_timing_benchmark_*.json 2>/dev/null | head -n1)

if [ -z "$LATEST_TIMING" ]; then
    echo "ERROR: No timing benchmark results file found!"
    exit 1
fi

echo "Latest timing benchmark file: $LATEST_TIMING"
echo ""

# Step 2: Run timing analysis
echo "[2/5] Running Timing Analysis..."
echo "----------------------------------"
python analyze_timing.py "$LATEST_TIMING"
if [ $? -ne 0 ]; then
    echo "ERROR: Timing analysis failed!"
    exit 1
fi

# Step 3: Run the comprehensive benchmark (all molecules)
echo "[3/5] Running Comprehensive PFAS Benchmark..."
echo "-----------------------------------------------"
echo "This will run the full enhanced benchmark with all molecules..."
echo "1" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "ERROR: Enhanced benchmark failed!"
    exit 1
fi

# Get the most recent comprehensive benchmark file
LATEST_BENCHMARK=$(ls -t pfas_enhanced_benchmark_*.json 2>/dev/null | head -n1)

if [ -z "$LATEST_BENCHMARK" ]; then
    echo "ERROR: No benchmark results file found!"
    exit 1
fi

echo "Latest comprehensive benchmark file: $LATEST_BENCHMARK"
echo ""

# Step 4: Run the enhanced analysis
echo "[4/5] Running Enhanced Analysis..."
echo "----------------------------------"
python enhanced_analysis.py "$LATEST_BENCHMARK"
if [ $? -ne 0 ]; then
    echo "ERROR: Enhanced analysis failed!"
    exit 1
fi

# Step 5: Generate summaries
echo "[5/5] Generating Summaries..."
echo "------------------------------"
python benchmark_summary.py "$LATEST_BENCHMARK"
if [ $? -ne 0 ]; then
    echo "ERROR: Summary generation failed!"
    exit 1
fi

python enhanced_summary.py "$LATEST_BENCHMARK"
if [ $? -ne 0 ]; then
    echo "ERROR: Enhanced summary generation failed!"
    exit 1
fi

# Get the most recent analysis files
LATEST_ANALYSIS=$(ls -t enhanced_pfas_analysis_*.html 2>/dev/null | head -n1)
LATEST_TIMING_ANALYSIS=$(ls -t timing_analysis_*.html 2>/dev/null | head -n1)

echo ""
echo "============================================"
echo "  ENHANCED BENCHMARK COMPLETE!"
echo "============================================"
echo ""
echo "Generated Files:"
echo "  📊 Timing Benchmark Data: $LATEST_TIMING"
echo "  📈 Timing Analysis Report: $LATEST_TIMING_ANALYSIS"
echo "  🧪 Comprehensive Benchmark Data: $LATEST_BENCHMARK"
echo "  📋 Enhanced Analysis Report: $LATEST_ANALYSIS"
echo "  🎯 Performance Heatmaps (PNG/SVG)"
echo "  🔗 Sankey Diagrams (PNG/SVG)"
echo "  ⏱️  Timing Visualizations (PNG/SVG)"
echo "  📊 Statistical Error Bar Charts"
echo ""
echo "Main Analysis Reports:"
echo "  - Timing Analysis: $LATEST_TIMING_ANALYSIS"
echo "  - Enhanced Analysis: $LATEST_ANALYSIS"
echo ""
echo "🎉 Complete pipeline executed successfully!"
echo "   ✅ Statistical timing analysis with 10x iterations"
echo "   ✅ Comprehensive molecule testing"
echo "   ✅ System specifications included"
echo "   ✅ Enhanced visualizations generated"