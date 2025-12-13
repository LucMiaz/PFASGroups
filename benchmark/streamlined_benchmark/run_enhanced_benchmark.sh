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

# Step 1: Run the enhanced benchmark
echo "[1/3] Running Enhanced PFAS Benchmark..."
echo "----------------------------------------"
python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "ERROR: Enhanced benchmark failed!"
    exit 1
fi

# Get the most recent benchmark file
LATEST_BENCHMARK=$(ls -t pfas_enhanced_benchmark_*.json 2>/dev/null | head -n1)

if [ -z "$LATEST_BENCHMARK" ]; then
    echo "ERROR: No benchmark results file found!"
    exit 1
fi

echo "Latest benchmark file: $LATEST_BENCHMARK"
echo ""

# Step 2: Run the enhanced analysis
echo "[2/3] Running Enhanced Analysis..."
echo "----------------------------------"
python enhanced_analysis.py "$LATEST_BENCHMARK"
if [ $? -ne 0 ]; then
    echo "ERROR: Enhanced analysis failed!"
    exit 1
fi

# Step 3: Show summary
echo "[3/3] Generating Summary..."
echo "----------------------------"
python benchmark_summary.py "$LATEST_BENCHMARK"
if [ $? -ne 0 ]; then
    echo "ERROR: Summary generation failed!"
    exit 1
fi

# Get the most recent analysis file
LATEST_ANALYSIS=$(ls -t enhanced_pfas_analysis_*.html 2>/dev/null | head -n1)

echo ""
echo "============================================"
echo "  ENHANCED BENCHMARK COMPLETE!"
echo "============================================"
echo ""
echo "Generated Files:"
echo "  - Benchmark Data: $LATEST_BENCHMARK"
echo "  - Analysis Report: $LATEST_ANALYSIS"
echo "  - Performance Heatmaps (PNG/SVG)"
echo "  - Sankey Diagrams (PNG/SVG)"
echo ""
echo "Analysis report saved as: $LATEST_ANALYSIS"
echo ""
echo "Pipeline completed successfully!"