#!/bin/bash
# Quick Unified Report Generator
# =============================
# Generates unified HTML report from existing benchmark results
# Usage: ./generate_report.sh

echo "🚀 QUICK UNIFIED REPORT GENERATOR"
echo "=================================="

# Check if we're in the right directory
if [ ! -f "generate_unified_report.py" ]; then
    echo "ERROR: generate_unified_report.py not found!"
    echo "Please run this script from the streamlined_benchmark directory"
    exit 1
fi

echo "📊 Generating unified HTML report from existing benchmark data..."
python generate_unified_report.py

if [ $? -eq 0 ]; then
    # Find the most recent unified report
    unified_report=$(ls -t unified_pfas_benchmark_report_*.html 2>/dev/null | head -n1)
    
    if [[ -n "$unified_report" ]]; then
        echo ""
        echo "🎉 SUCCESS: Unified report generated!"
        echo "📄 Report file: $unified_report"
        echo ""
        echo "🌐 Open this file in your web browser to view the complete analysis dashboard"
        echo "📊 The report includes all available benchmark results in a single view"
        exit 0
    else
        echo "❌ ERROR: Report file not found after generation"
        exit 1
    fi
else
    echo "❌ ERROR: Report generation failed"
    echo "📋 Make sure you have benchmark result files (*.json) in this directory"
    echo "📋 Run the full benchmark pipeline first if needed: ./run_enhanced_benchmark.sh"
    exit 1
fi