#!/bin/bash
# Unified PFAS Benchmark Runner (Linux/macOS)
# Runs all benchmark types and generates unified HTML report

echo "🚀 PFAS BENCHMARK SUITE - UNIFIED RUNNER"
echo "========================================="
echo ""

# Check if we're in the right directory
if [ ! -f "enhanced_pfas_benchmark.py" ]; then
    echo "❌ Error: enhanced_pfas_benchmark.py not found!"
    echo "   Please run this script from the benchmark/streamlined_benchmark directory"
    exit 1
fi

echo "📋 Running all PFAS benchmarks..."
echo ""

# Functional Groups Benchmark
echo "1️⃣ Running Functional Groups Benchmark..."
echo "1" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Functional Groups Benchmark failed"
    exit 1
fi
echo "✅ Functional Groups Benchmark completed"
echo ""

# OECD Validation Benchmark
echo "2️⃣ Running OECD Validation Benchmark..."
echo "2" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ OECD Validation Benchmark failed"
    exit 1
fi
echo "✅ OECD Validation Benchmark completed"
echo ""

# Timing Performance Benchmark
echo "3️⃣ Running Timing Performance Benchmark..."
echo "3" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Timing Performance Benchmark failed"
    exit 1
fi
echo "✅ Timing Performance Benchmark completed"
echo ""

# Non-Fluorinated Exclusion Benchmark
echo "4️⃣ Running Non-Fluorinated Exclusion Benchmark..."
echo "4" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Non-Fluorinated Exclusion Benchmark failed"
    exit 1
fi
echo "✅ Non-Fluorinated Exclusion Benchmark completed"
echo ""

# Complex Branched Structures Benchmark
echo "5️⃣ Running Complex Branched Structures Benchmark..."
echo "5" | python enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Complex Branched Structures Benchmark failed"
    exit 1
fi
echo "✅ Complex Branched Structures Benchmark completed"
echo ""

# Generate Unified Report
echo "📊 Generating Unified HTML Report..."
python generate_unified_report.py
if [ $? -ne 0 ]; then
    echo "❌ Unified Report generation failed"
    exit 1
fi

echo ""
echo "🎉 ALL BENCHMARKS COMPLETED SUCCESSFULLY!"
echo ""
echo "📋 Results Summary:"
echo "   • Functional Groups: Tests basic PFAS group detection"
echo "   • OECD Validation: Validates against OECD database"
echo "   • Timing Performance: Measures execution speed scaling"
echo "   • Non-Fluorinated: Ensures proper exclusion of non-PFAS"
echo "   • Complex Branched: Tests complex molecular structures"
echo ""
echo "📄 Unified Report: $(ls unified_pfas_benchmark_report_*.html | tail -1)"
echo "🌐 Open the HTML file in your browser to view detailed results"
echo ""
echo "✨ Benchmark Suite Complete!"