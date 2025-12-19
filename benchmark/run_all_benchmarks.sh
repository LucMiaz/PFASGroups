#!/bin/bash
# Unified PFAS Benchmark Runner (Linux/macOS)
# Runs all benchmark types and generates unified HTML report
# Now includes database integration for the review app

echo "🚀 PFAS BENCHMARK SUITE - UNIFIED RUNNER"
echo "========================================="
echo ""

# Check if we're in the right directory
if [ ! -f "enhanced_pfas_benchmark.py" ]; then
    echo "❌ Error: enhanced_pfas_benchmark.py not found!"
    echo "   Please run this script from the benchmark directory"
    exit 1
fi

# Create output directories
mkdir -p imgs html data

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

# Organize output files
echo "📁 Organizing output files..."
# Move data files
mv pfas_*_benchmark_*.json data/ 2>/dev/null || true
# Move HTML files
mv *.html html/ 2>/dev/null || true
# Move image files
mv *.png imgs/ 2>/dev/null || true
mv *.svg imgs/ 2>/dev/null || true

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
echo "📊 Updating Review App Database..."
if [ -f "enhanced_benchmark_with_db.py" ]; then
    python enhanced_benchmark_with_db.py
    echo "✅ Database updated with latest benchmark results"
else
    echo "⚠️  Database integration script not found - skipping database update"
fi
echo ""
echo "📄 Unified Report: html/$(ls html/unified_pfas_benchmark_report_*.html | tail -1 | xargs basename)"
echo "🌐 Open the HTML file in your browser to view detailed results"
echo "📁 Files organized in: data/ (JSON), html/ (reports), imgs/ (plots)"
echo ""
echo "🔬 Review App:"
echo "   • Run 'cd review-app && ./start-dev.sh' to start the review interface"
echo "   • Manual validation at http://localhost:3000"
echo "   • Database file: review-app/database/pfas_benchmark.db"
echo ""
echo "✨ Benchmark Suite Complete!"