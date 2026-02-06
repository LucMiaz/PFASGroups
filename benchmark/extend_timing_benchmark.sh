#!/bin/bash
# Extend Timing Benchmark
# Loads the most recent timing results and adds more molecules
# Usage: ./extend_timing_benchmark.sh [num_molecules] [iterations]

NUM_MOLECULES=${1:-2500}
ITERATIONS=${2:-5}

echo "🔄 EXTEND TIMING BENCHMARK"
echo "=========================="
echo "Adding ${NUM_MOLECULES} molecules with ${ITERATIONS} iterations each"
echo ""

# Check if we're in the right directory
if [ ! -f "scripts/enhanced_pfas_benchmark.py" ]; then
    echo "❌ Error: scripts/enhanced_pfas_benchmark.py not found!"
    echo "   Please run this script from the benchmark directory"
    exit 1
fi

# Run timing benchmark with reuse option
printf "3\ny\n" | python scripts/enhanced_pfas_benchmark.py

if [ $? -ne 0 ]; then
    echo "❌ Timing benchmark extension failed"
    exit 1
fi

echo ""
echo "✅ Timing benchmark extended successfully!"
echo ""
echo "📊 Next steps:"
echo "   1. Import to database: cd review-app/scripts && node import-benchmark-data.js"
echo "   2. Generate report: python scripts/generate_enhanced_timing_report.py"
echo "   3. Analyze complexity: python scripts/analyze_timing_with_complexity.py data/pfas_timing_benchmark_*.json"
