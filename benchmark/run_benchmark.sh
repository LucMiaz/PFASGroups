#!/bin/bash

# Runner script for PFAS-atlas benchmark
# This script activates the correct environment and runs the benchmark

set -e

echo "🚀 Starting PFAS-atlas benchmark..."
echo "================================="

# Check if mamba is available
if ! command -v mamba &> /dev/null; then
    echo "❌ mamba command not found. Please install mamba or modify script to use conda."
    exit 1
fi

# Run with mamba run to use the pfasatlas environment
echo "🔧 Using mamba run with pfasatlas environment..."

echo "✅ Mamba configured to use pfasatlas environment"

# Check if required files exist
SPECIFICITY_FILE="/home/luc/git/PFASGroups/PFASgroups/tests/results/specificity_test_results.csv"
PFAS_ATLAS_DIR="/home/luc/git/PFAS-atlas"

if [[ ! -f "$SPECIFICITY_FILE" ]]; then
    echo "❌ Specificity test file not found: $SPECIFICITY_FILE"
    exit 1
fi

if [[ ! -d "$PFAS_ATLAS_DIR" ]]; then
    echo "❌ PFAS-atlas directory not found: $PFAS_ATLAS_DIR"
    exit 1
fi

echo "✅ Input files validated"

# Run the benchmark
echo "🧪 Running benchmark comparison..."
cd "$(dirname "$0")"
mamba run -n pfasatlas python pfas_atlas_benchmark.py

echo "🎉 Benchmark complete!"
echo ""
echo "📁 Results saved in:"
echo "   - pfas_atlas_benchmark_report.md (detailed report)"
echo "   - pfas_atlas_benchmark_results.csv (full results)"
echo "   - pfas_atlas_benchmark_summary.json (summary statistics)"