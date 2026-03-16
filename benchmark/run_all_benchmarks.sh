#!/bin/bash
# Unified PFAS Benchmark Runner (Linux/macOS)
# Runs all benchmark types and generates unified HTML report
# Now includes database integration for the review app
# Usage: ./run_all_benchmarks.sh [--reuse-timing] [--quick]

# Parse command line arguments
REUSE_TIMING=false
QUICK_RUN=false
for arg in "$@"; do
    case $arg in
        --reuse-timing)
            REUSE_TIMING=true
            shift
            ;;
        --quick)
            QUICK_RUN=true
            shift
            ;;
        --help)
            echo "Usage: ./run_all_benchmarks.sh [--reuse-timing] [--quick]"
            echo ""
            echo "Options:"
            echo "  --reuse-timing    Load and extend previous timing benchmark results"
            echo "  --quick           Run a reduced-size benchmark for smoke testing"
            echo "  --help           Show this help message"
            exit 0
            ;;
    esac
done

echo "🚀 PFAS BENCHMARK SUITE - UNIFIED RUNNER"
echo "========================================="
if [ "$REUSE_TIMING" = true ]; then
    echo "♻️  Reusing previous timing results mode enabled"
fi
if [ "$QUICK_RUN" = true ]; then
    echo "⚡ Quick-run mode enabled (reduced sizes)"
    export PFAS_BENCH_QUICK=1
    export PFAS_BENCH_REPLICATES=3
    export PFAS_BENCH_TIMING_MOLECULES=80
    export PFAS_BENCH_TIMING_ITERATIONS=2
    export PFAS_BENCH_NONFLUOR=5
    export PFAS_BENCH_COMPLEX=5
    export PFAS_BENCH_OECD_LIMIT=200
    export PFAS_BENCH_BRANCHED_COMPONENTS=2
    export PFAS_BENCH_BRANCHED_GROUPS=6
    export PFAS_BENCH_BRANCHED_DISTANCE_MAX=1
fi
echo ""

# Ensure pfasatlas environment is active (mamba/conda)
if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "pfasatlas" ]; then
    if command -v mamba >/dev/null 2>&1; then
        eval "$(mamba shell.bash hook)"
        mamba activate pfasatlas
    elif command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate pfasatlas
    else
        echo "⚠️  mamba/conda not found; please activate 'pfasatlas' manually"
    fi
fi

# Check if we're in the right directory
if [ ! -f "scripts/enhanced_pfas_benchmark.py" ]; then
    echo "❌ Error: scripts/enhanced_pfas_benchmark.py not found!"
    echo "   Please run this script from the benchmark directory"
    exit 1
fi

# Create output directories
mkdir -p scripts reports data imgs

echo "📋 Running all PFAS benchmarks..."
echo ""

# Functional Groups Benchmark
echo "1️⃣ Running Functional Groups Benchmark..."
echo "1" | python scripts/enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Functional Groups Benchmark failed"
    exit 1
fi
echo "✅ Functional Groups Benchmark completed"
echo ""

# OECD Validation Benchmark
echo "2️⃣ Running OECD Validation Benchmark..."
echo "2" | python scripts/enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ OECD Validation Benchmark failed"
    exit 1
fi
echo "✅ OECD Validation Benchmark completed"
echo ""

# Timing Performance Benchmarks (3 profiles)
echo "3️⃣ Running Timing Performance Benchmarks (3 profiles)..."

echo "   • Profile: full (all features)"
if [ "$REUSE_TIMING" = true ]; then
    echo "♻️  Reusing previous timing results..."
    printf "3\ny\n" | env PFAS_BENCH_TIMING_PROFILE=full python scripts/enhanced_pfas_benchmark.py
else
    printf "3\nn\n" | env PFAS_BENCH_TIMING_PROFILE=full python scripts/enhanced_pfas_benchmark.py
fi
if [ $? -ne 0 ]; then
    echo "❌ Timing Performance Benchmark failed (full profile)"
    exit 1
fi

echo "   • Profile: no_resistance (skip effective graph resistance)"
if [ "$REUSE_TIMING" = true ]; then
    echo "♻️  Reusing previous timing results..."
    printf "3\ny\n" | env PFAS_BENCH_TIMING_PROFILE=no_resistance python scripts/enhanced_pfas_benchmark.py
else
    printf "3\nn\n" | env PFAS_BENCH_TIMING_PROFILE=no_resistance python scripts/enhanced_pfas_benchmark.py
fi
if [ $? -ne 0 ]; then
    echo "❌ Timing Performance Benchmark failed (no_resistance profile)"
    exit 1
fi

echo "   • Profile: no_metrics (skip component graph metrics)"
if [ "$REUSE_TIMING" = true ]; then
    echo "♻️  Reusing previous timing results..."
    printf "3\ny\n" | env PFAS_BENCH_TIMING_PROFILE=no_metrics python scripts/enhanced_pfas_benchmark.py
else
    printf "3\nn\n" | env PFAS_BENCH_TIMING_PROFILE=no_metrics python scripts/enhanced_pfas_benchmark.py
fi
if [ $? -ne 0 ]; then
    echo "❌ Timing Performance Benchmark failed (no_metrics profile)"
    exit 1
fi

echo "✅ Timing Performance Benchmarks completed"
echo ""

# Timing Profile Comparison (plots + LaTeX)
echo "🧪 Analyzing timing profiles (comparison + LaTeX)..."
python scripts/compare_timing_profiles.py
if [ $? -ne 0 ]; then
    echo "❌ Timing profile comparison failed"
    exit 1
fi
echo "✅ Timing profile comparison completed"
echo ""

# Non-Fluorinated Exclusion Benchmark
echo "4️⃣ Running Non-Fluorinated Exclusion Benchmark..."
echo "4" | python scripts/enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Non-Fluorinated Exclusion Benchmark failed"
    exit 1
fi
echo "✅ Non-Fluorinated Exclusion Benchmark completed"
echo ""

# Complex Branched Structures Benchmark
echo "5️⃣ Running Complex Branched Structures Benchmark..."
echo "5" | python scripts/enhanced_pfas_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ Complex Branched Structures Benchmark failed"
    exit 1
fi
echo "✅ Complex Branched Structures Benchmark completed"
echo ""

# Highly Branched Compounds Test
echo "6️⃣ Running Highly Branched Compounds Test..."
python scripts/test_highly_branched.py
if [ $? -ne 0 ]; then
    echo "❌ Highly Branched Compounds Test failed"
    exit 1
fi
echo "✅ Highly Branched Compounds Test completed"
echo ""

# Telomer Validation
echo "7️⃣ Running Telomer Detection Validation..."
python scripts/validate_telomers.py
if [ $? -ne 0 ]; then
    echo "❌ Telomer Validation failed"
    exit 1
fi
echo "✅ Telomer Validation completed"
echo ""

# PFAS Definitions Benchmark
echo "8️⃣ Running PFAS Definitions Benchmark..."
python scripts/benchmark_pfas_definitions.py
if [ $? -ne 0 ]; then
    echo "❌ PFAS Definitions Benchmark failed"
    exit 1
fi
echo "✅ PFAS Definitions Benchmark completed"
echo ""

# Analyze PFAS Definitions Results
echo "9️⃣ Analyzing PFAS Definitions Results..."
python scripts/analyze_definitions_benchmark.py
if [ $? -ne 0 ]; then
    echo "❌ PFAS Definitions Analysis failed"
    exit 1
fi
echo "✅ PFAS Definitions Analysis completed"
echo ""

# Add Test Metadata to JSON files
echo "🏷️  Adding test metadata to PFAS groups and definitions..."
if [ -f "scripts/add_test_metadata.py" ]; then
    python scripts/add_test_metadata.py
    if [ $? -ne 0 ]; then
        echo "⚠️  Test metadata addition failed (non-critical)"
    fi
else
    echo "⚠️  add_test_metadata.py not found (non-critical)"
fi
echo ""

# Generate Telomer Validation LaTeX
echo "📝 Generating telomer validation LaTeX content..."
python scripts/generate_telomer_latex.py
if [ $? -ne 0 ]; then
    echo "⚠️  Telomer LaTeX generation failed (non-critical)"
fi
echo ""

# Generate Unified Report
echo "📊 Generating Unified HTML Report..."
python scripts/generate_unified_report.py
if [ $? -ne 0 ]; then
    echo "❌ Unified Report generation failed"
    exit 1
fi
echo "✅ Unified Report generated"
echo ""

# Run analysis scripts BEFORE organizing files
echo "📊 Running Analysis Scripts..."
echo ""

# Timing Analysis
echo "⏱️  Analyzing timing performance..."
TIMING_FILE=$(ls data/pfas_timing_benchmark_*.json 2>/dev/null | tail -1)
if [ -f "$TIMING_FILE" ]; then
    python scripts/analyze_timing.py "$TIMING_FILE"
    if [ $? -eq 0 ]; then
        echo "✅ Timing analysis completed"
    else
        echo "⚠️  Timing analysis failed"
    fi
else
    echo "⚠️  No timing benchmark file found"
fi
echo ""

# Timing Models Analysis
echo "📈 Analyzing timing models (exponential fit)..."
if [ -f "$TIMING_FILE" ]; then
    python scripts/analyze_timing_models.py
    if [ $? -eq 0 ]; then
        echo "✅ Timing models analysis completed"
    else
        echo "⚠️  Timing models analysis failed"
    fi
else
    echo "⚠️  No timing benchmark file found"
fi
echo ""

# Definitions Benchmark Analysis
echo "📋 Analyzing PFAS definitions benchmark..."
DEFINITION_FILE=$(ls data/pfas_definitions_benchmark_*.json 2>/dev/null | tail -1)
if [ -f "$DEFINITION_FILE" ]; then
    python scripts/analyze_definitions_benchmark.py "$DEFINITION_FILE"
    if [ $? -eq 0 ]; then
        echo "✅ Definitions analysis completed"
    else
        echo "⚠️  Definitions analysis failed"
    fi
else
    echo "⚠️  No definitions benchmark file found (this is optional)"
fi
echo ""

# Complex Branched Analysis
echo "🧬 Analyzing complex branched structures..."
COMPLEX_FILE=$(ls data/pfas_complex_branched_benchmark_*.json 2>/dev/null | tail -1)
if [ -f "$COMPLEX_FILE" ]; then
    if [ -f "scripts/analyze_complex.py" ]; then
        python scripts/analyze_complex.py "$COMPLEX_FILE"
        if [ $? -eq 0 ]; then
            echo "✅ Complex branched analysis completed"
        else
            echo "⚠️  Complex branched analysis failed"
        fi
    else
        echo "⚠️  analyze_complex.py not found (non-critical)"
    fi
else
    echo "⚠️  No complex branched benchmark file found"
fi
echo ""

# Enhanced Analysis
echo "🔬 Analyzing enhanced functional groups and OECD..."
ENHANCED_FILE=$(ls data/pfas_enhanced_benchmark_*.json 2>/dev/null | tail -1)
OECD_FILE=$(ls data/pfas_oecd_benchmark_*.json 2>/dev/null | tail -1)
if [ -f "$ENHANCED_FILE" ] && [ -f "$OECD_FILE" ]; then
    python scripts/enhanced_analysis.py "$ENHANCED_FILE" "$OECD_FILE"
    if [ $? -eq 0 ]; then
        echo "✅ Enhanced analysis completed"
    else
        echo "⚠️  Enhanced analysis failed"
    fi
else
    echo "⚠️  Missing enhanced or OECD benchmark files"
    [ ! -f "$ENHANCED_FILE" ] && echo "    Missing: enhanced benchmark"
    [ ! -f "$OECD_FILE" ] && echo "    Missing: OECD benchmark"
fi
echo ""

# Organize output files
echo "📁 Organizing output files..."
# Keep data files in data/ directory (they're already there)
# Move HTML files
mv *.html reports/ 2>/dev/null || true
mv unified_pfas_benchmark_report_*.html reports/ 2>/dev/null || true
# Move image files
mv *.png imgs/ 2>/dev/null || true
mv *.svg imgs/ 2>/dev/null || true
# Move any analysis JSON files to data/
mv *_analysis.json data/ 2>/dev/null || true
mv *_results.json data/ 2>/dev/null || true
mv *_benchmark_*.json data/ 2>/dev/null || true
echo ""
echo "🎉 ALL BENCHMARKS COMPLETED SUCCESSFULLY!"
echo ""
echo "📋 Results Summary:"
echo "   • Functional Groups: Tests basic PFAS group detection"
echo "   • OECD Validation: Validates against OECD database"
echo "   • Timing Performance: Measures execution speed scaling"
echo "   • Non-Fluorinated: Ensures proper exclusion of non-PFAS"
echo "   • Complex Branched: Tests complex molecular structures"
echo "   • Highly Branched: Tests functional groups on perfluorinated components"
echo "   • Telomer Validation: Tests detection of fluorotelomers on PubChem dataset"
echo "   • PFAS Definitions: Benchmarks 5 PFAS definitions (OECD, EU, OPPT, UK, PFASSTRUCTv5)"
echo "   • Comprehensive Statistics: LaTeX tables and benchmark summary (reports/benchmark_summary.json)"
echo ""

# Telomer Validation Report
echo "🧪 Generating telomer validation report..."
if [ -f "data/telomer_validation_results.json" ]; then
    python scripts/generate_telomer_report.py
    if [ $? -eq 0 ]; then
        echo "✅ Telomer validation report completed"
    else
        echo "⚠️  Telomer validation report failed"
    fi
else
    echo "⚠️  No telomer validation data found"
fi
echo ""

# Comprehensive Benchmark Analysis
echo "📊 Generating comprehensive benchmark statistics..."
if [ -f "scripts/analyze_benchmarks_simple.py" ]; then
    python scripts/analyze_benchmarks_simple.py
    if [ $? -eq 0 ]; then
        echo "✅ Comprehensive benchmark analysis completed"
    else
        echo "⚠️  Comprehensive benchmark analysis failed"
    fi
else
    echo "⚠️  analyze_benchmarks_simple.py not found (non-critical)"
fi
echo ""

# LaTeX Tables Generation for Article
echo "📝 Generating LaTeX tables for article..."
if [ -f "scripts/generate_latex_tables.py" ]; then
    python scripts/generate_latex_tables.py
    if [ $? -eq 0 ]; then
        echo "✅ LaTeX tables generated successfully"
        echo "   • Main content: reports/pfasgroups_latex_results.tex"
        echo "   • Summary: reports/pfasgroups_latex_summary.tex"
    else
        echo "⚠️  LaTeX generation failed"
    fi
else
    echo "⚠️  generate_latex_tables.py not found"
fi
echo ""

# Organize analysis results
echo "📁 Organizing analysis results..."

echo "✅ Analysis results organized"
echo ""

# Import data to database
echo "💾 Importing data to Review App database..."
cd review-app

# Delete existing database if present (fresh schema on import)
if [ -f "database/pfas_benchmark.db" ]; then
    echo "🗑️  Deleting existing database: database/pfas_benchmark.db"
    rm -f database/pfas_benchmark.db
fi

# Import benchmark data
echo "📥 Importing benchmark data..."
node scripts/import-benchmark-data.js
if [ $? -eq 0 ]; then
    echo "✅ Data imported successfully"
    
    # Calculate molecular formulas
    echo "🧪 Calculating molecular formulas..."
    python scripts/calculate-formulas.py
    if [ $? -eq 0 ]; then
        echo "✅ Molecular formulas calculated"
    else
        echo "⚠️  Formula calculation failed (non-critical)"
    fi
else
    echo "❌ Data import failed"
    cd ..
    exit 1
fi

cd ..
echo ""

echo "📊 Database Update Complete!"
echo "   • Database: review-app/database/pfas_benchmark.db"
echo "   • Analysis reports: reports"
echo ""
REPORT_FILE=$(ls reports/unified_pfas_benchmark_report_*.html 2>/dev/null | tail -1)
if [ -n "$REPORT_FILE" ]; then
    echo "📄 Unified Report: $REPORT_FILE"
else
    echo "📄 Unified Report: reports/ (no report found)"
fi
echo "🌐 Open the HTML file in your browser to view detailed results"
echo "📁 Files organized in: data/ (JSON), html/ (reports), imgs/ (plots)"
echo ""
echo "🔬 Review App:"
echo "   • Navigate to: cd review-app"
echo "   • Start server: node server.js"
echo "   • Open browser: http://localhost:5000"
echo "   • View Analysis Reports tab for timing, complex, and enhanced analysis"
echo ""
echo "✨ Benchmark Suite Complete!"
