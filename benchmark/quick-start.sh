#!/bin/bash

echo "🚀 Quick Start: PFAS Benchmark Reviewer"
echo "======================================"

# Check if we're in the right directory
if [ ! -f "setup-review-app.sh" ]; then
    echo "❌ Error: Please run this from the benchmark directory"
    echo "   cd /home/luc/git/PFASGroups/benchmark"
    exit 1
fi

echo ""
echo "📋 Setup checklist:"
echo "  □ Install Node.js dependencies"
echo "  □ Setup SQLite database"  
echo "  □ Import existing benchmark data"
echo "  □ Start the review application"
echo ""

# Run setup
echo "1️⃣ Running setup script..."
./setup-review-app.sh

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Setup completed successfully!"
    echo ""
    
    # Generate timing analysis report
    echo "3️⃣ Generating timing analysis report..."
    if command -v conda &> /dev/null; then
        conda run -n chem python generate_timing_report.py benchmark_timing_report.html 2>/dev/null
    else
        python generate_timing_report.py benchmark_timing_report.html 2>/dev/null
    fi
    
    if [ -f "benchmark_timing_report.html" ]; then
        echo "✅ Timing report generated: benchmark_timing_report.html"
        echo "   View it by opening the file in your browser"
    else
        echo "⚠️  Timing report generation skipped (optional)"
    fi
    echo ""
    echo "🎯 What's next?"
    echo ""
    echo "Option 1: Development mode (recommended for reviewing)"
    echo "  cd review-app"
    echo "  ./start-dev.sh"
    echo "  Open http://localhost:3000"
    echo ""
    echo "Option 2: Production mode"
    echo "  cd review-app"  
    echo "  ./start-prod.sh"
    echo "  Open http://localhost:5000"
    echo ""
    echo "📊 Features available:"
    echo "  ✓ Prioritized display of misclassified molecules"
    echo "  ✓ Manual review with click buttons (✅❌🤷)"
    echo "  ✓ Enhanced performance metrics"
    echo "  ✓ Real-time accuracy computation"
    echo "  ✓ Export capabilities (JSON/CSV)"
    echo "  ✓ RDKit.js molecular visualization"
    echo ""
    echo "🔬 Pro tip: Enable 'Prioritize Misclassified' filter to review"
    echo "   the most important molecules first!"
    echo ""
else
    echo ""
    echo "❌ Setup failed. Please check the error messages above."
    echo "   You may need to install Node.js or resolve dependency issues."
    exit 1
fi