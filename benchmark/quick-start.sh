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