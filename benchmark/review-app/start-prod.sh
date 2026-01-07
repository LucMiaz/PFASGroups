#!/bin/bash
echo "🚀 Building and starting PFAS Benchmark Reviewer..."
cd client && npm run build && cd ..
echo "✓ React app built"
echo "🌐 Starting server at http://localhost:5000"
npm start
