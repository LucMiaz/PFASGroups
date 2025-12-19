#!/bin/bash
set -e

echo "🚀 Setting up PFAS Benchmark Reviewer Application"
echo "================================================"

# Change to review-app directory
cd "$(dirname "$0")/review-app"

echo "📦 Installing Node.js dependencies..."
npm install

echo "📦 Installing React app dependencies..."
# Check if the automated React app creation completed
if [ -d "client" ]; then
    echo "Using auto-generated React app..."
    cd client
    
    # Install additional dependencies
    npm install axios @rdkit/rdkit bootstrap react-bootstrap react-router-dom react-paginate react-select
    
    # Copy our custom React components
    if [ -d "../client-src/src" ]; then
        echo "📝 Copying custom React components..."
        cp -r ../client-src/src/* src/
        cp ../client-src/public/index.html public/
        cp ../client-src/package.json package-temp.json
        
        # Merge package.json dependencies
        node -e "
        const fs = require('fs');
        const existing = JSON.parse(fs.readFileSync('package.json', 'utf8'));
        const custom = JSON.parse(fs.readFileSync('package-temp.json', 'utf8'));
        
        existing.dependencies = { ...existing.dependencies, ...custom.dependencies };
        existing.proxy = custom.proxy;
        
        fs.writeFileSync('package.json', JSON.stringify(existing, null, 2));
        fs.unlinkSync('package-temp.json');
        console.log('✓ Package.json updated');
        "
        
        # Install updated dependencies
        npm install
    fi
    
    cd ..
else
    echo "Using pre-built React components..."
    # If auto-generation didn't complete, use our pre-built structure
    mv client-src client
    cd client
    npm install
    cd ..
fi

echo "🗄️  Initializing database..."
node -e "
const Database = require('./database/database.js');
const db = new Database();
console.log('✓ Database initialized');
"

echo "📊 Importing existing benchmark data..."
node scripts/import-benchmark-data.js || echo "⚠️  Data import failed - will continue without existing data"

echo "🔧 Creating startup scripts..."

# Create development startup script
cat > start-dev.sh << 'EOF'
#!/bin/bash
echo "🚀 Starting PFAS Benchmark Reviewer in development mode..."
echo "Server: http://localhost:5000"
echo "React Dev Server: http://localhost:3000 (if available)"
echo ""
npm run dev
EOF
chmod +x start-dev.sh

# Create production startup script
cat > start-prod.sh << 'EOF'
#!/bin/bash
echo "🚀 Building and starting PFAS Benchmark Reviewer..."
cd client && npm run build && cd ..
echo "✓ React app built"
echo "🌐 Starting server at http://localhost:5000"
npm start
EOF
chmod +x start-prod.sh

# Create data import script
cat > import-latest-data.sh << 'EOF'
#!/bin/bash
echo "📊 Importing latest benchmark data..."
node scripts/import-benchmark-data.js
echo "✅ Data import completed"
EOF
chmod +x import-latest-data.sh

echo ""
echo "✅ Setup completed successfully!"
echo ""
echo "🎯 Quick Start:"
echo "  Development mode: ./start-dev.sh"
echo "  Production mode:  ./start-prod.sh"
echo "  Import data:      ./import-latest-data.sh"
echo ""
echo "📝 Manual Steps:"
echo "  1. Run ./import-latest-data.sh to load existing benchmark data"
echo "  2. Run ./start-dev.sh to start the development server"
echo "  3. Open http://localhost:3000 (dev) or http://localhost:5000 (prod)"
echo ""
echo "🔧 Development:"
echo "  - Server runs on port 5000"
echo "  - React dev server on port 3000 (proxies API calls to 5000)"
echo "  - Database: SQLite file in database/pfas_benchmark.db"
echo ""
echo "📊 Features Available:"
echo "  ✓ Molecule visualization with RDKit.js"
echo "  ✓ Pagination and filtering"
echo "  ✓ Manual review interface"
echo "  ✓ Accuracy computation"
echo "  ✓ Data export (JSON/CSV)"
echo "  ✓ Dashboard with statistics"