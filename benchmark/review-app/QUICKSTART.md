# 🚀 PFAS Benchmark Reviewer - Quick Start Guide

## ✅ Setup Complete!

The application is now ready to use. Here's how to get started:

## Starting the Application

### Option 1: Development Mode (Recommended)
Best for active development and reviewing. Hot reload enabled.

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark\review-app
.\start-dev.ps1
```

This will:
- Start the backend server on **http://localhost:5000**
- Start the React dev server on **http://localhost:3000**
- Open **http://localhost:3000** in your browser

### Option 2: Production Mode
Optimized build for production use.

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark\review-app
.\start-prod.ps1
```

This will:
- Build the React app
- Start the server on **http://localhost:5000**
- Open **http://localhost:5000** in your browser

### Option 3: Server Only
Just the backend API (if you already built the React app).

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark\review-app
.\start-server-only.ps1
```

## 📊 Features Available

- ✅ Molecule visualization with RDKit.js
- ✅ Pagination and filtering by dataset
- ✅ Manual review interface with quick action buttons
- ✅ Prioritize misclassified molecules
- ✅ Real-time accuracy computation
- ✅ Dashboard with statistics
- ✅ Data export (JSON/CSV)

## 🗄️ Database

The application uses SQLite with the database file at:
```
c:\Users\luc\git\PFASGroups\benchmark\review-app\database\pfas_benchmark.db
```

## 📝 Next Steps

1. **Import benchmark data** (if you have existing data):
   ```powershell
   node scripts/import-benchmark-data.js
   ```

2. **Start the app** using one of the methods above

3. **Navigate to the web interface** and start reviewing!

## 🔧 Troubleshooting

- **Port already in use**: Kill the process using port 5000 or 3000
- **React build errors**: Try `cd client && npm install` again
- **Database errors**: Delete the database file to recreate it

## 🛠️ Development Commands

```powershell
# Install/reinstall dependencies
cd client
npm install

# Build React app manually
cd client
npm run build

# Run backend tests (if available)
npm test
```

## 📞 Need Help?

- Check the README.md in the review-app directory
- Review server.js for API endpoints
- Check browser console for React errors
