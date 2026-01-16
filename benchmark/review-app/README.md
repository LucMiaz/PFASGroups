# PFAS Benchmark Reviewer

A comprehensive web application for manually reviewing and validating PFAS classification results from PFASGroups and PFAS-Atlas algorithms.

## 🎯 Features

- **Interactive Molecule Visualization**: RDKit.js integration for 2D molecular structure display
- **Manual Classification Review**: Click-button interface to validate algorithm results
- **Pagination & Filtering**: Efficient browsing through large datasets
- **Multi-Dataset Support**: OECD, Enhanced, Timing, Complex Branched, and Non-Fluorinated datasets
- **Review Status Tracking**: Continue reviewing where you left off
- **Accuracy Computation**: Real-time accuracy metrics against manually reviewed entries
- **Data Export**: Export reviews in JSON and CSV formats
- **SQLite Database**: Persistent storage with efficient querying
- **Dashboard Analytics**: Overview of review progress and accuracy statistics

## 🚀 Quick Start

### Prerequisites
- Node.js (v14 or higher)
- Python 3.7+
- Existing PFAS benchmark data (JSON files)

### Setup and Installation

1. **Run the setup script:**
   ```bash
   cd /home/luc/git/PFASGroups/benchmark
   ./setup-review-app.sh
   ```

2. **Import existing benchmark data:**
   ```bash
   cd review-app
   ./import-latest-data.sh
   ```

3. **Start the application:**
   ```bash
   # Development mode (recommended)
   ./start-dev.sh
   
   # Or production mode
   ./start-prod.sh
   ```

4. **Open in browser:**
   - Development: http://localhost:3000
   - Production: http://localhost:5000

## 📁 Project Structure

```
review-app/
├── server.js                 # Express.js server with API endpoints
├── package.json              # Node.js dependencies
├── database/
│   ├── database.js           # SQLite database setup and helpers
│   └── pfas_benchmark.db     # SQLite database file
├── scripts/
│   └── import-benchmark-data.js  # Data import utilities
├── client/                   # React frontend application
│   ├── src/
│   │   ├── App.js           # Main React app component
│   │   ├── components/
│   │   │   ├── Dashboard.js          # Overview dashboard
│   │   │   ├── MoleculeReviewer.js   # Main review interface
│   │   │   ├── MoleculeViewer.js     # RDKit.js molecule display
│   │   │   └── AccuracyReport.js     # Accuracy analytics
│   │   └── App.css          # Application styles
│   └── public/
│       └── index.html       # HTML template
└── README.md                # This file
```

## 🗄️ Database Schema

### Tables

- **molecules**: Core molecular data (SMILES, properties, dataset info)
- **pfasgroups_results**: PFASGroups algorithm results (includes matched path types)
- **atlas_results**: PFAS-Atlas algorithm results  
- **manual_reviews**: Human reviewer validations

Note: Timing benchmark data is analyzed separately using Python scripts that read JSON files directly.

### Key Fields

- **Review Status**: `pfasgroups_correct`, `atlas_correct` (boolean or null)
- **Reviewer Notes**: Free-text annotations
- **Correct Classifications**: Manual ground truth data
- **Review Timestamps**: Track review progress over time

## 🔬 Review Interface

### Molecule Review Panel

Each molecule displays:
- **2D Structure**: Interactive RDKit.js visualization
- **Molecular Properties**: SMILES, molecular weight, chain length
- **Algorithm Results**: PFASGroups detected groups, PFAS-Atlas classifications
- **Execution Times**: Performance metrics for both algorithms
- **Review Buttons**: ✅ Correct, ❌ Incorrect, 🤷 Unclear for each algorithm

### Filtering Options

- **Dataset Type**: Filter by OECD, Enhanced, Timing, etc.
- **Review Status**: All, Reviewed, Unreviewed
- **Text Search**: Search by SMILES string or group names
- **Pagination**: Configurable page size (default: 10 molecules)

## 📊 Dashboard Features

### Overview Statistics
- Total molecules per dataset
- Review progress percentages
- Overall accuracy metrics
- Algorithm comparison charts

### Dataset Breakdown
- Per-dataset review completion
- Progress bars and status badges
- Drill-down capabilities

## 📈 Accuracy Reports

### Metrics Computed
- **Overall Accuracy**: Percentage of correct classifications
- **Per-Dataset Accuracy**: Breakdown by dataset type
- **Algorithm Comparison**: Side-by-side accuracy comparison
- **Confidence Intervals**: Statistical significance measures

### Export Capabilities
- **JSON Export**: Complete review data with metadata
- **CSV Export**: Tabular format for external analysis
- **Real-time Updates**: Metrics update as reviews are submitted

## 🔧 API Endpoints

### Core Endpoints
- `GET /api/molecules` - Paginated molecule list with filtering
- `POST /api/review` - Submit manual review
- `GET /api/stats` - Dashboard statistics
- `GET /api/accuracy` - Accuracy metrics
- `GET /api/export/reviews` - Export review data

### Query Parameters
- `page`, `limit` - Pagination
- `dataset` - Filter by dataset type
- `reviewStatus` - Filter by review completion
- `search` - Text search

## 🛠️ Development

### Adding New Datasets

1. **Update Database Schema**: Add new dataset types to the enum
2. **Import Script**: Modify `import-benchmark-data.js` for new data format
3. **Frontend Filters**: Add new dataset options to filter dropdowns
4. **Validation Logic**: Update review validation for dataset-specific requirements

### Customizing Review Interface

- **Review Fields**: Modify the `manual_reviews` table schema
- **UI Components**: Edit React components in `client/src/components/`
- **Validation Logic**: Update server-side validation in `server.js`

### Performance Optimization

- **Database Indexing**: Indexes already created on key columns
- **Pagination**: Configurable page sizes to manage memory
- **Lazy Loading**: Components load data on demand
- **Caching**: Consider Redis for session management in production

## 🔄 Integration with Benchmark Pipeline

### Automatic Data Import

The `run_all_benchmarks.sh` script now includes database integration:

```bash
# After running benchmarks, data is automatically imported
./run_all_benchmarks.sh
# Database is updated with latest results
```

### Manual Data Import

```bash
# Import specific JSON files
node scripts/import-benchmark-data.js

# Or use the helper script
./import-latest-data.sh
```

## 📝 Review Workflow

### Recommended Process

1. **Start with OECD Dataset**: Well-characterized reference molecules
2. **Filter by Unreviewed**: Focus on molecules needing validation
3. **Review in Batches**: Process 20-50 molecules per session
4. **Add Detailed Notes**: Document reasoning for classifications
5. **Export Progress**: Regular backups of review data

### Quality Control

- **Inter-reviewer Agreement**: Multiple reviewers for critical molecules
- **Review Timestamps**: Track review patterns and consistency
- **Audit Trail**: Complete history of review changes
- **Statistical Validation**: Confidence intervals and significance tests

## 🐛 Troubleshooting

### Common Issues

1. **RDKit.js Loading**: Molecule structures not displaying
   - Check browser console for RDKit.js errors
   - Fallback to SMILES string display

2. **Database Connection**: SQLite file permissions
   - Ensure write permissions for `database/pfas_benchmark.db`
   - Check file path in database configuration

3. **Memory Usage**: Large datasets causing slowdown
   - Reduce pagination limit
   - Filter to smaller dataset subsets
   - Monitor browser memory usage

4. **API Timeouts**: Slow database queries
   - Check database indexes
   - Optimize query filters
   - Consider database vacuum/analyze

### Debug Mode

```bash
# Enable detailed logging
DEBUG=* npm run dev

# Database debugging
sqlite3 database/pfas_benchmark.db ".schema"
```

## 📋 Production Deployment

### Requirements

- **Server**: Linux/macOS with Node.js
- **Memory**: 2GB RAM minimum (4GB recommended)
- **Storage**: 10GB for database and exports
- **Network**: HTTPS recommended for production

### Security Considerations

- **Input Validation**: All API inputs sanitized
- **SQL Injection**: Parameterized queries used throughout
- **File Access**: Restricted to designated directories
- **Session Management**: Stateless design with JWT tokens (if implementing auth)

### Performance Monitoring

```bash
# Monitor database size
du -h database/pfas_benchmark.db

# Check review completion rates
sqlite3 database/pfas_benchmark.db "SELECT dataset_type, COUNT(*) as total, SUM(CASE WHEN mr.id IS NOT NULL THEN 1 ELSE 0 END) as reviewed FROM molecules m LEFT JOIN manual_reviews mr ON m.id = mr.molecule_id GROUP BY dataset_type;"
```

## 🤝 Contributing

### Adding Features

1. **Database Changes**: Update schema in `database/database.js`
2. **API Endpoints**: Add routes in `server.js`
3. **Frontend Components**: Create React components in `client/src/components/`
4. **Documentation**: Update this README

### Code Style

- **Backend**: JavaScript ES6+ with async/await
- **Frontend**: React functional components with hooks
- **Database**: SQLite with prepared statements
- **Styling**: Bootstrap 5 with custom CSS

---

## 📞 Support

For issues, questions, or feature requests, please refer to the main PFASGroups repository or contact the development team.

**Version**: 1.0.0  
**Last Updated**: December 19, 2025  
**Maintainer**: Luc Miaz