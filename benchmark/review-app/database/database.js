const sqlite3 = require('sqlite3').verbose();
const path = require('path');

class Database {
    constructor() {
        const dbPath = path.join(__dirname, 'pfas_benchmark.db');
        this.db = new sqlite3.Database(dbPath, (err) => {
            if (err) {
                console.error('Error opening database:', err.message);
            } else {
                console.log('Connected to SQLite database');
                this.initTables();
            }
        });
    }

    initTables() {
        // Main molecules table
        this.db.run(`
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT NOT NULL,
                molecular_weight REAL,
                num_atoms INTEGER,
                num_bonds INTEGER,
                chain_length INTEGER,
                target_groups TEXT,
                generation_type TEXT,
                group_id INTEGER,
                group_name TEXT,
                dataset_type TEXT NOT NULL,
                benchmark_date TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )
        `);

        // PFASGroups results table
        this.db.run(`
            CREATE TABLE IF NOT EXISTS pfasgroups_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule_id INTEGER,
                detected_groups TEXT,
                success BOOLEAN,
                error_message TEXT,
                execution_time REAL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (molecule_id) REFERENCES molecules (id)
            )
        `);

        // PFAS-Atlas results table
        this.db.run(`
            CREATE TABLE IF NOT EXISTS atlas_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule_id INTEGER,
                first_class TEXT,
                second_class TEXT,
                success BOOLEAN,
                error_message TEXT,
                execution_time REAL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (molecule_id) REFERENCES molecules (id)
            )
        `);

        // Manual reviews table
        this.db.run(`
            CREATE TABLE IF NOT EXISTS manual_reviews (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule_id INTEGER,
                pfasgroups_correct BOOLEAN,
                atlas_correct BOOLEAN,
                reviewer_notes TEXT,
                reviewer_name TEXT,
                review_date DATETIME DEFAULT CURRENT_TIMESTAMP,
                is_pfas BOOLEAN,
                correct_groups TEXT,
                correct_classification TEXT,
                FOREIGN KEY (molecule_id) REFERENCES molecules (id)
            )
        `);

        // Timing benchmarks table
        this.db.run(`
            CREATE TABLE IF NOT EXISTS timing_benchmarks (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule_id INTEGER,
                iterations INTEGER,
                pfasgroups_time_avg REAL,
                pfasgroups_time_std REAL,
                pfasgroups_time_min REAL,
                pfasgroups_time_max REAL,
                atlas_time_avg REAL,
                atlas_time_std REAL,
                atlas_time_min REAL,
                atlas_time_max REAL,
                pfasgroups_success_rate REAL,
                atlas_success_rate REAL,
                benchmark_date TEXT,
                FOREIGN KEY (molecule_id) REFERENCES molecules (id)
            )
        `);

        // Create indexes for better performance
        this.db.run(`CREATE INDEX IF NOT EXISTS idx_molecules_dataset ON molecules (dataset_type)`);
        this.db.run(`CREATE INDEX IF NOT EXISTS idx_molecules_smiles ON molecules (smiles)`);
        this.db.run(`CREATE INDEX IF NOT EXISTS idx_reviews_molecule ON manual_reviews (molecule_id)`);
        this.db.run(`CREATE INDEX IF NOT EXISTS idx_pfasgroups_molecule ON pfasgroups_results (molecule_id)`);
        this.db.run(`CREATE INDEX IF NOT EXISTS idx_atlas_molecule ON atlas_results (molecule_id)`);
    }

    // Helper method to run queries with promises
    run(sql, params = []) {
        return new Promise((resolve, reject) => {
            this.db.run(sql, params, function(err) {
                if (err) {
                    reject(err);
                } else {
                    resolve({ id: this.lastID, changes: this.changes });
                }
            });
        });
    }

    // Helper method to get single row
    get(sql, params = []) {
        return new Promise((resolve, reject) => {
            this.db.get(sql, params, (err, row) => {
                if (err) {
                    reject(err);
                } else {
                    resolve(row);
                }
            });
        });
    }

    // Helper method to get multiple rows
    all(sql, params = []) {
        return new Promise((resolve, reject) => {
            this.db.all(sql, params, (err, rows) => {
                if (err) {
                    reject(err);
                } else {
                    resolve(rows);
                }
            });
        });
    }

    // Close database connection
    close() {
        return new Promise((resolve, reject) => {
            this.db.close((err) => {
                if (err) {
                    reject(err);
                } else {
                    resolve();
                }
            });
        });
    }
}

module.exports = Database;