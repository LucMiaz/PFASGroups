const initSqlJs = require('sql.js');
const fs = require('fs');
const path = require('path');

class DatabaseWrapper {
    constructor() {
        this.dbPath = path.join(__dirname, 'pfas_benchmark.db');
        this.db = null;
        this.isReady = false;
        this.initPromise = this.initialize();
        this.saveTimeout = null;
        this.isSaving = false;
    }

    async initialize() {
        try {
            const SQL = await initSqlJs();
            
            // Check if database file exists
            if (fs.existsSync(this.dbPath)) {
                const buffer = fs.readFileSync(this.dbPath);
                this.db = new SQL.Database(buffer);
                console.log('Connected to existing SQLite database');
            } else {
                this.db = new SQL.Database();
                console.log('Created new SQLite database');
            }
            
            this.initTables();
            this.isReady = true;
        } catch (err) {
            console.error('Error opening database:', err.message);
            throw err;
        }
    }

    async waitForReady() {
        await this.initPromise;
    }

    initTables() {
        // Main molecules table
        this.db.exec(`
            CREATE TABLE IF NOT EXISTS molecules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                smiles TEXT NOT NULL,
                molecular_formula TEXT,
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

        // PFASGroups results table (default flavor: bycomponent=False)
        this.db.exec(`
            CREATE TABLE IF NOT EXISTS pfasgroups_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule_id INTEGER,
                detected_groups TEXT,
                detected_definitions TEXT,
                success BOOLEAN,
                error_message TEXT,
                execution_time REAL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (molecule_id) REFERENCES molecules (id)
            )
        `);

        // PFASGroups bycomponent results table (bycomponent=True flavor)
        this.db.exec(`
            CREATE TABLE IF NOT EXISTS pfasgroups_results_bycomponent (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                molecule_id INTEGER,
                detected_groups TEXT,
                detected_definitions TEXT,
                success BOOLEAN,
                error_message TEXT,
                execution_time REAL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (molecule_id) REFERENCES molecules (id)
            )
        `);

        // PFAS-Atlas results table
        this.db.exec(`
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
        this.db.exec(`
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
        this.db.exec(`
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
        this.db.exec(`CREATE INDEX IF NOT EXISTS idx_molecules_dataset ON molecules (dataset_type)`);
        this.db.exec(`CREATE INDEX IF NOT EXISTS idx_molecules_smiles ON molecules (smiles)`);
        this.db.exec(`CREATE INDEX IF NOT EXISTS idx_reviews_molecule ON manual_reviews (molecule_id)`);
        this.db.exec(`CREATE INDEX IF NOT EXISTS idx_pfasgroups_molecule ON pfasgroups_results (molecule_id)`);
        this.db.exec(`CREATE INDEX IF NOT EXISTS idx_pfasgroups_bycomp_molecule ON pfasgroups_results_bycomponent (molecule_id)`);
        this.db.exec(`CREATE INDEX IF NOT EXISTS idx_atlas_molecule ON atlas_results (molecule_id)`);
        
        // Add molecular_formula column if it doesn't exist (migration)
        try {
            this.db.exec(`ALTER TABLE molecules ADD COLUMN molecular_formula TEXT`);
        } catch (e) {
            // Column already exists, ignore error
        }
        
        // Add detected_definitions column to pfasgroups_results if it doesn't exist (migration)
        try {
            this.db.exec(`ALTER TABLE pfasgroups_results ADD COLUMN detected_definitions TEXT`);
        } catch (e) {
            // Column already exists, ignore error
        }
        
        // Add detected_definitions column to pfasgroups_results_bycomponent if it doesn't exist (migration)
        try {
            this.db.exec(`ALTER TABLE pfasgroups_results_bycomponent ADD COLUMN detected_definitions TEXT`);
        } catch (e) {
            // Column already exists, ignore error
        }
        
        // Save after creating tables
        this.save();
    }

    save() {
        if (!this.db) return;
        
        // Debounce saves to avoid too frequent writes
        if (this.saveTimeout) {
            clearTimeout(this.saveTimeout);
        }
        
        this.saveTimeout = setTimeout(() => {
            if (this.isSaving) {
                console.log('Save already in progress, skipping...');
                return;
            }
            
            this.isSaving = true;
            
            try {
                const data = this.db.export();
                const buffer = Buffer.from(data);
                
                // Write to a temporary file first
                const tempPath = this.dbPath + '.tmp';
                fs.writeFileSync(tempPath, buffer);
                
                // Then rename it (atomic operation on most systems)
                if (fs.existsSync(this.dbPath)) {
                    // On Windows, we need to remove the target first
                    const backupPath = this.dbPath + '.backup';
                    if (fs.existsSync(backupPath)) {
                        try {
                            fs.unlinkSync(backupPath);
                        } catch (e) {
                            // Ignore if backup can't be removed
                        }
                    }
                    try {
                        fs.renameSync(this.dbPath, backupPath);
                    } catch (e) {
                        console.warn('Could not create backup:', e.message);
                    }
                }
                fs.renameSync(tempPath, this.dbPath);
                
                // Clean up old backup
                const backupPath = this.dbPath + '.backup';
                if (fs.existsSync(backupPath)) {
                    try {
                        fs.unlinkSync(backupPath);
                    } catch (e) {
                        // Ignore if backup can't be removed
                    }
                }
            } catch (err) {
                console.error('Error saving database:', err.message);
                // Try to recover from temp file if it exists
                const tempPath = this.dbPath + '.tmp';
                if (fs.existsSync(tempPath)) {
                    try {
                        fs.unlinkSync(tempPath);
                    } catch (e) {
                        // Ignore
                    }
                }
            } finally {
                this.isSaving = false;
            }
        }, 100); // Debounce by 100ms
    }

    // Helper method to run queries
    async run(sql, params = []) {
        await this.waitForReady();
        try {
            this.db.run(sql, params);
            const lastId = this.db.exec("SELECT last_insert_rowid() as id")[0]?.values[0]?.[0] || 0;
            
            // Save asynchronously to avoid blocking
            setImmediate(() => this.save());
            
            return { id: lastId, changes: this.db.getRowsModified() };
        } catch (err) {
            throw err;
        }
    }

    // Helper method to get single row
    async get(sql, params = []) {
        await this.waitForReady();
        try {
            const stmt = this.db.prepare(sql);
            stmt.bind(params);
            let row = null;
            if (stmt.step()) {
                const columns = stmt.getColumnNames();
                const values = stmt.get();
                row = {};
                columns.forEach((col, idx) => {
                    row[col] = values[idx];
                });
            }
            stmt.free();
            return row;
        } catch (err) {
            throw err;
        }
    }

    // Helper method to get multiple rows
    async all(sql, params = []) {
        await this.waitForReady();
        try {
            const stmt = this.db.prepare(sql);
            stmt.bind(params);
            const rows = [];
            while (stmt.step()) {
                const columns = stmt.getColumnNames();
                const values = stmt.get();
                const row = {};
                columns.forEach((col, idx) => {
                    row[col] = values[idx];
                });
                rows.push(row);
            }
            stmt.free();
            return rows;
        } catch (err) {
            throw err;
        }
    }

    // Close database connection
    async close() {
        await this.waitForReady();
        try {
            // Clear any pending save timeout
            if (this.saveTimeout) {
                clearTimeout(this.saveTimeout);
            }
            
            // Wait for any ongoing save to complete
            while (this.isSaving) {
                await new Promise(resolve => setTimeout(resolve, 50));
            }
            
            // Do final save synchronously
            if (this.db) {
                try {
                    const data = this.db.export();
                    const buffer = Buffer.from(data);
                    fs.writeFileSync(this.dbPath, buffer);
                    console.log('Database saved on close');
                } catch (err) {
                    console.error('Error saving database on close:', err.message);
                }
                this.db.close();
            }
        } catch (err) {
            console.error('Error closing database:', err.message);
            throw err;
        }
    }
}

module.exports = DatabaseWrapper;