const Database = require('../database/database');
const fs = require('fs-extra');
const path = require('path');

class DataImporter {
    constructor() {
        this.db = new Database();
    }

    async importBenchmarkData(dataDirectory = '../data') {
        const dataPath = path.join(__dirname, dataDirectory);
        
        try {
            // Get all JSON files in the data directory
            const files = await fs.readdir(dataPath);
            const jsonFiles = files.filter(file => file.endsWith('.json'));
            
            console.log(`Found ${jsonFiles.length} JSON files to import`);
            
            for (const file of jsonFiles) {
                const filePath = path.join(dataPath, file);
                console.log(`Importing ${file}...`);
                
                // Determine dataset type from filename
                const datasetType = this.determineDatasetType(file);
                
                if (datasetType === 'timing') {
                    await this.importTimingData(filePath, file);
                } else {
                    await this.importBenchmarkFile(filePath, datasetType, file);
                }
            }
            
            console.log('Data import completed successfully');
        } catch (error) {
            console.error('Error importing data:', error);
        }
    }

    determineDatasetType(filename) {
        if (filename.includes('timing')) return 'timing';
        if (filename.includes('oecd')) return 'oecd';
        if (filename.includes('enhanced')) return 'enhanced';
        if (filename.includes('non_fluorinated')) return 'non_fluorinated';
        if (filename.includes('complex_branched')) return 'complex_branched';
        return 'unknown';
    }

    async importBenchmarkFile(filePath, datasetType, filename) {
        const data = JSON.parse(await fs.readFile(filePath, 'utf8'));
        const benchmarkDate = this.extractDateFromFilename(filename);
        
        for (const record of data) {
            try {
                // Insert molecule data
                const moleculeData = record.molecule_data || record;
                const moleculeId = await this.insertMolecule(moleculeData, datasetType, benchmarkDate);
                
                // Insert PFASGroups result if exists
                if (record.pfasgroups_result) {
                    await this.insertPFASGroupsResult(moleculeId, record.pfasgroups_result);
                }
                
                // Insert Atlas result if exists
                if (record.atlas_result) {
                    await this.insertAtlasResult(moleculeId, record.atlas_result);
                }
            } catch (error) {
                console.error(`Error importing record in ${filename}:`, error);
                // Continue with next record
            }
        }
        
        console.log(`✓ Imported ${data.length} records from ${filename}`);
    }

    async importTimingData(filePath, filename) {
        const data = JSON.parse(await fs.readFile(filePath, 'utf8'));
        const benchmarkDate = this.extractDateFromFilename(filename);
        
        for (const record of data) {
            try {
                // Insert molecule data
                const moleculeId = await this.insertTimingMolecule(record, benchmarkDate);
                
                // Insert timing benchmark data
                await this.insertTimingBenchmark(moleculeId, record, benchmarkDate);
            } catch (error) {
                console.error(`Error importing timing record in ${filename}:`, error);
            }
        }
        
        console.log(`✓ Imported ${data.length} timing records from ${filename}`);
    }

    async insertMolecule(moleculeData, datasetType, benchmarkDate) {
        const result = await this.db.run(`
            INSERT INTO molecules (
                smiles, molecular_weight, num_atoms, num_bonds, 
                chain_length, target_groups, generation_type, 
                group_id, group_name, dataset_type, benchmark_date
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `, [
            moleculeData.smiles,
            moleculeData.molecular_weight || null,
            moleculeData.num_atoms || null,
            moleculeData.num_bonds || null,
            moleculeData.chain_length || null,
            JSON.stringify(moleculeData.target_groups || []),
            moleculeData.generation_type || null,
            moleculeData.group_id || null,
            moleculeData.group_name || null,
            datasetType,
            benchmarkDate
        ]);
        
        return result.id;
    }

    async insertTimingMolecule(record, benchmarkDate) {
        const result = await this.db.run(`
            INSERT INTO molecules (
                smiles, molecular_weight, num_atoms, num_bonds, 
                chain_length, dataset_type, benchmark_date
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
        `, [
            record.smiles,
            record.molecular_weight || null,
            record.num_atoms || null,
            record.num_bonds || null,
            record.chain_length || null,
            'timing',
            benchmarkDate
        ]);
        
        return result.id;
    }

    async insertPFASGroupsResult(moleculeId, result) {
        await this.db.run(`
            INSERT INTO pfasgroups_results (
                molecule_id, detected_groups, success, error_message, execution_time
            ) VALUES (?, ?, ?, ?, ?)
        `, [
            moleculeId,
            JSON.stringify(result.detected_groups || []),
            result.success || false,
            result.error || null,
            result.execution_time || null
        ]);
    }

    async insertAtlasResult(moleculeId, result) {
        await this.db.run(`
            INSERT INTO atlas_results (
                molecule_id, first_class, second_class, success, error_message, execution_time
            ) VALUES (?, ?, ?, ?, ?, ?)
        `, [
            moleculeId,
            result.first_class || null,
            result.second_class || null,
            result.success || false,
            result.error || null,
            result.execution_time || null
        ]);
    }

    async insertTimingBenchmark(moleculeId, record, benchmarkDate) {
        await this.db.run(`
            INSERT INTO timing_benchmarks (
                molecule_id, iterations, 
                pfasgroups_time_avg, pfasgroups_time_std, pfasgroups_time_min, pfasgroups_time_max,
                atlas_time_avg, atlas_time_std, atlas_time_min, atlas_time_max,
                pfasgroups_success_rate, atlas_success_rate, benchmark_date
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `, [
            moleculeId,
            record.iterations || 1,
            record.pfasgroups_time_avg || record.pfasgroups_time || null,
            record.pfasgroups_time_std || null,
            record.pfasgroups_time_min || null,
            record.pfasgroups_time_max || null,
            record.atlas_time_avg || record.atlas_time || null,
            record.atlas_time_std || null,
            record.atlas_time_min || null,
            record.atlas_time_max || null,
            record.pfasgroups_success_rate || null,
            record.atlas_success_rate || null,
            benchmarkDate
        ]);
    }

    extractDateFromFilename(filename) {
        // Extract date from filename pattern like "pfas_oecd_benchmark_20251219_085844.json"
        const match = filename.match(/_(\d{8})_\d{6}\.json$/);
        if (match) {
            const dateStr = match[1];
            return `${dateStr.substr(0,4)}-${dateStr.substr(4,2)}-${dateStr.substr(6,2)}`;
        }
        return new Date().toISOString().split('T')[0];
    }

    async getImportStats() {
        const stats = {};
        
        // Count by dataset type
        const datasets = await this.db.all(`
            SELECT dataset_type, COUNT(*) as count 
            FROM molecules 
            GROUP BY dataset_type
        `);
        
        stats.molecules = datasets;
        
        // Count reviews
        const reviewCount = await this.db.get(`
            SELECT COUNT(*) as count 
            FROM manual_reviews
        `);
        
        stats.reviews = reviewCount.count;
        
        return stats;
    }

    async close() {
        await this.db.close();
    }
}

// Run import if called directly
if (require.main === module) {
    const importer = new DataImporter();
    
    importer.importBenchmarkData()
        .then(async () => {
            const stats = await importer.getImportStats();
            console.log('\\nImport Statistics:');
            console.log(JSON.stringify(stats, null, 2));
            await importer.close();
        })
        .catch(async (error) => {
            console.error('Import failed:', error);
            await importer.close();
            process.exit(1);
        });
}

module.exports = DataImporter;