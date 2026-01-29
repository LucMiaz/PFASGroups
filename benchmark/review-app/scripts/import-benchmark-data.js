const Database = require('../database/database');
const fs = require('fs-extra');
const path = require('path');
const { execSync } = require('child_process');

/**
 * DataImporter - Imports benchmark data into the database
 * 
 * DEDUPLICATION STRATEGY:
 * ----------------------
 * This importer automatically removes duplicates during data import:
 * 
 * 1. Benchmark Data (OECD, Enhanced, etc.):
 *    - Duplicates are identified by matching BOTH:
 *      a) SMILES string (molecular structure)
 *      b) All classification results (PFASGroups, Atlas)
 *    - Keeps the first occurrence of each unique combination
 *    - Example: Same molecule with different classifications = kept as separate records
 * 
 * 2. Timing Data:
 *    - Duplicates are identified by SMILES only
 *    - Keeps the first occurrence of each unique SMILES
 *    - Timing measurements are averaged across the dataset
 * 
 * This ensures data quality while preserving legitimate variations in results.
 */
class DataImporter {
    constructor() {
        this.db = new Database();
    }

    async importBenchmarkData(dataDirectory = '../../data') {
        const dataPath = path.isAbsolute(dataDirectory) 
            ? dataDirectory 
            : path.join(__dirname, dataDirectory);
        
        console.log(`Looking for data in: ${dataPath}`);
        
        // Wait for database to be ready
        await this.db.waitForReady();
        
        try {
            // Get all JSON files in the data directory
            const files = await fs.readdir(dataPath);
            const jsonFiles = files.filter(file => file.endsWith('.json'));
            
            console.log(`Found ${jsonFiles.length} JSON files to import`);
            
            for (const file of jsonFiles) {
                const filePath = path.join(dataPath, file);
                
                // Skip analysis/summary/test files
                if (file.includes('validation') || file.includes('analysis') || file.includes('summary') || 
                    file.includes('test_api') || !file.startsWith('pfas_')) {
                    console.log(`Skipping non-benchmark file: ${file}`);
                    continue;
                }
                
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
        if (filename.includes('highly_branched')) return 'highly_branched';
        if (filename.includes('definitions')) return 'definitions';
        return 'unknown';
    }

    async importBenchmarkFile(filePath, datasetType, filename) {
        let data = JSON.parse(await fs.readFile(filePath, 'utf8'));
        const benchmarkDate = this.extractDateFromFilename(filename);
        
        // Handle different data structures
        
        // Highly branched benchmark has: { metadata: {...}, summary: {...}, details: [...] }
        if (data.metadata && data.details && Array.isArray(data.details)) {
            const allMolecules = [];
            
            // Transform details array into molecule records
            for (const test of data.details) {
                if (test.smiles) {
                    allMolecules.push({
                        molecule_data: {
                            smiles: test.smiles,
                            group_id: test.group_id,
                            group_name: test.group_name
                        },
                        pfasgroups_result: {
                            detected_groups: test.group_id ? [test.group_id] : [],
                            success: test.passed === true,
                            error: test.passed === false ? 'Test failed' : null,
                            execution_time: null
                        },
                        atlas_result: {
                            first_class: null,
                            second_class: null,
                            success: false,
                            error: false,
                            execution_time: null
                        }
                    });
                }
            }
            data = allMolecules;
        }
        // Definitions benchmark data has a nested structure: { metadata: {...}, benchmarks: {...} }
        else if (data.metadata && data.benchmarks) {
            const allMolecules = [];
            
            // Flatten all molecules from all categories and subcategories
            for (const [categoryName, categoryData] of Object.entries(data.benchmarks)) {
                if (categoryData.subcategories) {
                    for (const [subcategoryName, molecules] of Object.entries(categoryData.subcategories)) {
                        if (Array.isArray(molecules)) {
                            for (const mol of molecules) {
                                if (mol.smiles) {
                                    allMolecules.push({
                                        molecule_data: {
                                            smiles: mol.smiles,
                                            molecular_weight: mol.molecular_weight,
                                            num_atoms: mol.num_atoms
                                        },
                                        pfasgroups_result: {
                                            detected_groups: mol.pfasgroups_groups || [],
                                            detected_definitions: mol.pfasgroups_definitions || [],
                                            success: mol.pfasgroups_detected !== false,
                                            error: null,
                                            execution_time: mol.pfasgroups_execution_time
                                        },
                                        atlas_result: {
                                            first_class: mol.atlas_first_class,
                                            second_class: mol.atlas_second_class,
                                            success: mol.atlas_detected !== false,
                                            error: null,
                                            execution_time: mol.atlas_execution_time
                                        }
                                    });
                                }
                            }
                        }
                    }
                }
            }
            data = allMolecules;
        }
        // Complex branched data has a wrapper object with "molecules" array
        else if (Array.isArray(data) && data.length > 0 && data[0].molecules) {
            // Flatten the molecules from all test cases
            const allMolecules = [];
            for (const testCase of data) {
                if (testCase.molecules && Array.isArray(testCase.molecules)) {
                    // Transform the molecules to the expected format
                    for (const mol of testCase.molecules) {
                        allMolecules.push({
                            molecule_data: {
                                smiles: mol.smiles,
                                molecular_weight: mol.molecular_weight,
                                num_atoms: mol.num_atoms
                            },
                            pfasgroups_result: {
                                detected_groups: mol.pfasgroups_groups || [],
                                success: mol.pfasgroups_detected !== false,
                                error: null,
                                execution_time: mol.pfasgroups_execution_time
                            },
                            atlas_result: {
                                first_class: mol.atlas_first_class,
                                second_class: mol.atlas_second_class,
                                success: mol.atlas_detected !== false,
                                error: null,
                                execution_time: mol.atlas_execution_time
                            }
                        });
                    }
                }
            }
            data = allMolecules;
        }
        
        // Deduplicate records based on SMILES and classification results
        const originalCount = data.length;
        data = this.deduplicateRecords(data);
        const duplicatesRemoved = originalCount - data.length;
        if (duplicatesRemoved > 0) {
            console.log(`  Removed ${duplicatesRemoved} duplicate(s) from ${filename}`);
        }
        
        for (const record of data) {
            try {
                // Insert molecule data
                const moleculeData = record.molecule_data || record;
                
                // Calculate molecular formula if not present (skip for now - too slow)
                // Can be calculated on-demand or batch-processed later
                if (!moleculeData.molecular_formula && moleculeData.smiles) {
                    moleculeData.molecular_formula = null; // Placeholder
                }
                
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
        
        // DON'T deduplicate timing records - we want to keep triplicates for statistical analysis
        console.log(`  Keeping all ${data.length} timing records (triplicates preserved for analysis)`);
        
        for (const record of data) {
            try {
                // Calculate molecular formula if not present (skip for now - too slow)
                // Can be calculated on-demand or batch-processed later
                if (!record.molecular_formula && record.smiles) {
                    record.molecular_formula = null; // Placeholder
                }
                
                // Insert molecule data
                const moleculeId = await this.insertTimingMolecule(record, benchmarkDate);
                
                // Note: Timing benchmark data is analyzed separately via Python scripts,
                // not stored in database
            } catch (error) {
                console.error(`Error importing timing record in ${filename}:`, error);
            }
        }
        
        console.log(`✓ Imported ${data.length} timing records from ${filename}`);
    }

    /**
     * Deduplicate records based on SMILES and classification results
     * Keeps the first occurrence of each unique combination
     */
    deduplicateRecords(records) {
        const seen = new Map();
        const deduplicated = [];
        
        for (const record of records) {
            const moleculeData = record.molecule_data || record;
            const smiles = moleculeData.smiles;
            
            if (!smiles) {
                deduplicated.push(record);
                continue;
            }
            
            // Create a unique key based on SMILES and classification results
            const pfasgroupsGroups = JSON.stringify((record.pfasgroups_result?.detected_groups || []).sort());
            const atlasFirstClass = record.atlas_result?.first_class || '';
            const atlasSecondClass = record.atlas_result?.second_class || '';
            
            const key = `${smiles}|${pfasgroupsGroups}|${atlasFirstClass}|${atlasSecondClass}`;
            
            if (!seen.has(key)) {
                seen.set(key, true);
                deduplicated.push(record);
            }
        }
        
        return deduplicated;
    }

    /**
     * Deduplicate timing records based on SMILES only
     * For timing data, we keep the first occurrence of each SMILES
     */
    deduplicateTimingRecords(records) {
        const seen = new Set();
        const deduplicated = [];
        
        for (const record of records) {
            const smiles = record.smiles;
            
            if (!smiles) {
                deduplicated.push(record);
                continue;
            }
            
            if (!seen.has(smiles)) {
                seen.add(smiles);
                deduplicated.push(record);
            }
        }
        
        return deduplicated;
    }

    async insertMolecule(moleculeData, datasetType, benchmarkDate) {
        const result = await this.db.run(`
            INSERT INTO molecules (
                smiles, molecular_formula, molecular_weight, num_atoms, num_bonds, 
                chain_length, target_groups, generation_type, 
                group_id, group_name, dataset_type, benchmark_date
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        `, [
            moleculeData.smiles,
            moleculeData.molecular_formula || null,
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
                smiles, molecular_formula, molecular_weight, num_atoms, num_bonds, 
                chain_length, dataset_type, benchmark_date
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        `, [
            record.smiles,
            record.molecular_formula || null,
            record.molecular_weight || null,
            record.num_atoms || null,
            record.num_bonds || null,
            record.chain_length || null,
            'timing',
            benchmarkDate
        ]);
        
        return result.id;
    }

    /**
     * Extract matched path types from matches data
     * Returns a map of groupId -> pathTypes (e.g., "Perfluoroalkyl", "Polyfluoroalkyl", or "Perfluoroalkyl,Polyfluoroalkyl")
     */
    extractMatchedPathTypes(matches, detectedGroups) {
        const pathTypesMap = {};
        
        if (!matches || !Array.isArray(matches)) {
            return pathTypesMap;
        }
        
        // Process only group matches (not definition matches)
        for (const match of matches) {
            if (match.type === 'group' && match.id && match.components_types) {
                // Store all component types joined by comma (can be multiple per match)
                if (Array.isArray(match.components_types) && match.components_types.length > 0) {
                    pathTypesMap[match.id] = match.components_types.join(',');
                }
            }
        }
        
        return pathTypesMap;
    }

    async insertPFASGroupsResult(moleculeId, result) {
        // Extract matched path types from matches data
        const matchedPathTypes = this.extractMatchedPathTypes(result.matches, result.detected_groups);
        
        await this.db.run(`
            INSERT INTO pfasgroups_results (
                molecule_id, detected_groups, detected_definitions, matched_path_types, success, error_message, execution_time
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
        `, [
            moleculeId,
            JSON.stringify(result.detected_groups || []),
            JSON.stringify(result.detected_definitions || []),
            JSON.stringify(matchedPathTypes),
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

    /**
     * Calculate molecular formula from SMILES using RDKit
     */
    calculateMolecularFormula(smiles) {
        try {
            const pythonCode = `
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles('${smiles.replace(/'/g, "\\'")}')
if mol:
    print(Descriptors.rdMolDescriptors.CalcMolFormula(mol))
else:
    print('')
`;
            
            const result = execSync(`python -c "${pythonCode.replace(/"/g, '\\"')}"`, {
                encoding: 'utf8',
                timeout: 5000,
                windowsHide: true
            }).trim();
            
            return result || null;
        } catch (error) {
            return null;
        }
    }

    async close() {
        await this.db.close();
    }
}

// Run import if called directly
if (require.main === module) {
    const importer = new DataImporter();
    
    // Look for data in the benchmark/data directory (one level up from review-app)
    const dataDir = path.join(__dirname, '../../data');
    
    importer.importBenchmarkData(dataDir)
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
