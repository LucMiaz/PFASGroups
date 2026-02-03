const Database = require('./database/database');
const path = require('path');

async function quickReimport() {
    const db = new Database();
    await db.waitForReady();
    
    console.log('🗑️  Clearing highly_branched data...');
    
    // Delete only highly_branched data
    await db.run("DELETE FROM manual_reviews WHERE molecule_id IN (SELECT id FROM molecules WHERE dataset_type = 'highly_branched')");
    await db.run("DELETE FROM pfasgroups_results WHERE molecule_id IN (SELECT id FROM molecules WHERE dataset_type = 'highly_branched')");
    await db.run("DELETE FROM atlas_results WHERE molecule_id IN (SELECT id FROM molecules WHERE dataset_type = 'highly_branched')");
    await db.run("DELETE FROM molecules WHERE dataset_type = 'highly_branched'");
    
    console.log('✅ Cleared highly_branched data');
    console.log('📥 Re-importing highly_branched data...\n');
    
    // Now import just the highly branched data
    const DataImporter = require('./scripts/import-benchmark-data');
    const importer = new DataImporter();
    
    const dataPath = path.join(__dirname, '../data');
    const fs = require('fs');
    const files = fs.readdirSync(dataPath);
    const hbFile = files.find(f => f.includes('highly_branched') && f.endsWith('.json'));
    
    if (hbFile) {
        console.log(`Found file: ${hbFile}`);
        const filePath = path.join(dataPath, hbFile);
        await importer.importBenchmarkFile(filePath, 'highly_branched', hbFile);
        console.log('✅ Successfully re-imported highly_branched data');
        
        // Force save and wait for it to complete
        console.log('💾 Saving database...');
        await new Promise(resolve => {
            db.save();
            // Wait for save timeout + some buffer
            setTimeout(resolve, 2000);
        });
        console.log('✅ Database saved');
    } else {
        console.log('❌ No highly_branched file found');
    }
    
    process.exit(0);
}

quickReimport().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
