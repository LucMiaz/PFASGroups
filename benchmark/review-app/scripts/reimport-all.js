const Database = require('../database/database');
const fs = require('fs');
const path = require('path');

async function reimportAll() {
    const db = new Database();
    await db.waitForReady();
    
    console.log('🗑️  Clearing all existing data...');
    
    // Clear all tables
    await db.run('DELETE FROM manual_reviews');
    await db.run('DELETE FROM pfasgroups_results');
    await db.run('DELETE FROM pfasgroups_results_bycomponent');
    await db.run('DELETE FROM atlas_results');
    await db.run('DELETE FROM molecules');
    
    console.log('✅ Database cleared');
    console.log('📥 Re-importing data...\n');
    
    // Re-run the import script
    const importScript = path.join(__dirname, 'import-benchmark-data.js');
    require(importScript);
}

reimportAll().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
