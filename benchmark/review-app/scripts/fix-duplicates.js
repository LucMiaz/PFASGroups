const Database = require('../database/database');

async function fixDuplicates() {
    const db = new Database();
    await db.waitForReady();
    
    console.log('Analyzing duplicates...');
    
    // Get count before cleanup
    const beforeCount = await db.get('SELECT COUNT(*) as count FROM molecules');
    console.log(`Total molecules before cleanup: ${beforeCount.count}`);
    
    const uniqueCount = await db.get('SELECT COUNT(DISTINCT smiles) as count FROM molecules');
    console.log(`Unique SMILES: ${uniqueCount.count}`);
    console.log(`Duplicates to remove: ${beforeCount.count - uniqueCount.count}`);
    
    // For each duplicate SMILES, keep only the first occurrence (lowest id)
    // EXCEPT for timing dataset - keep all triplicates for timing analysis
    console.log('\nRemoving duplicates...');
    
    // Find all SMILES with duplicates (excluding timing dataset)
    const duplicates = await db.all(`
        SELECT smiles, COUNT(*) as cnt 
        FROM molecules 
        WHERE dataset_type != 'timing'
        GROUP BY smiles 
        HAVING cnt > 1
    `);
    
    console.log(`Found ${duplicates.length} SMILES with duplicates (excluding timing dataset)`);
    
    let totalRemoved = 0;
    for (const dup of duplicates) {
        // Get all IDs for this SMILES (only non-timing)
        const ids = await db.all(
            'SELECT id FROM molecules WHERE smiles = ? AND dataset_type != ? ORDER BY id',
            [dup.smiles, 'timing']
        );
        
        // Keep the first ID, delete the rest
        const keepId = ids[0].id;
        const deleteIds = ids.slice(1).map(row => row.id);
        
        for (const deleteId of deleteIds) {
            // Delete from all related tables
            await db.run('DELETE FROM pfasgroups_results WHERE molecule_id = ?', [deleteId]);
            await db.run('DELETE FROM pfasgroups_results_bycomponent WHERE molecule_id = ?', [deleteId]);
            await db.run('DELETE FROM atlas_results WHERE molecule_id = ?', [deleteId]);
            await db.run('DELETE FROM manual_reviews WHERE molecule_id = ?', [deleteId]);
            await db.run('DELETE FROM molecules WHERE id = ?', [deleteId]);
            totalRemoved++;
        }
        
        if (totalRemoved % 100 === 0) {
            console.log(`Removed ${totalRemoved} duplicates...`);
        }
    }
    
    console.log(`\nTotal duplicates removed: ${totalRemoved}`);
    
    // Get count after cleanup
    const afterCount = await db.get('SELECT COUNT(*) as count FROM molecules');
    console.log(`Total molecules after cleanup: ${afterCount.count}`);
    
    // Verify no duplicates remain
    const remainingDups = await db.get(`
        SELECT COUNT(*) as count 
        FROM (
            SELECT smiles, COUNT(*) as cnt 
            FROM molecules 
            GROUP BY smiles 
            HAVING cnt > 1
        )
    `);
    
    if (remainingDups.count === 0) {
        console.log('✅ All duplicates successfully removed!');
    } else {
        console.log(`⚠️ Warning: ${remainingDups.count} SMILES still have duplicates`);
    }
    
    // Show stats by dataset
    console.log('\nMolecules by dataset:');
    const stats = await db.all(`
        SELECT dataset_type, COUNT(*) as count 
        FROM molecules 
        GROUP BY dataset_type 
        ORDER BY count DESC
    `);
    stats.forEach(stat => {
        console.log(`  ${stat.dataset_type}: ${stat.count}`);
    });
    
    process.exit(0);
}

fixDuplicates().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
