const Database = require('./database/database');

const db = new Database();

db.waitForReady().then(async () => {
    console.log('Finding duplicate SMILES with different results...\n');
    
    const dupes = await db.all(`
        SELECT m.smiles, COUNT(*) as count, GROUP_CONCAT(m.id) as ids, GROUP_CONCAT(pg.success) as successes
        FROM molecules m
        LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
        GROUP BY m.smiles
        HAVING COUNT(*) > 1 AND successes LIKE '%0%'
        ORDER BY m.id
        LIMIT 50
    `);
    
    console.log(`Found ${dupes.length} duplicate SMILES with errors\n`);
    
    const idsToDelete = [];
    
    dupes.forEach(d => {
        const ids = d.ids.split(',').map(Number);
        const successes = d.successes.split(',').map(Number);
        
        console.log(`SMILES: ${d.smiles.substring(0, 60)}...`);
        console.log(`  IDs: ${ids.join(', ')} | Success: ${successes.join(', ')}`);
        
        // Find IDs where success = 0
        ids.forEach((id, idx) => {
            if (successes[idx] === 0) {
                idsToDelete.push(id);
                console.log(`  -> Will delete ID ${id} (has error)`);
            }
        });
        console.log();
    });
    
    if (idsToDelete.length > 0) {
        console.log(`\nDeleting ${idsToDelete.length} molecules with errors...`);
        
        const placeholders = idsToDelete.map(() => '?').join(',');
        await db.run(`DELETE FROM pfasgroups_results WHERE molecule_id IN (${placeholders})`, idsToDelete);
        await db.run(`DELETE FROM atlas_results WHERE molecule_id IN (${placeholders})`, idsToDelete);
        await db.run(`DELETE FROM molecules WHERE id IN (${placeholders})`, idsToDelete);
        
        console.log('Deleted!');
        
        const count = await db.get('SELECT COUNT(*) as count FROM molecules');
        console.log(`Total molecules now: ${count.count}`);
    } else {
        console.log('No duplicates with errors to delete.');
    }
    
    setTimeout(() => process.exit(0), 200);
}).catch(e => {
    console.error('Error:', e);
    process.exit(1);
});
