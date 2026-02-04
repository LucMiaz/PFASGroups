const Database = require('./database/database');

async function checkDatabase() {
    const db = new Database();
    await db.initialize();
    
    const smiles = 'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CC[Si](Cl)(Cl)Cl';
    
    const result = await db.get(`
        SELECT 
            m.id,
            m.smiles,
            pg.detected_groups
        FROM molecules m
        LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
        WHERE m.smiles = ?
    `, [smiles]);
    
    if (result) {
        console.log('Found in database:');
        console.log('Molecule ID:', result.id);
        console.log('SMILES:', result.smiles);
        console.log('Detected groups (raw):', result.detected_groups);
        
        if (result.detected_groups) {
            const groups = JSON.parse(result.detected_groups);
            console.log('Detected groups (parsed):', groups);
            
            // Check mapping
            const fs = require('fs');
            const groupsMap = JSON.parse(fs.readFileSync('./pfas_groups_map.json', 'utf8'));
            const groupsDict = {};
            groupsMap.forEach(g => groupsDict[g.id] = g);
            
            console.log('\nMapping for these groups:');
            groups.forEach(gid => {
                if (groupsDict[gid]) {
                    console.log(`  ID ${gid}: ${groupsDict[gid].name} (alias: ${groupsDict[gid].alias})`);
                } else {
                    console.log(`  ID ${gid}: NOT FOUND IN MAP`);
                }
            });
        }
    } else {
        console.log('Molecule not found in database');
    }
    
    await db.close();
}

checkDatabase().catch(console.error);
