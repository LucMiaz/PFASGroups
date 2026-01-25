const Database = require('./database/database');

const smiles_list = [
    'c1c(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)c(F)c(C(F)(F)C(=O)O)c(F)c1F',
    'C1(C(F)(F)F)(C(F)(F)F)C(C(F)(F)C(F)(F)F)C(C(F)(F)F)C(C(=O)O)C1(F)F'
];

const db = new Database();

db.waitForReady().then(async () => {
    console.log('Checking molecules in database:\n');
    
    for (const smiles of smiles_list) {
        console.log(`SMILES: ${smiles}`);
        
        const result = await db.all(`
            SELECT m.id, m.smiles, 
                   pg.detected_groups, pg.success as pfasgroups_success, pg.error_message as pfasgroups_error,
                   ar.first_class, ar.success as atlas_success, ar.error_message as atlas_error
            FROM molecules m
            LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
            LEFT JOIN atlas_results ar ON m.id = ar.molecule_id
            WHERE m.smiles = ?
        `, [smiles]);
        
        if (result.length === 0) {
            console.log('  NOT FOUND in database');
        } else {
            result.forEach(row => {
                console.log(`  ID: ${row.id}`);
                console.log(`  PFASGroups Success: ${row.pfasgroups_success}`);
                console.log(`  PFASGroups Error: ${row.pfasgroups_error}`);
                console.log(`  PFASGroups Groups: ${row.detected_groups}`);
                console.log(`  Atlas Success: ${row.atlas_success}`);
                console.log(`  Atlas Error: ${row.atlas_error}`);
            });
        }
        console.log('-'.repeat(80));
    }
    
    setTimeout(() => process.exit(0), 100);
}).catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
