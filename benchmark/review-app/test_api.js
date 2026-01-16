const Database = require('./database/database');
const fs = require('fs');
const path = require('path');

// Load PFAS groups mapping
let pfasGroupsMap = {};
let pfasGroupsData = [];
try {
    pfasGroupsData = JSON.parse(fs.readFileSync(path.join(__dirname, 'pfas_groups_map.json'), 'utf8'));
    pfasGroupsMap = Object.fromEntries(pfasGroupsData.map(g => [g.id, { name: g.name, alias: g.alias }]));
    console.log(`Loaded ${Object.keys(pfasGroupsMap).length} PFAS groups`);
} catch (error) {
    console.error('Warning: Could not load PFAS groups map:', error.message);
}

// Helper function to convert group IDs to names with path type information
function enrichGroupData(groupIds, matchedPathTypes = {}) {
    if (!Array.isArray(groupIds)) return [];
    return groupIds.map(id => {
        const groupInfo = pfasGroupsMap[id] || { name: `Group ${id}`, alias: `Group ${id}` };
        const pathType = matchedPathTypes[id] || null;
        
        // Map path types to abbreviations
        const pathTypeAbbrev = {
            'Perfluoroalkyl': 'per',
            'Polyfluoroalkyl': 'poly',
            'Perfluoro': 'per',
            'Polyfluoro': 'poly',
            'cyclic': 'cyc'
        };
        
        return {
            id: id,
            name: groupInfo.name,
            alias: groupInfo.alias,
            matchedPathType: pathType ? (pathTypeAbbrev[pathType] || pathType) : null,
            matchedPathTypeFull: pathType
        };
    });
}

(async () => {
    const db = new Database();
    await db.waitForReady();
    
    // Get a molecule from enhanced/oecd dataset (id >= 7)
    const mol = await db.get(`
        SELECT 
            m.*,
            pg.detected_groups as pfasgroups_detected,
            pg.matched_path_types as pfasgroups_matched_path_types
        FROM molecules m
        LEFT JOIN pfasgroups_results pg ON m.id = pg.molecule_id
        WHERE m.id = 7
    `);
    
    console.log('\n=== Raw database data ===');
    console.log('pfasgroups_detected:', mol.pfasgroups_detected);
    console.log('pfasgroups_matched_path_types:', mol.pfasgroups_matched_path_types);
    
    // Parse and enrich
    const matchedPathTypes = mol.pfasgroups_matched_path_types ? JSON.parse(mol.pfasgroups_matched_path_types) : {};
    const enrichedGroups = enrichGroupData(mol.pfasgroups_detected ? JSON.parse(mol.pfasgroups_detected) : [], matchedPathTypes);
    
    console.log('\n=== Enriched data ===');
    console.log(JSON.stringify(enrichedGroups, null, 2));
    
    await db.close();
    process.exit(0);
})();
