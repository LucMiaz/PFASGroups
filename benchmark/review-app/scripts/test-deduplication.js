const DataImporter = require('./import-benchmark-data');

// Create a test importer
const importer = new DataImporter();

console.log('Testing Deduplication Logic\n');
console.log('=' .repeat(60));

// Test 1: Benchmark data deduplication
console.log('\n📋 Test 1: Benchmark Data Deduplication');
console.log('-'.repeat(60));

const testBenchmarkData = [
    {
        molecule_data: { smiles: 'C(F)(F)F', molecular_weight: 100 },
        pfasgroups_result: { detected_groups: ['CF3'], success: true },
        atlas_result: { first_class: 'PFAS', second_class: 'Type1', success: true }
    },
    {
        molecule_data: { smiles: 'C(F)(F)F', molecular_weight: 100 },
        pfasgroups_result: { detected_groups: ['CF3'], success: true },
        atlas_result: { first_class: 'PFAS', second_class: 'Type1', success: true }
    },
    {
        molecule_data: { smiles: 'C(F)(F)F', molecular_weight: 100 },
        pfasgroups_result: { detected_groups: ['CF2'], success: true },
        atlas_result: { first_class: 'PFAS', second_class: 'Type2', success: true }
    },
    {
        molecule_data: { smiles: 'CC(F)(F)F', molecular_weight: 114 },
        pfasgroups_result: { detected_groups: ['CF3'], success: true },
        atlas_result: { first_class: 'PFAS', second_class: 'Type1', success: true }
    }
];

console.log(`Original count: ${testBenchmarkData.length}`);
const dedupBenchmark = importer.deduplicateRecords(testBenchmarkData);
console.log(`After deduplication: ${dedupBenchmark.length}`);
console.log(`Duplicates removed: ${testBenchmarkData.length - dedupBenchmark.length}`);

console.log('\nExpected: Remove 1 duplicate (first two records are identical)');
console.log('Keep: 3 records (first C(F)(F)F with CF3/Type1, second C(F)(F)F with CF2/Type2, CC(F)(F)F)');

// Test 2: Timing data deduplication
console.log('\n\n📊 Test 2: Timing Data Deduplication');
console.log('-'.repeat(60));

const testTimingData = [
    { smiles: 'C(F)(F)F', molecular_weight: 100, pfasgroups_time: 0.001 },
    { smiles: 'C(F)(F)F', molecular_weight: 100, pfasgroups_time: 0.002 },
    { smiles: 'CC(F)(F)F', molecular_weight: 114, pfasgroups_time: 0.001 },
    { smiles: 'C(F)(F)F', molecular_weight: 100, pfasgroups_time: 0.003 }
];

console.log(`Original count: ${testTimingData.length}`);
const dedupTiming = importer.deduplicateTimingRecords(testTimingData);
console.log(`After deduplication: ${dedupTiming.length}`);
console.log(`Duplicates removed: ${testTimingData.length - dedupTiming.length}`);

console.log('\nExpected: Remove 2 duplicates (SMILES-based only)');
console.log('Keep: 2 records (first C(F)(F)F, first CC(F)(F)F)');

// Test 3: Edge cases
console.log('\n\n🔧 Test 3: Edge Cases');
console.log('-'.repeat(60));

const edgeCases = [
    { molecule_data: { smiles: null }, pfasgroups_result: { detected_groups: [] } },
    { molecule_data: { smiles: 'CCF' }, pfasgroups_result: { detected_groups: [] } },
    { molecule_data: { smiles: 'CCF' }, pfasgroups_result: null },
    { molecule_data: { smiles: 'CCF' } }
];

console.log(`Original count: ${edgeCases.length}`);
const dedupEdge = importer.deduplicateRecords(edgeCases);
console.log(`After deduplication: ${dedupEdge.length}`);
console.log(`Duplicates removed: ${edgeCases.length - dedupEdge.length}`);

console.log('\nExpected: Handle null/undefined gracefully');
console.log('Keep: Records with null SMILES, and deduplicate based on available data');

console.log('\n' + '='.repeat(60));
console.log('✅ Deduplication tests completed!\n');

// Close database
importer.close().then(() => {
    console.log('Database connection closed.');
});
