import json

data = json.load(open('../data/pfas_highly_branched_benchmark_20260202_105121.json'))

print('First 5 records:')
for r in data['details'][:5]:
    print(f"  Group {r['group_id']}: pfasgroups_passed={r.get('pfasgroups_passed')}, atlas_passed={r.get('atlas_passed')}")

print(f"\nTotal records: {len(data['details'])}")
print(f"Metadata: {data['metadata']}")
