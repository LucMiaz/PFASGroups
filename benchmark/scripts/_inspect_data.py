"""Inspect benchmark data file structure."""
import json
import os

DATA = r"C:\Users\luc\git\PFASGroups\benchmark\data"

# 1. Timing analysis
print("=" * 60)
print("TIMING ANALYSIS")
print("=" * 60)
with open(os.path.join(DATA, "timing_analysis_20260315_065228.json")) as f:
    d = json.load(f)

s = d["summary"]
for k, v in s.items():
    if not isinstance(v, (dict, list)):
        print(f"  {k}: {v}")
    elif isinstance(v, dict):
        print(f"  {k}:")
        for kk, vv in list(v.items())[:10]:
            if not isinstance(vv, (dict, list)):
                print(f"    {kk}: {vv}")
            else:
                print(f"    {kk}: {type(vv).__name__}")

# Check for per-molecule data
if "figures" in d:
    print(f"\n  figures: list of {len(d['figures'])} items")
    if d["figures"]:
        print(f"  first figure keys: {list(d['figures'][0].keys()) if isinstance(d['figures'][0], dict) else type(d['figures'][0])}")

print(f"\n  system_specs: {d.get('system_specs', 'N/A')}")

# 2. Clinventory comparison
print("\n" + "=" * 60)
print("CLINVENTORY COMPARISON")
print("=" * 60)
with open(os.path.join(DATA, "clinventory_comparison_20260315T065332.json")) as f:
    c = json.load(f)

for k, v in c.items():
    if k == "molecules":
        print(f"  {k}: list of {len(v)} molecules")
        if v:
            print(f"    first mol keys: {list(v[0].keys())}")
    elif isinstance(v, dict):
        print(f"  {k}:")
        for kk, vv in v.items():
            if not isinstance(vv, (dict, list)):
                print(f"    {kk}: {vv}")
    elif isinstance(v, list):
        print(f"  {k}: list[{len(v)}]")
    else:
        print(f"  {k}: {v}")

# 3. Telomer validation
print("\n" + "=" * 60)
print("TELOMER VALIDATION")
print("=" * 60)
with open(os.path.join(DATA, "telomer_validation_results.json")) as f:
    t = json.load(f)

for k, v in t.items():
    if k == "molecules" or k == "results":
        print(f"  {k}: list/dict of {len(v)} items")
        if isinstance(v, list) and v:
            print(f"    first item keys: {list(v[0].keys()) if isinstance(v[0], dict) else type(v[0])}")
        elif isinstance(v, dict):
            first_key = list(v.keys())[0]
            print(f"    first key: {first_key}")
            print(f"    first value: {type(v[first_key])}")
    elif isinstance(v, dict):
        print(f"  {k}: dict with {len(v)} entries")
        for kk, vv in list(v.items())[:15]:
            if not isinstance(vv, (dict, list)):
                print(f"    {kk}: {vv}")
    elif isinstance(v, list):
        print(f"  {k}: list[{len(v)}]")
    else:
        print(f"  {k}: {v}")

# 4. toxcast results CSV
print("\n" + "=" * 60)
print("TOXCAST RESULTS CSV")
print("=" * 60)
import csv
with open(os.path.join(DATA, "toxcast_comparison_results.csv")) as f:
    reader = csv.DictReader(f)
    cols = reader.fieldnames
    print(f"  columns: {cols}")
    rows = list(reader)
    print(f"  rows: {len(rows)}")
    # Show unique experiments, feature_sets, models, endpoints
    exps = set(r.get("experiment", "") for r in rows)
    fsets = set(r.get("feature_set", "") for r in rows)
    models = set(r.get("model", "") for r in rows)
    endpoints = set(r.get("endpoint", "") for r in rows)
    print(f"  experiments: {sorted(exps)}")
    print(f"  feature_sets: {sorted(fsets)}")
    print(f"  models: {sorted(models)}")
    print(f"  endpoints ({len(endpoints)}): {sorted(endpoints)}")
