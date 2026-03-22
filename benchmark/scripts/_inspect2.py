import json

DATA = r"C:\Users\luc\git\PFASGroups\benchmark\data"

d = json.load(open(f"{DATA}/clinventory_comparison_20260315T065332.json"))

tb = d["timing_by_atom_bracket"]
print("TIMING BY BRACKET:")
for bracket, vals in tb.items():
    print(f"  {bracket}:")
    for k, v in vals.items():
        if not isinstance(v, (dict, list)):
            print(f"    {k}: {v}")

print()
print("TIMING OVERALL:")
to = d["timing_overall"]
for k, v in to.items():
    if isinstance(v, dict):
        print(f"  {k}:")
        for kk, vv in v.items():
            print(f"    {kk}: {vv}")
    else:
        print(f"  {k}: {v}")

print()
print("ATLAS CLASS1 DIST (top 10):")
ad = d.get("atlas_class1_distribution", {})
for k, v in list(ad.items())[:10]:
    print(f"  {k}: {v}")

print()
t = json.load(open(f"{DATA}/telomer_validation_results.json"))
print("TELOMER GROUP COUNTS:")
for item in t["group_counts"]:
    print(f"  {item}")
