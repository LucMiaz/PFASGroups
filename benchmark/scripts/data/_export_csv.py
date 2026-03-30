import json
import csv

with open('PFASGroups/data/Halogen_groups_smarts.json', encoding='utf-8') as f:
    data = json.load(f)


def fmt_smarts(smarts: dict) -> str:
    if not smarts:
        return ""
    parts = []
    for pattern, count in smarts.items():
        parts.append(f"{pattern} ×{count}" if count != 1 else pattern)
    return " | ".join(parts)


def fmt_formula_constraints(c: dict) -> str:
    if not c:
        return ""
    parts = []
    if "eq" in c:
        inner = ", ".join(f"{k}={v}" for k, v in c["eq"].items())
        parts.append(f"eq({inner})")
    if "gte" in c:
        inner = ", ".join(f"{k}≥{v}" for k, v in c["gte"].items())
        parts.append(f"gte({inner})")
    if "lte" in c:
        inner = ", ".join(f"{k}≤{v}" for k, v in c["lte"].items())
        if inner:
            parts.append(f"lte({inner})")
    if "only" in c:
        parts.append(f"only({','.join(c['only'])})")
    if "rel" in c:
        rel_parts = []
        for atom, spec in c["rel"].items():
            atoms_str = "+".join(spec.get("atoms", []))
            add = spec.get("add", 0)
            div = spec.get("div", 1)
            sign = "+" if add >= 0 else ""
            rel_parts.append(f"{atom}≤({atoms_str}{sign}{add})/{div}")
        parts.append("rel(" + ", ".join(rel_parts) + ")")
    return "; ".join(parts)


def fmt_component_constraints(g: dict) -> str:
    parts = []
    comp = g.get("componentSmarts")
    if comp:
        parts.append(f"component={comp}")
    sat = g.get("componentSaturation")
    if sat:
        parts.append(f"saturation={sat}")
    form = g.get("componentForm")
    if form:
        parts.append(f"form={form}")
    dist = g.get("max_dist_from_comp")
    if dist is not None:
        parts.append(f"max_dist={dist}")
    linker = g.get("linker_smarts")
    if linker:
        parts.append(f"linker={linker}")
    return "; ".join(parts)


rows = []
for g in data:
    rows.append({
        "category":              g["sorting"]["category"],
        "id":                    g["id"],
        "Group Name":            g["name"],
        "SMARTS Patterns":       fmt_smarts(g.get("smarts", {})),
        "Formula Constraints":   fmt_formula_constraints(g.get("constraints", {})),
        "Component Constraints": fmt_component_constraints(g),
    })

out_path = "PFASGroups/data/Halogen_groups_smarts.csv"
fieldnames = ["category", "id", "Group Name", "SMARTS Patterns",
              "Formula Constraints", "Component Constraints"]

with open(out_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)

print(f"Written {len(rows)} rows to {out_path}")
