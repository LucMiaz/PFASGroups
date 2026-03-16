#!/usr/bin/env python3
"""PFASGroups dependency + empirical co-occurrence network diagram.

Data sources
------------
pfasgroup_overlap.csv          - pairwise compound overlap from FluoroDB
specificity_test_groups.json   - structural dependency rules (per/poly/both/oneway)
Halogen_groups_smarts.json     - group metadata (category)

Visualisation (plotly, interactive HTML + PNG/PDF)
---------------------------------------------------
Nodes  sized by log^2(n_compounds), coloured by category (OECD=orange, generic=magenta,
       telomer=blue)
Rule edges  curved directed arrows coloured by type:
            per=orange, poly=magenta, both=dark-purple, oneway=grey
Residual edges  grey lines, linewidth proportional to jaccard
                (pairs not covered by any rule AND max(prop1,prop2) <= containment threshold)

Layout: networkx kamada_kawai_layout weighted by rules + high-jaccard co-occurrences.

Usage
-----
    python benchmark/scripts/group_dependencies.py \
        --overlap ../../zeropmdb/database/elements/data/pfasgroup_overlap.csv

    python benchmark/scripts/group_dependencies.py \
        --overlap /path/to/pfasgroup_overlap.csv \
        --min-jaccard 0.05 --top-residual 30 --outdir benchmark/imgs
"""
from __future__ import annotations
import argparse, json, math, sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    import plotly.graph_objects as go
    import plotly.io as pio
    PLOTLY = True
except ImportError:
    PLOTLY = False
    print("WARNING: plotly not available")

SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent
REPO_ROOT     = BENCHMARK_DIR.parent
GROUPS_JSON   = REPO_ROOT / "PFASGroups" / "data" / "Halogen_groups_smarts.json"
_RULES_CANDIDATES = [
    REPO_ROOT.parent / "zeropmdb" / "database" / "elements" / "data" / "specificity_test_groups.json",
    REPO_ROOT / "tests" / "test_data" / "specificity_test_groups.json",
]
RULES_JSON = next((p for p in _RULES_CANDIDATES if p.exists()),
                  _RULES_CANDIDATES[-1])

_C0 = "#E15D0B"; _C1 = "#306DBA"; _C2 = "#9D206C"; _C3 = "#51127C"
_GREY = "#888888"; _LGREY = "#BBBBBB"

_ETYPE_COLOUR: Dict[str, str] = {"per": _C0, "poly": _C2, "both": _C3, "oneway": _GREY}
_ETYPE_LABEL:  Dict[str, str] = {
    "per":    "per-fluorinated implication (source -> target)",
    "poly":   "poly-fluorinated implication (source -> target)",
    "both":   "per- + poly-fluorinated (source -> target)",
    "oneway": "saturation-agnostic (source -> target)",
}
_CAT_COLOUR: Dict[str, str] = {"OECD": _C0, "generic": _C2, "telomer": _C1, "other": _GREY}
_ALIAS_NORMALISE: Dict[str, str] = {
    "azole":                            "Heterocyclic azole",
    "azine":                            "Heterocyclic azine",
    "ethene":                           "alkene",
    "Polyhalogenated aryl compounds":   "Polyfluoroaryl compounds",
    "Perhalogenated aryl compounds":    "Perfluoroaryl compounds",
    "Polyhalogenated cyclic compounds": "Polyfluoro cyclic compounds",
    "Perhalogenated cyclic compounds":  "Perfluoro cyclic compounds",
}

def _rgba(h: str, a: float) -> str:
    h = h.lstrip("#")
    r, g, b = int(h[0:2],16), int(h[2:4],16), int(h[4:6],16)
    return f"rgba({r},{g},{b},{a:.2f})"

# ---------------------------------------------------------------------------
# Group metadata (no RDKit)
# ---------------------------------------------------------------------------
def _load_group_meta() -> Dict[int, Dict]:
    if not GROUPS_JSON.exists():
        return {}
    with GROUPS_JSON.open(encoding="utf-8") as f:
        raw = json.load(f)
    result: Dict[int, Dict] = {}
    for g in raw:
        name = g["name"]
        short = name[5:] if name.startswith("OECD ") else name
        result[g["id"]] = {"name": name, "short_name": short,
                           "category": g.get("test", {}).get("category", "other")}
    return result

_GROUP_META: Dict[int, Dict] = _load_group_meta()

def _build_name_index() -> Dict[str, int]:
    idx: Dict[str, int] = {}
    for gid, info in _GROUP_META.items():
        name = info["name"].lower().strip()
        idx[name] = gid
        if name.startswith("oecd "):
            idx[name[5:]] = gid
    return idx

_NAME_IDX: Dict[str, int] = _build_name_index()

def _resolve_name(raw_name: str) -> Optional[int]:
    name = _ALIAS_NORMALISE.get(raw_name, raw_name)
    key  = name.lower().strip()
    if key in _NAME_IDX:
        return _NAME_IDX[key]
    for k, gid in _NAME_IDX.items():
        if key in k and len(key) > 4:
            return gid
    for k, gid in _NAME_IDX.items():
        if k in key and len(k) > 4:
            return gid
    return None

# ---------------------------------------------------------------------------
# Load overlap CSV
# ---------------------------------------------------------------------------
def load_overlap(csv_path: Path):
    """Load pfasgroup_overlap.csv using stdlib csv only (no pandas needed for data loading).

    Returns group_data {gid: {alias, n_compounds, category}} and
    rows list of dicts with gid1, gid2, n_shared, jaccard, prop1, prop2.
    """
    import csv
    rows_raw = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows_raw.append(row)

    # Extract group sizes
    sizes: Dict[int, Dict] = {}
    for row in rows_raw:
        g1, g2 = int(row["group1_id"]), int(row["group2_id"])
        n1, n2 = int(row["n_group1"]),  int(row["n_group2"])
        alias1, alias2 = row["group1_alias"], row["group2_alias"]
        if g1 not in sizes or n1 > sizes[g1]["n"]:
            sizes[g1] = {"alias": alias1, "n": n1}
        if g2 not in sizes or n2 > sizes[g2]["n"]:
            sizes[g2] = {"alias": alias2, "n": n2}

    group_data: Dict[int, Dict] = {}
    for gid, s in sizes.items():
        meta = _GROUP_META.get(gid, {})
        group_data[gid] = {"alias": s["alias"], "n_compounds": s["n"],
                           "category": meta.get("category", "other")}

    pairs = []
    for row in rows_raw:
        pairs.append({
            "gid1":    int(row["group1_id"]),
            "gid2":    int(row["group2_id"]),
            "n_shared": int(row["n_shared"]),
            "jaccard": float(row["jaccard"]),
            "prop1":   float(row["prop_of_group1"]),
            "prop2":   float(row["prop_of_group2"]),
        })
    return group_data, pairs

# ---------------------------------------------------------------------------
# Load dependency rules
# ---------------------------------------------------------------------------
def load_rules(rules_path: Path, active_ids: set) -> List[Tuple[int, int, str]]:
    """Resolve rules to (src_gid, tgt_gid, edge_type). Deduplicates by keeping most specific type."""
    if not rules_path.exists():
        print(f"WARNING: rules file not found: {rules_path}")
        return []
    with rules_path.open(encoding="utf-8") as f:
        raw = json.load(f)
    _priority = {"per": 3, "poly": 3, "both": 2, "oneway": 1}
    best: Dict[Tuple[int, int], str] = {}
    skipped: set = set()
    for r in raw:
        src_id = _resolve_name(r["source"])
        tgt_id = _resolve_name(r["target"])
        if src_id is None: skipped.add(r["source"]); continue
        if tgt_id is None: skipped.add(r["target"]); continue
        if src_id not in active_ids or tgt_id not in active_ids: continue
        key = (src_id, tgt_id); et = r["edge_type"]
        if key not in best or _priority.get(et, 0) > _priority.get(best[key], 0):
            best[key] = et
    if skipped:
        print("  Rules: " + str(len(skipped)) + " unresolvable: "
              + ", ".join(f"'{n}'" for n in sorted(skipped)))
    edges = [(s, t, et) for (s, t), et in best.items()]
    print(f"  Rules: {len(raw)} raw -> {len(edges)} resolved edges, {len(skipped)} unresolvable")
    return edges

# ---------------------------------------------------------------------------
# Layout — networkx kamada_kawai (needs scipy) with pure-Python FR fallback
# ---------------------------------------------------------------------------
def compute_layout(
    group_ids: List[int],
    rule_edges: List[Tuple[int, int, str]],
    pairs: List[Dict],
    min_jaccard_layout: float = 0.05,
) -> Dict[int, Tuple[float, float]]:
    """Return {gid: (x,y)} using networkx kamada_kawai_layout (preferred) or
    a pure-Python Fruchterman-Reingold fallback if scipy is unavailable.

    Rule edges use weight 2.0; overlap pairs contribute their jaccard value.
    """
    try:
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(group_ids)
        for s, t, _ in rule_edges:
            if s in G.nodes and t in G.nodes:
                existing = G[s][t]["weight"] if G.has_edge(s, t) else 0.0
                G.add_edge(s, t, weight=max(existing, 2.0))
        for row in pairs:
            if row["jaccard"] < min_jaccard_layout:
                continue
            g1, g2 = row["gid1"], row["gid2"]
            if g1 in G.nodes and g2 in G.nodes:
                w = row["jaccard"]
                if G.has_edge(g1, g2):
                    G[g1][g2]["weight"] = max(G[g1][g2]["weight"], w)
                else:
                    G.add_edge(g1, g2, weight=w)
        pos_nx = nx.kamada_kawai_layout(G, weight="weight")
        print("  Layout: kamada_kawai")
        return {int(g): (float(x), float(y)) for g, (x, y) in pos_nx.items()}
    except Exception as exc:
        print(f"  kamada_kawai failed ({exc}); falling back to pure-Python FR layout")
        return _fr_layout(group_ids, rule_edges, pairs, min_jaccard_layout)


def _fr_layout(
    group_ids: List[int],
    rule_edges: List[Tuple[int, int, str]],
    pairs: List[Dict],
    min_jaccard_layout: float = 0.05,
    iterations: int = 500,
    seed: int = 42,
) -> Dict[int, Tuple[float, float]]:
    """Pure-Python Fruchterman-Reingold, no numpy/scipy required."""
    import random
    rng = random.Random(seed)
    n = max(len(group_ids), 1)
    pos: Dict[int, List[float]] = {}
    for i, gid in enumerate(group_ids):
        angle = 2.0 * math.pi * i / n
        r = 1.0 + rng.uniform(-0.04, 0.04)
        pos[gid] = [r * math.cos(angle), r * math.sin(angle)]
    adj: Dict[Tuple[int, int], float] = {}
    for s, t, _ in rule_edges:
        if s in pos and t in pos:
            key = (min(s, t), max(s, t))
            adj[key] = max(adj.get(key, 0.0), 2.0)
    for row in pairs:
        g1, g2 = row["gid1"], row["gid2"]
        if g1 in pos and g2 in pos and row["jaccard"] >= min_jaccard_layout:
            key = (min(g1, g2), max(g1, g2))
            adj[key] = max(adj.get(key, 0.0), row["jaccard"])
    gids_list = list(pos.keys())
    k_fr = math.sqrt(1.0 / n)
    for it in range(iterations):
        temp = 0.15 * (1.0 - it / iterations)
        disp: Dict[int, List[float]] = {gid: [0.0, 0.0] for gid in gids_list}
        for i in range(len(gids_list)):
            u = gids_list[i]
            for j in range(i + 1, len(gids_list)):
                v = gids_list[j]
                dx = pos[u][0] - pos[v][0]
                dy = pos[u][1] - pos[v][1]
                d  = math.sqrt(dx * dx + dy * dy) or 1e-4
                f  = k_fr * k_fr / d
                disp[u][0] += dx / d * f
                disp[u][1] += dy / d * f
                disp[v][0] -= dx / d * f
                disp[v][1] -= dy / d * f
        for (a, b), w in adj.items():
            dx = pos[a][0] - pos[b][0]
            dy = pos[a][1] - pos[b][1]
            d  = math.sqrt(dx * dx + dy * dy) or 1e-4
            f  = d * d / k_fr * w
            disp[a][0] -= dx / d * f
            disp[a][1] -= dy / d * f
            disp[b][0] += dx / d * f
            disp[b][1] += dy / d * f
        for gid in gids_list:
            disp[gid][0] -= 0.02 * pos[gid][0]
            disp[gid][1] -= 0.02 * pos[gid][1]
        for gid in gids_list:
            dx, dy = disp[gid]
            d = math.sqrt(dx * dx + dy * dy) or 1e-4
            scale = min(d, temp) / d
            pos[gid][0] += dx * scale
            pos[gid][1] += dy * scale
    return {gid: (pos[gid][0], pos[gid][1]) for gid in gids_list}

# ---------------------------------------------------------------------------
# Plot (plotly)
# ---------------------------------------------------------------------------
def plot_network(
    group_data: Dict[int, Dict],
    rule_edges: List[Tuple[int, int, str]],
    pairs: List[Dict],
    pos: Dict[int, Tuple[float, float]],
    *,
    min_jaccard: float = 0.05,
    containment_threshold: float = 0.90,
    top_residual: int = 25,
    title: str = "PFASGroups: structural dependencies + co-occurrence",
    output_base: Optional[Path] = None,
) -> None:
    if not PLOTLY:
        return

    rule_pair_set: Set[frozenset] = {frozenset([s, t]) for s, t, _ in rule_edges}
    active = set(pos)

    # Filter residual co-occurrence edges
    residual = []
    for row in pairs:
        g1, g2 = row["gid1"], row["gid2"]
        if g1 not in active or g2 not in active:
            continue
        if max(row["prop1"], row["prop2"]) > containment_threshold:
            continue  # containment
        if frozenset([g1, g2]) in rule_pair_set:
            continue  # already a rule
        if row["jaccard"] < min_jaccard:
            continue
        residual.append(row)
    residual.sort(key=lambda r: -r["jaccard"])
    residual = residual[:top_residual]

    fig = go.Figure()
    max_n = max(d["n_compounds"] for d in group_data.values())

    # 1. Residual co-occurrence edges (one trace per edge, grey, width prop. to jaccard)
    max_jac = residual[0]["jaccard"] if residual else 1.0
    for row in residual:
        g1, g2 = row["gid1"], row["gid2"]
        x0, y0 = pos[g1]; x1, y1 = pos[g2]
        jac = row["jaccard"]
        lw  = 0.8 + 5.0 * jac / max_jac
        alp = 0.20 + 0.55 * jac / max_jac
        alias1 = group_data.get(g1, {}).get("alias", str(g1))
        alias2 = group_data.get(g2, {}).get("alias", str(g2))
        n_shared = row["n_shared"]
        fig.add_trace(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            mode="lines",
            line=dict(color=_rgba(_LGREY, alp), width=lw),
            hovertext=(f"<b>[{g1}] {alias1}</b>  +  <b>[{g2}] {alias2}</b><br>"
                       f"Shared: {n_shared:,}<br>"
                       f"Jaccard: {jac:.3f}<br>"
                       f"prop([{g1}] in [{g2}]): {row['prop1']:.1%}<br>"
                       f"prop([{g2}] in [{g1}]): {row['prop2']:.1%}<br>"
                       "(non-rule co-occurrence)"),
            hoverinfo="text",
            showlegend=False,
        ))

    # 2. Rule edges (one trace per type for legend, arrows via annotations)
    annotations = []
    _rad_sign = {"per": +1, "poly": -1, "both": +1, "oneway": -1}
    legend_done: Set[str] = set()
    for s, t, et in rule_edges:
        if s not in pos or t not in pos:
            continue
        x0, y0 = pos[s]; x1, y1 = pos[t]
        col = _ETYPE_COLOUR.get(et, _GREY)
        show_leg = et not in legend_done
        if show_leg:
            legend_done.add(et)
        fig.add_trace(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            mode="lines",
            line=dict(color=_rgba(col, 0.55), width=1.5),
            name=_ETYPE_LABEL[et] if show_leg else "",
            legendgroup=f"rule_{et}",
            showlegend=show_leg,
            hovertext=(f"[{s}] -> [{t}]  [{et}]<br>"
                       f"{group_data.get(s,{}).get('alias','')} -> "
                       f"{group_data.get(t,{}).get('alias','')}"),
            hoverinfo="text",
        ))
        annotations.append(dict(
            x=x1, y=y1, ax=x0, ay=y0,
            xref="x", yref="y", axref="x", ayref="y",
            showarrow=True, arrowhead=2, arrowsize=0.9, arrowwidth=1.5,
            arrowcolor=_rgba(col, 0.75), text="",
        ))

    # 3. Residual legend entry
    if residual:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color=_rgba(_LGREY, 0.70), width=3),
            name=f"Non-rule co-occurrence (Jaccard>={min_jaccard}, top {top_residual};<br>"
                 f"width prop. to Jaccard; containment>{containment_threshold:.0%} excluded)",
            showlegend=True,
        ))

    # 4. Category legend
    for cat, col in _CAT_COLOUR.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=11, color=col, line=dict(color="white", width=1.5)),
            name=f"{cat} group", showlegend=True,
        ))

    # 5. Nodes
    gids     = sorted(pos)
    node_xs  = [pos[g][0] for g in gids]
    node_ys  = [pos[g][1] for g in gids]
    node_c   = [_CAT_COLOUR.get(group_data.get(g,{}).get("category","other"), _GREY)
                for g in gids]
    node_sz  = [10 + 30*(math.log1p(group_data.get(g,{}).get("n_compounds",1))/math.log1p(max_n))**2
                for g in gids]
    labels   = [f"[{g}] {group_data.get(g,{}).get('alias','')}"
                + f"<br>({group_data.get(g,{}).get('n_compounds',0):,})" for g in gids]
    hover    = []
    for g in gids:
        d = group_data.get(g, {})
        hover.append(f"<b>[{g}] {d.get('alias','')}</b><br>"
                     f"Category: {d.get('category','?')}<br>"
                     f"Compounds: {d.get('n_compounds',0):,}")
    fig.add_trace(go.Scatter(
        x=node_xs, y=node_ys, mode="markers+text",
        marker=dict(size=node_sz, color=node_c,
                    line=dict(color="white", width=1.8)),
        text=labels, textposition="middle right", textfont=dict(size=9),
        hovertext=hover, hoverinfo="text", showlegend=False, cliponaxis=False,
    ))

    n_nod  = len(pos)
    n_rule = len(rule_edges)
    n_res  = len(residual)
    fig.update_layout(
        annotations=annotations,
        title=dict(
            text=(f"<b>{title}</b><br>"
                  f"<sup>{n_nod} groups  .  {n_rule} rule edges  .  {n_res} residual edges  "
                  f"(containment>{containment_threshold:.0%} excluded; "
                  f"node size = log^2(n_compounds))</sup>"),
            x=0.5, font=dict(size=13),
        ),
        font=dict(family="Ubuntu, DejaVu Sans, sans-serif", size=10),
        paper_bgcolor="white", plot_bgcolor="#fafafa",
        showlegend=True,
        legend=dict(x=1.01, y=1.0, xanchor="left", bgcolor="white",
                    bordercolor=_LGREY, borderwidth=1, font=dict(size=9)),
        height=900, width=1400,
        xaxis=dict(visible=False), yaxis=dict(visible=False, scaleanchor="x", scaleratio=1),
        margin=dict(l=20, r=340, t=80, b=20),
        hovermode="closest",
    )

    if output_base is not None:
        html_path = output_base.with_suffix(".html")
        fig.write_html(str(html_path), include_plotlyjs="cdn")
        print(f"  Saved: {html_path}")
        for ext, kw in [("png", {"scale": 2}), ("pdf", {})]:
            try:
                p = output_base.with_suffix(f".{ext}")
                fig.write_image(str(p), **kw)
                print(f"  Saved: {p}")
            except Exception as exc:
                print(f"  {ext.upper()} skipped: {str(exc).split(chr(10))[0][:80]}")
    else:
        fig.show()

# ---------------------------------------------------------------------------
# Default CSV lookup
# ---------------------------------------------------------------------------
def _find_default_overlap() -> Optional[Path]:
    candidates = [
        REPO_ROOT.parent / "zeropmdb" / "database" / "elements" / "data" / "pfasgroup_overlap.csv",
        BENCHMARK_DIR / "data" / "pfasgroup_overlap.csv",
        Path("pfasgroup_overlap.csv"),
    ]
    for c in candidates:
        if c.exists():
            return c
    return None

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="PFASGroups dependency + co-occurrence network diagram",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--overlap", type=Path, default=_find_default_overlap(),
                   metavar="CSV", help="Path to pfasgroup_overlap.csv")
    p.add_argument("--min-compounds",         type=int,   default=1000,
                   help="Exclude groups with fewer than this many compounds")
    p.add_argument("--min-jaccard",           type=float, default=0.05)
    p.add_argument("--containment-threshold", type=float, default=0.90)
    p.add_argument("--top-residual",          type=int,   default=25)
    p.add_argument("--outdir", type=Path, default=BENCHMARK_DIR / "imgs")
    p.add_argument("--title",
                   default="PFASGroups: structural dependencies + empirical co-occurrence")
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    if args.overlap is None or not Path(args.overlap).exists():
        sys.exit(
            "ERROR: pfasgroup_overlap.csv not found.\n"
            "Supply via --overlap PATH\n"
            "e.g. --overlap .../zeropmdb/database/elements/data/pfasgroup_overlap.csv"
        )
    print(f"Reading overlap data from {args.overlap} ...")
    group_data, pairs = load_overlap(Path(args.overlap))
    print(f"  {len(group_data)} groups, {len(pairs):,} overlap pairs")

    # Filter small groups if requested
    if args.min_compounds > 0:
        before = len(group_data)
        group_data = {g: d for g, d in group_data.items()
                      if d["n_compounds"] >= args.min_compounds}
        pairs = [r for r in pairs
                 if r["gid1"] in group_data and r["gid2"] in group_data]
        print(f"  After --min-compounds={args.min_compounds}: "
              f"{len(group_data)} groups (removed {before - len(group_data)}), "
              f"{len(pairs):,} pairs")

    print(f"Loading dependency rules from {RULES_JSON} ...")
    rule_edges = load_rules(RULES_JSON, set(group_data))
    print("Computing layout (pure-Python FR, may take a few seconds) ...")
    pos = compute_layout(list(group_data), rule_edges, pairs)
    args.outdir.mkdir(exist_ok=True, parents=True)
    ts      = datetime.now().strftime("%Y%m%d_%H%M%S")
    outpath = args.outdir / f"group_dependencies_{ts}"
    print(f"Drawing -> {outpath}.html ...")
    plot_network(group_data=group_data, rule_edges=rule_edges, pairs=pairs, pos=pos,
                 min_jaccard=args.min_jaccard, containment_threshold=args.containment_threshold,
                 top_residual=args.top_residual, title=args.title, output_base=outpath)
    print("Done.")

if __name__ == "__main__":
    main()
