#!/usr/bin/env python3
"""Hierarchical Sankey of PFASGroups compound counts.

Visualises the full FluoroDB population flowing through the PFASGroups
classification hierarchy.  Layout:

  col 0  "FluoroDB"  (total)
  col 1  top-level generic groups  (polyfluoroalkyl / Side-chain aromatics /
          Perfluoroaryl / Polyfluoroaryl / Perfluoro cyclic / Polyfluoro cyclic …)
  col 2  functional-group children  (ether, amide, alcohol, amine…)
  col 3  more-specific children     (sulfonamide, sulfonic acid …)

Flow widths are proportional to n_compounds.
Containment edges (prop_of_child ≈ 1.0) are omitted as Sankey links; instead
the parent→child relationship determines the node column.

Usage (from benchmark/ directory):
    python scripts/group_hierarchy_sankey.py
    python scripts/group_hierarchy_sankey.py \\
        --overlap ../../zeropmdb/database/elements/data/pfasgroup_overlap.csv
    python scripts/group_hierarchy_sankey.py --min-compounds 5000
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    import plotly.graph_objects as go
    import plotly.io as pio
except ImportError:
    sys.exit("ERROR: plotly not available.  pip install plotly")
else:
    pio.templates.default = "plotly_white"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]
REPO_ROOT     = BENCHMARK_DIR.parent

_RULES_CANDIDATES = [
    REPO_ROOT.parent / "zeropmdb" / "database" / "elements" / "data" / "specificity_test_groups.json",
    REPO_ROOT / "tests" / "test_data" / "specificity_test_groups.json",
]
_RULES_DEFAULT = next((p for p in _RULES_CANDIDATES if p.exists()), _RULES_CANDIDATES[-1])

_OVERLAP_CANDIDATES = [
    REPO_ROOT.parent / "zeropmdb" / "database" / "elements" / "data" / "pfasgroup_overlap.csv",
    BENCHMARK_DIR / "data" / "pfasgroup_overlap.csv",
]
_OVERLAP_DEFAULT = next((p for p in _OVERLAP_CANDIDATES if p.exists()), None)

_GROUPS_JSON = REPO_ROOT / "PFASGroups" / "data" / "Halogen_groups_smarts.json"

# ---------------------------------------------------------------------------
# Colour helpers  (no numpy)
# ---------------------------------------------------------------------------
def _load_palette() -> List[str]:
    defaults = ["#E15D0B", "#306DBA", "#9D206C", "#51127C"]
    try:
        _p = REPO_ROOT / "PFASGroups" / "data" / "color_scheme.yaml"
        found = re.findall(r'"(#[0-9A-Fa-f]{6})"', _p.read_text())
        if len(found) >= 4:
            return found[:4]
    except Exception:
        pass
    return defaults

_C0, _C1, _C2, _C3 = _load_palette()
_GREY  = "#AAAAAA"
_LGREY = "#DDDDDD"


def _hex_to_rgb(h: str) -> Tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _lerp(c1: str, c2: str, t: float) -> str:
    r1, g1, b1 = _hex_to_rgb(c1)
    r2, g2, b2 = _hex_to_rgb(c2)
    return "#{:02X}{:02X}{:02X}".format(
        int(r1 + (r2 - r1) * t),
        int(g1 + (g2 - g1) * t),
        int(b1 + (b2 - b1) * t),
    )


def _rgba(h: str, a: float) -> str:
    r, g, b = _hex_to_rgb(h)
    return f"rgba({r},{g},{b},{a:.2f})"


# ---------------------------------------------------------------------------
# Alias normalisation  (rule names → overlap CSV aliases)
# ---------------------------------------------------------------------------
_ALIAS_MAP: Dict[str, str] = {
    "azole":                            "Heterocyclic azole",
    "azine":                            "Heterocyclic azine",
    "ethene":                           "alkene",
    "Polyhalogenated aryl compounds":   "Polyfluoroaryl compounds",
    "Perhalogenated aryl compounds":    "Perfluoroaryl compounds",
    "Polyhalogenated cyclic compounds": "Polyfluoro cyclic compounds",
    "Perhalogenated cyclic compounds":  "Perfluoro cyclic compounds",
    "fluorotelomer alcohols":           "alcohol",
    "iodide":                           "iodide",
}


def _normalise(name: str) -> str:
    return _ALIAS_MAP.get(name, name)


# ---------------------------------------------------------------------------
# Category colours
# ---------------------------------------------------------------------------
def _load_group_categories() -> Dict[int, str]:
    """Return {gid: category} from Halogen_groups_smarts.json."""
    if not _GROUPS_JSON.exists():
        return {}
    try:
        raw = json.loads(_GROUPS_JSON.read_text(encoding="utf-8"))
        return {g["id"]: g.get("test", {}).get("category", "other") for g in raw}
    except Exception:
        return {}


_CATEGORY_COLOUR: Dict[str, str] = {
    "OECD":    _C0,
    "generic": _C2,
    "telomer": _C1,
    "other":   _GREY,
}

_GID_CATEGORY: Dict[int, str] = _load_group_categories()


def _node_colour(gid: int) -> str:
    cat = _GID_CATEGORY.get(gid, "other")
    return _CATEGORY_COLOUR.get(cat, _GREY)


# ---------------------------------------------------------------------------
# Load overlap CSV
# ---------------------------------------------------------------------------
def load_overlap(
    csv_path: Path,
) -> Tuple[Dict[int, Dict], Dict[Tuple[int, int], Dict]]:
    """Return (group_data, overlap_lookup).

    group_data:      {gid: {alias, n_compounds}}
    overlap_lookup:  {(gid1, gid2): {n_shared, jaccard, prop1, prop2}}
                     keyed both ways.
    """
    raw_rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    group_data: Dict[int, Dict] = {}
    for r in raw_rows:
        g1, g2 = int(r["group1_id"]), int(r["group2_id"])
        n1, n2 = int(r["n_group1"]),  int(r["n_group2"])
        if g1 not in group_data:
            group_data[g1] = {"alias": r["group1_alias"], "n_compounds": n1}
        if g2 not in group_data:
            group_data[g2] = {"alias": r["group2_alias"], "n_compounds": n2}

    lookup: Dict[Tuple[int, int], Dict] = {}
    for r in raw_rows:
        g1, g2 = int(r["group1_id"]), int(r["group2_id"])
        d = {
            "n_shared": int(r["n_shared"]),
            "jaccard":  float(r["jaccard"]),
            "prop1":    float(r["prop_of_group1"]),
            "prop2":    float(r["prop_of_group2"]),
        }
        lookup[(g1, g2)] = d
        lookup[(g2, g1)] = {"n_shared": d["n_shared"], "jaccard": d["jaccard"],
                            "prop1": d["prop2"], "prop2": d["prop1"]}
    return group_data, lookup


# ---------------------------------------------------------------------------
# Load + resolve rules
# ---------------------------------------------------------------------------
def load_rules(
    rules_path: Path,
    group_data: Dict[int, Dict],
) -> List[Dict]:
    """Return list of {src_id, tgt_id, src_alias, tgt_alias, edge_type}."""
    if not rules_path.exists():
        print(f"WARNING: rules file not found: {rules_path}")
        return []

    alias_to_id: Dict[str, int] = {
        d["alias"].lower(): gid for gid, d in group_data.items()
    }

    raw = json.loads(rules_path.read_text(encoding="utf-8"))
    resolved = []
    for r in raw:
        src_norm = _normalise(r["source"])
        tgt_norm = _normalise(r["target"])
        sid = alias_to_id.get(src_norm.lower())
        tid = alias_to_id.get(tgt_norm.lower())
        if sid is None or tid is None:
            continue
        if sid == tid:
            continue
        resolved.append({
            "src_id":    sid,
            "tgt_id":    tid,
            "src_alias": group_data[sid]["alias"],
            "tgt_alias": group_data[tid]["alias"],
            "edge_type": r["edge_type"],
        })
    print(f"  Rules: {len(raw)} raw -> {len(resolved)} resolved")
    return resolved


# ---------------------------------------------------------------------------
# Build hierarchy  (assign column / depth to each node)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Backbone group IDs – pinned to fixed columns in the hierarchy.
# ---------------------------------------------------------------------------
# col 0 (display col 1) – broad poly classes, one per fluorocarbon type:
_POLY_ID:        int = 35  # polyfluoroalkyl               (3.8 M)
_POLY_ARYL_ID:   int = 38  # Polyfluoroaryl compounds      (1.01 M)
_POLY_CYCLIC_ID: int = 45  # Polyfluoro cyclic compounds   (388 K)

# col 1 (display col 2) – fully-fluorinated sub-classes:
_PERF_ID:        int = 34  # perfluoroalkyl                (3.27 M)
_PERF_ARYL_ID:   int = 37  # Perfluoroaryl compounds       (1.01 M, ≡ poly-aryl in data)
_PERF_CYCLIC_ID: int = 44  # Perfluoro cyclic compounds    (387 K,  ≈ poly-cyclic in data)

_BACKBONE_COL0: Set[int] = {_POLY_ID, _POLY_ARYL_ID, _POLY_CYCLIC_ID}
_BACKBONE_COL1: Set[int] = {_PERF_ID, _PERF_ARYL_ID, _PERF_CYCLIC_ID}
_BACKBONE:      Set[int] = _BACKBONE_COL0 | _BACKBONE_COL1

# Pairs (poly → perf) within each fluorocarbon class:
_POLY_TO_PERF: Dict[int, int] = {
    _POLY_ID:        _PERF_ID,
    _POLY_ARYL_ID:   _PERF_ARYL_ID,
    _POLY_CYCLIC_ID: _PERF_CYCLIC_ID,
}


def build_columns(
    group_data: Dict[int, Dict],
    rules: List[Dict],
    min_compounds: int,
) -> Tuple[Dict[int, int], Dict[int, Set[int]]]:
    """Return (col_of_node, children_of_node).

    Column assignment:
      col 0  polyfluoroalkyl, Polyfluoroaryl, Polyfluoro cyclic  (pinned)
      col 1  perfluoroalkyl, Perfluoroaryl, Perfluoro cyclic     (pinned)
      col 2  functional-group rule roots
      col 3  rule-children
      col 4  rule-grandchildren
    Backbone groups are placed by hard-coding because the rule file does not
    contain the backbone edges.  All other groups start at col 2+.
    """
    active: Set[int] = {
        gid for gid, d in group_data.items()
        if d["n_compounds"] >= min_compounds
    }

    # Build parent/child sets from rules, excluding ALL backbone as rule-parents.
    parents:  Dict[int, Set[int]] = {g: set() for g in active}
    children: Dict[int, Set[int]] = {g: set() for g in active}
    for r in rules:
        s, t = r["src_id"], r["tgt_id"]
        if s in active and t in active and s not in _BACKBONE:
            parents[t].add(s)
            children[s].add(t)

    # Assign backbone columns
    col: Dict[int, int] = {}
    for g in _BACKBONE_COL0:
        if g in active:
            col[g] = 0
    for g in _BACKBONE_COL1:
        if g in active:
            col[g] = 1

    # Rule-hierarchy roots (no rule-parent, not backbone) → col 2
    rule_roots = [g for g in active if g not in _BACKBONE and not parents[g]]
    for g in rule_roots:
        col[g] = 2

    # BFS for deeper rule levels (col 3, 4, …)
    visited: Set[int] = set(rule_roots) | _BACKBONE
    queue   = list(rule_roots)
    while queue:
        next_q = []
        for g in queue:
            for c in children[g]:
                new_col = col[g] + 1
                if c not in col or col[c] < new_col:
                    col[c] = new_col
                if c not in visited:
                    visited.add(c)
                    next_q.append(c)
        queue = next_q

    # Safety net: any active group not yet assigned → col 2
    for g in active:
        if g not in col:
            col[g] = 2

    return col, children


# ---------------------------------------------------------------------------
# Build Sankey flows
# ---------------------------------------------------------------------------

# Sentinel IDs for virtual nodes (all negative, cannot clash with real group IDs).
_FLUORODB_ID  = -1    # "FluoroDB" root                  display col 0
_POLY_ONLY_ID = -2    # "Polyfluoroalkyl only"            display col 2
# Category-aggregate sentinels for compact mode:
_OTHER_OECD_ID    = -10
_OTHER_TELOMER_ID = -11
_OTHER_GENERIC_ID = -12
_OTHER_REST_ID    = -13
_CAT_SENTINEL: Dict[str, int] = {
    "OECD":    _OTHER_OECD_ID,
    "telomer": _OTHER_TELOMER_ID,
    "generic": _OTHER_GENERIC_ID,
    "other":   _OTHER_REST_ID,
}


def build_sankey(
    group_data:         Dict[int, Dict],
    overlap:            Dict[Tuple[int, int], Dict],
    rules:              List[Dict],
    col_of:             Dict[int, int],
    children_of:        Dict[int, Set[int]],
    min_flow:           int = 100,
    min_root_compounds: int = 0,
    dataset_name:       str = "FluoroDB",
) -> go.Figure:
    """Build and return a plotly Sankey figure.

    Display layout (left → right):
      col 0  FluoroDB (virtual)
      col 1  [polyfluoroalkyl, Polyfluoroaryl, Polyfluoro cyclic]
              — backbone poly nodes, slightly staggered in x
      col 2  [perfluoroalkyl, Perfluoroaryl, Perfluoro cyclic, Polyfluoroalkyl-only]
              — backbone perf nodes + the partially-fluorinated alkyl remainder
      col 3  functional-group roots (col_of==2)
             OR "Other [category] groups" aggregates in compact mode
      col 4  rule-children
      col 5  rule-grandchildren

    arrangement='snap' lets Plotly auto-determine column layout from the DAG topology.
    Overlaps between backbone groups (poly-aryl ↔ poly-alkyl, etc.) are noted in
    hover-tooltips rather than as explicit cross-links, keeping flows conservative.
    """
    from collections import defaultdict

    active = set(col_of.keys())
    active_backbone = _BACKBONE & active

    # ─── Backbone flow values ────────────────────────────────────────────────
    def _n(gid: int) -> int:
        return group_data.get(gid, {}).get("n_compounds", 0)

    def _ovl(a: int, b: int) -> int:
        return overlap.get((a, b), {}).get("n_shared", 0)

    n_poly         = _n(_POLY_ID)
    n_perf         = _n(_PERF_ID)
    n_poly_aryl    = _n(_POLY_ARYL_ID)
    n_poly_cyclic  = _n(_POLY_CYCLIC_ID)
    n_perf_in_poly = _ovl(_POLY_ID, _PERF_ID)
    n_poly_only    = max(0, n_poly - n_perf_in_poly)

    # FluoroDB rough total (sum of the 3 poly backbones – direct pairwise overlaps;
    # ignores triple overlap, so this is an overestimate, noted in label):
    ov_poly_aryl   = _ovl(_POLY_ID, _POLY_ARYL_ID)
    ov_poly_cyclic = _ovl(_POLY_ID, _POLY_CYCLIC_ID)
    ov_aryl_cyclic = _ovl(_POLY_ARYL_ID, _POLY_CYCLIC_ID)
    n_fluorodb_est = (n_poly + n_poly_aryl + n_poly_cyclic
                      - ov_poly_aryl - ov_poly_cyclic - ov_aryl_cyclic)

    # ─── Functional-root filtering ───────────────────────────────────────────
    col2_all = [g for g in active if col_of[g] == 2]
    col2_shown = sorted(
        [g for g in col2_all
         if children_of.get(g) or
            group_data[g]["n_compounds"] >= min_root_compounds],
        key=lambda g: -group_data[g]["n_compounds"],
    )
    col2_hidden = [g for g in col2_all if g not in set(col2_shown)]

    # BFS: propagate shown status to rule-children/grandchildren
    shown_deep: Set[int] = set(col2_shown)
    bfs_q = list(col2_shown)
    while bfs_q:
        next_q = []
        for g in bfs_q:
            for c in children_of.get(g, set()):
                if c in active and c not in shown_deep:
                    shown_deep.add(c)
                    next_q.append(c)
        bfs_q = next_q

    all_real_shown: Set[int] = active_backbone | shown_deep

    # ─── Compact-mode category aggregates ───────────────────────────────────
    # Group the hidden col-2 roots by category and create virtual aggregate nodes.
    # Each aggregate node appears at display col 3, flow from perfluoroalkyl.
    cat_hidden: Dict[str, List[int]] = defaultdict(list)
    for g in col2_hidden:
        cat = _GID_CATEGORY.get(g, "other")
        cat_hidden[cat].append(g)

    # {sentinel_id: {dcol, n, label, colour, gids}}
    cat_aggregates: Dict[int, Dict] = {}
    cat_labels = {
        "OECD":    "Other OECD groups",
        "telomer": "Other telomer groups",
        "generic": "Other generic groups",
        "other":   "Other groups",
    }
    for cat, gids in cat_hidden.items():
        if not gids:
            continue
        sid = _CAT_SENTINEL[cat]
        total_n  = sum(group_data[g]["n_compounds"] for g in gids)
        cat_aggregates[sid] = {
            "dcol":   3,
            "n":      total_n,
            "label":  f"{cat_labels[cat]}<br>({len(gids)} groups, {total_n:,} cpds)",
            "colour": _CATEGORY_COLOUR.get(cat, _GREY),
            "gids":   gids,
        }

    # ─── Virtual nodes ───────────────────────────────────────────────────────
    virtual: Dict[int, Dict] = {
        _FLUORODB_ID: {
            "dcol":   0,
            "n":      n_fluorodb_est,
            "label":  f"{dataset_name}<br>(~{n_fluorodb_est:,} est.)",
            "colour": _LGREY,
        },
        _POLY_ONLY_ID: {
            "dcol":   2,
            "n":      n_poly_only,
            "label":  f"Polyfluoroalkyl only<br>({n_poly_only:,})",
            "colour": _lerp(_C1, "#FFFFFF", 0.4),
        },
    }
    virtual.update(cat_aggregates)

    # ─── Node enumeration ────────────────────────────────────────────────────
    def _dcol(gid: int) -> int:
        if gid in virtual:
            return virtual[gid]["dcol"]
        return col_of[gid] + 1   # real: col 0→dcol 1, col 1→dcol 2, col 2→dcol 3 …

    real_ordered = sorted(
        all_real_shown,
        key=lambda g: (col_of[g], -group_data[g]["n_compounds"]),
    )
    all_node_ids: List[int] = sorted(virtual.keys()) + real_ordered
    node_idx: Dict[int, int] = {nid: i for i, nid in enumerate(all_node_ids)}

    # ─── Group nodes by display column (for height/width computation) ───────
    by_dcol: Dict[int, List[int]] = defaultdict(list)
    for gid in all_node_ids:
        by_dcol[_dcol(gid)].append(gid)
    max_dcol = max(by_dcol.keys()) if by_dcol else 3

    def _n_of(gid: int) -> int:
        return virtual[gid]["n"] if gid in virtual else group_data[gid]["n_compounds"]

    for c in by_dcol:
        by_dcol[c].sort(key=lambda g: (0 if g in virtual else 1, -_n_of(g)))

    labels, colours = [], []
    for gid in all_node_ids:
        if gid in virtual:
            labels.append(virtual[gid]["label"])
            colours.append(virtual[gid]["colour"])
        else:
            d = group_data[gid]
            labels.append(f"{d['alias']}<br>({d['n_compounds']:,})")
            colours.append(_node_colour(gid))

    # ─── Flows ───────────────────────────────────────────────────────────────
    _ETYPE_COL: Dict[str, str] = {
        "per": _C0, "poly": _C2, "both": _C3, "oneway": _GREY,
    }
    flows: List[Dict] = []

    def _add(src: int, tgt: int, value: int, colour: str, label: str) -> None:
        if src not in node_idx or tgt not in node_idx:
            return
        if value < min_flow:
            return
        flows.append({"src": src, "tgt": tgt, "value": value,
                      "colour": colour, "label": label})

    # Dataset root → each poly backbone
    if _POLY_ID in active:
        _add(_FLUORODB_ID, _POLY_ID, n_poly,
             _rgba(_C2, 0.30), f"{dataset_name} \u2192 polyfluoroalkyl<br>{n_poly:,} compounds")
    if _POLY_ARYL_ID in active:
        _add(_FLUORODB_ID, _POLY_ARYL_ID, n_poly_aryl,
             _rgba(_C0, 0.30), f"{dataset_name} \u2192 Polyfluoroaryl<br>{n_poly_aryl:,} compounds")
    if _POLY_CYCLIC_ID in active:
        _add(_FLUORODB_ID, _POLY_CYCLIC_ID, n_poly_cyclic,
             _rgba(_C3, 0.30), f"{dataset_name} \u2192 Polyfluoro cyclic<br>{n_poly_cyclic:,} compounds")

    # Poly → Perf (within each fluorocarbon class)
    for poly_g, perf_g in _POLY_TO_PERF.items():
        if poly_g not in active or perf_g not in active:
            continue
        n = _ovl(poly_g, perf_g)
        name_poly = group_data[poly_g]["alias"]
        name_perf = group_data[perf_g]["alias"]
        _add(poly_g, perf_g, n,
             _rgba(_C0, 0.35),
             f"{name_poly} → {name_perf}<br>{n:,} fully-fluorinated compounds")

    # polyfluoroalkyl → "Polyfluoroalkyl only"
    if _POLY_ID in active:
        _add(_POLY_ID, _POLY_ONLY_ID, n_poly_only,
             _rgba(_C1, 0.35),
             f"Polyfluoroalkyl not perfluoroalkyl<br>{n_poly_only:,} partially-fluorinated")

    # perfluoroalkyl → shown col-2 functional roots
    if _PERF_ID in active:
        for gid in col2_shown:
            n = _ovl(_PERF_ID, gid)
            _add(_PERF_ID, gid, n,
                 _rgba(_C0, 0.18),
                 f"perfluoroalkyl → {group_data[gid]['alias']}<br>"
                 f"{n:,} compounds in both")

    # perfluoroalkyl → compact category aggregates
    if _PERF_ID in active:
        for sid, agg in cat_aggregates.items():
            n = sum(_ovl(_PERF_ID, g) for g in agg["gids"])
            _add(_PERF_ID, sid, n,
                 _rgba(_CATEGORY_COLOUR.get(next(
                     (cat for cat, s in _CAT_SENTINEL.items() if s == sid), "other"
                 ), _GREY), 0.18),
                 f"perfluoroalkyl → {agg['label']}<br>{n:,} total compounds in both")

    # Rule-based flows within shown functional groups
    _priority = {"per": 4, "poly": 3, "both": 2, "oneway": 1}
    sorted_rules = sorted(
        [r for r in rules
         if r["src_id"] in shown_deep and r["tgt_id"] in shown_deep
         and r["src_id"] not in _BACKBONE],
        key=lambda r: (
            _priority.get(r["edge_type"], 0),
            -overlap.get((r["src_id"], r["tgt_id"]), {}).get("n_shared", 0),
        ),
        reverse=True,
    )
    used_pairs: Set[Tuple[int, int]] = set()
    for r in sorted_rules:
        s, t = r["src_id"], r["tgt_id"]
        if (s, t) in used_pairs:
            continue
        ov_val = overlap.get((s, t), {}).get("n_shared", 0)
        used_pairs.add((s, t))
        base = _ETYPE_COL.get(r["edge_type"], _GREY)
        _add(s, t, ov_val, _rgba(base, 0.35),
             f"{group_data[s]['alias']} → {group_data[t]['alias']}<br>"
             f"{ov_val:,} shared compounds  ({r['edge_type']})")

    # ─── Build link arrays ────────────────────────────────────────────────────
    srcs, tgts, vals, lcols, llabs = [], [], [], [], []
    for f in flows:
        srcs.append(node_idx[f["src"]])
        tgts.append(node_idx[f["tgt"]])
        vals.append(f["value"])
        lcols.append(f["colour"])
        llabs.append(f["label"])

    # ─── Figure ──────────────────────────────────────────────────────────────
    n_col3 = len(by_dcol.get(3, []))
    height = max(800, n_col3 * 30 + 350)
    width  = max(1300, (max_dcol - 2) * 200 + 600)

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            label=labels, color=colours,
            pad=10, thickness=18,
            line=dict(color="white", width=0.8),
            hovertemplate="%{label}<extra></extra>",
        ),
        link=dict(
            source=srcs, target=tgts, value=vals,
            color=lcols, customdata=llabs,
            hovertemplate="%{customdata}<extra></extra>",
        ))
    )
    # Legend
    for leg_col, leg_name in [
        (_C2,   f"{dataset_name} → polyfluoroalkyl"),
        (_C0,   f"{dataset_name} → Polyfluoroaryl"),
        (_C3,   f"{dataset_name} → Polyfluoro cyclic"),
        (_C0,   "poly → perf (same class)"),
        (_C1,   "polyfluoroalkyl only (not perfluoro)"),
    ]:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color=leg_col, width=4),
            name=leg_name, showlegend=True,
        ))
    for et, ec in _ETYPE_COL.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color=ec, width=4),
            name=f"rule: {et}", showlegend=True,
        ))
    for cat, cc in _CATEGORY_COLOUR.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=cc),
            name=f"category: {cat}", showlegend=True,
        ))

    compact_note = (f"  ·  {sum(len(a['gids']) for a in cat_aggregates.values())} "
                    f"groups collapsed into categories"
                    if cat_aggregates else "")
    fig.update_layout(
        title=dict(
            text=(
                f"<b>PFASGroups hierarchy — compound counts from {dataset_name}</b><br>"
                f"<sup>polyfluoroalkyl: {n_poly:,}  ·  "
                f"perfluoroalkyl: {n_perf:,}  ·  "
                f"polyfluoroalkyl only: {n_poly_only:,}  ·  "
                f"Polyfluoroaryl: {n_poly_aryl:,}  ·  "
                f"Polyfluoro cyclic: {n_poly_cyclic:,}  ·  "
                f"{len(all_real_shown)} real groups{compact_note}</sup>"
            ),
            x=0.5, font=dict(size=12, family='Ubuntu'),

        ),
        font=dict(family="Ubuntu, DejaVu Sans, sans-serif", size=10),
        paper_bgcolor="white",
        height=height, width=width,
        margin=dict(l=20, r=250, t=100, b=20),
        showlegend=True,
        legend=dict(
            x=1.01, y=1.0, xanchor="left", bgcolor="white",
            bordercolor=_LGREY, borderwidth=1, font=dict(size=9),
        ),
    )
    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Hierarchical Sankey of PFASGroups compound counts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--overlap",       type=Path, default=_OVERLAP_DEFAULT, metavar="CSV")
    p.add_argument("--rules",          type=Path, default=_RULES_DEFAULT,   metavar="JSON")
    p.add_argument("--dataset-name",   default=None, metavar="NAME",
                   help="Label for the root node and diagram title (default: derived from "
                        "overlap filename, e.g. 'FluoroDB', 'CLinventory', 'OECD')")
    p.add_argument("--min-compounds", type=int, default=100,
                   help="Exclude groups with fewer compounds")
    p.add_argument("--min-root-compounds", type=int, default=50_000,
                   help="Min compounds for a functional root to be shown individually "
                        "(roots with rule-children are always shown; "
                        "hidden roots are aggregated by category in compact mode)")
    p.add_argument("--all", dest="show_all", action="store_true",
                   help="Show all groups (overrides --min-root-compounds to 0). "
                        "Produces a taller, more detailed diagram.")
    p.add_argument("--min-flow", type=int, default=100,
                   help="Exclude Sankey links with fewer shared compounds")
    p.add_argument("--outdir",   type=Path, default=BENCHMARK_DIR / "imgs")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    if not args.overlap or not args.overlap.exists():
        sys.exit("ERROR: pfasgroup_overlap.csv not found. Use --overlap PATH")

    print(f"Loading overlap: {args.overlap}")
    group_data, overlap = load_overlap(args.overlap)
    print(f"  {len(group_data)} groups")

    # Derive dataset name from filename if not given
    dataset_name = args.dataset_name
    if not dataset_name:
        stem = args.overlap.stem.lower()  # e.g. "pfasgroup_overlap" or "clinventory_overlap"
        if "clinventory" in stem or "cl_inventory" in stem:
            dataset_name = "CLinventory"
        elif "oecd" in stem:
            dataset_name = "OECD"
        else:
            # Default overlap file ("pfasgroup_overlap") or any FluoroDB-derived file
            dataset_name = "FluoroDB"
    file_slug = dataset_name.replace(" ", "_").lower()
    print(f"Loading rules:   {args.rules}")
    rules = load_rules(args.rules, group_data)

    print(f"Building columns (min_compounds={args.min_compounds}) ...")
    col_of, children_of = build_columns(group_data, rules, args.min_compounds)
    col_counts = {}
    for c in col_of.values():
        col_counts[c] = col_counts.get(c, 0) + 1
    print("  Nodes per column:", dict(sorted(col_counts.items())))

    min_root = 0 if args.show_all else args.min_root_compounds
    mode = "full (all groups)" if args.show_all else f"compact (threshold={min_root:,})"
    print(f"Building Sankey ({mode}) ...")
    fig = build_sankey(
        group_data, overlap, rules, col_of, children_of,
        min_flow=args.min_flow,
        min_root_compounds=min_root,
        dataset_name=dataset_name,
    )

    args.outdir.mkdir(exist_ok=True, parents=True)
    from datetime import datetime
    ts     = datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = "_all" if args.show_all else ""
    base   = args.outdir / f"group_hierarchy_sankey_{file_slug}{suffix}_{ts}"

    html = base.with_suffix(".html")
    fig.write_html(str(html), include_plotlyjs="cdn")
    print(f"  Saved: {html}")

    for ext, kw in [("png", {"scale": 2}), ("pdf", {})]:
        try:
            p = base.with_suffix(f".{ext}")
            fig.write_image(str(p), **kw)
            print(f"  Saved: {p}")
        except Exception as exc:
            print(f"  {ext.upper()} skipped: {str(exc).split(chr(10))[0][:80]}")

    print("Done.")


if __name__ == "__main__":
    main()
