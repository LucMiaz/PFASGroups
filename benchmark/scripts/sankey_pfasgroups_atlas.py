#!/usr/bin/env python3
"""
Sankey diagrams: PFASGroups groups → PFAS-Atlas categories.

Left side:  individual PFASGroups group IDs and names (format: "ID: name")
Right side: PFAS-Atlas atlas_class2 classification (falls back to atlas_class1)

Each side is labelled with the algorithm name.
Colors are shades of the project palette (color_scheme.yaml) + grey for
"not classified" / "error" nodes.

Usage (from benchmark/ directory):
    python scripts/sankey_pfasgroups_atlas.py
    python scripts/sankey_pfasgroups_atlas.py --input data/oecd_clinventory_timing_20260315_062059.json
    python scripts/sankey_pfasgroups_atlas.py --sample 5000           # quick test
    python scripts/sankey_pfasgroups_atlas.py --no-cache              # force re-classify
    python scripts/sankey_pfasgroups_atlas.py --dataset OECD          # only one dataset
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent
REPO_ROOT     = BENCHMARK_DIR.parent

for _p in (str(REPO_ROOT),):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError as _e:
    sys.exit(f"ERROR: rdkit not available: {_e}")

try:
    from PFASGroups import parse_mol
    from PFASGroups.core import rdkit_disable_log
    from PFASGroups.PFASEmbeddings import _get_group_info, _CLASSIFY_EXCLUDED_IDS
    rdkit_disable_log()
except ImportError as _e:
    sys.exit(f"ERROR: PFASGroups not available: {_e}")

try:
    import plotly.graph_objects as go
    PLOTLY = True
except ImportError:
    PLOTLY = False
    print("WARNING: plotly not available — Sankey output will be text only")

try:
    import numpy as np
    NUMPY = True
except ImportError:
    NUMPY = False

# ---------------------------------------------------------------------------
# Colour palette  (loaded from color_scheme.yaml)
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


_PALETTE = _load_palette()
_C0, _C1, _C2, _C3 = _PALETTE  # orange, blue, magenta, dark-purple

_GREY  = "#AAAAAA"
_LGREY = "#DDDDDD"


def _hex_to_rgb(h: str) -> Tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _rgb_to_hex(r: int, g: int, b: int) -> str:
    return f"#{r:02X}{g:02X}{b:02X}"


def _lerp_color(c1: str, c2: str, t: float) -> str:
    r1, g1, b1 = _hex_to_rgb(c1)
    r2, g2, b2 = _hex_to_rgb(c2)
    return _rgb_to_hex(
        int(r1 + (r2 - r1) * t),
        int(g1 + (g2 - g1) * t),
        int(b1 + (b2 - b1) * t),
    )


def _hex_to_rgba(hexcol: str, alpha: float = 0.5) -> str:
    r, g, b = _hex_to_rgb(hexcol)
    return f"rgba({r},{g},{b},{alpha:.2f})"


def _palette_for_n(n: int, c_start: str, c_end: str) -> List[str]:
    """Generate n colours interpolated from c_start to c_end."""
    if n <= 1:
        return [c_start]
    return [_lerp_color(c_start, c_end, i / (n - 1)) for i in range(n)]


# Fixed Atlas-side colours (right nodes) — covers both class1 and class2 values
_ATLAS_FIXED: Dict[str, str] = {
    # ── class1 values ───────────────────────────────────────────────────────
    "Not PFAS":                             _GREY,
    "Other PFASs":                          _C0,
    "Other PFASs, cyclic":                  _lerp_color(_C0, "#ffffff", 0.40),
    "PFAA precursors":                      _C1,
    "PFAA precursors, cyclic":              _lerp_color(_C1, "#ffffff", 0.40),
    "Polyfluoroalkyl acids":                _C2,
    "Polyfluoroalkyl acids, cyclic":        _lerp_color(_C2, "#ffffff", 0.40),
    "PFAAs":                                _C3,
    "PFAAs, cyclic":                        _lerp_color(_C3, "#ffffff", 0.40),
    "Complex structure":                    _lerp_color(_GREY, "#ffffff", 0.30),
    # ── class2 — PFAAs subtypes (dark-purple → blue shades) ─────────────────
    "PFCAs":                                _C3,
    "PFSAs":                                _lerp_color(_C3, _C1, 0.20),
    "PFSIAs":                               _lerp_color(_C3, _C1, 0.35),
    "PFECAs":                               _lerp_color(_C3, _C1, 0.50),
    "PFESAs":                               _lerp_color(_C3, _C1, 0.60),
    "PFPAs":                                _lerp_color(_C3, _C1, 0.70),
    "PFPIAs":                               _lerp_color(_C3, _C1, 0.80),
    "PFdiCAs":                              _lerp_color(_C3, "#ffffff", 0.30),
    "PFdiSAs":                              _lerp_color(_C3, "#ffffff", 0.45),
    "PFCA-anhydrides":                      _lerp_color(_C3, "#ffffff", 0.55),
    "PFCA-ester derivatives":               _lerp_color(_C3, "#ffffff", 0.65),
    "PFSA derivatives":                     _lerp_color(_C3, "#ffffff", 0.70),
    # ── class2 — PFAA precursors subtypes (blue shades) ─────────────────────
    "PASF-based substances":                _C1,
    "PASFs":                                _lerp_color(_C1, "#ffffff", 0.20),
    "PACFs":                                _lerp_color(_C1, "#ffffff", 0.30),
    "PFAIs":                                _lerp_color(_C1, "#ffffff", 0.40),
    "PFALs":                                _lerp_color(_C1, "#ffffff", 0.48),
    "SFAs":                                 _lerp_color(_C1, "#ffffff", 0.55),
    "HFCs":                                 _lerp_color(_C1, _C3, 0.20),
    "HFEs":                                 _lerp_color(_C1, _C3, 0.35),
    "HFOs":                                 _lerp_color(_C1, _C3, 0.50),
    "PFAenes":                              _lerp_color(_C1, _C3, 0.60),
    "PAECFs":                               _lerp_color(_C1, _C3, 0.70),
    "PFACs":                                _lerp_color(_C1, _C3, 0.80),
    "PFAKs":                                _lerp_color(_C1, "#2DAA5F", 0.30),
    "PolyFEAenes":                          _lerp_color(_C1, "#2DAA5F", 0.50),
    "PolyFEACs":                            _lerp_color(_C1, "#2DAA5F", 0.70),
    "Sulfonyl chloride":                    _lerp_color(_C1, "#888888", 0.40),
    "Acid chloride":                        _lerp_color(_C1, "#888888", 0.55),
    "n:2 fluorotelomer-based substances":   _lerp_color(_C1, "#ffffff", 0.60),
    "n:1 FTOHs":                            _lerp_color(_C1, "#ffffff", 0.68),
    "SFAenes":                              _lerp_color(_C1, "#ffffff", 0.74),
    "PolyFAenes":                           _lerp_color(_C1, "#ffffff", 0.78),
    "PolyFACs":                             _lerp_color(_C1, "#ffffff", 0.82),
    "PolyFAC derivatives":                  _lerp_color(_C1, "#ffffff", 0.86),
    # ── class2 — Polyfluoroalkyl acids subtypes (magenta shades) ────────────
    "PolyFCAs":                             _C2,
    "PolyFSAs":                             _lerp_color(_C2, _C1, 0.30),
    "PolyFECAs":                            _lerp_color(_C2, _C1, 0.55),
    "PolyFESAs":                            _lerp_color(_C2, _C3, 0.30),
    "PolyFSA derivatives":                  _lerp_color(_C2, "#ffffff", 0.35),
    "PolyFCA derivatives":                  _lerp_color(_C2, "#ffffff", 0.50),
    # ── class2 — Other PFASs subtypes (orange shades) ───────────────────────
    "PFAS derivatives":                     _C0,
    "Aromatic PFASs":                       _lerp_color(_C0, "#ffffff", 0.25),
    "side-chain fluorinated aromatics":     _lerp_color(_C0, "#ffffff", 0.40),
    "Perfluoroalkylethers":                 _lerp_color(_C0, "#ffffff", 0.50),
    "Amide derivatives":                    _lerp_color(_C0, "#ffffff", 0.58),
    "Perfluoroalkanes":                     _lerp_color(_C0, "#ffffff", 0.65),
    "Polyfluoroalkanes":                    _lerp_color(_C0, "#ffffff", 0.72),
    "Si PFASs":                             _lerp_color(_C0, _C2, 0.40),
    "(Hg,Sn,Ge,Sb,Se,B) PFASs":           _lerp_color(_C0, _C2, 0.60),
    "others":                               _lerp_color(_C0, "#ffffff", 0.80),
    # ── non-cyclic entries seen in practice but not listed above ─────────────
    "PFAK derivatives":                     _lerp_color(_C1, "#2DAA5F", 0.40),
    "PFPEs":                                _lerp_color(_C0, "#ffffff", 0.55),
    "Perfluoroalkyl-tert-amines":           _lerp_color(_C0, "#ffffff", 0.75),
    # ── fallback / error nodes ───────────────────────────────────────────────
    "Not PFAS by current definition":       _GREY,
    "error":                                _LGREY,
    "unavailable":                          _LGREY,
    "Unknown":                              _LGREY,
}

# ---------------------------------------------------------------------------
# PFASGroups metadata + group extraction
# ---------------------------------------------------------------------------
_GROUP_INFO = _get_group_info()

_OTHER_GROUPS   = "Other PFAS groups"


def _group_category(gid: int) -> str:
    """Return category string ('OECD', 'generic', 'telomer', 'other') for a group id."""
    return _GROUP_INFO.get(gid, {}).get("category", "other")


def _pg_groups(embedding) -> List[Tuple[int, str]]:
    """Return deduplicated list of (group_id, name) for all matched groups.

    Groups in _CLASSIFY_EXCLUDED_IDS (51, 52 — generic catch-alls) are skipped.
    """
    seen: set = set()
    result: List[Tuple[int, str]] = []
    for m in embedding.iter_group_matches():
        gid = m.group_id
        if gid is None or gid in seen or gid in _CLASSIFY_EXCLUDED_IDS:
            continue
        seen.add(gid)
        name = m.group_name or _GROUP_INFO.get(gid, {}).get("name", f"group_{gid}")
        result.append((gid, name))
    return result


def _pg_has_generic(embedding) -> bool:
    """Return True if a generic catch-all group (51 or 52) matched.

    These groups indicate a poly/per-fluorinated component was detected but
    no more specific PFASGroups group was assigned.
    """
    for m in embedding.iter_group_matches():
        if m.group_id in _CLASSIFY_EXCLUDED_IDS:
            return True
    return False


def _group_label(gid: int, name: str) -> str:
    """Format a group as 'name (id XX)'."""
    return f"{name} (id {gid})"


# ---------------------------------------------------------------------------
# Classification loop
# ---------------------------------------------------------------------------
@rdkit_disable_log()
def classify_records(
    molecules: List[dict],
    cache_path: Path,
    force: bool = False,
) -> List[dict]:
    """
    For each molecule record add pg_groups: list of {gid, name} dicts.

    Uses / writes a JSON cache at cache_path unless force=True.
    """
    if not force and cache_path.exists():
        print(f"  Loading cached groups from {cache_path.name} …")
        with cache_path.open(encoding="utf-8") as fh:
            raw_records = json.load(fh)
        # Detect stale cache: if any entry is missing pg_has_generic the cache
        # was built before that field was added and must be regenerated.
        if any("pg_has_generic" not in rec for rec in raw_records):
            print("  Cache is stale (missing pg_has_generic) — re-classifying …")
        else:
            cached = {
                rec["id"]: {
                    "pg_groups":      rec.get("pg_groups", []),
                    "pg_has_generic": rec["pg_has_generic"],
                }
                for rec in raw_records
            }
            n_found = sum(1 for r in molecules if str(r.get("id")) in cached)
            print(f"  Cache hit: {n_found:,}/{len(molecules):,} molecules")
            enriched = []
            for r in molecules:
                info = cached.get(str(r.get("id")), {"pg_groups": [], "pg_has_generic": False})
                enriched.append({**r, "pg_groups": info["pg_groups"],
                                 "pg_has_generic": info["pg_has_generic"]})
            return enriched

    print(f"  Classifying {len(molecules):,} molecules (no cache / --no-cache) …")
    n = len(molecules)
    report_every = max(1, n // 20)
    results = []
    t0 = time.perf_counter()

    for i, mol_info in enumerate(molecules, 1):
        if i % report_every == 0 or i == n:
            elapsed = time.perf_counter() - t0
            rate = i / elapsed if elapsed > 0 else 0
            eta = (n - i) / rate if rate > 0 else 0
            print(f"    {i:,}/{n:,}  ({rate:.0f} mol/s, ETA {eta:.0f}s) …", end="\r")

        smiles = mol_info.get("smiles", "")
        if not smiles:
            results.append({**mol_info, "pg_groups": []})
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({**mol_info, "pg_groups": []})
            continue

        try:
            from rdkit import RDLogger as _RL
            _RL.DisableLog("rdApp.*")
            embedding = parse_mol(mol, include_PFAS_definitions=True)
            groups      = _pg_groups(embedding)
            has_generic = _pg_has_generic(embedding)
        except Exception:
            groups      = []
            has_generic = False

        results.append({
            **mol_info,
            "pg_groups":      [{"gid": gid, "name": name} for gid, name in groups],
            "pg_has_generic": has_generic,
        })

    print()
    elapsed_total = time.perf_counter() - t0
    print(f"  Done in {elapsed_total:.1f}s  ({n / elapsed_total:.0f} mol/s)")

    # Save cache
    cache_path.parent.mkdir(exist_ok=True, parents=True)
    cache_records = [
        {
            "id":            str(r.get("id", "")),
            "pg_groups":     r.get("pg_groups", []),
            "pg_has_generic": r.get("pg_has_generic", False),
        }
        for r in results
    ]
    with cache_path.open("w", encoding="utf-8") as fh:
        json.dump(cache_records, fh, indent=2)
    print(f"  Cache saved: {cache_path}")
    return results


# ---------------------------------------------------------------------------
# Flow builder
# ---------------------------------------------------------------------------
_NOT_PFAS_LABELS = {"Not PFAS", "Not PFAS by current definition"}

_NOT_CLASSIFIED_GENERIC = "Not classified (fluorinated component found)"
_NOT_CLASSIFIED_NONE    = "Not classified (no fluorinated component)"


def build_flows(
    records: List[dict],
    right_key: str = "atlas_class2",
    top_n: int = 20,
    group_ids: Optional[set] = None,
    exclude_not_pfas: bool = False,
) -> List[Tuple[str, str, int]]:
    """
    Returns list of (left_label, right_label, count).

    Left:  'name (id XX)'  (or _NOT_CLASSIFIED_GENERIC / _NOT_CLASSIFIED_NONE if no groups matched).
    Right: atlas_class2 value (falls back to atlas_class1 when class2 is absent).
    Rare groups (outside top_n by molecule count) are merged into _OTHER_GROUPS.
    If group_ids is given, only groups whose id is in that set are kept;
    molecules with none of those groups are split into the two Not-classified buckets.
    If exclude_not_pfas is True, molecules that are both unmatched by PFASGroups
    and classified as non-PFAS by Atlas are silently dropped.
    """
    raw: Counter = Counter()
    for r in records:
        right = r.get(right_key) or r.get("atlas_class1") or "Unknown"
        if right in (None, "error", "unavailable"):
            right = "error"
        groups = r.get("pg_groups", [])
        if group_ids is not None:
            groups = [g for g in groups if g["gid"] in group_ids]
        if not groups:
            if exclude_not_pfas and right in _NOT_PFAS_LABELS:
                continue
            has_generic = r.get("pg_has_generic", False)
            nc_label = _NOT_CLASSIFIED_GENERIC if has_generic else _NOT_CLASSIFIED_NONE
            raw[(nc_label, right)] += 1
        else:
            for g in groups:
                label = _group_label(g["gid"], g["name"])
                raw[(label, right)] += 1

    # Identify top-N left labels by total count
    left_totals: Counter = Counter()
    for (left, _right), cnt in raw.items():
        left_totals[left] += cnt

    _NC_LABELS = {_NOT_CLASSIFIED_GENERIC, _NOT_CLASSIFIED_NONE}
    top_labels = {label for label, _ in left_totals.most_common(top_n)}
    top_labels |= _NC_LABELS  # always keep both Not-classified buckets

    collapsed: Dict[Tuple[str, str], int] = defaultdict(int)
    for (left, right), cnt in raw.items():
        eff_left = left if left in top_labels else _OTHER_GROUPS
        collapsed[(eff_left, right)] += cnt

    # Sort flows: specials (Not classified, Other PFAS groups) last, then by count desc
    _SPECIAL_LEFT = {_NOT_CLASSIFIED_GENERIC, _NOT_CLASSIFIED_NONE, _OTHER_GROUPS}
    left_totals_collapsed: Counter = Counter()
    for (left, _), cnt in collapsed.items():
        left_totals_collapsed[left] += cnt

    return sorted(
        [(l, r, c) for (l, r), c in collapsed.items()],
        key=lambda x: (1 if x[0] in _SPECIAL_LEFT else 0, -left_totals_collapsed[x[0]], x[1]),
    )


# ---------------------------------------------------------------------------
# Color assignment
# ---------------------------------------------------------------------------
def _atlas_node_color(label: str) -> str:
    """Look up a color for an Atlas class2 node label.

    If the label is not in ``_ATLAS_FIXED`` but ends with '', cyclic', the
    color is derived automatically from the non-cyclic base (lightened toward
    white by 35 %).
    """
    if label in _ATLAS_FIXED:
        return _ATLAS_FIXED[label]
    if label.endswith(", cyclic"):
        base = label[: -len(", cyclic")]
        base_color = _ATLAS_FIXED.get(base, _C0)
        return _lerp_color(base_color, "#ffffff", 0.35)
    return _C0


def assign_colors(
    left_labels: List[str],
    right_labels: List[str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (left_color_map, right_color_map) keyed by label string."""
    # Right side — Atlas fixed palette (with auto cyclic fallback)
    right_colors: Dict[str, str] = {
        lab: _atlas_node_color(lab) for lab in right_labels
    }

    # Left side — separate by group category
    oecd_labels:    List[str] = []
    generic_labels: List[str] = []
    telomer_labels: List[str] = []

    _NC_LABELS = {_NOT_CLASSIFIED_GENERIC, _NOT_CLASSIFIED_NONE}
    for lab in left_labels:
        if lab in _NC_LABELS or lab == _OTHER_GROUPS:
            continue
        m = re.search(r'\(id (\d+)\)', lab)
        if not m:
            continue
        gid = int(m.group(1))
        cat = _group_category(gid)
        if cat == "OECD":
            oecd_labels.append(lab)
        elif cat == "generic":
            generic_labels.append(lab)
        elif cat == "telomer":
            telomer_labels.append(lab)
        else:
            generic_labels.append(lab)   # fallback: treat as generic

    left_colors: Dict[str, str] = {}
    # OECD groups: orange → blue
    for lab, col in zip(oecd_labels, _palette_for_n(len(oecd_labels) or 1, _C0, _C1)):
        left_colors[lab] = col
    # Generic groups: magenta → dark-purple
    for lab, col in zip(generic_labels, _palette_for_n(len(generic_labels) or 1, _C2, _C3)):
        left_colors[lab] = col
    # Telomeric groups: blue → teal
    _teal = "#2DAA5F"
    for lab, col in zip(telomer_labels, _palette_for_n(len(telomer_labels) or 1, _lerp_color(_C1, _teal, 0.5), _teal)):
        left_colors[lab] = col
    # Special nodes
    left_colors[_NOT_CLASSIFIED_GENERIC] = _lerp_color(_GREY, _C0, 0.30)   # warm grey — fluorinated found
    left_colors[_NOT_CLASSIFIED_NONE]    = _GREY                             # plain grey — nothing found
    left_colors[_OTHER_GROUPS]           = _lerp_color(_C0, "#ffffff", 0.30)

    return left_colors, right_colors


# ---------------------------------------------------------------------------
# Plotly Sankey
# ---------------------------------------------------------------------------
def plot_sankey(
    flows: List[Tuple[str, str, int]],
    title: str,
    output_base: Path,
    min_flow: int = 5,
) -> None:
    if not PLOTLY:
        print("  [plotly unavailable — skipping]")
        return

    flows = [(l, r, c) for l, r, c in flows if c >= min_flow]
    if not flows:
        print("  No flows above threshold.")
        return

    # Order nodes: left by total descending, right by total descending
    left_totals:  Counter = Counter()
    right_totals: Counter = Counter()
    for l, r, c in flows:
        left_totals[l]  += c
        right_totals[r] += c

    # Left: flow-descending, specials (Not classified, Other PFAS groups) at bottom
    _LEFT_BOTTOM  = {_NOT_CLASSIFIED_GENERIC, _NOT_CLASSIFIED_NONE, _OTHER_GROUPS}
    left_main   = [l for l, _ in left_totals.most_common()  if l not in _LEFT_BOTTOM]
    left_bottom = [l for l, _ in left_totals.most_common()  if l in  _LEFT_BOTTOM]
    left_labels_ordered = left_main + left_bottom

    # Right: flow-descending, non-PFAS / error nodes at bottom
    _RIGHT_BOTTOM = {"Not PFAS", "Not PFAS by current definition",
                     "error", "unavailable", "Unknown"}
    right_main   = [r for r, _ in right_totals.most_common() if r not in _RIGHT_BOTTOM]
    right_bottom = [r for r, _ in right_totals.most_common() if r in  _RIGHT_BOTTOM]
    right_labels_ordered = right_main + right_bottom

    all_labels = left_labels_ordered + right_labels_ordered
    label_idx  = {lab: i for i, lab in enumerate(all_labels)}

    left_colors, right_colors = assign_colors(left_labels_ordered, right_labels_ordered)

    node_colors: List[str] = []
    for lab in left_labels_ordered:
        node_colors.append(left_colors.get(lab, _C0))
    for lab in right_labels_ordered:
        node_colors.append(right_colors.get(lab, _C0))

    sources, targets, values, link_colors = [], [], [], []
    for l_lab, r_lab, cnt in flows:
        if l_lab not in label_idx or r_lab not in label_idx:
            continue
        sources.append(label_idx[l_lab])
        targets.append(label_idx[r_lab])
        values.append(cnt)
        link_colors.append(_hex_to_rgba(left_colors.get(l_lab, _C0), 0.40))

    # ── Layout parameters ────────────────────────────────────────────────
    n_left    = len(left_labels_ordered)
    n_right   = len(right_labels_ordered)
    NODE_THICKNESS = 20
    height_px = max(920, n_left * 50 + 200)
    width_px  = 1640
    margin_t  = 110
    margin_b  = 40
    margin_l  = 20
    margin_r  = 20

    # Full label for each node (name + count), shown by Plotly automatically
    all_totals = {**left_totals, **right_totals}
    full_labels = [
        f"{lab}  ({all_totals.get(lab, 0):,})" for lab in all_labels
    ]

    node_x = [0.001] * n_left + [0.999] * n_right

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=8,
            thickness=NODE_THICKNESS,
            line=dict(color="white", width=0.5),
            label=full_labels,
            color=node_colors,
            x=node_x,
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
        ),
    ))

    n_pairs = sum(c for _, _, c in flows)

    fig.update_layout(
        title_text=(
            f"<b>{title}</b>"
            f"<br><sup>n = {n_pairs:,} molecule–group assignments</sup>"
        ),
        title_x=0.5,
        font=dict(family="Ubuntu, sans-serif", size=12),
        paper_bgcolor="white",
        plot_bgcolor="white",
        height=height_px,
        width=width_px,
        margin=dict(l=margin_l, r=margin_r, t=margin_t, b=margin_b),
        annotations=[
            dict(
                text="<b>PFASGroups</b>",
                x=0.01, y=1.06,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=14, color="#333333"),
                xanchor="left",
            ),
            dict(
                text="<b>PFAS-Atlas</b>",
                x=0.99, y=1.06,
                xref="paper", yref="paper",
                showarrow=False,
                font=dict(size=14, color="#333333"),
                xanchor="right",
            ),
        ],
    )

    html_path = output_base.with_suffix(".html")
    fig.write_html(str(html_path), include_plotlyjs="cdn")
    print(f"  Saved: {html_path.relative_to(BENCHMARK_DIR)}")

    for ext, kwargs in [("png", {"scale": 2}), ("pdf", {})]:
        try:
            img_path = output_base.with_suffix(f".{ext}")
            fig.write_image(str(img_path), **kwargs)
            print(f"  Saved: {img_path.relative_to(BENCHMARK_DIR)}")
        except Exception as exc:
            short = str(exc).split("\n")[0][:80]
            print(f"  {ext.upper()} export skipped: {short}")
            print(f"  → Open the .html in a browser to view/export the Sankey")


# ---------------------------------------------------------------------------
# Text summary (always printed)
# ---------------------------------------------------------------------------
def print_flows(flows: List[Tuple[str, str, int]], title: str) -> None:
    print()
    print(f"  {title}")
    print("  " + "-" * 74)
    totals: Counter = Counter()
    for l, r, c in flows:
        totals[l] += c
    grand = sum(totals.values())
    for l_lab in sorted(totals, key=lambda x: -totals[x]):
        sub = [(r, c) for ll, r, c in flows if ll == l_lab]
        sub.sort(key=lambda x: -x[1])
        pct = 100 * totals[l_lab] / grand if grand else 0
        print(f"  {l_lab:<55} n={totals[l_lab]:>7,}  ({pct:5.1f}%)")
        for r_lab, cnt in sub[:5]:
            p2 = 100 * cnt / totals[l_lab]
            print(f"      → {r_lab:<40} {cnt:>6,}  ({p2:.1f}%)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate Sankey diagrams: PFASGroups groups → PFAS-Atlas"
    )
    p.add_argument("--input",    type=Path, default=None,
                   help="Timing JSON file (auto-detects latest if omitted)")
    p.add_argument("--sample",   type=int,  default=None,
                   help="Use only first N molecules (for quick testing)")
    p.add_argument("--no-cache", action="store_true",
                   help="Force re-classification even if cache exists")
    p.add_argument("--top-n",    type=int,  default=20,
                   help="Max distinct left-side group labels per Sankey (default 20)")
    p.add_argument("--min-flow", type=int,  default=5,
                   help="Minimum flow count to include in Sankey (default 5)")
    p.add_argument("--dataset",  default=None,
                   help="Filter to one dataset name, e.g. 'OECD' or 'clinventory'")
    p.add_argument("--outdir",   type=Path,
                   default=BENCHMARK_DIR / "imgs",
                   help="Output directory for figures (default imgs/)")
    p.add_argument("--groups",   default=None,
                   help="Filter to a group-ID range, e.g. '1-28' or '29-115'")
    p.add_argument("--exclude-not-pfas", action="store_true",
                   help="Drop molecules that are unmatched by PFASGroups AND "
                        "classified as non-PFAS by Atlas (removes Not classified → Not PFAS flows)")
    p.add_argument("--right-key", default="atlas_class2",
                   choices=["atlas_class1", "atlas_class2"],
                   help="Atlas field to use for the right-side nodes "
                        "(default: atlas_class2 — detailed; atlas_class1 — broad)")
    p.add_argument("--combine-datasets", action="store_true",
                   help="Combine all datasets into a single Sankey instead of "
                        "one per dataset")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")

    # ── Locate timing JSON ────────────────────────────────────────────────
    input_path = args.input
    if input_path is None:
        candidates = sorted(
            (BENCHMARK_DIR / "data").glob("oecd_clinventory_timing_*.json"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if not candidates:
            sys.exit("ERROR: no oecd_clinventory_timing_*.json found in data/")
        input_path = candidates[0]
    print(f"Input: {input_path}")

    # ── Load molecules ────────────────────────────────────────────────────
    with input_path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    molecules: List[dict] = data.get("molecules", [])
    print(f"Loaded {len(molecules):,} molecules")

    if args.dataset:
        molecules = [m for m in molecules if m.get("dataset") == args.dataset]
        print(f"Filtered to dataset '{args.dataset}': {len(molecules):,} molecules")

    if args.sample:
        molecules = molecules[: args.sample]
        print(f"Sampled: {args.sample:,}")

    # ── Classify: get PFASGroups group matches per molecule ───────────────
    cache_path = BENCHMARK_DIR / "data" / f"sankey_groups_{input_path.stem}.json"
    molecules = classify_records(molecules, cache_path, force=args.no_cache)

    # ── Group by dataset ──────────────────────────────────────────────────
    by_dataset: Dict[str, List[dict]] = defaultdict(list)
    for m in molecules:
        by_dataset[m.get("dataset", "unknown")].append(m)

    print(f"\nDatasets: { {k: len(v) for k, v in by_dataset.items()} }")

    # ── Parse group-ID filter ─────────────────────────────────────────────
    group_ids: Optional[set] = None
    groups_suffix = ""
    if args.groups:
        m = re.match(r'^(\d+)-(\d+)$', args.groups.strip())
        if not m:
            sys.exit(f"ERROR: --groups must be a range like '1-28', got: {args.groups!r}")
        lo, hi = int(m.group(1)), int(m.group(2))
        group_ids = set(range(lo, hi + 1))
        groups_suffix = f"_groups{lo}-{hi}"
        print(f"Group filter: IDs {lo}–{hi} ({len(group_ids)} IDs)")

    right_key    = args.right_key
    class_label  = " (detailed categories)" if right_key == "atlas_class2" else ""
    right_suffix = "_class2" if right_key == "atlas_class2" else "_class1"

    # ── Build dataset iterator ────────────────────────────────────────────
    if args.combine_datasets:
        ds_items: List[Tuple[str, List[dict]]] = [("OECD + CLinventory", molecules)]
    else:
        ds_items = list(by_dataset.items())

    # ── Generate Sankeys ──────────────────────────────────────────────────
    args.outdir.mkdir(exist_ok=True, parents=True)
    print()

    for ds_name, ds_records in ds_items:
        safe_ds = re.sub(r'_+', '_', ds_name.replace(" ", "_").replace("(", "").replace(")", "").replace("+", "")).strip("_")
        print(f"{'='*70}")
        print(f"Dataset: {ds_name}  ({len(ds_records):,} molecules)")
        print(f"{'='*70}")

        flows = build_flows(ds_records, right_key=right_key, top_n=args.top_n,
                            group_ids=group_ids, exclude_not_pfas=args.exclude_not_pfas)
        title = f"PFASGroups → PFAS-Atlas{class_label}  |  {ds_name}"
        if groups_suffix:
            title += f"  |  groups {lo}–{hi}"
        print_flows(flows, title)

        stem = f"sankey_{safe_ds}{groups_suffix}{right_suffix}_{ts}"
        out_base = args.outdir / stem
        print(f"\n  Writing Sankey: {stem}.*")
        plot_sankey(flows, title, out_base, min_flow=args.min_flow)

        # ── CSV of unclassified-by-PFASGroups but PFAS-Atlas-classified ──
        if group_ids is not None:
            _pfas_not = {"Not PFAS", "Unknown", "error", "unavailable", "", None}
            unclassified = [
                r for r in ds_records
                if not any(g["gid"] in group_ids for g in r.get("pg_groups", []))
                and (r.get("atlas_class2") or r.get("atlas_class1")) not in _pfas_not
            ]
            csv_stem = f"unclassified_{safe_ds}{groups_suffix}_{ts}"
            csv_path = args.outdir / f"{csv_stem}.csv"
            csv_fields = [
                "id", "name", "smiles", "formula",
                "atlas_class1", "atlas_class2",
                "pg_all_groups",
            ]
            with csv_path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
                writer.writeheader()
                for r in unclassified:
                    all_grps = "; ".join(
                        f"{g['gid']}: {g['name']}" for g in r.get("pg_groups", [])
                    )
                    writer.writerow({
                        "id":           r.get("id", ""),
                        "name":         r.get("name", ""),
                        "smiles":       r.get("smiles", ""),
                        "formula":      r.get("formula", ""),
                        "atlas_class1": r.get("atlas_class1", ""),
                        "atlas_class2": r.get("atlas_class2", ""),
                        "pg_all_groups": all_grps,
                    })
            print(f"  Unclassified CSV ({len(unclassified):,} molecules): "
                  f"{csv_path.relative_to(BENCHMARK_DIR)}")
        print()

    print("Done.")


if __name__ == "__main__":
    main()
