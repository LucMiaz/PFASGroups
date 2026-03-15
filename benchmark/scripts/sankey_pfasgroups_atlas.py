#!/usr/bin/env python3
"""
Sankey diagrams: PFASGroups groups → PFAS-Atlas categories.

Left side:  individual PFASGroups group IDs and names (format: "ID: name")
Right side: PFAS-Atlas atlas_class1 classification

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
import json
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
    from PFASGroups.results_model import _get_group_info, _CLASSIFY_EXCLUDED_IDS
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


# Fixed Atlas-side colours (right nodes)
_ATLAS_FIXED: Dict[str, str] = {
    "Not PFAS":                      _GREY,
    "Other PFASs":                   _C0,
    "PFAA precursors":               _C1,
    "Other PFASs, cyclic":           _lerp_color(_C0, "#ffffff", 0.40),
    "Polyfluoroalkyl acids":         _C2,
    "Polyfluoroalkyl acids, cyclic": _lerp_color(_C2, "#ffffff", 0.40),
    "PFAA precursors, cyclic":       _lerp_color(_C1, "#ffffff", 0.40),
    "PFAAs":                         _C3,
    "PFAAs, cyclic":                 _lerp_color(_C3, "#ffffff", 0.40),
    "error":                         _LGREY,
    "unavailable":                   _LGREY,
    "Unknown":                       _LGREY,
}

# ---------------------------------------------------------------------------
# PFASGroups metadata + group extraction
# ---------------------------------------------------------------------------
_GROUP_INFO = _get_group_info()

_NOT_CLASSIFIED = "Not classified"
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


def _group_label(gid: int, name: str) -> str:
    """Format a group as 'ID: name'."""
    return f"{gid}: {name}"


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
            cached = {r["id"]: r.get("pg_groups", []) for r in json.load(fh)}
        n_found = sum(1 for r in molecules if str(r.get("id")) in cached)
        print(f"  Cache hit: {n_found:,}/{len(molecules):,} molecules")
        enriched = []
        for r in molecules:
            groups = cached.get(str(r.get("id")), [])
            enriched.append({**r, "pg_groups": groups})
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
            groups = _pg_groups(embedding)
        except Exception:
            groups = []

        results.append({
            **mol_info,
            "pg_groups": [{"gid": gid, "name": name} for gid, name in groups],
        })

    print()
    elapsed_total = time.perf_counter() - t0
    print(f"  Done in {elapsed_total:.1f}s  ({n / elapsed_total:.0f} mol/s)")

    # Save cache
    cache_path.parent.mkdir(exist_ok=True, parents=True)
    cache_records = [
        {"id": str(r.get("id", "")), "pg_groups": r.get("pg_groups", [])}
        for r in results
    ]
    with cache_path.open("w", encoding="utf-8") as fh:
        json.dump(cache_records, fh, indent=2)
    print(f"  Cache saved: {cache_path}")
    return results


# ---------------------------------------------------------------------------
# Flow builder
# ---------------------------------------------------------------------------
def build_flows(
    records: List[dict],
    right_key: str = "atlas_class1",
    top_n: int = 20,
) -> List[Tuple[str, str, int]]:
    """
    Returns list of (left_label, right_label, count).

    Left:  'ID: group name'  (or _NOT_CLASSIFIED if no groups matched).
    Right: atlas_class1 value.
    Rare groups (outside top_n by molecule count) are merged into _OTHER_GROUPS.
    """
    raw: Counter = Counter()
    for r in records:
        right = r.get(right_key) or "Unknown"
        if right in (None, "error", "unavailable"):
            right = "error"
        groups = r.get("pg_groups", [])
        if not groups:
            raw[(_NOT_CLASSIFIED, right)] += 1
        else:
            for g in groups:
                label = _group_label(g["gid"], g["name"])
                raw[(label, right)] += 1

    # Identify top-N left labels by total count
    left_totals: Counter = Counter()
    for (left, _right), cnt in raw.items():
        left_totals[left] += cnt

    top_labels = {label for label, _ in left_totals.most_common(top_n)}
    top_labels.add(_NOT_CLASSIFIED)

    collapsed: Dict[Tuple[str, str], int] = defaultdict(int)
    for (left, right), cnt in raw.items():
        eff_left = left if left in top_labels else _OTHER_GROUPS
        collapsed[(eff_left, right)] += cnt

    return [(l, r, c) for (l, r), c in sorted(collapsed.items())]


# ---------------------------------------------------------------------------
# Color assignment
# ---------------------------------------------------------------------------
def assign_colors(
    left_labels: List[str],
    right_labels: List[str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (left_color_map, right_color_map) keyed by label string."""
    # Right side — Atlas fixed palette
    right_colors: Dict[str, str] = {
        lab: _ATLAS_FIXED.get(lab, _C0) for lab in right_labels
    }

    # Left side — separate by group category
    oecd_labels:    List[str] = []
    generic_labels: List[str] = []
    telomer_labels: List[str] = []

    for lab in left_labels:
        if lab in (_NOT_CLASSIFIED, _OTHER_GROUPS):
            continue
        m = re.match(r'^(\d+):', lab)
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
    left_colors[_NOT_CLASSIFIED] = _GREY
    left_colors[_OTHER_GROUPS]   = _lerp_color(_C0, "#ffffff", 0.30)

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

    left_labels_ordered  = [l for l, _ in left_totals.most_common()]
    right_labels_ordered = [r for r, _ in right_totals.most_common()]

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

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="white", width=0.5),
            label=all_labels,
            color=node_colors,
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
        height=820,
        width=1200,
        margin=dict(l=20, r=20, t=110, b=20),
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

    # ── Generate Sankeys ──────────────────────────────────────────────────
    args.outdir.mkdir(exist_ok=True, parents=True)
    print()

    for ds_name, ds_records in by_dataset.items():
        safe_ds = ds_name.replace(" ", "_").replace("(", "").replace(")", "")
        print(f"{'='*70}")
        print(f"Dataset: {ds_name}  ({len(ds_records):,} molecules)")
        print(f"{'='*70}")

        flows = build_flows(ds_records, right_key="atlas_class1", top_n=args.top_n)
        title = f"PFASGroups → PFAS-Atlas  |  {ds_name}"
        print_flows(flows, title)

        stem = f"sankey_{safe_ds}_{ts}"
        out_base = args.outdir / stem
        print(f"\n  Writing Sankey: {stem}.*")
        plot_sankey(flows, title, out_base, min_flow=args.min_flow)
        print()

    print("Done.")


if __name__ == "__main__":
    main()
