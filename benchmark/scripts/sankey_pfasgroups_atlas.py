#!/usr/bin/env python3
"""
Sankey diagrams: PFASGroups categories → PFAS-Atlas categories.

Three versions of the PFASGroups left-side labelling:
  oecd          : OECD-category group names (groups 1–28)
  generic       : generic-category group names (functional-group types, 29–73)
  generic_tel   : generic + telomeric group names, telomers as an explicit bucket

For each version, the right-side labels are the PFAS-Atlas ``atlas_class1``
values stored in the timing-benchmark JSON.

Usage (from benchmark/ directory):
    python scripts/sankey_pfasgroups_atlas.py
    python scripts/sankey_pfasgroups_atlas.py --input data/oecd_clinventory_timing_20260315_010107.json
    python scripts/sankey_pfasgroups_atlas.py --sample 5000           # quick test
    python scripts/sankey_pfasgroups_atlas.py --no-cache              # force re-classify
    python scripts/sankey_pfasgroups_atlas.py --dataset OECD          # only one dataset
    python scripts/sankey_pfasgroups_atlas.py --dataset "clinventory (F)"
"""
from __future__ import annotations

import argparse
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
# Colour palette
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

_NOT_PFAS_COLOR = "#AAAAAA"
_ERROR_COLOR    = "#DDDDDD"


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
    """Generate n colors interpolated from c_start to c_end."""
    if n <= 1:
        return [c_start]
    return [_lerp_color(c_start, c_end, i / (n - 1)) for i in range(n)]


# Fixed Atlas-side colours
_ATLAS_FIXED: Dict[str, str] = {
    "Not PFAS":                      _NOT_PFAS_COLOR,
    "Other PFASs":                   _C0,
    "PFAA precursors":               _C1,
    "Other PFASs, cyclic":           _lerp_color(_C0, "#FFFFFF", 0.35),
    "Polyfluoroalkyl acids":         _C2,
    "Polyfluoroalkyl acids, cyclic": _lerp_color(_C2, "#FFFFFF", 0.35),
    "PFAA precursors, cyclic":       _lerp_color(_C1, "#FFFFFF", 0.35),
    "PFAAs":                         _C3,
    "PFAAs, cyclic":                 _lerp_color(_C3, "#FFFFFF", 0.35),
    "error":                         _ERROR_COLOR,
    "unavailable":                   _ERROR_COLOR,
}

# ---------------------------------------------------------------------------
# PFASGroups category extraction
# ---------------------------------------------------------------------------
_GROUP_INFO = _get_group_info()


def _pg_labels(embedding) -> Dict[str, List[str]]:
    """Extract group names by category type from a PFASEmbedding result."""
    oecd_names:    List[str] = []
    generic_names: List[str] = []
    telomer_names: List[str] = []
    seen_oecd:    set = set()
    seen_generic: set = set()
    seen_telomer: set = set()

    for m in embedding.iter_group_matches():
        gid = m.group_id
        if gid is None:
            continue
        info = _GROUP_INFO.get(gid, {})
        cat  = info.get("category", "other")
        name = m.group_name or info.get("name", f"group_{gid}")

        if cat == "OECD" and name not in seen_oecd:
            seen_oecd.add(name)
            oecd_names.append(name)
        elif cat == "generic" and gid not in _CLASSIFY_EXCLUDED_IDS and name not in seen_generic:
            seen_generic.add(name)
            generic_names.append(name)
        elif cat == "telomer" and name not in seen_telomer:
            seen_telomer.add(name)
            telomer_names.append(name)

    return {"oecd": oecd_names, "generic": generic_names, "telomer": telomer_names}


def _first_or_two(names: List[str]) -> str:
    """Return first name, or 'A, B' if two short names."""
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    if len(names[0]) + len(names[1]) < 40:
        return f"{names[0]}, {names[1]}"
    return names[0]


def label_oecd(labels: Dict[str, List[str]]) -> str:
    """Left label for the OECD version."""
    if labels["oecd"]:
        return _first_or_two(labels["oecd"])
    return "Not PFAS"


def label_generic(labels: Dict[str, List[str]]) -> str:
    """Left label for the generic functional-group version."""
    if labels["generic"]:
        return _first_or_two(labels["generic"])
    return "Not PFAS"


def label_generic_tel(labels: Dict[str, List[str]]) -> str:
    """Left label for the generic+telomers version (telomers explicit)."""
    if labels["telomer"]:
        return "Telomers"
    if labels["generic"]:
        return _first_or_two(labels["generic"])
    return "Not PFAS"


_LABEL_FUNS = {
    "oecd":        label_oecd,
    "generic":     label_generic,
    "generic_tel": label_generic_tel,
}

_VERSION_TITLES = {
    "oecd":        "OECD groups",
    "generic":     "Generic functional groups",
    "generic_tel": "Generic + Telomers",
}

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
    For each molecule record add pg_label_{oecd,generic,generic_tel}.

    Uses / writes a JSON cache at cache_path unless force=True.
    """
    if not force and cache_path.exists():
        print(f"  Loading cached categories from {cache_path.name} …")
        with cache_path.open(encoding="utf-8") as fh:
            cached = {r["id"]: r for r in json.load(fh)}
        n_found = sum(1 for r in molecules if str(r.get("id")) in cached)
        print(f"  Cache hit: {n_found:,}/{len(molecules):,} molecules")
        enriched = []
        for r in molecules:
            c = cached.get(str(r.get("id")))
            if c:
                r = {**r, **{k: c[k] for k in ("pg_label_oecd", "pg_label_generic", "pg_label_generic_tel")}}
            else:
                r = {**r, "pg_label_oecd": "Not PFAS",
                     "pg_label_generic": "Not PFAS",
                     "pg_label_generic_tel": "Not PFAS"}
            enriched.append(r)
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
        base_entry = {
            "id":   str(mol_info.get("id", "")),
            "pg_label_oecd":        "Not PFAS",
            "pg_label_generic":     "Not PFAS",
            "pg_label_generic_tel": "Not PFAS",
        }

        if not smiles:
            results.append({**mol_info, **base_entry})
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({**mol_info, **base_entry})
            continue

        try:
            from rdkit import RDLogger as _RL
            _RL.DisableLog("rdApp.*")
            embedding = parse_mol(mol, include_PFAS_definitions=True)
            labels = _pg_labels(embedding)
        except Exception:
            results.append({**mol_info, **base_entry})
            continue

        entry = {
            "id":                   str(mol_info.get("id", "")),
            "pg_label_oecd":        label_oecd(labels),
            "pg_label_generic":     label_generic(labels),
            "pg_label_generic_tel": label_generic_tel(labels),
        }
        results.append({**mol_info, **entry})

    print()
    elapsed_total = time.perf_counter() - t0
    print(f"  Done in {elapsed_total:.1f}s  ({n / elapsed_total:.0f} mol/s)")

    # Save cache
    cache_path.parent.mkdir(exist_ok=True, parents=True)
    cache_records = [{k: r[k] for k in ("id", "pg_label_oecd", "pg_label_generic", "pg_label_generic_tel")} for r in results]
    with cache_path.open("w", encoding="utf-8") as fh:
        json.dump(cache_records, fh, indent=2)
    print(f"  Cache saved: {cache_path}")
    return results


# ---------------------------------------------------------------------------
# Flow builder
# ---------------------------------------------------------------------------
def build_flows(
    records: List[dict],
    left_key: str,
    right_key: str = "atlas_class1",
    top_n: int = 12,
    other_pfas_label: str = "Other PFAS",
) -> List[Tuple[str, str, int]]:
    """
    Returns list of (left_label, right_label, count), aggregating rare left
    labels into other_pfas_label (preserving "Not PFAS").
    """
    raw: Counter = Counter()
    for r in records:
        left  = r.get(left_key) or "Not PFAS"
        right = r.get(right_key) or "Unknown"
        if right in (None, "error", "unavailable"):
            right = "error"
        raw[(left, right)] += 1

    # Identify top-N left labels by total count (keeping "Not PFAS" always)
    left_totals: Counter = Counter()
    for (left, right), cnt in raw.items():
        left_totals[left] += cnt

    top_labels = {label for label, _ in left_totals.most_common(top_n)}
    top_labels.add("Not PFAS")

    # Build collapsed flows
    collapsed: Dict[Tuple[str, str], int] = defaultdict(int)
    for (left, right), cnt in raw.items():
        eff_left = left if left in top_labels else other_pfas_label
        collapsed[(eff_left, right)] += cnt

    return [(l, r, c) for (l, r), c in sorted(collapsed.items())]


# ---------------------------------------------------------------------------
# Color assignment
# ---------------------------------------------------------------------------
def assign_colors(
    left_labels: List[str],
    right_labels: List[str],
    version: str,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (left_color_map, right_color_map)."""
    left_colors: Dict[str, str] = {}
    right_colors: Dict[str, str] = {}

    # Right side — Atlas fixed palette
    for lab in right_labels:
        right_colors[lab] = _ATLAS_FIXED.get(lab, _C0)

    # Left side
    pfas_labels = [l for l in left_labels if l != "Not PFAS"]
    n_pfas = len(pfas_labels)

    if version == "oecd":
        palette = _palette_for_n(max(n_pfas, 1), _C0, _C1)
    elif version == "generic":
        palette = _palette_for_n(max(n_pfas, 1), _C2, _C3)
    else:  # generic_tel
        palette = _palette_for_n(max(n_pfas, 1), "#2DAA5F", _C1)

    for i, lab in enumerate(pfas_labels):
        left_colors[lab] = palette[i]
    left_colors["Not PFAS"] = _NOT_PFAS_COLOR

    # "Other PFAS" bucket
    for lab in left_labels:
        if lab not in left_colors:
            left_colors[lab] = _lerp_color(_C0, "#FFFFFF", 0.25)

    return left_colors, right_colors


# ---------------------------------------------------------------------------
# Plotly Sankey
# ---------------------------------------------------------------------------
def plot_sankey(
    flows: List[Tuple[str, str, int]],
    title: str,
    version: str,
    output_base: Path,
    min_flow: int = 5,
) -> None:
    if not PLOTLY:
        print("  [plotly unavailable — skipping]")
        return

    # Filter tiny flows for readability
    flows = [(l, r, c) for l, r, c in flows if c >= min_flow]
    if not flows:
        print("  No flows above threshold.")
        return

    left_labels_ordered  = list(dict.fromkeys(l for l, _, _ in flows))
    right_labels_ordered = list(dict.fromkeys(r for _, r, _ in flows))

    all_labels = left_labels_ordered + right_labels_ordered
    label_idx  = {lab: i for i, lab in enumerate(all_labels)}

    left_colors, right_colors = assign_colors(
        left_labels_ordered, right_labels_ordered, version
    )

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
        link_colors.append(_hex_to_rgba(left_colors.get(l_lab, _C0), 0.42))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=18,
            thickness=22,
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

    total = sum(c for _, _, c in flows)
    fig.update_layout(
        title_text=f"<b>{title}</b><br><sup>n = {total:,} molecules</sup>",
        title_x=0.5,
        font=dict(family="Ubuntu, sans-serif", size=12),
        paper_bgcolor="white",
        plot_bgcolor="white",
        height=720,
        width=1100,
        margin=dict(l=30, r=30, t=80, b=30),
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
            print(f"  {ext.upper()} export skipped (kaleido/display issue): {short}")
            print(f"  → Open the .html file in a browser to view/export the Sankey")


# ---------------------------------------------------------------------------
# Text summary (fallback / always printed)
# ---------------------------------------------------------------------------
def print_flows(flows: List[Tuple[str, str, int]], title: str) -> None:
    print()
    print(f"  {title}")
    print("  " + "-" * 70)
    totals: Counter = Counter()
    for l, r, c in flows:
        totals[l] += c
    grand = sum(totals.values())
    for l_lab in sorted(totals, key=lambda x: -totals[x]):
        sub = [(r, c) for ll, r, c in flows if ll == l_lab]
        sub.sort(key=lambda x: -x[1])
        pct = 100 * totals[l_lab] / grand if grand else 0
        print(f"  {l_lab:<45} n={totals[l_lab]:>7,}  ({pct:5.1f}%)")
        for r_lab, cnt in sub[:5]:
            p2 = 100 * cnt / totals[l_lab]
            print(f"      → {r_lab:<40} {cnt:>6,} ({p2:.1f}%)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate Sankey diagrams: PFASGroups categories → PFAS-Atlas"
    )
    p.add_argument("--input",    type=Path, default=None,
                   help="Timing JSON file (auto-detects latest if omitted)")
    p.add_argument("--sample",   type=int,  default=None,
                   help="Use only first N molecules (for quick testing)")
    p.add_argument("--no-cache", action="store_true",
                   help="Force re-classification even if cache exists")
    p.add_argument("--top-n",    type=int,  default=12,
                   help="Max distinct left-side labels per Sankey (default 12)")
    p.add_argument("--min-flow", type=int,  default=5,
                   help="Minimum flow count to include in Sankey (default 5)")
    p.add_argument("--dataset",  default=None,
                   help="Filter to one dataset name, e.g. 'OECD' or 'clinventory (F)'")
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

    # Optionally filter by dataset name
    if args.dataset:
        molecules = [m for m in molecules if m.get("dataset") == args.dataset]
        print(f"Filtered to dataset '{args.dataset}': {len(molecules):,} molecules")

    # F-only subset for clinventory
    aug: List[dict] = []
    for m in molecules:
        if m.get("dataset") == "clinventory" and m.get("n_fluorine", 0) > 0:
            aug.append({**m, "dataset": "clinventory (F)"})
    molecules = molecules + aug
    print(f"Total records (incl. clinventory-F subset): {len(molecules):,}")

    if args.sample:
        molecules = molecules[: args.sample]
        print(f"Sampled: {args.sample:,}")

    # ── Re-classify for PFASGroups categories ────────────────────────────
    cache_path = (BENCHMARK_DIR / "data" /
                  f"sankey_categories_{input_path.stem}.json")
    molecules = classify_records(molecules, cache_path, force=args.no_cache)

    # ── Filter to F-subset after enrichment for re-use ───────────────────
    # (clinventory (F) records are duplicated entries — ok since cache keyed by id+dataset)
    # Re-label F subset
    augmented: List[dict] = []
    for m in molecules:
        if m.get("dataset") == "clinventory" and m.get("n_fluorine", 0) > 0:
            augmented.append({**m, "dataset": "clinventory (F)"})
    molecules = molecules + augmented

    # De-duplicate (keep first occurrence per id+dataset)
    seen: set = set()
    deduped: List[dict] = []
    for m in molecules:
        key = (m.get("id"), m.get("dataset"))
        if key not in seen:
            seen.add(key)
            deduped.append(m)
    molecules = deduped

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

        for version, left_key in [
            ("oecd",        "pg_label_oecd"),
            ("generic",     "pg_label_generic"),
            ("generic_tel", "pg_label_generic_tel"),
        ]:
            other_label = f"Other PFAS ({_VERSION_TITLES[version]})"
            flows = build_flows(
                ds_records,
                left_key=left_key,
                right_key="atlas_class1",
                top_n=args.top_n,
                other_pfas_label=other_label,
            )
            title = (
                f"PFASGroups [{_VERSION_TITLES[version]}] → PFAS-Atlas  |  "
                f"{ds_name}"
            )
            print_flows(flows, title)

            stem = f"sankey_{safe_ds}_{version}_{ts}"
            out_base = args.outdir / stem
            print(f"\n  Writing Sankey: {stem}.*")
            plot_sankey(flows, title, version, out_base, min_flow=args.min_flow)
            print()

    print("Done.")


if __name__ == "__main__":
    main()
