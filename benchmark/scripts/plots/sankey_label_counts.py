#!/usr/bin/env python3
"""Sankey diagrams from PFASGroups label-assignment counts.

Reads the per-compound classification JSON produced by the timing pipeline
(fields: id, pg_label_oecd, pg_label_generic, pg_label_generic_tel) and
generates one Sankey per label field showing compound distribution across
classification groups.

Layout:
  col 0  Dataset source node (total compound count)
  col 1  Each unique primary classification label (sorted by count desc,
          "Not PFAS" always last)

Outputs HTML + PDF for every label field processed.

Usage (from benchmark/ directory):
    python scripts/sankey_label_counts.py \\
        --timing-json data/sankey_categories_oecd_clinventory_timing_20260315_062059.json

    python scripts/sankey_label_counts.py \\
        --timing-json data/sankey_categories_oecd_clinventory_timing_20260315_062059.json \\
        --label-field oecd --outdir imgs/
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import plotly.graph_objects as go
    import plotly.io as pio
except ImportError:
    sys.exit("ERROR: plotly not available.  pip install plotly")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR    = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[1]

LABEL_FIELDS = {
    "oecd":        "pg_label_oecd",
    "generic":     "pg_label_generic",
    "generic_tel": "pg_label_generic_tel",
}

FIELD_TITLES = {
    "oecd":        "OECD classification",
    "generic":     "Generic groups",
    "generic_tel": "Generic groups (with telomers)",
}

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
_NOT_PFAS = "Not PFAS"

# Source node: teal
_SRC_COL   = "rgba(44, 160, 101, 0.85)"
# PFAS group nodes: blue-green gradient (will cycle)
_PFAS_COLS = [
    "rgba(31, 119, 180, 0.80)",
    "rgba(255, 127, 14, 0.80)",
    "rgba(44, 160, 44, 0.80)",
    "rgba(214, 39, 40, 0.80)",
    "rgba(148, 103, 189, 0.80)",
    "rgba(140, 86, 75, 0.80)",
    "rgba(227, 119, 194, 0.80)",
    "rgba(127, 127, 127, 0.65)",
    "rgba(188, 189, 34, 0.80)",
    "rgba(23, 190, 207, 0.80)",
]
# Not-PFAS sink: grey
_NOT_PFAS_COL = "rgba(160, 160, 160, 0.65)"
# Link colour (semi-transparent version of target colour)
_LINK_ALPHA = 0.35


def _link_colour(hex_or_rgba: str) -> str:
    """Return a more transparent version for links."""
    if hex_or_rgba.startswith("rgba("):
        # Replace last alpha value with _LINK_ALPHA
        parts = hex_or_rgba.rstrip(")").split(",")
        parts[-1] = f" {_LINK_ALPHA})"
        return ",".join(parts)
    return hex_or_rgba


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_timing_json(path: Path) -> List[dict]:
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if isinstance(data, dict):
        # might be wrapped: {"compounds": [...]}
        for key in ("compounds", "data", "results"):
            if key in data and isinstance(data[key], list):
                return data[key]
    if isinstance(data, list):
        return data
    raise ValueError(f"Unexpected JSON structure in {path}")


def primary_label(raw: Optional[str]) -> str:
    """Return the first (highest-priority) label from a comma-separated string."""
    if not raw or not str(raw).strip():
        return "unclassified"
    return str(raw).split(",")[0].strip()


def count_labels(records: List[dict], field: str) -> Counter:
    """Count compounds by primary label for the given field name."""
    c: Counter = Counter()
    for rec in records:
        lbl = primary_label(rec.get(field))
        c[lbl] += 1
    return c


# ---------------------------------------------------------------------------
# Sankey builder
# ---------------------------------------------------------------------------

def build_sankey_figure(
    counts: Counter,
    dataset_name: str,
    field_key: str,
    min_count: int = 1,
) -> go.Figure:
    """Build and return a Sankey figure for one label field."""

    total = sum(counts.values())
    title_field = FIELD_TITLES.get(field_key, field_key)

    # Sort: PFAS groups by count desc, then "Not PFAS" last
    pfas_items = sorted(
        [(lbl, n) for lbl, n in counts.items() if lbl != _NOT_PFAS and n >= min_count],
        key=lambda x: -x[1],
    )
    not_pfas_count = counts.get(_NOT_PFAS, 0)
    if not_pfas_count >= min_count:
        pfas_items.append((_NOT_PFAS, not_pfas_count))

    # Nodes: 0 = source, 1..N = label targets
    node_labels:  List[str] = [f"{dataset_name}<br>({total:,})"]
    node_colours: List[str] = [_SRC_COL]

    link_sources: List[int] = []
    link_targets: List[int] = []
    link_values:  List[int] = []
    link_colours: List[str] = []
    link_labels:  List[str] = []

    for i, (lbl, n) in enumerate(pfas_items):
        node_idx = i + 1
        pct = 100 * n / total if total else 0
        node_labels.append(f"{lbl}<br>({n:,}, {pct:.1f}%)")
        if lbl == _NOT_PFAS:
            col = _NOT_PFAS_COL
        else:
            col = _PFAS_COLS[i % len(_PFAS_COLS)]
        node_colours.append(col)

        link_sources.append(0)
        link_targets.append(node_idx)
        link_values.append(n)
        link_colours.append(_link_colour(col))
        link_labels.append(f"{lbl}: {n:,} ({pct:.1f}%)")

    n_targets = len(pfas_items)
    height = max(600, n_targets * 30 + 200)

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            label=node_labels,
            color=node_colours,
            pad=12,
            thickness=18,
            hovertemplate="%{label}<extra></extra>",
        ),
        link=dict(
            source=link_sources,
            target=link_targets,
            value=link_values,
            color=link_colours,
            label=link_labels,
            hovertemplate="%{label}<extra></extra>",
        ),
    ))

    fig.update_layout(
        title=dict(
            text=(
                f"<b>PFASGroups — {title_field}</b><br>"
                f"<sup>{dataset_name} · {total:,} compounds</sup>"
            ),
            x=0.01,
            font=dict(size=16),
        ),
        font=dict(family="Arial, sans-serif", size=12),
        paper_bgcolor="white",
        height=height,
        width=900,
        margin=dict(l=20, r=20, t=80, b=20),
    )
    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--timing-json",
        required=True,
        metavar="FILE",
        help="Path to the timing JSON with pg_label_* fields",
    )
    p.add_argument(
        "--label-field",
        choices=list(LABEL_FIELDS.keys()),
        default=None,
        metavar="{oecd,generic,generic_tel}",
        help="Which label field to process (default: all three)",
    )
    p.add_argument(
        "--dataset-name",
        default=None,
        metavar="NAME",
        help="Display name for the dataset (default: derived from filename)",
    )
    p.add_argument(
        "--min-count",
        type=int,
        default=1,
        metavar="N",
        help="Minimum compound count to include a label node (default: 1)",
    )
    p.add_argument(
        "--outdir",
        default=None,
        metavar="DIR",
        help="Output directory (default: benchmark/imgs/)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    timing_path = Path(args.timing_json)
    if not timing_path.exists():
        sys.exit(f"ERROR: timing JSON not found: {timing_path}")

    outdir = Path(args.outdir) if args.outdir else BENCHMARK_DIR / "imgs"
    outdir.mkdir(parents=True, exist_ok=True)

    # Derive a dataset name from the filename if not given
    dataset_name = args.dataset_name or timing_path.stem.split("_timing_")[0].replace("_", " ").title()

    print(f"Loading {timing_path} …")
    records = load_timing_json(timing_path)
    print(f"  {len(records):,} records loaded.")

    fields_to_run = (
        {args.label_field: LABEL_FIELDS[args.label_field]}
        if args.label_field
        else LABEL_FIELDS
    )

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    for field_key, field_name in fields_to_run.items():
        counts = count_labels(records, field_name)
        n_labels = sum(1 for n in counts.values() if n >= args.min_count)
        print(f"\n[{field_key}] field='{field_name}' → {n_labels} unique labels")
        for lbl, n in sorted(counts.items(), key=lambda x: -x[1])[:10]:
            pct = 100 * n / len(records) if records else 0
            print(f"  {n:>8,}  ({pct:5.1f}%)  {lbl}")
        if len(counts) > 10:
            print(f"  … and {len(counts)-10} more")

        fig = build_sankey_figure(
            counts=counts,
            dataset_name=dataset_name,
            field_key=field_key,
            min_count=args.min_count,
        )

        stem = f"sankey_labels_{field_key}_{ts}"
        html_path = outdir / f"{stem}.html"
        pdf_path  = outdir / f"{stem}.pdf"

        fig.write_html(str(html_path), include_plotlyjs="cdn")
        print(f"  HTML → {html_path}")

        try:
            fig.write_image(str(pdf_path), format="pdf")
            print(f"  PDF  → {pdf_path}")
        except Exception as exc:
            print(f"  PDF  skipped ({exc}). Install kaleido: pip install kaleido")

    print("\nDone.")


if __name__ == "__main__":
    main()
