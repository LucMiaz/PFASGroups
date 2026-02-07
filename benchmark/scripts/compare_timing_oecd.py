#!/usr/bin/env python3
"""
Compare timing benchmark vs OECD benchmark timings.

Creates a side-by-side comparison for PFASGroups and PFAS-Atlas:
- Timing benchmark (synthetic coverage over wide size range)
- OECD benchmark (real dataset)
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - RDKit should be available in this repo
    Chem = None


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def latest_file(pattern: str) -> Optional[Path]:
    files = sorted(Path("benchmark/data").glob(pattern), key=lambda x: x.stat().st_mtime, reverse=True)
    return files[0] if files else None


def atoms_from_smiles(smiles: str) -> Optional[int]:
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol.GetNumAtoms()


def extract_timing_benchmark(path: Path) -> Dict[str, List[Tuple[int, float]]]:
    data = load_json(path)
    pfas = []
    atlas = []
    for row in data:
        atoms = row.get("num_atoms")
        if atoms is None:
            continue
        if "pfasgroups_time_avg" in row:
            p_time = row.get("pfasgroups_time_avg")
            a_time = row.get("atlas_time_avg")
        else:
            p_time = row.get("pfasgroups_time")
            a_time = row.get("atlas_time")
        if p_time is not None:
            pfas.append((atoms, p_time * 1000))
        if a_time is not None:
            atlas.append((atoms, a_time * 1000))
    return {"pfasgroups": pfas, "atlas": atlas}


def extract_oecd_benchmark(path: Path) -> Dict[str, List[Tuple[int, float]]]:
    data = load_json(path)
    pfas = []
    atlas = []
    for row in data:
        smiles = row.get("molecule_data", {}).get("smiles")
        atoms = atoms_from_smiles(smiles) if smiles else None
        if atoms is None:
            continue
        p_time = row.get("pfasgroups_result", {}).get("execution_time")
        a_time = row.get("atlas_result", {}).get("execution_time")
        if p_time is not None:
            pfas.append((atoms, p_time * 1000))
        if a_time is not None:
            atlas.append((atoms, a_time * 1000))
    return {"pfasgroups": pfas, "atlas": atlas}


def percentile_band(data: List[Tuple[int, float]], percentiles=(10, 50, 90)):
    if not data:
        return None
    atoms = np.array([d[0] for d in data])
    times = np.array([d[1] for d in data])
    # Bin by atom counts (20 bins over observed range)
    min_atoms = atoms.min()
    max_atoms = atoms.max()
    if min_atoms == max_atoms:
        return None
    bins = np.linspace(min_atoms, max_atoms, 21)
    centers = 0.5 * (bins[:-1] + bins[1:])
    band = {"x": [], "p10": [], "p50": [], "p90": []}
    for i in range(len(bins) - 1):
        mask = (atoms >= bins[i]) & (atoms < bins[i + 1])
        if not np.any(mask):
            continue
        values = times[mask]
        band["x"].append(float(centers[i]))
        band["p10"].append(float(np.percentile(values, percentiles[0])))
        band["p50"].append(float(np.percentile(values, percentiles[1])))
        band["p90"].append(float(np.percentile(values, percentiles[2])))
    return band


def add_dataset_traces(fig, row, col, timing_data, oecd_data, title):
    timing_atoms = [x for x, _ in timing_data]
    timing_times = [y for _, y in timing_data]
    oecd_atoms = [x for x, _ in oecd_data]
    oecd_times = [y for _, y in oecd_data]

    fig.add_trace(
        go.Scatter(
            x=timing_atoms,
            y=timing_times,
            mode="markers",
            name="Timing benchmark",
            marker=dict(color="#1f77b4", size=6, opacity=0.5),
            hovertemplate="atoms=%{x}, time=%{y:.2f}ms<extra>Timing benchmark</extra>",
        ),
        row=row,
        col=col,
    )

    fig.add_trace(
        go.Scatter(
            x=oecd_atoms,
            y=oecd_times,
            mode="markers",
            name="OECD benchmark",
            marker=dict(color="#ff7f0e", size=6, opacity=0.6),
            hovertemplate="atoms=%{x}, time=%{y:.2f}ms<extra>OECD benchmark</extra>",
        ),
        row=row,
        col=col,
    )

    band = percentile_band(timing_data)
    if band:
        fig.add_trace(
            go.Scatter(
                x=band["x"],
                y=band["p50"],
                mode="lines",
                name="Timing median",
                line=dict(color="#1f77b4", width=2),
                hovertemplate="atoms=%{x:.0f}, median=%{y:.2f}ms<extra></extra>",
                showlegend=False,
            ),
            row=row,
            col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=band["x"] + band["x"][::-1],
                y=band["p90"] + band["p10"][::-1],
                fill="toself",
                fillcolor="rgba(31,119,180,0.15)",
                line=dict(color="rgba(255,255,255,0)"),
                name="Timing P10-P90",
                hoverinfo="skip",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    fig.update_xaxes(title_text="Number of atoms", row=row, col=col)
    fig.update_yaxes(title_text="Execution time (ms)", row=row, col=col)
    fig.update_layout(annotations=fig.layout.annotations + (dict(text=title, x=0.5, y=1.08, xref="x domain", yref="y domain", showarrow=False),))


def main() -> int:
    timing_file = latest_file("pfas_timing_benchmark_*.json")
    oecd_file = latest_file("pfas_oecd_benchmark_*.json")

    if timing_file is None:
        print("❌ No timing benchmark file found in benchmark/data")
        return 1
    if oecd_file is None:
        print("❌ No OECD benchmark file found in benchmark/data")
        return 1

    timing_data = extract_timing_benchmark(timing_file)
    oecd_data = extract_oecd_benchmark(oecd_file)

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("PFASGroups", "PFAS-Atlas"),
        horizontal_spacing=0.08,
    )

    add_dataset_traces(
        fig,
        1,
        1,
        timing_data["pfasgroups"],
        oecd_data["pfasgroups"],
        "PFASGroups",
    )
    add_dataset_traces(
        fig,
        1,
        2,
        timing_data["atlas"],
        oecd_data["atlas"],
        "PFAS-Atlas",
    )

    fig.update_layout(
        title=(
            "Timing Benchmark vs OECD Benchmark\n"
            "<sub>Theoretical timing range vs real dataset timing</sub>"
        ),
        width=1400,
        height=650,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    reports_dir = Path("benchmark/reports")
    imgs_dir = Path("benchmark/imgs")
    reports_dir.mkdir(parents=True, exist_ok=True)
    imgs_dir.mkdir(parents=True, exist_ok=True)

    html_path = reports_dir / f"timing_oecd_comparison_{timestamp}.html"
    png_path = imgs_dir / f"timing_oecd_comparison_{timestamp}.png"

    fig.write_html(html_path)
    fig.write_image(png_path, width=1400, height=650)

    print("Comparison report generated:")
    print(f"  HTML: {html_path}")
    print(f"  PNG:  {png_path}")
    print(f"  Timing data: {timing_file}")
    print(f"  OECD data:   {oecd_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
