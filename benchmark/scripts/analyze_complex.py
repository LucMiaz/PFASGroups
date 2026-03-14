#!/usr/bin/env python3
"""
Analyze complex branched benchmark results.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path


def load_json(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def analyze_complex(data):
    total = 0
    pfas_detected = 0
    atlas_detected = 0
    pfas_correct = 0
    case_stats = []

    for test_case in data:
        molecules = test_case.get("molecules", []) or []
        total_case = len(molecules)
        total += total_case

        case_pfas = sum(1 for mol in molecules if mol.get("PFASGroups_detected") is True)
        case_atlas = sum(1 for mol in molecules if mol.get("atlas_detected") is True)
        case_correct = sum(1 for mol in molecules if mol.get("PFASGroups_correct") is True)

        pfas_detected += case_pfas
        atlas_detected += case_atlas
        pfas_correct += case_correct

        case_stats.append({
            "test_name": test_case.get("test_name"),
            "description": test_case.get("description"),
            "complexity": test_case.get("complexity"),
            "molecules": total_case,
            "PFASGroups_detected": case_pfas,
            "atlas_detected": case_atlas,
            "PFASGroups_correct": case_correct,
            "PFASGroups_rate": (case_pfas / total_case) * 100 if total_case else 0.0,
            "atlas_rate": (case_atlas / total_case) * 100 if total_case else 0.0,
            "PFASGroups_correct_rate": (case_correct / total_case) * 100 if total_case else 0.0,
        })

    summary = {
        "total_molecules": total,
        "PFASGroups_detected": pfas_detected,
        "atlas_detected": atlas_detected,
        "PFASGroups_correct": pfas_correct,
        "PFASGroups_detection_rate": (pfas_detected / total) * 100 if total else 0.0,
        "atlas_detection_rate": (atlas_detected / total) * 100 if total else 0.0,
        "PFASGroups_correct_rate": (pfas_correct / total) * 100 if total else 0.0,
    }

    return summary, case_stats


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python analyze_complex.py <complex_branched_benchmark.json>")
        return 1

    input_path = Path(sys.argv[1]).resolve()
    if not input_path.exists():
        print(f"Error: file not found: {input_path}")
        return 1

    data = load_json(input_path)
    if not isinstance(data, list):
        print("Error: expected a list of test cases")
        return 1

    summary, case_stats = analyze_complex(data)

    output_dir = Path("reports")
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = output_dir / f"complex_branched_analysis_{timestamp}.json"

    payload = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "source_file": str(input_path),
            "analysis_type": "complex_branched",
        },
        "summary": summary,
        "cases": case_stats,
    }

    with output_file.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    print("Complex branched analysis complete")
    print(f"  Total molecules: {summary['total_molecules']}")
    print(f"  PFASGroups detection rate: {summary['PFASGroups_detection_rate']:.1f}%")
    print(f"  Atlas detection rate: {summary['atlas_detection_rate']:.1f}%")
    print(f"  Output: {output_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
