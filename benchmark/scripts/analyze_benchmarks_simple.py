#!/usr/bin/env python3
"""
Generate a lightweight benchmark summary for reports/benchmark_summary.json.
"""

from __future__ import annotations

import json
import statistics
from datetime import datetime
from pathlib import Path
import glob


def get_latest_file(pattern: str) -> str | None:
    files = sorted(glob.glob(pattern), key=lambda x: Path(x).stat().st_mtime, reverse=True)
    return files[0] if files else None


def load_json(path: str):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def summarize_timing(data):
    if not data:
        return {}
    times = [row.get("HalogenGroups_time_avg") for row in data if row.get("HalogenGroups_time_avg") is not None]
    if not times:
        return {"n_molecules": len(data)}
    return {
        "n_molecules": len(data),
        "mean_time_ms": statistics.mean(times) * 1000,
        "median_time_ms": statistics.median(times) * 1000,
        "min_time_ms": min(times) * 1000,
        "max_time_ms": max(times) * 1000,
    }


def summarize_enhanced(data):
    if not data:
        return {}
    detected = sum(1 for row in data if row.get("HalogenGroups_result", {}).get("detected_groups"))
    return {
        "total": len(data),
        "detected": detected,
        "detection_rate": (detected / len(data)) * 100 if data else 0.0,
    }


def summarize_complex(data):
    if not data:
        return {}
    total = 0
    pfas_detected = 0
    atlas_detected = 0
    for test_case in data:
        molecules = test_case.get("molecules", []) or []
        total += len(molecules)
        pfas_detected += sum(1 for mol in molecules if mol.get("HalogenGroups_detected") is True)
        atlas_detected += sum(1 for mol in molecules if mol.get("atlas_detected") is True)
    return {
        "total": total,
        "HalogenGroups_detected": pfas_detected,
        "atlas_detected": atlas_detected,
        "HalogenGroups_detection_rate": (pfas_detected / total) * 100 if total else 0.0,
        "atlas_detection_rate": (atlas_detected / total) * 100 if total else 0.0,
    }


def summarize_highly_branched(data):
    details = []
    if isinstance(data, dict):
        details = data.get("details", []) or []
    elif isinstance(data, list):
        details = data

    if not details:
        return {}

    total = len(details)
    passed = sum(1 for row in details if row.get("HalogenGroups_passed") is True)
    atlas_passed = sum(1 for row in details if row.get("atlas_passed") is True)

    return {
        "total": total,
        "HalogenGroups_passed": passed,
        "atlas_passed": atlas_passed,
        "HalogenGroups_pass_rate": (passed / total) * 100 if total else 0.0,
        "atlas_pass_rate": (atlas_passed / total) * 100 if total else 0.0,
    }


def summarize_telomer(data):
    if not data:
        return {}

    if isinstance(data, dict):
        total = data.get("total_molecules", 0)
        detected = data.get("telomer_detected", 0)
        return {
            "total": total,
            "detected": detected,
            "detection_rate": (detected / total) * 100 if total else 0.0,
        }

    if isinstance(data, list):
        total = len(data)
        detected = sum(1 for row in data if row.get("detected") is True)
        return {
            "total": total,
            "detected": detected,
            "detection_rate": (detected / total) * 100 if total else 0.0,
        }

    return {}


def main() -> int:
    data_dir = Path("data")
    reports_dir = Path("reports")
    reports_dir.mkdir(parents=True, exist_ok=True)

    timing_file = get_latest_file(str(data_dir / "pfas_timing_benchmark_*.json"))
    enhanced_file = get_latest_file(str(data_dir / "pfas_enhanced_benchmark_*.json"))
    complex_file = get_latest_file(str(data_dir / "pfas_complex_branched_benchmark_*.json"))
    highly_file = get_latest_file(str(data_dir / "pfas_highly_branched_benchmark_*.json"))
    telomer_file = data_dir / "telomer_validation_results.json"

    timing_data = load_json(timing_file) if timing_file else []
    enhanced_data = load_json(enhanced_file) if enhanced_file else []
    complex_data = load_json(complex_file) if complex_file else []
    highly_data = load_json(highly_file) if highly_file else []
    telomer_data = load_json(str(telomer_file)) if telomer_file.exists() else []

    summary = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "timing_file": timing_file,
            "enhanced_file": enhanced_file,
            "complex_file": complex_file,
            "highly_file": highly_file,
            "telomer_file": str(telomer_file) if telomer_file.exists() else None,
        },
        "timing": summarize_timing(timing_data),
        "enhanced": summarize_enhanced(enhanced_data),
        "complex_branched": summarize_complex(complex_data),
        "highly_branched": summarize_highly_branched(highly_data),
        "telomer_validation": summarize_telomer(telomer_data),
    }

    output_file = reports_dir / "benchmark_summary.json"
    with output_file.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Summary saved to: {output_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
