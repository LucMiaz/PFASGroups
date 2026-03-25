#!/usr/bin/env python3
"""
Post-process benchmark outputs to ensure metadata files are discoverable.

This script keeps dataset structures unchanged and focuses on:
- Copying legacy benchmark outputs from benchmark/scripts/data to benchmark/data
- Adding minimal metadata to definition benchmark files if missing
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path


def copy_if_newer(src: Path, dest: Path) -> bool:
    if not dest.exists():
        shutil.copy2(src, dest)
        return True
    if src.stat().st_mtime > dest.stat().st_mtime:
        shutil.copy2(src, dest)
        return True
    return False


def ensure_definitions_metadata(path: Path) -> bool:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except Exception:
        return False

    if not isinstance(data, dict):
        return False

    if "metadata" in data:
        return False

    data["metadata"] = {
        "timestamp": datetime.now().isoformat(),
        "source": "add_test_metadata.py",
        "note": "metadata added to existing definitions benchmark file",
    }

    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)

    return True


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    benchmark_dir = script_dir.parents[1]
    data_dir = benchmark_dir / "data"
    legacy_dir = script_dir / "data"

    data_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    if legacy_dir.exists():
        for src in legacy_dir.glob("*.json"):
            if not (src.name.startswith("pfas_") or "telomer_validation" in src.name):
                continue
            dest = data_dir / src.name
            if copy_if_newer(src, dest):
                copied += 1

    updated_meta = 0
    for defs_file in data_dir.glob("pfas_definitions_benchmark_*.json"):
        if ensure_definitions_metadata(defs_file):
            updated_meta += 1

    print("Metadata post-processing complete")
    print(f"  Legacy files copied: {copied}")
    print(f"  Definitions metadata updated: {updated_meta}")


if __name__ == "__main__":
    main()
