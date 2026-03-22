#!/usr/bin/env bash
# run_toxcast_workflow.sh
# ──────────────────────────────────────────────────────────────────────────
# Run the full ToxCast fingerprint comparison workflow.
#
# Steps:
#   1. (Optional) Rebuild the ToxCast dataset parquet from MariaDB.
#   2. Run the fingerprint comparison: nested CV, Random Forests,
#      Gradient Boosting, PFASGroups vs Richard2023.
#   3. Print a summary of output files.
#
# Usage:
#   bash run_toxcast_workflow.sh              # full run
#   bash run_toxcast_workflow.sh --skip-build # skip dataset rebuild
#   CONDA_ENV=myenv bash run_toxcast_workflow.sh
# ──────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
CONDA_ENV="${CONDA_ENV:-chem}"
SKIP_BUILD=0

for arg in "$@"; do
    case "$arg" in
        --skip-build) SKIP_BUILD=1 ;;
        *) echo "[WARN] Unknown argument: $arg" ;;
    esac
done

# ---------------------------------------------------------------------------
# Locate benchmark directory (works whether script is sourced or run directly)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
PARQUET="$DATA_DIR/toxcast_dataset.parquet"

echo ""
echo "Benchmark dir : $SCRIPT_DIR"
echo "Conda env     : $CONDA_ENV"

# ---------------------------------------------------------------------------
# Set up conda
# ---------------------------------------------------------------------------
# Support both Anaconda/Miniconda and Miniforge layouts
for _conda_sh in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh"; do
    if [ -f "$_conda_sh" ]; then
        # shellcheck source=/dev/null
        source "$_conda_sh"
        break
    fi
done

conda activate "$CONDA_ENV"

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
run_step() {
    local label="$1"
    local script="$2"
    echo ""
    echo "=== $label ==="
    local t0
    t0=$(date +%s)
    python "$script"
    local code=$?
    local t1
    t1=$(date +%s)
    local elapsed=$(( t1 - t0 ))
    if [ $code -ne 0 ]; then
        echo "[ERROR] '$script' exited with code $code" >&2
        exit $code
    fi
    printf "[OK] finished in %dm%02ds\n" $(( elapsed / 60 )) $(( elapsed % 60 ))
}

# ---------------------------------------------------------------------------
# Step 1 – dataset build
# ---------------------------------------------------------------------------
if [ "$SKIP_BUILD" -eq 1 ]; then
    echo ""
    echo "=== Skipping dataset build (--skip-build) ==="
    if [ ! -f "$PARQUET" ]; then
        echo "[WARN] $PARQUET not found — comparison step may fail" >&2
    fi
elif [ -f "$PARQUET" ]; then
    echo ""
    echo "Dataset already exists: $PARQUET"
    read -r -p "Rebuild it? [y/N] " answer
    if [[ "$answer" =~ ^[Yy] ]]; then
        run_step "Building ToxCast dataset" "$SCRIPT_DIR/scripts/build_toxcast_dataset.py"
    else
        echo "=== Using existing dataset ==="
    fi
else
    run_step "Building ToxCast dataset" "$SCRIPT_DIR/scripts/build_toxcast_dataset.py"
fi

# ---------------------------------------------------------------------------
# Step 2 – fingerprint comparison
# ---------------------------------------------------------------------------
run_step "Running fingerprint comparison (nested CV)" \
         "$SCRIPT_DIR/scripts/compare_fingerprints_toxcast.py"

# ---------------------------------------------------------------------------
# Step 3 – summary
# ---------------------------------------------------------------------------
echo ""
echo "=== Output files ==="
if [ -d "$DATA_DIR" ]; then
    find "$DATA_DIR" -maxdepth 1 -name "toxcast_comparison_*" | sort | while read -r f; do
        size=$(du -h "$f" 2>/dev/null | cut -f1)
        echo "  $(basename "$f")  ($size)"
    done
else
    echo "[WARN] data directory not found: $DATA_DIR" >&2
fi

echo ""
echo "Workflow complete."
