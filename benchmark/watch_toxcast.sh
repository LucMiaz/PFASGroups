#!/usr/bin/env bash
# watch_toxcast.sh
# ─────────────────────────────────────────────────────────────────────────────
# Wait for compare_fingerprints_toxcast.py to finish, then run update_toxcast_si.py.
#
# Usage from the benchmark directory:
#   bash watch_toxcast.sh &
#   # or:
#   bash watch_toxcast.sh 2>&1 | tee watch_toxcast.log
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="${CONDA_ENV:-chem}"
LOCK_FILE="$SCRIPT_DIR/watch_toxcast.lock"
UPDATE_SCRIPT="$SCRIPT_DIR/scripts/update_toxcast_si.py"
RESULTS_CSV="$SCRIPT_DIR/data/toxcast_comparison_results.csv"
SUMMARY_CSV="$SCRIPT_DIR/data/toxcast_comparison_summary.csv"
LOG_FILE="$SCRIPT_DIR/watch_toxcast.log"

# Prevent double-run
if [ -f "$LOCK_FILE" ]; then
    echo "[WARN] Another watcher may be running (lock file exists: $LOCK_FILE). Exiting."
    exit 1
fi
touch "$LOCK_FILE"
trap 'rm -f "$LOCK_FILE"' EXIT

# Setup conda
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

echo ""
echo "$(date): watch_toxcast.sh started. Waiting for compare_fingerprints_toxcast.py to finish..."

# Record modification time of the results CSV before we start waiting
MTIME_BEFORE=0
if [ -f "$RESULTS_CSV" ]; then
    MTIME_BEFORE=$(stat -c %Y "$RESULTS_CSV" 2>/dev/null || echo 0)
fi

# Poll every 60 seconds for the process to end
while true; do
    # Check if the script is still running
    if pgrep -f "compare_fingerprints_toxcast.py" > /dev/null 2>&1; then
        sleep 60
        continue
    fi

    # Process ended — check that the CSV was actually updated (new results written)
    MTIME_AFTER=0
    if [ -f "$RESULTS_CSV" ]; then
        MTIME_AFTER=$(stat -c %Y "$RESULTS_CSV" 2>/dev/null || echo 0)
    fi

    if [ "$MTIME_AFTER" -gt "$MTIME_BEFORE" ] || [ ! -f "$RESULTS_CSV" ]; then
        # CSV was updated (or doesn't exist — run anyway and let the script report errors)
        break
    fi

    # Process ended but CSV not updated — might have crashed; still try
    echo "$(date): Process ended but CSV not updated. Proceeding anyway."
    break
done

echo ""
echo "$(date): compare_fingerprints_toxcast.py finished. Running update_toxcast_si.py ..."
echo ""

python "$UPDATE_SCRIPT" 2>&1 | tee -a "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "$(date): [SUCCESS] update_toxcast_si.py completed."
else
    echo "$(date): [ERROR] update_toxcast_si.py exited with code $EXIT_CODE."
    exit "$EXIT_CODE"
fi
