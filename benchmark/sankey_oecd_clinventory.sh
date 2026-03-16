#!/usr/bin/env bash
# Generate all Sankey diagrams for the OECD / CLinventory benchmark.
# Run from the benchmark/ directory:
#   cd benchmark && bash sankey_oecd_clinventory.sh
#
# Prerequisites — export overlap CSVs from the database (adjust list surnames if needed):
#   cd <zeropmdb>/database
#   python manage.py export_pfasgroup_overlap --list FluoroDB     --output <benchmark>/data/pfasgroup_overlap.csv
#   python manage.py export_pfasgroup_overlap --list CLinventory  --output <benchmark>/data/clinventory_overlap.csv
#   python manage.py export_pfasgroup_overlap --list OECD         --output <benchmark>/data/oecd_overlap.csv
set -uo pipefail

ATLAS_SCRIPT="scripts/sankey_pfasgroups_atlas.py"
HIERARCHY_SCRIPT="scripts/group_hierarchy_sankey.py"

# Adjust if your database list surnames differ:
FLUORODB_OVERLAP="data/pfasgroup_overlap.csv"
CLINVENTORY_OVERLAP="data/clinventory_overlap.csv"
OECD_OVERLAP="data/oecd_overlap.csv"

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'
BOLD='\033[1m'; RESET='\033[0m'

ok=0; fail=0
declare -a failed_jobs=()

run_job() {
    local desc="$1"; local script="$2"; shift 2
    echo -e "\n${CYAN}${BOLD}▶ ${desc}${RESET}"
    echo -e "  ${BOLD}Command:${RESET} python ${script} $*"
    local t0=$SECONDS
    if python "$script" "$@"; then
        local elapsed=$(( SECONDS - t0 ))
        echo -e "  ${GREEN}✔ Done${RESET}  (${elapsed}s)"
        (( ++ok ))
    else
        local elapsed=$(( SECONDS - t0 ))
        echo -e "  ${RED}✘ FAILED${RESET}  (${elapsed}s)"
        (( ++fail ))
        failed_jobs+=("$desc")
    fi
}

skip_missing() {
    local file="$1" desc="$2"
    if [[ ! -f "$file" ]]; then
        echo -e "\n${RED}${BOLD}⚠ SKIP${RESET} ${desc}"
        echo -e "  Overlap file not found: $file"
        echo -e "  Export it first:  python manage.py export_pfasgroup_overlap --list <SURNAME> --output $file"
        (( ++fail ))
        failed_jobs+=("SKIP: $desc (missing $file)")
        return 1
    fi
    return 0
}

total_t0=$SECONDS

# ── Section 1: Hierarchy Sankeys (backbone → functional groups) ───────────────
echo -e "\n${BOLD}═══ Hierarchy Sankeys ═══${RESET}"

# FluoroDB
skip_missing "$FLUORODB_OVERLAP" "FluoroDB hierarchy" \
  && run_job "FluoroDB — hierarchy Sankey (compact)" "$HIERARCHY_SCRIPT" \
       --overlap "$FLUORODB_OVERLAP"

skip_missing "$FLUORODB_OVERLAP" "FluoroDB hierarchy (all groups)" \
  && run_job "FluoroDB — hierarchy Sankey (all groups)" "$HIERARCHY_SCRIPT" \
       --overlap "$FLUORODB_OVERLAP" --all

# CLinventory
skip_missing "$CLINVENTORY_OVERLAP" "CLinventory hierarchy" \
  && run_job "CLinventory — hierarchy Sankey (compact)" "$HIERARCHY_SCRIPT" \
       --overlap "$CLINVENTORY_OVERLAP"

skip_missing "$CLINVENTORY_OVERLAP" "CLinventory hierarchy (all groups)" \
  && run_job "CLinventory — hierarchy Sankey (all groups)" "$HIERARCHY_SCRIPT" \
       --overlap "$CLINVENTORY_OVERLAP" --all

# OECD
skip_missing "$OECD_OVERLAP" "OECD hierarchy" \
  && run_job "OECD — hierarchy Sankey (compact)" "$HIERARCHY_SCRIPT" \
       --overlap "$OECD_OVERLAP"

skip_missing "$OECD_OVERLAP" "OECD hierarchy (all groups)" \
  && run_job "OECD — hierarchy Sankey (all groups)" "$HIERARCHY_SCRIPT" \
       --overlap "$OECD_OVERLAP" --all

# ── Section 2: Atlas Sankeys (PFASGroups → PFAS-Atlas categories) ────────────
echo -e "\n${BOLD}═══ Atlas Sankeys ═══${RESET}"

run_job "Groups 1-28 — top-28 labels, exclude non-PFAS" "$ATLAS_SCRIPT" \
    --top-n 28 --groups 1-28 --exclude-not-pfas

run_job "Groups 1-28 — top-12 labels, exclude non-PFAS" "$ATLAS_SCRIPT" \
    --top-n 12 --groups 1-28 --exclude-not-pfas

run_job "Groups 29-115 — top-12 labels, exclude non-PFAS" "$ATLAS_SCRIPT" \
    --top-n 12 --groups 29-115 --exclude-not-pfas

run_job "Groups 29-115 — top-25 labels, exclude non-PFAS" "$ATLAS_SCRIPT" \
    --top-n 25 --groups 29-115 --exclude-not-pfas

run_job "All groups — top-12 labels, exclude non-PFAS" "$ATLAS_SCRIPT" \
    --top-n 12 --exclude-not-pfas

# ── Summary ───────────────────────────────────────────────────────────────────
total_elapsed=$(( SECONDS - total_t0 ))
total=$(( ok + fail ))
echo -e "\n${BOLD}════════════════════════════════════════${RESET}"
echo -e "${BOLD}Summary${RESET}  (${total_elapsed}s total)"
echo -e "  ${GREEN}${ok}/${total} jobs succeeded${RESET}"
if (( fail > 0 )); then
    echo -e "  ${RED}${fail} job(s) failed:${RESET}"
    for j in "${failed_jobs[@]}"; do
        echo -e "    ${RED}✘${RESET}  $j"
    done
    exit 1
fi
echo -e "${BOLD}════════════════════════════════════════${RESET}"