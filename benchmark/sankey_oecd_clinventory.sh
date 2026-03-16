#!/usr/bin/env bash
# Generate all Sankey diagrams for the OECD / CLinventory benchmark.
# Run from the benchmark/ directory:
#   cd benchmark && bash sankey_oecd_clinventory.sh
set -uo pipefail

SCRIPT="scripts/sankey_pfasgroups_atlas.py"
ENRICH="scripts/enrich_atlas_class2.py"

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'
BOLD='\033[1m'; RESET='\033[0m'

ok=0; fail=0
declare -a failed_jobs=()

run_job() {
    local desc="$1"; shift
    echo -e "\n${CYAN}${BOLD}▶ ${desc}${RESET}"
    echo -e "  ${BOLD}Command:${RESET} python ${SCRIPT} $*"
    local t0=$SECONDS
    if python "$SCRIPT" "$@"; then
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

total_t0=$SECONDS

# ── Enrich timing JSON with atlas_class2 if missing ──────────────────────────
echo -e "\n${CYAN}${BOLD}▶ Enriching timing JSON with atlas_class2 (if needed)${RESET}"
conda run -n pfasatlas python "$ENRICH"

# ── Jobs ──────────────────────────────────────────────────────────────────────
# 1 & 2 — long diagrams: all groups visible, detailed class 2 on the right
run_job "Groups 1-28 — all labels, Atlas class 2 (detailed)" \
    --groups 1-28 --top-n 28 --right-key atlas_class2 --exclude-not-pfas

run_job "Groups 29-115 — all labels, Atlas class 2 (detailed)" \
    --groups 29-115 --top-n 87 --right-key atlas_class2 --exclude-not-pfas

# 3 & 4 — compact diagrams: top-12 groups, broad class 1 on the right
run_job "Groups 1-28 — top-12 labels, Atlas class 1 (broad)" \
    --groups 1-28 --top-n 12 --right-key atlas_class1 --exclude-not-pfas

run_job "Groups 29-115 — top-12 labels, Atlas class 1 (broad)" \
    --groups 29-115 --top-n 12 --right-key atlas_class1 --exclude-not-pfas

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