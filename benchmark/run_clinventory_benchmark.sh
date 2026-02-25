#!/bin/bash
# Clinventory Classification Benchmark Runner
# Classifies all molecules in the clinventory database using HalogenGroups
# and PFAS-Atlas, then compares timing and PFAS-detection performance.
#
# Usage:
#   ./run_clinventory_benchmark.sh [OPTIONS]
#
# Options:
#   --limit N            Process only the first N molecules (for quick tests)
#   --all-molecules      Include non-halogenated molecules
#   --db-password PASS   DB password (prompted if omitted)
#   --skip-classify      Skip classification steps (use existing JSON files)
#   --skip-hg            Skip HalogenGroups classification
#   --skip-atlas         Skip PFAS-Atlas classification
#   --quick              Alias for --limit 1000
#   --help               Show this help

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─── Parse arguments ──────────────────────────────────────────────────────────
DB_PASSWORD="${DATABASE_PASSWORD:-}"
LIMIT=""
ALL_MOL_FLAG=""
SKIP_CLASSIFY=false
SKIP_HG=false
SKIP_ATLAS=false
SKIP_PFASGROUPS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --limit)              LIMIT="$2";       shift 2 ;;
        --limit=*)            LIMIT="${1#--limit=}"; shift ;;
        --all-molecules)      ALL_MOL_FLAG="--all-molecules"; shift ;;
        --db-password)        DB_PASSWORD="$2"; shift 2 ;;
        --db-password=*)      DB_PASSWORD="${1#--db-password=}"; shift ;;
        --skip-classify)      SKIP_CLASSIFY=true; shift ;;
        --skip-hg)            SKIP_HG=true;       shift ;;
        --skip-atlas)         SKIP_ATLAS=true;     shift ;;
        --skip-pfasgroups)    SKIP_PFASGROUPS=true; shift ;;
        --quick)              LIMIT=1000;          shift ;;
        --help)
            echo "Usage: ./run_clinventory_benchmark.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --limit N            Process only the first N molecules"
            echo "  --all-molecules      Include non-halogenated molecules"
            echo "  --db-password PASS   DB password (prompted if omitted)"
            echo "  --skip-classify      Skip both classification steps"
            echo "  --skip-hg            Skip HalogenGroups classification"
            echo "  --skip-atlas         Skip PFAS-Atlas classification"
            echo "  --skip-pfasgroups    Skip PFASGroups (F-only) classification"
            echo "  --quick              Alias for --limit 1000"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "🔬 CLINVENTORY CLASSIFICATION BENCHMARK"
echo "======================================="
if [ -n "$LIMIT" ]; then
    echo "⚡ Limit: ${LIMIT} molecules"
fi
echo ""

# ─── Activate pfasatlas environment ──────────────────────────────────────────
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" != "pfasatlas" ]; then
    if command -v mamba >/dev/null 2>&1; then
        eval "$(mamba shell.bash hook)"
        mamba activate pfasatlas
    elif command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate pfasatlas
    else
        echo "⚠️  mamba/conda not found; please activate 'pfasatlas' manually before running."
        echo "   mamba activate pfasatlas"
        exit 1
    fi
fi
echo "✅ Environment: $CONDA_DEFAULT_ENV"

# ─── Directory check ──────────────────────────────────────────────────────────
if [ ! -f "scripts/classify_halogengroups_clinventory.py" ]; then
    echo "❌ Error: scripts/classify_halogengroups_clinventory.py not found!"
    echo "   Please run this script from the benchmark/ directory."
    exit 1
fi

mkdir -p data imgs reports

# ─── DB password handling ─────────────────────────────────────────────────────
if [ -z "$DB_PASSWORD" ]; then
    read -s -p "Enter clinventory DB password (user 'luc'): " DB_PASSWORD
    echo ""
fi
export DATABASE_PASSWORD="$DB_PASSWORD"
export DATABASE_USER="${DATABASE_USER:-luc}"
export DATABASE_NAME="${DATABASE_NAME:-clinventory}"
export DATABASE_HOST="${DATABASE_HOST:-localhost}"
export DATABASE_PORT="${DATABASE_PORT:-5432}"

# ─── Build optional flags ─────────────────────────────────────────────────────
LIMIT_FLAG=""
[ -n "$LIMIT" ] && LIMIT_FLAG="--limit $LIMIT"

# ─── 1. HalogenGroups classification ─────────────────────────────────────────
if [ "$SKIP_CLASSIFY" = false ] && [ "$SKIP_HG" = false ]; then
    echo ""
    echo "1️⃣  Running HalogenGroups classification …"
    python scripts/classify_halogengroups_clinventory.py \
        --db-password "$DB_PASSWORD" \
        $LIMIT_FLAG $ALL_MOL_FLAG
    if [ $? -ne 0 ]; then
        echo "❌ HalogenGroups classification failed"
        exit 1
    fi
    echo "✅ HalogenGroups classification completed"
else
    echo "⏭️  Skipping HalogenGroups classification"
fi

# ─── 2. PFAS-Atlas classification ─────────────────────────────────────────────
if [ "$SKIP_CLASSIFY" = false ] && [ "$SKIP_ATLAS" = false ]; then
    echo ""
    echo "2️⃣  Running PFAS-Atlas classification …"
    python scripts/classify_pfasatlas_clinventory.py \
        --db-password "$DB_PASSWORD" \
        $LIMIT_FLAG
    if [ $? -ne 0 ]; then
        echo "❌ PFAS-Atlas classification failed"
        exit 1
    fi
    echo "✅ PFAS-Atlas classification completed"
else
    echo "⏭️  Skipping PFAS-Atlas classification"
fi

# ─── 2.5. PFASGroups (F-only) classification ─────────────────────────────────
if [ "$SKIP_CLASSIFY" = false ] && [ "$SKIP_PFASGROUPS" = false ]; then
    echo ""
    echo "2½️ Running PFASGroups (F-only) classification …"
    python scripts/classify_pfasgroups_clinventory.py \
        --db-password "$DB_PASSWORD" \
        $LIMIT_FLAG
    if [ $? -ne 0 ]; then
        echo "❌ PFASGroups classification failed"
        exit 1
    fi
    echo "✅ PFASGroups (F-only) classification completed"
else
    echo "⏭️  Skipping PFASGroups classification"
fi

# ─── 3. Compare classifiers ───────────────────────────────────────────────────
echo ""
echo "3️⃣  Comparing classifiers …"
HG_FILE=$(ls data/halogengroups_clinventory_*.json 2>/dev/null | tail -1)
ATLAS_FILE=$(ls data/pfasatlas_clinventory_*.json  2>/dev/null | tail -1)
PFG_FILE=$(ls data/pfasgroups_clinventory_*.json   2>/dev/null | tail -1)

if [ -z "$HG_FILE" ]; then
    echo "❌ No HalogenGroups result file found in data/"
    exit 1
fi
if [ -z "$ATLAS_FILE" ]; then
    echo "❌ No PFAS-Atlas result file found in data/"
    exit 1
fi

echo "   HalogenGroups : $HG_FILE"
echo "   PFAS-Atlas    : $ATLAS_FILE"
PFG_FLAG=""
if [ -n "$PFG_FILE" ]; then
    echo "   PFASGroups(F) : $PFG_FILE"
    PFG_FLAG="--pfasgroups-file $PFG_FILE"
fi

python scripts/compare_clinventory_classifiers.py \
    --hg-file  "$HG_FILE" \
    --atlas-file "$ATLAS_FILE" \
    ${PFG_FLAG}
if [ $? -ne 0 ]; then
    echo "❌ Comparison failed"
    exit 1
fi
echo "✅ Comparison plots and data generated"

# ─── 4. Generate LaTeX ────────────────────────────────────────────────────────
echo ""
echo "4️⃣  Generating LaTeX tables and text …"
CMP_FILE=$(ls data/clinventory_comparison_*.json 2>/dev/null | tail -1)

if [ -z "$CMP_FILE" ]; then
    echo "❌ No comparison JSON found in data/"
    exit 1
fi

python scripts/generate_clinventory_latex.py \
    --hg-file    "$HG_FILE" \
    --atlas-file "$ATLAS_FILE" \
    --cmp-file   "$CMP_FILE" \
    ${PFG_FLAG}
if [ $? -ne 0 ]; then
    echo "❌ LaTeX generation failed"
    exit 1
fi
echo "✅ LaTeX generated"

# ─── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "🎉 CLINVENTORY BENCHMARK COMPLETE"
echo "================================="
echo ""
echo "📊 Data files (data/):"
ls data/halogengroups_clinventory_*.json 2>/dev/null | xargs -I{} echo "   {}" || true
ls data/pfasatlas_clinventory_*.json     2>/dev/null | xargs -I{} echo "   {}" || true
ls data/pfasgroups_clinventory_*.json    2>/dev/null | xargs -I{} echo "   {}" || true
ls data/clinventory_comparison_*.json   2>/dev/null | xargs -I{} echo "   {}" || true
echo ""
echo "📈 Plots (imgs/):"
ls imgs/clinventory_*.pdf 2>/dev/null | xargs -I{} echo "   {}" || echo "   (none yet)"
echo ""
echo "📝 LaTeX (reports/):"
ls reports/clinventory_comparison_*.tex 2>/dev/null | xargs -I{} echo "   {}" || true
echo ""
echo "Tables generated:"
echo "   • clinventory_comparison_tables.tex  (Tables 1-6)"
echo "   • clinventory_comparison_text.tex    (Narrative section)"
echo "   • clinventory_comparison_figures.tex (\\includegraphics snippets)"
echo ""
echo "Plots generated:"
echo "   • clinventory_timing_box              — overall timing box plots"
echo "   • clinventory_timing_cdf              — timing CDF (log scale)"
echo "   • clinventory_timing_by_bracket       — median timing by atom-count bracket"
echo "   • clinventory_timing_by_halogen       — median timing by halogen class"
echo "   • clinventory_timing_ratio_by_bracket — HG/Atlas timing ratio by bracket"
echo "   • clinventory_classification_agreement — 2×2 agreement bar chart"
echo "   • clinventory_atlas_class_distribution — PFAS-Atlas class histogram"
echo "   • clinventory_hg_top_groups            — top-20 HalogenGroups groups"
