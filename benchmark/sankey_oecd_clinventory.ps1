#Requires -Version 5.1
<#
.SYNOPSIS
    Generate all Sankey diagrams for the OECD / CLinventory benchmark.
.DESCRIPTION
    Run from the benchmark/ directory:
        cd benchmark; .\sankey_oecd_clinventory.ps1

    Prerequisites — export overlap CSVs from the database (adjust list surnames if needed):
        cd <zeropmdb>\database
        python manage.py export_pfasgroup_overlap --list FluoroDB    --output <benchmark>\data\pfasgroup_overlap.csv
        python manage.py export_pfasgroup_overlap --list CLinventory --output <benchmark>\data\clinventory_overlap.csv
        python manage.py export_pfasgroup_overlap --list OECD        --output <benchmark>\data\oecd_overlap.csv
#>
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$AtlasScript     = 'scripts/plots/sankey_pfasgroups_atlas.py'
$HierarchyScript = 'scripts/plots/group_hierarchy_sankey.py'

# Overlap CSV paths (adjust list surnames in export command if needed):
$FluoroDBOverlap     = 'data/pfasgroup_overlap.csv'
$CLinventoryOverlap  = 'data/clinventory_overlap.csv'
$OecdOverlap         = 'data/oecd_overlap.csv'

$ok         = 0
$fail       = 0
$failedJobs = [System.Collections.Generic.List[string]]::new()

function Invoke-SankeyJob {
    param(
        [string]   $Description,
        [string]   $ScriptPath,
        [string[]] $JobArgs
    )
    Write-Host ""
    Write-Host "▶ $Description" -ForegroundColor Cyan
    Write-Host "  Command: python $ScriptPath $($JobArgs -join ' ')"
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    try {
        python $ScriptPath @JobArgs
        if ($LASTEXITCODE -ne 0) { throw "exit code $LASTEXITCODE" }
        $sw.Stop()
        Write-Host "  ✔ Done  ($([int]$sw.Elapsed.TotalSeconds)s)" -ForegroundColor Green
        $script:ok++
    }
    catch {
        $sw.Stop()
        Write-Host "  ✘ FAILED  ($([int]$sw.Elapsed.TotalSeconds)s): $_" -ForegroundColor Red
        $script:fail++
        $script:failedJobs.Add($Description)
    }
}

function Skip-Missing {
    param([string] $FilePath, [string] $Description)
    if (-not (Test-Path $FilePath)) {
        Write-Host ""
        Write-Host "⚠ SKIP  $Description" -ForegroundColor Red
        Write-Host "  Overlap file not found: $FilePath"
        Write-Host "  Export it first:  python manage.py export_pfasgroup_overlap --list <SURNAME> --output $FilePath"
        $script:fail++
        $script:failedJobs.Add("SKIP: $Description (missing $FilePath)")
        return $false
    }
    return $true
}

$totalSw = [System.Diagnostics.Stopwatch]::StartNew()

# ── Section 1: Hierarchy Sankeys (backbone → functional groups) ───────────────
Write-Host "`n═══ Hierarchy Sankeys ═══" -ForegroundColor Cyan

# FluoroDB
if (Skip-Missing $FluoroDBOverlap 'FluoroDB hierarchy') {
    Invoke-SankeyJob 'FluoroDB — hierarchy Sankey (compact)' $HierarchyScript @(
        '--overlap', $FluoroDBOverlap
    )
    Invoke-SankeyJob 'FluoroDB — hierarchy Sankey (all groups)' $HierarchyScript @(
        '--overlap', $FluoroDBOverlap, '--all'
    )
}

# CLinventory
if (Skip-Missing $CLinventoryOverlap 'CLinventory hierarchy') {
    Invoke-SankeyJob 'CLinventory — hierarchy Sankey (compact)' $HierarchyScript @(
        '--overlap', $CLinventoryOverlap
    )
    Invoke-SankeyJob 'CLinventory — hierarchy Sankey (all groups)' $HierarchyScript @(
        '--overlap', $CLinventoryOverlap, '--all'
    )
}

# OECD
if (Skip-Missing $OecdOverlap 'OECD hierarchy') {
    Invoke-SankeyJob 'OECD — hierarchy Sankey (compact)' $HierarchyScript @(
        '--overlap', $OecdOverlap
    )
    Invoke-SankeyJob 'OECD — hierarchy Sankey (all groups)' $HierarchyScript @(
        '--overlap', $OecdOverlap, '--all'
    )
}

# ── Section 2: Atlas Sankeys (PFASGroups → PFAS-Atlas categories) ─────────────
Write-Host "`n═══ Atlas Sankeys ═══" -ForegroundColor Cyan

Invoke-SankeyJob 'Groups 1-28 — top-28 labels, exclude non-PFAS' $AtlasScript @(
    '--top-n', '28', '--groups', '1-28', '--exclude-not-pfas'
)

Invoke-SankeyJob 'Groups 1-28 — top-12 labels, exclude non-PFAS' $AtlasScript @(
    '--top-n', '12', '--groups', '1-28', '--exclude-not-pfas'
)

Invoke-SankeyJob 'Groups 29-115 — top-12 labels, exclude non-PFAS' $AtlasScript @(
    '--top-n', '12', '--groups', '29-115', '--exclude-not-pfas'
)

Invoke-SankeyJob 'Groups 29-115 — top-25 labels, exclude non-PFAS' $AtlasScript @(
    '--top-n', '25', '--groups', '29-115', '--exclude-not-pfas'
)

Invoke-SankeyJob 'All groups — top-12 labels, exclude non-PFAS' $AtlasScript @(
    '--top-n', '12', '--exclude-not-pfas'
)

# ── Summary ───────────────────────────────────────────────────────────────────
$totalSw.Stop()
$total = $ok + $fail
Write-Host ""
Write-Host "════════════════════════════════════════"
Write-Host "Summary  ($([int]$totalSw.Elapsed.TotalSeconds)s total)"
Write-Host "  $ok/$total jobs succeeded" -ForegroundColor Green
if ($fail -gt 0) {
    Write-Host "  $fail job(s) failed:" -ForegroundColor Red
    foreach ($j in $failedJobs) {
        Write-Host "    ✘  $j" -ForegroundColor Red
    }
    exit 1
}
Write-Host "════════════════════════════════════════"
