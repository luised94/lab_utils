#!/usr/bin/env bash
# =============================================================================
# COPY RESOLVED PDFs TO TARGET DIRECTORY
# =============================================================================
# Reads the paths TSV from resolve-paths.sh, deduplicates by file path,
# generates a temporary PowerShell script, and executes it once to:
#   1. Create the target directory
#   2. Hydrate each Dropbox cloud-only file
#   3. Copy each PDF to the target directory
#   4. Validate each copy exists at destination
#
# WHY POWERSHELL INSTEAD OF CP:
#   Dropbox smart-sync stores files as cloud-only placeholders on disk.
#   These appear as normal files to `ls` but contain no actual data.
#   PowerShell file operations (Copy-Item, Get-Content) trigger Dropbox
#   to download ("hydrate") the real content automatically. WSL's `cp`
#   does not trigger this hydration and would copy empty stubs.
#
# WHY A SINGLE .PS1 SCRIPT:
#   Each `powershell.exe` invocation from WSL takes ~1-2 seconds of startup
#   overhead. For 100 files, that's 2-3 minutes of pure overhead. Generating
#   a .ps1 and executing it once eliminates this.
#
# WHY -LiteralPath:
#   PowerShell's default -Path parameter interprets wildcards and special
#   characters (brackets, backticks, dollar signs). -LiteralPath treats the
#   path as a literal string. Academic filenames frequently contain Unicode
#   (”, ú, -, ?) and special characters that would break -Path.
#
# HYDRATION CHECK (0x00400000):
#   The Windows file attribute RECALL_ON_DATA_ACCESS (0x00400000) indicates
#   a Dropbox cloud-only placeholder. The script checks this attribute before
#   copying, hydrates if needed, and validates the attribute is cleared.
#   Pattern taken from Backup-ZoteroStorage.ps1.
#
# Uses -LiteralPath throughout and escapes filenames for PowerShell safety.
#
# Usage:
#   bash collect-reference-pdfs_3_copy-pdfs.sh <paths_tsv> <target_directory> [--dry-run]
#
# Output:
#   PDFs copied to <target_directory>/
#   Results TSV written to <target_directory>/_copy_results.tsv
#   Summary printed to stderr
#
# Prerequisites:
#   - powershell.exe available (WSL interop)
#   - paths TSV from collect-reference-pdfs_2_resolve-paths.sh
# =============================================================================
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

DRY_RUN=false

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [ $# -lt 2 ]; then
    echo "Usage: bash collect-reference-pdfs_3_copy-pdfs.sh <paths_tsv> <target_directory> [--dry-run]" >&2
    exit 1
fi

paths_tsv_file="$1"
target_directory="$2"

if [ "${3:-}" = "--dry-run" ]; then
    DRY_RUN=true
fi

if [ ! -f "$paths_tsv_file" ]; then
    echo "Error: Paths file not found: ${paths_tsv_file}" >&2
    exit 1
fi

if ! command -v powershell.exe &>/dev/null; then
    echo "Error: powershell.exe not available (WSL interop required)" >&2
    exit 1
fi

# Count data rows (subtract 1 for header)
total_rows=$(( $(wc -l < "$paths_tsv_file") - 1 ))
if [ "$total_rows" -le 0 ]; then
    echo "Error: No data rows in paths file" >&2
    exit 1
fi

echo "Paths file:       ${paths_tsv_file} (${total_rows} rows)" >&2
echo "Target directory: ${target_directory}" >&2
echo "Dry run:          ${DRY_RUN}" >&2
echo "" >&2

# =============================================================================
# CONVERT TARGET DIRECTORY TO WINDOWS PATH
# =============================================================================

if [[ "$target_directory" != /mnt/[a-z]/* ]]; then
    echo "Error: Target directory must be under /mnt/<drive>/... for PowerShell access" >&2
    echo "  Got: ${target_directory}" >&2
    exit 1
fi

drive_letter=$(echo "$target_directory" | cut -c6 | tr '[:lower:]' '[:upper:]')
remaining_path=$(echo "$target_directory" | cut -c7-)
target_directory_windows="${drive_letter}:${remaining_path//\//\\}"

echo "Target (Windows): ${target_directory_windows}" >&2
echo "" >&2

# =============================================================================
# DEDUPLICATE PATHS
# =============================================================================

deduplicated_paths_file=$(mktemp /tmp/zotero_dedup_XXXXXX.tsv)

# Skip header, extract absolute_path (column 3), sort unique
tail -n +2 "$paths_tsv_file" | cut -f3 | sort -u > "$deduplicated_paths_file"

unique_file_count=$(wc -l < "$deduplicated_paths_file")
duplicate_count=$(( total_rows - unique_file_count ))

echo "Unique files to copy: ${unique_file_count}" >&2
if [ "$duplicate_count" -gt 0 ]; then
    echo "Duplicates skipped:   ${duplicate_count}" >&2
fi
echo "" >&2

# =============================================================================
# CONVERT WSL PATHS TO WINDOWS PATHS AND ESCAPE FOR POWERSHELL
# =============================================================================

windows_paths_file=$(mktemp /tmp/zotero_winpaths_XXXXXX.txt)

while IFS= read -r absolute_path; do
    if [[ "$absolute_path" == /mnt/[a-z]/* ]]; then
        source_drive=$(echo "$absolute_path" | cut -c6 | tr '[:lower:]' '[:upper:]')
        source_remaining=$(echo "$absolute_path" | cut -c7-)
        windows_path="${source_drive}:${source_remaining//\//\\}"
        echo "$windows_path" >> "$windows_paths_file"
    else
        echo "WARNING: Skipping non-/mnt/ path: ${absolute_path}" >&2
    fi
done < "$deduplicated_paths_file"

windows_path_count=$(wc -l < "$windows_paths_file")
echo "Converted ${windows_path_count} paths to Windows format" >&2
echo "" >&2

# =============================================================================
# GENERATE POWERSHELL SCRIPT
# =============================================================================

powershell_script_file=$(mktemp /tmp/zotero_copy_XXXXXX.ps1)
results_filename="_copy_results.tsv"

# Escape a string for use inside PowerShell single quotes: ' becomes ''
escape_for_powershell_single_quote() {
    echo "$1" | sed "s/'/''/g"
}

cat > "$powershell_script_file" << 'PS_HEADER'
param(
    [string]$TargetDir,
    [string]$SourcePathsFile,
    [string]$ResultsFile,
    [switch]$DryRun
)
$ErrorActionPreference = "Continue"

# --- Create target directory ---
if (-not $DryRun) {
    if (-not (Test-Path -LiteralPath $TargetDir)) {
        New-Item -Path $TargetDir -ItemType Directory -Force | Out-Null
        Write-Host "[INFO]  Created target directory: $TargetDir"
    } else {
        Write-Host "[INFO]  Target directory exists: $TargetDir"
    }
}

# --- Read source paths ---
$SourcePaths = Get-Content -LiteralPath $SourcePathsFile -Encoding UTF8
$TotalFiles = $SourcePaths.Count
Write-Host "[INFO]  Processing $TotalFiles files..."

# --- Results tracking ---
$Results = [System.Collections.Generic.List[PSCustomObject]]::new()
$CopiedCount = 0
$SkippedCount = 0
$HydrateFailCount = 0
$CopyFailCount = 0
$ValidateFailCount = 0
$Current = 0

foreach ($SourcePath in $SourcePaths) {
    $Current++
    $FileName = Split-Path $SourcePath -Leaf
    $DestPath = Join-Path $TargetDir $FileName

    # --- Skip if already exists at destination ---
    if (Test-Path -LiteralPath $DestPath) {
        Write-Host "  [$Current/$TotalFiles] SKIP (exists): $FileName"
        $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="SKIPPED"; Error="" })
        $SkippedCount++
        continue
    }

    if ($DryRun) {
        Write-Host "  [$Current/$TotalFiles] DRY RUN: $FileName"
        $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="DRY_RUN"; Error="" })
        $CopiedCount++
        continue
    }

    # --- Check source exists ---
    if (-not (Test-Path -LiteralPath $SourcePath)) {
        Write-Host "  [$Current/$TotalFiles] FAIL (not found): $FileName" -ForegroundColor Yellow
        $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="NOT_FOUND"; Error="Source file does not exist" })
        $CopyFailCount++
        continue
    }

    # --- Hydrate if Dropbox cloud-only placeholder ---
    try {
        $Attrs = [System.IO.File]::GetAttributes($SourcePath)
        $IsPlaceholder = ($Attrs -band 0x00400000) -ne 0
    } catch {
        $IsPlaceholder = $false
    }

    if ($IsPlaceholder) {
        try {
            Get-Content -LiteralPath $SourcePath -TotalCount 1 -ErrorAction Stop | Out-Null
        } catch {
            Write-Host "  [$Current/$TotalFiles] FAIL (hydrate): $FileName" -ForegroundColor Yellow
            $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="HYDRATE_FAIL"; Error=$_.Exception.Message })
            $HydrateFailCount++
            continue
        }

        # Validate hydration completed
        try {
            $AttrsPost = [System.IO.File]::GetAttributes($SourcePath)
            if ($AttrsPost -band 0x00400000) {
                Write-Host "  [$Current/$TotalFiles] FAIL (still placeholder): $FileName" -ForegroundColor Yellow
                $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="HYDRATE_FAIL"; Error="Still placeholder after hydration" })
                $HydrateFailCount++
                continue
            }
        } catch {
            # Attribute check failed but hydration may have worked - proceed with copy
        }
    }

    # --- Copy ---
    try {
        Copy-Item -LiteralPath $SourcePath -Destination $DestPath -Force -ErrorAction Stop
    } catch {
        Write-Host "  [$Current/$TotalFiles] FAIL (copy): $FileName" -ForegroundColor Yellow
        $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="COPY_FAIL"; Error=$_.Exception.Message })
        $CopyFailCount++
        continue
    }

    # --- Validate copy exists at destination ---
    if (-not (Test-Path -LiteralPath $DestPath)) {
        Write-Host "  [$Current/$TotalFiles] FAIL (validate): $FileName" -ForegroundColor Yellow
        $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="VALIDATE_FAIL"; Error="File not found at destination after copy" })
        $ValidateFailCount++
        continue
    }

    Write-Host "  [$Current/$TotalFiles] OK: $FileName"
    $Results.Add([PSCustomObject]@{ Path=$SourcePath; Status="OK"; Error="" })
    $CopiedCount++
}

# --- Write results TSV ---
$ResultLines = [System.Collections.Generic.List[string]]::new()
$ResultLines.Add("source_path`tstatus`terror")
foreach ($R in $Results) {
    $ResultLines.Add("$($R.Path)`t$($R.Status)`t$($R.Error)")
}
[System.IO.File]::WriteAllLines($ResultsFile, $ResultLines, [System.Text.Encoding]::UTF8)

# --- Summary ---
Write-Host ""
Write-Host "--- Summary ---"
Write-Host "Total files:       $TotalFiles"
Write-Host "Copied:            $CopiedCount"
Write-Host "Already existed:   $SkippedCount"
Write-Host "Hydration failed:  $HydrateFailCount"
Write-Host "Copy failed:       $CopyFailCount"
Write-Host "Validation failed: $ValidateFailCount"

if ($DryRun) {
    Write-Host ""
    Write-Host "(DRY RUN - no files were actually copied)"
}

# Exit with error if any failures
$FailTotal = $HydrateFailCount + $CopyFailCount + $ValidateFailCount
if ($FailTotal -gt 0) { exit 1 }
exit 0
PS_HEADER

echo "Generated PowerShell script: ${powershell_script_file}" >&2

# =============================================================================
# CONVERT PATHS FILE TO WINDOWS PATH FOR POWERSHELL ACCESS
# =============================================================================

# PowerShell needs Windows paths for -LiteralPath on the paths file and results file
windows_paths_file_win=$(wslpath -w "$windows_paths_file")
results_file_wsl="${target_directory}/${results_filename}"

if [[ "$results_file_wsl" == /mnt/[a-z]/* ]]; then
    results_drive=$(echo "$results_file_wsl" | cut -c6 | tr '[:lower:]' '[:upper:]')
    results_remaining=$(echo "$results_file_wsl" | cut -c7-)
    results_file_windows="${results_drive}:${results_remaining//\//\\}"
else
    results_file_windows="$results_file_wsl"
fi

# =============================================================================
# EXECUTE POWERSHELL SCRIPT
# =============================================================================

echo "Executing PowerShell copy script..." >&2
echo "" >&2

powershell_script_file_win=$(wslpath -w "$powershell_script_file")

dry_run_flag=""
if [ "$DRY_RUN" = true ]; then
    dry_run_flag="-DryRun"
fi

# Create target directory from bash side first (needed for results file on dry run)
mkdir -p "$target_directory"

powershell.exe -NoProfile -ExecutionPolicy Bypass \
    -File "$powershell_script_file_win" \
    -TargetDir "$target_directory_windows" \
    -SourcePathsFile "$windows_paths_file_win" \
    -ResultsFile "$results_file_windows" \
    $dry_run_flag

powershell_exit_code=$?

echo "" >&2

# =============================================================================
# PARSE RESULTS
# =============================================================================

if [ -f "$results_file_wsl" ]; then
    failed_lines=$(tail -n +2 "$results_file_wsl" | grep -v -E "OK|SKIPPED|DRY_RUN" || true)
    if [ -n "$failed_lines" ]; then
        echo "Failed files:" >&2
        echo "$failed_lines" | while IFS=$'\t' read -r path status error; do
            echo "  [${status}] $(basename "$path"): ${error}" >&2
        done
    fi
    echo "" >&2
    echo "Results file: ${results_file_wsl}" >&2
else
    echo "WARNING: Results file not created at ${results_file_wsl}" >&2
fi

# =============================================================================
# CLEANUP
# =============================================================================

rm -f "$deduplicated_paths_file"
rm -f "$windows_paths_file"
rm -f "$powershell_script_file"

# =============================================================================
# EXIT
# =============================================================================

if [ "$powershell_exit_code" -ne 0 ]; then
    echo "PowerShell exited with code ${powershell_exit_code} (some files failed)" >&2
    exit 1
fi

exit 0
