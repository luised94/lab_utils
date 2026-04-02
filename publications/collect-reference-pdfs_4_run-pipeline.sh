#!/usr/bin/env bash
# =============================================================================
# COLLECT REFERENCE PDFs - FULL PIPELINE
# =============================================================================
# Orchestrates the three pipeline scripts:
#   1. Extract Zotero item keys from .docx
#   2. Resolve keys to PDF file paths via SQLite
#   3. Copy PDFs to target directory via PowerShell
#
# Usage:
#   bash collect-reference-pdfs_4_run-pipeline.sh <docx_path> <target_directory> [--dry-run]
#
# Example:
#   bash collect-reference-pdfs_4_run-pipeline.sh \
#       "/mnt/c/.../manuscript.docx" \
#       "/mnt/c/.../reference_pdfs"
#
# Prerequisites:
#   - python3, sqlite3, unzip, powershell.exe
#   - Zotero must be closed
#   - All three pipeline scripts in the same directory as this script
# =============================================================================
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

DRY_RUN_FLAG=""

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [ $# -lt 2 ]; then
    echo "Usage: bash collect-reference-pdfs_4_run-pipeline.sh <docx_path> <target_directory> [--dry-run]" >&2
    exit 1
fi

docx_path="$1"
target_directory="$2"

if [ "${3:-}" = "--dry-run" ]; then
    DRY_RUN_FLAG="--dry-run"
fi

if [ ! -f "$docx_path" ]; then
    echo "Error: .docx file not found: ${docx_path}" >&2
    exit 1
fi

# Resolve script directory (all scripts must be co-located)
script_directory="$(cd "$(dirname "$0")" && pwd)"
extract_script="${script_directory}/collect-reference-pdfs_1_extract-keys.py"
resolve_script="${script_directory}/collect-reference-pdfs_2_resolve-paths.sh"
copy_script="${script_directory}/collect-reference-pdfs_3_copy-pdfs.sh"

for script_file in "$extract_script" "$resolve_script" "$copy_script"; do
    if [ ! -f "$script_file" ]; then
        echo "Error: Pipeline script not found: ${script_file}" >&2
        exit 1
    fi
done

# Intermediate files go in a temp directory
working_directory=$(mktemp -d /tmp/zotero_pipeline_XXXXXX)
docx_stem=$(basename "$docx_path" .docx)

echo "==========================================" >&2
echo "COLLECT REFERENCE PDFs - PIPELINE" >&2
echo "==========================================" >&2
echo "Document:   ${docx_path}" >&2
echo "Target:     ${target_directory}" >&2
echo "Working:    ${working_directory}" >&2
echo "Dry run:    ${DRY_RUN_FLAG:-false}" >&2
echo "" >&2

# =============================================================================
# STEP 1: EXTRACT KEYS
# =============================================================================

echo "==========================================" >&2
echo "STEP 1/3: Extract citation keys from .docx" >&2
echo "==========================================" >&2

python3 "$extract_script" "$docx_path" "$working_directory"

keys_file="${working_directory}/${docx_stem}_keys.txt"
if [ ! -f "$keys_file" ]; then
    echo "Error: Keys file not created: ${keys_file}" >&2
    exit 1
fi

key_count=$(wc -l < "$keys_file")
echo "" >&2

# =============================================================================
# STEP 2: RESOLVE PATHS
# =============================================================================

echo "==========================================" >&2
echo "STEP 2/3: Resolve keys to PDF paths" >&2
echo "==========================================" >&2

bash "$resolve_script" "$keys_file" "$working_directory"

paths_file="${working_directory}/${docx_stem}_keys_paths.tsv"
unresolved_file="${working_directory}/${docx_stem}_keys_unresolved.tsv"

if [ ! -f "$paths_file" ]; then
    echo "Error: Paths file not created: ${paths_file}" >&2
    exit 1
fi

resolved_count=$(( $(wc -l < "$paths_file") - 1 ))
unresolved_count=0
if [ -f "$unresolved_file" ]; then
    unresolved_count=$(( $(wc -l < "$unresolved_file") - 1 ))
fi

echo "" >&2

# =============================================================================
# STEP 3: COPY PDFs
# =============================================================================

echo "==========================================" >&2
echo "STEP 3/3: Copy PDFs to target directory" >&2
echo "==========================================" >&2

bash "$copy_script" "$paths_file" "$target_directory" $DRY_RUN_FLAG

echo "" >&2

# =============================================================================
# PIPELINE SUMMARY
# =============================================================================

echo "==========================================" >&2
echo "PIPELINE COMPLETE" >&2
echo "==========================================" >&2
echo "Citation keys extracted:  ${key_count}" >&2
echo "Resolved to PDF paths:   ${resolved_count}" >&2
echo "Unresolved (no PDF):     ${unresolved_count}" >&2
echo "Target directory:        ${target_directory}" >&2
echo "" >&2

if [ -f "$unresolved_file" ] && [ "$unresolved_count" -gt 0 ]; then
    echo "Unresolved items:" >&2
    tail -n +2 "$unresolved_file" | while IFS=$'\t' read -r key item_type title; do
        echo "  ${key}  ${item_type}  ${title}" >&2
    done
    echo "" >&2
fi

echo "Working files:  ${working_directory}" >&2
echo "Cleanup:        rm -rf ${working_directory}" >&2
