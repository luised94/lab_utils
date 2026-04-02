#!/usr/bin/env bash
# =============================================================================
# RESOLVE ZOTERO ITEM KEYS TO PDF FILE PATHS
# =============================================================================
# Takes a keys file (one key per line, from extract-keys.py) and queries
# zotero.sqlite to find the absolute path to each linked PDF attachment.
#
# Usage:
#   bash collect-reference-pdfs_2_resolve-paths.sh <keys_file> [output_directory]
#
# Output:
#   <output_directory>/<keys_stem>_paths.tsv       (key \t relative_path \t absolute_path)
#   <output_directory>/<keys_stem>_unresolved.tsv   (key \t item_type \t title)
#
# Prerequisites:
#   - Zotero must be closed (script checks and exits if running)
#   - sqlite3 must be installed
#
# Path resolution:
#   Attanger/ZotFile store linked files with paths like:
#       attachments:AuthorLastName/Author_Year_Title.pdf
#   This script strips the "attachments:" prefix and prepends the
#   baseAttachmentPath from Zotero's prefs.js.
# =============================================================================
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

ZOTERO_SQLITE_PATH="/mnt/c/Users/Luised94/Zotero/zotero.sqlite"
ZOTERO_PREFS_GLOB="/mnt/c/Users/Luised94/AppData/Roaming/Zotero/Zotero/Profiles/*/prefs.js"
ATTACHMENTS_PREFIX="attachments:"

# Fallback if prefs.js cannot be read or does not contain baseAttachmentPath
FALLBACK_BASE_ATTACHMENT_PATH="/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/zotero-storage"

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [ $# -lt 1 ]; then
    echo "Usage: bash collect-reference-pdfs_2_resolve-paths.sh <keys_file> [output_directory]" >&2
    exit 1
fi

keys_file="$1"
output_directory="${2:-.}"

if [ ! -f "$keys_file" ]; then
    echo "Error: Keys file not found: ${keys_file}" >&2
    exit 1
fi

if [ ! -d "$output_directory" ]; then
    echo "Error: Output directory not found: ${output_directory}" >&2
    exit 1
fi

if ! command -v sqlite3 &>/dev/null; then
    echo "Error: sqlite3 not installed (sudo apt install sqlite3)" >&2
    exit 1
fi

key_count=$(wc -l < "$keys_file")
if [ "$key_count" -eq 0 ]; then
    echo "Error: Keys file is empty: ${keys_file}" >&2
    exit 1
fi

keys_stem=$(basename "$keys_file" .txt)
output_paths_file="${output_directory}/${keys_stem}_paths.tsv"
output_unresolved_file="${output_directory}/${keys_stem}_unresolved.tsv"

echo "Keys file:   ${keys_file} (${key_count} keys)" >&2
echo "Output dir:  ${output_directory}" >&2
echo "" >&2

# =============================================================================
# CHECK ZOTERO IS NOT RUNNING
# =============================================================================

zotero_pid=$(powershell.exe -NoProfile -Command \
    "Get-Process zotero -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Id" \
    2>/dev/null | tr -d '\r' || true)

if [ -n "$zotero_pid" ]; then
    echo "Error: Zotero is running (PID: ${zotero_pid}). Close it first." >&2
    exit 1
fi

echo "Zotero process check: not running" >&2

# =============================================================================
# CHECK SQLITE FILE EXISTS
# =============================================================================

if [ ! -f "$ZOTERO_SQLITE_PATH" ]; then
    echo "Error: zotero.sqlite not found at: ${ZOTERO_SQLITE_PATH}" >&2
    exit 1
fi

# =============================================================================
# RESOLVE baseAttachmentPath FROM prefs.js
# =============================================================================

base_attachment_path=""

# shellcheck disable=SC2086
prefs_file=$(ls $ZOTERO_PREFS_GLOB 2>/dev/null | head -1 || true)

if [ -n "$prefs_file" ] && [ -f "$prefs_file" ]; then
    echo "Found prefs.js: ${prefs_file}" >&2

    # prefs.js stores the path as:
    #   user_pref("extensions.zotero.baseAttachmentPath", "C:\\Users\\...");
    # Line format: user_pref("extensions.zotero.baseAttachmentPath", "C:\\Users\\...");
    # Extract the second quoted string (the value) from the matching line
    raw_pref_value=$(grep 'extensions\.zotero\.baseAttachmentPath' "$prefs_file" 2>/dev/null \
        | sed 's/.*", "//;s/".*//' || true)

    if [ -n "$raw_pref_value" ]; then
        # Convert Windows path (with escaped backslashes) to WSL path
        # prefs.js stores: C:\\Users\\... - first unescape, then convert
        windows_path=$(echo "$raw_pref_value" | sed 's|\\\\|\\|g')
        # Convert C:\path\to\dir to /mnt/c/path/to/dir
        drive_letter=$(echo "$windows_path" | head -c 1 | tr '[:upper:]' '[:lower:]')
        remaining_path=$(echo "$windows_path" | cut -c3- | sed 's|\\|/|g')
        base_attachment_path="/mnt/${drive_letter}${remaining_path}"
        echo "baseAttachmentPath from prefs: ${base_attachment_path}" >&2
    else
        echo "baseAttachmentPath not found in prefs.js, using fallback" >&2
    fi
else
    echo "prefs.js not found at glob pattern, using fallback" >&2
fi

if [ -z "$base_attachment_path" ]; then
    base_attachment_path="$FALLBACK_BASE_ATTACHMENT_PATH"
    echo "Using fallback path: ${base_attachment_path}" >&2
fi

if [ ! -d "$base_attachment_path" ]; then
    echo "Error: Base attachment path does not exist: ${base_attachment_path}" >&2
    exit 1
fi

echo "Resolved base path: ${base_attachment_path}" >&2
echo "" >&2

# =============================================================================
# COPY SQLITE TO TEMP (AVOID WAL/LOCK ISSUES)
# =============================================================================

temporary_db=$(mktemp /tmp/zotero_resolve_XXXXXX.sqlite)
cp "$ZOTERO_SQLITE_PATH" "$temporary_db"
[ -f "${ZOTERO_SQLITE_PATH}-wal" ] && cp "${ZOTERO_SQLITE_PATH}-wal" "${temporary_db}-wal"
[ -f "${ZOTERO_SQLITE_PATH}-shm" ] && cp "${ZOTERO_SQLITE_PATH}-shm" "${temporary_db}-shm"

echo "Copied database to: ${temporary_db}" >&2

# =============================================================================
# BUILD SQL IN CLAUSE FROM KEYS FILE
# =============================================================================

sql_in_clause=$(sed "s/.*/'&'/" "$keys_file" | paste -sd,)

# =============================================================================
# QUERY: RESOLVED KEYS (have linkMode=2 PDF attachments)
# =============================================================================

echo "Querying for resolved keys..." >&2

sqlite3 -separator $'\t' "$temporary_db" \
"SELECT
    pi.key,
    ia.path
FROM items pi
JOIN itemAttachments ia ON ia.parentItemID = pi.itemID
WHERE pi.key IN (${sql_in_clause})
  AND ia.linkMode = 2
  AND ia.contentType = 'application/pdf'
ORDER BY pi.key;" > "${output_directory}/.resolved_raw.tmp"

# =============================================================================
# RESOLVE RELATIVE PATHS TO ABSOLUTE PATHS
# =============================================================================

resolved_count=0
missing_on_disk_count=0

# Write header
printf "key\trelative_path\tabsolute_path\texists_on_disk\n" > "$output_paths_file"

while IFS=$'\t' read -r item_key raw_path; do
    # Strip "attachments:" prefix if present
    if [[ "$raw_path" == "${ATTACHMENTS_PREFIX}"* ]]; then
        relative_path="${raw_path#"${ATTACHMENTS_PREFIX}"}"
        absolute_path="${base_attachment_path}/${relative_path}"
    else
        # Absolute path or other format - use as-is after converting to WSL
        relative_path="$raw_path"
        # Try Windows-to-WSL conversion if it looks like a Windows path
        if [[ "$raw_path" =~ ^[A-Za-z]:\\ ]]; then
            drive_letter=$(echo "$raw_path" | head -c 1 | tr '[:upper:]' '[:lower:]')
            remaining=$(echo "$raw_path" | cut -c3- | sed 's|\\|/|g')
            absolute_path="/mnt/${drive_letter}${remaining}"
        else
            absolute_path="$raw_path"
        fi
    fi

    exists_on_disk="yes"
    if [ ! -f "$absolute_path" ]; then
        exists_on_disk="no"
        missing_on_disk_count=$((missing_on_disk_count + 1))
    fi

    printf "%s\t%s\t%s\t%s\n" "$item_key" "$relative_path" "$absolute_path" "$exists_on_disk" >> "$output_paths_file"
    resolved_count=$((resolved_count + 1))

done < "${output_directory}/.resolved_raw.tmp"

rm -f "${output_directory}/.resolved_raw.tmp"

echo "Resolved ${resolved_count} keys to PDF paths" >&2

# =============================================================================
# QUERY: UNRESOLVED KEYS (no linkMode=2 PDF attachment)
# =============================================================================

echo "Querying for unresolved keys..." >&2

# Write header
printf "key\titem_type\ttitle\n" > "$output_unresolved_file"

sqlite3 -separator $'\t' "$temporary_db" \
"SELECT
    pi.key,
    COALESCE(it.typeName, 'unknown'),
    COALESCE(
        (SELECT iv.value
         FROM itemData id2
         JOIN itemDataValues iv ON id2.valueID = iv.valueID
         JOIN fields f ON id2.fieldID = f.fieldID
         WHERE id2.itemID = pi.itemID AND f.fieldName = 'title'),
        'no title')
FROM items pi
LEFT JOIN itemTypes it ON pi.itemTypeID = it.itemTypeID
WHERE pi.key IN (${sql_in_clause})
  AND pi.key NOT IN (
      SELECT pi2.key
      FROM items pi2
      JOIN itemAttachments ia ON ia.parentItemID = pi2.itemID
      WHERE pi2.key IN (${sql_in_clause})
        AND ia.linkMode = 2
        AND ia.contentType = 'application/pdf'
  )
ORDER BY pi.key;" >> "$output_unresolved_file"

# Count unresolved (subtract 1 for header)
unresolved_count=$(( $(wc -l < "$output_unresolved_file") - 1 ))

# =============================================================================
# CLEANUP TEMP DB
# =============================================================================

rm -f "$temporary_db" "${temporary_db}-wal" "${temporary_db}-shm"

# =============================================================================
# SUMMARY
# =============================================================================

echo "" >&2
echo "--- Summary ---" >&2
echo "Input keys:          ${key_count}" >&2
echo "Resolved to PDF:     ${resolved_count}" >&2
echo "Unresolved:          ${unresolved_count}" >&2
echo "Missing on disk:     ${missing_on_disk_count}" >&2
echo "" >&2
echo "Paths file:          ${output_paths_file}" >&2
echo "Unresolved file:     ${output_unresolved_file}" >&2

if [ "$missing_on_disk_count" -gt 0 ]; then
    echo "" >&2
    echo "WARNING: ${missing_on_disk_count} resolved paths do not exist on disk." >&2
    echo "  These may be Dropbox cloud-only files (will hydrate during copy step)." >&2
    echo "  Or the files may have been moved/deleted." >&2
fi

if [ "$unresolved_count" -gt 0 ]; then
    echo "" >&2
    echo "Unresolved items (no linked PDF):" >&2
    tail -n +2 "$output_unresolved_file" | while IFS=$'\t' read -r key item_type title; do
        echo "  ${key}  ${item_type}  ${title}" >&2
    done
fi
