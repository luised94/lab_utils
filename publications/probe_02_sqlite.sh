#!/usr/bin/env bash
# =============================================================================
# PROBE 02: ZOTERO SQLITE SCHEMA AND LINKED FILE PATHS
# =============================================================================
# Destructive operations: NONE.
#   - Copies zotero.sqlite to /tmp, queries the copy. Original untouched.
#   - All output to stdout.
#
# Prerequisite: Zotero must be closed (probe_00 checks this).
#
# Assumptions tested:
#   A5. items.key matches the key from Word URIs (probe_01 output).
#   A6. itemAttachments.path holds the linked file path for Attanger files.
#   A7. linkMode=2 indicates linked files.
#
# Expected outcome: We see the schema, path format (absolute vs relative),
#   and can verify a key-to-path lookup works end to end.
# =============================================================================
set -euo pipefail

ZOTERO_DB="/mnt/c/Users/Luised94/Zotero/zotero.sqlite"
PREFS_FILE="/mnt/c/Users/Luised94/Zotero/prefs.js"

# Accept probe_01 temp dir as argument, or find most recent one
if [ -n "${1:-}" ] && [ -d "$1" ]; then
    PROBE_DIR="$1"
else
    PROBE_DIR=$(ls -dt /tmp/zotero_probe_docx_* 2>/dev/null | head -1 || true)
fi

KEYS_FILE="${PROBE_DIR}/zotero_item_keys.txt"
DB_COPY=$(mktemp /tmp/zotero_probe_db_XXXXXX.sqlite)

echo "========================================"
echo "PROBE 02: ZOTERO SQLITE"
echo "========================================"
echo ""

# --- Check Zotero closed ---
ZOTERO_PID=$(powershell.exe -NoProfile -Command \
    "Get-Process zotero -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Id" \
    2>/dev/null | tr -d '\r' || true)

if [ -n "$ZOTERO_PID" ]; then
    echo "FAIL: Zotero is running (PID: ${ZOTERO_PID}). Close it first."
    exit 1
fi

# --- Copy DB ---
echo "--- 2a. Copying database ---"
cp "$ZOTERO_DB" "$DB_COPY"
[ -f "${ZOTERO_DB}-wal" ] && cp "${ZOTERO_DB}-wal" "${DB_COPY}-wal"
[ -f "${ZOTERO_DB}-shm" ] && cp "${ZOTERO_DB}-shm" "${DB_COPY}-shm"
echo "  Copy: ${DB_COPY}"
echo ""

# --- Schema ---
echo "--- 2b. Relevant table schemas ---"
echo ""
echo ">> items:"
sqlite3 "$DB_COPY" ".schema items"
echo ""
echo ">> itemAttachments:"
sqlite3 "$DB_COPY" ".schema itemAttachments"
echo ""

# --- linkMode distribution ---
echo "--- 2c. Attachment linkMode distribution ---"
echo "  (0=imported_file, 1=imported_url, 2=linked_file, 3=linked_url)"
sqlite3 -header -column "$DB_COPY" \
    "SELECT linkMode, contentType, COUNT(*) as count
     FROM itemAttachments
     GROUP BY linkMode, contentType
     ORDER BY count DESC
     LIMIT 15;"
echo ""

# --- Sample linked PDFs ---
echo "--- 2d. Sample linked PDF attachments (linkMode=2) ---"
sqlite3 -header -column "$DB_COPY" \
    "SELECT
         pi.key AS parent_key,
         ci.key AS attach_key,
         ia.path,
         ia.linkMode
     FROM itemAttachments ia
     JOIN items ci ON ia.itemID = ci.itemID
     JOIN items pi ON ia.parentItemID = pi.itemID
     WHERE ia.linkMode = 2
       AND ia.contentType = 'application/pdf'
     LIMIT 5;"
echo ""

# --- Path format ---
echo "--- 2e. Path format analysis ---"
sqlite3 -header -column "$DB_COPY" \
    "SELECT
         CASE
             WHEN path LIKE 'attachments:%' THEN 'attachments: prefix'
             WHEN path LIKE 'storage:%' THEN 'storage: prefix'
             WHEN path LIKE '/%' THEN 'absolute unix'
             WHEN path LIKE '_:\\%' THEN 'absolute windows'
             ELSE 'other: ' || SUBSTR(path, 1, 40)
         END AS format,
         COUNT(*) as count
     FROM itemAttachments
     WHERE linkMode = 2 AND contentType = 'application/pdf'
     GROUP BY format;"
echo ""

# --- Zotero prefs (baseAttachmentPath) ---
echo "--- 2f. Attachment path preferences ---"
if [ -f "$PREFS_FILE" ]; then
    grep -iE "baseAttachmentPath|attanger|linkedAttachment" "$PREFS_FILE" 2>/dev/null || echo "  No matching prefs found"
else
    echo "  prefs.js not found"
fi
echo ""

# --- Key-to-path test using probe_01 keys ---
echo "--- 2g. Key-to-path lookup test ---"
if [ -f "$KEYS_FILE" ]; then
    KEY_COUNT=$(wc -l < "$KEYS_FILE")
    echo "  Found ${KEY_COUNT} keys from probe_01: ${KEYS_FILE}"
    SAMPLE_KEYS=$(head -3 "$KEYS_FILE" | tr '\n' ',' | sed 's/,$//')
    echo "  Testing with first 3: ${SAMPLE_KEYS}"
    echo ""

    # Build SQL IN clause
    IN_CLAUSE=$(head -3 "$KEYS_FILE" | sed "s/.*/'&'/" | paste -sd,)

    sqlite3 -header -column "$DB_COPY" \
        "SELECT
             pi.key AS item_key,
             ia.path,
             ia.linkMode
         FROM items pi
         JOIN itemAttachments ia ON ia.parentItemID = pi.itemID
         WHERE pi.key IN (${IN_CLAUSE})
           AND ia.linkMode = 2
           AND ia.contentType = 'application/pdf';"

    # Count how many keys resolve
    RESOLVED=$(sqlite3 "$DB_COPY" \
        "SELECT COUNT(DISTINCT pi.key)
         FROM items pi
         JOIN itemAttachments ia ON ia.parentItemID = pi.itemID
         WHERE pi.key IN ($(sed "s/.*/'&'/" "$KEYS_FILE" | paste -sd,))
           AND ia.linkMode = 2
           AND ia.contentType = 'application/pdf';")

    echo ""
    echo "  Keys from docx: ${KEY_COUNT}"
    echo "  Keys resolved to PDFs: ${RESOLVED}"
    if [ "$RESOLVED" -lt "$KEY_COUNT" ]; then
        echo "  WARNING: $(( KEY_COUNT - RESOLVED )) keys did not resolve."
        echo "  Possible reasons: no PDF attached, different linkMode, or key mismatch."
    fi
else
    echo "  No keys file found. Run probe_01 first."
    echo "  Expected: ${KEYS_FILE}"
fi
echo ""

# --- Cleanup note ---
echo "========================================"
echo "PROBE 02 COMPLETE"
echo "========================================"
echo "Cleanup: rm -f ${DB_COPY} ${DB_COPY}-wal ${DB_COPY}-shm"
echo "Probe 01 temp: rm -rf ${PROBE_DIR}"
