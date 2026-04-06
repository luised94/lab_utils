#!/usr/bin/env bash
# =============================================================================
# PROBE 00: PREFLIGHT - dependencies, files, Zotero process, paths
# =============================================================================
# Destructive operations: NONE. All read-only checks.
# Writes: Nothing.
# =============================================================================
set -euo pipefail

DOCX_PATH="/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/Lab/publications-and-presentations/lemr_publication_bypass_orc4r/manuscript_orc4r_bypass_v16.docx"
ZOTERO_DB="/mnt/c/Users/Luised94/Zotero/zotero.sqlite"
STORAGE_PATH="/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/zotero-storage"
TARGET_PATH="/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/Lab/publications-and-presentations/lemr_publication_bypass_orc4r/reference_pdfs"

echo "========================================"
echo "PROBE 00: PREFLIGHT"
echo "========================================"
echo ""

# --- Dependencies ---
echo "--- Dependencies ---"
MISSING=()
for cmd in unzip sqlite3 python3 grep powershell.exe; do
    if command -v "$cmd" &>/dev/null; then
        echo "  PASS: ${cmd}"
    else
        echo "  FAIL: ${cmd}"
        MISSING+=("$cmd")
    fi
done

if [ ${#MISSING[@]} -gt 0 ]; then
    echo ""
    echo "Install missing tools (likely: sudo apt install ${MISSING[*]})"
    exit 1
fi
echo ""

# --- Files ---
echo "--- File existence ---"
for label_path in "docx:${DOCX_PATH}" "sqlite:${ZOTERO_DB}"; do
    label="${label_path%%:*}"
    fpath="${label_path#*:}"
    if [ -f "$fpath" ]; then
        echo "  PASS: ${label} ($(stat --printf='%s' "$fpath") bytes)"
    else
        echo "  FAIL: ${label} not found at ${fpath}"
    fi
done
echo ""

# --- Zotero process ---
echo "--- Zotero process ---"
ZOTERO_PID=$(powershell.exe -NoProfile -Command \
    "Get-Process zotero -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Id" \
    2>/dev/null | tr -d '\r' || true)

if [ -n "$ZOTERO_PID" ]; then
    echo "  FAIL: Zotero running (PID: ${ZOTERO_PID}). Close it before probe_02."
else
    echo "  PASS: Zotero not running"
fi
echo ""

# --- Paths ---
echo "--- Dropbox paths ---"
if [ -d "$STORAGE_PATH" ]; then
    PDF_COUNT=$(find "$STORAGE_PATH" -name '*.pdf' -maxdepth 2 2>/dev/null | wc -l)
    echo "  PASS: zotero-storage exists (${PDF_COUNT} PDFs at depth<=2)"
    echo "  Sample files:"
    find "$STORAGE_PATH" -name '*.pdf' -maxdepth 2 2>/dev/null | head -3 | while read -r f; do
        echo "    $(basename "$f")"
    done
else
    echo "  FAIL: zotero-storage not found at ${STORAGE_PATH}"
fi

if [ -d "$TARGET_PATH" ]; then
    echo "  Target dir exists: ${TARGET_PATH}"
else
    echo "  Target dir does not exist yet (will create during pipeline)"
fi
echo ""
echo "Preflight complete. Proceed to probe_01_docx.sh"
