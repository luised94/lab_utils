#!/usr/bin/env bash
# ============================================================
# lw_backup.sh -- additive sync of lab directory to Dropbox
# ============================================================
# Requires LW_LAB_ROOT, LW_DB, and LW_BACKUP_DIR from lw.bash.
# No fallbacks: refuses to run if variables are missing or paths
# are wrong. Additive only -- never deletes from destination.
#
# Usage: bash lw_backup.sh [--dry-run]
# ============================================================
set -euo pipefail

# --- validate environment ---
_bail() {
    echo "lw_backup: error: $1" >&2
    exit 1
}

[ -n "${LW_LAB_ROOT:-}" ]  || _bail "LW_LAB_ROOT not set. Source lw.bash first."
[ -n "${LW_DB:-}"        ]  || _bail "LW_DB not set. Source lw.bash first."
[ -n "${LW_BACKUP_DIR:-}" ] || _bail "LW_BACKUP_DIR not set. Source lw.bash first."
[ -d "${LW_LAB_ROOT}"     ] || _bail "Lab root does not exist: ${LW_LAB_ROOT}"
[ -f "${LW_DB}"            ] || _bail "Database not found: ${LW_DB}"

command -v sqlite3 >/dev/null 2>&1 || _bail "sqlite3 not found."
command -v rsync   >/dev/null 2>&1 || _bail "rsync not found."

# --- parse arguments ---
DRY_RUN=""
if [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN="--dry-run"
    echo "lw_backup: dry run (no files will be copied)"
fi

# --- safe database snapshot ---
# sqlite3 .backup uses the backup API -- produces a consistent
# copy even if another process is writing. The temporary copy
# lives in LAB_ROOT and is cleaned up after sync.
DB_SNAPSHOT="${LW_LAB_ROOT}/lab_backup.db"

if [ -z "${DRY_RUN}" ]; then
    sqlite3 "${LW_DB}" ".backup '${DB_SNAPSHOT}'"
    echo "lw_backup: database snapshot created"
else
    echo "lw_backup: would create database snapshot"
fi

# --- sync files (additive, no deletes) ---
# --no-perms --no-owner --no-group: avoids WSL permission
# issues on the Windows filesystem.
# --exclude lab_backup.db: temporary snapshot, synced separately.
# --exclude lab.db: replaced by the safe snapshot.
rsync -av --no-perms --no-owner --no-group \
    --exclude='lab_backup.db' \
    --exclude='lab.db' \
    ${DRY_RUN} \
    "${LW_LAB_ROOT}/" "${LW_BACKUP_DIR}/"

# --- copy safe database snapshot ---
if [ -z "${DRY_RUN}" ]; then
    cp "${DB_SNAPSHOT}" "${LW_BACKUP_DIR}/lab.db"
    rm "${DB_SNAPSHOT}"
    echo "lw_backup: clean database copy synced"
else
    echo "lw_backup: would copy database snapshot to ${LW_BACKUP_DIR}/lab.db"
fi

# --- summary ---
echo "lw_backup: complete -> ${LW_BACKUP_DIR}"
