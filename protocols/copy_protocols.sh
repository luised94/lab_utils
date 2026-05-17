#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# ENVIRONMENT CHECKS
# =============================================================================

if [ -z "${MC_WINDOWS_USER:-}" ]; then
    echo "ERROR: MC_WINDOWS_USER is not set."
    echo "       export MC_WINDOWS_USER=YourWindowsUsername"
    echo "       Add to ~/.bashrc to persist across sessions."
    exit 1
fi

# =============================================================================
# CONFIGURATION
# =============================================================================

SOURCE_DIR="/mnt/c/Users/$MC_WINDOWS_USER/MIT Dropbox/Luis Martinez/Lab/Protocols"
DEST_DIR="$HOME/lab/protocols"

# =============================================================================
# DERIVED / PREFLIGHT
# =============================================================================

if [ ! -d "$SOURCE_DIR" ]; then
    echo "ERROR: Source directory not found: $SOURCE_DIR"
    exit 1
fi

mkdir -p "$DEST_DIR"

echo "INFO: Source : $SOURCE_DIR"
echo "INFO: Dest   : $DEST_DIR"
echo "---"

# =============================================================================
# MAIN
# =============================================================================

copied=0
skipped_temp=0
skipped_collision=0

while IFS= read -r -d '' src; do

    filename=$(basename "$src")

    # Skip system artifacts and Office temp/lock files
    case "$filename" in
        ._*|"~$"*|.~lock.*|.DS_Store|desktop.ini|Thumbs.db)
            echo "SKIP (system): $filename"
            skipped_temp=$(( skipped_temp + 1 ))
            continue
            ;;
        *.mp4|*.m)
            echo "SKIP (wrong extension): $filename"
            skipped_temp=$(( skipped_temp + 1))
            continue
            ;;
    esac

    dest="$DEST_DIR/$filename"

    # Collision: prefix with immediate parent directory name
    if [ -e "$dest" ]; then
        parent=$(basename "$(dirname "$src")")
        dest="$DEST_DIR/${parent}__${filename}"

        # Second collision: skip and warn - manual resolution needed
        if [ -e "$dest" ]; then
            echo "WARN (collision): $filename already exists as ${parent}__${filename}, skipping."
            skipped_collision=$(( skipped_collision + 1 ))
            continue
        fi

        echo "RENAMED : $filename -> ${parent}__${filename}"
    fi

    cp "$src" "$dest"
    echo "OK: $filename"
    copied=$(( copied + 1 ))

done < <(find "$SOURCE_DIR" -type f -print0)

# =============================================================================
# SUMMARY
# =============================================================================

echo "---"
echo "Done. Copied: $copied | Temp skipped: $skipped_temp | Collision skipped: $skipped_collision"
