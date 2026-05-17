# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "pymupdf",
#   "python-docx",
#   "Pillow",
# ]
# ///

import os
import sys
from pathlib import Path
from collections import defaultdict

import fitz
import docx
from PIL import Image

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_DIR          = Path.home() / "lab" / "protocols"
SUPPORTED_EXTS     = {".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".docx"}
WARN_EXTS          = {".doc"}
TEMP_PREFIXES      = ("._", "~$", ".~lock.")
SYSTEM_FILES       = {".DS_Store", "desktop.ini", "Thumbs.db"}
LARGE_FILE_WARN_MB = 50

# Magic bytes: extension -> (offset, expected_bytes)
MAGIC = {
    ".pdf":  (0, b"%PDF"),
    ".docx": (0, b"PK\x03\x04"),
    ".png":  (0, b"\x89PNG\r\n\x1a\n"),
    ".jpg":  (0, b"\xff\xd8\xff"),
    ".jpeg": (0, b"\xff\xd8\xff"),
    ".tif":  (0, None),   # two valid headers, checked separately
    ".tiff": (0, None),
}
TIFF_HEADERS = (b"II*\x00", b"MM\x00*")

# =============================================================================
# PREFLIGHT
# =============================================================================

if not INPUT_DIR.exists():
    print(f"ERROR: Input directory not found: {INPUT_DIR}")
    sys.exit(1)

if not INPUT_DIR.is_dir():
    print(f"ERROR: {INPUT_DIR} is not a directory.")
    sys.exit(1)

# =============================================================================
# COLLECT FILES
# =============================================================================

all_files = sorted(f for f in INPUT_DIR.iterdir() if f.is_file())

print(f"INFO: Inspecting {INPUT_DIR}")
print(f"INFO: {len(all_files)} files found")
print()

# =============================================================================
# MAIN INSPECTION LOOP
# =============================================================================

# Per-extension accumulators
ext_counts = defaultdict(int)
ext_bytes  = defaultdict(int)

# Issue tracking
issues = []   # (filename, issue_description)

for f in all_files:

    name = f.name
    ext  = f.suffix.lower()
    size = f.stat().st_size

    ext_counts[ext or "(no extension)"] += 1
    ext_bytes[ext or "(no extension)"]  += size

    # --- Temp / system artifact check ---
    if any(name.startswith(p) for p in TEMP_PREFIXES) or name in SYSTEM_FILES:
        issues.append((name, "system artifact - safe to delete"))
        print(f"  SYSTEM  : {name}")
        continue

    # --- Zero-byte check ---
    if size == 0:
        issues.append((name, "zero bytes - likely a Dropbox placeholder not downloaded"))
        print(f"  EMPTY   : {name}")
        continue

    # --- Large file warning ---
    size_mb = size / (1024 * 1024)
    if size_mb > LARGE_FILE_WARN_MB:
        issues.append((name, f"large file ({size_mb:.1f} MB) - OCR may be slow"))
        print(f"  LARGE   : {name}  ({size_mb:.1f} MB)")

    # --- Unsupported / warned extensions ---
    if ext in WARN_EXTS:
        issues.append((name, ".doc format - not supported, needs manual conversion"))
        print(f"  DOC     : {name}")
        continue

    if ext not in SUPPORTED_EXTS:
        print(f"  UNKNOWN : {name}  ({ext})")
        continue

    # --- Magic byte check ---
    magic_ok = True
    try:
        with open(f, "rb") as fh:
            header = fh.read(8)

        if ext in (".tif", ".tiff"):
            if not any(header.startswith(h) for h in TIFF_HEADERS):
                issues.append((name, f"magic mismatch - header {header[:4]!r} not a valid TIFF"))
                print(f"  MAGIC?  : {name}")
                magic_ok = False
        else:
            expected = MAGIC[ext][1]
            if not header.startswith(expected):
                issues.append((name, f"magic mismatch - header {header[:4]!r} does not match {ext}"))
                print(f"  MAGIC?  : {name}")
                magic_ok = False

    except OSError as e:
        issues.append((name, f"could not read file: {e}"))
        print(f"  UNREAD  : {name}")
        magic_ok = False

    if not magic_ok:
        continue

    # --- Structural integrity check ---
    ok = True
    detail = ""

    try:
        if ext == ".pdf":
            doc = fitz.open(f)
            pages = doc.page_count
            doc.close()
            if pages == 0:
                ok = False
                detail = "opened but has 0 pages"
            else:
                detail = f"{pages}p"

        elif ext == ".docx":
            docx.Document(f)
            detail = "ok"

        elif ext in (".png", ".jpg", ".jpeg", ".tif", ".tiff"):
            img = Image.open(f)
            img.verify()
            detail = "ok"

    except Exception as e:
        ok = False
        detail = str(e)

    if ok:
        size_kb = size / 1024
        print(f"  OK      : {name:<50}  {size_kb:>8.1f} KB  {detail}")
    else:
        issues.append((name, f"structural check failed: {detail}"))
        print(f"  CORRUPT?: {name}  - {detail}")

# =============================================================================
# DUPLICATE FILENAME CHECK
# =============================================================================

seen_names = defaultdict(list)
for f in all_files:
    seen_names[f.name].append(f)

duplicates = {n: paths for n, paths in seen_names.items() if len(paths) > 1}
if duplicates:
    print()
    print("WARN: Duplicate filenames detected:")
    for name, paths in duplicates.items():
        for p in paths:
            issues.append((name, "duplicate filename in directory"))
        print(f"  {name}  ({len(paths)} copies)")

# =============================================================================
# SUMMARY
# =============================================================================

print()
print("-" * 60)
print("File breakdown by extension:")
for ext, count in sorted(ext_counts.items()):
    total_kb = ext_bytes[ext] / 1024
    print(f"  {ext:<12} {count:>4} file(s)   {total_kb:>8.1f} KB")

total_size_kb = sum(ext_bytes.values()) / 1024
print(f"  {'TOTAL':<12} {len(all_files):>4} file(s)   {total_size_kb:>8.1f} KB")
print()

if issues:
    print(f"Issues found: {len(issues)}")
    for name, desc in issues:
        print(f"  ! {name}: {desc}")
else:
    print("No issues found. Directory looks clean.")

print("-" * 60)
