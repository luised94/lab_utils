# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "pymupdf",
#   "python-docx",
#   "pytesseract",
#   "Pillow",
# ]
# ///

import os
import sys
import shutil
import time
from pathlib import Path

import fitz
import docx
import pytesseract
from PIL import Image

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_DIR              = Path.home() / "lab" / "protocols"
OUTPUT_DIR_SUFFIX      = "_extracted"
SCANNED_PAGE_THRESHOLD = 50
TESSERACT_LANG         = "eng"
SUPPORTED_EXTENSIONS   = {".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".docx"}
WARN_EXTENSIONS        = {".doc"}

# =============================================================================
# DERIVED VALUES
# =============================================================================

output_dir = INPUT_DIR.parent / (INPUT_DIR.name + OUTPUT_DIR_SUFFIX)

# =============================================================================
# PREFLIGHT
# =============================================================================

if shutil.which("tesseract") is None:
    print("ERROR: tesseract-ocr not found on PATH.")
    print("       Install with: sudo apt install -y tesseract-ocr")
    sys.exit(1)

if not INPUT_DIR.exists():
    print(f"ERROR: Input directory not found: {INPUT_DIR}")
    sys.exit(1)

if not INPUT_DIR.is_dir():
    print(f"ERROR: {INPUT_DIR} exists but is not a directory.")
    sys.exit(1)

if not os.access(INPUT_DIR, os.R_OK):
    print(f"ERROR: No read permission on {INPUT_DIR}")
    sys.exit(1)

if not os.access(output_dir.parent, os.W_OK):
    print(f"ERROR: No write permission in {output_dir.parent}")
    print(f"       Cannot create output directory.")
    sys.exit(1)

output_dir.mkdir(exist_ok=True)

print(f"INFO: Input  : {INPUT_DIR}")
print(f"INFO: Output : {output_dir}")
