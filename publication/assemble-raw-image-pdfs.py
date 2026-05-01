# /// script
# dependencies = [
#   "pypdf",
# ]
# requires-python = ">=3.9"
# ///

"""
Assemble annotated raw image PDFs into single compiled PDF for PLOS ONE submission.

Purpose: Merge individual exported raw image PDFs (from Illustrator) into
S1_raw_images.pdf for PLOS ONE Supporting Information and Zenodo repository submission.

Only blot and gel figures are exported and assembled. Figures that are not blots or gels
(e.g., Fig2, S2, S3) are intentionally absent from the input set.

Usage:
    uv run assemble_raw_images.py

Requirements:
    - Input: pdf_export/ directory with PDFs matching pattern {FigureID}_src-panel_raw-image_{NN}.pdf
    - Output: S1_raw_images.pdf in current working directory

Author: Luis Martinez
Date: April 30, 2026
Project: LEMR Publication Bypass ORC4R Dataset
"""

import hashlib
import pathlib
import re
import sys

import pypdf

# ===========================================================================
# Configuration
# ===========================================================================

INPUT_DIRECTORY = pathlib.Path("pdf_export")
OUTPUT_FILENAME = pathlib.Path("S1_raw_images.pdf")
EXPECTED_COUNT = 33
FILENAME_PATTERN = re.compile(r"^(Fig\d+|S\d+)_src-panel_raw-image_\d{2}\.pdf$")
EXPECTED_FIGURE_IDS = {"Fig1", "Fig3", "Fig4", "Fig5", "Fig6", "S1", "S4", "S5", "S6"}
SIZE_WARNING_MEGABYTES = 20
MANIFEST_FILENAME = pathlib.Path("S1_raw_images_manifest.txt")
DRY_RUN = False

# ===========================================================================
# Phase 1: Scan input directory
# ===========================================================================

print("=" * 60)
print("Phase 1: Scanning input directory")
print("=" * 60)

if not INPUT_DIRECTORY.exists():
    print(f"ERROR: Input directory '{INPUT_DIRECTORY}' does not exist.")
    print("Create the directory and export PDFs from Illustrator before running this script.")
    sys.exit(1)

if not INPUT_DIRECTORY.is_dir():
    print(f"ERROR: '{INPUT_DIRECTORY}' exists but is not a directory.")
    sys.exit(1)

pdf_file_paths = list(INPUT_DIRECTORY.glob("*.pdf"))
pdf_file_count = len(pdf_file_paths)

print(f"Found {pdf_file_count} PDF file(s) in '{INPUT_DIRECTORY}/'.")

if pdf_file_count == 0:
    print("ERROR: No PDF files found. Export PDFs from Illustrator first.")
    sys.exit(1)

# ===========================================================================
# Phase 2: Validate naming pattern
# ===========================================================================

print()
print("=" * 60)
print("Phase 2: Validating filename patterns")
print("=" * 60)

non_matching_filenames = []
found_figure_ids = set()

for pdf_file_path in pdf_file_paths:
    filename = pdf_file_path.name
    match_result = FILENAME_PATTERN.match(filename)
    if match_result is None:
        non_matching_filenames.append(filename)
    else:
        figure_id = match_result.group(1)
        found_figure_ids.add(figure_id)


if len(non_matching_filenames) > 0:
    print(f"NOTE: {len(non_matching_filenames)} file(s) do not match the expected pattern and will be skipped:")
    print(f"  Pattern: {{FigureID}}_src-panel_raw-image_{{NN}}.pdf")
    for skipped_filename in non_matching_filenames:
        print(f"  - {skipped_filename}")

# Keep only the files that matched the pattern
pdf_file_paths = [
    pdf_file_path
    for pdf_file_path in pdf_file_paths
    if FILENAME_PATTERN.match(pdf_file_path.name) is not None
]
pdf_file_count = len(pdf_file_paths)

if pdf_file_count == 0:
    print("ERROR: No files matched the expected naming pattern.")
    sys.exit(1)

print(f"{pdf_file_count} file(s) match the expected pattern.")

# Check for missing figure IDs (warn, do not halt)
missing_figure_ids = EXPECTED_FIGURE_IDS - found_figure_ids
if len(missing_figure_ids) > 0:
    missing_sorted = sorted(missing_figure_ids)
    print(f"WARNING: No PDFs found for figure ID(s): {', '.join(missing_sorted)}")
    print("  Note: Only blot/gel figures are expected. This may be intentional.")

# Check for unexpected figure IDs (warn, do not halt)
unexpected_figure_ids = found_figure_ids - EXPECTED_FIGURE_IDS
if len(unexpected_figure_ids) > 0:
    unexpected_sorted = sorted(unexpected_figure_ids)
    print(f"WARNING: Unexpected figure ID(s) found: {', '.join(unexpected_sorted)}")

if len(missing_figure_ids) == 0 and len(unexpected_figure_ids) == 0:
    print(f"All expected figure IDs present: {', '.join(sorted(EXPECTED_FIGURE_IDS))}")

# ===========================================================================
# Phase 3: Validate count
# ===========================================================================

print()
print("=" * 60)
print("Phase 3: Validating file count")
print("=" * 60)

if pdf_file_count != EXPECTED_COUNT:
    print(f"ERROR: Expected {EXPECTED_COUNT} PDFs but found {pdf_file_count}.")
    print("Check that all Illustrator files have been exported to pdf_export/.")
    sys.exit(1)

print(f"File count matches expected: {pdf_file_count}/{EXPECTED_COUNT}.")

# ===========================================================================
# Phase 4: Validate readability and page count
# ===========================================================================

print()
print("=" * 60)
print("Phase 4: Validating PDF readability and page counts")
print("=" * 60)

for pdf_file_path in pdf_file_paths:
    try:
        reader = pypdf.PdfReader(pdf_file_path)
        page_count = len(reader.pages)
    except Exception as read_error:
        print(f"ERROR: Cannot read '{pdf_file_path.name}': {read_error}")
        print("Re-export this file from Illustrator and try again.")
        sys.exit(1)

    if page_count != 1:
        print(f"ERROR: '{pdf_file_path.name}' has {page_count} pages (expected 1).")
        print("Each Illustrator file should have a single artboard exported as a single-page PDF.")
        sys.exit(1)

print("All PDFs are readable and contain exactly 1 page each.")

# ===========================================================================
# Phase 5: Sort
# ===========================================================================

print()
print("=" * 60)
print("Phase 5: Sorting files (lexicographic order)")
print("=" * 60)

sorted_pdf_file_paths = sorted(pdf_file_paths, key=lambda path: path.name)

print("Merge order:")
for page_number, pdf_file_path in enumerate(sorted_pdf_file_paths, start=1):
    print(f"  Page {page_number:2d}: {pdf_file_path.name}")

# ===========================================================================
# Phase 6: Merge
# ===========================================================================

print()
print("=" * 60)
print("Phase 6: Merging PDFs")
print("=" * 60)

if DRY_RUN:
    print("DRY RUN: Validation and sort order shown above. No output file written.")
    sys.exit(0)

writer = pypdf.PdfWriter()

for page_index, pdf_file_path in enumerate(sorted_pdf_file_paths):
    reader = pypdf.PdfReader(pdf_file_path)
    page = reader.pages[0]
    writer.add_page(page)
    current_page_number = page_index + 1
    if current_page_number % 10 == 0 or current_page_number == len(sorted_pdf_file_paths):
        print(f"  Merged {current_page_number}/{len(sorted_pdf_file_paths)} pages...")

with open(OUTPUT_FILENAME, "wb") as output_file:
    writer.write(output_file)

print(f"Output written to '{OUTPUT_FILENAME}'.")

# ===========================================================================
# Phase 7: Report
# ===========================================================================

print()
print("=" * 60)
print("Phase 7: Summary report")
print("=" * 60)

output_file_size_bytes = OUTPUT_FILENAME.stat().st_size
output_file_size_megabytes = output_file_size_bytes / (1024 * 1024)
total_pages_merged = len(sorted_pdf_file_paths)

print(f"Pages merged: {total_pages_merged}")
print(f"Output file:  {OUTPUT_FILENAME}")
print(f"File size:    {output_file_size_megabytes:.2f} MB")

if output_file_size_megabytes > SIZE_WARNING_MEGABYTES:
    print()
    print(f"WARNING: Output file exceeds {SIZE_WARNING_MEGABYTES} MB.")
    print("Upload to Zenodo and reference DOI in PLOS ONE Data Availability Statement.")

# ===========================================================================
# Phase 8: Generate manifest with SHA-256 hashes
# ===========================================================================

if not DRY_RUN:
    print()
    print("=" * 60)
    print("Phase 8: Generating manifest with SHA-256 hashes")
    print("=" * 60)

    manifest_lines = []
    manifest_lines.append("# S1_raw_images.pdf - Page-to-source manifest")
    manifest_lines.append(f"# Generated: {OUTPUT_FILENAME}")
    manifest_lines.append(f"# Total pages: {total_pages_merged}")
    manifest_lines.append("#")
    manifest_lines.append("# Page | Source Filename | SHA-256")
    manifest_lines.append("#" + "-" * 59)

    for page_index, pdf_file_path in enumerate(sorted_pdf_file_paths):
        page_number = page_index + 1

        sha256_hash = hashlib.sha256()
        with open(pdf_file_path, "rb") as hash_input_file:
            while True:
                chunk = hash_input_file.read(8192)
                if len(chunk) == 0:
                    break
                sha256_hash.update(chunk)
        hex_digest = sha256_hash.hexdigest()

        manifest_lines.append(f"{page_number:2d} | {pdf_file_path.name} | {hex_digest}")

    manifest_text = "\n".join(manifest_lines) + "\n"

    with open(MANIFEST_FILENAME, "w") as manifest_output_file:
        manifest_output_file.write(manifest_text)

    print(f"Manifest written to '{MANIFEST_FILENAME}'.")
    print(f"Contains {total_pages_merged} entries with SHA-256 hashes.")

print()
print("Done.")
