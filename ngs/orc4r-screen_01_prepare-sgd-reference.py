# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
# ]
# ///
"""
Download and prepare the yeast gene coordinate reference from SGD.

Downloads SGD_features.tab from the Saccharomyces Genome Database, filters for
gene-like features on nuclear chromosomes, maps chromosome identifiers to the
chrI-chrXVI format used in BAM files, and saves as a TSV with a companion
metadata file recording the source URL, download date, and SHA256 hash.

The output files are used by analyze_mutations.py for gene annotation.

Usage:
    uv run prepare_reference.py            (download and process)
    uv run prepare_reference.py --force    (re-download even if files exist)
    uv run prepare_reference.py --help

Output:
    reference/gene_coordinates.tsv
    reference/gene_coordinates.meta.json
"""

import argparse
import datetime
import hashlib
import json
import pathlib
import sys
import urllib.request

import pandas


# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIRECTORY = pathlib.Path(__file__).resolve().parent
REFERENCE_DIRECTORY = SCRIPT_DIRECTORY / "reference"

GENE_COORDINATES_FILE = REFERENCE_DIRECTORY / "gene_coordinates.tsv"
GENE_COORDINATES_METADATA_FILE = REFERENCE_DIRECTORY / "gene_coordinates.meta.json"

SGD_FEATURES_URL = (
    "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
)

# SGD_features.tab is headerless. These are the column names in positional order.
# Reference: https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.README
SGD_COLUMN_NAMES = [
    "sgdid",
    "feature_type",
    "feature_qualifier",
    "feature_name",  # systematic name, e.g. YAL001C
    "standard_name",  # common gene name, e.g. TFC3
    "alias",  # pipe-separated aliases
    "parent_feature_name",
    "secondary_sgdid",
    "chromosome",  # 1-16, Mito, or 2-micron
    "start",
    "stop",
    "strand",  # W (Watson/+) or C (Crick/-)
    "genetic_position",
    "coordinate_version",
    "sequence_version",
    "description",
]

# Map SGD numeric chromosome identifiers to BAM-style roman numeral names
CHROMOSOME_NUMBER_TO_NAME = {
    "1": "chrI",
    "2": "chrII",
    "3": "chrIII",
    "4": "chrIV",
    "5": "chrV",
    "6": "chrVI",
    "7": "chrVII",
    "8": "chrVIII",
    "9": "chrIX",
    "10": "chrX",
    "11": "chrXI",
    "12": "chrXII",
    "13": "chrXIII",
    "14": "chrXIV",
    "15": "chrXV",
    "16": "chrXVI",
}

# Feature types that represent genes or gene-like genomic elements.
# These are the features useful for annotating point mutations.
INCLUDED_FEATURE_TYPES = {
    "ORF",
    "tRNA gene",
    "rRNA gene",
    "snRNA gene",
    "snoRNA gene",
    "ncRNA gene",
    "pseudogene",
    "transposable_element_gene",
}

# Columns to keep in the output TSV, in this order
OUTPUT_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "gene_symbol",
    "gene_name",
    "gene_alias",
    "feature_type",
    "feature_qualifier",
    "description",
]


# =============================================================================
# Argument parsing
# =============================================================================

argument_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
argument_parser.add_argument(
    "--force",
    action="store_true",
    help="Re-download and overwrite existing reference files.",
)
arguments = argument_parser.parse_args()


# =============================================================================
# Check if reference files already exist
# =============================================================================

if (
    not arguments.force
    and GENE_COORDINATES_FILE.exists()
    and GENE_COORDINATES_METADATA_FILE.exists()
):
    print(f"Reference files already exist:")
    print(f"  {GENE_COORDINATES_FILE}")
    print(f"  {GENE_COORDINATES_METADATA_FILE}")
    print(f"Run with --force to re-download and overwrite.")
    sys.exit(0)


# =============================================================================
# Step 1: Download SGD_features.tab
# =============================================================================

print("=" * 60)
print("STEP 1: DOWNLOAD SGD_features.tab")
print("=" * 60)
print(f"  Source: {SGD_FEATURES_URL}")

download_date = datetime.date.today().isoformat()

try:
    response = urllib.request.urlopen(SGD_FEATURES_URL, timeout=60)
    raw_bytes = response.read()
    raw_text = raw_bytes.decode("utf-8")
except urllib.error.URLError as download_error:
    sys.exit(f"[ERROR] Download failed: {download_error}")

line_count = raw_text.count("\n")
print(f"  Downloaded {len(raw_bytes):,} bytes ({line_count:,} lines)")
print(f"  Date: {download_date}")
print()


# =============================================================================
# Step 2: Parse and filter
# =============================================================================

print("=" * 60)
print("STEP 2: PARSE AND FILTER")
print("=" * 60)

# Parse the headerless tab-separated file
# StringIO avoids writing a temp file to disk
import io

all_features = pandas.read_csv(
    io.StringIO(raw_text),
    sep="\t",
    header=None,
    names=SGD_COLUMN_NAMES,
    dtype=str,
    na_values=[""],
    keep_default_na=False,
)
print(f"  Total features in SGD_features.tab: {len(all_features):,}")

# Filter to nuclear chromosomes only (exclude Mito, 2-micron, and any blanks)
nuclear_chromosomes = set(CHROMOSOME_NUMBER_TO_NAME.keys())
is_nuclear = all_features["chromosome"].isin(nuclear_chromosomes)
nuclear_features = all_features[is_nuclear].copy()
print(f"  Nuclear chromosome features: {len(nuclear_features):,}")

# Filter to gene-like feature types
is_gene_like = nuclear_features["feature_type"].isin(INCLUDED_FEATURE_TYPES)
gene_features = nuclear_features[is_gene_like].copy()
print(f"  Gene-like features retained: {len(gene_features):,}")

# Report feature type breakdown
feature_type_counts = gene_features["feature_type"].value_counts()
for feature_type, count in feature_type_counts.items():
    print(f"    {feature_type:<30} {count:>5}")

assert len(gene_features) > 0, (
    "No gene features found after filtering. Check SGD file format."
)
print()


# =============================================================================
# Step 3: Transform to output format
# =============================================================================

print("=" * 60)
print("STEP 3: TRANSFORM TO OUTPUT FORMAT")
print("=" * 60)

# Map chromosome numbers to chrI-chrXVI names
gene_features["chromosome"] = gene_features["chromosome"].map(CHROMOSOME_NUMBER_TO_NAME)
unmapped_chromosomes = gene_features["chromosome"].isna().sum()
assert unmapped_chromosomes == 0, (
    f"{unmapped_chromosomes} rows have unmapped chromosome values"
)
print(f"  [OK] Chromosome identifiers mapped to chrI-chrXVI format")

# Convert coordinates to integers
gene_features["start"] = gene_features["start"].astype(int)
gene_features["stop"] = gene_features["stop"].astype(int)

# SGD coordinates: some features have start > stop (Crick strand).
# Normalize so start is always the lower coordinate for interval lookups.
needs_swap = gene_features["start"] > gene_features["stop"]
swap_count = needs_swap.sum()
if swap_count > 0:
    gene_features.loc[needs_swap, ["start", "stop"]] = gene_features.loc[
        needs_swap, ["stop", "start"]
    ].values
    print(f"  [OK] Swapped start/stop for {swap_count:,} Crick-strand features")

# Rename columns to match what analyze_mutations.py expects
gene_features = gene_features.rename(
    columns={
        "stop": "end",
        "standard_name": "gene_symbol",
        "feature_name": "gene_name",
        "alias": "gene_alias",
    }
)

# Select and order output columns
gene_coordinates = gene_features[OUTPUT_COLUMNS].copy()

# Sort by chromosome (in biological order) then by start position
chromosome_sort_order = {
    name: index for index, name in enumerate(CHROMOSOME_NUMBER_TO_NAME.values())
}
gene_coordinates = gene_coordinates.sort_values(
    by=["chromosome", "start"],
    key=lambda column: (
        column.map(chromosome_sort_order) if column.name == "chromosome" else column
    ),
    ignore_index=True,
)

print(f"  [OK] {len(gene_coordinates):,} gene annotations ready")
print()


# =============================================================================
# Step 4: Save reference files
# =============================================================================

print("=" * 60)
print("STEP 4: SAVE REFERENCE FILES")
print("=" * 60)

REFERENCE_DIRECTORY.mkdir(exist_ok=True)

# Write TSV
gene_coordinates.to_csv(
    str(GENE_COORDINATES_FILE),
    sep="\t",
    index=False,
)
print(f"  Saved: {GENE_COORDINATES_FILE}")

# Compute SHA256 hash of the written file
file_bytes = GENE_COORDINATES_FILE.read_bytes()
file_sha256 = hashlib.sha256(file_bytes).hexdigest()

# Write metadata
metadata = {
    "source_url": SGD_FEATURES_URL,
    "download_date": download_date,
    "sha256": file_sha256,
    "feature_types_included": sorted(INCLUDED_FEATURE_TYPES),
    "gene_count": len(gene_coordinates),
    "note": (
        "Generated by prepare_reference.py from SGD_features.tab. "
        "See https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.README "
        "for column definitions."
    ),
}

GENE_COORDINATES_METADATA_FILE.write_text(json.dumps(metadata, indent=2) + "\n")
print(f"  Saved: {GENE_COORDINATES_METADATA_FILE}")
print(f"  SHA256: {file_sha256}")
print()

print("Done. Reference files are ready for analyze_mutations.py.")
