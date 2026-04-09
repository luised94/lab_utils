# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pysam",
#     "pandas",
#     "numpy",
# ]
# ///
"""
Identify potential causative mutations from yeast whole-genome sequencing data.

Compares suppressor ("hit") samples against wild-type ("ctrl") samples across
multiple experiments. Candidate mutations are identified where:
  1) More than 10 reads map to the base (both samples)
  2) More than 90% of reads agree on the nucleotide (both samples)
  3) The nucleotide differs between wild-type and suppressor

Candidate positions are then annotated with gene names using a local copy
of the SGD gene coordinate reference (see reference/ directory).

Usage:
    uv run analyze_mutations.py
    uv run analyze_mutations.py --help

Prerequisites:
    uv run prepare_reference.py   (one-time setup for gene annotations)

Expected directory structure:
    analyze_mutations.py
    prepare_reference.py
    reference/
        gene_coordinates.tsv       (from prepare_reference.py)
        gene_coordinates.meta.json (from prepare_reference.py)
    data/
        Exp_1_hit_sorted.bam       (+ .bai index, auto-generated if missing)
        Exp_1_ctrl_sorted.bam
        ...through Exp_9...
    output/                        (created automatically)
        pickles/                   (intermediate cached data)
        candidate_mutations_with_genes.csv
"""

import argparse
import collections
import hashlib
import json
import pathlib
import sys

import numpy
import pandas
import pysam


# =============================================================================
# Configuration - all thresholds, paths, and constants in one place
# =============================================================================

SCRIPT_DIRECTORY = pathlib.Path(__file__).resolve().parent

DATA_DIRECTORY = SCRIPT_DIRECTORY / "data"
OUTPUT_DIRECTORY = SCRIPT_DIRECTORY / "output"
PICKLE_DIRECTORY = OUTPUT_DIRECTORY / "pickles"
REFERENCE_DIRECTORY = SCRIPT_DIRECTORY / "reference"

GENE_COORDINATES_FILE = REFERENCE_DIRECTORY / "gene_coordinates.tsv"
GENE_COORDINATES_METADATA_FILE = REFERENCE_DIRECTORY / "gene_coordinates.meta.json"

EXPERIMENT_NUMBERS = list(range(1, 10))  # experiments 1 through 9

HIT_BAM_TEMPLATE = "Exp_{experiment_number}_hit_sorted.bam"
CTRL_BAM_TEMPLATE = "Exp_{experiment_number}_ctrl_sorted.bam"

# Canonical chromosome names used in DataFrames and the gene reference file.
# BAM files may use a different convention; this is resolved during preflight.
REFERENCE_CHROMOSOMES = [
    "chrI",
    "chrII",
    "chrIII",
    "chrIV",
    "chrV",
    "chrVI",
    "chrVII",
    "chrVIII",
    "chrIX",
    "chrX",
    "chrXI",
    "chrXII",
    "chrXIII",
    "chrXIV",
    "chrXV",
    "chrXVI",
]

ROMAN_NUMERALS = [
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "X",
    "XI",
    "XII",
    "XIII",
    "XIV",
    "XV",
    "XVI",
]

# Known BAM chromosome naming conventions and their mapping to REFERENCE_CHROMOSOMES.
# Checked in order during preflight if the BAM files don't use chrI-chrXVI directly.
KNOWN_CHROMOSOME_CONVENTIONS = {
    "numeric": {
        str(i): f"chr{roman}" for i, roman in enumerate(ROMAN_NUMERALS, start=1)
    },
    "roman_no_prefix": {roman: f"chr{roman}" for roman in ROMAN_NUMERALS},
    "chr_numeric": {
        f"chr{i}": f"chr{roman}" for i, roman in enumerate(ROMAN_NUMERALS, start=1)
    },
}

# Filtering thresholds (must match manuscript methods section)
MINIMUM_READ_COVERAGE = 10  # reads mapped to the base
MINIMUM_CONFIDENCE_PERCENT = 90  # percent of reads agreeing on the nucleotide

OUTPUT_CSV_FILENAME = "candidate_mutations_with_genes.csv"


# =============================================================================
# Argument parsing (--help only, no arguments required)
# =============================================================================

argument_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
argument_parser.parse_args()


# =============================================================================
# Preflight checks
# =============================================================================

print("=" * 60)
print("PREFLIGHT CHECKS")
print("=" * 60)

# -- Check data directory exists --
if not DATA_DIRECTORY.is_dir():
    sys.exit(
        f"[ERROR] Data directory not found: {DATA_DIRECTORY}\n"
        f"  Create it and place sorted BAM files inside."
    )
print(f"  [OK] Data directory: {DATA_DIRECTORY}")

# -- Check all BAM files exist --
missing_bam_files = []
for experiment_number in EXPERIMENT_NUMBERS:
    for template in [HIT_BAM_TEMPLATE, CTRL_BAM_TEMPLATE]:
        bam_path = DATA_DIRECTORY / template.format(experiment_number=experiment_number)
        if not bam_path.exists():
            missing_bam_files.append(bam_path.name)

if len(missing_bam_files) > 0:
    sys.exit(
        f"[ERROR] Missing BAM files in {DATA_DIRECTORY}:\n"
        + "\n".join(f"  - {filename}" for filename in missing_bam_files)
    )
print(f"  [OK] All {len(EXPERIMENT_NUMBERS) * 2} BAM files found")

# -- Index BAM files if .bai is missing (pysam requires indexes for pileup) --
newly_indexed_count = 0
for experiment_number in EXPERIMENT_NUMBERS:
    for template in [HIT_BAM_TEMPLATE, CTRL_BAM_TEMPLATE]:
        bam_path = DATA_DIRECTORY / template.format(experiment_number=experiment_number)
        bai_path = pathlib.Path(str(bam_path) + ".bai")
        if not bai_path.exists():
            print(f"  [INDEXING] {bam_path.name} ...")
            pysam.index(str(bam_path))
            assert bai_path.exists(), f"Indexing failed for {bam_path.name}"
            newly_indexed_count += 1

if newly_indexed_count > 0:
    print(f"  [OK] Indexed {newly_indexed_count} BAM files")
else:
    print(f"  [OK] All BAM index files (.bai) present")

# -- Detect BAM chromosome naming convention --
# Open one BAM header and compare its contig names against our expected format.
# BAM files commonly use chrI, 1, I, or chr1 for yeast chromosomes.
sample_bam_path = DATA_DIRECTORY / HIT_BAM_TEMPLATE.format(
    experiment_number=EXPERIMENT_NUMBERS[0]
)
sample_bam_file = pysam.AlignmentFile(str(sample_bam_path), mode="rb")
bam_contig_names = set(sample_bam_file.references)
sample_bam_file.close()

if set(REFERENCE_CHROMOSOMES).issubset(bam_contig_names):
    # BAM already uses chrI-chrXVI - identity mapping
    bam_to_reference_chromosome_map = {
        chromosome: chromosome for chromosome in REFERENCE_CHROMOSOMES
    }
    print(f"  [OK] BAM chromosome names match reference format (chrI-chrXVI)")
else:
    # Try each known alternative convention
    bam_to_reference_chromosome_map = None
    matched_convention_name = None

    for convention_name, convention_mapping in KNOWN_CHROMOSOME_CONVENTIONS.items():
        if set(convention_mapping.keys()).issubset(bam_contig_names):
            bam_to_reference_chromosome_map = convention_mapping
            matched_convention_name = convention_name
            break

    if bam_to_reference_chromosome_map is None:
        # Show what we found to help the user diagnose
        sample_contigs = sorted(bam_contig_names)[:20]
        sys.exit(
            f"[ERROR] Unrecognized chromosome names in BAM files.\n"
            f"  Found in BAM header: {sample_contigs}\n"
            f"  Expected one of these conventions:\n"
            f"    chrI, chrII, ... chrXVI\n"
            f"    1, 2, ... 16\n"
            f"    I, II, ... XVI\n"
            f"    chr1, chr2, ... chr16\n"
            f"  Check that the BAM files are aligned to a yeast reference genome."
        )

    print(
        f"  [OK] BAM uses '{matched_convention_name}' chromosome names "
        f"- auto-mapped to chrI-chrXVI"
    )

# Build the reverse mapping: reference name  BAM contig name (for pileup calls)
reference_to_bam_chromosome_map = {
    reference_name: bam_name
    for bam_name, reference_name in bam_to_reference_chromosome_map.items()
}

# -- Check gene coordinate reference file and verify integrity --
if not GENE_COORDINATES_FILE.exists():
    sys.exit(
        f"[ERROR] Gene coordinate reference not found: {GENE_COORDINATES_FILE}\n"
        f"  Run 'uv run prepare_reference.py' first to download from SGD."
    )
if not GENE_COORDINATES_METADATA_FILE.exists():
    sys.exit(
        f"[ERROR] Reference metadata not found: {GENE_COORDINATES_METADATA_FILE}\n"
        f"  Run 'uv run prepare_reference.py' to regenerate."
    )

reference_metadata = json.loads(GENE_COORDINATES_METADATA_FILE.read_text())
expected_sha256 = reference_metadata["sha256"]
actual_sha256 = hashlib.sha256(GENE_COORDINATES_FILE.read_bytes()).hexdigest()

if actual_sha256 != expected_sha256:
    sys.exit(
        f"[ERROR] Gene coordinate file integrity check failed.\n"
        f"  Expected SHA256: {expected_sha256}\n"
        f"  Actual SHA256:   {actual_sha256}\n"
        f"  The file may have been modified. Run 'uv run prepare_reference.py' to regenerate."
    )
print(
    f"  [OK] Gene coordinates verified (source: {reference_metadata['source_url']}, "
    f"downloaded: {reference_metadata['download_date']})"
)

# -- Create output directories --
OUTPUT_DIRECTORY.mkdir(exist_ok=True)
PICKLE_DIRECTORY.mkdir(exist_ok=True)
print(f"  [OK] Output directory: {OUTPUT_DIRECTORY}")
print()


# =============================================================================
# Load gene coordinate reference
#
# Pre-group by chromosome so annotation lookups in Step 3 are fast.
# =============================================================================

gene_coordinates = pandas.read_csv(
    GENE_COORDINATES_FILE,
    sep="\t",
    dtype={"chromosome": str, "start": int, "end": int},
)

assert len(gene_coordinates) > 0, "Gene coordinate file is empty."
print(f"  Loaded {len(gene_coordinates):,} gene annotations from reference file")

gene_coordinates_by_chromosome = {
    chromosome: chromosome_genes.reset_index(drop=True)
    for chromosome, chromosome_genes in gene_coordinates.groupby("chromosome")
}
print()


# =============================================================================
# Step 1: Extract sequencing data from BAM files
#
# For each experiment, read the hit and ctrl BAM files, extract per-position
# consensus base calls and coverage, merge on shared positions, and cache
# the result as a pickle for fast re-runs.
#
# The BAM chromosome names may differ from the reference format (e.g. "1" vs
# "chrI"). We use reference_to_bam_chromosome_map for pileup calls and store
# the reference-format name in all DataFrames so downstream steps are consistent.
# =============================================================================

print("=" * 60)
print("STEP 1: EXTRACT SEQUENCING DATA FROM BAM FILES")
print("=" * 60)

for experiment_number in EXPERIMENT_NUMBERS:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"

    if pickle_path.exists():
        print(f"  [SKIP] Experiment {experiment_number}: cached pickle exists")
        continue

    print(f"  [PROCESSING] Experiment {experiment_number}")

    # One dataframe per sample (hit, ctrl), then merge.
    sample_dataframes = {}

    for suffix, bam_template in [
        ("hit", HIT_BAM_TEMPLATE),
        ("ctrl", CTRL_BAM_TEMPLATE),
    ]:
        bam_path = DATA_DIRECTORY / bam_template.format(
            experiment_number=experiment_number
        )

        # Accumulator lists - one entry per genomic position
        chromosome_names = []
        base_positions = []
        read_coverages = []
        consensus_bases = []
        consensus_confidences = []

        bam_file = pysam.AlignmentFile(str(bam_path), mode="rb")

        for reference_chromosome in REFERENCE_CHROMOSOMES:
            bam_chromosome = reference_to_bam_chromosome_map[reference_chromosome]
            print(f"    {suffix:<4} | {reference_chromosome} (BAM: {bam_chromosome})")

            for pileup_column in bam_file.pileup(contig=bam_chromosome, truncate=True):
                # Gather the base call from every read overlapping this position.
                # query_position is None when the read has a deletion here.
                base_calls_at_position = []
                for pileup_read in pileup_column.pileups:
                    if pileup_read.query_position is not None:
                        base_calls_at_position.append(
                            pileup_read.alignment.query_sequence[
                                pileup_read.query_position
                            ]
                        )
                    else:
                        base_calls_at_position.append("N")

                # Consensus: most frequent base call and its proportion
                if len(base_calls_at_position) > 0:
                    base_counter = collections.Counter(base_calls_at_position)
                    most_common_base, most_common_count = base_counter.most_common(1)[0]
                    confidence_percent = round(
                        (most_common_count / len(base_calls_at_position)) * 100, 2
                    )
                else:
                    most_common_base = "N"
                    confidence_percent = numpy.nan

                # Always store the reference-format chromosome name
                chromosome_names.append(reference_chromosome)
                base_positions.append(pileup_column.pos)
                read_coverages.append(pileup_column.n)
                consensus_bases.append(most_common_base)
                consensus_confidences.append(confidence_percent)

        bam_file.close()

        sample_dataframe = pandas.DataFrame(
            {
                "chromosome": chromosome_names,
                "base_position": base_positions,
                f"base_value_{suffix}": consensus_bases,
                f"coverage_{suffix}": read_coverages,
                f"confidence_{suffix}": consensus_confidences,
            }
        )

        assert len(sample_dataframe) > 0, (
            f"No data extracted from {bam_path.name}. "
            f"Verify the BAM file is sorted by coordinate and properly indexed."
        )
        print(f"    {suffix:<4} | {len(sample_dataframe):,} positions extracted")

        sample_dataframes[suffix] = sample_dataframe

    # -- Merge hit and ctrl on shared genomic positions --
    # Inner join: only positions sequenced in BOTH samples are useful for comparison.
    hit_position_count = len(sample_dataframes["hit"])
    ctrl_position_count = len(sample_dataframes["ctrl"])

    merged_dataframe = pandas.merge(
        left=sample_dataframes["hit"],
        right=sample_dataframes["ctrl"],
        on=["chromosome", "base_position"],
        how="inner",
    )

    merged_position_count = len(merged_dataframe)
    print(
        f"    Merged: {hit_position_count:,} (hit) + {ctrl_position_count:,} (ctrl) "
        f" {merged_position_count:,} shared positions"
    )

    # Sanity check: large drop suggests misaligned references or truncated files
    larger_input_count = max(hit_position_count, ctrl_position_count)
    if merged_position_count < larger_input_count * 0.9:
        print(
            f"    [WARNING] Merge retained <90% of positions. "
            f"Verify both BAM files are aligned to the same reference genome."
        )

    merged_dataframe.to_pickle(str(pickle_path))
    print(f"    Cached: {pickle_path.name}")

print()


# =============================================================================
# Step 2: Filter for candidate mutations
#
# Apply the three published criteria symmetrically to both samples:
#   - Coverage must exceed MINIMUM_READ_COVERAGE in both hit and ctrl
#   - Confidence must exceed MINIMUM_CONFIDENCE_PERCENT in both hit and ctrl
#   - The consensus base must differ between hit and ctrl
# =============================================================================

print("=" * 60)
print("STEP 2: FILTER FOR CANDIDATE MUTATIONS")
print("=" * 60)
print(
    f"  Criteria: coverage > {MINIMUM_READ_COVERAGE}, "
    f"confidence > {MINIMUM_CONFIDENCE_PERCENT}%, "
    f"base differs between hit and ctrl"
)

all_candidate_dataframes = []

for experiment_number in EXPERIMENT_NUMBERS:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    passes_coverage_hit = merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE
    passes_coverage_ctrl = merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE
    passes_confidence_hit = (
        merged_dataframe["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT
    )
    passes_confidence_ctrl = (
        merged_dataframe["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT
    )
    bases_differ = (
        merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]
    )

    is_candidate_mutation = (
        passes_coverage_hit
        & passes_coverage_ctrl
        & passes_confidence_hit
        & passes_confidence_ctrl
        & bases_differ
    )

    candidate_dataframe = merged_dataframe[is_candidate_mutation].copy()
    candidate_dataframe["experiment"] = f"Exp_{experiment_number}"

    print(
        f"  Experiment {experiment_number}: {len(candidate_dataframe):,} candidate positions"
    )
    all_candidate_dataframes.append(candidate_dataframe)

all_candidates = pandas.concat(all_candidate_dataframes, ignore_index=True)
total_candidate_count = len(all_candidates)
print(f"  Total candidates across all experiments: {total_candidate_count:,}")

if total_candidate_count == 0:
    print("  [NOTE] No candidates found. Adjust thresholds or check data quality.")
    print("  Writing empty output file.")
    all_candidates.to_csv(str(OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME), index=False)
    sys.exit(0)

print()


# =============================================================================
# Step 3: Annotate candidate positions with gene names
#
# For each candidate position, find which gene(s) in the SGD reference overlap
# that chromosomal coordinate. A gene overlaps when: gene_start <= position < gene_end.
#
# Positions falling in intergenic regions will have no gene annotation.
# When multiple genes overlap a position, the first (by start coordinate) is kept.
# =============================================================================

print("=" * 60)
print("STEP 3: ANNOTATE WITH GENE NAMES")
print("=" * 60)

all_candidates["gene_symbol"] = pandas.NA
all_candidates["gene_name"] = pandas.NA
all_candidates["gene_alias"] = pandas.NA
all_candidates["gene_feature_type"] = pandas.NA
all_candidates["gene_description"] = pandas.NA

intergenic_count = 0

for candidate_index in range(total_candidate_count):
    chromosome = all_candidates.at[candidate_index, "chromosome"]
    base_position = all_candidates.at[candidate_index, "base_position"]
    experiment = all_candidates.at[candidate_index, "experiment"]

    print(
        f"  [{candidate_index + 1}/{total_candidate_count}] "
        f"{experiment} | {chromosome}:{base_position}",
        end="",
    )

    # Look up genes overlapping this position on this chromosome
    chromosome_genes = gene_coordinates_by_chromosome.get(
        chromosome, pandas.DataFrame()
    )

    if len(chromosome_genes) == 0:
        print(f"  [WARNING] no genes loaded for {chromosome}")
        intergenic_count += 1
        continue

    overlapping_genes = chromosome_genes[
        (chromosome_genes["start"] <= base_position)
        & (chromosome_genes["end"] > base_position)
    ]

    if len(overlapping_genes) == 0:
        print("  (intergenic)")
        intergenic_count += 1
        continue

    # Take the first overlapping gene (sorted by start coordinate in the reference file)
    first_gene = overlapping_genes.iloc[0]
    all_candidates.at[candidate_index, "gene_symbol"] = first_gene.get(
        "gene_symbol", pandas.NA
    )
    all_candidates.at[candidate_index, "gene_name"] = first_gene.get(
        "gene_name", pandas.NA
    )
    all_candidates.at[candidate_index, "gene_alias"] = first_gene.get(
        "gene_alias", pandas.NA
    )
    all_candidates.at[candidate_index, "gene_feature_type"] = first_gene.get(
        "feature_type", pandas.NA
    )
    all_candidates.at[candidate_index, "gene_description"] = first_gene.get(
        "description", pandas.NA
    )

    display_name = (
        first_gene.get("gene_symbol")
        or first_gene.get("gene_name")
        or first_gene.get("gene_alias")
        or "unnamed"
    )
    print(f"  {display_name}")

if intergenic_count > 0:
    print(f"  [NOTE] {intergenic_count} positions in intergenic regions")

print()


# =============================================================================
# Step 4: Save results
#
# Sort by chromosome and position for readability.
# Column order: experiment context  genomic location  sequencing data  annotation.
# =============================================================================

print("=" * 60)
print("STEP 4: SAVE RESULTS")
print("=" * 60)

output_column_order = [
    "experiment",
    "chromosome",
    "base_position",
    "base_value_hit",
    "base_value_ctrl",
    "coverage_hit",
    "coverage_ctrl",
    "confidence_hit",
    "confidence_ctrl",
    "gene_symbol",
    "gene_name",
    "gene_alias",
    "gene_feature_type",
    "gene_description",
]

chromosome_sort_order = {
    chromosome: index for index, chromosome in enumerate(REFERENCE_CHROMOSOMES)
}

all_candidates_sorted = all_candidates[output_column_order].sort_values(
    by=["experiment", "chromosome", "base_position"],
    key=lambda column: (
        column.map(chromosome_sort_order) if column.name == "chromosome" else column
    ),
    ignore_index=True,
)

output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
all_candidates_sorted.to_csv(str(output_csv_path), index=False)
print(f"  Saved {total_candidate_count:,} candidates to: {output_csv_path}")
print()


# =============================================================================
# Summary: genes appearing across multiple experiments
#
# Genes mutated in more than one independent experiment are the strongest
# candidates for causative mutations.
# =============================================================================

print("=" * 60)
print("SUMMARY: RECURRENT GENES ACROSS EXPERIMENTS")
print("=" * 60)

annotated_candidates = all_candidates_sorted[
    all_candidates_sorted["gene_symbol"].notna()
]

if len(annotated_candidates) == 0:
    print("  No annotated genes to summarize.")
else:
    gene_experiment_counts = (
        annotated_candidates.groupby("gene_symbol")["experiment"]
        .nunique()
        .sort_values(ascending=False)
    )

    recurrent_genes = gene_experiment_counts[gene_experiment_counts > 1]

    if len(recurrent_genes) == 0:
        print("  No genes found in more than one experiment.")
    else:
        print(f"  {len(recurrent_genes)} gene(s) mutated in multiple experiments:\n")
        for gene_symbol, experiment_count in recurrent_genes.items():
            gene_rows = annotated_candidates[
                annotated_candidates["gene_symbol"] == gene_symbol
            ]
            experiments_list = ", ".join(sorted(gene_rows["experiment"].unique()))
            print(
                f"    {gene_symbol:<15} {experiment_count} experiments ({experiments_list})"
            )

print()
print("Done.")
