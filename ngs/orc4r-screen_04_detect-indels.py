# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pysam",
#     "pandas",
#     "numpy",
# ]
# ///
"""
Detect insertion/deletion (indel) mutations from yeast whole-genome sequencing data.

Companion to orc4r-screen_02_analyze-bam.py, which detects single nucleotide
variants (SNVs). This script targets experiments where the suppressor mutation
is an indel rather than a point mutation - identified by the coverage diagnostic
(orc4r-screen_03_coverage-diagnostic.py) as experiments with high-quality data
but zero SNV candidates and ~57% confidence at differing positions.

For each genomic position, counts reads showing insertions or deletions via
pysam's pileup interface, then identifies positions where the hit (suppressor)
sample has a significant indel signal absent in the ctrl (wild-type) sample.

Usage:
    uv run orc4r-screen_04_detect-indels.py
    uv run orc4r-screen_04_detect-indels.py --experiments 1,7,9
    uv run orc4r-screen_04_detect-indels.py --clean
    uv run orc4r-screen_04_detect-indels.py --help

Prerequisites:
    uv run orc4r-screen_01_prepare-reference.py   (gene annotations)
    BAM files in data/ directory (same as orc4r-screen_02_analyze-bam.py)
"""

import argparse
import collections
import hashlib
import json
import pathlib
import sys
import time

import numpy
import pandas
import pysam


# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIRECTORY = pathlib.Path(__file__).resolve().parent

DATA_DIRECTORY = SCRIPT_DIRECTORY / "data"
OUTPUT_DIRECTORY = SCRIPT_DIRECTORY / "output"
INDEL_PICKLE_DIRECTORY = OUTPUT_DIRECTORY / "pickles_indel"
REFERENCE_DIRECTORY = SCRIPT_DIRECTORY / "reference"

GENE_COORDINATES_FILE = REFERENCE_DIRECTORY / "gene_coordinates.tsv"
GENE_COORDINATES_METADATA_FILE = REFERENCE_DIRECTORY / "gene_coordinates.meta.json"

DEFAULT_EXPERIMENT_NUMBERS = list(range(1, 10))

HIT_BAM_TEMPLATE = "Exp_{experiment_number}_hit_sorted.bam"
CTRL_BAM_TEMPLATE = "Exp_{experiment_number}_ctrl_sorted.bam"

REFERENCE_CHROMOSOMES = [
    "chrI",    "chrII",   "chrIII",  "chrIV",
    "chrV",    "chrVI",   "chrVII",  "chrVIII",
    "chrIX",   "chrX",    "chrXI",   "chrXII",
    "chrXIII", "chrXIV",  "chrXV",   "chrXVI",
]

ROMAN_NUMERALS = [
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII",
    "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI",
]

KNOWN_CHROMOSOME_CONVENTIONS = {
    "numeric": {
        str(i): f"chr{roman}"
        for i, roman in enumerate(ROMAN_NUMERALS, start=1)
    },
    "roman_no_prefix": {
        roman: f"chr{roman}"
        for roman in ROMAN_NUMERALS
    },
    "chr_numeric": {
        f"chr{i}": f"chr{roman}"
        for i, roman in enumerate(ROMAN_NUMERALS, start=1)
    },
}

# Filtering thresholds for indel detection
MINIMUM_READ_COVERAGE = 10
MINIMUM_INDEL_FRACTION = 0.20   # at least 20% of reads must show the indel
MAXIMUM_BACKGROUND_FRACTION = 0.05  # noise ceiling in the other sample

OUTPUT_CSV_FILENAME = "candidate_indels_with_genes.csv"


# =============================================================================
# Argument parsing
# =============================================================================

argument_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
argument_parser.add_argument(
    "--experiments",
    type=str,
    default=None,
    help="Comma-separated experiment numbers to process (default: 1-9). Example: --experiments 1,7,9",
)
argument_parser.add_argument(
    "--clean",
    action="store_true",
    help="Delete cached indel pickle files before running.",
)
arguments = argument_parser.parse_args()

if arguments.experiments is not None:
    experiment_numbers = [int(x.strip()) for x in arguments.experiments.split(",")]
    print(f"  Running experiments: {experiment_numbers}")
else:
    experiment_numbers = DEFAULT_EXPERIMENT_NUMBERS

start_time = time.perf_counter()


# =============================================================================
# Preflight checks
# =============================================================================

print("=" * 60)
print("PREFLIGHT CHECKS")
print("=" * 60)

# -- Handle --clean flag --
if arguments.clean and INDEL_PICKLE_DIRECTORY.exists():
    import shutil
    shutil.rmtree(INDEL_PICKLE_DIRECTORY)
    print(f"  [CLEAN] Deleted cached indel pickles: {INDEL_PICKLE_DIRECTORY}")

# -- Check data directory --
if not DATA_DIRECTORY.is_dir():
    sys.exit(
        f"[ERROR] Data directory not found: {DATA_DIRECTORY}\n"
        f"  Create it and place sorted BAM files inside."
    )
print(f"  [OK] Data directory: {DATA_DIRECTORY}")

# -- Check BAM files --
missing_bam_files = []
for experiment_number in experiment_numbers:
    for template in [HIT_BAM_TEMPLATE, CTRL_BAM_TEMPLATE]:
        bam_path = DATA_DIRECTORY / template.format(experiment_number=experiment_number)
        if not bam_path.exists():
            missing_bam_files.append(bam_path.name)

if len(missing_bam_files) > 0:
    sys.exit(
        f"[ERROR] Missing BAM files in {DATA_DIRECTORY}:\n"
        + "\n".join(f"  - {filename}" for filename in missing_bam_files)
    )
print(f"  [OK] All {len(experiment_numbers) * 2} BAM files found")

# -- Validate BAM sort order --
sample_bam_path = DATA_DIRECTORY / HIT_BAM_TEMPLATE.format(
    experiment_number=experiment_numbers[0]
)
sample_bam_file = pysam.AlignmentFile(str(sample_bam_path), mode="rb")
bam_header = sample_bam_file.header.to_dict()

sort_order = bam_header.get("HD", {}).get("SO", "unknown")
if sort_order != "coordinate":
    print(
        f"  [WARNING] BAM sort order is '{sort_order}', expected 'coordinate'.\n"
        f"    Sort with: samtools sort -o sorted.bam input.bam"
    )
else:
    print(f"  [OK] BAM sort order: coordinate")

# -- Index BAM files if needed --
newly_indexed_count = 0
for experiment_number in experiment_numbers:
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

# -- Detect chromosome naming convention --
bam_contig_names = set(sample_bam_file.references)
sample_bam_file.close()

if set(REFERENCE_CHROMOSOMES).issubset(bam_contig_names):
    bam_to_reference_chromosome_map = {
        chromosome: chromosome for chromosome in REFERENCE_CHROMOSOMES
    }
    print(f"  [OK] BAM chromosome names match reference format (chrI-chrXVI)")
else:
    bam_to_reference_chromosome_map = None
    matched_convention_name = None

    for convention_name, convention_mapping in KNOWN_CHROMOSOME_CONVENTIONS.items():
        if set(convention_mapping.keys()).issubset(bam_contig_names):
            bam_to_reference_chromosome_map = convention_mapping
            matched_convention_name = convention_name
            break

    if bam_to_reference_chromosome_map is None:
        sample_contigs = sorted(bam_contig_names)[:20]
        sys.exit(
            f"[ERROR] Unrecognized chromosome names in BAM files.\n"
            f"  Found: {sample_contigs}\n"
            f"  Expected: chrI-chrXVI, 1-16, I-XVI, or chr1-chr16"
        )

    print(
        f"  [OK] BAM uses '{matched_convention_name}' chromosome names "
        f"- auto-mapped to chrI-chrXVI"
    )

reference_to_bam_chromosome_map = {
    reference_name: bam_name
    for bam_name, reference_name in bam_to_reference_chromosome_map.items()
}

# -- Verify gene coordinate reference --
if not GENE_COORDINATES_FILE.exists():
    sys.exit(
        f"[ERROR] Gene coordinate reference not found: {GENE_COORDINATES_FILE}\n"
        f"  Run 'uv run orc4r-screen_01_prepare-reference.py' first."
    )
if not GENE_COORDINATES_METADATA_FILE.exists():
    sys.exit(
        f"[ERROR] Reference metadata not found: {GENE_COORDINATES_METADATA_FILE}\n"
        f"  Run 'uv run orc4r-screen_01_prepare-reference.py' to regenerate."
    )

reference_metadata = json.loads(GENE_COORDINATES_METADATA_FILE.read_text())
expected_sha256 = reference_metadata["sha256"]
actual_sha256 = hashlib.sha256(GENE_COORDINATES_FILE.read_bytes()).hexdigest()

if actual_sha256 != expected_sha256:
    sys.exit(
        f"[ERROR] Gene coordinate file integrity check failed.\n"
        f"  Run 'uv run orc4r-screen_01_prepare-reference.py' to regenerate."
    )
print(
    f"  [OK] Gene coordinates verified "
    f"(downloaded: {reference_metadata['download_date']})"
)

# -- Create output directories --
OUTPUT_DIRECTORY.mkdir(exist_ok=True)
INDEL_PICKLE_DIRECTORY.mkdir(exist_ok=True)
print(f"  [OK] Output directory: {OUTPUT_DIRECTORY}")
print()


# =============================================================================
# Load gene coordinate reference
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
    for chromosome, chromosome_genes
    in gene_coordinates.groupby("chromosome")
}
print()


# =============================================================================
# Step 1: Extract indel data from BAM files
#
# For each position in the pileup, count reads showing insertions or deletions.
#
# pysam pileup read fields used:
#   pileup_read.indel > 0  : insertion of this many bases AFTER this position
#   pileup_read.indel < 0  : deletion of |indel| bases STARTING after this position
#   pileup_read.is_del     : this position falls WITHIN a deletion in this read
#
# We capture both insertion and deletion signals separately, plus the most
# common indel size for characterizing the variant.
# =============================================================================

print("=" * 60)
print("STEP 1: EXTRACT INDEL DATA FROM BAM FILES")
print("=" * 60)

try:
    for experiment_number in experiment_numbers:
        pickle_path = INDEL_PICKLE_DIRECTORY / f"df_indel_exp{experiment_number}.pickle"

        if pickle_path.exists():
            print(f"  [SKIP] Experiment {experiment_number}: cached pickle exists")
            continue

        print(f"  [PROCESSING] Experiment {experiment_number}")

        sample_dataframes = {}

        for suffix, bam_template in [("hit", HIT_BAM_TEMPLATE), ("ctrl", CTRL_BAM_TEMPLATE)]:
            bam_path = DATA_DIRECTORY / bam_template.format(experiment_number=experiment_number)

            chromosome_names = []
            base_positions = []
            read_coverages = []
            insertion_counts = []
            deletion_counts = []
            dominant_indel_sizes = []

            bam_file = pysam.AlignmentFile(str(bam_path), mode="rb")

            for reference_chromosome in REFERENCE_CHROMOSOMES:
                bam_chromosome = reference_to_bam_chromosome_map[reference_chromosome]

                if reference_chromosome == bam_chromosome:
                    print(f"    {suffix:<4} | {reference_chromosome}")
                else:
                    print(f"    {suffix:<4} | {reference_chromosome} (BAM: {bam_chromosome})")

                for pileup_column in bam_file.pileup(contig=bam_chromosome, truncate=True):

                    position_insertion_count = 0
                    position_deletion_count = 0
                    indel_sizes_at_position = []

                    for pileup_read in pileup_column.pileups:
                        # Insertion: extra bases inserted after this position
                        if pileup_read.indel > 0:
                            position_insertion_count += 1
                            indel_sizes_at_position.append(pileup_read.indel)

                        # Deletion start: bases deleted starting after this position
                        if pileup_read.indel < 0:
                            indel_sizes_at_position.append(pileup_read.indel)

                        # Position falls within an existing deletion in this read
                        if pileup_read.is_del:
                            position_deletion_count += 1

                    # Determine the most common indel size at this position
                    # Positive = insertion, negative = deletion
                    if len(indel_sizes_at_position) > 0:
                        size_counter = collections.Counter(indel_sizes_at_position)
                        dominant_size = size_counter.most_common(1)[0][0]
                    else:
                        dominant_size = 0

                    chromosome_names.append(reference_chromosome)
                    base_positions.append(pileup_column.pos)
                    read_coverages.append(pileup_column.n)
                    insertion_counts.append(position_insertion_count)
                    deletion_counts.append(position_deletion_count)
                    dominant_indel_sizes.append(dominant_size)

            bam_file.close()

            sample_dataframe = pandas.DataFrame({
                "chromosome": chromosome_names,
                "base_position": base_positions,
                f"coverage_{suffix}": read_coverages,
                f"insertion_count_{suffix}": insertion_counts,
                f"deletion_count_{suffix}": deletion_counts,
                f"dominant_indel_size_{suffix}": dominant_indel_sizes,
            })

            assert len(sample_dataframe) > 0, (
                f"No data extracted from {bam_path.name}."
            )
            print(f"    {suffix:<4} | {len(sample_dataframe):,} positions extracted")

            sample_dataframes[suffix] = sample_dataframe

        # -- Merge hit and ctrl --
        hit_count = len(sample_dataframes["hit"])
        ctrl_count = len(sample_dataframes["ctrl"])

        merged_dataframe = pandas.merge(
            left=sample_dataframes["hit"],
            right=sample_dataframes["ctrl"],
            on=["chromosome", "base_position"],
            how="inner",
        )

        merged_count = len(merged_dataframe)
        print(
            f"    Merged: {hit_count:,} (hit) + {ctrl_count:,} (ctrl) "
            f"-> {merged_count:,} shared positions"
        )

        larger_input_count = max(hit_count, ctrl_count)
        if merged_count < larger_input_count * 0.9:
            print(f"    [WARNING] Merge retained <90% of positions.")

        merged_dataframe.to_pickle(str(pickle_path))
        print(f"    Cached: {pickle_path.name}")

except KeyboardInterrupt:
    print("\n  [INTERRUPTED] Received Ctrl+C during BAM extraction.")
    print("  Any fully processed experiments are cached. Re-run to resume.")
    sys.exit(1)

print()


# =============================================================================
# Step 2: Filter for candidate indels
#
# A candidate indel position meets ALL of these criteria:
#   1) Coverage > MINIMUM_READ_COVERAGE in both samples
#   2) Indel fraction > MINIMUM_INDEL_FRACTION in one sample (the "signal")
#   3) Indel fraction < MAXIMUM_BACKGROUND_FRACTION in the other (the "background")
#
# We check both directions: indel in hit but not ctrl (suppressor gained an
# indel) and indel in ctrl but not hit (less likely but checked for completeness).
# =============================================================================

print("=" * 60)
print("STEP 2: FILTER FOR CANDIDATE INDELS")
print("=" * 60)
print(
    f"  Criteria: coverage > {MINIMUM_READ_COVERAGE}, "
    f"indel fraction > {MINIMUM_INDEL_FRACTION:.0%} in one sample, "
    f"< {MAXIMUM_BACKGROUND_FRACTION:.0%} in the other"
)

all_candidate_dataframes = []

for experiment_number in experiment_numbers:
    pickle_path = INDEL_PICKLE_DIRECTORY / f"df_indel_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    # Compute indel fractions (insertions + deletions combined)
    # Guard against division by zero at positions with zero coverage
    merged_dataframe["indel_fraction_hit"] = numpy.where(
        merged_dataframe["coverage_hit"] > 0,
        (merged_dataframe["insertion_count_hit"] + merged_dataframe["deletion_count_hit"])
        / merged_dataframe["coverage_hit"],
        0.0,
    )
    merged_dataframe["indel_fraction_ctrl"] = numpy.where(
        merged_dataframe["coverage_ctrl"] > 0,
        (merged_dataframe["insertion_count_ctrl"] + merged_dataframe["deletion_count_ctrl"])
        / merged_dataframe["coverage_ctrl"],
        0.0,
    )

    # Coverage filter: both samples must be well-sequenced
    passes_coverage = (
        (merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE)
        & (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
    )

    # Direction 1: indel present in hit (suppressor), absent in ctrl (wild-type)
    indel_in_hit = (
        passes_coverage
        & (merged_dataframe["indel_fraction_hit"] > MINIMUM_INDEL_FRACTION)
        & (merged_dataframe["indel_fraction_ctrl"] < MAXIMUM_BACKGROUND_FRACTION)
    )

    # Direction 2: indel present in ctrl, absent in hit (less expected but checked)
    indel_in_ctrl = (
        passes_coverage
        & (merged_dataframe["indel_fraction_ctrl"] > MINIMUM_INDEL_FRACTION)
        & (merged_dataframe["indel_fraction_hit"] < MAXIMUM_BACKGROUND_FRACTION)
    )

    is_candidate = indel_in_hit | indel_in_ctrl

    candidate_dataframe = merged_dataframe[is_candidate].copy()
    candidate_dataframe["experiment"] = f"Exp_{experiment_number}"

    # Label which sample carries the indel
    candidate_dataframe["indel_sample"] = numpy.where(
        indel_in_hit[is_candidate].values, "hit", "ctrl"
    )

    # Determine indel type based on which count is dominant
    candidate_dataframe["indel_type"] = "unknown"
    for row_index in candidate_dataframe.index:
        sample = candidate_dataframe.at[row_index, "indel_sample"]
        ins_count = candidate_dataframe.at[row_index, f"insertion_count_{sample}"]
        del_count = candidate_dataframe.at[row_index, f"deletion_count_{sample}"]
        if ins_count > del_count:
            candidate_dataframe.at[row_index, "indel_type"] = "insertion"
        elif del_count > ins_count:
            candidate_dataframe.at[row_index, "indel_type"] = "deletion"
        else:
            candidate_dataframe.at[row_index, "indel_type"] = "mixed"

    # Copy the dominant indel size from the signal sample
    candidate_dataframe["indel_size"] = numpy.where(
        candidate_dataframe["indel_sample"] == "hit",
        candidate_dataframe["dominant_indel_size_hit"],
        candidate_dataframe["dominant_indel_size_ctrl"],
    )

    candidate_count = len(candidate_dataframe)
    position_word = "position" if candidate_count == 1 else "positions"
    print(f"  Experiment {experiment_number}: {candidate_count:,} candidate {position_word}")

    if candidate_count > 0:
        hit_count = (candidate_dataframe["indel_sample"] == "hit").sum()
        ctrl_count = (candidate_dataframe["indel_sample"] == "ctrl").sum()
        ins_count = (candidate_dataframe["indel_type"] == "insertion").sum()
        del_count = (candidate_dataframe["indel_type"] == "deletion").sum()
        print(f"    In hit: {hit_count}, in ctrl: {ctrl_count} | Insertions: {ins_count}, Deletions: {del_count}")

    all_candidate_dataframes.append(candidate_dataframe)

all_candidates = pandas.concat(all_candidate_dataframes, ignore_index=True)
total_candidate_count = len(all_candidates)
print(f"  Total indel candidates across all experiments: {total_candidate_count:,}")

if total_candidate_count == 0:
    print("  [NOTE] No indel candidates found. Consider adjusting thresholds.")
    all_candidates.to_csv(str(OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME), index=False)
    elapsed_seconds = time.perf_counter() - start_time
    print(f"\nTotal runtime: {elapsed_seconds:.1f} seconds")
    sys.exit(0)

print()


# =============================================================================
# Step 3: Annotate candidate positions with gene names
#
# Same interval lookup as orc4r-screen_02_analyze-bam.py.
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

    chromosome_genes = gene_coordinates_by_chromosome.get(chromosome, pandas.DataFrame())

    if len(chromosome_genes) == 0:
        print(f" -> [WARNING] no genes loaded for {chromosome}")
        intergenic_count += 1
        continue

    overlapping_genes = chromosome_genes[
        (chromosome_genes["start"] <= base_position)
        & (chromosome_genes["end"] > base_position)
    ]

    if len(overlapping_genes) == 0:
        print(" -> (intergenic)")
        intergenic_count += 1
        continue

    first_gene = overlapping_genes.iloc[0]
    all_candidates.at[candidate_index, "gene_symbol"] = first_gene.get("gene_symbol", pandas.NA)
    all_candidates.at[candidate_index, "gene_name"] = first_gene.get("gene_name", pandas.NA)
    all_candidates.at[candidate_index, "gene_alias"] = first_gene.get("gene_alias", pandas.NA)
    all_candidates.at[candidate_index, "gene_feature_type"] = first_gene.get("feature_type", pandas.NA)
    all_candidates.at[candidate_index, "gene_description"] = first_gene.get("description", pandas.NA)

    display_name = (
        first_gene.get("gene_symbol")
        or first_gene.get("gene_name")
        or first_gene.get("gene_alias")
        or "unnamed"
    )
    print(f" -> {display_name}")

if intergenic_count > 0:
    print(f"  [NOTE] {intergenic_count} position(s) in intergenic regions")

print()


# =============================================================================
# Step 4: Save results
# =============================================================================

print("=" * 60)
print("STEP 4: SAVE RESULTS")
print("=" * 60)

output_column_order = [
    "experiment",
    "chromosome",
    "base_position",
    "indel_type",
    "indel_size",
    "indel_sample",
    "indel_fraction_hit",
    "indel_fraction_ctrl",
    "coverage_hit",
    "coverage_ctrl",
    "insertion_count_hit",
    "insertion_count_ctrl",
    "deletion_count_hit",
    "deletion_count_ctrl",
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

# Round fractions for readability
all_candidates_sorted["indel_fraction_hit"] = all_candidates_sorted["indel_fraction_hit"].round(3)
all_candidates_sorted["indel_fraction_ctrl"] = all_candidates_sorted["indel_fraction_ctrl"].round(3)

output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
all_candidates_sorted.to_csv(str(output_csv_path), index=False)
print(f"  Saved {total_candidate_count:,} candidates to: {output_csv_path}")
print()


# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("SUMMARY")
print("=" * 60)

# Per-experiment breakdown
for experiment_number in experiment_numbers:
    exp_label = f"Exp_{experiment_number}"
    exp_candidates = all_candidates_sorted[all_candidates_sorted["experiment"] == exp_label]

    if len(exp_candidates) == 0:
        print(f"  Experiment {experiment_number}: no indel candidates")
        continue

    print(f"\n  Experiment {experiment_number}: {len(exp_candidates)} candidate(s)")

    for _, row in exp_candidates.iterrows():
        gene_display = row["gene_symbol"] if pandas.notna(row["gene_symbol"]) else "(intergenic)"
        print(
            f"    {row['chromosome']}:{row['base_position']} "
            f"{row['indel_type']} (size={row['indel_size']:+d}) "
            f"in {row['indel_sample']} "
            f"[hit={row['indel_fraction_hit']:.1%}, ctrl={row['indel_fraction_ctrl']:.1%}] "
            f"-> {gene_display}"
        )

# Recurrent genes
annotated = all_candidates_sorted[all_candidates_sorted["gene_symbol"].notna()]
if len(annotated) > 0:
    gene_counts = (
        annotated.groupby("gene_symbol")["experiment"]
        .nunique()
        .sort_values(ascending=False)
    )
    recurrent = gene_counts[gene_counts > 1]
    if len(recurrent) > 0:
        print(f"\n  {len(recurrent)} gene(s) with indels in multiple experiments:")
        for gene_symbol, count in recurrent.items():
            experiments = ", ".join(sorted(
                annotated[annotated["gene_symbol"] == gene_symbol]["experiment"].unique()
            ))
            print(f"    {gene_symbol:<15} {count} experiments ({experiments})")

print()
elapsed_seconds = time.perf_counter() - start_time
print(f"Total runtime: {elapsed_seconds:.1f} seconds")
print("Done.")
