# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pysam",
#     "pandas",
#     "numpy",
#     "intermine",
# ]
# ///
"""
Identify potential causative mutations from yeast whole-genome sequencing data.

Compares suppressor ("hit") samples against wild-type ("ctrl") samples across
multiple experiments. Candidate mutations are identified where:
  1) More than 10 reads map to the base (both samples)
  2) More than 90% of reads agree on the nucleotide (both samples)
  3) The nucleotide differs between wild-type and suppressor

Candidate positions are then annotated with gene names via YeastMine (SGD).

Usage:
    uv run analyze_mutations.py
    uv run analyze_mutations.py --help

Expected directory structure:
    analyze_mutations.py
    data/
        Exp_1_hit_sorted.bam   (+ .bai index, auto-generated if missing)
        Exp_1_ctrl_sorted.bam
        ...through Exp_9...
    output/                    (created automatically)
        pickles/               (intermediate cached data)
        candidate_mutations_with_genes.csv
"""

import argparse
import collections
import pathlib
import sys

import numpy as np
import pandas as pd
import pysam
from intermine.webservice import Service


# =============================================================================
# Configuration - all thresholds, paths, and constants in one place
# =============================================================================

SCRIPT_DIRECTORY = pathlib.Path(__file__).resolve().parent
DATA_DIRECTORY = SCRIPT_DIRECTORY / "data"
OUTPUT_DIRECTORY = SCRIPT_DIRECTORY / "output"
PICKLE_DIRECTORY = OUTPUT_DIRECTORY / "pickles"

EXPERIMENT_NUMBERS = list(range(1, 10))  # experiments 1 through 9

HIT_BAM_TEMPLATE = "Exp_{experiment_number}_hit_sorted.bam"
CTRL_BAM_TEMPLATE = "Exp_{experiment_number}_ctrl_sorted.bam"

YEAST_CHROMOSOMES = [
    "chrI",   "chrII",   "chrIII",  "chrIV",
    "chrV",   "chrVI",   "chrVII",  "chrVIII",
    "chrIX",  "chrX",    "chrXI",   "chrXII",
    "chrXIII","chrXIV",  "chrXV",   "chrXVI",
]

# Filtering thresholds (must match manuscript methods section)
MINIMUM_READ_COVERAGE = 10     # reads mapped to the base
MINIMUM_CONFIDENCE_PERCENT = 90  # percent of reads agreeing on the nucleotide

YEASTMINE_SERVICE_URL = "https://yeastmine.yeastgenome.org/yeastmine/service"

OUTPUT_CSV_FILENAME = "candidate_mutations_with_genes.csv"

SAMPLE_SUFFIXES = ["hit", "ctrl"]


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
print(f"  [OK] Data directory exists: {DATA_DIRECTORY}")

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

# -- Create output directories --
OUTPUT_DIRECTORY.mkdir(exist_ok=True)
PICKLE_DIRECTORY.mkdir(exist_ok=True)
print(f"  [OK] Output directory ready: {OUTPUT_DIRECTORY}")
print()


# =============================================================================
# Step 1: Extract sequencing data from BAM files
#
# For each experiment, read the hit and ctrl BAM files, extract per-position
# consensus base calls and coverage, merge on shared positions, and cache
# the result as a pickle for fast re-runs.
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

    # We build one dataframe per sample (hit, ctrl), then merge them.
    # Storing in a dict keyed by suffix keeps the loop body identical.
    sample_dataframes = {}

    for suffix, bam_template in [("hit", HIT_BAM_TEMPLATE), ("ctrl", CTRL_BAM_TEMPLATE)]:
        bam_path = DATA_DIRECTORY / bam_template.format(experiment_number=experiment_number)

        # Accumulator lists - one entry per genomic position
        chromosome_names = []
        base_positions = []
        read_coverages = []
        consensus_bases = []
        consensus_confidences = []

        bam_file = pysam.AlignmentFile(str(bam_path), mode="rb")

        for chromosome in YEAST_CHROMOSOMES:
            print(f"    {suffix:<4} | {chromosome}")

            for pileup_column in bam_file.pileup(contig=chromosome, truncate=True):

                # Gather the base call from every read overlapping this position.
                # query_position is None when the read has a deletion here.
                base_calls_at_position = []
                for pileup_read in pileup_column.pileups:
                    if pileup_read.query_position is not None:
                        base_calls_at_position.append(
                            pileup_read.alignment.query_sequence[pileup_read.query_position]
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
                    confidence_percent = np.nan

                chromosome_names.append(chromosome)
                base_positions.append(pileup_column.pos)
                read_coverages.append(pileup_column.n)
                consensus_bases.append(most_common_base)
                consensus_confidences.append(confidence_percent)

        bam_file.close()

        sample_dataframe = pd.DataFrame({
            "chromosome": chromosome_names,
            "base_position": base_positions,
            f"base_value_{suffix}": consensus_bases,
            f"coverage_{suffix}": read_coverages,
            f"confidence_{suffix}": consensus_confidences,
        })

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

    merged_dataframe = pd.merge(
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

    # Sanity check: if merge drops a large fraction, the BAM files may be misaligned
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
    merged_dataframe = pd.read_pickle(str(pickle_path))

    passes_coverage_hit = merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE
    passes_coverage_ctrl = merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE
    passes_confidence_hit = merged_dataframe["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT
    passes_confidence_ctrl = merged_dataframe["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT
    bases_differ = merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]

    is_candidate_mutation = (
        passes_coverage_hit
        & passes_coverage_ctrl
        & passes_confidence_hit
        & passes_confidence_ctrl
        & bases_differ
    )

    candidate_dataframe = merged_dataframe[is_candidate_mutation].copy()
    candidate_dataframe["experiment"] = f"Exp_{experiment_number}"

    print(f"  Experiment {experiment_number}: {len(candidate_dataframe):,} candidate positions")
    all_candidate_dataframes.append(candidate_dataframe)

all_candidates = pd.concat(all_candidate_dataframes, ignore_index=True)
total_candidate_count = len(all_candidates)
print(f"  Total candidates across all experiments: {total_candidate_count:,}")

if total_candidate_count == 0:
    print("  [NOTE] No candidates found. Adjust thresholds or check data quality.")
    print("  Writing empty output file.")
    all_candidates.to_csv(str(OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME), index=False)
    sys.exit(0)

print()


# =============================================================================
# Step 3: Annotate candidate positions with gene names via YeastMine
#
# For each candidate position, query the SGD YeastMine service to find which
# gene (if any) overlaps that chromosomal coordinate. Positions falling in
# intergenic regions will have no gene annotation.
#
# Note: YeastMine is an external web service maintained by SGD. If this step
# fails due to connectivity issues, the candidate CSV from Step 2 is still
# valid - only gene names will be missing.
#
# When multiple genes overlap a position, the first match (sorted by
# primaryIdentifier) is kept.
# =============================================================================

print("=" * 60)
print("STEP 3: ANNOTATE WITH GENE NAMES (YeastMine)")
print("=" * 60)

all_candidates["gene_alias"] = pd.NA
all_candidates["gene_name"] = pd.NA
all_candidates["gene_symbol"] = pd.NA

try:
    yeastmine_service = Service(YEASTMINE_SERVICE_URL)
except Exception as connection_error:
    print(f"  [ERROR] Could not connect to YeastMine: {connection_error}")
    print(f"  Skipping gene annotation. Candidate positions are still saved.")
    yeastmine_service = None

annotation_failure_count = 0
intergenic_count = 0

if yeastmine_service is not None:
    for candidate_index in range(total_candidate_count):
        chromosome = all_candidates.at[candidate_index, "chromosome"]
        base_position = all_candidates.at[candidate_index, "base_position"]
        experiment = all_candidates.at[candidate_index, "experiment"]

        print(
            f"  [{candidate_index + 1}/{total_candidate_count}] "
            f"{experiment} | {chromosome}:{base_position}",
            end="",
        )

        try:
            query = yeastmine_service.new_query("Gene")
            query.add_view(
                "chromosome.primaryIdentifier",
                "primaryIdentifier",
                "secondaryIdentifier",
                "featureType",
                "symbol",
                "name",
                "sgdAlias",
                "organism.shortName",
                "qualifier",
                "chromosomeLocation.start",
                "chromosomeLocation.end",
                "chromosomeLocation.strand",
                "description",
            )
            query.add_sort_order("Gene.primaryIdentifier", "ASC")
            # Gene overlaps this base if: gene_start <= base_position < gene_end
            query.add_constraint("chromosome.primaryIdentifier", "=", chromosome)
            query.add_constraint("chromosomeLocation.start", "<=", str(base_position))
            query.add_constraint("chromosomeLocation.end", ">", str(base_position))

            gene_found = False
            for gene_row in query.rows():
                all_candidates.at[candidate_index, "gene_alias"] = gene_row["sgdAlias"]
                all_candidates.at[candidate_index, "gene_name"] = gene_row["name"]
                all_candidates.at[candidate_index, "gene_symbol"] = gene_row["symbol"]
                gene_found = True
                display_name = (
                    gene_row["symbol"] or gene_row["name"] or gene_row["sgdAlias"]
                )
                print(f"  {display_name}")
                break  # take first overlapping gene; see note above

            if not gene_found:
                print("  (intergenic)")
                intergenic_count += 1

        except Exception as query_error:
            print(f"  [ERROR] {query_error}")
            annotation_failure_count += 1

if annotation_failure_count > 0:
    print(f"  [WARNING] {annotation_failure_count} annotation lookups failed")
if intergenic_count > 0:
    print(f"  [NOTE] {intergenic_count} positions in intergenic regions")

print()


# =============================================================================
# Step 4: Save results
# =============================================================================

print("=" * 60)
print("STEP 4: SAVE RESULTS")
print("=" * 60)

output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
all_candidates.to_csv(str(output_csv_path), index=False)
print(f"  Saved {total_candidate_count:,} candidates to: {output_csv_path}")
print()
print("Done.")
