# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "numpy",
# ]
# ///
"""
Detect structural variants using coverage imbalance and mismatch clustering.

Uses the existing SNV pickles from orc4r-screen_02_analyze-bam.py - no BAM
reprocessing required. Targets experiments where neither SNV detection (02) nor
small indel detection (04) found candidates.

Detection strategies:
  1) Coverage imbalance: regions where one sample has dramatically lower
     coverage, indicating large deletions or insertions that prevent read
     mapping.
  2) Low-confidence mismatch clustering: spatial clusters of positions where
     bases differ at low confidence (~50-60%), the hallmark of misaligned
     reads flanking a structural breakpoint.
  3) Positional gap detection: unexplained gaps in the merged position data
     where one sample may be missing reads entirely.

After detection, a cross-experiment filter removes regions that also appear in
experiments with known SNV hits (background experiments). Regions shared across
many experiments are artifacts (repetitive regions, GC-biased coverage dips).
Only regions unique to the target experiments are reported.

Usage:
    uv run orc4r-screen_05_structural-variants.py
    uv run orc4r-screen_05_structural-variants.py --experiments 1,7,9
    uv run orc4r-screen_05_structural-variants.py --help

Prerequisites:
    uv run orc4r-screen_02_analyze-bam.py   (generates SNV pickles)
    uv run orc4r-screen_01_prepare-reference.py   (gene annotations)
"""

import argparse
import hashlib
import json
import pathlib
import sys
import time

import numpy
import pandas


# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIRECTORY = pathlib.Path(__file__).resolve().parent

PICKLE_DIRECTORY = SCRIPT_DIRECTORY / "output" / "pickles"
OUTPUT_DIRECTORY = SCRIPT_DIRECTORY / "output"
REFERENCE_DIRECTORY = SCRIPT_DIRECTORY / "reference"

GENE_COORDINATES_FILE = REFERENCE_DIRECTORY / "gene_coordinates.tsv"
GENE_COORDINATES_METADATA_FILE = REFERENCE_DIRECTORY / "gene_coordinates.meta.json"

DEFAULT_EXPERIMENT_NUMBERS = list(range(1, 10))

REFERENCE_CHROMOSOMES = [
    "chrI",    "chrII",   "chrIII",  "chrIV",
    "chrV",    "chrVI",   "chrVII",  "chrVIII",
    "chrIX",   "chrX",    "chrXI",   "chrXII",
    "chrXIII", "chrXIV",  "chrXV",   "chrXVI",
]

# Must match thresholds in orc4r-screen_02_analyze-bam.py
MINIMUM_READ_COVERAGE = 10

# Coverage imbalance: flag positions where one sample has <25% of the other's coverage
COVERAGE_IMBALANCE_RATIO = 0.25

# Minimum contiguous region length (bp) to report as a coverage drop
MINIMUM_REGION_LENGTH_BP = 50

# Maximum gap between positions to consider them part of the same region
MAX_REGION_GAP_BP = 10

# Mismatch clustering: low-confidence differing positions
MAXIMUM_CLUSTER_CONFIDENCE = 80  # positions with confidence below this are "low quality mismatches"
MAX_CLUSTER_GAP_BP = 500         # maximum gap between mismatches to form a cluster
MIN_CLUSTER_SIZE = 3             # minimum positions in a reportable cluster

# Positional gap detection: flag unexplained gaps in merged position data
MIN_SUSPICIOUS_GAP_BP = 200      # gaps larger than this are reported
MIN_FLANKING_COVERAGE = 15       # positions flanking the gap must be well-covered

# Cross-experiment filtering
# Experiments with known SNV hits serve as negative controls. Any region also
# found in these experiments is a shared artifact, not a suppressor-specific variant.
BACKGROUND_EXPERIMENT_NUMBERS = [2, 3, 4, 5, 6, 8]
TARGET_EXPERIMENT_NUMBERS = [1, 7, 9]
REGION_PROXIMITY_BP = 1000  # regions within this distance on same chromosome are considered shared

OUTPUT_CSV_FILENAME = "candidate_structural_variants.csv"


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
    help=(
        "Comma-separated experiment numbers to analyze (default: all 1-9). "
        "For effective background filtering, include both target and background experiments."
    ),
)
arguments = argument_parser.parse_args()

if arguments.experiments is not None:
    experiment_numbers = [int(x.strip()) for x in arguments.experiments.split(",")]
else:
    experiment_numbers = DEFAULT_EXPERIMENT_NUMBERS

# Determine which of the requested experiments are targets vs background
active_target_experiments = [n for n in experiment_numbers if n in TARGET_EXPERIMENT_NUMBERS]
active_background_experiments = [n for n in experiment_numbers if n in BACKGROUND_EXPERIMENT_NUMBERS]

start_time = time.perf_counter()


# =============================================================================
# Preflight checks
# =============================================================================

print("=" * 60)
print("PREFLIGHT CHECKS")
print("=" * 60)

if not PICKLE_DIRECTORY.is_dir():
    sys.exit(
        f"[ERROR] Pickle directory not found: {PICKLE_DIRECTORY}\n"
        f"  Run 'uv run orc4r-screen_02_analyze-bam.py' first."
    )

missing_pickles = []
for experiment_number in experiment_numbers:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    if not pickle_path.exists():
        missing_pickles.append(pickle_path.name)

if len(missing_pickles) > 0:
    sys.exit(
        f"[ERROR] Missing pickle files:\n"
        + "\n".join(f"  - {f}" for f in missing_pickles)
        + f"\n  Run 'uv run orc4r-screen_02_analyze-bam.py' first."
    )
print(f"  [OK] All {len(experiment_numbers)} SNV pickle files found")

if len(active_background_experiments) == 0:
    print(
        f"  [WARNING] No background experiments included. Cross-experiment filtering\n"
        f"    will be skipped. For best results, run all 9 experiments."
    )
else:
    print(f"  [OK] Background experiments (known SNVs): {active_background_experiments}")
    print(f"  [OK] Target experiments (no SNVs): {active_target_experiments}")

# -- Verify gene coordinates --
if not GENE_COORDINATES_FILE.exists() or not GENE_COORDINATES_METADATA_FILE.exists():
    sys.exit(
        f"[ERROR] Gene coordinate reference not found.\n"
        f"  Run 'uv run orc4r-screen_01_prepare-reference.py' first."
    )

reference_metadata = json.loads(GENE_COORDINATES_METADATA_FILE.read_text())
expected_sha256 = reference_metadata["sha256"]
actual_sha256 = hashlib.sha256(GENE_COORDINATES_FILE.read_bytes()).hexdigest()

if actual_sha256 != expected_sha256:
    sys.exit(
        f"[ERROR] Gene coordinate file integrity check failed.\n"
        f"  Run 'uv run orc4r-screen_01_prepare-reference.py' to regenerate."
    )
print(f"  [OK] Gene coordinates verified")
print()


# =============================================================================
# Load gene coordinate reference
# =============================================================================

gene_coordinates = pandas.read_csv(
    GENE_COORDINATES_FILE,
    sep="\t",
    dtype={"chromosome": str, "start": int, "end": int},
)

print(f"  Loaded {len(gene_coordinates):,} gene annotations")

gene_coordinates_by_chromosome = {
    chromosome: chromosome_genes.reset_index(drop=True)
    for chromosome, chromosome_genes
    in gene_coordinates.groupby("chromosome")
}
print()


# =============================================================================
# Helper: find genes overlapping a genomic region
# =============================================================================

def find_overlapping_genes(chromosome, region_start, region_end):
    """Return list of gene symbols overlapping the given region."""
    chromosome_genes = gene_coordinates_by_chromosome.get(chromosome, pandas.DataFrame())
    if len(chromosome_genes) == 0:
        return []
    overlapping = chromosome_genes[
        (chromosome_genes["start"] < region_end)
        & (chromosome_genes["end"] > region_start)
    ]
    symbols = []
    for _, gene in overlapping.iterrows():
        name = gene.get("gene_symbol") or gene.get("gene_name") or gene.get("gene_alias")
        if pandas.notna(name):
            symbols.append(name)
    return symbols


# =============================================================================
# Helper: cluster sorted positions into contiguous regions
# =============================================================================

def cluster_positions(positions, max_gap):
    """
    Group a sorted list of integer positions into clusters.
    Consecutive positions within max_gap of each other are merged.
    Returns list of (start, end, count) tuples.
    """
    if len(positions) == 0:
        return []

    clusters = []
    cluster_start = positions[0]
    cluster_end = positions[0]
    cluster_count = 1

    for i in range(1, len(positions)):
        if positions[i] - cluster_end <= max_gap:
            cluster_end = positions[i]
            cluster_count += 1
        else:
            clusters.append((cluster_start, cluster_end, cluster_count))
            cluster_start = positions[i]
            cluster_end = positions[i]
            cluster_count = 1

    clusters.append((cluster_start, cluster_end, cluster_count))
    return clusters


# =============================================================================
# Helper: check if a region overlaps any region in a background set
# =============================================================================

def overlaps_any_background_region(chromosome, region_start, region_end, background_regions, proximity):
    """
    Check if a region on a given chromosome overlaps (within proximity) any
    region in the background set. Returns True if shared, False if unique.
    """
    same_chromosome = background_regions[background_regions["chromosome"] == chromosome]
    for _, background_row in same_chromosome.iterrows():
        # Two regions overlap within proximity if:
        # region_start - proximity <= bg_end AND bg_start - proximity <= region_end
        if (region_start - proximity <= background_row["region_end"]
                and background_row["region_start"] - proximity <= region_end):
            return True
    return False


# =============================================================================
# Main analysis loop - process each experiment
# =============================================================================

all_candidate_regions = []

for experiment_number in experiment_numbers:
    print("=" * 60)
    print(f"EXPERIMENT {experiment_number}")
    experiment_label = "TARGET" if experiment_number in TARGET_EXPERIMENT_NUMBERS else "BACKGROUND"
    print(f"  ({experiment_label})")
    print("=" * 60)

    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    median_coverage_hit = merged_dataframe["coverage_hit"].median()
    median_coverage_ctrl = merged_dataframe["coverage_ctrl"].median()
    print(
        f"  Median coverage: hit={median_coverage_hit:.0f}x, "
        f"ctrl={median_coverage_ctrl:.0f}x"
    )
    print()


    # =================================================================
    # Section 1: Coverage imbalance
    # =================================================================

    print(f"  --- Coverage Imbalance (ratio < {COVERAGE_IMBALANCE_RATIO}) ---")

    drop_in_hit = (
        (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        & (merged_dataframe["coverage_hit"] < merged_dataframe["coverage_ctrl"] * COVERAGE_IMBALANCE_RATIO)
    )

    drop_in_ctrl = (
        (merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE)
        & (merged_dataframe["coverage_ctrl"] < merged_dataframe["coverage_hit"] * COVERAGE_IMBALANCE_RATIO)
    )

    for direction_label, direction_mask in [("hit", drop_in_hit), ("ctrl", drop_in_ctrl)]:
        flagged_positions = merged_dataframe[direction_mask].copy()

        if len(flagged_positions) == 0:
            print(f"    Drop in {direction_label}: 0 positions flagged")
            continue

        print(f"    Drop in {direction_label}: {len(flagged_positions):,} positions flagged")

        for chromosome in REFERENCE_CHROMOSOMES:
            chrom_positions = flagged_positions[
                flagged_positions["chromosome"] == chromosome
            ]["base_position"].sort_values().tolist()

            if len(chrom_positions) == 0:
                continue

            regions = cluster_positions(chrom_positions, max_gap=MAX_REGION_GAP_BP)

            for region_start, region_end, position_count in regions:
                region_length = region_end - region_start + 1

                if region_length < MINIMUM_REGION_LENGTH_BP:
                    continue

                region_data = merged_dataframe[
                    (merged_dataframe["chromosome"] == chromosome)
                    & (merged_dataframe["base_position"] >= region_start)
                    & (merged_dataframe["base_position"] <= region_end)
                ]
                mean_cov_hit = region_data["coverage_hit"].mean()
                mean_cov_ctrl = region_data["coverage_ctrl"].mean()

                genes = find_overlapping_genes(chromosome, region_start, region_end)
                gene_display = ", ".join(genes) if genes else "(intergenic)"

                all_candidate_regions.append({
                    "experiment": f"Exp_{experiment_number}",
                    "chromosome": chromosome,
                    "region_start": region_start,
                    "region_end": region_end,
                    "region_length_bp": region_length,
                    "position_count": position_count,
                    "detection_method": "coverage_imbalance",
                    "affected_sample": direction_label,
                    "mean_coverage_hit": round(mean_cov_hit, 1),
                    "mean_coverage_ctrl": round(mean_cov_ctrl, 1),
                    "genes": gene_display,
                })

    print()


    # =================================================================
    # Section 2: Low-confidence mismatch clustering
    # =================================================================

    print(f"  --- Low-Confidence Mismatch Clusters (conf < {MAXIMUM_CLUSTER_CONFIDENCE}%) ---")

    bases_differ = merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]
    low_confidence = (
        (merged_dataframe["confidence_hit"] < MAXIMUM_CLUSTER_CONFIDENCE)
        | (merged_dataframe["confidence_ctrl"] < MAXIMUM_CLUSTER_CONFIDENCE)
    )
    has_coverage = (
        (merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE)
        & (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
    )

    mismatch_positions = merged_dataframe[bases_differ & low_confidence & has_coverage].copy()

    if len(mismatch_positions) == 0:
        print(f"    0 low-confidence mismatch positions found")
    else:
        print(f"    {len(mismatch_positions):,} low-confidence mismatch positions found")

        for chromosome in REFERENCE_CHROMOSOMES:
            chrom_mismatches = mismatch_positions[
                mismatch_positions["chromosome"] == chromosome
            ].sort_values("base_position")

            if len(chrom_mismatches) == 0:
                continue

            positions = chrom_mismatches["base_position"].tolist()
            clusters = cluster_positions(positions, max_gap=MAX_CLUSTER_GAP_BP)

            for cluster_start, cluster_end, cluster_count in clusters:
                if cluster_count < MIN_CLUSTER_SIZE:
                    continue

                cluster_span = cluster_end - cluster_start + 1

                cluster_data = chrom_mismatches[
                    (chrom_mismatches["base_position"] >= cluster_start)
                    & (chrom_mismatches["base_position"] <= cluster_end)
                ]
                mean_conf_hit = cluster_data["confidence_hit"].mean()
                mean_conf_ctrl = cluster_data["confidence_ctrl"].mean()
                mean_cov_hit = cluster_data["coverage_hit"].mean()
                mean_cov_ctrl = cluster_data["coverage_ctrl"].mean()

                genes = find_overlapping_genes(chromosome, cluster_start, cluster_end)
                gene_display = ", ".join(genes) if genes else "(intergenic)"

                all_candidate_regions.append({
                    "experiment": f"Exp_{experiment_number}",
                    "chromosome": chromosome,
                    "region_start": cluster_start,
                    "region_end": cluster_end,
                    "region_length_bp": cluster_span,
                    "position_count": cluster_count,
                    "detection_method": "mismatch_cluster",
                    "affected_sample": "hit",
                    "mean_coverage_hit": round(mean_cov_hit, 1),
                    "mean_coverage_ctrl": round(mean_cov_ctrl, 1),
                    "genes": gene_display,
                })

    print()


    # =================================================================
    # Section 3: Positional gap detection
    # =================================================================

    print(f"  --- Positional Gaps (> {MIN_SUSPICIOUS_GAP_BP} bp) ---")

    gap_count = 0

    for chromosome in REFERENCE_CHROMOSOMES:
        chrom_data = merged_dataframe[
            merged_dataframe["chromosome"] == chromosome
        ].sort_values("base_position")

        if len(chrom_data) < 2:
            continue

        positions = chrom_data["base_position"].values
        coverages_hit = chrom_data["coverage_hit"].values
        coverages_ctrl = chrom_data["coverage_ctrl"].values

        for i in range(1, len(positions)):
            gap_size = positions[i] - positions[i - 1]

            if gap_size <= MIN_SUSPICIOUS_GAP_BP:
                continue

            flank_before_cov_hit = coverages_hit[i - 1]
            flank_before_cov_ctrl = coverages_ctrl[i - 1]
            flank_after_cov_hit = coverages_hit[i]
            flank_after_cov_ctrl = coverages_ctrl[i]

            before_well_covered = (
                flank_before_cov_hit > MIN_FLANKING_COVERAGE
                or flank_before_cov_ctrl > MIN_FLANKING_COVERAGE
            )
            after_well_covered = (
                flank_after_cov_hit > MIN_FLANKING_COVERAGE
                or flank_after_cov_ctrl > MIN_FLANKING_COVERAGE
            )

            if not (before_well_covered and after_well_covered):
                continue

            gap_start = int(positions[i - 1])
            gap_end = int(positions[i])

            genes = find_overlapping_genes(chromosome, gap_start, gap_end)
            gene_display = ", ".join(genes) if genes else "(intergenic)"

            avg_flank_hit = (flank_before_cov_hit + flank_after_cov_hit) / 2
            avg_flank_ctrl = (flank_before_cov_ctrl + flank_after_cov_ctrl) / 2

            if avg_flank_hit < avg_flank_ctrl * 0.5:
                likely_sample = "hit"
            elif avg_flank_ctrl < avg_flank_hit * 0.5:
                likely_sample = "ctrl"
            else:
                likely_sample = "unclear"

            gap_count += 1

            all_candidate_regions.append({
                "experiment": f"Exp_{experiment_number}",
                "chromosome": chromosome,
                "region_start": gap_start,
                "region_end": gap_end,
                "region_length_bp": gap_size,
                "position_count": 0,
                "detection_method": "positional_gap",
                "affected_sample": likely_sample,
                "mean_coverage_hit": round(avg_flank_hit, 1),
                "mean_coverage_ctrl": round(avg_flank_ctrl, 1),
                "genes": gene_display,
            })

    if gap_count == 0:
        print(f"    No suspicious gaps found")

    print()


# =============================================================================
# Cross-experiment filtering
#
# Regions found in background experiments (those with known SNV hits) are
# shared genomic artifacts - repetitive regions, GC-biased coverage dips,
# reference assembly gaps, etc. These are NOT suppressor-specific variants.
#
# For each candidate region in a target experiment, check whether a similar
# region (same chromosome, overlapping within REGION_PROXIMITY_BP) exists in
# ANY background experiment. If so, discard it.
# =============================================================================

print("=" * 60)
print("CROSS-EXPERIMENT FILTERING")
print("=" * 60)

results_dataframe = pandas.DataFrame(all_candidate_regions)

if len(results_dataframe) == 0:
    print("  No candidate regions to filter.")
    elapsed_seconds = time.perf_counter() - start_time
    print(f"\nTotal runtime: {elapsed_seconds:.1f} seconds")
    sys.exit(0)

target_labels = [f"Exp_{n}" for n in active_target_experiments]
background_labels = [f"Exp_{n}" for n in active_background_experiments]

target_regions = results_dataframe[results_dataframe["experiment"].isin(target_labels)].copy()
background_regions = results_dataframe[results_dataframe["experiment"].isin(background_labels)].copy()

total_raw_target = len(target_regions)
print(f"  Target experiments: {active_target_experiments}")
print(f"  Background experiments: {active_background_experiments}")
print(f"  Raw candidate regions in target experiments: {total_raw_target:,}")
print(f"  Regions in background experiments: {len(background_regions):,}")
print(f"  Proximity threshold: {REGION_PROXIMITY_BP:,} bp")

if len(active_background_experiments) == 0 or len(background_regions) == 0:
    print("  [SKIP] No background regions available for filtering.")
    filtered_regions = target_regions.copy()
else:
    is_shared = []
    for _, target_row in target_regions.iterrows():
        shared = overlaps_any_background_region(
            chromosome=target_row["chromosome"],
            region_start=target_row["region_start"],
            region_end=target_row["region_end"],
            background_regions=background_regions,
            proximity=REGION_PROXIMITY_BP,
        )
        is_shared.append(shared)

    shared_count = sum(is_shared)
    unique_count = total_raw_target - shared_count

    print(f"\n  Shared with background (removed): {shared_count:,}")
    print(f"  Unique to target experiments: {unique_count:,}")

    filtered_regions = target_regions[~pandas.Series(is_shared, index=target_regions.index)].copy()

# Per-method breakdown of what survived
if len(filtered_regions) > 0:
    print("\n  Surviving regions by detection method:")
    for method in ["coverage_imbalance", "mismatch_cluster", "positional_gap"]:
        method_count = (filtered_regions["detection_method"] == method).sum()
        if method_count > 0:
            print(f"    {method.replace('_', ' ').title()}: {method_count}")

print()


# =============================================================================
# Save results
# =============================================================================

print("=" * 60)
print("SAVE RESULTS")
print("=" * 60)

if len(filtered_regions) == 0:
    print("  No unique structural variant candidates after filtering.")
    output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
    pandas.DataFrame(columns=results_dataframe.columns).to_csv(
        str(output_csv_path), index=False
    )
    print(f"  Saved empty results to: {output_csv_path}")
    elapsed_seconds = time.perf_counter() - start_time
    print(f"\nTotal runtime: {elapsed_seconds:.1f} seconds")
    sys.exit(0)

chromosome_sort_order = {c: i for i, c in enumerate(REFERENCE_CHROMOSOMES)}

filtered_regions_sorted = filtered_regions.sort_values(
    by=["experiment", "chromosome", "region_start"],
    key=lambda column: (
        column.map(chromosome_sort_order) if column.name == "chromosome" else column
    ),
    ignore_index=True,
)

output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
filtered_regions_sorted.to_csv(str(output_csv_path), index=False)
print(f"  Saved {len(filtered_regions_sorted)} unique candidate regions to: {output_csv_path}")
print()


# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("SUMMARY: UNIQUE STRUCTURAL VARIANT CANDIDATES")
print("=" * 60)

for experiment_number in active_target_experiments:
    exp_label = f"Exp_{experiment_number}"
    exp_regions = filtered_regions_sorted[filtered_regions_sorted["experiment"] == exp_label]

    if len(exp_regions) == 0:
        print(f"\n  Experiment {experiment_number}: no unique candidates")
        continue

    print(f"\n  Experiment {experiment_number}: {len(exp_regions)} unique candidate(s)")

    for _, row in exp_regions.iterrows():
        print(
            f"    [{row['detection_method'].replace('_', ' ')}] "
            f"{row['chromosome']}:{row['region_start']}-{row['region_end']} "
            f"({row['region_length_bp']:,} bp) "
            f"in {row['affected_sample']} "
            f"cov: hit={row['mean_coverage_hit']:.0f}, ctrl={row['mean_coverage_ctrl']:.0f} "
            f"-> {row['genes']}"
        )

# List unique genes
all_genes = set()
for genes_str in filtered_regions_sorted["genes"]:
    if genes_str != "(intergenic)":
        for gene in genes_str.split(", "):
            all_genes.add(gene)

if len(all_genes) > 0:
    print(f"\n  Genes in unique candidate regions: {', '.join(sorted(all_genes))}")

print()
elapsed_seconds = time.perf_counter() - start_time
print(f"Total runtime: {elapsed_seconds:.1f} seconds")
print("Done.")
