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

# Maximum gap between positions to consider them part of the same region/cluster
MAX_REGION_GAP_BP = 10

# Mismatch clustering: low-confidence differing positions
MAXIMUM_CLUSTER_CONFIDENCE = 80  # positions with confidence below this are "low quality mismatches"
MAX_CLUSTER_GAP_BP = 500         # maximum gap between mismatches to form a cluster
MIN_CLUSTER_SIZE = 3             # minimum positions in a reportable cluster

# Positional gap detection: flag unexplained gaps in merged position data
MIN_SUSPICIOUS_GAP_BP = 200      # gaps larger than this are reported
MIN_FLANKING_COVERAGE = 15       # positions flanking the gap must be well-covered

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
    help="Comma-separated experiment numbers to analyze (default: 1-9).",
)
arguments = argument_parser.parse_args()

if arguments.experiments is not None:
    experiment_numbers = [int(x.strip()) for x in arguments.experiments.split(",")]
else:
    experiment_numbers = DEFAULT_EXPERIMENT_NUMBERS

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
# Main analysis loop - process each experiment
# =============================================================================

all_candidate_regions = []

for experiment_number in experiment_numbers:
    print("=" * 60)
    print(f"EXPERIMENT {experiment_number}")
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
    #
    # Find positions where one sample has <25% of the other's coverage,
    # then cluster consecutive positions into regions.
    # =================================================================

    print(f"  --- Coverage Imbalance (ratio < {COVERAGE_IMBALANCE_RATIO}) ---")

    # Drop in hit (possible deletion/insertion in suppressor)
    drop_in_hit = (
        (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        & (merged_dataframe["coverage_hit"] < merged_dataframe["coverage_ctrl"] * COVERAGE_IMBALANCE_RATIO)
    )

    # Drop in ctrl (possible deletion/insertion in wild-type - less expected)
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

                # Get coverage stats for this region
                region_data = merged_dataframe[
                    (merged_dataframe["chromosome"] == chromosome)
                    & (merged_dataframe["base_position"] >= region_start)
                    & (merged_dataframe["base_position"] <= region_end)
                ]
                mean_cov_hit = region_data["coverage_hit"].mean()
                mean_cov_ctrl = region_data["coverage_ctrl"].mean()

                genes = find_overlapping_genes(chromosome, region_start, region_end)
                gene_display = ", ".join(genes) if genes else "(intergenic)"

                print(
                    f"      {chromosome}:{region_start}-{region_end} "
                    f"({region_length:,} bp, {position_count} pos) "
                    f"cov: hit={mean_cov_hit:.0f}, ctrl={mean_cov_ctrl:.0f} "
                    f"-> {gene_display}"
                )

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
    #
    # Find positions where bases differ between hit and ctrl but
    # confidence is low (reads split ~50/50). Cluster spatially close
    # positions - clusters indicate breakpoint regions.
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

                # Get confidence stats for this cluster
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

                print(
                    f"      {chromosome}:{cluster_start}-{cluster_end} "
                    f"({cluster_span:,} bp span, {cluster_count} mismatches) "
                    f"conf: hit={mean_conf_hit:.0f}%, ctrl={mean_conf_ctrl:.0f}% "
                    f"-> {gene_display}"
                )

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
    #
    # The SNV pipeline uses an inner merge, so positions present in one
    # sample but not the other are dropped. Large gaps in the merged
    # data that have well-covered flanking positions suggest a region
    # deleted (or unmappable) in one sample.
    #
    # Note: we cannot determine WHICH sample is missing reads in the
    # gap. The flanking coverage provides hints.
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

            # Check flanking coverage: positions just before and after the gap
            flank_before_cov_hit = coverages_hit[i - 1]
            flank_before_cov_ctrl = coverages_ctrl[i - 1]
            flank_after_cov_hit = coverages_hit[i]
            flank_after_cov_ctrl = coverages_ctrl[i]

            # Only report if flanking positions are well-covered in at least one sample
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

            # Determine which sample likely has the gap
            # If flanking hit coverage is lower, the gap is probably in hit
            avg_flank_hit = (flank_before_cov_hit + flank_after_cov_hit) / 2
            avg_flank_ctrl = (flank_before_cov_ctrl + flank_after_cov_ctrl) / 2

            if avg_flank_hit < avg_flank_ctrl * 0.5:
                likely_sample = "hit"
            elif avg_flank_ctrl < avg_flank_hit * 0.5:
                likely_sample = "ctrl"
            else:
                likely_sample = "unclear"

            print(
                f"      {chromosome}:{gap_start}-{gap_end} "
                f"({gap_size:,} bp gap) "
                f"flank cov: hit={avg_flank_hit:.0f}, ctrl={avg_flank_ctrl:.0f} "
                f"likely in: {likely_sample} "
                f"-> {gene_display}"
            )

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
# Save results
# =============================================================================

print("=" * 60)
print("SAVE RESULTS")
print("=" * 60)

if len(all_candidate_regions) == 0:
    print("  No structural variant candidates found across any experiment.")
    elapsed_seconds = time.perf_counter() - start_time
    print(f"\nTotal runtime: {elapsed_seconds:.1f} seconds")
    sys.exit(0)

results_dataframe = pandas.DataFrame(all_candidate_regions)

chromosome_sort_order = {c: i for i, c in enumerate(REFERENCE_CHROMOSOMES)}
results_dataframe = results_dataframe.sort_values(
    by=["experiment", "chromosome", "region_start"],
    key=lambda column: (
        column.map(chromosome_sort_order) if column.name == "chromosome" else column
    ),
    ignore_index=True,
)

output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
results_dataframe.to_csv(str(output_csv_path), index=False)
print(f"  Saved {len(results_dataframe)} candidate regions to: {output_csv_path}")
print()


# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("SUMMARY")
print("=" * 60)

for method in ["coverage_imbalance", "mismatch_cluster", "positional_gap"]:
    method_results = results_dataframe[results_dataframe["detection_method"] == method]
    if len(method_results) > 0:
        print(f"\n  {method.replace('_', ' ').title()}: {len(method_results)} region(s)")
        for _, row in method_results.iterrows():
            print(
                f"    {row['experiment']} | {row['chromosome']}:{row['region_start']}-{row['region_end']} "
                f"({row['region_length_bp']:,} bp) "
                f"in {row['affected_sample']} -> {row['genes']}"
            )

# Cross-reference: genes appearing across multiple experiments or detection methods
all_genes_mentioned = set()
for genes_str in results_dataframe["genes"]:
    if genes_str != "(intergenic)":
        for gene in genes_str.split(", "):
            all_genes_mentioned.add(gene)

if len(all_genes_mentioned) > 0:
    print(f"\n  Genes in candidate regions: {', '.join(sorted(all_genes_mentioned))}")

print()
elapsed_seconds = time.perf_counter() - start_time
print(f"Total runtime: {elapsed_seconds:.1f} seconds")
print("Done.")
