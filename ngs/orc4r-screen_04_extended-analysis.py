# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pysam",
#     "pandas",
#     "numpy",
# ]
# ///
"""
Extended analysis: indel detection and structural variant search.

Documents the hypotheses tested after the SNV pipeline (02) and diagnostic (03)
identified experiments with zero SNV candidates. Both analyses were ultimately
ruled out - the suppressor mechanism was chromosome X disomy (found by 03).

Part A - Indel Detection:
  Re-processes BAM files to capture per-position insertion/deletion counts.
  Filters for positions with significant indel signal in one sample absent
  in the other. Result: zero indel candidates found.

Part B - Structural Variant Analysis:
  Uses existing SNV pickles to detect coverage imbalances, low-confidence
  mismatch clusters, and positional gaps. Applies cross-experiment background
  filtering. Result: no experiment-specific structural variants after filtering.

Usage:
    uv run orc4r-screen_04_extended-analysis.py                    (run both)
    uv run orc4r-screen_04_extended-analysis.py --mode indel       (Part A only)
    uv run orc4r-screen_04_extended-analysis.py --mode structural  (Part B only)
    uv run orc4r-screen_04_extended-analysis.py --experiments 1,7,9
    uv run orc4r-screen_04_extended-analysis.py --clean
    uv run orc4r-screen_04_extended-analysis.py --help

Prerequisites:
    uv run orc4r-screen_01_prepare-reference.py   (gene annotations)
    uv run orc4r-screen_02_analyze-bam.py          (SNV pickles, needed for Part B)
    BAM files in data/ directory                   (needed for Part A)
"""

import argparse
import collections
import hashlib
import json
import pathlib
import shutil
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
SNV_PICKLE_DIRECTORY = OUTPUT_DIRECTORY / "pickles"
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

# -- Indel detection thresholds (Part A) --
MINIMUM_READ_COVERAGE = 10
MINIMUM_INDEL_FRACTION = 0.20
MAXIMUM_BACKGROUND_FRACTION = 0.05

# -- Structural variant thresholds (Part B) --
COVERAGE_IMBALANCE_RATIO = 0.25
MINIMUM_REGION_LENGTH_BP = 200
MAX_REGION_GAP_BP = 10
MAXIMUM_CLUSTER_CONFIDENCE = 80
MAX_CLUSTER_GAP_BP = 500
MIN_CLUSTER_SIZE = 3
MIN_SUSPICIOUS_GAP_BP = 200
MIN_FLANKING_COVERAGE = 15

# -- Cross-experiment filtering (Part B) --
BACKGROUND_EXPERIMENT_NUMBERS = [2, 3, 4, 5, 6, 8, 9]
TARGET_EXPERIMENT_NUMBERS = [1, 7]
REGION_PROXIMITY_BP = 1000

INDEL_CSV_FILENAME = "candidate_indels_with_genes.csv"
STRUCTURAL_CSV_FILENAME = "candidate_structural_variants.csv"


# =============================================================================
# Argument parsing
# =============================================================================

argument_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
argument_parser.add_argument(
    "--mode",
    choices=["indel", "structural", "all"],
    default="all",
    help="Which analysis to run: 'indel' (Part A), 'structural' (Part B), or 'all' (default).",
)
argument_parser.add_argument(
    "--experiments",
    type=str,
    default=None,
    help="Comma-separated experiment numbers (default: 1-9).",
)
argument_parser.add_argument(
    "--clean",
    action="store_true",
    help="Delete cached indel pickle files before running.",
)
arguments = argument_parser.parse_args()

if arguments.experiments is not None:
    experiment_numbers = [int(x.strip()) for x in arguments.experiments.split(",")]
else:
    experiment_numbers = DEFAULT_EXPERIMENT_NUMBERS

run_indels = arguments.mode in ("indel", "all")
run_structural = arguments.mode in ("structural", "all")

start_time = time.perf_counter()


# =============================================================================
# Preflight checks
# =============================================================================

print("=" * 80)
print("PREFLIGHT CHECKS")
print("=" * 80)
print(f"  Mode: {arguments.mode}")

# -- Handle --clean --
if arguments.clean and INDEL_PICKLE_DIRECTORY.exists():
    shutil.rmtree(INDEL_PICKLE_DIRECTORY)
    print(f"  [CLEAN] Deleted indel pickles: {INDEL_PICKLE_DIRECTORY}")

# -- Check BAM files (Part A only) --
reference_to_bam_chromosome_map = None

if run_indels:
    if not DATA_DIRECTORY.is_dir():
        sys.exit(f"[ERROR] Data directory not found: {DATA_DIRECTORY}")

    missing_bam_files = []
    for experiment_number in experiment_numbers:
        for template in [HIT_BAM_TEMPLATE, CTRL_BAM_TEMPLATE]:
            bam_path = DATA_DIRECTORY / template.format(experiment_number=experiment_number)
            if not bam_path.exists():
                missing_bam_files.append(bam_path.name)

    if len(missing_bam_files) > 0:
        sys.exit(
            f"[ERROR] Missing BAM files:\n"
            + "\n".join(f"  - {f}" for f in missing_bam_files)
        )
    print(f"  [OK] All {len(experiment_numbers) * 2} BAM files found")

    # Index BAMs if needed
    newly_indexed = 0
    for experiment_number in experiment_numbers:
        for template in [HIT_BAM_TEMPLATE, CTRL_BAM_TEMPLATE]:
            bam_path = DATA_DIRECTORY / template.format(experiment_number=experiment_number)
            bai_path = pathlib.Path(str(bam_path) + ".bai")
            if not bai_path.exists():
                print(f"  [INDEXING] {bam_path.name} ...")
                pysam.index(str(bam_path))
                assert bai_path.exists(), f"Indexing failed for {bam_path.name}"
                newly_indexed += 1

    if newly_indexed > 0:
        print(f"  [OK] Indexed {newly_indexed} BAM files")
    else:
        print(f"  [OK] All BAM index files present")

    # Detect chromosome naming convention
    sample_bam_path = DATA_DIRECTORY / HIT_BAM_TEMPLATE.format(
        experiment_number=experiment_numbers[0]
    )
    sample_bam_file = pysam.AlignmentFile(str(sample_bam_path), mode="rb")
    bam_contig_names = set(sample_bam_file.references)
    sample_bam_file.close()

    if set(REFERENCE_CHROMOSOMES).issubset(bam_contig_names):
        bam_to_reference_map = {c: c for c in REFERENCE_CHROMOSOMES}
        print(f"  [OK] BAM chromosome names match reference format")
    else:
        bam_to_reference_map = None
        for conv_name, conv_map in KNOWN_CHROMOSOME_CONVENTIONS.items():
            if set(conv_map.keys()).issubset(bam_contig_names):
                bam_to_reference_map = conv_map
                print(f"  [OK] BAM uses '{conv_name}' names - auto-mapped")
                break
        if bam_to_reference_map is None:
            sys.exit(f"[ERROR] Unrecognized BAM chromosome names: {sorted(bam_contig_names)[:10]}")

    reference_to_bam_chromosome_map = {v: k for k, v in bam_to_reference_map.items()}

# -- Check SNV pickles (Part B only) --
if run_structural:
    if not SNV_PICKLE_DIRECTORY.is_dir():
        sys.exit(
            f"[ERROR] SNV pickle directory not found: {SNV_PICKLE_DIRECTORY}\n"
            f"  Run 'uv run orc4r-screen_02_analyze-bam.py' first."
        )
    missing_snv_pickles = []
    for experiment_number in experiment_numbers:
        p = SNV_PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
        if not p.exists():
            missing_snv_pickles.append(p.name)
    if len(missing_snv_pickles) > 0:
        sys.exit(
            f"[ERROR] Missing SNV pickles:\n"
            + "\n".join(f"  - {f}" for f in missing_snv_pickles)
        )
    print(f"  [OK] All {len(experiment_numbers)} SNV pickle files found")

# -- Check gene coordinates --
if not GENE_COORDINATES_FILE.exists() or not GENE_COORDINATES_METADATA_FILE.exists():
    sys.exit(
        f"[ERROR] Gene coordinate reference not found.\n"
        f"  Run 'uv run orc4r-screen_01_prepare-reference.py' first."
    )

reference_metadata = json.loads(GENE_COORDINATES_METADATA_FILE.read_text())
actual_sha256 = hashlib.sha256(GENE_COORDINATES_FILE.read_bytes()).hexdigest()
if actual_sha256 != reference_metadata["sha256"]:
    sys.exit(f"[ERROR] Gene coordinate integrity check failed.")
print(f"  [OK] Gene coordinates verified")

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
print(f"  Loaded {len(gene_coordinates):,} gene annotations")

gene_coordinates_by_chromosome = {
    chromosome: group.reset_index(drop=True)
    for chromosome, group in gene_coordinates.groupby("chromosome")
}
print()


# =============================================================================
# Shared helpers
# =============================================================================

def find_overlapping_genes(chromosome, region_start, region_end):
    """Return list of gene symbols overlapping a genomic region."""
    chrom_genes = gene_coordinates_by_chromosome.get(chromosome, pandas.DataFrame())
    if len(chrom_genes) == 0:
        return []
    overlapping = chrom_genes[
        (chrom_genes["start"] < region_end)
        & (chrom_genes["end"] > region_start)
    ]
    symbols = []
    for _, gene in overlapping.iterrows():
        name = gene.get("gene_symbol") or gene.get("gene_name") or gene.get("gene_alias")
        if pandas.notna(name):
            symbols.append(name)
    return symbols


def find_gene_at_position(chromosome, position):
    """Return first gene overlapping a single position, or None."""
    chrom_genes = gene_coordinates_by_chromosome.get(chromosome, pandas.DataFrame())
    if len(chrom_genes) == 0:
        return None
    overlapping = chrom_genes[
        (chrom_genes["start"] <= position) & (chrom_genes["end"] > position)
    ]
    if len(overlapping) == 0:
        return None
    first = overlapping.iloc[0]
    return first.get("gene_symbol") or first.get("gene_name") or first.get("gene_alias") or None


def cluster_positions(positions, max_gap):
    """Group sorted positions into clusters by max_gap. Returns [(start, end, count)]."""
    if len(positions) == 0:
        return []
    clusters = []
    start = positions[0]
    end = positions[0]
    count = 1
    for i in range(1, len(positions)):
        if positions[i] - end <= max_gap:
            end = positions[i]
            count += 1
        else:
            clusters.append((start, end, count))
            start = positions[i]
            end = positions[i]
            count = 1
    clusters.append((start, end, count))
    return clusters


# #############################################################################
#
#  PART A: INDEL DETECTION
#
# #############################################################################

if run_indels:
    print("#" * 80)
    print("# PART A: INDEL DETECTION")
    print("#" * 80)
    print()

    # =========================================================================
    # A1: Extract indel data from BAM files
    # =========================================================================

    print("=" * 80)
    print("A1: EXTRACT INDEL DATA FROM BAM FILES")
    print("=" * 80)

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

                for ref_chrom in REFERENCE_CHROMOSOMES:
                    bam_chrom = reference_to_bam_chromosome_map[ref_chrom]

                    if ref_chrom == bam_chrom:
                        print(f"    {suffix:<4} | {ref_chrom}")
                    else:
                        print(f"    {suffix:<4} | {ref_chrom} (BAM: {bam_chrom})")

                    for pileup_column in bam_file.pileup(contig=bam_chrom, truncate=True):
                        ins_count = 0
                        del_count = 0
                        indel_sizes = []

                        for pileup_read in pileup_column.pileups:
                            if pileup_read.indel > 0:
                                ins_count += 1
                                indel_sizes.append(pileup_read.indel)
                            if pileup_read.indel < 0:
                                indel_sizes.append(pileup_read.indel)
                            if pileup_read.is_del:
                                del_count += 1

                        if len(indel_sizes) > 0:
                            size_counter = collections.Counter(indel_sizes)
                            dominant_size = size_counter.most_common(1)[0][0]
                        else:
                            dominant_size = 0

                        chromosome_names.append(ref_chrom)
                        base_positions.append(pileup_column.pos)
                        read_coverages.append(pileup_column.n)
                        insertion_counts.append(ins_count)
                        deletion_counts.append(del_count)
                        dominant_indel_sizes.append(dominant_size)

                bam_file.close()

                sample_df = pandas.DataFrame({
                    "chromosome": chromosome_names,
                    "base_position": base_positions,
                    f"coverage_{suffix}": read_coverages,
                    f"insertion_count_{suffix}": insertion_counts,
                    f"deletion_count_{suffix}": deletion_counts,
                    f"dominant_indel_size_{suffix}": dominant_indel_sizes,
                })

                assert len(sample_df) > 0, f"No data from {bam_path.name}."
                print(f"    {suffix:<4} | {len(sample_df):,} positions extracted")
                sample_dataframes[suffix] = sample_df

            merged = pandas.merge(
                left=sample_dataframes["hit"],
                right=sample_dataframes["ctrl"],
                on=["chromosome", "base_position"],
                how="inner",
            )
            print(f"    Merged: {len(merged):,} shared positions")
            merged.to_pickle(str(pickle_path))
            print(f"    Cached: {pickle_path.name}")

    except KeyboardInterrupt:
        print("\n  [INTERRUPTED] Partial results cached. Re-run to resume.")
        sys.exit(1)

    print()

    # =========================================================================
    # A2: Filter for candidate indels
    # =========================================================================

    print("=" * 80)
    print("A2: FILTER FOR CANDIDATE INDELS")
    print("=" * 80)
    print(
        f"  Criteria: coverage > {MINIMUM_READ_COVERAGE}, "
        f"indel fraction > {MINIMUM_INDEL_FRACTION:.0%} in one sample, "
        f"< {MAXIMUM_BACKGROUND_FRACTION:.0%} in the other"
    )

    all_indel_candidates = []

    for experiment_number in experiment_numbers:
        pickle_path = INDEL_PICKLE_DIRECTORY / f"df_indel_exp{experiment_number}.pickle"
        merged = pandas.read_pickle(str(pickle_path))

        merged["indel_fraction_hit"] = numpy.where(
            merged["coverage_hit"] > 0,
            (merged["insertion_count_hit"] + merged["deletion_count_hit"]) / merged["coverage_hit"],
            0.0,
        )
        merged["indel_fraction_ctrl"] = numpy.where(
            merged["coverage_ctrl"] > 0,
            (merged["insertion_count_ctrl"] + merged["deletion_count_ctrl"]) / merged["coverage_ctrl"],
            0.0,
        )

        passes_cov = (
            (merged["coverage_hit"] > MINIMUM_READ_COVERAGE)
            & (merged["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        )

        indel_in_hit = (
            passes_cov
            & (merged["indel_fraction_hit"] > MINIMUM_INDEL_FRACTION)
            & (merged["indel_fraction_ctrl"] < MAXIMUM_BACKGROUND_FRACTION)
        )
        indel_in_ctrl = (
            passes_cov
            & (merged["indel_fraction_ctrl"] > MINIMUM_INDEL_FRACTION)
            & (merged["indel_fraction_hit"] < MAXIMUM_BACKGROUND_FRACTION)
        )

        is_candidate = indel_in_hit | indel_in_ctrl
        candidate_df = merged[is_candidate].copy()
        candidate_df["experiment"] = f"Exp_{experiment_number}"

        candidate_df["indel_sample"] = numpy.where(
            indel_in_hit[is_candidate].values, "hit", "ctrl"
        )

        candidate_df["indel_type"] = "unknown"
        for idx in candidate_df.index:
            sample = candidate_df.at[idx, "indel_sample"]
            ins = candidate_df.at[idx, f"insertion_count_{sample}"]
            dels = candidate_df.at[idx, f"deletion_count_{sample}"]
            if ins > dels:
                candidate_df.at[idx, "indel_type"] = "insertion"
            elif dels > ins:
                candidate_df.at[idx, "indel_type"] = "deletion"
            else:
                candidate_df.at[idx, "indel_type"] = "mixed"

        candidate_df["indel_size"] = numpy.where(
            candidate_df["indel_sample"] == "hit",
            candidate_df["dominant_indel_size_hit"],
            candidate_df["dominant_indel_size_ctrl"],
        )

        count = len(candidate_df)
        word = "position" if count == 1 else "positions"
        print(f"  Experiment {experiment_number}: {count:,} candidate {word}")
        all_indel_candidates.append(candidate_df)

    all_indel_df = pandas.concat(all_indel_candidates, ignore_index=True)
    total_indels = len(all_indel_df)
    print(f"  Total indel candidates: {total_indels:,}")

    # =========================================================================
    # A3: Annotate and save indel results
    # =========================================================================

    if total_indels > 0:
        print()
        print("=" * 80)
        print("A3: ANNOTATE INDEL CANDIDATES")
        print("=" * 80)

        all_indel_df["gene_symbol"] = pandas.NA
        for idx in range(total_indels):
            chrom = all_indel_df.at[idx, "chromosome"]
            pos = all_indel_df.at[idx, "base_position"]
            gene = find_gene_at_position(chrom, pos)
            all_indel_df.at[idx, "gene_symbol"] = gene or pandas.NA

        indel_output_cols = [
            "experiment", "chromosome", "base_position", "indel_type", "indel_size",
            "indel_sample", "indel_fraction_hit", "indel_fraction_ctrl",
            "coverage_hit", "coverage_ctrl", "gene_symbol",
        ]
        all_indel_df["indel_fraction_hit"] = all_indel_df["indel_fraction_hit"].round(3)
        all_indel_df["indel_fraction_ctrl"] = all_indel_df["indel_fraction_ctrl"].round(3)
        indel_csv_path = OUTPUT_DIRECTORY / INDEL_CSV_FILENAME
        all_indel_df[indel_output_cols].to_csv(str(indel_csv_path), index=False)
        print(f"  Saved {total_indels} candidates to: {indel_csv_path}")
    else:
        print("  [NOTE] No indel candidates found.")
        indel_csv_path = OUTPUT_DIRECTORY / INDEL_CSV_FILENAME
        pandas.DataFrame().to_csv(str(indel_csv_path), index=False)
        print(f"  Saved empty results to: {indel_csv_path}")

    print()


# #############################################################################
#
#  PART B: STRUCTURAL VARIANT ANALYSIS
#
# #############################################################################

if run_structural:
    print("#" * 80)
    print("# PART B: STRUCTURAL VARIANT ANALYSIS")
    print("#" * 80)
    print()

    active_targets = [n for n in experiment_numbers if n in TARGET_EXPERIMENT_NUMBERS]
    active_background = [n for n in experiment_numbers if n in BACKGROUND_EXPERIMENT_NUMBERS]

    if len(active_background) == 0:
        print(
            "  [WARNING] No background experiments included. Cross-experiment filtering\n"
            "    will be skipped. For best results, run all 9 experiments."
        )
    else:
        print(f"  Background experiments (known SNVs): {active_background}")
        print(f"  Target experiments (no SNVs): {active_targets}")

    print()

    # =========================================================================
    # B1: Detect candidate regions per experiment
    # =========================================================================

    print("=" * 80)
    print("B1: DETECT CANDIDATE REGIONS")
    print("=" * 80)

    all_candidate_regions = []

    for experiment_number in experiment_numbers:
        print(f"\n  --- Experiment {experiment_number} ---")
        exp_label = "TARGET" if experiment_number in TARGET_EXPERIMENT_NUMBERS else "BACKGROUND"
        print(f"  ({exp_label})")

        merged = pandas.read_pickle(
            str(SNV_PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle")
        )

        # -- Coverage imbalance --
        drop_in_hit = (
            (merged["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
            & (merged["coverage_hit"] < merged["coverage_ctrl"] * COVERAGE_IMBALANCE_RATIO)
        )
        drop_in_ctrl = (
            (merged["coverage_hit"] > MINIMUM_READ_COVERAGE)
            & (merged["coverage_ctrl"] < merged["coverage_hit"] * COVERAGE_IMBALANCE_RATIO)
        )

        for direction, mask in [("hit", drop_in_hit), ("ctrl", drop_in_ctrl)]:
            flagged = merged[mask]
            for chromosome in REFERENCE_CHROMOSOMES:
                positions = flagged[flagged["chromosome"] == chromosome]["base_position"].sort_values().tolist()
                if len(positions) == 0:
                    continue
                for reg_start, reg_end, pos_count in cluster_positions(positions, MAX_REGION_GAP_BP):
                    reg_length = reg_end - reg_start + 1
                    if reg_length < MINIMUM_REGION_LENGTH_BP:
                        continue
                    region_data = merged[
                        (merged["chromosome"] == chromosome)
                        & (merged["base_position"] >= reg_start)
                        & (merged["base_position"] <= reg_end)
                    ]
                    genes = find_overlapping_genes(chromosome, reg_start, reg_end)
                    all_candidate_regions.append({
                        "experiment": f"Exp_{experiment_number}",
                        "chromosome": chromosome,
                        "region_start": reg_start,
                        "region_end": reg_end,
                        "region_length_bp": reg_length,
                        "position_count": pos_count,
                        "detection_method": "coverage_imbalance",
                        "affected_sample": direction,
                        "mean_coverage_hit": round(region_data["coverage_hit"].mean(), 1),
                        "mean_coverage_ctrl": round(region_data["coverage_ctrl"].mean(), 1),
                        "genes": ", ".join(genes) if genes else "(intergenic)",
                    })

        # -- Mismatch clustering --
        bases_differ = merged["base_value_hit"] != merged["base_value_ctrl"]
        low_conf = (
            (merged["confidence_hit"] < MAXIMUM_CLUSTER_CONFIDENCE)
            | (merged["confidence_ctrl"] < MAXIMUM_CLUSTER_CONFIDENCE)
        )
        has_cov = (
            (merged["coverage_hit"] > MINIMUM_READ_COVERAGE)
            & (merged["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        )
        mismatches = merged[bases_differ & low_conf & has_cov]

        for chromosome in REFERENCE_CHROMOSOMES:
            chrom_mm = mismatches[mismatches["chromosome"] == chromosome].sort_values("base_position")
            if len(chrom_mm) == 0:
                continue
            positions = chrom_mm["base_position"].tolist()
            for cl_start, cl_end, cl_count in cluster_positions(positions, MAX_CLUSTER_GAP_BP):
                if cl_count < MIN_CLUSTER_SIZE:
                    continue
                cluster_data = chrom_mm[
                    (chrom_mm["base_position"] >= cl_start)
                    & (chrom_mm["base_position"] <= cl_end)
                ]
                genes = find_overlapping_genes(chromosome, cl_start, cl_end)
                all_candidate_regions.append({
                    "experiment": f"Exp_{experiment_number}",
                    "chromosome": chromosome,
                    "region_start": cl_start,
                    "region_end": cl_end,
                    "region_length_bp": cl_end - cl_start + 1,
                    "position_count": cl_count,
                    "detection_method": "mismatch_cluster",
                    "affected_sample": "hit",
                    "mean_coverage_hit": round(cluster_data["coverage_hit"].mean(), 1),
                    "mean_coverage_ctrl": round(cluster_data["coverage_ctrl"].mean(), 1),
                    "genes": ", ".join(genes) if genes else "(intergenic)",
                })

        # -- Positional gaps --
        for chromosome in REFERENCE_CHROMOSOMES:
            chrom_data = merged[merged["chromosome"] == chromosome].sort_values("base_position")
            if len(chrom_data) < 2:
                continue
            positions = chrom_data["base_position"].values
            cov_hit = chrom_data["coverage_hit"].values
            cov_ctrl = chrom_data["coverage_ctrl"].values

            for i in range(1, len(positions)):
                gap = positions[i] - positions[i - 1]
                if gap <= MIN_SUSPICIOUS_GAP_BP:
                    continue
                before_ok = cov_hit[i-1] > MIN_FLANKING_COVERAGE or cov_ctrl[i-1] > MIN_FLANKING_COVERAGE
                after_ok = cov_hit[i] > MIN_FLANKING_COVERAGE or cov_ctrl[i] > MIN_FLANKING_COVERAGE
                if not (before_ok and after_ok):
                    continue

                avg_hit = (cov_hit[i-1] + cov_hit[i]) / 2
                avg_ctrl = (cov_ctrl[i-1] + cov_ctrl[i]) / 2
                if avg_hit < avg_ctrl * 0.5:
                    likely = "hit"
                elif avg_ctrl < avg_hit * 0.5:
                    likely = "ctrl"
                else:
                    likely = "unclear"

                genes = find_overlapping_genes(chromosome, int(positions[i-1]), int(positions[i]))
                all_candidate_regions.append({
                    "experiment": f"Exp_{experiment_number}",
                    "chromosome": chromosome,
                    "region_start": int(positions[i-1]),
                    "region_end": int(positions[i]),
                    "region_length_bp": gap,
                    "position_count": 0,
                    "detection_method": "positional_gap",
                    "affected_sample": likely,
                    "mean_coverage_hit": round(avg_hit, 1),
                    "mean_coverage_ctrl": round(avg_ctrl, 1),
                    "genes": ", ".join(genes) if genes else "(intergenic)",
                })

    print()

    # =========================================================================
    # B2: Cross-experiment filtering
    # =========================================================================

    print("=" * 80)
    print("B2: CROSS-EXPERIMENT FILTERING")
    print("=" * 80)

    results_df = pandas.DataFrame(all_candidate_regions)

    if len(results_df) == 0:
        print("  No candidate regions to filter.")
        filtered_regions = pandas.DataFrame()
    else:
        target_labels = [f"Exp_{n}" for n in active_targets]
        background_labels = [f"Exp_{n}" for n in active_background]

        target_df = results_df[results_df["experiment"].isin(target_labels)].copy()
        background_df = results_df[results_df["experiment"].isin(background_labels)].copy()

        total_raw = len(target_df)
        print(f"  Raw candidate regions in target experiments: {total_raw:,}")
        print(f"  Regions in background experiments: {len(background_df):,}")

        if len(active_background) == 0 or len(background_df) == 0:
            print("  [SKIP] No background for filtering.")
            filtered_regions = target_df.copy()
        else:
            is_shared = []
            for _, row in target_df.iterrows():
                same_chrom = background_df[background_df["chromosome"] == row["chromosome"]]
                shared = False
                for _, bg in same_chrom.iterrows():
                    if (row["region_start"] - REGION_PROXIMITY_BP <= bg["region_end"]
                            and bg["region_start"] - REGION_PROXIMITY_BP <= row["region_end"]):
                        shared = True
                        break
                is_shared.append(shared)

            shared_count = sum(is_shared)
            print(f"  Shared with background (removed): {shared_count:,}")
            print(f"  Unique to target experiments: {total_raw - shared_count:,}")

            filtered_regions = target_df[
                ~pandas.Series(is_shared, index=target_df.index)
            ].copy()

    print()

    # =========================================================================
    # B3: Save structural variant results
    # =========================================================================

    print("=" * 80)
    print("B3: SAVE STRUCTURAL VARIANT RESULTS")
    print("=" * 80)

    sv_csv_path = OUTPUT_DIRECTORY / STRUCTURAL_CSV_FILENAME

    if len(filtered_regions) == 0:
        print("  No unique structural variant candidates after filtering.")
        pandas.DataFrame().to_csv(str(sv_csv_path), index=False)
    else:
        chrom_sort = {c: i for i, c in enumerate(REFERENCE_CHROMOSOMES)}
        filtered_sorted = filtered_regions.sort_values(
            by=["experiment", "chromosome", "region_start"],
            key=lambda col: col.map(chrom_sort) if col.name == "chromosome" else col,
            ignore_index=True,
        )
        filtered_sorted.to_csv(str(sv_csv_path), index=False)
        print(f"  Saved {len(filtered_sorted)} candidate regions to: {sv_csv_path}")

        print()
        print("  --- Summary ---")
        for exp_number in active_targets:
            exp_label = f"Exp_{exp_number}"
            exp_regions = filtered_sorted[filtered_sorted["experiment"] == exp_label]
            if len(exp_regions) == 0:
                print(f"  Experiment {exp_number}: no unique candidates")
            else:
                print(f"  Experiment {exp_number}: {len(exp_regions)} unique candidate(s)")
                for _, row in exp_regions.iterrows():
                    print(
                        f"    [{row['detection_method'].replace('_', ' ')}] "
                        f"{row['chromosome']}:{row['region_start']}-{row['region_end']} "
                        f"({row['region_length_bp']:,} bp) "
                        f"-> {row['genes']}"
                    )

    print()


# =============================================================================
# Final summary
# =============================================================================

print("=" * 80)
print("COMPLETE")
print("=" * 80)

elapsed_seconds = time.perf_counter() - start_time
print(f"Total runtime: {elapsed_seconds:.1f} seconds")
print("Done.")
