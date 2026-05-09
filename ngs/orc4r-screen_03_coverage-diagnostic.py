# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "numpy",
# ]
# ///
"""
Comprehensive diagnostic for the ORC4 suppressor screen.

Reads cached pickle files from orc4r-screen_02_analyze-bam.py and performs:
  1) Per-experiment coverage summary
  2) Filter funnel analysis - where candidates are lost
  3) Near-miss analysis - characterizes positions that almost passed filters
  4) Aneuploidy check - per-chromosome coverage ratios across all experiments
  5) Gene catalog - lists genes on any chromosome showing gain or loss

This diagnostic traces the analytical path from "why do some experiments
produce zero SNV candidates?" through to the discovery of chromosome X
disomy as an alternative suppressor mechanism.

Usage:
    uv run orc4r-screen_03_diagnostic.py
    uv run orc4r-screen_03_diagnostic.py --experiments 1,7,9
    uv run orc4r-screen_03_diagnostic.py --help

Prerequisites:
    uv run orc4r-screen_02_analyze-bam.py   (generates cached pickles)
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

# Must match thresholds in orc4r-screen_02_analyze-bam.py
MINIMUM_READ_COVERAGE = 10
MINIMUM_CONFIDENCE_PERCENT = 90

# Aneuploidy detection thresholds
# Normalized ratio = (chromosome hit/ctrl ratio) / (genome-wide hit/ctrl ratio)
# Values near 2.0 suggest chromosome gain; near 0.5 suggest loss.
ANEUPLOIDY_GAIN_THRESHOLD = 1.5
ANEUPLOIDY_LOSS_THRESHOLD = 0.67

# Replication-related search terms for gene description filtering
REPLICATION_SEARCH_TERMS = [
    "orc ",
    "origin recognition",
    "replication",
    "mcm",
    "cdc6",
    "helicase",
    "licensing",
    "prereplicat",
    "minichromosome",
]

# Specific genes of interest for the ORC4 suppressor screen
GENES_OF_INTEREST = [
    "ORC1",
    "ORC2",
    "ORC3",
    "ORC4",
    "ORC5",
    "ORC6",
    "CDC6",
    "CDT1",
    "TAH11",
    "MCM2",
    "MCM3",
    "MCM4",
    "MCM5",
    "MCM6",
    "MCM7",
    "DPB11",
    "PSF2",
    "POL31",
    "POL32",
    "RFC2",
    "RFA3",
]

OUTPUT_CSV_FILENAME = "coverage_diagnostic.csv"
OUTPUT_REPORT_FILENAME = "diagnostic_report.txt"


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

# Accumulate report text for file output
report_lines = []


def out(text=""):
    """Print to stdout and accumulate for report file."""
    print(text)
    report_lines.append(text)


# =============================================================================
# Preflight checks
# =============================================================================

out("=" * 80)
out("PREFLIGHT CHECKS")
out("=" * 80)

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
        f"[ERROR] Missing pickle files in {PICKLE_DIRECTORY}:\n"
        + "\n".join(f"  - {filename}" for filename in missing_pickles)
        + f"\n  Run 'uv run orc4r-screen_02_analyze-bam.py' first."
    )
out(f"  [OK] All {len(experiment_numbers)} pickle files found")

# Check gene coordinate reference
has_reference = False
if GENE_COORDINATES_FILE.exists() and GENE_COORDINATES_METADATA_FILE.exists():
    reference_metadata = json.loads(GENE_COORDINATES_METADATA_FILE.read_text())
    expected_sha256 = reference_metadata["sha256"]
    actual_sha256 = hashlib.sha256(GENE_COORDINATES_FILE.read_bytes()).hexdigest()
    if actual_sha256 == expected_sha256:
        has_reference = True
        out(f"  [OK] Gene coordinates verified")
    else:
        out(
            f"  [WARNING] Gene coordinate file integrity check failed. Gene catalog will be skipped."
        )
else:
    out(f"  [NOTE] Gene coordinates not found. Gene catalog will be skipped.")
    out(f"    Run 'uv run orc4r-screen_01_prepare-reference.py' to enable.")

out()


# =============================================================================
# Section 1: Per-experiment coverage summary
# =============================================================================

out("=" * 80)
out("SECTION 1: COVERAGE SUMMARY PER EXPERIMENT")
out("=" * 80)

experiment_summaries = []

for experiment_number in experiment_numbers:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    total_positions = len(merged_dataframe)

    summary = {
        "experiment": experiment_number,
        "total_positions": total_positions,
        "median_coverage_hit": merged_dataframe["coverage_hit"].median(),
        "median_coverage_ctrl": merged_dataframe["coverage_ctrl"].median(),
        "mean_coverage_hit": round(merged_dataframe["coverage_hit"].mean(), 1),
        "mean_coverage_ctrl": round(merged_dataframe["coverage_ctrl"].mean(), 1),
        "pct_above_cov_hit": round(
            (merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE).sum()
            / total_positions
            * 100,
            1,
        ),
        "pct_above_cov_ctrl": round(
            (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE).sum()
            / total_positions
            * 100,
            1,
        ),
        "pct_above_conf_hit": round(
            (merged_dataframe["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT).sum()
            / total_positions
            * 100,
            1,
        ),
        "pct_above_conf_ctrl": round(
            (merged_dataframe["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT).sum()
            / total_positions
            * 100,
            1,
        ),
    }
    experiment_summaries.append(summary)

out(
    f"\n  {'Exp':<5} {'Positions':>12} {'Med Cov Hit':>12} {'Med Cov Ctrl':>13} "
    f"{'%>Cov Hit':>10} {'%>Cov Ctrl':>11} {'%>Conf Hit':>11} {'%>Conf Ctrl':>12}"
)
out(
    f"  {'---':<5} {'-' * 10:>12} {'-' * 10:>12} {'-' * 10:>13} "
    f"{'-' * 8:>10} {'-' * 8:>11} {'-' * 8:>11} {'-' * 8:>12}"
)

for summary in experiment_summaries:
    out(
        f"  {summary['experiment']:<5} "
        f"{summary['total_positions']:>12,} "
        f"{summary['median_coverage_hit']:>12.0f} "
        f"{summary['median_coverage_ctrl']:>13.0f} "
        f"{summary['pct_above_cov_hit']:>9.1f}% "
        f"{summary['pct_above_cov_ctrl']:>10.1f}% "
        f"{summary['pct_above_conf_hit']:>10.1f}% "
        f"{summary['pct_above_conf_ctrl']:>11.1f}%"
    )

out()


# =============================================================================
# Section 2: Filter funnel analysis
#
# For each experiment, show how many positions survive each successive filter.
# This reveals WHERE candidates are lost for zero-hit experiments.
# =============================================================================

out("=" * 80)
out("SECTION 2: FILTER FUNNEL ANALYSIS")
out("=" * 80)
out(
    f"  Thresholds: coverage > {MINIMUM_READ_COVERAGE}, "
    f"confidence > {MINIMUM_CONFIDENCE_PERCENT}%"
)

funnel_rows = []

for experiment_number in experiment_numbers:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    total_positions = len(merged_dataframe)

    passes_coverage_both = (
        merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE
    ) & (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
    after_coverage = passes_coverage_both.sum()

    passes_confidence_both = (
        passes_coverage_both
        & (merged_dataframe["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT)
        & (merged_dataframe["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT)
    )
    after_confidence = passes_confidence_both.sum()

    bases_differ = (
        merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]
    )
    after_differ = (passes_confidence_both & bases_differ).sum()

    funnel_rows.append(
        {
            "experiment": experiment_number,
            "total_positions": total_positions,
            "pass_coverage": after_coverage,
            "pass_confidence": after_confidence,
            "pass_differ": after_differ,
        }
    )

out(
    f"\n  {'Exp':<5} {'Total':>12} {'Pass Cov':>12} {'Pass Conf':>12} {'Bases Differ':>13}"
)
out(f"  {'---':<5} {'-' * 10:>12} {'-' * 10:>12} {'-' * 10:>12} {'-' * 10:>13}")

for row in funnel_rows:
    marker = " <--" if row["pass_differ"] == 0 else ""
    out(
        f"  {row['experiment']:<5} "
        f"{row['total_positions']:>12,} "
        f"{row['pass_coverage']:>12,} "
        f"{row['pass_confidence']:>12,} "
        f"{row['pass_differ']:>13,}"
        f"{marker}"
    )

zero_hit_experiments = [
    row["experiment"] for row in funnel_rows if row["pass_differ"] == 0
]
has_hit_experiments = [
    row["experiment"] for row in funnel_rows if row["pass_differ"] > 0
]

out()
if len(zero_hit_experiments) > 0:
    out(f"  Experiments with zero SNV candidates: {zero_hit_experiments}")
    out(f"  Experiments with SNV candidates: {has_hit_experiments}")
else:
    out("  All experiments produced SNV candidates.")

out()


# =============================================================================
# Section 3: Near-miss analysis (zero-hit experiments)
#
# For experiments with zero candidates, characterize positions where bases
# differ but failed quality filters. The confidence distribution at these
# positions reveals whether the mutations are SNVs with poor data or
# structural changes causing read misalignment.
# =============================================================================

if len(zero_hit_experiments) > 0:
    out("=" * 80)
    out("SECTION 3: NEAR-MISS ANALYSIS (ZERO-HIT EXPERIMENTS)")
    out("=" * 80)

    for experiment_number in zero_hit_experiments:
        pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
        merged_dataframe = pandas.read_pickle(str(pickle_path))

        bases_differ = (
            merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]
        )
        differing_positions = merged_dataframe[bases_differ].copy()
        total_differing = len(differing_positions)

        out(
            f"\n  Experiment {experiment_number}: {total_differing:,} positions with different bases (before quality filters)"
        )

        if total_differing == 0:
            out("    No base differences at any quality level.")
            continue

        has_cov_both = (differing_positions["coverage_hit"] > MINIMUM_READ_COVERAGE) & (
            differing_positions["coverage_ctrl"] > MINIMUM_READ_COVERAGE
        )
        covered_diffs = has_cov_both.sum()

        has_conf_hit = has_cov_both & (
            differing_positions["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT
        )
        has_conf_both = has_conf_hit & (
            differing_positions["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT
        )

        out(f"    Differing positions:                  {total_differing:>8,}")
        out(
            f"    With coverage > {MINIMUM_READ_COVERAGE} (both):             {covered_diffs:>8,}"
        )
        out(
            f"    + confidence > {MINIMUM_CONFIDENCE_PERCENT}% (hit):            {has_conf_hit.sum():>8,}"
        )
        out(
            f"    + confidence > {MINIMUM_CONFIDENCE_PERCENT}% (both):           {has_conf_both.sum():>8,}"
        )

        if has_conf_both.sum() == 0 and covered_diffs > 0:
            low_conf_diffs = differing_positions[has_cov_both]
            median_conf_hit = low_conf_diffs["confidence_hit"].median()
            median_conf_ctrl = low_conf_diffs["confidence_ctrl"].median()
            out(
                f"\n    Covered but low-confidence: median confidence at differing positions: "
                f"hit={median_conf_hit:.0f}%, ctrl={median_conf_ctrl:.0f}%"
            )
            out(
                f"    ~50-60% confidence indicates reads split roughly equally between bases."
            )
            out(
                f"    This pattern is consistent with structural variants causing read misalignment,"
            )
            out(f"    or aneuploidy affecting coverage ratios genome-wide.")

    out()


# =============================================================================
# Section 4: Aneuploidy check (all experiments)
#
# For each experiment, compute the per-chromosome coverage ratio (hit/ctrl)
# normalized by the genome-wide ratio. A chromosome with normalized ratio
# near 2.0 suggests gain (disomy); near 0.5 suggests loss.
#
# This check runs on ALL experiments so that any chromosome gains unique to
# zero-hit experiments can be identified by contrast.
# =============================================================================

out("=" * 80)
out("SECTION 4: ANEUPLOIDY CHECK (ALL EXPERIMENTS)")
out("=" * 80)
out()
out(f"  Gain threshold: normalized ratio > {ANEUPLOIDY_GAIN_THRESHOLD}")
out(f"  Loss threshold: normalized ratio < {ANEUPLOIDY_LOSS_THRESHOLD}")
out()

all_normalized_ratios = {}

for experiment_number in experiment_numbers:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    global_med_hit = merged_dataframe["coverage_hit"].median()
    global_med_ctrl = merged_dataframe["coverage_ctrl"].median()
    global_ratio = global_med_hit / max(global_med_ctrl, 1)

    out(
        f"  --- Experiment {experiment_number} (genome-wide hit/ctrl ratio: {global_ratio:.2f}) ---"
    )
    out(
        f"    {'Chromosome':<12} {'Med Hit':>8} {'Med Ctrl':>9} {'Ratio':>6} {'Norm':>6} {'Flag':>8}"
    )
    out(
        f"    {'-' * 11:<12} {'-' * 7:>8} {'-' * 8:>9} {'-' * 5:>6} {'-' * 5:>6} {'-' * 7:>8}"
    )

    exp_ratios = {}
    for chromosome in REFERENCE_CHROMOSOMES:
        chrom_data = merged_dataframe[merged_dataframe["chromosome"] == chromosome]
        med_hit = chrom_data["coverage_hit"].median()
        med_ctrl = chrom_data["coverage_ctrl"].median()
        chrom_ratio = med_hit / max(med_ctrl, 1)
        normalized = chrom_ratio / max(global_ratio, 0.01)

        flag = ""
        if normalized > ANEUPLOIDY_GAIN_THRESHOLD:
            flag = "GAIN?"
        elif normalized < ANEUPLOIDY_LOSS_THRESHOLD:
            flag = "LOSS?"

        exp_ratios[chromosome] = normalized
        out(
            f"    {chromosome:<12} {med_hit:>8.0f} {med_ctrl:>9.0f} {chrom_ratio:>6.2f} {normalized:>6.2f} {flag:>8}"
        )

    all_normalized_ratios[experiment_number] = exp_ratios
    out()

# -- Summary: flagged chromosomes across all experiments --
out("  --- SUMMARY: Chromosomes with GAIN or LOSS flags ---")

any_flags = False
for chromosome in REFERENCE_CHROMOSOMES:
    flagged = []
    for experiment_number in experiment_numbers:
        norm = all_normalized_ratios[experiment_number][chromosome]
        if norm > ANEUPLOIDY_GAIN_THRESHOLD:
            flagged.append(f"Exp {experiment_number} (GAIN, {norm:.2f})")
        elif norm < ANEUPLOIDY_LOSS_THRESHOLD:
            flagged.append(f"Exp {experiment_number} (LOSS, {norm:.2f})")
    if flagged:
        any_flags = True
        out(f"    {chromosome:<12} {', '.join(flagged)}")

if not any_flags:
    out("    No chromosomes flagged in any experiment.")

out()

# -- Summary table: correlation between aneuploidy and SNV results --
flagged_chromosomes = set()
for chromosome in REFERENCE_CHROMOSOMES:
    for experiment_number in experiment_numbers:
        norm = all_normalized_ratios[experiment_number][chromosome]
        if norm > ANEUPLOIDY_GAIN_THRESHOLD or norm < ANEUPLOIDY_LOSS_THRESHOLD:
            flagged_chromosomes.add(chromosome)

if len(flagged_chromosomes) > 0:
    out("  --- Correlation: aneuploidy vs SNV detection ---")
    header_chroms = "  ".join(f"{c:>6}" for c in sorted(flagged_chromosomes))
    out(f"    {'Exp':>4}  {'SNV hits':>9}  {header_chroms}")

    for experiment_number in experiment_numbers:
        snv_status = "yes" if experiment_number in has_hit_experiments else "no"
        chrom_values = "  ".join(
            f"{all_normalized_ratios[experiment_number][c]:>6.2f}"
            for c in sorted(flagged_chromosomes)
        )
        out(f"    {experiment_number:>4}  {snv_status:>9}  {chrom_values}")

    out()


# =============================================================================
# Section 5: Gene catalog for flagged chromosomes
#
# For any chromosome showing gain or loss in the zero-hit experiments,
# list all genes and highlight replication-related genes that could
# explain suppression of the ORC4-ts phenotype.
# =============================================================================

# Identify chromosomes gained/lost specifically in zero-hit experiments
zero_hit_flagged_chromosomes = {}
for chromosome in REFERENCE_CHROMOSOMES:
    for experiment_number in zero_hit_experiments:
        norm = all_normalized_ratios[experiment_number][chromosome]
        if norm > ANEUPLOIDY_GAIN_THRESHOLD:
            zero_hit_flagged_chromosomes.setdefault(chromosome, []).append(
                (experiment_number, "GAIN", norm)
            )
        elif norm < ANEUPLOIDY_LOSS_THRESHOLD:
            zero_hit_flagged_chromosomes.setdefault(chromosome, []).append(
                (experiment_number, "LOSS", norm)
            )

if has_reference and len(zero_hit_flagged_chromosomes) > 0:
    out("=" * 80)
    out("SECTION 5: GENE CATALOG FOR FLAGGED CHROMOSOMES")
    out("=" * 80)
    out()

    gene_coordinates = pandas.read_csv(GENE_COORDINATES_FILE, sep="\t", dtype=str)

    for chromosome, flags in sorted(zero_hit_flagged_chromosomes.items()):
        flag_summary = ", ".join(
            f"Exp {exp} ({direction}, {norm:.2f})" for exp, direction, norm in flags
        )
        out(f"  --- {chromosome}: {flag_summary} ---")
        out()

        chrom_genes = gene_coordinates[
            gene_coordinates["chromosome"] == chromosome
        ].copy()
        chrom_genes["start"] = chrom_genes["start"].astype(int)
        chrom_genes = chrom_genes.sort_values("start")

        out(f"    Total genes on {chromosome}: {len(chrom_genes)}")
        out()

        # 5a: Replication-related genes on this chromosome
        out(f"    Replication-related genes on {chromosome}:")
        out(f"    {'Symbol':<12} {'Name':<12} {'Start':>8} {'End':>8}  Description")
        out(f"    {'-' * 11:<12} {'-' * 11:<12} {'-' * 7:>8} {'-' * 7:>8}  -----------")

        replication_gene_count = 0
        for _, gene in chrom_genes.iterrows():
            desc = str(gene.get("description", "") or "").lower()
            symbol = str(gene.get("gene_symbol", "") or "").lower()
            name = str(gene.get("gene_name", "") or "").lower()
            searchable = f"{desc} {symbol} {name}"

            is_replication = any(
                term in searchable for term in REPLICATION_SEARCH_TERMS
            )
            if is_replication:
                replication_gene_count += 1
                sym_display = gene.get("gene_symbol", "") or ""
                name_display = gene.get("gene_name", "") or ""
                desc_display = str(gene.get("description", "") or "")
                desc_short = (
                    desc_display[:65] + "..."
                    if len(desc_display) > 65
                    else desc_display
                )
                out(
                    f"    {sym_display:<12} {name_display:<12} {gene['start']:>8} {gene['end']:>8}  {desc_short}"
                )

        if replication_gene_count == 0:
            out(f"    (none found)")

        out()

        # 5b: Specific genes of interest on this chromosome
        out(f"    Genes of interest on {chromosome}:")
        out(f"    {'Gene':<12} {'Location':>20}  Description")
        out(f"    {'-' * 11:<12} {'-' * 19:>20}  -----------")

        found_any = False
        for target in GENES_OF_INTEREST:
            matches = chrom_genes[
                (chrom_genes["gene_symbol"].fillna("").str.upper() == target)
                | (chrom_genes["gene_name"].fillna("").str.upper() == target)
            ]
            for _, gene in matches.iterrows():
                found_any = True
                sym_display = gene.get("gene_symbol", "") or gene.get("gene_name", "")
                desc_display = str(gene.get("description", "") or "")
                desc_short = (
                    desc_display[:60] + "..."
                    if len(desc_display) > 60
                    else desc_display
                )
                out(
                    f"    {sym_display:<12} {gene['start']:>8}-{gene['end']:<8}  {desc_short}"
                )

        if not found_any:
            out(f"    (none found)")

        out()

    # 5c: All genes of interest across all chromosomes (for reference)
    out("  --- All genes of interest: chromosomal locations ---")
    out(
        f"    {'Gene':<12} {'Chrom':<8} {'Start':>8} {'End':>8}  {'On flagged chrom?':>18}"
    )
    out(f"    {'-' * 11:<12} {'-' * 7:<8} {'-' * 7:>8} {'-' * 7:>8}  {'-' * 17:>18}")

    for target in GENES_OF_INTEREST:
        matches = gene_coordinates[
            (gene_coordinates["gene_symbol"].fillna("").str.upper() == target)
            | (gene_coordinates["gene_name"].fillna("").str.upper() == target)
        ]
        if len(matches) > 0:
            for _, gene in matches.iterrows():
                sym_display = gene.get("gene_symbol", "") or gene.get("gene_name", "")
                chrom = gene.get("chromosome", "")
                on_flagged = "<<<" if chrom in zero_hit_flagged_chromosomes else ""
                out(
                    f"    {sym_display:<12} {chrom:<8} {gene['start']:>8} {gene['end']:>8}  {on_flagged:>18}"
                )
        else:
            out(f"    {target:<12} NOT FOUND")

    out()

elif len(zero_hit_flagged_chromosomes) > 0 and not has_reference:
    out("  [NOTE] Gene catalog skipped - reference file not available.")
    out()


# =============================================================================
# Save results
# =============================================================================

out("=" * 80)
out("SAVE RESULTS")
out("=" * 80)

# Save coverage diagnostic CSV
summary_dataframe = pandas.DataFrame(experiment_summaries)
funnel_dataframe = pandas.DataFrame(funnel_rows)
diagnostic_dataframe = pandas.merge(
    left=summary_dataframe,
    right=funnel_dataframe,
    on="experiment",
    how="inner",
    suffixes=("_summary", "_funnel"),
)

if "total_positions_summary" in diagnostic_dataframe.columns:
    diagnostic_dataframe = diagnostic_dataframe.rename(
        columns={"total_positions_summary": "total_positions"}
    )
    diagnostic_dataframe = diagnostic_dataframe.drop(
        columns=["total_positions_funnel"], errors="ignore"
    )

csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
diagnostic_dataframe.to_csv(str(csv_path), index=False)
out(f"  Saved diagnostic CSV: {csv_path}")

# Save full text report
report_path = OUTPUT_DIRECTORY / OUTPUT_REPORT_FILENAME
report_path.write_text("\n".join(report_lines) + "\n")
out(f"  Saved diagnostic report: {report_path}")

out()
elapsed_seconds = time.perf_counter() - start_time
out(f"Total runtime: {elapsed_seconds:.1f} seconds")
out("Done.")
