# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas",
#     "numpy",
# ]
# ///
"""
Coverage and filter diagnostic for the ORC4 suppressor screen.

Reads cached pickle files from orc4r-screen_02_analyze-bam.py and reports
per-experiment coverage statistics, filter funnel analysis, and a comparison
between experiments that yielded candidates and those that did not.

This diagnostic helps determine why some experiments (e.g., 1, 7, 9) produce
zero SNV candidates, guiding whether to pursue indel detection or investigate
data quality issues.

Usage:
    uv run orc4r-screen_03_coverage-diagnostic.py
    uv run orc4r-screen_03_coverage-diagnostic.py --help
    uv run orc4r-screen_03_coverage-diagnostic.py --experiments 1,7,9

Prerequisites:
    uv run orc4r-screen_02_analyze-bam.py   (generates cached pickles)
"""

import argparse
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

DEFAULT_EXPERIMENT_NUMBERS = list(range(1, 10))

REFERENCE_CHROMOSOMES = [
    "chrI",    "chrII",   "chrIII",  "chrIV",
    "chrV",    "chrVI",   "chrVII",  "chrVIII",
    "chrIX",   "chrX",    "chrXI",   "chrXII",
    "chrXIII", "chrXIV",  "chrXV",   "chrXVI",
]

# Must match thresholds in orc4r-screen_02_analyze-bam.py
MINIMUM_READ_COVERAGE = 10
MINIMUM_CONFIDENCE_PERCENT = 90

OUTPUT_CSV_FILENAME = "coverage_diagnostic.csv"


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
        f"[ERROR] Missing pickle files in {PICKLE_DIRECTORY}:\n"
        + "\n".join(f"  - {filename}" for filename in missing_pickles)
        + f"\n  Run 'uv run orc4r-screen_02_analyze-bam.py' first."
    )

print(f"  [OK] All {len(experiment_numbers)} pickle files found")
print()


# =============================================================================
# Section 1: Per-experiment coverage summary
#
# For each experiment, report total positions, median coverage, and the
# fraction of positions meeting the coverage and confidence thresholds.
# =============================================================================

print("=" * 60)
print("SECTION 1: COVERAGE SUMMARY PER EXPERIMENT")
print("=" * 60)

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
            / total_positions * 100, 1
        ),
        "pct_above_cov_ctrl": round(
            (merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE).sum()
            / total_positions * 100, 1
        ),
        "pct_above_conf_hit": round(
            (merged_dataframe["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT).sum()
            / total_positions * 100, 1
        ),
        "pct_above_conf_ctrl": round(
            (merged_dataframe["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT).sum()
            / total_positions * 100, 1
        ),
    }
    experiment_summaries.append(summary)

summary_dataframe = pandas.DataFrame(experiment_summaries)

print(f"\n  {'Exp':<5} {'Positions':>12} {'Med Cov Hit':>12} {'Med Cov Ctrl':>13} "
      f"{'%>Cov Hit':>10} {'%>Cov Ctrl':>11} {'%>Conf Hit':>11} {'%>Conf Ctrl':>12}")
print(f"  {'-'*4:<5} {'-'*10:>12} {'-'*10:>12} {'-'*10:>13} "
      f"{'-'*8:>10} {'-'*8:>11} {'-'*8:>11} {'-'*8:>12}")

for summary in experiment_summaries:
    print(
        f"  {summary['experiment']:<5} "
        f"{summary['total_positions']:>12,} "
        f"{summary['median_coverage_hit']:>12.0f} "
        f"{summary['median_coverage_ctrl']:>13.0f} "
        f"{summary['pct_above_cov_hit']:>9.1f}% "
        f"{summary['pct_above_cov_ctrl']:>10.1f}% "
        f"{summary['pct_above_conf_hit']:>10.1f}% "
        f"{summary['pct_above_conf_ctrl']:>11.1f}%"
    )

print()


# =============================================================================
# Section 2: Filter funnel analysis
#
# For each experiment, show how many positions survive each successive filter
# criterion. This reveals WHERE candidates are lost for experiments that
# produce zero hits.
#
# Funnel stages:
#   All positions
#   -> Pass coverage (both > 10)
#   -> Pass confidence (both > 90%)
#   -> Bases differ between hit and ctrl
# =============================================================================

print("=" * 60)
print("SECTION 2: FILTER FUNNEL ANALYSIS")
print("=" * 60)
print(
    f"  Thresholds: coverage > {MINIMUM_READ_COVERAGE}, "
    f"confidence > {MINIMUM_CONFIDENCE_PERCENT}%"
)

funnel_rows = []

for experiment_number in experiment_numbers:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    total_positions = len(merged_dataframe)

    passes_coverage_hit = merged_dataframe["coverage_hit"] > MINIMUM_READ_COVERAGE
    passes_coverage_ctrl = merged_dataframe["coverage_ctrl"] > MINIMUM_READ_COVERAGE
    passes_coverage_both = passes_coverage_hit & passes_coverage_ctrl
    after_coverage = passes_coverage_both.sum()

    passes_confidence_hit = merged_dataframe["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT
    passes_confidence_ctrl = merged_dataframe["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT
    after_confidence = (passes_coverage_both & passes_confidence_hit & passes_confidence_ctrl).sum()

    bases_differ = merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]
    after_differ = (
        passes_coverage_both & passes_confidence_hit & passes_confidence_ctrl & bases_differ
    ).sum()

    funnel_rows.append({
        "experiment": experiment_number,
        "total_positions": total_positions,
        "pass_coverage": after_coverage,
        "pass_confidence": after_confidence,
        "pass_differ": after_differ,
    })

print(f"\n  {'Exp':<5} {'Total':>12} {'Pass Cov':>12} {'Pass Conf':>12} {'Bases Differ':>13}")
print(f"  {'-'*4:<5} {'-'*10:>12} {'-'*10:>12} {'-'*10:>12} {'-'*10:>13}")

for row in funnel_rows:
    differ_marker = " <--" if row["pass_differ"] == 0 else ""
    print(
        f"  {row['experiment']:<5} "
        f"{row['total_positions']:>12,} "
        f"{row['pass_coverage']:>12,} "
        f"{row['pass_confidence']:>12,} "
        f"{row['pass_differ']:>13,}"
        f"{differ_marker}"
    )

print()

# Identify experiments with zero candidates
zero_hit_experiments = [row["experiment"] for row in funnel_rows if row["pass_differ"] == 0]
has_candidates_experiments = [row["experiment"] for row in funnel_rows if row["pass_differ"] > 0]

if len(zero_hit_experiments) > 0:
    print(f"  Experiments with zero SNV candidates: {zero_hit_experiments}")

    # Diagnose: do the zero-hit experiments still have positions passing coverage + confidence?
    zero_hit_pass_confidence = [
        row for row in funnel_rows
        if row["experiment"] in zero_hit_experiments and row["pass_confidence"] > 0
    ]

    if len(zero_hit_pass_confidence) == len(zero_hit_experiments):
        print(
            "  All zero-hit experiments have well-covered, high-confidence positions\n"
            "  but no base differences. This suggests the suppressor mutations may be\n"
            "  insertions or deletions (indels) rather than single nucleotide variants."
        )
    else:
        low_quality_experiments = [
            row["experiment"] for row in funnel_rows
            if row["experiment"] in zero_hit_experiments and row["pass_confidence"] == 0
        ]
        if len(low_quality_experiments) > 0:
            print(
                f"  Experiments {low_quality_experiments} have no positions passing both\n"
                f"  coverage and confidence filters. Data quality may be insufficient."
            )

print()


# =============================================================================
# Section 3: Per-chromosome coverage for zero-hit experiments
#
# Break down coverage by chromosome for experiments that yielded no candidates.
# This can reveal if specific chromosomes have unusually low coverage, suggesting
# alignment issues or library preparation problems.
# =============================================================================

if len(zero_hit_experiments) > 0:
    print("=" * 60)
    print("SECTION 3: PER-CHROMOSOME DETAIL (ZERO-HIT EXPERIMENTS)")
    print("=" * 60)

    for experiment_number in zero_hit_experiments:
        pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
        merged_dataframe = pandas.read_pickle(str(pickle_path))

        print(f"\n  Experiment {experiment_number}:")
        print(f"  {'Chromosome':<10} {'Positions':>10} {'Med Cov Hit':>12} {'Med Cov Ctrl':>13} "
              f"{'%>Cov Both':>11} {'%>Conf Both':>12} {'Bases Diff':>11}")
        print(f"  {'-'*9:<10} {'-'*9:>10} {'-'*10:>12} {'-'*10:>13} "
              f"{'-'*9:>11} {'-'*9:>12} {'-'*9:>11}")

        chromosome_sort_order = {c: i for i, c in enumerate(REFERENCE_CHROMOSOMES)}

        for chromosome in sorted(
            merged_dataframe["chromosome"].unique(),
            key=lambda c: chromosome_sort_order.get(c, 999),
        ):
            chrom_data = merged_dataframe[merged_dataframe["chromosome"] == chromosome]
            chrom_count = len(chrom_data)

            passes_cov = (
                (chrom_data["coverage_hit"] > MINIMUM_READ_COVERAGE)
                & (chrom_data["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
            )
            passes_conf = (
                passes_cov
                & (chrom_data["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT)
                & (chrom_data["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT)
            )
            bases_diff = (
                passes_conf
                & (chrom_data["base_value_hit"] != chrom_data["base_value_ctrl"])
            )

            pct_cov = round(passes_cov.sum() / chrom_count * 100, 1) if chrom_count > 0 else 0
            pct_conf = round(passes_conf.sum() / chrom_count * 100, 1) if chrom_count > 0 else 0

            print(
                f"  {chromosome:<10} "
                f"{chrom_count:>10,} "
                f"{chrom_data['coverage_hit'].median():>12.0f} "
                f"{chrom_data['coverage_ctrl'].median():>13.0f} "
                f"{pct_cov:>10.1f}% "
                f"{pct_conf:>11.1f}% "
                f"{bases_diff.sum():>11,}"
            )

    print()


# =============================================================================
# Section 4: Base composition at high-quality differing positions
#
# For experiments WITH candidates, show what types of mutations were found.
# For experiments WITHOUT candidates, show how close they came - positions
# where bases differ but failed coverage or confidence filters.
# =============================================================================

print("=" * 60)
print("SECTION 4: NEAR-MISS ANALYSIS (ZERO-HIT EXPERIMENTS)")
print("=" * 60)

for experiment_number in zero_hit_experiments:
    pickle_path = PICKLE_DIRECTORY / f"df_merged_exp{experiment_number}.pickle"
    merged_dataframe = pandas.read_pickle(str(pickle_path))

    # Positions where bases differ regardless of quality filters
    bases_differ = merged_dataframe["base_value_hit"] != merged_dataframe["base_value_ctrl"]
    differing_positions = merged_dataframe[bases_differ].copy()

    total_differing = len(differing_positions)

    print(f"\n  Experiment {experiment_number}: {total_differing:,} positions with different bases (before quality filters)")

    if total_differing == 0:
        print("    No base differences found at any quality level.")
        print("    Suppressor is very likely an indel or structural variant.")
        continue

    # Break down why they were filtered out
    has_cov_hit = (differing_positions["coverage_hit"] > MINIMUM_READ_COVERAGE).sum()
    has_cov_ctrl = (differing_positions["coverage_ctrl"] > MINIMUM_READ_COVERAGE).sum()
    has_cov_both = (
        (differing_positions["coverage_hit"] > MINIMUM_READ_COVERAGE)
        & (differing_positions["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
    ).sum()
    has_conf_hit = (
        (differing_positions["coverage_hit"] > MINIMUM_READ_COVERAGE)
        & (differing_positions["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        & (differing_positions["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT)
    ).sum()
    has_conf_both = (
        (differing_positions["coverage_hit"] > MINIMUM_READ_COVERAGE)
        & (differing_positions["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        & (differing_positions["confidence_hit"] > MINIMUM_CONFIDENCE_PERCENT)
        & (differing_positions["confidence_ctrl"] > MINIMUM_CONFIDENCE_PERCENT)
    ).sum()

    print(f"    Differing positions:                  {total_differing:>8,}")
    print(f"    With coverage > {MINIMUM_READ_COVERAGE} (hit):              {has_cov_hit:>8,}")
    print(f"    With coverage > {MINIMUM_READ_COVERAGE} (ctrl):             {has_cov_ctrl:>8,}")
    print(f"    With coverage > {MINIMUM_READ_COVERAGE} (both):             {has_cov_both:>8,}")
    print(f"    + confidence > {MINIMUM_CONFIDENCE_PERCENT}% (hit):            {has_conf_hit:>8,}")
    print(f"    + confidence > {MINIMUM_CONFIDENCE_PERCENT}% (both):           {has_conf_both:>8,}")

    if has_conf_both == 0 and has_cov_both > 0:
        # Positions differ and have coverage, but confidence is low
        # This suggests mixed reads - possible heterozygous indels causing misalignment
        low_conf_diffs = differing_positions[
            (differing_positions["coverage_hit"] > MINIMUM_READ_COVERAGE)
            & (differing_positions["coverage_ctrl"] > MINIMUM_READ_COVERAGE)
        ]
        median_conf_hit = low_conf_diffs["confidence_hit"].median()
        median_conf_ctrl = low_conf_diffs["confidence_ctrl"].median()
        print(
            f"\n    Covered but low-confidence positions suggest mixed base calls.\n"
            f"    Median confidence at differing, covered positions: "
            f"hit={median_conf_hit:.0f}%, ctrl={median_conf_ctrl:.0f}%\n"
            f"    This pattern is consistent with indels causing read misalignment."
        )

    if total_differing > 0 and has_cov_both == 0:
        print(
            f"\n    Base differences exist but all at low coverage.\n"
            f"    Sequencing depth may be insufficient for this experiment."
        )

print()


# =============================================================================
# Save diagnostic summary to CSV
# =============================================================================

print("=" * 60)
print("SAVE DIAGNOSTIC SUMMARY")
print("=" * 60)

funnel_dataframe = pandas.DataFrame(funnel_rows)
diagnostic_dataframe = pandas.merge(
    left=summary_dataframe,
    right=funnel_dataframe,
    on="experiment",
    how="inner",
    suffixes=("_summary", "_funnel"),
)

# Resolve the duplicate total_positions columns from the merge
if "total_positions_summary" in diagnostic_dataframe.columns:
    diagnostic_dataframe = diagnostic_dataframe.rename(
        columns={"total_positions_summary": "total_positions"}
    )
    diagnostic_dataframe = diagnostic_dataframe.drop(
        columns=["total_positions_funnel"], errors="ignore"
    )

output_csv_path = OUTPUT_DIRECTORY / OUTPUT_CSV_FILENAME
diagnostic_dataframe.to_csv(str(output_csv_path), index=False)
print(f"  Saved diagnostic summary to: {output_csv_path}")

print()
elapsed_seconds = time.perf_counter() - start_time
print(f"Total runtime: {elapsed_seconds:.1f} seconds")
print("Done.")
