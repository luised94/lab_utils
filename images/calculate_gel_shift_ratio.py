# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Turn a gel shift gel's two region extractions into a per-lane bound fraction.

    uv run gel_shift_ratio.py '<gel.tif>'
    # or  uv run gel_shift_ratio.py '<stem>_gel_analysis'

A gel shift lane carries signal in two migration bands: a slower (higher on the
gel, smaller migration millimetres) BOUND band and a faster (lower, larger
millimetres) FREE band. You extract each band as its own region with
extract_lane_values.py --region, producing two extract_region_*.csv in the gel
analysis directory. This script pairs those two by lane and computes, per lane:

    bound_fraction = bound / (bound + free)

That fraction, not a plain bound/free ratio, is the reported quantity because it
stays bounded in [0, 1] when either band goes to zero: an all-bound lane reads
1.0, an all-free lane reads 0.0, and only a lane with NO signal in either band
(the denominator vanishes) is undefined -- that lane is emitted blank and flagged,
never invented as a 0 or a 1.

TOP/BOTTOM is assigned by migration position, not by any label: the region with
the smaller region_start_millimetres is the bound (top) band. The two regions must
be different windows; two identical windows are a hard failure, since then there is
no top and bottom to pair.

Identity and the biological factor columns come from sample_sheet.csv (the extract
CSVs carry lane identity but not the per-experiment factors), joined by lane_index,
and are carried through verbatim so the downstream stats/plot step can group by any
of them.

Output: gel_shift_ratio_<bound-window>_over_<free-window>.csv in the gel analysis
directory, one row per lane the sample sheet describes, plus a _checks.json with
the run's verdicts. Flat and procedural; emit_message / die / record_check are the
only helpers, duplicated verbatim per the house convention (no cross-script import).
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys

GEL_ANALYSIS_DIRECTORY_SUFFIX = "_gel_analysis"
SAMPLE_SHEET_FILENAME = "sample_sheet.csv"
EXTRACT_REGION_GLOB = "extract_region_*.csv"
# A lane is undefined only when NEITHER band carries signal: the extractor writes a
# blank value for a lane empty in a window, and a present-but-zero value sums the
# same. The fraction is computed whenever the summed signal exceeds this floor.
# Kept at exactly zero (no magic threshold): a lane is excluded only when both
# bands are blank/zero. Raise this if a real signal floor is ever wanted.
BOUND_PLUS_FREE_SIGNAL_FLOOR = 0.0

ACCUMULATED_RUN_LOG_LINES = []
accumulated_check_records = []
# output_stem is referenced by record_check on a hard failure that fires before the
# stem can be derived from the two windows, so it has a placeholder until then.
output_stem = "gel_shift_ratio"


def emit_message(source_tag, message_text):
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    ACCUMULATED_RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


# record_check appends a verdict and, on a HARD failure, writes the checks file
# before exiting so a failed run is diagnosable from the artifact and not only the
# terminal. gel_analysis_directory and output_stem it needs are set before the
# checks that can hard-fail through it run.
def record_check(check_name, severity, passed, detail_text):
    accumulated_check_records.append(
        {
            "check": check_name,
            "severity": severity,
            "passed": bool(passed),
            "detail": detail_text,
        }
    )
    if severity == "hard" and not passed:
        hard_failure_checks_path = gel_analysis_directory / (
            output_stem + "_checks.json"
        )
        hard_failure_checks_path.write_text(
            json.dumps(accumulated_check_records, indent=2), encoding="utf-8"
        )
        emit_message("check", "HARD FAIL: " + check_name + " -- " + detail_text)
        sys.exit(1)
    emit_message(
        "check",
        ("PASS " if passed else "soft ") + check_name + " -- " + detail_text,
    )


# =============================================================================
# Arguments and TIFF-addressed path resolution (the pipeline's single scheme: the
# gel .tif is the entry point, and <stem>_gel_analysis beside it holds the data).
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="gel_shift_ratio.py",
    description="Pair a gel's two region extractions into a per-lane bound fraction.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "gel_path",
    help="the gel analysis directory, or the gel .tif beside it (stem maps to "
    "<stem>_gel_analysis)",
)
parsed_arguments = argument_parser.parse_args()

given_path = pathlib.Path(parsed_arguments.gel_path)
if given_path.is_dir():
    gel_analysis_directory = given_path
elif given_path.is_file():
    if given_path.suffix.lower() != ".tif":
        die(
            "input",
            "file entry point must be the gel .tif (stem maps to <stem>"
            + GEL_ANALYSIS_DIRECTORY_SUFFIX
            + "); got: "
            + str(given_path),
        )
    gel_analysis_directory = given_path.parent / (
        given_path.stem + GEL_ANALYSIS_DIRECTORY_SUFFIX
    )
    if not gel_analysis_directory.is_dir():
        die(
            "input",
            "derived analysis directory does not exist: "
            + str(gel_analysis_directory)
            + " (expected beside the TIFF; run the ImageJ export first)",
        )
else:
    die("input", "not a file or directory: " + str(given_path))
emit_message("input", "gel analysis directory: " + str(gel_analysis_directory))

sample_sheet_path = gel_analysis_directory / SAMPLE_SHEET_FILENAME
if not sample_sheet_path.is_file():
    die(
        "input",
        "sample sheet not found at "
        + str(sample_sheet_path)
        + " (author it in the gel analysis directory)",
    )

# =============================================================================
# Find EXACTLY two region extractions. Two is the whole contract: a gel shift lane
# has a bound and a free band, extracted as two windows. Zero/one means the
# extractions have not been produced yet; three or more means a stale extract is
# present and the pairing would be ambiguous -- refuse rather than guess which two.
# =============================================================================

region_extract_paths = sorted(gel_analysis_directory.glob(EXTRACT_REGION_GLOB))
if len(region_extract_paths) != 2:
    die(
        "input",
        "expected exactly two %s (one bound band, one free band); found %d: %s"
        % (
            EXTRACT_REGION_GLOB,
            len(region_extract_paths),
            ", ".join(path.name for path in region_extract_paths) or "none",
        ),
    )


# Read one region extraction: its per-lane value keyed by lane_index, and the
# region window (constant across the file) used to order bound vs free. utf-8-sig
# strips an Excel BOM; every cell is trimmed so a hand-touched CSV cannot break the
# lane join or the numeric parse invisibly. The value field is blank for a lane with
# no signal in the window; blank parses to None here and is treated as zero signal
# only where the fraction arithmetic makes that explicit.
def read_region_extract(region_extract_path):
    value_by_lane_index = {}
    region_start_millimetres_seen = set()
    region_end_millimetres_seen = set()
    with region_extract_path.open(newline="", encoding="utf-8-sig") as region_handle:
        region_reader = csv.DictReader(region_handle)
        for raw_row in region_reader:
            row = {
                column_name: (
                    cell_value.strip() if isinstance(cell_value, str) else cell_value
                )
                for column_name, cell_value in raw_row.items()
            }
            lane_index = int(row["lane_index"])
            value_cell = row.get("value", "")
            value_by_lane_index[lane_index] = (
                None if value_cell == "" else float(value_cell)
            )
            region_start_millimetres_seen.add(float(row["region_start_millimetres"]))
            region_end_millimetres_seen.add(float(row["region_end_millimetres"]))
    # One extraction is one window, so the start/end are constant; if a file somehow
    # mixes windows the pairing has no single window to name, so refuse it.
    if len(region_start_millimetres_seen) != 1 or len(region_end_millimetres_seen) != 1:
        die(
            "input",
            "extract has more than one region window, cannot pair: "
            + str(region_extract_path),
        )
    return (
        value_by_lane_index,
        next(iter(region_start_millimetres_seen)),
        next(iter(region_end_millimetres_seen)),
    )


first_values, first_start_mm, first_end_mm = read_region_extract(
    region_extract_paths[0]
)
second_values, second_start_mm, second_end_mm = read_region_extract(
    region_extract_paths[1]
)

# Bound is the band higher on the gel: the smaller region_start_millimetres. Equal
# windows have no top/bottom to distinguish, so that is a hard failure.
if first_start_mm == second_start_mm:
    record_check(
        "the two region windows differ",
        "hard",
        False,
        "both regions start at %g mm; bound and free must be different windows"
        % first_start_mm,
    )
if first_start_mm < second_start_mm:
    bound_values, bound_start_mm, bound_end_mm = (
        first_values,
        first_start_mm,
        first_end_mm,
    )
    free_values, free_start_mm, free_end_mm = (
        second_values,
        second_start_mm,
        second_end_mm,
    )
else:
    bound_values, bound_start_mm, bound_end_mm = (
        second_values,
        second_start_mm,
        second_end_mm,
    )
    free_values, free_start_mm, free_end_mm = (
        first_values,
        first_start_mm,
        first_end_mm,
    )
emit_message(
    "input",
    "bound band [%g, %g] mm; free band [%g, %g] mm"
    % (bound_start_mm, bound_end_mm, free_start_mm, free_end_mm),
)

# The output name records both windows and their roles, so a directory that later
# holds several ratio files stays self-describing.
output_stem = "gel_shift_ratio_%g-%gmm_bound_over_%g-%gmm_free" % (
    bound_start_mm,
    bound_end_mm,
    free_start_mm,
    free_end_mm,
)

# The two extractions are the same gel, so their lane sets must match; a mismatch
# means they were not produced from the same profile and cannot be paired.
bound_lane_indices = set(bound_values)
free_lane_indices = set(free_values)
record_check(
    "the two extractions cover the same lanes",
    "hard",
    bound_lane_indices == free_lane_indices,
    "bound lanes %s vs free lanes %s"
    % (sorted(bound_lane_indices), sorted(free_lane_indices)),
)

# =============================================================================
# Read the sample sheet for identity and the per-experiment factor columns, keyed
# by lane_index. Every column is carried through verbatim (minus the ones this
# script computes) so the stats/plot step can group by any factor.
# =============================================================================

with sample_sheet_path.open(newline="", encoding="utf-8-sig") as sample_sheet_handle:
    sample_sheet_reader = csv.DictReader(sample_sheet_handle)
    sample_sheet_column_names = sample_sheet_reader.fieldnames or []
    sample_sheet_row_by_lane_index = {}
    for raw_row in sample_sheet_reader:
        row = {
            column_name: (
                cell_value.strip() if isinstance(cell_value, str) else cell_value
            )
            for column_name, cell_value in raw_row.items()
        }
        lane_index_cell = row.get("lane_index", "")
        if lane_index_cell == "":
            continue  # a trailing blank template row carries no lane
        sample_sheet_row_by_lane_index[int(lane_index_cell)] = row

# Columns this script writes itself; a sample sheet must not silently override them.
COMPUTED_COLUMN_NAMES = [
    "bound_value",
    "free_value",
    "bound_fraction",
    "bound_region_start_millimetres",
    "bound_region_end_millimetres",
    "free_region_start_millimetres",
    "free_region_end_millimetres",
]
carried_sample_sheet_columns = [
    column_name
    for column_name in sample_sheet_column_names
    if column_name not in COMPUTED_COLUMN_NAMES
]

# =============================================================================
# Compute the bound fraction per lane. Emit one row per lane the sample sheet
# describes, in lane order, carrying the sheet's columns plus the computed values.
# =============================================================================

output_rows = []
lanes_undefined_both_bands = []
lanes_missing_from_sheet = []
for lane_index in sorted(bound_lane_indices):
    sheet_row = sample_sheet_row_by_lane_index.get(lane_index)
    if sheet_row is None:
        # A lane measured in the extract but absent from the sheet (a ladder/empty
        # lane the sheet omits): it has no identity or factors to report, so it is
        # not an analyte row. Note it and skip.
        lanes_missing_from_sheet.append(lane_index)
        continue

    bound_value = bound_values[lane_index]
    free_value = free_values[lane_index]
    # A blank (no signal in the window) contributes zero to the sum; the fraction is
    # undefined only when BOTH are blank/zero, i.e. the lane has no signal at all.
    bound_signal = 0.0 if bound_value is None else bound_value
    free_signal = 0.0 if free_value is None else free_value
    total_signal = bound_signal + free_signal
    if total_signal <= BOUND_PLUS_FREE_SIGNAL_FLOOR:
        bound_fraction_cell = ""  # undefined; never invented as 0 or 1
        lanes_undefined_both_bands.append(lane_index)
    else:
        bound_fraction_cell = round(bound_signal / total_signal, 6)

    output_row = {
        column_name: sheet_row.get(column_name, "")
        for column_name in carried_sample_sheet_columns
    }
    output_row["bound_value"] = "" if bound_value is None else round(bound_value, 4)
    output_row["free_value"] = "" if free_value is None else round(free_value, 4)
    output_row["bound_fraction"] = bound_fraction_cell
    output_row["bound_region_start_millimetres"] = bound_start_mm
    output_row["bound_region_end_millimetres"] = bound_end_mm
    output_row["free_region_start_millimetres"] = free_start_mm
    output_row["free_region_end_millimetres"] = free_end_mm
    output_rows.append(output_row)

record_check(
    "at least one lane produced a bound fraction",
    "hard",
    any(output_row["bound_fraction"] != "" for output_row in output_rows),
    "%d lane(s) with a defined bound fraction"
    % sum(1 for output_row in output_rows if output_row["bound_fraction"] != ""),
)
record_check(
    "no analyte lane is undefined in both bands",
    "soft",
    len(lanes_undefined_both_bands) == 0,
    "lanes with no signal in either band (emitted blank): %s"
    % (sorted(lanes_undefined_both_bands) if lanes_undefined_both_bands else "none"),
)
record_check(
    "every extracted lane is described by the sample sheet",
    "soft",
    len(lanes_missing_from_sheet) == 0,
    "extracted lanes absent from the sheet (skipped): %s"
    % (sorted(lanes_missing_from_sheet) if lanes_missing_from_sheet else "none"),
)

# =============================================================================
# Write the ratio CSV and the checks file beside the extractions.
# =============================================================================

output_column_names = carried_sample_sheet_columns + COMPUTED_COLUMN_NAMES
output_csv_path = gel_analysis_directory / (output_stem + ".csv")
with output_csv_path.open("w", newline="", encoding="utf-8") as output_csv_handle:
    csv_writer = csv.DictWriter(output_csv_handle, fieldnames=output_column_names)
    csv_writer.writeheader()
    csv_writer.writerows(output_rows)
emit_message("output", "wrote " + str(output_csv_path))

checks_path = gel_analysis_directory / (output_stem + "_checks.json")
checks_path.write_text(
    json.dumps(accumulated_check_records, indent=2), encoding="utf-8"
)
emit_message("output", "wrote " + str(checks_path))

# Console summary to stderr (stdout stays empty so the CSV path can be diffed).
defined_fractions = [
    float(output_row["bound_fraction"])
    for output_row in output_rows
    if output_row["bound_fraction"] != ""
]
if defined_fractions:
    emit_message(
        "summary",
        "n=%d bound_fraction min=%.4f max=%.4f"
        % (len(defined_fractions), min(defined_fractions), max(defined_fractions)),
    )
print(str(output_csv_path))
