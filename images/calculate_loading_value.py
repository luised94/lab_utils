# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Turn a loading gel's single region extraction into a blank-corrected per-lane
value, carrying every sample-sheet factor column through for downstream grouping.

    uv run calculate_loading_value.py '<gel_analysis_dir>/extract_region_<window>.csv'
    # or  uv run calculate_loading_value.py '<gel_analysis_dir>'   (one extract only)
    # or  uv run calculate_loading_value.py '<gel.tif>'            (one extract only)

A loading lane carries a single quantity: the summed signal over one migration
window (a band, or a group of bands, chosen when you draw the region with
extract_lane_values.py --region). Unlike gel shift there is no second band and no
top/total ratio to form here; the loading ratio is lane-to-REFERENCE-lane, and that
normalization is deliberately left to the aggregate/plot step (plot_loading.R),
which needs every gel's reference lane to divide against. This script therefore
emits the measured per-lane value plus the blank it removed, and carries the
analysis_role column so the reference lane is findable downstream.

Entry point is the extract CSV itself, by preference. A loading gel directory can
hold several extract_region_*.csv as you iterate on which window to sum, so globbing
the directory and demanding exactly one would fight that iteration. Naming the
extract picks the window unambiguously. A directory (or the gel .tif beside it) is
still accepted for convenience, but only when exactly one extract_region_*.csv is
present; two or more is a hard failure that asks you to name the one you mean. This
mirrors plot_single_experiment.py's front-end, not the gel-shift "exactly two" rule.

Blank correction follows the loading precedent (plot_single_experiment.py), not the
gel-shift decision to skip it: a region sum over a fixed window includes whatever
residual streak sits in that window even in an empty lane, so
    blank_baseline = mean extract value over lanes whose lane_content is empty
    value_corrected = max(0, value - blank_baseline)
Both the raw value and value_corrected are emitted, with blank_baseline on every
row, so the durable artifact stays honest about what was subtracted and the plot
step can normalize on either column.

Output: loading_value_<window>.csv in the gel analysis directory, one row per lane
the sample sheet describes, plus a _checks.json with the run's verdicts. Flat and
procedural; emit_message / die / record_check are the only helpers, duplicated
verbatim per the house convention (no cross-script import).
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
# The extract filename prefix; the output window token is the extract stem with this
# prefix removed, so the loading_value file names the same window it summed.
EXTRACT_FILENAME_PREFIX = "extract_"

# lane_content vocabulary from the sample sheet schema. Empty lanes define the blank;
# every described lane is emitted (the reference lane is a sample and must survive
# into the output so plot_loading.R can normalize against it).
LANE_CONTENT_EMPTY = "empty"

# A blank-subtracted value cannot go below zero: a lane reading under the blank is
# background, not negative protein. Clips are counted and surfaced as a soft check.
ZERO_FLOOR = 0.0

ACCUMULATED_RUN_LOG_LINES = []
accumulated_check_records = []
# output_stem is referenced by record_check on a hard failure that can fire before
# the stem is derived from the extract window, so it holds a placeholder until then.
output_stem = "loading_value"


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
# Arguments and extract-addressed path resolution. The extract CSV is the primary
# entry point (it names the window); a directory or the gel .tif beside it is
# accepted only when exactly one extract_region_*.csv is present.
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="calculate_loading_value.py",
    description="Blank-correct one gel's single region extraction into a per-lane "
    "loading value, carrying the sample-sheet factor columns through.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "gel_path",
    help="the extract_region_*.csv (preferred; names the window), or the gel "
    "analysis directory / the gel .tif beside it when exactly one extract is present",
)
parsed_arguments = argument_parser.parse_args()

given_path = pathlib.Path(parsed_arguments.gel_path)

# Resolve to a single extract_region_*.csv and its analysis directory. Three shapes:
# the extract CSV directly; a directory holding exactly one extract; the gel .tif
# whose sibling <stem>_gel_analysis holds exactly one extract. Ambiguity (2+
# extracts) hard-fails asking for the specific file, rather than guessing a window.
if given_path.is_file() and given_path.name.startswith("extract_region_"):
    region_extract_path = given_path
    gel_analysis_directory = given_path.parent
elif given_path.is_dir():
    gel_analysis_directory = given_path
    candidate_extract_paths = sorted(gel_analysis_directory.glob(EXTRACT_REGION_GLOB))
    if len(candidate_extract_paths) != 1:
        die(
            "input",
            "expected exactly one %s in %s to auto-select; found %d: %s -- name the "
            "specific extract CSV instead"
            % (
                EXTRACT_REGION_GLOB,
                str(gel_analysis_directory),
                len(candidate_extract_paths),
                ", ".join(path.name for path in candidate_extract_paths) or "none",
            ),
        )
    region_extract_path = candidate_extract_paths[0]
elif given_path.is_file():
    # A .tif entry point: derive <stem>_gel_analysis beside it, then require exactly
    # one extract inside. A non-.tif file that is not an extract cannot be resolved.
    if given_path.suffix.lower() != ".tif":
        die(
            "input",
            "file entry point must be an extract_region_*.csv or the gel .tif; got: "
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
    candidate_extract_paths = sorted(gel_analysis_directory.glob(EXTRACT_REGION_GLOB))
    if len(candidate_extract_paths) != 1:
        die(
            "input",
            "expected exactly one %s in %s to auto-select; found %d: %s -- name the "
            "specific extract CSV instead"
            % (
                EXTRACT_REGION_GLOB,
                str(gel_analysis_directory),
                len(candidate_extract_paths),
                ", ".join(path.name for path in candidate_extract_paths) or "none",
            ),
        )
    region_extract_path = candidate_extract_paths[0]
else:
    die("input", "not a file or directory: " + str(given_path))
emit_message("input", "gel analysis directory: " + str(gel_analysis_directory))
emit_message("input", "region extract: " + str(region_extract_path))

sample_sheet_path = gel_analysis_directory / SAMPLE_SHEET_FILENAME
if not sample_sheet_path.is_file():
    die(
        "input",
        "sample sheet not found at "
        + str(sample_sheet_path)
        + " (author it in the gel analysis directory)",
    )

# The output window token names the same window the extract summed. The extract stem
# is "extract_region_<window>"; strip the leading "extract_" so the loading_value
# file sits beside its extract and names the same selection.
extract_filename_stem = region_extract_path.stem
if not extract_filename_stem.startswith(EXTRACT_FILENAME_PREFIX):
    die(
        "input",
        "extract filename does not start with '"
        + EXTRACT_FILENAME_PREFIX
        + "': "
        + region_extract_path.name,
    )
selector = extract_filename_stem[len(EXTRACT_FILENAME_PREFIX) :]
output_stem = "loading_value_" + selector

# =============================================================================
# Read the single region extraction: its per-lane value keyed by lane_index.
# utf-8-sig strips an Excel BOM and every cell is trimmed so a hand-touched CSV
# cannot break the lane join or the numeric parse invisibly. The value field is
# blank for a lane with no signal in the window; blank parses to None here and is
# treated as zero signal only where the arithmetic makes that explicit.
# =============================================================================

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

# One extraction is one window, so the start/end are constant; a file that somehow
# mixes windows has no single window to name, so refuse it.
if len(region_start_millimetres_seen) != 1 or len(region_end_millimetres_seen) != 1:
    die(
        "input",
        "extract has more than one region window, cannot summarize: "
        + str(region_extract_path),
    )
region_start_millimetres = next(iter(region_start_millimetres_seen))
region_end_millimetres = next(iter(region_end_millimetres_seen))
emit_message(
    "input",
    "region window [%g, %g] mm" % (region_start_millimetres, region_end_millimetres),
)

# =============================================================================
# Read the sample sheet for identity, lane disposition, and the per-experiment
# factor columns, keyed by lane_index. Every column is carried through verbatim
# (minus the ones this script computes) so plot_loading.R can group by any factor.
# lane_content is needed here to identify the empty lanes that define the blank.
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

if "lane_content" not in sample_sheet_column_names:
    die(
        "input",
        "sample sheet has no lane_content column; cannot identify the empty lanes "
        "that define the blank baseline",
    )

# The two extraction/sheet lane sets must match; a mismatch means the extract was
# produced from a different profile than the sheet describes and cannot be joined.
extract_lane_indices = set(value_by_lane_index)
sheet_lane_indices = set(sample_sheet_row_by_lane_index)
record_check(
    "the extract and the sample sheet cover the same lanes",
    "hard",
    extract_lane_indices == sheet_lane_indices,
    "extract lanes %s vs sheet lanes %s"
    % (sorted(extract_lane_indices), sorted(sheet_lane_indices)),
)

# =============================================================================
# Blank baseline from the empty lanes, then blank-correct every lane. The blank is
# the residual local background/streak a region sum picks up even where nothing was
# loaded; it is measured on lane_content == empty lanes.
# =============================================================================

empty_lane_indices = sorted(
    lane_index
    for lane_index in value_by_lane_index
    if sample_sheet_row_by_lane_index[lane_index].get("lane_content", "")
    == LANE_CONTENT_EMPTY
)
empty_lane_values = [
    value_by_lane_index[lane_index]
    for lane_index in empty_lane_indices
    if value_by_lane_index[lane_index] is not None
]
if empty_lane_values:
    blank_baseline = sum(empty_lane_values) / len(empty_lane_values)
else:
    # No empty lane to measure the residual background: correction is a no-op. This
    # is legal (a gel with no blank lane) but worth surfacing, so it is a soft check.
    blank_baseline = ZERO_FLOOR
record_check(
    "at least one empty lane defines the blank baseline",
    "soft",
    len(empty_lane_values) > 0,
    "%d empty lane(s) with signal used; blank_baseline=%.4f"
    % (len(empty_lane_values), blank_baseline),
)

# Columns this script writes itself; a sample sheet must not silently override them.
# analysis_role is NOT computed here -- it is carried from the sheet so plot_loading.R
# can find the reference lane -- so it is deliberately absent from this list.
COMPUTED_COLUMN_NAMES = [
    "value",
    "blank_baseline",
    "value_corrected",
    "region_start_millimetres",
    "region_end_millimetres",
]
carried_sample_sheet_columns = [
    column_name
    for column_name in sample_sheet_column_names
    if column_name not in COMPUTED_COLUMN_NAMES
]

# =============================================================================
# Compute the blank-corrected value per lane. Emit one row per lane the sample sheet
# describes, in lane order, carrying the sheet's columns plus the computed values.
# Every described lane is emitted (including the reference lane and, unlike
# plot_single_experiment.py, the empty/ladder lanes) so the downstream step has the
# full gel: it needs the reference to normalize against and can filter the rest.
# =============================================================================

output_rows = []
clipped_to_zero_lane_count = 0
missing_value_lane_count = 0
for lane_index in sorted(extract_lane_indices):
    sheet_row = sample_sheet_row_by_lane_index[lane_index]
    raw_value = value_by_lane_index[lane_index]
    if raw_value is None:
        # A lane with no signal in the window: correct to zero and flag it, rather
        # than emitting a blank the downstream numeric join would trip on.
        missing_value_lane_count += 1
        raw_value_for_output = ""
        value_corrected = ZERO_FLOOR
    else:
        raw_value_for_output = round(raw_value, 4)
        value_corrected = raw_value - blank_baseline
        if value_corrected < ZERO_FLOOR:
            value_corrected = ZERO_FLOOR
            clipped_to_zero_lane_count += 1

    output_row = {
        column_name: sheet_row.get(column_name, "")
        for column_name in carried_sample_sheet_columns
    }
    output_row["value"] = raw_value_for_output
    output_row["blank_baseline"] = round(blank_baseline, 4)
    output_row["value_corrected"] = round(value_corrected, 4)
    output_row["region_start_millimetres"] = region_start_millimetres
    output_row["region_end_millimetres"] = region_end_millimetres
    output_rows.append(output_row)

record_check(
    "at least one lane produced a value",
    "hard",
    len(output_rows) > 0,
    "%d lane(s) emitted" % len(output_rows),
)
record_check(
    "no lane clipped to zero by blank subtraction",
    "soft",
    clipped_to_zero_lane_count == 0,
    "%d lane(s) fell below the blank and were clipped to zero"
    % clipped_to_zero_lane_count,
)
record_check(
    "every lane has a value in the window",
    "soft",
    missing_value_lane_count == 0,
    "%d lane(s) had a blank value and were corrected to zero"
    % missing_value_lane_count,
)

# The downstream reference normalization needs exactly one analysis_role == reference
# lane per gel; surface a wrong count here (soft: the calc still emits a valid CSV,
# but plot_loading.R will hard-fail, so warn early where it is cheap to notice).
if "analysis_role" in sample_sheet_column_names:
    reference_lane_count = sum(
        1
        for output_row in output_rows
        if output_row.get("analysis_role", "") == "reference"
    )
    record_check(
        "exactly one reference lane for downstream normalization",
        "soft",
        reference_lane_count == 1,
        "%d lane(s) with analysis_role=reference" % reference_lane_count,
    )
else:
    record_check(
        "sample sheet carries analysis_role for downstream normalization",
        "soft",
        False,
        "no analysis_role column; plot_loading.R cannot find the reference lane",
    )

# =============================================================================
# Write the loading value CSV and the checks file beside the extraction.
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
defined_values = [
    float(output_row["value_corrected"])
    for output_row in output_rows
    if output_row["value_corrected"] != ""
]
if defined_values:
    emit_message(
        "summary",
        "n=%d value_corrected min=%.4f max=%.4f blank=%.4f"
        % (len(defined_values), min(defined_values), max(defined_values), blank_baseline),
    )
print(str(output_csv_path))
