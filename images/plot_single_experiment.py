# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib>=3.11.1",
# ]
# ///
r"""
Turn one gel's per-lane extract into blank-corrected per-sample numbers and a bar
chart. This is the single-experiment step: it consumes an extract_region_*.csv or
extract_band_*.csv (written by extract_lane_values.py), subtracts the empty-lane
background that the global plate median leaves behind, drops the non-analyte lanes,
and writes exactly the samples an experiment reports, with the blank it removed
carried alongside every number.

What it does, in order:

  - Read the extract CSV for the per-lane raw value column, joined by lane_index to
    sample_sheet.csv for lane_content and condition_type (the extract does not read
    the sheet, so empty-vs-ladder-vs-sample is only known here).
  - blank_baseline = mean of the extract value over lanes whose lane_content is
    empty. This is the residual local background/streak the plate median cannot
    remove (see EXPLANATION_and_decisions.md, "Background and blank lanes").
  - value_corrected = value - blank_baseline, clipped at zero (a sample below the
    blank cannot carry negative protein; a clip is flagged as a soft check).
  - Exclude ladder and empty lanes; write one row and one bar per sample lane.

Addressing follows the rest of the pipeline: point at the gel analysis directory,
or at the extract CSV itself, and the standard filenames resolve. --extract-csv and
--sample-sheet override when the convention does not hold.

Outputs, written into the gel directory so they travel with the analysis and are
checked every run:
    single_experiment_<selector>.csv          one row per sample lane
    single_experiment_<selector>.png          one bar per sample lane (NOT viewed)
    single_experiment_<selector>_checks.json  inputs, blank, per-sample values, verdicts

<selector> is the extract filename with the leading "extract_" removed, so the
single-experiment files sit next to the extract they came from and name the same
selection (region+baseline or band+quantity). Flat and procedural; emit_message /
die / record_check are the only helpers, duplicated verbatim per the house
convention.
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys

# Never render interactively; this script only writes a PNG file and must never
# open a display. Set the non-interactive backend before importing pyplot.
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot

# Standard filenames resolved inside the gel analysis directory.
SAMPLE_SHEET_FILENAME = "sample_sheet.csv"

# The two extract families this step consumes; the selector is the filename with
# this prefix removed, so single_experiment files name the same selection.
EXTRACT_FILENAME_PREFIX = "extract_"
EXTRACT_REGION_GLOB = "extract_region_*.csv"
EXTRACT_BAND_GLOB = "extract_band_*.csv"

# lane_content vocabulary from the sample sheet schema. Only samples are plotted;
# empty lanes define the blank; ladder lanes are dropped.
LANE_CONTENT_SAMPLE = "sample"
LANE_CONTENT_EMPTY = "empty"
LANE_CONTENT_LADDER = "ladder"

# A blank-subtracted value cannot go below zero: a sample lane reading under the
# blank is background, not negative protein. Clips are counted and surfaced.
ZERO_FLOOR = 0.0

# Distinct condition_type values are coloured from this fixed palette (assigned in
# sorted order so the same condition gets the same colour run to run). ASCII hex.
CONDITION_COLOUR_PALETTE = [
    "#4C72B0",
    "#DD8452",
    "#55A868",
    "#C44E52",
    "#8172B3",
    "#937860",
    "#DA8BC3",
    "#8C8C8C",
]
FALLBACK_CONDITION_COLOUR = "#333333"

PLOT_FIGURE_WIDTH_INCHES = 11.0
PLOT_FIGURE_HEIGHT_INCHES = 6.0
PLOT_DOTS_PER_INCH = 150

ACCUMULATED_RUN_LOG_LINES = []


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
# terminal. output_directory, output_stem, gel_id and selector are set before any
# check runs.
accumulated_check_records = []


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
        hard_failure_checks_path = output_directory / (output_stem + "_checks.json")
        hard_failure_checks_path.write_text(
            json.dumps(
                {
                    "schema_version": "plot_single_experiment_1",
                    "generated_at": datetime.datetime.now()
                    .astimezone()
                    .isoformat(timespec="seconds"),
                    "gel_id": gel_id,
                    "selector": selector,
                    "checks": accumulated_check_records,
                    "failed": {"check": check_name, "detail": detail_text},
                },
                indent=2,
            )
            + "\n"
        )
        die("checks", "hard check failed: " + check_name + " (" + detail_text + ")")


# =============================================================================
# Arguments
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="plot_single_experiment.py",
    description="Blank-correct and plot one gel's per-lane extract values.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "gel_path",
    help="the gel analysis directory, or the extract CSV inside it",
)
argument_parser.add_argument(
    "--extract-csv",
    dest="extract_csv_override",
    default=None,
    help="path to the extract_region_*/extract_band_* CSV when auto-find is ambiguous",
)
argument_parser.add_argument(
    "--sample-sheet", dest="sample_sheet_override", default=None
)
parsed_arguments = argument_parser.parse_args()

# =============================================================================
# Resolve the extract CSV and, from it, the gel directory and selector. These
# failures happen before the checks file can be named (the selector IS the name),
# so they exit plainly rather than through record_check.
# =============================================================================

given_path = pathlib.Path(parsed_arguments.gel_path)
if parsed_arguments.extract_csv_override:
    extract_csv_path = pathlib.Path(parsed_arguments.extract_csv_override)
    if not extract_csv_path.is_file():
        print(
            "[input] ERROR: --extract-csv not found: " + str(extract_csv_path),
            file=sys.stderr,
        )
        sys.exit(1)
elif given_path.is_file() and (
    given_path.name.startswith("extract_region_")
    or given_path.name.startswith("extract_band_")
):
    # Pointed straight at an extract CSV: use it.
    extract_csv_path = given_path
else:
    # Pointed at the directory (or any other file in it): auto-find the extract.
    if given_path.is_dir():
        search_directory = given_path
    elif given_path.is_file():
        search_directory = given_path.parent
    else:
        print("[input] ERROR: path does not exist: " + str(given_path), file=sys.stderr)
        sys.exit(1)
    candidate_extract_paths = sorted(
        list(search_directory.glob(EXTRACT_REGION_GLOB))
        + list(search_directory.glob(EXTRACT_BAND_GLOB))
    )
    if len(candidate_extract_paths) == 0:
        print(
            "[input] ERROR: no extract_region_*.csv or extract_band_*.csv in "
            + str(search_directory),
            file=sys.stderr,
        )
        sys.exit(1)
    if len(candidate_extract_paths) > 1:
        # Do not guess which extraction the author meant; make them name it.
        print(
            "[input] ERROR: multiple extract CSVs in "
            + str(search_directory)
            + "; pass --extract-csv to choose one of: "
            + ", ".join(path.name for path in candidate_extract_paths),
            file=sys.stderr,
        )
        sys.exit(1)
    extract_csv_path = candidate_extract_paths[0]

output_directory = extract_csv_path.parent
extract_filename_stem = extract_csv_path.stem
if not extract_filename_stem.startswith(EXTRACT_FILENAME_PREFIX):
    print(
        "[input] ERROR: extract filename does not start with '"
        + EXTRACT_FILENAME_PREFIX
        + "': "
        + extract_csv_path.name,
        file=sys.stderr,
    )
    sys.exit(1)
selector = extract_filename_stem[len(EXTRACT_FILENAME_PREFIX) :]
output_stem = "single_experiment_" + selector

sample_sheet_path = (
    pathlib.Path(parsed_arguments.sample_sheet_override)
    if parsed_arguments.sample_sheet_override
    else output_directory / SAMPLE_SHEET_FILENAME
)

# gel_id is referenced by record_check's hard-failure writer; set a placeholder now
# and fill it from the extract's first row once read.
gel_id = ""

# =============================================================================
# Read the extract CSV. utf-8-sig strips an Excel BOM and every cell is trimmed,
# because a hand-exported CSV carries a BOM and stray whitespace that would break
# the exact lane joins and the value parsing invisibly.
# =============================================================================

with extract_csv_path.open(newline="", encoding="utf-8-sig") as extract_csv_handle:
    extract_rows = [
        {
            column_name: (
                cell_value.strip() if isinstance(cell_value, str) else cell_value
            )
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in csv.DictReader(extract_csv_handle)
    ]
if extract_rows:
    gel_id = extract_rows[0].get("gel_id", "")
record_check(
    "extract csv has rows",
    "hard",
    len(extract_rows) > 0,
    str(extract_csv_path),
)

# The extract value column is the per-lane number (raw_sum-derived in region mode,
# reported area/apex in band mode). It may be blank when a lane had no signal in a
# region window; blank parses to None and is excluded from the blank mean.
value_by_lane_index = {}
for extract_row in extract_rows:
    lane_index = int(extract_row["lane_index"])
    value_cell = extract_row.get("value", "")
    value_by_lane_index[lane_index] = (
        None if value_cell == "" or value_cell is None else float(value_cell)
    )
# Identity as the extract recorded it, kept only to cross-check the sheet below.
extract_identity_by_lane_index = {
    int(extract_row["lane_index"]): {
        "well_number": extract_row.get("well_number", ""),
        "sample_label": extract_row.get("sample_label", ""),
    }
    for extract_row in extract_rows
}

# =============================================================================
# Read the sample sheet. It is authoritative for lane disposition (lane_content),
# design class (condition_type) and identity (well_number, sample_label); the
# extract cannot know empty-vs-ladder because it never reads the sheet.
# =============================================================================

record_check(
    "sample sheet present",
    "hard",
    sample_sheet_path.is_file(),
    str(sample_sheet_path),
)
with sample_sheet_path.open(newline="", encoding="utf-8-sig") as sample_sheet_handle:
    sample_sheet_rows = [
        {
            column_name: (
                cell_value.strip() if isinstance(cell_value, str) else cell_value
            )
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in csv.DictReader(sample_sheet_handle)
    ]
record_check(
    "sample sheet has rows",
    "hard",
    len(sample_sheet_rows) > 0,
    str(sample_sheet_path),
)

sheet_by_lane_index = {}
for sample_sheet_row in sample_sheet_rows:
    lane_index = int(sample_sheet_row["lane_index"])
    sheet_by_lane_index[lane_index] = {
        "well_number": sample_sheet_row["well_number"],
        "lane_content": sample_sheet_row["lane_content"],
        "condition_type": sample_sheet_row["condition_type"],
        "sample_label": sample_sheet_row["sample_label"],
    }

# Every lane in the extract must be described by the sheet, or the blank mean and
# the sample/ladder/empty split are computed against an incomplete map.
extract_lanes_missing_from_sheet = sorted(
    lane_index
    for lane_index in value_by_lane_index
    if lane_index not in sheet_by_lane_index
)
record_check(
    "every extract lane is in the sample sheet",
    "hard",
    len(extract_lanes_missing_from_sheet) == 0,
    "lanes missing from sheet: "
    + (
        ", ".join(str(lane_index) for lane_index in extract_lanes_missing_from_sheet)
        if extract_lanes_missing_from_sheet
        else "none"
    ),
)

# =============================================================================
# Blank baseline from the empty lanes, then blank-correct the sample lanes.
# =============================================================================

empty_lane_indices = sorted(
    lane_index
    for lane_index in value_by_lane_index
    if sheet_by_lane_index[lane_index]["lane_content"] == LANE_CONTENT_EMPTY
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
    # is legal (a gel with no blank) but worth surfacing, so it is a soft check.
    blank_baseline = ZERO_FLOOR
record_check(
    "at least one empty lane defines the blank baseline",
    "soft",
    len(empty_lane_values) > 0,
    "%d empty lane(s) with signal used; blank_baseline=%.4f"
    % (len(empty_lane_values), blank_baseline),
)

sample_lane_indices = sorted(
    lane_index
    for lane_index in value_by_lane_index
    if sheet_by_lane_index[lane_index]["lane_content"] == LANE_CONTENT_SAMPLE
)
record_check(
    "at least one sample lane to plot",
    "hard",
    len(sample_lane_indices) > 0,
    "%d sample lane(s)" % len(sample_lane_indices),
)

output_rows = []
corrected_values = []
clipped_to_zero_lane_count = 0
identity_mismatch_lane_count = 0
missing_value_sample_lane_count = 0
for lane_index in sample_lane_indices:
    sheet_entry = sheet_by_lane_index[lane_index]
    raw_value = value_by_lane_index[lane_index]
    if raw_value is None:
        # A sample lane with no signal in the window: correct to zero and flag it,
        # rather than dropping the sample silently.
        missing_value_sample_lane_count += 1
        raw_value_for_output = ""
        value_corrected = ZERO_FLOOR
    else:
        raw_value_for_output = round(raw_value, 4)
        value_corrected = raw_value - blank_baseline
        if value_corrected < ZERO_FLOOR:
            value_corrected = ZERO_FLOOR
            clipped_to_zero_lane_count += 1
    corrected_values.append(value_corrected)

    # Cross-check the extract's recorded identity against the sheet; a disagreement
    # means a stale extract sitting next to a newer sheet (or vice versa).
    extract_identity = extract_identity_by_lane_index.get(
        lane_index, {"well_number": "", "sample_label": ""}
    )
    if (
        extract_identity["well_number"] != sheet_entry["well_number"]
        or extract_identity["sample_label"] != sheet_entry["sample_label"]
    ):
        identity_mismatch_lane_count += 1

    output_rows.append(
        {
            "lane_index": lane_index,
            "well_number": sheet_entry["well_number"],
            "sample_label": sheet_entry["sample_label"],
            "condition_type": sheet_entry["condition_type"],
            "value": raw_value_for_output,
            "blank_baseline": round(blank_baseline, 4),
            "value_corrected": round(value_corrected, 4),
        }
    )

record_check(
    "no sample lane clipped to zero by blank subtraction",
    "soft",
    clipped_to_zero_lane_count == 0,
    "%d of %d sample lanes fell below the blank and were clipped to zero"
    % (clipped_to_zero_lane_count, len(sample_lane_indices)),
)
record_check(
    "every sample lane has a value in the window",
    "soft",
    missing_value_sample_lane_count == 0,
    "%d of %d sample lanes had a blank value and were corrected to zero"
    % (missing_value_sample_lane_count, len(sample_lane_indices)),
)
record_check(
    "extract identity agrees with the sample sheet",
    "soft",
    identity_mismatch_lane_count == 0,
    "%d of %d sample lanes disagree on well_number/sample_label between extract and sheet"
    % (identity_mismatch_lane_count, len(sample_lane_indices)),
)

# =============================================================================
# Summary statistics over the blank-corrected values (the spread is what "is this
# consistent?" reduces to once the blank is removed).
# =============================================================================

ordered_values = sorted(corrected_values)
value_count = len(ordered_values)
mean_value = sum(ordered_values) / value_count
if value_count % 2:
    median_value = ordered_values[value_count // 2]
else:
    median_value = 0.5 * (
        ordered_values[value_count // 2 - 1] + ordered_values[value_count // 2]
    )
if value_count > 1 and mean_value != 0:
    sample_variance = sum((value - mean_value) ** 2 for value in ordered_values) / (
        value_count - 1
    )
    coefficient_of_variation = (sample_variance**0.5) / abs(mean_value)
else:
    coefficient_of_variation = None
value_summary = {
    "n": value_count,
    "minimum": min(ordered_values),
    "median": median_value,
    "maximum": max(ordered_values),
    "coefficient_of_variation": (
        round(coefficient_of_variation, 4)
        if coefficient_of_variation is not None
        else None
    ),
}

# =============================================================================
# Write the per-sample CSV.
# =============================================================================

output_csv_path = output_directory / (output_stem + ".csv")
output_column_names = list(output_rows[0].keys())
with output_csv_path.open("w", newline="", encoding="utf-8") as output_csv_handle:
    csv_writer = csv.DictWriter(output_csv_handle, fieldnames=output_column_names)
    csv_writer.writeheader()
    csv_writer.writerows(output_rows)
emit_message("output", "wrote " + str(output_csv_path))

# =============================================================================
# Write the bar chart. One bar per sample lane, height value_corrected, coloured by
# condition_type, x labelled by sample_label and well_number. The PNG is written
# and never opened by this script.
# =============================================================================

condition_types_in_order = sorted(
    {output_row["condition_type"] for output_row in output_rows}
)
colour_by_condition_type = {}
for condition_ordinal, condition_type_value in enumerate(condition_types_in_order):
    if condition_ordinal < len(CONDITION_COLOUR_PALETTE):
        colour_by_condition_type[condition_type_value] = CONDITION_COLOUR_PALETTE[
            condition_ordinal
        ]
    else:
        colour_by_condition_type[condition_type_value] = FALLBACK_CONDITION_COLOUR

bar_positions = list(range(len(output_rows)))
bar_heights = [output_row["value_corrected"] for output_row in output_rows]
bar_colours = [
    colour_by_condition_type[output_row["condition_type"]] for output_row in output_rows
]
bar_tick_labels = [
    output_row["sample_label"] + "\nwell " + str(output_row["well_number"])
    for output_row in output_rows
]

figure_object, axes_object = matplotlib.pyplot.subplots(
    figsize=(PLOT_FIGURE_WIDTH_INCHES, PLOT_FIGURE_HEIGHT_INCHES)
)
axes_object.bar(bar_positions, bar_heights, color=bar_colours)
axes_object.set_xticks(bar_positions)
axes_object.set_xticklabels(bar_tick_labels, rotation=45, ha="right", fontsize=8)
axes_object.set_ylabel("value_corrected (blank-subtracted)")
axes_object.set_title(gel_id + "\n" + selector)
legend_handles = [
    matplotlib.pyplot.Rectangle(
        (0, 0), 1, 1, color=colour_by_condition_type[condition_type_value]
    )
    for condition_type_value in condition_types_in_order
]
axes_object.legend(
    legend_handles, condition_types_in_order, title="condition_type", fontsize=8
)
figure_object.tight_layout()
output_png_path = output_directory / (output_stem + ".png")
figure_object.savefig(output_png_path, dpi=PLOT_DOTS_PER_INCH)
matplotlib.pyplot.close(figure_object)
emit_message("output", "wrote " + str(output_png_path))

# =============================================================================
# Write the checks JSON.
# =============================================================================

checks_path = output_directory / (output_stem + "_checks.json")
checks_path.write_text(
    json.dumps(
        {
            "schema_version": "plot_single_experiment_1",
            "generated_at": datetime.datetime.now()
            .astimezone()
            .isoformat(timespec="seconds"),
            "gel_id": gel_id,
            "selector": selector,
            "inputs": {
                "extract_csv": str(extract_csv_path),
                "sample_sheet": str(sample_sheet_path),
            },
            "provenance": {
                "blank_baseline": round(blank_baseline, 4),
                "empty_lane_indices": empty_lane_indices,
                "empty_lane_count_used": len(empty_lane_values),
                "sample_lane_count": len(sample_lane_indices),
                "clipped_to_zero_lane_count": clipped_to_zero_lane_count,
                "note": (
                    "blank_baseline is the mean extract value over empty lanes "
                    "(residual local background/streak the plate median leaves); "
                    "value_corrected = value - blank_baseline, clipped at zero. "
                    "ladder and empty lanes are excluded from the output"
                ),
            },
            "value_summary_over_sample_lanes": value_summary,
            "outputs": {
                "single_experiment_csv": str(output_csv_path),
                "single_experiment_png": str(output_png_path),
            },
            "checks": accumulated_check_records,
        },
        indent=2,
    )
    + "\n"
)
emit_message("output", "wrote " + str(checks_path))

# Console summary to stderr (stdout stays empty so paths can be diffed): the blank,
# the per-sample corrected values and their spread.
emit_message("summary", "gel " + gel_id)
emit_message("summary", "selector: " + selector)
emit_message(
    "summary",
    "blank_baseline=%.4f from %d empty lane(s) %s"
    % (blank_baseline, len(empty_lane_values), empty_lane_indices),
)
for output_row in output_rows:
    emit_message(
        "value",
        "lane %2s well %2s %-14s %-16s value %-14s corrected %s"
        % (
            output_row["lane_index"],
            output_row["well_number"],
            output_row["condition_type"],
            output_row["sample_label"],
            output_row["value"],
            output_row["value_corrected"],
        ),
    )
emit_message(
    "summary",
    "n=%s min=%s median=%s max=%s CV=%s"
    % (
        value_summary.get("n"),
        value_summary.get("minimum"),
        value_summary.get("median"),
        value_summary.get("maximum"),
        value_summary.get("coefficient_of_variation"),
    ),
)
for check_record in accumulated_check_records:
    if check_record["severity"] == "soft" and not check_record["passed"]:
        emit_message(
            "check", "WARNING " + check_record["check"] + ": " + check_record["detail"]
        )
