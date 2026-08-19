# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "scipy",
#     "matplotlib",
# ]
# ///
r"""
Stage 3a of the manual gel path: read the ImageJ export (manual_lane_profiles.csv),
assert the invariants that must hold for a well-formed export, add the same
detrended_value that measure_gel.py computes per lane, and plot every lane so the
peaks can be checked by eye. This stage calls no bands and integrates nothing; it
exists to make the export trustworthy and legible before anything is built on it.

Reuses only what the pipeline already depends on (numpy, scipy, matplotlib). Every
non-standard-library call is written at its call site with the full namespace, per
the project convention (numpy.x, scipy.ndimage.x, matplotlib.pyplot.x).

Outputs, written next to the input CSV so they travel with the analysis and are
checked every run rather than discarded:
    profile_checks.json            the invariant checks and their verdicts
    lane_profiles_with_detrend.csv the input plus a detrended_value column
    lane_profiles_grid.png         one small-multiple panel per lane

Usage:
    uv run analyze_profiles.py '/mnt/c/.../manual_lane_profiles.csv'
"""

import argparse
import csv
import datetime
import json
import math
import os
import pathlib
import sys

import numpy
import scipy.ndimage
import matplotlib

matplotlib.use("Agg")  # no display in WSL; write files only
import matplotlib.pyplot

# =============================================================================
# Configuration. Named so a reader sees the constraint, not a bare literal.
# =============================================================================

# 16-bit container ceiling. raw_value is a sum of (pixel - background) over a lane
# width, so its own ceiling is width * (CONTAINER_MAXIMUM - background); a value
# above that is impossible and means a wrong background, orientation, or overflow.
CONTAINER_MAXIMUM_VALUE_16_BIT = 65535

# Plate background sanity band for these scanners. The two known files sit near
# 1067 (fluorescence) and 441 (phosphor); outside this range points at a bad crop
# or a miscalibrated container, so it is a warning rather than a hard stop.
PLATE_BACKGROUND_PLAUSIBLE_MINIMUM = 50.0
PLATE_BACKGROUND_PLAUSIBLE_MAXIMUM = 8000.0

# A lane whose peak barely rises above its own median carries no measurable band.
# Flagged (not failed): empty and ladder lanes legitimately look like this.
EMPTY_LANE_PEAK_OVER_MEDIAN_RATIO = 3.0

# Detrend window, copied verbatim from measure_gel.py's per-lane migration profile
# so detrended_value here is identical to what the automated pipeline would write.
DETREND_WINDOW_PROFILE_FRACTION = 2.0
DETREND_WINDOW_MINIMUM_PIXELS = 3

CHECKS_JSON_FILENAME = "profile_checks.json"
DETRENDED_CSV_FILENAME = "lane_profiles_with_detrend.csv"
GRID_PLOT_FILENAME = "lane_profiles_grid.png"
FIGURE_DOTS_PER_INCH = 150

# =============================================================================
# The two permitted helpers, matching the example scripts' logging contract.
# =============================================================================

ACCUMULATED_RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    ACCUMULATED_RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


# =============================================================================
# Arguments and input path
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="analyze_profiles.py",
    description="Validate and detrend an ImageJ lane-profile export; plot each lane.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "input_profile_csv_path", help="Path to manual_lane_profiles.csv"
)
parsed_arguments = argument_parser.parse_args()

input_profile_csv_path = pathlib.Path(
    os.path.abspath(parsed_arguments.input_profile_csv_path)
)
if not input_profile_csv_path.is_file():
    die("input", "not a file: " + str(input_profile_csv_path))
output_directory_path = input_profile_csv_path.parent
emit_message("input", "reading " + str(input_profile_csv_path))

# =============================================================================
# Read the CSV. Flat: one pass into typed column lists keyed by lane.
# =============================================================================

EXPECTED_COLUMNS = [
    "lane_index",
    "drawn_order",
    "roi_name",
    "lane_detection_status",
    "prediction_span",
    "comb_well_count",
    "roi_x",
    "roi_y",
    "roi_w",
    "roi_h",
    "plate_background_median",
    "migration_position_pixels",
    "migration_position_millimetres",
    "raw_value",
]

with open(input_profile_csv_path, newline="") as input_file_handle:
    csv_reader = csv.DictReader(input_file_handle)
    csv_header_fields = csv_reader.fieldnames
    all_rows = list(csv_reader)

missing_columns = [
    name for name in EXPECTED_COLUMNS if name not in (csv_header_fields or [])
]
if missing_columns:
    die(
        "schema",
        "CSV is missing required columns: "
        + ", ".join(missing_columns)
        + ". Re-export with the current export_lane_profiles.ijm.",
    )
if len(all_rows) == 0:
    die("schema", "CSV has a header but no data rows.")

# Group rows by lane while preserving the per-lane position order as written.
rows_by_lane_index = {}
for row in all_rows:
    lane_index_value = int(row["lane_index"])
    rows_by_lane_index.setdefault(lane_index_value, []).append(row)
sorted_lane_indices = sorted(rows_by_lane_index.keys())

# Scalars that must be constant across the whole file; read from the first row and
# asserted below rather than trusted.
comb_well_count = int(all_rows[0]["comb_well_count"])
plate_background_median = float(all_rows[0]["plate_background_median"])
first_roi_width = int(all_rows[0]["roi_w"])

emit_message(
    "read",
    "%d rows across %d lanes; comb_well_count=%d, plate_background=%.1f"
    % (
        len(all_rows),
        len(sorted_lane_indices),
        comb_well_count,
        plate_background_median,
    ),
)

# =============================================================================
# Invariants. Each check appends a verdict; hard failures stop after the report is
# written, so a failure is always diagnosable from profile_checks.json.
# =============================================================================

check_records = []


def record_check(check_name, severity, passed, detail_text):
    # severity is "hard" (a failure invalidates the export) or "soft" (a warning).
    check_records.append(
        {
            "check_name": check_name,
            "severity": severity,
            "passed": bool(passed),
            "detail": detail_text,
        }
    )


# The theoretical per-row ceiling for raw_value given this background and width.
raw_value_ceiling = first_roi_width * (
    CONTAINER_MAXIMUM_VALUE_16_BIT - plate_background_median
)

all_raw_values = numpy.array([float(row["raw_value"]) for row in all_rows], dtype=float)

record_check(
    "raw_value_non_negative",
    "hard",
    bool(numpy.all(all_raw_values >= 0.0)),
    "min raw_value = %.3f" % float(all_raw_values.min()),
)
record_check(
    "raw_value_under_ceiling",
    "hard",
    bool(numpy.all(all_raw_values <= raw_value_ceiling + 1.0)),
    "max raw_value = %.1f, ceiling width*(65535-bg) = %.1f"
    % (float(all_raw_values.max()), raw_value_ceiling),
)

# Positions must be a contiguous 0..n-1 run within every lane, or the migration
# axis is scrambled and no downstream position mapping can be trusted.
positions_contiguous = True
for lane_index_value in sorted_lane_indices:
    lane_positions = [
        int(r["migration_position_pixels"])
        for r in rows_by_lane_index[lane_index_value]
    ]
    if lane_positions != list(range(len(lane_positions))):
        positions_contiguous = False
        break
record_check(
    "positions_contiguous_per_lane",
    "hard",
    positions_contiguous,
    "each lane's migration_position_pixels runs 0..n-1",
)

# comb_well_count and plate_background_median must be single-valued across the file.
comb_values_distinct = {r["comb_well_count"] for r in all_rows}
background_values_distinct = {r["plate_background_median"] for r in all_rows}
record_check(
    "comb_well_count_constant",
    "hard",
    len(comb_values_distinct) == 1,
    "distinct comb_well_count values: " + ", ".join(sorted(comb_values_distinct)),
)
record_check(
    "plate_background_constant",
    "hard",
    len(background_values_distinct) == 1,
    "distinct plate_background_median values: " + str(len(background_values_distinct)),
)

# Lanes drawn must equal the comb size (the macro asserts this too; re-checked here
# because the CSV may outlive the macro run that made it).
record_check(
    "lane_count_matches_comb",
    "hard",
    len(sorted_lane_indices) == comb_well_count,
    "lanes=%d, comb_well_count=%d" % (len(sorted_lane_indices), comb_well_count),
)

# mm-per-pixel implied by the data must be one constant ratio.
implied_ratios = set()
for row in all_rows:
    position_pixels = int(row["migration_position_pixels"])
    if position_pixels > 0:
        implied_ratios.add(
            round(float(row["migration_position_millimetres"]) / position_pixels, 6)
        )
record_check(
    "millimetres_per_pixel_constant",
    "hard",
    len(implied_ratios) <= 1,
    "distinct mm/pixel ratios: " + str(sorted(implied_ratios)),
)

# --- soft checks ---
record_check(
    "plate_background_plausible",
    "soft",
    PLATE_BACKGROUND_PLAUSIBLE_MINIMUM
    <= plate_background_median
    <= PLATE_BACKGROUND_PLAUSIBLE_MAXIMUM,
    "plate_background_median = %.1f, plausible band [%.0f, %.0f]"
    % (
        plate_background_median,
        PLATE_BACKGROUND_PLAUSIBLE_MINIMUM,
        PLATE_BACKGROUND_PLAUSIBLE_MAXIMUM,
    ),
)

# Per-lane "carries signal": peak well above the lane's own median. Empty and
# ladder lanes are expected to fail this; it is information, not an error.
lanes_without_signal = []
for lane_index_value in sorted_lane_indices:
    lane_raw = numpy.array(
        [float(r["raw_value"]) for r in rows_by_lane_index[lane_index_value]],
        dtype=float,
    )
    lane_median = float(numpy.median(lane_raw))
    lane_peak = float(lane_raw.max())
    carries_signal = (
        lane_median > 0.0
        and lane_peak >= EMPTY_LANE_PEAK_OVER_MEDIAN_RATIO * lane_median
    )
    if not carries_signal:
        lanes_without_signal.append(lane_index_value)
record_check(
    "all_lanes_carry_signal",
    "soft",
    len(lanes_without_signal) == 0,
    "lanes below peak/median ratio %.1f (likely empty or ladder): %s"
    % (
        EMPTY_LANE_PEAK_OVER_MEDIAN_RATIO,
        ", ".join(str(i) for i in lanes_without_signal)
        if lanes_without_signal
        else "none",
    ),
)

hard_failures = [
    c for c in check_records if c["severity"] == "hard" and not c["passed"]
]
soft_warnings = [
    c for c in check_records if c["severity"] == "soft" and not c["passed"]
]

checks_report = {
    "schema_version": "profile_checks_1",
    "generated_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "input_profile_csv": str(input_profile_csv_path),
    "lane_count": len(sorted_lane_indices),
    "comb_well_count": comb_well_count,
    "plate_background_median": plate_background_median,
    "raw_value_ceiling": raw_value_ceiling,
    "overall_status": "fail"
    if hard_failures
    else ("warning" if soft_warnings else "pass"),
    "hard_failure_count": len(hard_failures),
    "warning_count": len(soft_warnings),
    "checks": check_records,
}
checks_json_path = output_directory_path / CHECKS_JSON_FILENAME
checks_json_path.write_text(json.dumps(checks_report, indent=2) + "\n")
emit_message(
    "checks",
    "wrote " + str(checks_json_path) + "; overall " + checks_report["overall_status"],
)

for warning_check in soft_warnings:
    emit_message(
        "checks",
        "WARNING " + warning_check["check_name"] + ": " + warning_check["detail"],
    )
if hard_failures:
    for failed_check in hard_failures:
        emit_message(
            "checks",
            "FAIL " + failed_check["check_name"] + ": " + failed_check["detail"],
        )
    die(
        "checks",
        "%d hard invariant(s) failed; see %s. Nothing downstream should run."
        % (len(hard_failures), str(checks_json_path)),
    )

# =============================================================================
# Detrend, identical to measure_gel.py's per-lane migration profile, and write the
# augmented CSV.
# =============================================================================

detrended_value_by_row_id = {}
for lane_index_value in sorted_lane_indices:
    lane_rows = rows_by_lane_index[lane_index_value]
    raw_lane_profile = numpy.array(
        [float(r["raw_value"]) for r in lane_rows], dtype=float
    )
    profile_length = raw_lane_profile.size
    profile_smoothing_window = max(
        DETREND_WINDOW_MINIMUM_PIXELS,
        int(
            round(
                DETREND_WINDOW_PROFILE_FRACTION
                * profile_length
                / float(comb_well_count)
            )
        ),
    )
    profile_baseline = scipy.ndimage.uniform_filter1d(
        raw_lane_profile, size=profile_smoothing_window, mode="nearest"
    )
    detrended_lane_profile = raw_lane_profile - profile_baseline
    detrended_lane_profile = detrended_lane_profile - detrended_lane_profile.mean()
    for position_index in range(profile_length):
        detrended_value_by_row_id[id(lane_rows[position_index])] = float(
            detrended_lane_profile[position_index]
        )

detrended_csv_path = output_directory_path / DETRENDED_CSV_FILENAME
augmented_header = list(csv_header_fields) + ["detrended_value"]
with open(detrended_csv_path, "w", newline="") as output_file_handle:
    csv_writer = csv.DictWriter(output_file_handle, fieldnames=augmented_header)
    csv_writer.writeheader()
    for row in all_rows:
        output_row = dict(row)
        output_row["detrended_value"] = "%.6f" % detrended_value_by_row_id[id(row)]
        csv_writer.writerow(output_row)
emit_message(
    "detrend",
    "wrote "
    + str(detrended_csv_path)
    + " (detrend window per lane = max(%d, round(%.1f*n/%d)))"
    % (DETREND_WINDOW_MINIMUM_PIXELS, DETREND_WINDOW_PROFILE_FRACTION, comb_well_count),
)

# =============================================================================
# Small-multiples plot: one panel per lane, raw and detrended, on a shared x in mm.
# =============================================================================

lane_panel_count = len(sorted_lane_indices)
grid_column_count = int(math.ceil(math.sqrt(lane_panel_count)))
grid_row_count = int(math.ceil(lane_panel_count / float(grid_column_count)))
grid_figure, grid_axes = matplotlib.pyplot.subplots(
    grid_row_count,
    grid_column_count,
    figsize=(grid_column_count * 2.6, grid_row_count * 2.0),
    squeeze=False,
)
for panel_index in range(grid_row_count * grid_column_count):
    panel_axis = grid_axes[panel_index // grid_column_count][
        panel_index % grid_column_count
    ]
    if panel_index >= lane_panel_count:
        panel_axis.axis("off")
        continue
    lane_index_value = sorted_lane_indices[panel_index]
    lane_rows = rows_by_lane_index[lane_index_value]
    position_millimetres = numpy.array(
        [float(r["migration_position_millimetres"]) for r in lane_rows]
    )
    raw_lane_profile = numpy.array([float(r["raw_value"]) for r in lane_rows])
    detrended_lane_profile = numpy.array(
        [detrended_value_by_row_id[id(r)] for r in lane_rows]
    )
    panel_axis.plot(
        position_millimetres, raw_lane_profile, color="0.6", linewidth=0.6, label="raw"
    )
    panel_axis.plot(
        position_millimetres,
        detrended_lane_profile,
        color="black",
        linewidth=0.8,
        label="detrended",
    )
    panel_axis.axhline(0.0, color="red", linewidth=0.4)
    lane_note = "" if lane_index_value not in lanes_without_signal else " (low)"
    panel_axis.set_title("lane %d%s" % (lane_index_value, lane_note), fontsize="small")
    panel_axis.tick_params(labelsize=6)
grid_figure.supxlabel("migration position (mm)", fontsize="small")
grid_figure.suptitle(
    "Lane profiles: raw (grey) and detrended (black). Check peaks by eye.",
    fontsize="small",
)
grid_figure.tight_layout()
grid_plot_path = output_directory_path / GRID_PLOT_FILENAME
grid_figure.savefig(grid_plot_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(grid_figure)
emit_message("plot", "wrote " + str(grid_plot_path))

emit_message(
    "done",
    "overall %s; %d warning(s). Inspect %s and %s before continuing."
    % (
        checks_report["overall_status"],
        len(soft_warnings),
        GRID_PLOT_FILENAME,
        CHECKS_JSON_FILENAME,
    ),
)
sys.exit(0)
