# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy",
#     "scipy",
#     "matplotlib",
# ]
# ///
r"""
The manual gel analysis: read the ImageJ export (manual_lane_profiles.csv) plus the
validated sample sheet, assert the export invariants, add the detrended_value that
measure_gel.py computes per lane, detect bands (per-lane and consensus), then
integrate every measured lane at the consensus windows to produce band_measurements
.csv, reusing measure_gel.py's valley-to-valley and rolling-ball baselines verbatim.
Reports both a region-net AREA and an apex HEIGHT per band, flags saturated bands
from the at_ceiling_count column, and scores the designated reference lane so a bad
divisor is caught before normalization. Every stage writes a checkable artifact.

Reuses only what the pipeline already depends on (numpy, scipy, matplotlib). Every
non-standard-library call is written at its call site with the full namespace, per
the project convention (numpy.x, scipy.ndimage.x, matplotlib.pyplot.x).

Outputs, written into the gel analysis directory so they travel with the analysis
and are checked every run rather than discarded:
    profile_checks.json            the export invariant checks and verdicts
    lane_profiles_with_detrend.csv the input plus a detrended_value column
    lane_profiles_grid.png         one small-multiple panel per lane
    band_candidates_per_lane.csv   per-lane detected peaks
    consensus_bands.csv            consensus band positions, windows, occupancy
    band_model_diagnostics.json    consensus-vs-per-lane agreement metrics
    band_detection_overlay.png     peaks, windows, apex exclusion, by lane
    band_measurements.csv          per (well, band) area, height, baselines, flags
    reference_quality.json         the reference lane score and per-input ranking
    band_measurement_overlay.png   windows, baselines, and both quantities by lane

Usage (directory-addressed, like the manifest and the sheet validator):
    uv run analyze_gel.py '<gel_analysis_dir>'
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
import scipy.signal
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

# --- Band detection, all values copied from measure_gel.py so the manual path
# calls bands with the same behaviour as the automated pipeline. ---
# Smoothing applied to the raw profile before peak finding (mm of migration).
BAND_MIGRATION_PROFILE_SMOOTHING_MILLIMETRES = 0.4
# A peak must rise this fraction of the lane's own maximum to be called.
BAND_PEAK_PROMINENCE_FRACTION = 0.12
# Two called peaks closer than this are one band (mm). Also the consensus cluster
# tolerance, tied together exactly as measure_gel.py ties them.
BAND_MINIMUM_SEPARATION_MILLIMETRES = 1.0
# Fixed integration-window width centred on each consensus band (mm), before it is
# clipped at the midpoint to each neighbour so windows never overlap.
BAND_REGION_WINDOW_WIDTH_MILLIMETRES = 2.0
# Top-of-lane zone excluded from DETECTION only (not from the data or the plot).
# Justified on this gel: the apex feature sits within ~1 mm of the ROI top in
# every lane and the nearest real band is >15 mm below it, so this cannot clip a
# real band. The plot hatches this zone so the claim is checkable every run.
APEX_EXCLUSION_MILLIMETRES = 5.0
# A consensus cluster is "well supported" if this fraction of lanes contributed a
# peak; used only to score whether consensus is trustworthy, never to drop bands.
CONSENSUS_WELL_SUPPORTED_LANE_FRACTION = 0.5

BAND_CANDIDATES_CSV_FILENAME = "band_candidates_per_lane.csv"
CONSENSUS_BANDS_CSV_FILENAME = "consensus_bands.csv"
BAND_MODEL_DIAGNOSTICS_FILENAME = "band_model_diagnostics.json"
BAND_DETECTION_PLOT_FILENAME = "band_detection_overlay.png"

# --- Integration, copied from measure_gel.py so the manual path integrates bands
# the same way as the automated pipeline. ---
# How far outside each window the valley-to-valley baseline searches for a flanking
# minimum (mm), so neighbouring bands sit on a consistent local baseline.
BASELINE_FLANK_SEARCH_MILLIMETRES = 1.5
# Width of the flat rolling-ball opening element (mm), the cross-check baseline.
# Must exceed the widest single band or the baseline rides up into the peak.
ROLLING_BALL_WIDTH_MILLIMETRES = 4.0
# A band with at least this many at-ceiling pixels anywhere in its window is
# saturated: its area is a floor (clipped pixels lost their true height), so the
# saturation basis wins the reported-value contract.
SATURATION_AT_CEILING_PIXEL_FLOOR = 1
# A cell where the two baselines disagree by more than this fraction of the lane's
# signal scale is fragile; there the cluster-inclusive opening net is reported.
BASELINE_FRAGILE_DISAGREEMENT_FRACTION = 0.25
BAND_MEASUREMENTS_CSV_FILENAME = "band_measurements.csv"
REFERENCE_QUALITY_FILENAME = "reference_quality.json"
BAND_MEASUREMENT_PLOT_FILENAME = "band_measurement_overlay.png"

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
    prog="analyze_gel.py",
    description="Validate, detrend, detect bands, and measure a manual gel export.",
    allow_abbrev=False,
)
# Directory-addressed, matching the sheet validator and the manifest: point at the
# gel analysis directory (or any file in it) and the standard files are resolved
# inside it. The sample sheet is REQUIRED, not optional: this stage measures only
# the lanes the sheet marks measured/reference, so running without it is undefined.
argument_parser.add_argument(
    "gel_directory_or_file", help="The gel analysis directory, or any file inside it"
)
argument_parser.add_argument(
    "--profile-csv", default=None, help="Override profile CSV path"
)
argument_parser.add_argument(
    "--sample-sheet", default=None, help="Override sample sheet path"
)
parsed_arguments = argument_parser.parse_args()

STANDARD_PROFILE_FILENAME = "manual_lane_profiles.csv"
STANDARD_SAMPLE_SHEET_FILENAME = "sample_sheet.csv"

gel_path_argument = pathlib.Path(
    os.path.abspath(parsed_arguments.gel_directory_or_file)
)
if gel_path_argument.is_dir():
    gel_directory_path = gel_path_argument
elif gel_path_argument.is_file():
    gel_directory_path = gel_path_argument.parent
else:
    die("input", "not a file or directory: " + str(gel_path_argument))
output_directory_path = gel_directory_path
emit_message("input", "gel analysis directory: " + str(gel_directory_path))

if parsed_arguments.profile_csv is not None:
    input_profile_csv_path = pathlib.Path(os.path.abspath(parsed_arguments.profile_csv))
else:
    input_profile_csv_path = gel_directory_path / STANDARD_PROFILE_FILENAME
if parsed_arguments.sample_sheet is not None:
    sample_sheet_csv_path = pathlib.Path(os.path.abspath(parsed_arguments.sample_sheet))
else:
    sample_sheet_csv_path = gel_directory_path / STANDARD_SAMPLE_SHEET_FILENAME

if not input_profile_csv_path.is_file():
    die(
        "input",
        "profile CSV not found at "
        + str(input_profile_csv_path)
        + ". Run the ImageJ export macro first, or pass --profile-csv.",
    )
if not sample_sheet_csv_path.is_file():
    die(
        "input",
        "sample sheet not found at "
        + str(sample_sheet_csv_path)
        + ". Author and validate it first, or pass --sample-sheet.",
    )
emit_message("input", "profile CSV: " + str(input_profile_csv_path))
emit_message("input", "sample sheet: " + str(sample_sheet_csv_path))

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
    "at_ceiling_count",
]

with open(
    input_profile_csv_path, newline="", encoding="utf-8-sig"
) as input_file_handle:
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
        + ". Re-export with the current export_lane_profiles.ijm (needs at_ceiling_count).",
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
# Read the sample sheet. It is required and assumed already validated by
# validate_sample_sheet.py; a few load-bearing invariants are re-checked here
# (lanes match, exactly one reference) because this CSV may outlive that run, and
# a mismatch would silently measure the wrong lanes.
# =============================================================================

with open(sample_sheet_csv_path, newline="", encoding="utf-8-sig") as sheet_file_handle:
    sheet_reader = csv.DictReader(sheet_file_handle)
    sheet_rows = []
    for raw_sheet_row in sheet_reader:
        sheet_rows.append(
            {key: (value or "").strip() for key, value in raw_sheet_row.items()}
        )

for required_sheet_column in (
    "lane_index",
    "well_number",
    "lane_content",
    "analysis_role",
):
    if not sheet_rows or required_sheet_column not in sheet_rows[0]:
        die(
            "sheet",
            "sample sheet missing column %r; validate it first."
            % required_sheet_column,
        )

well_number_by_lane_index = {}
analysis_role_by_lane_index = {}
lane_content_by_lane_index = {}
sample_label_by_lane_index = {}
for sheet_row in sheet_rows:
    sheet_lane_index = int(sheet_row["lane_index"])
    well_number_by_lane_index[sheet_lane_index] = int(sheet_row["well_number"])
    analysis_role_by_lane_index[sheet_lane_index] = sheet_row["analysis_role"]
    lane_content_by_lane_index[sheet_lane_index] = sheet_row["lane_content"]
    sample_label_by_lane_index[sheet_lane_index] = sheet_row.get("sample_label", "")

if sorted(well_number_by_lane_index.keys()) != sorted_lane_indices:
    die(
        "sheet",
        "sample sheet lanes %s do not match profile lanes %s; wrong sheet?"
        % (sorted(well_number_by_lane_index.keys()), sorted_lane_indices),
    )
reference_lane_indices = [
    lane for lane, role in analysis_role_by_lane_index.items() if role == "reference"
]
if len(reference_lane_indices) != 1:
    die(
        "sheet",
        "need exactly one analysis_role=reference; found %d. Re-validate the sheet."
        % len(reference_lane_indices),
    )
reference_lane_index = reference_lane_indices[0]

# Lanes to measure: reference and measured. Excluded (ladder/empty) are skipped for
# integration but stay in detection and plots.
measured_lane_indices = sorted(
    lane
    for lane, role in analysis_role_by_lane_index.items()
    if role in ("reference", "measured")
)
emit_message(
    "sheet",
    "%d measured lane(s), reference is lane %d (well %d, %s)"
    % (
        len(measured_lane_indices),
        reference_lane_index,
        well_number_by_lane_index[reference_lane_index],
        sample_label_by_lane_index.get(reference_lane_index, ""),
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

# =============================================================================
# Band detection. Two models computed side by side so the choice between them is
# made from the diagnostics below, not assumed: per-lane (each lane's own peaks)
# and consensus (peaks that line up across lanes, clustered as measure_gel.py
# clusters them). Integration is deliberately NOT done here; this stage proposes
# bands and reports how well the two models agree, then stops for a human gate.
# =============================================================================

# The single mm-per-pixel already asserted constant above; recover it from any
# row with a non-zero position so windows can be reported in mm.
millimetres_per_pixel = 0.0
for row in all_rows:
    position_pixels = int(row["migration_position_pixels"])
    if position_pixels > 0:
        millimetres_per_pixel = (
            float(row["migration_position_millimetres"]) / position_pixels
        )
        break
if millimetres_per_pixel <= 0.0:
    die(
        "bands",
        "cannot recover mm-per-pixel; migration_position_millimetres looks unset.",
    )

band_smoothing_pixels = max(
    1, int(round(BAND_MIGRATION_PROFILE_SMOOTHING_MILLIMETRES / millimetres_per_pixel))
)
band_minimum_separation_pixels = max(
    1, int(round(BAND_MINIMUM_SEPARATION_MILLIMETRES / millimetres_per_pixel))
)
band_region_window_pixels = BAND_REGION_WINDOW_WIDTH_MILLIMETRES / millimetres_per_pixel
apex_exclusion_pixels = int(round(APEX_EXCLUSION_MILLIMETRES / millimetres_per_pixel))

# Per-lane detection. find_peaks on the smoothed raw profile with lane-scaled
# prominence, dropping any peak inside the apex zone. Kept per lane so both the
# per-lane model and the consensus pool are built from the same peaks.
peaks_by_lane_index = {}
band_candidate_rows = []
for lane_index_value in sorted_lane_indices:
    lane_rows = rows_by_lane_index[lane_index_value]
    raw_lane_profile = numpy.array(
        [float(r["raw_value"]) for r in lane_rows], dtype=float
    )
    smoothed_lane_profile = scipy.ndimage.uniform_filter1d(
        raw_lane_profile, size=band_smoothing_pixels, mode="nearest"
    )
    lane_profile_maximum = float(smoothed_lane_profile.max())
    detected_peak_positions, _ = scipy.signal.find_peaks(
        smoothed_lane_profile,
        distance=band_minimum_separation_pixels,
        prominence=BAND_PEAK_PROMINENCE_FRACTION * lane_profile_maximum,
    )
    kept_peak_positions = [
        int(p) for p in detected_peak_positions if int(p) >= apex_exclusion_pixels
    ]
    peaks_by_lane_index[lane_index_value] = kept_peak_positions
    for peak_rank, peak_position in enumerate(detected_peak_positions):
        band_candidate_rows.append(
            {
                "lane_index": lane_index_value,
                "peak_position_pixels": int(peak_position),
                "peak_position_millimetres": "%.4f"
                % (int(peak_position) * millimetres_per_pixel),
                "peak_smoothed_value": "%.2f"
                % float(smoothed_lane_profile[int(peak_position)]),
                "in_apex_zone": "yes"
                if int(peak_position) < apex_exclusion_pixels
                else "no",
            }
        )

band_candidates_csv_path = output_directory_path / BAND_CANDIDATES_CSV_FILENAME
with open(band_candidates_csv_path, "w", newline="") as candidate_file_handle:
    candidate_writer = csv.DictWriter(
        candidate_file_handle,
        fieldnames=[
            "lane_index",
            "peak_position_pixels",
            "peak_position_millimetres",
            "peak_smoothed_value",
            "in_apex_zone",
        ],
    )
    candidate_writer.writeheader()
    for candidate_row in band_candidate_rows:
        candidate_writer.writerow(candidate_row)
emit_message("bands", "wrote " + str(band_candidates_csv_path))

# Consensus clustering, verbatim in spirit from measure_gel.py: pool every kept
# peak with its lane, sort, then single-linkage cluster within the separation
# tolerance with a same-lane exclusion so one lane cannot appear twice in a band.
pooled_peak_records = []
for lane_index_value in sorted_lane_indices:
    for peak_position in peaks_by_lane_index[lane_index_value]:
        pooled_peak_records.append((peak_position, lane_index_value))
pooled_peak_records.sort()

consensus_clusters = []
if pooled_peak_records:
    current_cluster = [pooled_peak_records[0]]
    for peak_position, lane_index_value in pooled_peak_records[1:]:
        lanes_in_current_cluster = {member_lane for _, member_lane in current_cluster}
        if (
            peak_position - current_cluster[-1][0] <= band_minimum_separation_pixels
            and lane_index_value not in lanes_in_current_cluster
        ):
            current_cluster.append((peak_position, lane_index_value))
        else:
            consensus_clusters.append(current_cluster)
            current_cluster = [(peak_position, lane_index_value)]
    consensus_clusters.append(current_cluster)

# One canonical position per cluster (median of members), with occupancy and the
# cross-lane spread that says whether the peaks actually line up.
consensus_band_records = []
for cluster in consensus_clusters:
    member_positions = [position for position, _ in cluster]
    member_lanes = sorted({member_lane for _, member_lane in cluster})
    cross_lane_spread_pixels = int(max(member_positions) - min(member_positions))
    consensus_band_records.append(
        {
            "canonical_position_pixels": float(numpy.median(member_positions)),
            "member_lanes": member_lanes,
            "occupancy": len(member_lanes),
            "cross_lane_spread_pixels": cross_lane_spread_pixels,
        }
    )
consensus_band_records.sort(key=lambda record: record["canonical_position_pixels"])

# Fixed-width windows clipped at neighbour midpoints so no two windows overlap.
# The per-band decision flags (well-supported, spread-fits-window) are computed
# here, once the window exists, so the CSV and the diagnostics report the same
# flags rather than one preceding the other.
lane_count_for_support = len(sorted_lane_indices)
for band_position_index in range(len(consensus_band_records)):
    canonical_position = consensus_band_records[band_position_index][
        "canonical_position_pixels"
    ]
    window_start = canonical_position - band_region_window_pixels / 2.0
    window_end = canonical_position + band_region_window_pixels / 2.0
    if band_position_index > 0:
        left_midpoint = 0.5 * (
            consensus_band_records[band_position_index - 1]["canonical_position_pixels"]
            + canonical_position
        )
        window_start = max(window_start, left_midpoint)
    if band_position_index < len(consensus_band_records) - 1:
        right_midpoint = 0.5 * (
            canonical_position
            + consensus_band_records[band_position_index + 1][
                "canonical_position_pixels"
            ]
        )
        window_end = min(window_end, right_midpoint)
    clipped_window_start = max(0, int(round(window_start)))
    clipped_window_end = int(round(window_end))
    consensus_band_records[band_position_index]["window_start_pixels"] = (
        clipped_window_start
    )
    consensus_band_records[band_position_index]["window_end_pixels"] = (
        clipped_window_end
    )
    # A band is safe for consensus integration when its cross-lane spread fits
    # inside the window every member is measured in; smile wider than the window
    # is what would actually corrupt the measurement.
    band_window_width_pixels = clipped_window_end - clipped_window_start
    consensus_band_records[band_position_index]["spread_fits_window"] = (
        consensus_band_records[band_position_index]["cross_lane_spread_pixels"]
        <= band_window_width_pixels
    )
    consensus_band_records[band_position_index]["is_well_supported"] = (
        consensus_band_records[band_position_index]["occupancy"]
        >= CONSENSUS_WELL_SUPPORTED_LANE_FRACTION * lane_count_for_support
    )

consensus_bands_csv_path = output_directory_path / CONSENSUS_BANDS_CSV_FILENAME
with open(consensus_bands_csv_path, "w", newline="") as consensus_file_handle:
    consensus_writer = csv.DictWriter(
        consensus_file_handle,
        fieldnames=[
            "consensus_band_index",
            "canonical_position_pixels",
            "canonical_position_millimetres",
            "window_start_pixels",
            "window_end_pixels",
            "occupancy",
            "member_lanes",
            "cross_lane_spread_pixels",
            "cross_lane_spread_millimetres",
            "is_well_supported",
            "spread_fits_window",
        ],
    )
    consensus_writer.writeheader()
    for consensus_band_index, consensus_band in enumerate(consensus_band_records):
        consensus_writer.writerow(
            {
                "consensus_band_index": consensus_band_index,
                "canonical_position_pixels": "%.1f"
                % consensus_band["canonical_position_pixels"],
                "canonical_position_millimetres": "%.4f"
                % (consensus_band["canonical_position_pixels"] * millimetres_per_pixel),
                "window_start_pixels": consensus_band["window_start_pixels"],
                "window_end_pixels": consensus_band["window_end_pixels"],
                "occupancy": consensus_band["occupancy"],
                "member_lanes": " ".join(
                    str(member_lane) for member_lane in consensus_band["member_lanes"]
                ),
                "cross_lane_spread_pixels": consensus_band["cross_lane_spread_pixels"],
                "cross_lane_spread_millimetres": "%.4f"
                % (consensus_band["cross_lane_spread_pixels"] * millimetres_per_pixel),
                "is_well_supported": "yes"
                if consensus_band["is_well_supported"]
                else "no",
                "spread_fits_window": "yes"
                if consensus_band["spread_fits_window"]
                else "no",
            }
        )
emit_message("bands", "wrote " + str(consensus_bands_csv_path))

# Agreement diagnostics: the numbers that decide consensus vs per-lane. A peak is
# "consensus-explained" if it sits in a cluster that both spans enough lanes and
# is tight; the fraction of peaks so explained is the headline score.
total_kept_peaks = sum(len(peaks_by_lane_index[lane]) for lane in sorted_lane_indices)
# Score consensus by what actually breaks integration: a well-supported band is
# safe if its cross-lane spread fits inside its own window (flags already set per
# band above). Absolute mm thresholds punish honest smile; window-fit does not.
well_supported_band_count = 0
peaks_in_well_supported_bands = 0
peaks_in_safe_bands = 0
singleton_band_count = 0
max_cross_lane_spread_millimetres = 0.0
well_supported_spread_values_millimetres = []
for consensus_band in consensus_band_records:
    spread_millimetres = (
        consensus_band["cross_lane_spread_pixels"] * millimetres_per_pixel
    )
    max_cross_lane_spread_millimetres = max(
        max_cross_lane_spread_millimetres, spread_millimetres
    )
    if consensus_band["is_well_supported"]:
        well_supported_band_count += 1
        peaks_in_well_supported_bands += consensus_band["occupancy"]
        well_supported_spread_values_millimetres.append(spread_millimetres)
        if consensus_band["spread_fits_window"]:
            peaks_in_safe_bands += consensus_band["occupancy"]
    if consensus_band["occupancy"] == 1:
        singleton_band_count += 1
# Headline: of the peaks that landed in well-supported bands, what fraction sit in
# bands whose spread still fits the window. Near 1.0 means consensus is safe here.
consensus_safe_fraction = (
    peaks_in_safe_bands / float(peaks_in_well_supported_bands)
    if peaks_in_well_supported_bands
    else 0.0
)
median_well_supported_spread_millimetres = (
    float(numpy.median(well_supported_spread_values_millimetres))
    if well_supported_spread_values_millimetres
    else 0.0
)

band_model_diagnostics = {
    "schema_version": "band_model_diagnostics_1",
    "generated_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "millimetres_per_pixel": millimetres_per_pixel,
    "apex_exclusion_millimetres": APEX_EXCLUSION_MILLIMETRES,
    "lane_count": lane_count_for_support,
    "per_lane_peak_counts": {
        str(lane): len(peaks_by_lane_index[lane]) for lane in sorted_lane_indices
    },
    "total_kept_peaks": total_kept_peaks,
    "consensus_band_count": len(consensus_band_records),
    "well_supported_band_count": well_supported_band_count,
    "singleton_band_count": singleton_band_count,
    "consensus_safe_fraction": round(consensus_safe_fraction, 4),
    "median_well_supported_cross_lane_spread_millimetres": round(
        median_well_supported_spread_millimetres, 4
    ),
    "max_cross_lane_spread_millimetres": round(max_cross_lane_spread_millimetres, 4),
    "region_window_width_millimetres": BAND_REGION_WINDOW_WIDTH_MILLIMETRES,
    "interpretation": (
        "consensus_safe_fraction near 1.0 means the well-supported bands' cross-lane "
        "spread fits inside their integration windows, so consensus is safe. A large "
        "median spread relative to the window is lane smile: widen the window or "
        "prefer a reference template. Many singletons favour per-lane."
    ),
}
band_model_diagnostics_path = output_directory_path / BAND_MODEL_DIAGNOSTICS_FILENAME
band_model_diagnostics_path.write_text(
    json.dumps(band_model_diagnostics, indent=2) + "\n"
)
emit_message(
    "bands",
    "wrote "
    + str(band_model_diagnostics_path)
    + "; consensus safe fraction = %.2f, median smile = %.2f mm over %d consensus bands"
    % (
        consensus_safe_fraction,
        median_well_supported_spread_millimetres,
        len(consensus_band_records),
    ),
)

# Gate plot: every lane profile with its own detected peaks (dots), the consensus
# windows shaded across all lanes, and the apex exclusion zone hatched.
band_figure, band_axes = matplotlib.pyplot.subplots(
    grid_row_count,
    grid_column_count,
    figsize=(grid_column_count * 2.6, grid_row_count * 2.0),
    squeeze=False,
)
apex_zone_millimetres = apex_exclusion_pixels * millimetres_per_pixel
for panel_index in range(grid_row_count * grid_column_count):
    panel_axis = band_axes[panel_index // grid_column_count][
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
    panel_axis.plot(
        position_millimetres, raw_lane_profile, color="black", linewidth=0.7
    )
    for consensus_band in consensus_band_records:
        panel_axis.axvspan(
            consensus_band["window_start_pixels"] * millimetres_per_pixel,
            consensus_band["window_end_pixels"] * millimetres_per_pixel,
            color="tab:blue",
            alpha=0.10,
            linewidth=0,
        )
    panel_axis.axvspan(
        0.0,
        apex_zone_millimetres,
        facecolor="none",
        edgecolor="0.5",
        hatch="///",
        linewidth=0,
    )
    for peak_position in peaks_by_lane_index[lane_index_value]:
        panel_axis.plot(
            peak_position * millimetres_per_pixel,
            raw_lane_profile[peak_position],
            marker="v",
            color="tab:red",
            markersize=4,
        )
    panel_axis.set_title(
        "lane %d (%d peaks)"
        % (lane_index_value, len(peaks_by_lane_index[lane_index_value])),
        fontsize="small",
    )
    panel_axis.tick_params(labelsize=6)
band_figure.supxlabel("migration position (mm)", fontsize="small")
band_figure.suptitle(
    "Band detection: red = per-lane peaks, blue = consensus windows, hatched = apex excluded",
    fontsize="small",
)
band_figure.tight_layout()
band_detection_plot_path = output_directory_path / BAND_DETECTION_PLOT_FILENAME
band_figure.savefig(band_detection_plot_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(band_figure)
emit_message("bands", "wrote " + str(band_detection_plot_path))

# =============================================================================
# Measurement. Integrate each MEASURED lane at the consensus windows, reusing
# measure_gel.py's math: a valley-to-valley baseline (flanking minima, sloped line,
# region net, negative-clip) and a rolling-ball opening cross-check, reporting both
# a region-net AREA and an apex HEIGHT per band. Saturation (from at_ceiling_count)
# and baseline fragility select the reported basis exactly as measure_gel.py does.
# Excluded lanes (ladder/empty) are not measured.
# =============================================================================

baseline_flank_search_pixels = max(
    1, int(round(BASELINE_FLANK_SEARCH_MILLIMETRES / millimetres_per_pixel))
)
rolling_ball_width_pixels = max(
    1, int(round(ROLLING_BALL_WIDTH_MILLIMETRES / millimetres_per_pixel))
)

# Per-lane raw profile and per-row at-ceiling count, keyed by lane, for the lanes
# we measure. Also a lane signal scale (the profile's own spread) used to judge
# baseline fragility, matching measure_gel.py's use of a lane-scaled tolerance.
raw_profile_by_lane = {}
at_ceiling_by_lane = {}
lane_signal_scale_by_lane = {}
for lane_index_value in measured_lane_indices:
    lane_rows = rows_by_lane_index[lane_index_value]
    raw_profile_by_lane[lane_index_value] = numpy.array(
        [float(r["raw_value"]) for r in lane_rows], dtype=float
    )
    at_ceiling_by_lane[lane_index_value] = numpy.array(
        [int(r["at_ceiling_count"]) for r in lane_rows], dtype=int
    )
    lane_profile = raw_profile_by_lane[lane_index_value]
    lane_signal_scale_by_lane[lane_index_value] = float(
        numpy.percentile(lane_profile, 99.0) - numpy.median(lane_profile)
    )

migration_extent_pixels = raw_profile_by_lane[measured_lane_indices[0]].size

band_measurement_rows = []
for lane_index_value in measured_lane_indices:
    lane_profile = raw_profile_by_lane[lane_index_value]
    lane_ceiling_counts = at_ceiling_by_lane[lane_index_value]
    lane_signal_scale = lane_signal_scale_by_lane[lane_index_value]
    # Rolling-ball opening baseline for the whole lane once (the cross-check).
    opening_odd_width = rolling_ball_width_pixels + (1 - rolling_ball_width_pixels % 2)
    opening_baseline = scipy.ndimage.grey_opening(
        lane_profile, size=opening_odd_width, mode="nearest"
    )

    for consensus_band_index, consensus_band in enumerate(consensus_band_records):
        window_start = consensus_band["window_start_pixels"]
        window_end = consensus_band["window_end_pixels"]

        # --- valley-to-valley baseline: flanking minima just outside the window ---
        left_search_start = max(0, window_start - baseline_flank_search_pixels)
        left_flank = lane_profile[left_search_start : window_start + 1]
        if left_flank.size:
            left_reference_column = left_search_start + int(numpy.argmin(left_flank))
            left_reference_value = float(lane_profile[left_reference_column])
        else:
            left_reference_column = window_start
            left_reference_value = float(lane_profile[window_start])
        right_search_end = min(
            migration_extent_pixels, window_end + baseline_flank_search_pixels
        )
        right_flank = lane_profile[max(window_start, window_end - 1) : right_search_end]
        if right_flank.size:
            right_reference_column = max(window_start, window_end - 1) + int(
                numpy.argmin(right_flank)
            )
            right_reference_value = float(lane_profile[right_reference_column])
        else:
            right_reference_column = max(window_start, window_end - 1)
            right_reference_value = float(lane_profile[right_reference_column])

        window_columns = numpy.arange(window_start, window_end)
        if window_columns.size == 0:
            valley_net_area = 0.0
            opening_net_area = 0.0
            apex_height_above_valley = 0.0
            apex_raw_value = 0.0
            apex_migration_pixels = int(window_start)
            measurement_status = "empty_window"
        else:
            if right_reference_column != left_reference_column:
                valley_baseline = left_reference_value + (
                    (right_reference_value - left_reference_value)
                    * (window_columns - left_reference_column)
                    / float(right_reference_column - left_reference_column)
                )
            else:
                valley_baseline = numpy.full(window_columns.shape, left_reference_value)
            profile_over_window = lane_profile[window_start:window_end]
            valley_net_unclipped = float((profile_over_window - valley_baseline).sum())
            opening_net_area = float(
                (profile_over_window - opening_baseline[window_start:window_end]).sum()
            )
            if valley_net_unclipped < 0.0:
                valley_net_area = 0.0
                if window_start == 0 or window_end == migration_extent_pixels:
                    measurement_status = "window_truncated_at_crop_edge"
                else:
                    measurement_status = "baseline_above_signal_clipped_to_zero"
            else:
                valley_net_area = valley_net_unclipped
                if window_start == 0 or window_end == migration_extent_pixels:
                    measurement_status = "window_truncated_at_crop_edge"
                else:
                    measurement_status = "ok"
            apex_offset = int(numpy.argmax(profile_over_window))
            apex_raw_value = float(profile_over_window[apex_offset])
            apex_migration_pixels = int(window_start + apex_offset)
            apex_height_above_valley = float(
                (profile_over_window - valley_baseline).max()
            )

        # --- saturation, from the at-ceiling counts inside the window ---
        window_at_ceiling_total = int(
            lane_ceiling_counts[window_start:window_end].sum()
        )
        is_saturated = window_at_ceiling_total >= SATURATION_AT_CEILING_PIXEL_FLOOR

        # --- baseline fragility: the two nets disagree relative to lane scale ---
        if lane_signal_scale > 0.0:
            baseline_disagreement_fraction = (
                abs(valley_net_area - opening_net_area) / lane_signal_scale
            )
        else:
            baseline_disagreement_fraction = 0.0
        baseline_is_fragile = (
            baseline_disagreement_fraction >= BASELINE_FRAGILE_DISAGREEMENT_FRACTION
        )

        # --- reported-value contract, verbatim from measure_gel.py: saturation is a
        # floor and wins first; a fragile cell reports the cluster-inclusive opening
        # net; otherwise the valley-to-valley region net. Basis recorded either way.
        if is_saturated:
            reported_area = valley_net_area
            reported_value_basis = "saturated_floor"
        elif baseline_is_fragile:
            reported_area = opening_net_area
            reported_value_basis = "opening_net_cluster_inclusive"
        else:
            reported_area = valley_net_area
            reported_value_basis = "region_net"

        band_measurement_rows.append(
            {
                "gel_id": input_profile_csv_path.parent.name,
                "well_number": well_number_by_lane_index[lane_index_value],
                "lane_index": lane_index_value,
                "sample_label": sample_label_by_lane_index.get(lane_index_value, ""),
                "analysis_role": analysis_role_by_lane_index[lane_index_value],
                "canonical_band_index": consensus_band_index,
                "canonical_position_millimetres": round(
                    consensus_band["canonical_position_pixels"] * millimetres_per_pixel,
                    4,
                ),
                "window_start_pixels": window_start,
                "window_end_pixels": window_end,
                "region_net_area": round(valley_net_area, 2),
                "opening_net_area": round(opening_net_area, 2),
                "reported_area": round(reported_area, 2),
                "reported_value_basis": reported_value_basis,
                "apex_height_above_baseline": round(apex_height_above_valley, 2),
                "apex_raw_value": round(apex_raw_value, 2),
                "apex_migration_pixels": apex_migration_pixels,
                "cross_lane_spread_millimetres": round(
                    consensus_band["cross_lane_spread_pixels"] * millimetres_per_pixel,
                    4,
                ),
                "window_at_ceiling_pixels": window_at_ceiling_total,
                "saturation_status": "saturated" if is_saturated else "clean",
                "baseline_disagreement_fraction": round(
                    baseline_disagreement_fraction, 5
                ),
                "baseline_agreement_status": "fragile"
                if baseline_is_fragile
                else "agree",
                "measurement_status": measurement_status,
                "band_occupancy": consensus_band["occupancy"],
                "band_is_well_supported": "yes"
                if consensus_band["is_well_supported"]
                else "no",
            }
        )

band_measurements_csv_path = output_directory_path / BAND_MEASUREMENTS_CSV_FILENAME
band_measurement_field_order = [
    "gel_id",
    "well_number",
    "lane_index",
    "sample_label",
    "analysis_role",
    "canonical_band_index",
    "canonical_position_millimetres",
    "window_start_pixels",
    "window_end_pixels",
    "region_net_area",
    "opening_net_area",
    "reported_area",
    "reported_value_basis",
    "apex_height_above_baseline",
    "apex_raw_value",
    "apex_migration_pixels",
    "cross_lane_spread_millimetres",
    "window_at_ceiling_pixels",
    "saturation_status",
    "baseline_disagreement_fraction",
    "baseline_agreement_status",
    "measurement_status",
    "band_occupancy",
    "band_is_well_supported",
]
with open(band_measurements_csv_path, "w", newline="") as measurement_file_handle:
    measurement_writer = csv.DictWriter(
        measurement_file_handle, fieldnames=band_measurement_field_order
    )
    measurement_writer.writeheader()
    for measurement_row in band_measurement_rows:
        measurement_writer.writerow(measurement_row)
emit_message(
    "measure",
    "wrote "
    + str(band_measurements_csv_path)
    + " (%d bands x %d measured lanes)"
    % (len(consensus_band_records), len(measured_lane_indices)),
)

# =============================================================================
# Reference quality. The reference lane's control band divides every normalized
# value, so score the designated reference and every input lane on the four
# properties that make a divisor trustworthy, and warn if the pick is weak. The
# control band is the earliest consensus band (canonical_band_index 0) by default,
# matching the loading-control convention.
# =============================================================================

CONTROL_BAND_CANONICAL_INDEX = 0
measurements_by_lane_and_band = {}
for measurement_row in band_measurement_rows:
    measurements_by_lane_and_band[
        (measurement_row["lane_index"], measurement_row["canonical_band_index"])
    ] = measurement_row

# Total lane signal (for the "typical loading" property) over the measured lanes.
total_signal_by_lane = {
    lane: float(raw_profile_by_lane[lane].sum()) for lane in measured_lane_indices
}
median_total_signal = float(numpy.median(list(total_signal_by_lane.values())))

reference_quality_by_lane = {}
for lane_index_value in measured_lane_indices:
    control_measurement = measurements_by_lane_and_band.get(
        (lane_index_value, CONTROL_BAND_CANONICAL_INDEX)
    )
    if control_measurement is None:
        continue
    control_not_saturated = control_measurement["saturation_status"] == "clean"
    lane_scale = lane_signal_scale_by_lane[lane_index_value]
    control_snr = (
        (control_measurement["apex_height_above_baseline"] / lane_scale)
        if lane_scale > 0
        else 0.0
    )
    # Alignment: this lane's deviation of its control apex from the consensus centre.
    control_band_record = consensus_band_records[CONTROL_BAND_CANONICAL_INDEX]
    consensus_centre_millimetres = (
        control_band_record["canonical_position_pixels"] * millimetres_per_pixel
    )
    control_alignment_deviation_millimetres = abs(
        control_measurement["apex_migration_pixels"] * millimetres_per_pixel
        - consensus_centre_millimetres
    )
    loading_ratio = (
        (total_signal_by_lane[lane_index_value] / median_total_signal)
        if median_total_signal > 0
        else 0.0
    )
    reference_quality_by_lane[lane_index_value] = {
        "well_number": well_number_by_lane_index[lane_index_value],
        "sample_label": sample_label_by_lane_index.get(lane_index_value, ""),
        "control_not_saturated": control_not_saturated,
        "control_band_snr": round(control_snr, 3),
        "control_alignment_deviation_mm": round(
            control_alignment_deviation_millimetres, 3
        ),
        "loading_ratio_vs_median": round(loading_ratio, 3),
    }


# Rank candidate references: disqualify saturated control bands, then prefer high
# SNR, tight alignment, and loading near the median. A single scalar for ordering.
def reference_candidate_score(lane_index_value):
    quality = reference_quality_by_lane[lane_index_value]
    if not quality["control_not_saturated"]:
        return -1.0
    alignment_penalty = 1.0 / (1.0 + quality["control_alignment_deviation_mm"])
    loading_penalty = 1.0 / (1.0 + abs(quality["loading_ratio_vs_median"] - 1.0))
    return quality["control_band_snr"] * alignment_penalty * loading_penalty


ranked_reference_candidates = sorted(
    reference_quality_by_lane.keys(), key=reference_candidate_score, reverse=True
)
best_candidate_lane = (
    ranked_reference_candidates[0] if ranked_reference_candidates else None
)
designated_score = reference_candidate_score(reference_lane_index)
best_score = (
    reference_candidate_score(best_candidate_lane)
    if best_candidate_lane is not None
    else 0.0
)

reference_warnings = []
if not reference_quality_by_lane.get(reference_lane_index, {}).get(
    "control_not_saturated", True
):
    reference_warnings.append(
        "designated reference lane %d has a SATURATED control band"
        % reference_lane_index
    )
if (
    best_candidate_lane is not None
    and best_candidate_lane != reference_lane_index
    and best_score > 1.5 * max(designated_score, 1e-9)
):
    reference_warnings.append(
        "lane %d (well %d) scores clearly higher as reference than your pick, lane %d (well %d)"
        % (
            best_candidate_lane,
            well_number_by_lane_index[best_candidate_lane],
            reference_lane_index,
            well_number_by_lane_index[reference_lane_index],
        )
    )

reference_quality_report = {
    "schema_version": "reference_quality_1",
    "generated_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "designated_reference_lane_index": reference_lane_index,
    "designated_reference_well_number": well_number_by_lane_index[reference_lane_index],
    "control_band_canonical_index": CONTROL_BAND_CANONICAL_INDEX,
    "best_candidate_lane_index": best_candidate_lane,
    "per_lane_quality": {
        str(lane): reference_quality_by_lane[lane]
        for lane in sorted(reference_quality_by_lane)
    },
    "warnings": reference_warnings,
}
reference_quality_path = output_directory_path / REFERENCE_QUALITY_FILENAME
reference_quality_path.write_text(json.dumps(reference_quality_report, indent=2) + "\n")
emit_message("reference", "wrote " + str(reference_quality_path))
for reference_warning in reference_warnings:
    emit_message("reference", "WARNING " + reference_warning)

# =============================================================================
# Measurement overlay: per measured lane, the profile with consensus windows, the
# apex marked, and saturated bands flagged, so the numbers have a visual gate.
# =============================================================================

measured_panel_count = len(measured_lane_indices)
measure_grid_columns = int(math.ceil(math.sqrt(measured_panel_count)))
measure_grid_rows = int(math.ceil(measured_panel_count / float(measure_grid_columns)))
measure_figure, measure_axes = matplotlib.pyplot.subplots(
    measure_grid_rows,
    measure_grid_columns,
    figsize=(measure_grid_columns * 2.6, measure_grid_rows * 2.0),
    squeeze=False,
)
for panel_index in range(measure_grid_rows * measure_grid_columns):
    panel_axis = measure_axes[panel_index // measure_grid_columns][
        panel_index % measure_grid_columns
    ]
    if panel_index >= measured_panel_count:
        panel_axis.axis("off")
        continue
    lane_index_value = measured_lane_indices[panel_index]
    lane_profile = raw_profile_by_lane[lane_index_value]
    position_millimetres = numpy.arange(lane_profile.size) * millimetres_per_pixel
    panel_axis.plot(position_millimetres, lane_profile, color="black", linewidth=0.7)
    for consensus_band_index, consensus_band in enumerate(consensus_band_records):
        measurement_row = measurements_by_lane_and_band.get(
            (lane_index_value, consensus_band_index)
        )
        window_color = (
            "tab:red"
            if (measurement_row and measurement_row["saturation_status"] == "saturated")
            else "tab:blue"
        )
        panel_axis.axvspan(
            consensus_band["window_start_pixels"] * millimetres_per_pixel,
            consensus_band["window_end_pixels"] * millimetres_per_pixel,
            color=window_color,
            alpha=0.10,
            linewidth=0,
        )
        if measurement_row:
            panel_axis.plot(
                measurement_row["apex_migration_pixels"] * millimetres_per_pixel,
                measurement_row["apex_raw_value"],
                marker="v",
                color="tab:green",
                markersize=3,
            )
    role_tag = "REF" if lane_index_value == reference_lane_index else ""
    panel_axis.set_title(
        "lane %d w%d %s"
        % (lane_index_value, well_number_by_lane_index[lane_index_value], role_tag),
        fontsize="small",
    )
    panel_axis.tick_params(labelsize=6)
measure_figure.supxlabel("migration position (mm)", fontsize="small")
measure_figure.suptitle(
    "Measurement: blue=window, red=saturated window, green=apex. REF=reference lane.",
    fontsize="small",
)
measure_figure.tight_layout()
band_measurement_plot_path = output_directory_path / BAND_MEASUREMENT_PLOT_FILENAME
measure_figure.savefig(band_measurement_plot_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(measure_figure)
emit_message("measure", "wrote " + str(band_measurement_plot_path))

saturated_measurement_count = sum(
    1 for row in band_measurement_rows if row["saturation_status"] == "saturated"
)
fragile_measurement_count = sum(
    1 for row in band_measurement_rows if row["baseline_agreement_status"] == "fragile"
)
emit_message(
    "done",
    "measured %d bands x %d lanes; %d saturated, %d fragile; %d reference warning(s). "
    "Inspect %s, %s, and %s."
    % (
        len(consensus_band_records),
        len(measured_lane_indices),
        saturated_measurement_count,
        fragile_measurement_count,
        len(reference_warnings),
        BAND_MEASUREMENT_PLOT_FILENAME,
        BAND_MEASUREMENTS_CSV_FILENAME,
        REFERENCE_QUALITY_FILENAME,
    ),
)
sys.exit(0)
