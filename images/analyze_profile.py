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

emit_message(
    "done",
    "overall %s; %d warning(s). Inspect %s, %s, and %s before choosing a band model."
    % (
        checks_report["overall_status"],
        len(soft_warnings),
        BAND_DETECTION_PLOT_FILENAME,
        BAND_MODEL_DIAGNOSTICS_FILENAME,
        CONSENSUS_BANDS_CSV_FILENAME,
    ),
)
sys.exit(0)
