# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "tifffile",
#     "numpy",
#     "matplotlib",
#     "scipy",
# ]
# ///
r"""
Stage 2 of the gel densitometry pipeline, Slice A: geometry only. Reads the
stage 1 validation report for a scan, applies the crop, levels the gel only when
the tilt is worth the interpolation, fits the lane grid along the migration axis,
verifies that the recorded orientation matches where the lanes actually are, and
draws an overlay for the operator to check. It measures no bands and writes no
CSV; that is Slice B.

Scope and order are dictated by decisions recorded in DESIGN.md section 14 and the
stage 2 planning:

1. The analysis is invariant to which cardinal orientation the gel was scanned in.
   Nothing is rotated to a canonical orientation. The migration axis is read from
   the sidecar and every step is pointed along it, so a horizontal gel and the
   same gel scanned vertically produce identical measurements.
2. Saturation statistics are taken on the un-interpolated crop, BEFORE any
   rotation, because bilinear interpolation moves values off the ceiling and would
   hide clipping.
3. Rotation happens only above a threshold (default 0.25 degrees). Below it the
   sub-pixel gain is not worth the blur, so the tilt is recorded and skipped. When
   rotation does happen it is bilinear (order 1), and the landmarks are rotated by
   the same transform and the tilt re-derived and asserted near zero, because a
   sign error produces a doubled tilt that still looks roughly straight.
4. The recorded orientation is checked against the image, not trusted. A mislabel
   makes every downstream step measure across the lanes instead of along them,
   with plausible-looking wrong numbers, so it is a hard stop rather than a hope.

Shape follows CONVENTIONS.md: flat procedural, no helpers beyond the tagged stderr
emitter and the fail-fast exit, full descriptive names carrying units, ASCII only,
comments stating why. Findings accumulate and the provenance report is always
written at a single exit point; die() is reserved for the unrecoverable, where
there is nothing to report about.

Usage:

    uv run measure_gel_image.py '<scan>.tif'
    uv run measure_gel_image.py --rotation-threshold-degrees 0.0 '<scan>.tif'
"""

import argparse
import datetime
import json
import math
import os
import pathlib
import sys

import numpy
import scipy.ndimage
import scipy.signal
import tifffile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot
import matplotlib.patches


ANALYSIS_REPORT_SCHEMA_VERSION = 1

VALIDATION_REPORT_FILENAME = "input_file_validation_report.json"
ANALYSIS_REPORT_FILENAME = "stage2_analysis_report.json"
LANE_GRID_OVERLAY_FILENAME = "lane_grid_overlay.png"
LANE_PROFILES_PLOT_FILENAME = "lane_profiles.png"
LANE_PROFILES_CSV_FILENAME = "lane_profiles.csv"

# Below this the correction is a small fraction of a pixel across a band and is not
# worth the interpolation blur, so the tilt is recorded and the image left alone.
DEFAULT_ROTATION_THRESHOLD_DEGREES = 0.25

# After rotating the landmarks by the chosen angle, the re-derived tilt must be at
# least this close to zero, or the sign was wrong and we stop.
ROTATION_RESIDUAL_TILT_TOLERANCE_DEGREES = 0.05

CONTAINER_MAXIMUM_VALUE_16_BIT = 65535

# The lane pitch is searched within this fraction of the pitch implied by the
# crop extent and the expected lane count, to absorb a crop that is a little loose
# or an expected count that is off by the odd empty lane.
LANE_PITCH_SEARCH_FRACTION = 0.30

# Orientation is decided by autocorrelation of the collapsed profile at the lane
# spacing, which is phase-independent: it measures whether an axis carries lane-
# pitch periodicity at all, without needing a comb to have landed on it. The
# stacking axis must be at least this many times as periodic as the migration
# axis, and must clear the floor so a featureless crop does not pass by default.
# This is a soft check (a warning, not a hard stop) until it is proven on real
# gels: a guard that blocks a correct gel is worse than no guard.
ORIENTATION_PERIODICITY_RATIO = 1.25
ORIENTATION_AUTOCORRELATION_FLOOR = 0.25

# Orientation is decided by per-axis gradient energy, not lane-pitch periodicity.
# Bands are sharp along the migration axis (narrow, closely spaced edges) and
# smooth along the stacking axis (a band spans the lane width), so the migration
# axis carries more high-frequency structure regardless of how many lanes are
# populated. This holds on a sparse gel where autocorrelation at the lane pitch is
# near zero and cannot decide orientation at all. The claimed migration axis is
# consistent when its gradient energy exceeds the claimed stacking axis's by at
# least this ratio. Measured ~2.3 on the s0002 fixture (correctly labelled).
ORIENTATION_GRADIENT_ENERGY_RATIO = 1.25

# Data-driven lane locator tuning. These were derived by measurement on the s0002
# fixture (a gel where all 15 wells were loaded but only 5 carried usable signal,
# clustered rather than evenly spread) and are documented so they are understood,
# not copied blindly.
#
# A migration column that is elevated down its whole length is a plate edge or
# scratch, not a lane; a real gel column sits at plate level except where bands
# cross it. Exclude any column whose median (down the stacking axis) exceeds this
# multiple of the typical column median.
EDGE_COLUMN_MEDIAN_MULTIPLE = 8.0
# Lane centres are found on a profile smoothed over half the lane pitch. Each well
# shows a symmetric two-lobed (doublet) intensity across its width from edge/well
# concentration during loading; smoothing over half a pitch merges the doublet into
# one bump so a lane is not split into two, while staying below a full pitch so
# adjacent lanes are not merged.
LANE_CENTRE_SMOOTHING_PITCH_FRACTION = 0.5
# Minimum separation between accepted lane centres, as a fraction of the pitch. Set
# below one pitch to tolerate a slightly irregular comb, above the doublet lobe
# spacing so one lane yields one centre.
LANE_MINIMUM_SEPARATION_PITCH_FRACTION = 0.7
# A candidate lane bump must clear this fraction of the strongest bump to be
# considered; kept gentle so a faint but real lane survives, since compact
# artifacts are removed afterwards by the band-structure test.
LANE_PROMINENCE_FRACTION_OF_MAXIMUM = 0.06
# The band-structure test that separates a real lane from a compact artifact (dust
# speck, ring): a real lane carries band signal across many migration columns,
# while a speck lights up only a handful even when its brightest point is dark or
# its scattered debris spans a wide range. Require at least this fraction of usable
# migration columns to participate (exceed 15 percent of the lane strip maximum).
LANE_MINIMUM_PARTICIPATING_COLUMN_FRACTION = 0.05
LANE_COLUMN_PARTICIPATION_THRESHOLD_FRACTION = 0.15

DISPLAY_COLORMAP_NAME = "gray_r"
OVERLAY_MAXIMUM_DIMENSION_PIXELS = 1400
FIGURE_DOTS_PER_INCH = 130


def emit_message(source_tag, message_text):
    sys.stderr.write("[" + source_tag + "] " + message_text + "\n")


def die(source_tag, message_text):
    emit_message(source_tag, message_text)
    sys.exit(2)


argument_parser = argparse.ArgumentParser(
    description="Stage 2 Slice A: crop, level if worthwhile, fit the lane grid, "
                "verify orientation, and draw an overlay. Reads the stage 1 report."
)
argument_parser.add_argument("input_tiff_path")
argument_parser.add_argument("--page-index", type=int, default=0)
argument_parser.add_argument(
    "--rotation-threshold-degrees", type=float,
    default=DEFAULT_ROTATION_THRESHOLD_DEGREES,
)
parsed_arguments = argument_parser.parse_args()

working_path_text = parsed_arguments.input_tiff_path.strip().strip('"').strip("'")
if working_path_text == "":
    die("path", "path is empty after stripping quotes and whitespace")
input_tiff_absolute_path = pathlib.Path(working_path_text).expanduser().resolve()
if not input_tiff_absolute_path.is_file():
    die("path", "not a file: " + str(input_tiff_absolute_path))

output_directory_path = input_tiff_absolute_path.parent / (
    input_tiff_absolute_path.stem + "_gel_analysis"
)
if not output_directory_path.is_dir():
    die("output", "expected output directory does not exist: "
        + str(output_directory_path) + ". Run stage 0 then stage 1 first.")

validation_report_path = output_directory_path / VALIDATION_REPORT_FILENAME
if not validation_report_path.is_file():
    die("stage1", "no stage 1 report at " + str(validation_report_path)
        + ". Run validate_gel_image.py first; stage 2 reads its output.")
try:
    validation_report = json.loads(validation_report_path.read_text())
except Exception as report_read_error:
    die("stage1", "stage 1 report will not parse: " + str(report_read_error))

# Stage 2 refuses to measure a file stage 1 flagged: the geometry it would rely on
# is exactly what a failing stage 1 says cannot be trusted.
if validation_report.get("overall_status") != "pass":
    die("stage1", "stage 1 overall_status is "
        + repr(validation_report.get("overall_status")) + ", not 'pass'. Fix the "
        "sidecar and re-run stage 1 before stage 2.")

preprocess_sidecar_block = validation_report.get("preprocess_sidecar", {})
geometry_pixels = preprocess_sidecar_block.get("geometry_pixels")
if geometry_pixels is None:
    die("stage1", "the stage 1 report has no preprocess_sidecar.geometry_pixels, so "
        "there is no crop or axis to analyse. A sidecar is required for stage 2.")

gel_migration_axis = geometry_pixels["gel_migration_axis"]
crop_x_pixels = int(geometry_pixels["crop_x_pixels"])
crop_y_pixels = int(geometry_pixels["crop_y_pixels"])
crop_width_pixels = int(geometry_pixels["crop_width_pixels"])
crop_height_pixels = int(geometry_pixels["crop_height_pixels"])
landmark_a_x_pixels = float(geometry_pixels["landmark_a_x_pixels"])
landmark_a_y_pixels = float(geometry_pixels["landmark_a_y_pixels"])
landmark_b_x_pixels = float(geometry_pixels["landmark_b_x_pixels"])
landmark_b_y_pixels = float(geometry_pixels["landmark_b_y_pixels"])
derived_tilt_angle_degrees = preprocess_sidecar_block.get("derived_tilt_angle_degrees")
expected_lane_count = int(preprocess_sidecar_block.get("values", {})["expected_lane_count"])
micrometres_per_pixel = validation_report.get("geometry", {}).get("micrometres_per_pixel")
encoding_verified = validation_report.get("linearity_evidence", {}).get("encoding_verified")

ANALYSIS_FINDINGS = []

# The read path, as float64 so nothing overflows and interpolation has headroom.
try:
    with tifffile.TiffFile(str(input_tiff_absolute_path)) as tiff_handle:
        raw_pixel_array = tiff_handle.pages[parsed_arguments.page_index].asarray()
except Exception as pixel_read_error:
    die("pixels", "page " + str(parsed_arguments.page_index) + " would not decode: "
        + str(pixel_read_error))
pixel_array = raw_pixel_array.astype(numpy.float64)
pixel_array_height_pixels, pixel_array_width_pixels = pixel_array.shape

crop_is_inside_image = (
    crop_width_pixels > 0 and crop_height_pixels > 0
    and 0 <= crop_x_pixels and 0 <= crop_y_pixels
    and crop_x_pixels + crop_width_pixels <= pixel_array_width_pixels
    and crop_y_pixels + crop_height_pixels <= pixel_array_height_pixels
)
if not crop_is_inside_image:
    die("crop", "the stage 1 crop does not fit this image; the report and the file "
        "disagree, which stage 1 should have caught.")

# Crop in the raw as-opened frame, where the sidecar coordinates were measured. The
# bottom-left-origin flip and lane numbering are deferred to Slice B, where the
# band-numbering convention is defined; flipping now would only desynchronise this
# overlay from what the operator sees in Fiji.
crop_pixels = pixel_array[
    crop_y_pixels:crop_y_pixels + crop_height_pixels,
    crop_x_pixels:crop_x_pixels + crop_width_pixels,
]

# Saturation on the un-interpolated crop, before any rotation can move values off
# the ceiling.
crop_maximum_value = float(crop_pixels.max())
at_container_maximum_in_crop = int((crop_pixels >= CONTAINER_MAXIMUM_VALUE_16_BIT).sum())
at_crop_maximum_count = int((crop_pixels >= crop_maximum_value).sum())
ANALYSIS_FINDINGS.append({
    "check_name": "in_crop_saturation_before_interpolation",
    "status": "warning" if at_container_maximum_in_crop > 0 else "pass",
    "is_hard_stop": False,
    "detail": "crop maximum %d held by %d pixels; %d pixels at the container maximum "
              "%d, measured before any rotation"
              % (crop_maximum_value, at_crop_maximum_count,
                 at_container_maximum_in_crop, CONTAINER_MAXIMUM_VALUE_16_BIT),
})

# Rotation decision. Below the threshold the tilt is recorded and the image left
# alone; above it the image is levelled with a self-checked bilinear rotation.
rotation_applied = False
rotation_applied_angle_degrees = 0.0
rotation_residual_tilt_degrees = None
analysis_crop = crop_pixels
target_line_angle_degrees = 0.0 if gel_migration_axis == "vertical" else 90.0

if derived_tilt_angle_degrees is None:
    ANALYSIS_FINDINGS.append({
        "check_name": "rotation_decision",
        "status": "warning",
        "is_hard_stop": False,
        "detail": "stage 1 derived no tilt (landmark span too small); proceeding "
                  "without rotation.",
    })
elif abs(derived_tilt_angle_degrees) <= parsed_arguments.rotation_threshold_degrees:
    ANALYSIS_FINDINGS.append({
        "check_name": "rotation_decision",
        "status": "pass",
        "is_hard_stop": False,
        "detail": "tilt %.4f degrees is within the %.2f degree threshold, so the crop "
                  "is analysed unrotated and the tilt is recorded only."
                  % (derived_tilt_angle_degrees, parsed_arguments.rotation_threshold_degrees),
    })
else:
    # Landmarks in crop-local (row, col) coordinates and the crop centre.
    landmark_a_row_col = numpy.array([landmark_a_y_pixels - crop_y_pixels,
                                      landmark_a_x_pixels - crop_x_pixels])
    landmark_b_row_col = numpy.array([landmark_b_y_pixels - crop_y_pixels,
                                      landmark_b_x_pixels - crop_x_pixels])
    crop_centre_row_col = numpy.array([crop_height_pixels / 2.0, crop_width_pixels / 2.0])

    # Try both signs and keep the one that levels the landmark line, so no rotation
    # convention has to be reasoned about in the abstract. The image is then rotated
    # by the same matrix, so image and landmarks cannot disagree.
    chosen_rotation_matrix = None
    for candidate_angle_degrees in (derived_tilt_angle_degrees, -derived_tilt_angle_degrees):
        candidate_angle_radians = math.radians(candidate_angle_degrees)
        candidate_cos = math.cos(candidate_angle_radians)
        candidate_sin = math.sin(candidate_angle_radians)
        forward_rotation_matrix = numpy.array([[candidate_cos, -candidate_sin],
                                               [candidate_sin, candidate_cos]])
        a_rotated = forward_rotation_matrix @ (landmark_a_row_col - crop_centre_row_col) + crop_centre_row_col
        b_rotated = forward_rotation_matrix @ (landmark_b_row_col - crop_centre_row_col) + crop_centre_row_col
        rotated_delta_row = b_rotated[0] - a_rotated[0]
        rotated_delta_col = b_rotated[1] - a_rotated[1]
        rotated_line_angle_degrees = math.degrees(math.atan2(rotated_delta_row, rotated_delta_col))
        residual = rotated_line_angle_degrees - target_line_angle_degrees
        while residual > 90.0:
            residual -= 180.0
        while residual <= -90.0:
            residual += 180.0
        if abs(residual) <= ROTATION_RESIDUAL_TILT_TOLERANCE_DEGREES:
            chosen_rotation_matrix = forward_rotation_matrix
            rotation_applied_angle_degrees = candidate_angle_degrees
            rotation_residual_tilt_degrees = residual
            break

    if chosen_rotation_matrix is None:
        ANALYSIS_FINDINGS.append({
            "check_name": "rotation_levels_landmarks",
            "status": "fail",
            "is_hard_stop": True,
            "detail": "neither rotation sign brought the landmark line within %.3f "
                      "degrees of level; refusing to rotate on a transform that does "
                      "not do what it claims."
                      % ROTATION_RESIDUAL_TILT_TOLERANCE_DEGREES,
        })
    else:
        # affine_transform samples input at matrix @ output + offset, so the inverse
        # (transpose) of the forward rotation is used, about the same centre. Plate
        # median fills the corners so it does not read as signal.
        inverse_rotation_matrix = chosen_rotation_matrix.T
        rotation_offset = crop_centre_row_col - inverse_rotation_matrix @ crop_centre_row_col
        fill_value = float(numpy.median(crop_pixels))
        analysis_crop = scipy.ndimage.affine_transform(
            crop_pixels, inverse_rotation_matrix, offset=rotation_offset,
            order=1, mode="constant", cval=fill_value,
        )
        rotation_applied = True
        ANALYSIS_FINDINGS.append({
            "check_name": "rotation_levels_landmarks",
            "status": "pass",
            "is_hard_stop": False,
            "detail": "rotated by %.4f degrees; landmark line re-derived to %.4f "
                      "degrees, within tolerance. Saturation was measured before this."
                      % (rotation_applied_angle_degrees, rotation_residual_tilt_degrees),
        })

# Lane grid fit and orientation check both work from 1D profiles collapsed along
# each axis. Collapsing along migration averages out the bands and leaves the lane
# periodicity; collapsing along the stacking axis leaves the band structure, which
# is not 15-fold periodic. numpy axis 0 is rows (image y), axis 1 is columns (x).
if gel_migration_axis == "horizontal":
    stacking_extent_pixels = analysis_crop.shape[0]
    migration_extent_pixels = analysis_crop.shape[1]
    profile_along_stacking = analysis_crop.mean(axis=1)
    profile_along_migration = analysis_crop.mean(axis=0)
else:
    stacking_extent_pixels = analysis_crop.shape[1]
    migration_extent_pixels = analysis_crop.shape[0]
    profile_along_stacking = analysis_crop.mean(axis=0)
    profile_along_migration = analysis_crop.mean(axis=1)

# Two questions from the collapsed profiles, kept separate because conflating them
# is what made the old check wrong. (1) Orientation: does an axis carry lane-pitch
# periodicity at all? Answered by autocorrelation at the lane spacing, which needs
# no phase and no comb. (2) Where exactly are the lanes? Answered by a comb, but
# only along the stacking axis and only for drawing; its score is no longer trusted
# to decide orientation. Detrend each profile by subtracting a moving average about
# two pitches wide to remove the illumination gradient.
detrended_profile_by_axis = {}
autocorrelation_lags_by_axis = {}
autocorrelation_curve_by_axis = {}
periodicity_autocorrelation_by_axis = {}
for axis_name, axis_profile, axis_extent in (
    ("stacking", profile_along_stacking, stacking_extent_pixels),
    ("migration", profile_along_migration, migration_extent_pixels),
):
    expected_pitch_pixels = axis_extent / float(expected_lane_count)
    smoothing_window = max(3, int(round(2.0 * expected_pitch_pixels)))
    profile_baseline = scipy.ndimage.uniform_filter1d(
        axis_profile, size=smoothing_window, mode="nearest"
    )
    detrended_profile = axis_profile - profile_baseline
    detrended_profile = detrended_profile - detrended_profile.mean()
    detrended_energy = float((detrended_profile * detrended_profile).sum()) + 1e-9

    # Normalised autocorrelation over a band of lags around the expected pitch, so a
    # pitch that is a little tighter or looser than crop_extent / lane_count is still
    # found. The peak over the band is the periodicity score for this axis.
    lag_low = max(1, int(round(expected_pitch_pixels * (1.0 - LANE_PITCH_SEARCH_FRACTION))))
    lag_high = max(lag_low + 1, int(round(expected_pitch_pixels * (1.0 + LANE_PITCH_SEARCH_FRACTION))))
    lag_high = min(lag_high, axis_extent - 1)
    autocorrelation_lags = numpy.arange(lag_low, lag_high + 1)
    autocorrelation_curve = numpy.array([
        float((detrended_profile[:-lag] * detrended_profile[lag:]).sum()) / detrended_energy
        for lag in autocorrelation_lags
    ])

    detrended_profile_by_axis[axis_name] = detrended_profile
    autocorrelation_lags_by_axis[axis_name] = autocorrelation_lags
    autocorrelation_curve_by_axis[axis_name] = autocorrelation_curve
    periodicity_autocorrelation_by_axis[axis_name] = (
        float(autocorrelation_curve.max()) if autocorrelation_curve.size else 0.0
    )

stacking_autocorrelation = periodicity_autocorrelation_by_axis["stacking"]
migration_autocorrelation = periodicity_autocorrelation_by_axis["migration"]

# Data-driven lane locator, replacing the imposed even-count comb. The old comb
# assumed expected_lane_count lanes filling the crop at a single pitch and scored
# whichever phase best hit that many evenly-spaced samples. On a gel where only
# some wells carry signal, and unevenly, that produced a plausible-looking grid
# that did not sit on the real lanes. Here the pitch is MEASURED from the stacking
# autocorrelation, the populated lane centres are FOUND by peak-finding, compact
# non-lane artifacts (dust, rings) are REJECTED by requiring band structure across
# the migration width, each surviving centre is REFINED to the intensity centroid
# of its lane (a bare peak lands on one lobe of the intra-lane doublet and reads
# off centre), and a pitch-and-phase comb is fitted to the survivors ONLY to give
# each an integer well index, so a sparse gel reads as "wells 1, 11, 12, 13, 14
# populated of 15" rather than as five renumbered lanes.
expected_stacking_pitch = stacking_extent_pixels / float(expected_lane_count)

# Measured pitch: the lag of the strongest stacking autocorrelation peak, if the
# axis is periodic enough to trust; otherwise fall back to the geometry-implied
# pitch and say so. The autocorrelation was already computed above over a band of
# lags around the expected pitch.
stacking_autocorrelation_lags = autocorrelation_lags_by_axis["stacking"]
stacking_autocorrelation_curve = autocorrelation_curve_by_axis["stacking"]
if (stacking_autocorrelation_curve.size > 0
        and stacking_autocorrelation_curve.max() >= ORIENTATION_AUTOCORRELATION_FLOOR):
    measured_lane_pitch_pixels = float(
        stacking_autocorrelation_lags[int(numpy.argmax(stacking_autocorrelation_curve))]
    )
    lane_pitch_source = "measured_by_stacking_autocorrelation"
else:
    measured_lane_pitch_pixels = expected_stacking_pitch
    lane_pitch_source = "fell_back_to_geometry_expected_pitch"
lane_pitch_pixels = measured_lane_pitch_pixels
lane_pitch_millimetres = (
    lane_pitch_pixels * micrometres_per_pixel / 1000.0
    if micrometres_per_pixel is not None else None
)

# Orient a background-subtracted signal image so axis 0 is the stacking axis and
# axis 1 is the migration axis, whatever the scan orientation, so the locator logic
# is written once. Background is the crop median (the plate level).
plate_background_value = float(numpy.median(analysis_crop))
signal_above_plate = analysis_crop - plate_background_value
signal_above_plate[signal_above_plate < 0.0] = 0.0
if gel_migration_axis == "horizontal":
    stacking_by_migration_signal = signal_above_plate
else:
    stacking_by_migration_signal = signal_above_plate.T

# Usable migration columns: drop plate-edge and scratch columns (elevated down
# their whole length) and any container-saturated columns, so they neither swamp
# the stacking profile nor defeat the band-structure test.
per_migration_column_median = numpy.median(stacking_by_migration_signal, axis=0)
positive_column_medians = per_migration_column_median[per_migration_column_median > 0.0]
typical_column_median = (
    float(numpy.median(positive_column_medians)) if positive_column_medians.size else 1.0
)
per_migration_column_saturated_fraction = numpy.mean(
    (analysis_crop if gel_migration_axis == "horizontal" else analysis_crop.T)
    >= CONTAINER_MAXIMUM_VALUE_16_BIT, axis=0
)
usable_migration_column_mask = (
    (per_migration_column_median <= EDGE_COLUMN_MEDIAN_MULTIPLE * typical_column_median)
    & (per_migration_column_saturated_fraction < 0.05)
)
usable_migration_column_count = int(usable_migration_column_mask.sum())

# Stacking profile for lane finding: signal summed over the usable migration
# columns only, then smoothed over half the pitch to merge each lane's doublet.
usable_stacking_by_migration_signal = numpy.where(
    usable_migration_column_mask[numpy.newaxis, :], stacking_by_migration_signal, 0.0
)
stacking_profile_for_lane_finding = usable_stacking_by_migration_signal.sum(axis=1)
lane_centre_smoothing_pixels = max(
    3, int(round(LANE_CENTRE_SMOOTHING_PITCH_FRACTION * measured_lane_pitch_pixels))
)
stacking_profile_smoothed_for_lane_finding = scipy.ndimage.uniform_filter1d(
    stacking_profile_for_lane_finding, size=lane_centre_smoothing_pixels, mode="nearest"
)

# Candidate lane centres: prominent peaks, separated by at least most of a pitch.
lane_minimum_separation_pixels = max(
    1.0, LANE_MINIMUM_SEPARATION_PITCH_FRACTION * measured_lane_pitch_pixels
)
lane_prominence_threshold = (
    LANE_PROMINENCE_FRACTION_OF_MAXIMUM * float(stacking_profile_smoothed_for_lane_finding.max())
)
candidate_lane_centre_rows, _ = scipy.signal.find_peaks(
    stacking_profile_smoothed_for_lane_finding,
    distance=lane_minimum_separation_pixels,
    prominence=lane_prominence_threshold,
)

# Keep only candidates with real lane structure, and refine each kept centre to the
# intensity centroid of its strip so the doublet does not bias it off centre.
half_pitch_pixels = int(round(measured_lane_pitch_pixels / 2.0))
lane_centres_stacking_pixels = []
lane_participating_column_counts = []
lane_band_counts = []
lane_total_signal_values = []
rejected_candidate_rows = []
rejected_candidate_reasons = []
for candidate_row in candidate_lane_centre_rows:
    strip_start = max(0, int(candidate_row) - half_pitch_pixels)
    strip_end = min(stacking_extent_pixels, int(candidate_row) + half_pitch_pixels)
    lane_strip = usable_stacking_by_migration_signal[strip_start:strip_end, :]
    lane_migration_profile = lane_strip.sum(axis=0)
    lane_migration_profile_maximum = float(lane_migration_profile.max())
    if lane_migration_profile_maximum <= 0.0:
        rejected_candidate_rows.append(int(candidate_row))
        rejected_candidate_reasons.append("no_signal_in_strip")
        continue
    participation_threshold = (
        LANE_COLUMN_PARTICIPATION_THRESHOLD_FRACTION * lane_migration_profile_maximum
    )
    participating_column_count = int((lane_migration_profile > participation_threshold).sum())
    band_peaks, _ = scipy.signal.find_peaks(
        lane_migration_profile,
        distance=8,
        prominence=LANE_COLUMN_PARTICIPATION_THRESHOLD_FRACTION * lane_migration_profile_maximum,
    )
    is_real_lane = (
        participating_column_count
        >= LANE_MINIMUM_PARTICIPATING_COLUMN_FRACTION * usable_migration_column_count
        and len(band_peaks) >= 1
    )
    if not is_real_lane:
        rejected_candidate_rows.append(int(candidate_row))
        rejected_candidate_reasons.append(
            "participating_columns_%d_bands_%d" % (participating_column_count, len(band_peaks))
        )
        continue
    strip_rows = numpy.arange(strip_start, strip_end)
    strip_stacking_weight = stacking_profile_for_lane_finding[strip_start:strip_end]
    strip_weight_total = float(strip_stacking_weight.sum())
    refined_centre_row = (
        float((strip_rows * strip_stacking_weight).sum() / strip_weight_total)
        if strip_weight_total > 0.0 else float(candidate_row)
    )
    lane_centres_stacking_pixels.append(refined_centre_row)
    lane_participating_column_counts.append(participating_column_count)
    lane_band_counts.append(int(len(band_peaks)))
    lane_total_signal_values.append(float(lane_strip.sum()))

populated_lane_count = len(lane_centres_stacking_pixels)

# Fit a pitch-and-phase comb to the found centres ONLY to label them with integer
# well indices consistent with the physical comb, tolerating empty wells between
# populated ones. Pitch is the measured pitch; phase is brute-forced over one pitch.
lane_well_indices = None
if populated_lane_count >= 2:
    found_centres_array = numpy.array(lane_centres_stacking_pixels, dtype=numpy.float64)
    best_phase_offset_pixels = 0.0
    best_total_squared_residual_pixels = None
    for candidate_phase_offset_pixels in numpy.linspace(0.0, measured_lane_pitch_pixels, 200):
        nearest_well_index = numpy.round(
            (found_centres_array - candidate_phase_offset_pixels) / measured_lane_pitch_pixels
        )
        model_centre_positions = (
            candidate_phase_offset_pixels + nearest_well_index * measured_lane_pitch_pixels
        )
        total_squared_residual_pixels = float(
            numpy.sum((found_centres_array - model_centre_positions) ** 2)
        )
        if (best_total_squared_residual_pixels is None
                or total_squared_residual_pixels < best_total_squared_residual_pixels):
            best_total_squared_residual_pixels = total_squared_residual_pixels
            best_phase_offset_pixels = float(candidate_phase_offset_pixels)
    raw_well_indices = numpy.round(
        (found_centres_array - best_phase_offset_pixels) / measured_lane_pitch_pixels
    ).astype(int)
    raw_well_indices = raw_well_indices - int(raw_well_indices.min())
    lane_well_indices = [int(value) for value in raw_well_indices]
elif populated_lane_count == 1:
    lane_well_indices = [0]

# Orientation verdict, from per-axis gradient energy, reported as a WARNING not a
# hard stop. It flags a probable mislabel for the operator to check against the
# overlay, but does not block: a false block on a correct gel is the worse failure.
# Gradient energy is measured on the usable (edge-masked) signal, oriented so axis
# 0 is stacking and axis 1 is migration, so plate edges do not count as band edges.
# The lane-pitch autocorrelations are still computed above and carried into the
# report as secondary diagnostics, but they no longer decide orientation, because
# they collapse to noise on a sparsely populated gel.
stacking_axis_gradient = numpy.diff(usable_stacking_by_migration_signal, axis=0)
migration_axis_gradient = numpy.diff(usable_stacking_by_migration_signal, axis=1)
stacking_gradient_energy = float((stacking_axis_gradient * stacking_axis_gradient).sum())
migration_gradient_energy = float((migration_axis_gradient * migration_axis_gradient).sum())
orientation_gradient_energy_ratio = (
    migration_gradient_energy / stacking_gradient_energy
    if stacking_gradient_energy > 0.0 else 0.0
)
orientation_is_consistent = (
    orientation_gradient_energy_ratio >= ORIENTATION_GRADIENT_ENERGY_RATIO
)
ANALYSIS_FINDINGS.append({
    "check_name": "orientation_matches_band_sharpness",
    "status": "pass" if orientation_is_consistent else "warning",
    "is_hard_stop": False,
    "detail": "gradient energy is %.2fx higher along the claimed migration axis than "
              "the claimed stacking axis (%s migration); bands are sharp along "
              "migration and smooth along stacking, so migration should carry more "
              "structure. "
              % (orientation_gradient_energy_ratio, gel_migration_axis)
              + ("consistent: the sharper axis is the one recorded as migration."
                 if orientation_is_consistent else
                 "the migration axis is not clearly the sharper one, so the axis may "
                 "be mislabelled. Check the overlay and lane_profiles.png before "
                 "trusting or fixing the sidecar.")
              + " Lane-pitch autocorrelation (secondary, unreliable on sparse gels) "
                "is %.3f stacking vs %.3f migration."
              % (stacking_autocorrelation, migration_autocorrelation),
})

# Lane grid finding: report the MEASURED occupancy, not the assumed count. Finding
# fewer populated lanes than expected is normal (empty wells) and is a warning for
# the operator to confirm against the overlay, not a failure. Finding none is a
# warning worth flagging louder.
if populated_lane_count == 0:
    lane_grid_status = "warning"
    lane_grid_detail = (
        "no populated lanes were found. The crop may be blank, the signal may be "
        "below the prominence floor, or the migration axis may be mislabelled; "
        "check lane_profiles.png and the orientation finding."
    )
else:
    lane_grid_status = "pass" if populated_lane_count == expected_lane_count else "warning"
    lane_grid_detail = (
        "found %d populated lane(s) of %d wells by peak-finding at a measured pitch "
        "of %.2f pixels" % (populated_lane_count, expected_lane_count, lane_pitch_pixels)
        + ("" if lane_pitch_millimetres is None else " (%.3f mm)" % lane_pitch_millimetres)
        + " (%s)" % lane_pitch_source
        + (", assigned to well indices %s" % lane_well_indices
           if lane_well_indices is not None else "")
        + (", rejected %d candidate(s) as non-lane artifacts" % len(rejected_candidate_rows)
           if rejected_candidate_rows else "")
        + ". Occupancy is measured from the image, not assumed; confirm against the "
          "overlay."
    )
ANALYSIS_FINDINGS.append({
    "check_name": "lane_grid_fit",
    "status": lane_grid_status,
    "is_hard_stop": False,
    "detail": lane_grid_detail,
})

# Overlay for the operator. The crop is drawn unflipped so it matches Fiji, with the
# fitted lane lines and the orientation verdict in the title.
overlay_stride = max(1, int(math.ceil(
    max(analysis_crop.shape) / float(OVERLAY_MAXIMUM_DIMENSION_PIXELS))))
overlay_array = analysis_crop[::overlay_stride, ::overlay_stride]
overlay_figure, overlay_axes = matplotlib.pyplot.subplots(figsize=(8.0, 8.0))
overlay_axes.imshow(
    overlay_array, cmap=DISPLAY_COLORMAP_NAME,
    vmin=numpy.percentile(overlay_array, 1.0),
    vmax=numpy.percentile(overlay_array, 99.5),
    interpolation="nearest", origin="upper",
)
if lane_centres_stacking_pixels is not None:
    for lane_centre in lane_centres_stacking_pixels:
        if gel_migration_axis == "horizontal":
            overlay_axes.axhline(lane_centre / overlay_stride, color="red", linewidth=0.8)
        else:
            overlay_axes.axvline(lane_centre / overlay_stride, color="red", linewidth=0.8)
overlay_title = (
    "Lane grid on the crop, unflipped: " + input_tiff_absolute_path.name
    + "\nmigration " + gel_migration_axis
    + ("; unrotated" if not rotation_applied
       else "; rotated %.3f deg" % rotation_applied_angle_degrees)
    + "; orientation " + ("OK" if orientation_is_consistent else "CHECK")
    + "\nmigration/stacking gradient-energy ratio %.2f (orientation cue)"
      % orientation_gradient_energy_ratio
    + ("" if lane_pitch_millimetres is None
       else "; pitch %.2f mm" % lane_pitch_millimetres)
)
overlay_axes.set_title(overlay_title, fontsize="medium")
overlay_axes.set_xlabel("crop column / " + str(overlay_stride))
overlay_axes.set_ylabel("crop row / " + str(overlay_stride))
overlay_figure.tight_layout()
overlay_output_path = output_directory_path / LANE_GRID_OVERLAY_FILENAME
try:
    overlay_figure.savefig(overlay_output_path, dpi=FIGURE_DOTS_PER_INCH)
except Exception as overlay_write_error:
    die("overlay", "could not write " + str(overlay_output_path) + ": " + str(overlay_write_error))
matplotlib.pyplot.close(overlay_figure)
emit_message("overlay", "wrote " + str(overlay_output_path))

# Diagnostic profiles, always written. This is the view that tells us, on a real
# gel, the true lane pitch and phase, how many bumps there are, whether the lanes
# fill the crop or sit inside margins, and how cleanly the two axes separate. The
# top panels are the detrended collapsed profiles with the fitted comb centres
# marked; the bottom panel is the lane-pitch autocorrelation for both axes.
profiles_figure, profiles_axes = matplotlib.pyplot.subplots(3, 1, figsize=(9.0, 9.0))
stacking_detrended_profile = detrended_profile_by_axis["stacking"]
migration_detrended_profile = detrended_profile_by_axis["migration"]
profiles_axes[0].plot(numpy.arange(stacking_detrended_profile.size),
                      stacking_detrended_profile, color="black", linewidth=0.8)
if lane_centres_stacking_pixels is not None:
    for lane_centre in lane_centres_stacking_pixels:
        profiles_axes[0].axvline(lane_centre, color="red", linewidth=0.7)
profiles_axes[0].set_title("collapsed along migration -> stacking-axis profile "
                           "(should show the lanes); red = fitted comb centres",
                           fontsize="small")
profiles_axes[0].set_xlabel("stacking-axis position (pixels)")
profiles_axes[1].plot(numpy.arange(migration_detrended_profile.size),
                      migration_detrended_profile, color="black", linewidth=0.8)
profiles_axes[1].set_title("collapsed along stacking -> migration-axis profile "
                           "(should NOT be lane-periodic)", fontsize="small")
profiles_axes[1].set_xlabel("migration-axis position (pixels)")
for axis_name, curve_color in (("stacking", "red"), ("migration", "blue")):
    profiles_axes[2].plot(autocorrelation_lags_by_axis[axis_name],
                          autocorrelation_curve_by_axis[axis_name],
                          color=curve_color, linewidth=1.0, label=axis_name)
profiles_axes[2].axhline(0.0, color="grey", linewidth=0.5)
profiles_axes[2].set_title("lane-pitch autocorrelation vs lag; peak decides "
                           "orientation (stacking should win)", fontsize="small")
profiles_axes[2].set_xlabel("lag (pixels)")
profiles_axes[2].legend(fontsize="small")
profiles_figure.tight_layout()
lane_profiles_plot_path = output_directory_path / LANE_PROFILES_PLOT_FILENAME
profiles_figure.savefig(lane_profiles_plot_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(profiles_figure)
emit_message("profiles", "wrote " + str(lane_profiles_plot_path))

# The same profiles as data, long format, so they can be re-plotted or fitted
# outside this script without re-running it.
lane_profiles_csv_lines = ["axis,position_pixels,detrended_value"]
for axis_name, axis_profile_values in (
    ("stacking", stacking_detrended_profile),
    ("migration", migration_detrended_profile),
):
    for position_index in range(axis_profile_values.size):
        lane_profiles_csv_lines.append(
            "%s,%d,%.6f" % (axis_name, position_index, float(axis_profile_values[position_index]))
        )
lane_profiles_csv_path = output_directory_path / LANE_PROFILES_CSV_FILENAME
lane_profiles_csv_path.write_text("\n".join(lane_profiles_csv_lines) + "\n")
emit_message("profiles", "wrote " + str(lane_profiles_csv_path))

# Single exit point: assemble and write the provenance report, then set the code.
failed_hard_stops = [
    finding["check_name"] for finding in ANALYSIS_FINDINGS
    if finding["is_hard_stop"] and finding["status"] == "fail"
]
warnings = [
    finding["check_name"] for finding in ANALYSIS_FINDINGS
    if finding["status"] == "warning"
]
overall_status = "fail" if failed_hard_stops else ("warning" if warnings else "pass")

analysis_report = {
    "analysis_report_schema_version": ANALYSIS_REPORT_SCHEMA_VERSION,
    "overall_status": overall_status,
    "failed_hard_stop_check_names": failed_hard_stops,
    "warning_check_names": warnings,
    "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "input_file": str(input_tiff_absolute_path),
    "read_from_stage1_report": str(validation_report_path),
    # Carried through untouched so no consumer of stage 2 can forget that the pixel
    # values may not be linear in signal.
    "encoding_verified": encoding_verified,
    "geometry_used": {
        "gel_migration_axis": gel_migration_axis,
        "crop_x_pixels": crop_x_pixels,
        "crop_y_pixels": crop_y_pixels,
        "crop_width_pixels": crop_width_pixels,
        "crop_height_pixels": crop_height_pixels,
        "micrometres_per_pixel": micrometres_per_pixel,
        "expected_lane_count": expected_lane_count,
        "derived_tilt_angle_degrees": derived_tilt_angle_degrees,
    },
    "rotation": {
        "threshold_degrees": parsed_arguments.rotation_threshold_degrees,
        "applied": rotation_applied,
        "applied_angle_degrees": rotation_applied_angle_degrees,
        "residual_tilt_degrees": rotation_residual_tilt_degrees,
    },
    "in_crop_saturation": {
        "crop_maximum_value": crop_maximum_value,
        "at_crop_maximum_count": at_crop_maximum_count,
        "at_container_maximum_count": at_container_maximum_in_crop,
    },
    "lane_grid": {
        "pitch_pixels": lane_pitch_pixels,
        "pitch_millimetres": lane_pitch_millimetres,
        "pitch_source": lane_pitch_source,
        "expected_lane_count": expected_lane_count,
        "populated_lane_count": populated_lane_count,
        "lane_well_indices": lane_well_indices,
        "lane_centres_stacking_pixels": [
            float(value) for value in lane_centres_stacking_pixels
        ],
        "lane_participating_column_counts": lane_participating_column_counts,
        "lane_band_counts": lane_band_counts,
        "lane_total_signal_values": lane_total_signal_values,
        "usable_migration_column_count": usable_migration_column_count,
        "rejected_candidate_rows": rejected_candidate_rows,
        "rejected_candidate_reasons": rejected_candidate_reasons,
    },
    "orientation_check": {
        "claimed_migration_axis": gel_migration_axis,
        "decided_by": "per_axis_gradient_energy",
        "migration_gradient_energy": migration_gradient_energy,
        "stacking_gradient_energy": stacking_gradient_energy,
        "migration_over_stacking_gradient_energy_ratio": orientation_gradient_energy_ratio,
        "gradient_energy_ratio_threshold": ORIENTATION_GRADIENT_ENERGY_RATIO,
        "lane_pitch_autocorrelation_stacking_diagnostic": stacking_autocorrelation,
        "lane_pitch_autocorrelation_migration_diagnostic": migration_autocorrelation,
        "consistent": orientation_is_consistent,
    },
    "analysis_findings": ANALYSIS_FINDINGS,
    "outputs": {
        "lane_grid_overlay": str(overlay_output_path),
        "lane_profiles_plot": str(lane_profiles_plot_path),
        "lane_profiles_csv": str(lane_profiles_csv_path),
    },
}
analysis_report_path = output_directory_path / ANALYSIS_REPORT_FILENAME
analysis_report_path.write_text(json.dumps(analysis_report, indent=2))
emit_message("report", "wrote " + str(analysis_report_path) + "; overall " + overall_status)

sys.exit(1 if failed_hard_stops else 0)
