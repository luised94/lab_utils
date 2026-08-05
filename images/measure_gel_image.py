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

# The comb, along the stacking axis only, to place the lane centres for the overlay.
# This is the part that assumes a fixed count filling the crop, and is the next
# thing to make robust once the profiles show how the real lanes sit.
stacking_detrended = detrended_profile_by_axis["stacking"]
stacking_scale = float(stacking_detrended.std()) + 1e-9
expected_stacking_pitch = stacking_extent_pixels / float(expected_lane_count)
best_comb_score = -1e18
lane_centres_stacking_pixels = None
lane_pitch_pixels = expected_stacking_pitch
pitch_low = expected_stacking_pitch * (1.0 - LANE_PITCH_SEARCH_FRACTION)
pitch_high = expected_stacking_pitch * (1.0 + LANE_PITCH_SEARCH_FRACTION)
for candidate_pitch in numpy.linspace(pitch_low, pitch_high, 61):
    first_centre_step = max(1.0, candidate_pitch / 20.0)
    first_centre = 0.0
    while first_centre < candidate_pitch:
        centres = first_centre + candidate_pitch * numpy.arange(expected_lane_count)
        sample_indices = numpy.rint(centres).astype(int)
        if (sample_indices >= 0).all() and (sample_indices < stacking_extent_pixels).all():
            score = float(stacking_detrended[sample_indices].mean()) / stacking_scale
            if score > best_comb_score:
                best_comb_score = score
                lane_centres_stacking_pixels = centres
                lane_pitch_pixels = candidate_pitch
        first_centre += first_centre_step
lane_pitch_millimetres = (
    lane_pitch_pixels * micrometres_per_pixel / 1000.0
    if micrometres_per_pixel is not None else None
)

# Orientation verdict, from autocorrelation, reported as a WARNING not a hard stop.
# It flags a probable mislabel for the operator to check against the overlay, but
# does not block, because the statistic is not yet proven on real gels and a false
# block on a correct gel is the worse failure.
orientation_is_consistent = (
    stacking_autocorrelation >= ORIENTATION_AUTOCORRELATION_FLOOR
    and stacking_autocorrelation >= ORIENTATION_PERIODICITY_RATIO * migration_autocorrelation
)
ANALYSIS_FINDINGS.append({
    "check_name": "orientation_matches_lane_periodicity",
    "status": "pass" if orientation_is_consistent else "warning",
    "is_hard_stop": False,
    "detail": "lane-pitch autocorrelation is %.3f along the claimed stacking axis and "
              "%.3f along the claimed migration axis (%s migration). "
              % (stacking_autocorrelation, migration_autocorrelation, gel_migration_axis)
              + ("consistent: the lanes stack across the axis recorded."
                 if orientation_is_consistent else
                 "the stacking axis is not clearly the more periodic one. Either the "
                 "axis is mislabelled, or the lanes are faint or unevenly spaced and "
                 "the periodicity is weak. Check the overlay and lane_profiles.png "
                 "before trusting or fixing the sidecar."),
})

lane_count_found_plausible = (
    lane_centres_stacking_pixels is not None
    and len(lane_centres_stacking_pixels) == expected_lane_count
)
ANALYSIS_FINDINGS.append({
    "check_name": "lane_grid_fit",
    "status": "pass" if lane_count_found_plausible else "warning",
    "is_hard_stop": False,
    "detail": ("placed %d evenly-spaced lane centres at a pitch of %.2f pixels"
               % (expected_lane_count, lane_pitch_pixels)
               + ("" if lane_pitch_millimetres is None
                  else " (%.3f mm)" % lane_pitch_millimetres)
               + ". Note the count is assumed from expected_lane_count, not counted "
                 "from peaks; confirm against the overlay.")
              if lane_count_found_plausible else
              "could not place the expected number of lane centres inside the crop.",
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
    + "\nlane-pitch autocorr stacking %.2f vs migration %.2f"
      % (stacking_autocorrelation, migration_autocorrelation)
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
        "lane_count": expected_lane_count,
        "lane_centres_stacking_pixels": (
            [float(value) for value in lane_centres_stacking_pixels]
            if lane_centres_stacking_pixels is not None else None
        ),
    },
    "orientation_check": {
        "claimed_migration_axis": gel_migration_axis,
        "lane_pitch_autocorrelation_stacking": stacking_autocorrelation,
        "lane_pitch_autocorrelation_migration": migration_autocorrelation,
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
