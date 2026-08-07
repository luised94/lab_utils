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
Stage 2 of the gel densitometry pipeline: the single measurement stage. Reads one
already-validated TIFF plus the stage 1 report, crops to the sidecar geometry,
fits the lane grid, checks orientation, detects bands per lane, builds consensus
band windows across lanes, integrates every cell above a local baseline, and
writes one provenance report (band_detection_report.json) plus the overlays,
profiles, and CSVs.

This file is the union of the former Slice A (measure_gel_image.py: crop, lane
grid, orientation) and Slice B (measure_gel_bands.py: per-lane detection,
consensus, integration, figures). They were separate scripts that handed geometry
between them through stage2_analysis_report.json on disk; that handoff is gone.
The image is read once and cropped once, and the band stage consumes the lane
geometry from memory, so the two halves cannot disagree about the crop.

Tilt is measured but never applied. Stage 1 derives the tilt angle from the two
landmark clicks and this stage reports it and warns above a threshold; the gel is
run essentially straight on the imager, and the rare tilted gel is straightened or
rotated in Fiji before cropping. Removing the in-line resampler removes the
interpolation blur and the saturation-ordering constraint that a rotation would
have imposed: with no interpolation, the saturation measurement no longer has to
precede anything, and the crop the locator sees is the crop the integrator sees.

The pixel values may not be linear in signal until a dilution series is measured;
encoding_verified is carried from stage 1 into this report so no consumer can
forget that every quantitative number here is provisional.
"""

import argparse
import datetime
import json
import math
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

# The merged report carries its own top-level schema version, and retains the two
# section schema versions from the source reports so each section still declares
# the schema it was written against.
GEL_MEASUREMENT_REPORT_SCHEMA_VERSION = 1
ANALYSIS_REPORT_SCHEMA_VERSION = 1
BAND_DETECTION_REPORT_SCHEMA_VERSION = 1

# Filenames inherited from Slice A. The unified report is written to Slice B's
# BAND_DETECTION_REPORT_FILENAME (defined below); Slice A's stage2_analysis_report
# and Slice B's stage2-report constants are gone with the disk handoff.
VALIDATION_REPORT_FILENAME = "input_file_validation_report.json"
LANE_GRID_OVERLAY_FILENAME = "lane_grid_overlay.png"
LANE_PROFILES_PLOT_FILENAME = "lane_profiles.png"
LANE_PROFILES_CSV_FILENAME = "lane_profiles.csv"

# Tilt above this many degrees raises a warning (straighten and re-scan, or rotate
# in Fiji before cropping); the gel is measured un-rotated either way. This was
# Slice A's rotation threshold, now a warning-only threshold since nothing rotates.
DEFAULT_TILT_WARNING_THRESHOLD_DEGREES = 0.25

BAND_DETECTION_REPORT_FILENAME = "band_detection_report.json"
BAND_WINDOWS_OVERLAY_FILENAME = "overlay_band_windows_labeled.png"
LANE_MIGRATION_PROFILES_FILENAME = "lane_migration_profiles.png"
BAND_CENTERS_CSV_FILENAME = "band_centers.csv"
# Slice B2 outputs.
BAND_DEFINITIONS_CSV_FILENAME = "band_definitions.csv"
BAND_MEASUREMENTS_CSV_FILENAME = "band_measurements.csv"
CONSENSUS_CENTER_MAP_FILENAME = "consensus_center_map.png"
CONSENSUS_WINDOWS_OVERLAY_FILENAME = "overlay_consensus_windows_by_occupancy.png"
# Slice B3 output.
SATURATION_DERIVATION_FILENAME = "saturation_derivation.png"
BASELINE_COMPARISON_FILENAME = "baseline_comparison.png"
STRIP_TRACE_FILENAME = "figure_strip_trace.png"
METRIC_COMPARE_FILENAME = "figure_metric_compare.png"
NET_VS_APEX_FILENAME = "figure_net_vs_apex.png"
REPORTED_VALUES_OVERLAY_FILENAME = "overlay_reported_values.png"

CONTAINER_MAXIMUM_VALUE_16_BIT = 65535

# Slice B3 saturation, option 1: the effective ceiling is DERIVED per image, not
# baked, because the clip level is a sensor/processing property that moves between
# scans. On s0002 (Amersham Imager 680, Auto High-dynamic-range exposure) the high
# tail decays smoothly to the 16-bit container with no sub-ceiling pile-up, so the
# effective ceiling IS the container and the measured crop has zero saturated
# pixels; the derivation must report that honestly rather than assert a plateau
# that is not there. On a differently-exposed gel the sensor may clip below the
# container, which shows as an over-populated DN with almost nothing above it, and
# then that DN is the effective ceiling. The knee search below detects that case
# and otherwise defaults to the container.
#
# A sub-ceiling knee is the LOWEST DN, at or above this fraction of the container,
# whose pixel count is a strong spike over the local tail baseline just under it
# AND above which essentially no pixels exist. The search only looks below the
# container; a pile-up AT the container is the container clip, not a sensor knee.
SUB_CEILING_KNEE_SEARCH_FLOOR_FRACTION = 0.6
SUB_CEILING_KNEE_SPIKE_FACTOR = 8.0
SUB_CEILING_KNEE_LOCAL_BASELINE_SPAN_VALUES = 200
SUB_CEILING_KNEE_ABOVE_TOLERANCE_FRACTION = 0.001

# A cell whose window-by-strip patch has more than this fraction of pixels at or
# above the effective ceiling is flagged: its integrated intensity is a floor, not
# a measure, because clipped pixels lost their true height. A handful of hot or
# dust pixels should not condemn a cell, so this is a small fraction, not any-pixel.
DEFAULT_SATURATION_WARN_FRACTION = 0.005

# Slice B3.2 rolling-ball cross-check. Width pinned to 4.0 mm on the s0002 profiles:
# it exceeds the doublet (~40 px) and single-band support (15-40 px) so an isolated
# band sits fully above the opened baseline, but stays below the ~5 mm point where
# the plate begins bridging genuinely separate close bands into one. Converted
# through the report pixel size so a differently-scanned gel does not silently shift.
DEFAULT_ROLLING_BALL_WIDTH_MILLIMETRES = 4.0
DEFAULT_STERNBERG_BALL_RADIUS_MILLIMETRES = 4.0
DEFAULT_STERNBERG_BALL_VERTICAL_SCALE = 1.0

# A cell is baseline-fragile when the two baselines' nets differ by more than this
# fraction of the lane's signal scale (largest opened net in the lane). Normalising
# by the lane scale, not by the cell's own net, is what stops near-empty cells from
# reading as maximally fragile.
BASELINE_FRAGILE_DISAGREEMENT_FRACTION = 0.15
# A lane whose signal scale is below this fraction of the strongest lane's scale is
# treated as having no signal, so its cells read no_signal rather than fragile.
LANE_SIGNAL_SCALE_NO_SIGNAL_FRACTION = 0.05
# Width-adequacy: if even the best-case lane's opened baseline still rides this far
# up into its tallest band (as a fraction of peak), the width is too small globally.
WIDTH_ADEQUACY_CLIMB_WARN_FRACTION = 0.15
# A band whose opened baseline sits this far above the lane's opened-baseline floor
# (as a fraction of the lane's peak-above-floor) is being bridged by a neighbour, so
# its opened net is a cluster-shoulder quantity, not a clean per-band background.
OPENING_BASELINE_ELEVATED_FRACTION = 0.10

# A migration column elevated down its whole length is a plate edge or scratch,
# not a lane; exclude columns whose median exceeds this multiple of the typical
# column median. Used by both the lane locator above and the band-integration
# column mask below, one definition now that both run in the same script.
EDGE_COLUMN_MEDIAN_MULTIPLE = 8.0

# Config, lengths in millimetres per DESIGN section 4, converted to pixels using
# the pixel size read from the Slice A report so a differently-scanned gel does not
# silently change the analysis. All overridable on the command line.
DEFAULT_REGION_WINDOW_WIDTH_MILLIMETRES = 2.0
# Minimum separation between called band centres. Set from the s0002 measurement:
# the mid Mcm2-7 cluster band sits ~1.4 mm from its neighbour with ~0.45 relative
# prominence, so a 1.5 mm floor dropped a real, strong band; 1.0 mm recovers it and
# the two other close-but-real peaks, and 0.8 mm gives an identical result, so 1.0
# is a stable choice rather than a knife-edge. Peaks closer than this that are also
# faint are still gated by prominence.
DEFAULT_MINIMUM_BAND_SEPARATION_MILLIMETRES = 1.0
DEFAULT_MIGRATION_PROFILE_SMOOTHING_MILLIMETRES = 0.4
DEFAULT_LANE_STRIP_FIXED_HALF_HEIGHT_MILLIMETRES = 3.5

# Dimensionless. The band peak must clear this fraction of the lane's own profile
# maximum; kept as a fraction of the lane's own signal so a faint lane is judged on
# its own scale, not the strongest lane's.
DEFAULT_BAND_PEAK_PROMINENCE_FRACTION = 0.12
# The lane strip runs to where the lane's stacking cross-section falls to this
# fraction of its own maximum, capturing the doublet.
DEFAULT_LANE_STRIP_CROSS_SECTION_FRACTION = 0.25

# Slice B2 consensus and integration config.
# The cluster tolerance is NOT derived from a gap histogram: on real gels the
# pooled adjacent-gap distribution is not reliably bimodal (verified on s0002,
# where it is a continuous ramp with no clean valley), so a valley-based tolerance
# over-splits. Instead the tolerance defaults to the minimum band separation: two
# centres from DIFFERENT lanes closer than the tightest within-lane band spacing
# are the same band shifted by cross-lane misalignment. The same-lane exclusion in
# the clustering (never two centres from one lane in one canonical band) is the
# hard constraint that keeps this from over-merging, so the tolerance is a bounded
# choice rather than a delicate knob. Overridable on the command line.
# A canonical band whose cross-lane spread reaches this fraction of the minimum
# separation is flagged: its members may not be the same species, and the merge is
# uncertain. It is a warning with a suggested action, not a hard stop.
CANONICAL_SPREAD_WARN_FRACTION = 0.8
# Local valley-to-valley baseline: the flanking-gel minima are searched this far
# outside each window; where windows abut (clipped at a midpoint) the search runs
# into the shared valley, so neighbouring bands sit on a consistent baseline.
DEFAULT_BASELINE_FLANK_SEARCH_MILLIMETRES = 1.0

DISPLAY_COLORMAP_NAME = "gray_r"
FIGURE_DOTS_PER_INCH = 130


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
ORIENTATION_PERIODICITY_RATIO = 1.25  # inherited from Slice A; defined but unused
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

# The lane-grid overlay is downsampled so its longest side is at most this many
# pixels, to keep the written PNG a reasonable size.
OVERLAY_MAXIMUM_DIMENSION_PIXELS = 1400


def emit_message(source_tag, message_text):
    sys.stderr.write("[" + source_tag + "] " + message_text + "\n")


def die(source_tag, message_text):
    emit_message(source_tag, message_text)
    sys.exit(2)


def rolling_ball_opening_baseline(profile_values, width_pixels):
    """Flat morphological opening as the DESIGN 5.4 rolling-ball cross-check.

    A flat structuring element of the given odd width is slid beneath the 1-D
    migration profile (grey erosion then dilation). This is the DEFAULT cross-check
    because it has a single, unambiguous parameter (the width) and no vertical
    intensity-to-pixel scaling to choose. The opened baseline lies at or below the
    profile by construction, so net-above-it is always non-negative and can never
    produce the baseline-above-signal clip that valley-to-valley does. The width
    must exceed the widest single-band feature or the plate rides up into bands;
    that is the width-adequacy check downstream.
    """
    odd_width = max(1, int(round(width_pixels)))
    if odd_width % 2 == 0:
        odd_width += 1
    return scipy.ndimage.grey_opening(profile_values, size=odd_width, mode="nearest")


def sternberg_ball_baseline(profile_values, radius_pixels, vertical_scale):
    """Sternberg rolling ball, KEPT AS A SELECTABLE REFERENCE, not the default.

    A parabolic (ball-cap) element is rolled beneath the profile. Unlike the flat
    opening this element has a vertical extent, so it needs a scale relating
    intensity units to pixel units; that scale is an unpinned degree of freedom
    that drives the baseline as much as the radius does (this is exactly why the
    flat opening is the default). It is retained because DESIGN 5.4 names the
    rolling ball and the comparison is a useful reference; treat vertical_scale as
    experimental and confirm any run against the baseline_comparison panel.
    """
    ball_radius = max(1, int(round(radius_pixels)))
    offsets = numpy.arange(-ball_radius, ball_radius + 1)
    # Concave-down cap; curvature set by the radius, height set by vertical_scale.
    element_cap = vertical_scale * (
        numpy.sqrt(numpy.maximum(0.0, ball_radius**2 - offsets**2)) - ball_radius
    )
    padded = numpy.pad(profile_values, ball_radius, mode="edge")
    eroded = numpy.empty_like(profile_values, dtype=float)
    for index in range(profile_values.size):
        eroded[index] = numpy.min(
            padded[index : index + 2 * ball_radius + 1] - element_cap
        )
    padded_eroded = numpy.pad(eroded, ball_radius, mode="edge")
    opened = numpy.empty_like(profile_values, dtype=float)
    for index in range(profile_values.size):
        opened[index] = numpy.max(
            padded_eroded[index : index + 2 * ball_radius + 1] + element_cap
        )
    return opened


def cross_check_baseline(
    profile_values, method_name, width_pixels, radius_pixels, vertical_scale
):
    """Dispatch to the selected cross-check baseline. Opening is the default."""
    if method_name == "sternberg_ball":
        return sternberg_ball_baseline(profile_values, radius_pixels, vertical_scale)
    return rolling_ball_opening_baseline(profile_values, width_pixels)


argument_parser = argparse.ArgumentParser(
    description="Stage 3 Slice B1: per-lane band-centre detection and windows, "
    "drawn and dumped. Reads the Slice A geometry report. No integration."
)
argument_parser.add_argument("input_tiff_path")
argument_parser.add_argument("--page-index", type=int, default=0)
argument_parser.add_argument(
    "--region-window-width-millimetres",
    type=float,
    default=DEFAULT_REGION_WINDOW_WIDTH_MILLIMETRES,
)
argument_parser.add_argument(
    "--minimum-band-separation-millimetres",
    type=float,
    default=DEFAULT_MINIMUM_BAND_SEPARATION_MILLIMETRES,
)
argument_parser.add_argument(
    "--band-peak-prominence-fraction",
    type=float,
    default=DEFAULT_BAND_PEAK_PROMINENCE_FRACTION,
)
argument_parser.add_argument(
    "--lane-strip-cross-section-fraction",
    type=float,
    default=DEFAULT_LANE_STRIP_CROSS_SECTION_FRACTION,
)
argument_parser.add_argument(
    "--use-fixed-strip-height",
    action="store_true",
    help="use one fixed strip height for every lane instead of the per-lane "
    "measured doublet width",
)
argument_parser.add_argument(
    "--lane-strip-fixed-half-height-millimetres",
    type=float,
    default=DEFAULT_LANE_STRIP_FIXED_HALF_HEIGHT_MILLIMETRES,
)
argument_parser.add_argument(
    "--band-cluster-tolerance-millimetres",
    type=float,
    default=None,
    help="tolerance for grouping per-lane band centres into one canonical band "
    "across lanes; defaults to the minimum band separation when omitted",
)
argument_parser.add_argument(
    "--baseline-flank-search-millimetres",
    type=float,
    default=DEFAULT_BASELINE_FLANK_SEARCH_MILLIMETRES,
)
argument_parser.add_argument(
    "--saturation-warn-fraction",
    type=float,
    default=DEFAULT_SATURATION_WARN_FRACTION,
    help="a cell is flagged saturated when this fraction of its window-by-strip "
    "pixels sit at or above the derived effective ceiling",
)
argument_parser.add_argument(
    "--baseline-cross-check-method",
    choices=["opening", "sternberg_ball"],
    default="opening",
    help="rolling-ball-family baseline used to cross-check valley-to-valley; the "
    "flat opening is the default, the Sternberg ball is a reference option",
)
argument_parser.add_argument(
    "--rolling-ball-width-millimetres",
    type=float,
    default=DEFAULT_ROLLING_BALL_WIDTH_MILLIMETRES,
    help="flat-opening structuring-element width for the cross-check baseline",
)
argument_parser.add_argument(
    "--sternberg-ball-radius-millimetres",
    type=float,
    default=DEFAULT_STERNBERG_BALL_RADIUS_MILLIMETRES,
    help="reference-only: Sternberg ball radius when that method is selected",
)
argument_parser.add_argument(
    "--sternberg-ball-vertical-scale",
    type=float,
    default=DEFAULT_STERNBERG_BALL_VERTICAL_SCALE,
    help="reference-only, experimental: Sternberg ball intensity-to-pixel scale; "
    "this is the unpinned degree of freedom that motivated the flat default",
)
argument_parser.add_argument(
    "--use-existing-band-definitions",
    action="store_true",
    help="read canonical band positions and windows from an existing "
    "band_definitions.csv (operator override) instead of deriving them",
)

argument_parser.add_argument(
    "--tilt-warning-threshold-degrees",
    type=float,
    default=DEFAULT_TILT_WARNING_THRESHOLD_DEGREES,
    help="tilt above this many degrees raises a warning; the gel is measured "
    "un-rotated regardless (straighten and re-scan, or rotate in Fiji before "
    "cropping)",
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
    die(
        "output",
        "expected output directory does not exist: "
        + str(output_directory_path)
        + ". Run stage 0 then stage 1 first.",
    )

validation_report_path = output_directory_path / VALIDATION_REPORT_FILENAME
if not validation_report_path.is_file():
    die(
        "stage1",
        "no stage 1 report at "
        + str(validation_report_path)
        + ". Run validate_gel_image.py first; the measurement stage reads its output.",
    )
try:
    validation_report = json.loads(validation_report_path.read_text())
except Exception as report_read_error:
    die("stage1", "stage 1 report will not parse: " + str(report_read_error))

# The measurement stage refuses to measure a file stage 1 flagged: the geometry it
# would rely on is exactly what a failing stage 1 says cannot be trusted.
if validation_report.get("overall_status") != "pass":
    die(
        "stage1",
        "stage 1 overall_status is "
        + repr(validation_report.get("overall_status"))
        + ", not 'pass'. Fix the "
        "sidecar and re-run stage 1 first.",
    )

preprocess_sidecar_block = validation_report.get("preprocess_sidecar", {})
geometry_pixels = preprocess_sidecar_block.get("geometry_pixels")
if geometry_pixels is None:
    die(
        "stage1",
        "the stage 1 report has no preprocess_sidecar.geometry_pixels, so "
        "there is no crop or axis to analyse. A sidecar is required.",
    )

gel_migration_axis = geometry_pixels["gel_migration_axis"]
crop_x_pixels = int(geometry_pixels["crop_x_pixels"])
crop_y_pixels = int(geometry_pixels["crop_y_pixels"])
crop_width_pixels = int(geometry_pixels["crop_width_pixels"])
crop_height_pixels = int(geometry_pixels["crop_height_pixels"])
# The tilt angle is the only landmark-derived value used downstream now that the
# resampler is gone; the raw landmark pixel coordinates were consumed solely by the
# rotation block, so they are no longer read here.
derived_tilt_angle_degrees = preprocess_sidecar_block.get("derived_tilt_angle_degrees")
expected_lane_count = int(
    preprocess_sidecar_block.get("values", {})["expected_lane_count"]
)
micrometres_per_pixel = validation_report.get("geometry", {}).get(
    "micrometres_per_pixel"
)
encoding_verified = validation_report.get("linearity_evidence", {}).get(
    "encoding_verified"
)

ANALYSIS_FINDINGS = []

# The read path, as float64 so nothing overflows. The image is read exactly once
# for the whole measurement stage; the raw uint16 array is kept for the effective-
# ceiling derivation in the band stage, and the ImageDescription is captured here
# so the band stage can read the instrument's dark offset and exposure mode without
# opening the file a second time.
raw_image_description_text = ""
try:
    with tifffile.TiffFile(str(input_tiff_absolute_path)) as tiff_handle:
        tiff_page = tiff_handle.pages[parsed_arguments.page_index]
        raw_pixel_array = tiff_page.asarray()
        image_description_tag = tiff_page.tags.get("ImageDescription")
        if image_description_tag is not None:
            raw_image_description_text = str(image_description_tag.value)
except Exception as pixel_read_error:
    die(
        "pixels",
        "page "
        + str(parsed_arguments.page_index)
        + " would not decode: "
        + str(pixel_read_error),
    )
pixel_array = raw_pixel_array.astype(numpy.float64)
pixel_array_height_pixels, pixel_array_width_pixels = pixel_array.shape

crop_is_inside_image = (
    crop_width_pixels > 0
    and crop_height_pixels > 0
    and 0 <= crop_x_pixels
    and 0 <= crop_y_pixels
    and crop_x_pixels + crop_width_pixels <= pixel_array_width_pixels
    and crop_y_pixels + crop_height_pixels <= pixel_array_height_pixels
)
if not crop_is_inside_image:
    die(
        "crop",
        "the stage 1 crop does not fit this image; the report and the file "
        "disagree, which stage 1 should have caught.",
    )

# Crop in the raw as-opened frame, where the sidecar coordinates were measured. The
# bottom-left-origin flip and lane numbering are applied in the band stage, where
# the band-numbering convention is defined; flipping now would only desynchronise
# this overlay from what the operator sees in Fiji. This one crop is what both the
# lane locator and the band integrator work from.
crop_pixels = pixel_array[
    crop_y_pixels : crop_y_pixels + crop_height_pixels,
    crop_x_pixels : crop_x_pixels + crop_width_pixels,
]

# Saturation on the crop. With no rotation there is no interpolation to move values
# off the ceiling, so this no longer has to be measured before any resampling step.
crop_maximum_value = float(crop_pixels.max())
at_container_maximum_in_crop = int(
    (crop_pixels >= CONTAINER_MAXIMUM_VALUE_16_BIT).sum()
)
at_crop_maximum_count = int((crop_pixels >= crop_maximum_value).sum())
ANALYSIS_FINDINGS.append(
    {
        "check_name": "in_crop_saturation",
        "status": "warning" if at_container_maximum_in_crop > 0 else "pass",
        "is_hard_stop": False,
        "detail": "crop maximum %d held by %d pixels; %d pixels at the container maximum %d"
        % (
            crop_maximum_value,
            at_crop_maximum_count,
            at_container_maximum_in_crop,
            CONTAINER_MAXIMUM_VALUE_16_BIT,
        ),
    }
)

# Tilt is measured, reported, and warned on, but never applied. Stage 1 derived the
# angle from the two landmark clicks; here it only decides whether to warn. The gel
# is measured on the un-rotated crop in every branch.
analysis_crop = crop_pixels
if derived_tilt_angle_degrees is None:
    ANALYSIS_FINDINGS.append(
        {
            "check_name": "tilt_measurement",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "stage 1 derived no tilt angle (landmark span below its floor), so "
            "tilt could not be checked; measured un-rotated.",
        }
    )
elif abs(derived_tilt_angle_degrees) <= parsed_arguments.tilt_warning_threshold_degrees:
    ANALYSIS_FINDINGS.append(
        {
            "check_name": "tilt_measurement",
            "status": "pass",
            "is_hard_stop": False,
            "detail": "tilt %.4f deg within the %.4f deg threshold; measured un-rotated."
            % (
                derived_tilt_angle_degrees,
                parsed_arguments.tilt_warning_threshold_degrees,
            ),
        }
    )
else:
    ANALYSIS_FINDINGS.append(
        {
            "check_name": "tilt_measurement",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "tilt %.4f deg exceeds the %.4f deg threshold; straighten and re-scan, "
            "or rotate in Fiji before cropping. Measured un-rotated."
            % (
                derived_tilt_angle_degrees,
                parsed_arguments.tilt_warning_threshold_degrees,
            ),
        }
    )

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
    lag_low = max(
        1, int(round(expected_pitch_pixels * (1.0 - LANE_PITCH_SEARCH_FRACTION)))
    )
    lag_high = max(
        lag_low + 1,
        int(round(expected_pitch_pixels * (1.0 + LANE_PITCH_SEARCH_FRACTION))),
    )
    lag_high = min(lag_high, axis_extent - 1)
    autocorrelation_lags = numpy.arange(lag_low, lag_high + 1)
    autocorrelation_curve = numpy.array(
        [
            float((detrended_profile[:-lag] * detrended_profile[lag:]).sum())
            / detrended_energy
            for lag in autocorrelation_lags
        ]
    )

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
if (
    stacking_autocorrelation_curve.size > 0
    and stacking_autocorrelation_curve.max() >= ORIENTATION_AUTOCORRELATION_FLOOR
):
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
    if micrometres_per_pixel is not None
    else None
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
    float(numpy.median(positive_column_medians))
    if positive_column_medians.size
    else 1.0
)
per_migration_column_saturated_fraction = numpy.mean(
    (analysis_crop if gel_migration_axis == "horizontal" else analysis_crop.T)
    >= CONTAINER_MAXIMUM_VALUE_16_BIT,
    axis=0,
)
usable_migration_column_mask = (
    per_migration_column_median <= EDGE_COLUMN_MEDIAN_MULTIPLE * typical_column_median
) & (per_migration_column_saturated_fraction < 0.05)
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
lane_prominence_threshold = LANE_PROMINENCE_FRACTION_OF_MAXIMUM * float(
    stacking_profile_smoothed_for_lane_finding.max()
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
    participating_column_count = int(
        (lane_migration_profile > participation_threshold).sum()
    )
    band_peaks, _ = scipy.signal.find_peaks(
        lane_migration_profile,
        distance=8,
        prominence=LANE_COLUMN_PARTICIPATION_THRESHOLD_FRACTION
        * lane_migration_profile_maximum,
    )
    is_real_lane = (
        participating_column_count
        >= LANE_MINIMUM_PARTICIPATING_COLUMN_FRACTION * usable_migration_column_count
        and len(band_peaks) >= 1
    )
    if not is_real_lane:
        rejected_candidate_rows.append(int(candidate_row))
        rejected_candidate_reasons.append(
            "participating_columns_%d_bands_%d"
            % (participating_column_count, len(band_peaks))
        )
        continue
    strip_rows = numpy.arange(strip_start, strip_end)
    strip_stacking_weight = stacking_profile_for_lane_finding[strip_start:strip_end]
    strip_weight_total = float(strip_stacking_weight.sum())
    refined_centre_row = (
        float((strip_rows * strip_stacking_weight).sum() / strip_weight_total)
        if strip_weight_total > 0.0
        else float(candidate_row)
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
    for candidate_phase_offset_pixels in numpy.linspace(
        0.0, measured_lane_pitch_pixels, 200
    ):
        nearest_well_index = numpy.round(
            (found_centres_array - candidate_phase_offset_pixels)
            / measured_lane_pitch_pixels
        )
        model_centre_positions = (
            candidate_phase_offset_pixels
            + nearest_well_index * measured_lane_pitch_pixels
        )
        total_squared_residual_pixels = float(
            numpy.sum((found_centres_array - model_centre_positions) ** 2)
        )
        if (
            best_total_squared_residual_pixels is None
            or total_squared_residual_pixels < best_total_squared_residual_pixels
        ):
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
stacking_gradient_energy = float(
    (stacking_axis_gradient * stacking_axis_gradient).sum()
)
migration_gradient_energy = float(
    (migration_axis_gradient * migration_axis_gradient).sum()
)
orientation_gradient_energy_ratio = (
    migration_gradient_energy / stacking_gradient_energy
    if stacking_gradient_energy > 0.0
    else 0.0
)
orientation_is_consistent = (
    orientation_gradient_energy_ratio >= ORIENTATION_GRADIENT_ENERGY_RATIO
)
ANALYSIS_FINDINGS.append(
    {
        "check_name": "orientation_matches_band_sharpness",
        "status": "pass" if orientation_is_consistent else "warning",
        "is_hard_stop": False,
        "detail": "gradient energy is %.2fx higher along the claimed migration axis than "
        "the claimed stacking axis (%s migration); bands are sharp along "
        "migration and smooth along stacking, so migration should carry more "
        "structure. "
        % (orientation_gradient_energy_ratio, gel_migration_axis)
        + (
            "consistent: the sharper axis is the one recorded as migration."
            if orientation_is_consistent
            else "the migration axis is not clearly the sharper one, so the axis may "
            "be mislabelled. Check the overlay and lane_profiles.png before "
            "trusting or fixing the sidecar."
        )
        + " Lane-pitch autocorrelation (secondary, unreliable on sparse gels) "
        "is %.3f stacking vs %.3f migration."
        % (stacking_autocorrelation, migration_autocorrelation),
    }
)

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
    lane_grid_status = (
        "pass" if populated_lane_count == expected_lane_count else "warning"
    )
    lane_grid_detail = (
        "found %d populated lane(s) of %d wells by peak-finding at a measured pitch "
        "of %.2f pixels"
        % (populated_lane_count, expected_lane_count, lane_pitch_pixels)
        + (
            ""
            if lane_pitch_millimetres is None
            else " (%.3f mm)" % lane_pitch_millimetres
        )
        + " (%s)" % lane_pitch_source
        + (
            ", assigned to well indices %s" % lane_well_indices
            if lane_well_indices is not None
            else ""
        )
        + (
            ", rejected %d candidate(s) as non-lane artifacts"
            % len(rejected_candidate_rows)
            if rejected_candidate_rows
            else ""
        )
        + ". Occupancy is measured from the image, not assumed; confirm against the "
        "overlay."
    )
ANALYSIS_FINDINGS.append(
    {
        "check_name": "lane_grid_fit",
        "status": lane_grid_status,
        "is_hard_stop": False,
        "detail": lane_grid_detail,
    }
)

# Overlay for the operator. The crop is drawn unflipped so it matches Fiji, with the
# fitted lane lines and the orientation verdict in the title.
overlay_stride = max(
    1,
    int(math.ceil(max(analysis_crop.shape) / float(OVERLAY_MAXIMUM_DIMENSION_PIXELS))),
)
overlay_array = analysis_crop[::overlay_stride, ::overlay_stride]
overlay_figure, overlay_axes = matplotlib.pyplot.subplots(figsize=(8.0, 8.0))
overlay_axes.imshow(
    overlay_array,
    cmap=DISPLAY_COLORMAP_NAME,
    vmin=numpy.percentile(overlay_array, 1.0),
    vmax=numpy.percentile(overlay_array, 99.5),
    interpolation="nearest",
    origin="upper",
)
if lane_centres_stacking_pixels is not None:
    for lane_centre in lane_centres_stacking_pixels:
        if gel_migration_axis == "horizontal":
            overlay_axes.axhline(
                lane_centre / overlay_stride, color="red", linewidth=0.8
            )
        else:
            overlay_axes.axvline(
                lane_centre / overlay_stride, color="red", linewidth=0.8
            )
overlay_title = (
    "Lane grid on the crop, unflipped: "
    + input_tiff_absolute_path.name
    + "\nmigration "
    + gel_migration_axis
    + (
        "; tilt %.4f deg" % derived_tilt_angle_degrees
        if derived_tilt_angle_degrees is not None
        else "; tilt unavailable"
    )
    + "; orientation "
    + ("OK" if orientation_is_consistent else "CHECK")
    + "\nmigration/stacking gradient-energy ratio %.2f (orientation cue)"
    % orientation_gradient_energy_ratio
    + (
        ""
        if lane_pitch_millimetres is None
        else "; pitch %.2f mm" % lane_pitch_millimetres
    )
)
overlay_axes.set_title(overlay_title, fontsize="medium")
overlay_axes.set_xlabel("crop column / " + str(overlay_stride))
overlay_axes.set_ylabel("crop row / " + str(overlay_stride))
overlay_figure.tight_layout()
overlay_output_path = output_directory_path / LANE_GRID_OVERLAY_FILENAME
try:
    overlay_figure.savefig(overlay_output_path, dpi=FIGURE_DOTS_PER_INCH)
except Exception as overlay_write_error:
    die(
        "overlay",
        "could not write " + str(overlay_output_path) + ": " + str(overlay_write_error),
    )
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
profiles_axes[0].plot(
    numpy.arange(stacking_detrended_profile.size),
    stacking_detrended_profile,
    color="black",
    linewidth=0.8,
)
if lane_centres_stacking_pixels is not None:
    for lane_centre in lane_centres_stacking_pixels:
        profiles_axes[0].axvline(lane_centre, color="red", linewidth=0.7)
profiles_axes[0].set_title(
    "collapsed along migration -> stacking-axis profile "
    "(should show the lanes); red = fitted comb centres",
    fontsize="small",
)
profiles_axes[0].set_xlabel("stacking-axis position (pixels)")
profiles_axes[1].plot(
    numpy.arange(migration_detrended_profile.size),
    migration_detrended_profile,
    color="black",
    linewidth=0.8,
)
profiles_axes[1].set_title(
    "collapsed along stacking -> migration-axis profile (should NOT be lane-periodic)",
    fontsize="small",
)
profiles_axes[1].set_xlabel("migration-axis position (pixels)")
for axis_name, curve_color in (("stacking", "red"), ("migration", "blue")):
    profiles_axes[2].plot(
        autocorrelation_lags_by_axis[axis_name],
        autocorrelation_curve_by_axis[axis_name],
        color=curve_color,
        linewidth=1.0,
        label=axis_name,
    )
profiles_axes[2].axhline(0.0, color="grey", linewidth=0.5)
profiles_axes[2].set_title(
    "lane-pitch autocorrelation vs lag; peak decides orientation (stacking should win)",
    fontsize="small",
)
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
            "%s,%d,%.6f"
            % (axis_name, position_index, float(axis_profile_values[position_index]))
        )
lane_profiles_csv_path = output_directory_path / LANE_PROFILES_CSV_FILENAME
lane_profiles_csv_path.write_text("\n".join(lane_profiles_csv_lines) + "\n")
emit_message("profiles", "wrote " + str(lane_profiles_csv_path))


# ---- Seam between the lane locator (above) and the band stage (below). The band
# stage consumes the lane centres, well indices, and pitch from memory rather than
# re-reading a stage 2 report from disk, and works from the same crop_pixels the
# locator used, so the two halves cannot disagree about the geometry or the crop.
if not lane_centres_stacking_pixels:
    die(
        "lanes",
        "no populated lanes were found, so there is nothing to measure "
        "bands in. Check the lane_grid_overlay before proceeding.",
    )
if lane_well_indices is None or len(lane_well_indices) != len(
    lane_centres_stacking_pixels
):
    # Not fatal: well indices are provisional and only used as a label here.
    lane_well_indices = list(range(len(lane_centres_stacking_pixels)))

# Convert the millimetre config to pixels once, using the report's pixel size.
millimetres_per_pixel = micrometres_per_pixel / 1000.0
saturation_warn_fraction = parsed_arguments.saturation_warn_fraction
rolling_ball_width_pixels = (
    parsed_arguments.rolling_ball_width_millimetres / millimetres_per_pixel
)
sternberg_ball_radius_pixels = (
    parsed_arguments.sternberg_ball_radius_millimetres / millimetres_per_pixel
)
region_window_width_pixels = (
    parsed_arguments.region_window_width_millimetres / millimetres_per_pixel
)
minimum_band_separation_pixels = max(
    1.0, parsed_arguments.minimum_band_separation_millimetres / millimetres_per_pixel
)
migration_profile_smoothing_pixels = max(
    1,
    int(round(DEFAULT_MIGRATION_PROFILE_SMOOTHING_MILLIMETRES / millimetres_per_pixel)),
)
lane_strip_fixed_half_height_pixels = int(
    round(
        parsed_arguments.lane_strip_fixed_half_height_millimetres
        / millimetres_per_pixel
    )
)

BAND_DETECTION_FINDINGS = []

# Parse the "Key = Value" lines of the instrument description for the two fields
# the saturation block records. Absent or unparseable fields degrade to None
# rather than fail: they are provenance, not control.
image_description_fields = {}
for description_line in raw_image_description_text.splitlines():
    if "=" in description_line:
        key_text, _, value_text = description_line.partition("=")
        image_description_fields[key_text.strip()] = value_text.strip()
exposure_mode_text = image_description_fields.get("Exposure mode")
dark_offset_text = image_description_fields.get("Dark offset")
dark_offset_value = None
if dark_offset_text is not None:
    try:
        dark_offset_value = int(round(float(dark_offset_text)))
    except ValueError:
        dark_offset_value = None

# Derive the effective ceiling from the whole page (more clipped pixels than the
# crop, so a knee, if any, is better evidenced). Default to the 16-bit container;
# lower it only if a real sub-ceiling pile-up is found. See the constants block for
# the definition of a knee. On s0002 no knee is found and this stays 65535.
page_value_counts = numpy.bincount(
    raw_pixel_array.astype(numpy.int64).ravel(),
    minlength=CONTAINER_MAXIMUM_VALUE_16_BIT + 1,
)
knee_search_floor_value = int(
    SUB_CEILING_KNEE_SEARCH_FLOOR_FRACTION * CONTAINER_MAXIMUM_VALUE_16_BIT
)
pixel_count_above_search_floor = int(page_value_counts[knee_search_floor_value:].sum())
sub_ceiling_knee_value = None
for candidate_value in range(
    CONTAINER_MAXIMUM_VALUE_16_BIT - 1, knee_search_floor_value, -1
):
    candidate_count = int(page_value_counts[candidate_value])
    if candidate_count == 0:
        continue
    local_baseline_start = max(
        knee_search_floor_value,
        candidate_value - SUB_CEILING_KNEE_LOCAL_BASELINE_SPAN_VALUES,
    )
    local_baseline_counts = page_value_counts[local_baseline_start:candidate_value]
    local_baseline_count = (
        float(numpy.median(local_baseline_counts))
        if local_baseline_counts.size
        else 0.0
    )
    pixels_strictly_above_candidate = int(
        page_value_counts[candidate_value + 1 :].sum()
    )
    above_fraction = pixels_strictly_above_candidate / float(
        max(1, pixel_count_above_search_floor)
    )
    candidate_is_spike = candidate_count >= SUB_CEILING_KNEE_SPIKE_FACTOR * max(
        1.0, local_baseline_count
    )
    candidate_is_empty_above = (
        above_fraction <= SUB_CEILING_KNEE_ABOVE_TOLERANCE_FRACTION
    )
    if candidate_is_spike and candidate_is_empty_above:
        sub_ceiling_knee_value = candidate_value
        break
if sub_ceiling_knee_value is not None:
    effective_ceiling_value = sub_ceiling_knee_value
    sub_ceiling_knee_detected = True
else:
    effective_ceiling_value = CONTAINER_MAXIMUM_VALUE_16_BIT
    sub_ceiling_knee_detected = False
page_pixels_at_effective_ceiling = int(
    page_value_counts[effective_ceiling_value:].sum()
)
crop_pixels_at_effective_ceiling = int((crop_pixels >= effective_ceiling_value).sum())

# Re-derive background and the edge-column mask exactly as Slice A did, rather than
# trust a passed array, so the two scripts cannot silently diverge on preprocessing.
plate_background_value = float(numpy.median(crop_pixels))
signal_above_plate = crop_pixels - plate_background_value
signal_above_plate[signal_above_plate < 0.0] = 0.0

# Orient so axis 0 is the stacking axis (lanes) and axis 1 is the migration axis
# (bands), whatever the scan orientation, so the per-lane logic is written once.
if gel_migration_axis == "horizontal":
    stacking_by_migration_signal = signal_above_plate
    oriented_raw_stacking_by_migration = crop_pixels
else:
    stacking_by_migration_signal = signal_above_plate.T
    oriented_raw_stacking_by_migration = crop_pixels.T
stacking_extent_pixels = stacking_by_migration_signal.shape[0]
migration_extent_pixels = stacking_by_migration_signal.shape[1]

per_migration_column_median = numpy.median(stacking_by_migration_signal, axis=0)
positive_column_medians = per_migration_column_median[per_migration_column_median > 0.0]
typical_column_median = (
    float(numpy.median(positive_column_medians))
    if positive_column_medians.size
    else 1.0
)
per_migration_column_saturated_fraction = numpy.mean(
    (crop_pixels if gel_migration_axis == "horizontal" else crop_pixels.T)
    >= effective_ceiling_value,
    axis=0,
)
usable_migration_column_mask = (
    per_migration_column_median <= EDGE_COLUMN_MEDIAN_MULTIPLE * typical_column_median
) & (per_migration_column_saturated_fraction < 0.05)
usable_stacking_by_migration_signal = numpy.where(
    usable_migration_column_mask[numpy.newaxis, :], stacking_by_migration_signal, 0.0
)

# Sort lanes by stacking position so neighbour midpoints are well defined, keeping
# each lane's well index attached.
lane_order = numpy.argsort(numpy.array(lane_centres_stacking_pixels))
ordered_lane_centres = [float(lane_centres_stacking_pixels[i]) for i in lane_order]
ordered_well_indices = [int(lane_well_indices[i]) for i in lane_order]

# Per-lane band detection.
per_lane_records = []
band_centers_csv_rows = []
for lane_position in range(len(ordered_lane_centres)):
    lane_centre_stacking = ordered_lane_centres[lane_position]
    well_index = ordered_well_indices[lane_position]

    # Bound this lane's strip at the midpoints to its neighbours so wide doublets
    # cannot bleed across lanes or be double-counted. Outermost lanes are capped at
    # one pitch from the centre so an isolated lane does not swallow empty gel.
    if lane_position > 0:
        lower_bound_stacking = 0.5 * (
            ordered_lane_centres[lane_position - 1] + lane_centre_stacking
        )
    else:
        lower_bound_stacking = lane_centre_stacking - lane_pitch_pixels
    if lane_position < len(ordered_lane_centres) - 1:
        upper_bound_stacking = 0.5 * (
            lane_centre_stacking + ordered_lane_centres[lane_position + 1]
        )
    else:
        upper_bound_stacking = lane_centre_stacking + lane_pitch_pixels
    search_top = max(0, int(math.floor(lower_bound_stacking)))
    search_bottom = min(stacking_extent_pixels, int(math.ceil(upper_bound_stacking)))

    # The lane's stacking cross-section within its bounded territory, collapsed over
    # the usable migration columns. Its extent above a fraction of its own maximum
    # is the measured doublet strip.
    territory_signal = usable_stacking_by_migration_signal[search_top:search_bottom, :]
    lane_cross_section = territory_signal.sum(axis=1)
    cross_section_maximum = (
        float(lane_cross_section.max()) if lane_cross_section.size else 0.0
    )
    if cross_section_maximum > 0.0:
        cross_section_threshold = (
            parsed_arguments.lane_strip_cross_section_fraction * cross_section_maximum
        )
        rows_above_threshold = numpy.where(
            lane_cross_section > cross_section_threshold
        )[0]
    else:
        rows_above_threshold = numpy.array([], dtype=int)
    if rows_above_threshold.size:
        measured_strip_top = search_top + int(rows_above_threshold.min())
        measured_strip_bottom = search_top + int(rows_above_threshold.max()) + 1
    else:
        measured_strip_top = search_top
        measured_strip_bottom = search_bottom
    measured_strip_height_pixels = measured_strip_bottom - measured_strip_top

    # The fixed-height alternative, always computed for audit, and used instead when
    # requested. Clamped to the bounded territory so it too cannot bleed.
    fixed_strip_top = max(
        search_top,
        int(round(lane_centre_stacking)) - lane_strip_fixed_half_height_pixels,
    )
    fixed_strip_bottom = min(
        search_bottom,
        int(round(lane_centre_stacking)) + lane_strip_fixed_half_height_pixels,
    )
    fixed_strip_height_pixels = fixed_strip_bottom - fixed_strip_top

    if parsed_arguments.use_fixed_strip_height:
        strip_top = fixed_strip_top
        strip_bottom = fixed_strip_bottom
    else:
        strip_top = measured_strip_top
        strip_bottom = measured_strip_bottom

    # The lane's migration profile: collapse the strip across stacking over usable
    # columns, lightly smoothed. Peaks are band centres.
    lane_strip_signal = usable_stacking_by_migration_signal[strip_top:strip_bottom, :]
    lane_migration_profile = lane_strip_signal.sum(axis=0)
    lane_migration_profile_smoothed = scipy.ndimage.uniform_filter1d(
        lane_migration_profile, size=migration_profile_smoothing_pixels, mode="nearest"
    )
    lane_profile_maximum = float(lane_migration_profile_smoothed.max())
    band_prominence_threshold = (
        parsed_arguments.band_peak_prominence_fraction * lane_profile_maximum
    )
    band_centre_columns, _ = scipy.signal.find_peaks(
        lane_migration_profile_smoothed,
        distance=minimum_band_separation_pixels,
        prominence=band_prominence_threshold,
    )

    # merged_center: two called centres closer than a window width, so the overlap
    # case is visible rather than silently split or merged.
    half_window_pixels = region_window_width_pixels / 2.0
    merged_center_flags = [False] * len(band_centre_columns)
    for band_position in range(len(band_centre_columns)):
        for other_position in range(len(band_centre_columns)):
            if (
                other_position != band_position
                and abs(
                    int(band_centre_columns[band_position])
                    - int(band_centre_columns[other_position])
                )
                < region_window_width_pixels
            ):
                merged_center_flags[band_position] = True

    lane_band_records = []
    for band_position in range(len(band_centre_columns)):
        band_centre_migration = int(band_centre_columns[band_position])
        window_start_migration = max(
            0, int(round(band_centre_migration - half_window_pixels))
        )
        window_end_migration = min(
            migration_extent_pixels,
            int(round(band_centre_migration + half_window_pixels)),
        )
        band_record = {
            "well_index": well_index,
            "lane_centre_stacking_pixels": round(lane_centre_stacking, 2),
            "band_index_in_lane": band_position,
            "band_centre_migration_pixels": band_centre_migration,
            "band_centre_migration_millimetres": round(
                band_centre_migration * millimetres_per_pixel, 4
            ),
            "window_start_migration_pixels": window_start_migration,
            "window_end_migration_pixels": window_end_migration,
            "lane_strip_top_stacking_pixels": strip_top,
            "lane_strip_bottom_stacking_pixels": strip_bottom,
            "lane_strip_measured_height_pixels": measured_strip_height_pixels,
            "lane_strip_fixed_height_pixels": fixed_strip_height_pixels,
            "peak_profile_value": round(
                float(lane_migration_profile_smoothed[band_centre_migration]), 2
            ),
            "merged_center_flag": bool(merged_center_flags[band_position]),
        }
        lane_band_records.append(band_record)
        band_centers_csv_rows.append(band_record)

    if len(band_centre_columns) == 0:
        BAND_DETECTION_FINDINGS.append(
            {
                "check_name": "lane_has_detectable_bands",
                "status": "warning",
                "is_hard_stop": False,
                "detail": "lane at well index %d (stacking %d) yielded no band peaks above "
                "the prominence floor; it may be empty or faint. Check the "
                "profile panel." % (well_index, int(round(lane_centre_stacking))),
            }
        )

    per_lane_records.append(
        {
            "well_index": well_index,
            "lane_centre_stacking_pixels": round(lane_centre_stacking, 2),
            "strip_top_stacking_pixels": strip_top,
            "strip_bottom_stacking_pixels": strip_bottom,
            "measured_strip_height_pixels": measured_strip_height_pixels,
            "fixed_strip_height_pixels": fixed_strip_height_pixels,
            "band_count": len(band_centre_columns),
            "band_centre_migration_pixels": [
                int(value) for value in band_centre_columns
            ],
            "merged_center_flags": merged_center_flags,
            "profile_smoothed": lane_migration_profile_smoothed,
        }
    )

total_band_count = len(band_centers_csv_rows)

# Display array, computed once and shared by every overlay. Separate from the
# measured signal per DESIGN 5.9: contrast scaling is display-only.
display_source = signal_above_plate
display_normalized = numpy.clip(
    display_source / numpy.percentile(display_source, 99.5), 0.0, 1.0
)

# ==========================================================================
# Slice B2: consensus band definition across lanes, then region integration.
# ==========================================================================
baseline_flank_search_pixels = max(
    1,
    int(
        round(
            parsed_arguments.baseline_flank_search_millimetres / millimetres_per_pixel
        )
    ),
)
if parsed_arguments.band_cluster_tolerance_millimetres is None:
    cluster_tolerance_pixels = minimum_band_separation_pixels
    cluster_tolerance_source = "tied_to_minimum_band_separation"
else:
    cluster_tolerance_pixels = (
        parsed_arguments.band_cluster_tolerance_millimetres / millimetres_per_pixel
    )
    cluster_tolerance_source = "operator_override"

band_definitions_csv_path = output_directory_path / BAND_DEFINITIONS_CSV_FILENAME

# Recompute each lane's migration profile once, keyed by well index, so both the
# integration and the baseline read the same profile the detection used.
lane_profile_by_well_index = {}
lane_strip_bounds_by_well_index = {}
lane_detected_centres_by_well_index = {}
for lane_record in per_lane_records:
    strip_top = lane_record["strip_top_stacking_pixels"]
    strip_bottom = lane_record["strip_bottom_stacking_pixels"]
    lane_profile_by_well_index[lane_record["well_index"]] = (
        usable_stacking_by_migration_signal[strip_top:strip_bottom, :].sum(axis=0)
    )
    lane_strip_bounds_by_well_index[lane_record["well_index"]] = (
        strip_top,
        strip_bottom,
    )
    lane_detected_centres_by_well_index[lane_record["well_index"]] = list(
        lane_record["band_centre_migration_pixels"]
    )

canonical_band_records = []
if (
    parsed_arguments.use_existing_band_definitions
    and band_definitions_csv_path.is_file()
):
    # Operator override: read canonical positions and windows verbatim. The operator
    # owns the boundaries; the tool applies exactly what the file says.
    definition_lines = band_definitions_csv_path.read_text().strip().splitlines()
    definition_header = definition_lines[0].split(",")
    for definition_line in definition_lines[1:]:
        definition_values = dict(zip(definition_header, definition_line.split(",")))
        canonical_band_records.append(
            {
                "canonical_position_migration_pixels": float(
                    definition_values["canonical_position_migration_pixels"]
                ),
                "window_start_migration_pixels": int(
                    definition_values["window_start_migration_pixels"]
                ),
                "window_end_migration_pixels": int(
                    definition_values["window_end_migration_pixels"]
                ),
                "member_wells": [
                    int(w)
                    for w in definition_values["member_wells"].split(";")
                    if w != ""
                ],
                "cross_lane_spread_pixels": int(
                    definition_values.get("cross_lane_spread_pixels", "0")
                ),
                "spread_warns": definition_values.get("spread_warns", "False")
                == "True",
            }
        )
    canonical_band_records.sort(
        key=lambda record: record["canonical_position_migration_pixels"]
    )
    band_definition_source = "operator_band_definitions_file"
else:
    # Derive canonical bands. Pool every detected centre with its well index, sorted
    # by migration position, then single-linkage cluster with same-lane exclusion.
    pooled_centre_records = []
    for lane_record in per_lane_records:
        for band_centre in lane_record["band_centre_migration_pixels"]:
            pooled_centre_records.append((int(band_centre), lane_record["well_index"]))
    pooled_centre_records.sort()

    canonical_clusters = []
    if pooled_centre_records:
        current_cluster = [pooled_centre_records[0]]
        for centre_column, well_index in pooled_centre_records[1:]:
            wells_in_current_cluster = {
                member_well for _, member_well in current_cluster
            }
            if (
                centre_column - current_cluster[-1][0] <= cluster_tolerance_pixels
                and well_index not in wells_in_current_cluster
            ):
                current_cluster.append((centre_column, well_index))
            else:
                canonical_clusters.append(current_cluster)
                current_cluster = [(centre_column, well_index)]
        canonical_clusters.append(current_cluster)

    for cluster in canonical_clusters:
        member_columns = [column for column, _ in cluster]
        member_wells = sorted({member_well for _, member_well in cluster})
        cross_lane_spread_pixels = int(max(member_columns) - min(member_columns))
        canonical_band_records.append(
            {
                "canonical_position_migration_pixels": float(
                    numpy.median(member_columns)
                ),
                "member_columns": member_columns,
                "member_wells": member_wells,
                "cross_lane_spread_pixels": cross_lane_spread_pixels,
                "spread_warns": cross_lane_spread_pixels
                >= CANONICAL_SPREAD_WARN_FRACTION * minimum_band_separation_pixels,
            }
        )
    canonical_band_records.sort(
        key=lambda record: record["canonical_position_migration_pixels"]
    )

    # Fixed-width windows centred on each canonical position, clipped at the midpoint
    # to each neighbour so windows never overlap and the valley between two crowded
    # bands is assigned to the nearer band, not counted twice.
    half_window_pixels = region_window_width_pixels / 2.0
    for band_position in range(len(canonical_band_records)):
        canonical_position = canonical_band_records[band_position][
            "canonical_position_migration_pixels"
        ]
        window_start = canonical_position - half_window_pixels
        window_end = canonical_position + half_window_pixels
        if band_position > 0:
            left_midpoint = 0.5 * (
                canonical_band_records[band_position - 1][
                    "canonical_position_migration_pixels"
                ]
                + canonical_position
            )
            window_start = max(window_start, left_midpoint)
        if band_position < len(canonical_band_records) - 1:
            right_midpoint = 0.5 * (
                canonical_position
                + canonical_band_records[band_position + 1][
                    "canonical_position_migration_pixels"
                ]
            )
            window_end = min(window_end, right_midpoint)
        canonical_band_records[band_position]["window_start_migration_pixels"] = max(
            0, int(round(window_start))
        )
        canonical_band_records[band_position]["window_end_migration_pixels"] = min(
            migration_extent_pixels, int(round(window_end))
        )
    band_definition_source = "derived_by_consensus"

# Findings for canonical bands whose cross-lane spread reaches the separation limit:
# the merge is uncertain, so warn with the value and a suggested action, and keep
# going, because the complete table is still wanted.
for canonical_position_index in range(len(canonical_band_records)):
    canonical_band = canonical_band_records[canonical_position_index]
    if canonical_band.get("spread_warns"):
        BAND_DETECTION_FINDINGS.append(
            {
                "check_name": "canonical_band_spread_within_separation",
                "status": "warning",
                "is_hard_stop": False,
                "detail": "canonical band near migration %d has a cross-lane spread of %d "
                "pixels, reaching the %.0f%% of minimum-separation flag; its "
                "members may not be the same species. Suggested action: check "
                "this band on the center map, re-crop if the lanes are "
                "misaligned, or set --band-cluster-tolerance-millimetres "
                "explicitly."
                % (
                    int(round(canonical_band["canonical_position_migration_pixels"])),
                    canonical_band["cross_lane_spread_pixels"],
                    100.0 * CANONICAL_SPREAD_WARN_FRACTION,
                ),
            }
        )

# Write the editable band-definition file (unless we just read it back), so the
# operator can override the boundaries and re-run with --use-existing-band-definitions.
if band_definition_source != "operator_band_definitions_file":
    definition_lines = [
        "canonical_band_index,canonical_position_migration_pixels,"
        "canonical_position_migration_millimetres,window_start_migration_pixels,"
        "window_end_migration_pixels,member_well_count,member_wells,"
        "cross_lane_spread_pixels,spread_warns"
    ]
    for canonical_position_index in range(len(canonical_band_records)):
        canonical_band = canonical_band_records[canonical_position_index]
        definition_lines.append(
            "%d,%.2f,%.4f,%d,%d,%d,%s,%d,%s"
            % (
                canonical_position_index,
                canonical_band["canonical_position_migration_pixels"],
                canonical_band["canonical_position_migration_pixels"]
                * millimetres_per_pixel,
                canonical_band["window_start_migration_pixels"],
                canonical_band["window_end_migration_pixels"],
                len(canonical_band["member_wells"]),
                ";".join(str(well) for well in canonical_band["member_wells"]),
                canonical_band["cross_lane_spread_pixels"],
                canonical_band["spread_warns"],
            )
        )
    band_definitions_csv_path.write_text("\n".join(definition_lines) + "\n")
    emit_message("bands", "wrote " + str(band_definitions_csv_path))

# Occupancy binding: each detected band in a lane belongs to its NEAREST canonical
# band only, so one detected band cannot be claimed as locally detected by two
# adjacent canonical positions.
locally_detected_well_canonical_pairs = set()
for lane_record in per_lane_records:
    for detected_centre in lane_record["band_centre_migration_pixels"]:
        nearest_canonical_index = None
        nearest_distance = None
        for canonical_position_index in range(len(canonical_band_records)):
            distance = abs(
                detected_centre
                - canonical_band_records[canonical_position_index][
                    "canonical_position_migration_pixels"
                ]
            )
            if nearest_distance is None or distance < nearest_distance:
                nearest_distance = distance
                nearest_canonical_index = canonical_position_index
        if nearest_canonical_index is not None:
            locally_detected_well_canonical_pairs.add(
                (lane_record["well_index"], nearest_canonical_index)
            )

# Region integration: apply the same canonical windows to every lane (complete
# rectangular table), integrate above a local valley-to-valley baseline, tag
# occupancy so a real band and background-at-a-shared-position are distinguishable.
# Cross-check pre-pass: the opened (or Sternberg-ball) baseline is one curve per
# lane, so compute it and every window's net-above-it up front, before the main
# loop, because the per-cell fragility flag must normalise by the lane signal scale
# (the largest opened net in the lane) and the no-signal gate by the strongest lane.
# The lane's opened-baseline floor and peak are stored too: a band whose opened
# baseline sits well above the lane floor is being bridged by a neighbour (a
# cluster shoulder), which is annotated per cell rather than passed off as a clean
# per-band baseline.
cross_check_baseline_by_well_index = {}
cross_check_net_by_well_and_canonical = {}
lane_signal_scale_by_well_index = {}
lane_climb_into_peak_by_well_index = {}
lane_baseline_floor_by_well_index = {}
lane_peak_value_by_well_index = {}
for lane_record in per_lane_records:
    well_index = lane_record["well_index"]
    lane_migration_profile = lane_profile_by_well_index[well_index]
    lane_cross_check_baseline = cross_check_baseline(
        lane_migration_profile,
        parsed_arguments.baseline_cross_check_method,
        rolling_ball_width_pixels,
        sternberg_ball_radius_pixels,
        parsed_arguments.sternberg_ball_vertical_scale,
    )
    cross_check_baseline_by_well_index[well_index] = lane_cross_check_baseline
    lane_net_by_canonical = {}
    for canonical_position_index in range(len(canonical_band_records)):
        window_start = canonical_band_records[canonical_position_index][
            "window_start_migration_pixels"
        ]
        window_end = canonical_band_records[canonical_position_index][
            "window_end_migration_pixels"
        ]
        if window_end > window_start:
            net_segment = (
                lane_migration_profile[window_start:window_end]
                - lane_cross_check_baseline[window_start:window_end]
            )
            net_segment[net_segment < 0.0] = 0.0
            lane_net_by_canonical[canonical_position_index] = float(net_segment.sum())
        else:
            lane_net_by_canonical[canonical_position_index] = 0.0
    cross_check_net_by_well_and_canonical[well_index] = lane_net_by_canonical
    lane_signal_scale_by_well_index[well_index] = (
        max(lane_net_by_canonical.values()) if lane_net_by_canonical else 0.0
    )
    lane_baseline_floor_by_well_index[well_index] = (
        float(lane_cross_check_baseline.min())
        if lane_cross_check_baseline.size
        else 0.0
    )
    # Width-adequacy climb: how far the baseline rides up at the lane's tallest
    # column, as a fraction of that peak. Small on an isolated band well-cleared by
    # the plate; large only if the width is too small (or the lane is a cluster).
    lane_peak_value = (
        float(lane_migration_profile.max()) if lane_migration_profile.size else 0.0
    )
    lane_peak_column = (
        int(numpy.argmax(lane_migration_profile)) if lane_migration_profile.size else 0
    )
    lane_peak_value_by_well_index[well_index] = lane_peak_value
    lane_climb_into_peak_by_well_index[well_index] = (
        float(lane_cross_check_baseline[lane_peak_column]) / lane_peak_value
        if lane_peak_value > 0.0
        else 0.0
    )

global_max_lane_signal_scale = (
    max(lane_signal_scale_by_well_index.values())
    if lane_signal_scale_by_well_index
    else 0.0
)
no_signal_scale_floor = (
    LANE_SIGNAL_SCALE_NO_SIGNAL_FRACTION * global_max_lane_signal_scale
)
# The best-case lane sets the width-adequacy verdict: if even it rides up, the
# width is too small for every lane, not just the clusters.
minimum_lane_climb_into_peak = (
    min(lane_climb_into_peak_by_well_index.values())
    if lane_climb_into_peak_by_well_index
    else 0.0
)

band_measurement_rows = []
valley_to_valley_display_segments_by_well = {
    lane_record["well_index"]: [] for lane_record in per_lane_records
}
locally_detected_cell_count = 0
consensus_only_cell_count = 0
saturated_cell_count = 0
baseline_fragile_cell_count = 0
multi_maximum_window_count = 0
for lane_record in per_lane_records:
    well_index = lane_record["well_index"]
    lane_migration_profile = lane_profile_by_well_index[well_index]
    strip_top, strip_bottom = lane_strip_bounds_by_well_index[well_index]
    for canonical_position_index in range(len(canonical_band_records)):
        canonical_band = canonical_band_records[canonical_position_index]
        window_start = canonical_band["window_start_migration_pixels"]
        window_end = canonical_band["window_end_migration_pixels"]
        canonical_position = canonical_band["canonical_position_migration_pixels"]

        # Local valley-to-valley baseline: flanking minima just outside the window.
        left_search_start = max(0, window_start - baseline_flank_search_pixels)
        left_flank = lane_migration_profile[left_search_start : window_start + 1]
        if left_flank.size:
            left_reference_offset = int(numpy.argmin(left_flank))
            left_reference_column = left_search_start + left_reference_offset
            left_reference_value = float(left_flank[left_reference_offset])
        else:
            left_reference_column = window_start
            left_reference_value = float(lane_migration_profile[window_start])
        right_search_end = min(
            migration_extent_pixels, window_end + baseline_flank_search_pixels
        )
        right_flank = lane_migration_profile[
            max(window_start, window_end - 1) : right_search_end
        ]
        if right_flank.size:
            right_reference_offset = int(numpy.argmin(right_flank))
            right_reference_column = (
                max(window_start, window_end - 1) + right_reference_offset
            )
            right_reference_value = float(right_flank[right_reference_offset])
        else:
            right_reference_column = max(window_start, window_end - 1)
            right_reference_value = float(
                lane_migration_profile[right_reference_column]
            )

        window_columns = numpy.arange(window_start, window_end)
        # Peak mode and multiplicity, initialised so the empty-window branch is safe.
        # Apex is a height, so it is independent of the window width that makes the
        # region net a partial-peak integral; multiplicity counts maxima inside the
        # window (lane-scaled prominence) to mark where the single number conflates
        # sub-bands. Both are profile properties, shared by the cell's two rows.
        peak_apex_raw_value = 0.0
        peak_apex_migration_pixels = int(window_start)
        peak_apex_above_valley_to_valley = 0.0
        sub_peak_count = 0
        if window_columns.size == 0:
            integrated_intensity = 0.0
            raw_window_sum = 0.0
            baseline_mean_value = 0.0
            net_above_baseline = 0.0
            measurement_status = "empty_window"
        else:
            if right_reference_column != left_reference_column:
                baseline_over_window = left_reference_value + (
                    (right_reference_value - left_reference_value)
                    * (window_columns - left_reference_column)
                    / float(right_reference_column - left_reference_column)
                )
            else:
                baseline_over_window = numpy.full(
                    window_columns.shape, left_reference_value
                )
            # Keep the drawn valley-to-valley baseline for the comparison panel: the
            # sloped line actually used, as a display array separate from the measured
            # net, so the panel shows both baselines rather than the opening alone.
            valley_to_valley_display_segments_by_well[well_index].append(
                (
                    int(window_start),
                    int(window_end),
                    float(left_reference_value),
                    float(right_reference_value),
                    int(left_reference_column),
                    int(right_reference_column),
                )
            )
            profile_over_window = lane_migration_profile[window_start:window_end]
            raw_window_sum = float(profile_over_window.sum())
            net_above_baseline = float(
                (profile_over_window - baseline_over_window).sum()
            )
            baseline_mean_value = float(baseline_over_window.mean())
            # Clip the reported intensity at zero: a window on the tail of a strong
            # neighbouring band can put the sloped valley baseline above the profile,
            # giving a spurious large negative. The honest value there is no band. The
            # unclipped net is kept for audit, and a flag marks where clipping bit,
            # which is exactly where B3's rolling-ball baseline is the cross-check.
            if net_above_baseline < 0.0:
                integrated_intensity = 0.0
                if window_start == 0 or window_end == migration_extent_pixels:
                    measurement_status = "window_truncated_at_crop_edge"
                else:
                    measurement_status = "baseline_above_signal_clipped_to_zero"
            else:
                integrated_intensity = net_above_baseline
                if window_start == 0 or window_end == migration_extent_pixels:
                    measurement_status = "window_truncated_at_crop_edge"
                else:
                    measurement_status = "ok"

            # Peak apex within the window, raw and above the valley-to-valley
            # baseline, plus the column where it falls. The apex is window-width
            # independent, so it should track band strength more stably than the
            # partial-peak region net; whether it actually is more consistent here
            # is a claim to read off the numbers, not assume.
            apex_offset = int(numpy.argmax(profile_over_window))
            peak_apex_raw_value = float(profile_over_window[apex_offset])
            peak_apex_migration_pixels = int(window_start + apex_offset)
            peak_apex_above_valley_to_valley = float(
                (profile_over_window - baseline_over_window).max()
            )
            # Multiplicity: maxima inside the window above a lane-scaled prominence,
            # reusing the same prominence fraction as B1 detection. More than one
            # means the window holds sub-structure and its single integral is not a
            # per-band number; that is where deconvolution would earn its cost.
            window_sub_peak_columns, _ = scipy.signal.find_peaks(
                profile_over_window,
                prominence=parsed_arguments.band_peak_prominence_fraction
                * lane_peak_value_by_well_index[well_index],
            )
            sub_peak_count = int(window_sub_peak_columns.size)

        # Per-cell saturation, counted on the raw-DN 2-D patch (strip by window),
        # not on the collapsed profile, because a clipped pixel is a 2-D pixel that
        # summing to a profile would hide. Saturation is a raw-value property, so it
        # is judged against the raw crop and the derived effective ceiling, entirely
        # independent of the background subtraction and baseline. A saturated cell's
        # integrated_intensity is a floor: the clipped pixels lost their true height.
        saturation_patch = oriented_raw_stacking_by_migration[
            strip_top:strip_bottom, window_start:window_end
        ]
        saturation_patch_pixel_count = int(saturation_patch.size)
        saturated_pixel_count = int((saturation_patch >= effective_ceiling_value).sum())
        saturated_pixel_fraction = (
            saturated_pixel_count / float(saturation_patch_pixel_count)
            if saturation_patch_pixel_count
            else 0.0
        )
        if saturated_pixel_fraction > saturation_warn_fraction:
            saturation_status = "saturated_intensity_is_floor"
            saturated_cell_count += 1
        else:
            saturation_status = "clean"

        occupancy = (
            "locally_detected"
            if (well_index, canonical_position_index)
            in locally_detected_well_canonical_pairs
            else "consensus_only"
        )
        if occupancy == "locally_detected":
            locally_detected_cell_count += 1
        else:
            consensus_only_cell_count += 1

        # Cross-baseline agreement. The opened net is >= 0 by construction, so the
        # comparison uses the UNCLIPPED valley-to-valley net; that keeps a sign flip
        # (vtv drove the clipped intensity to zero while the opening still sees a
        # band) visible instead of hidden by the clip. Disagreement is normalised by
        # the lane signal scale so near-empty cells cannot read as maximally fragile,
        # and it is only judged where the lane carries signal. A canonical band that
        # shares the opening width with a neighbour is bridged by the plate, so its
        # opened net is a per-cluster quantity and is labelled as such rather than
        # offered as a clean per-band disagreement.
        lane_signal_scale = lane_signal_scale_by_well_index[well_index]
        cross_check_net = cross_check_net_by_well_and_canonical[well_index][
            canonical_position_index
        ]
        cross_check_baseline_mean_value = (
            float(
                cross_check_baseline_by_well_index[well_index][
                    window_start:window_end
                ].mean()
            )
            if window_end > window_start
            else 0.0
        )
        if lane_signal_scale > 0.0:
            disagreement_fraction = (
                abs(net_above_baseline - cross_check_net) / lane_signal_scale
            )
        else:
            disagreement_fraction = 0.0
        # A sign flip only counts when the opening sees a NON-trivial band where vtv
        # clipped to nothing; otherwise both baselines are agreeing on ~zero and the
        # bare inequality (0.0 vs a hair above 0) is noise, not disagreement.
        sign_flip_vtv_clipped_opening_positive = bool(
            integrated_intensity == 0.0 and cross_check_net > no_signal_scale_floor
        )
        # Elevation annotation: is the opened baseline under this band well above the
        # lane's own baseline floor? If so a neighbour is holding the plate up and the
        # opened net is a cluster-shoulder quantity. Kept separate from the verdict so
        # it explains, rather than overrides, agree/fragile.
        lane_baseline_floor = lane_baseline_floor_by_well_index[well_index]
        lane_peak_value = lane_peak_value_by_well_index[well_index]
        lane_peak_above_floor = lane_peak_value - lane_baseline_floor
        opening_baseline_elevated = bool(
            lane_peak_above_floor > 0.0
            and (cross_check_baseline_mean_value - lane_baseline_floor)
            > OPENING_BASELINE_ELEVATED_FRACTION * lane_peak_above_floor
        )
        if lane_signal_scale < no_signal_scale_floor:
            baseline_agreement_status = "no_signal"
        elif disagreement_fraction >= BASELINE_FRAGILE_DISAGREEMENT_FRACTION:
            baseline_agreement_status = "fragile"
            baseline_fragile_cell_count += 1
        else:
            baseline_agreement_status = "agree"

        # Reported-intensity bracket: the two methods give two numbers for the same
        # cell (valley-to-valley clipped = the per-band lower estimate, opening =
        # the cluster-inclusive upper estimate). On a fragile cluster cell the honest
        # statement is the interval, not either endpoint. Kept as min/max so it is a
        # valid interval regardless of which method is higher.
        reported_intensity_lower = min(integrated_intensity, cross_check_net)
        reported_intensity_upper = max(integrated_intensity, cross_check_net)
        reported_intensity_bracket_width = (
            reported_intensity_upper - reported_intensity_lower
        )
        # Apex above the opening baseline for the opening row (>= 0 by construction).
        if window_end > window_start:
            peak_apex_above_opening = float(
                (
                    lane_migration_profile[window_start:window_end]
                    - cross_check_baseline_by_well_index[well_index][
                        window_start:window_end
                    ]
                ).max()
            )
        else:
            peak_apex_above_opening = 0.0
        if sub_peak_count > 1:
            multi_maximum_window_count += 1

        # The reported-value contract, decided from the s0002 comparison figures and
        # meant to adapt per cell rather than impose one global metric. It stays a
        # single consistent quantity: an AREA (summed intensity), never mixed with a
        # height, so the column is comparable across cells. Saturation caps everything
        # (the intensity is a floor), so it wins first. A fragile cell's per-band
        # valley-to-valley area is baseline-sensitive and can clip to zero, so the
        # cluster-inclusive opening area (the upper end of the bracket) is reported
        # instead, with the bracket width carrying the uncertainty and the apex
        # columns available separately for anyone who wants a height. Everywhere else
        # the standard valley-to-valley region net is the number. The basis is
        # recorded so the published value is always traceable to how it was chosen.
        if saturation_status != "clean":
            reported_value = integrated_intensity
            reported_value_basis = "saturated_floor"
        elif baseline_agreement_status == "fragile":
            reported_value = cross_check_net
            reported_value_basis = "opening_net_cluster_inclusive"
        else:
            reported_value = integrated_intensity
            reported_value_basis = "region_net"

        # Fields shared by both baseline rows of this cell. The agreement fields are
        # properties of the cell comparison, so they ride on both rows identically.
        shared_cell_fields = {
            "well_index": well_index,
            "lane_centre_stacking_pixels": round(
                lane_record["lane_centre_stacking_pixels"], 2
            ),
            "canonical_band_index": canonical_position_index,
            "canonical_position_migration_pixels": int(round(canonical_position)),
            "canonical_position_migration_millimetres": round(
                canonical_position * millimetres_per_pixel, 4
            ),
            "window_start_migration_pixels": window_start,
            "window_end_migration_pixels": window_end,
            "lane_strip_top_stacking_pixels": strip_top,
            "lane_strip_bottom_stacking_pixels": strip_bottom,
            "integration_method": "region",
            "occupancy": occupancy,
            "canonical_member_well_count": len(canonical_band["member_wells"]),
            "cross_lane_spread_pixels": canonical_band["cross_lane_spread_pixels"],
            "measurement_status": measurement_status,
            "saturated_pixel_count": saturated_pixel_count,
            "saturated_pixel_fraction": round(saturated_pixel_fraction, 5),
            "saturation_status": saturation_status,
            "baseline_disagreement_fraction": round(disagreement_fraction, 5),
            "baseline_agreement_status": baseline_agreement_status,
            "sign_flip_vtv_clipped_opening_positive": sign_flip_vtv_clipped_opening_positive,
            "opening_baseline_elevated": opening_baseline_elevated,
            "peak_apex_raw_value": round(peak_apex_raw_value, 2),
            "peak_apex_migration_pixels": peak_apex_migration_pixels,
            "sub_peak_count": sub_peak_count,
            "reported_intensity_bracket_width": round(
                reported_intensity_bracket_width, 2
            ),
            "reported_value": round(reported_value, 2),
            "reported_value_basis": reported_value_basis,
        }
        valley_to_valley_row = dict(shared_cell_fields)
        valley_to_valley_row.update(
            {
                "baseline_method": "valley_to_valley",
                "integrated_intensity": round(integrated_intensity, 2),
                "net_above_baseline_unclipped": round(net_above_baseline, 2),
                "raw_window_sum": round(raw_window_sum, 2),
                "baseline_mean_value": round(baseline_mean_value, 2),
                "peak_apex_above_baseline": round(peak_apex_above_valley_to_valley, 2),
            }
        )
        opening_row = dict(shared_cell_fields)
        opening_row.update(
            {
                "baseline_method": parsed_arguments.baseline_cross_check_method,
                "integrated_intensity": round(cross_check_net, 2),
                "net_above_baseline_unclipped": round(cross_check_net, 2),
                "raw_window_sum": round(raw_window_sum, 2),
                "baseline_mean_value": round(cross_check_baseline_mean_value, 2),
                "peak_apex_above_baseline": round(peak_apex_above_opening, 2),
            }
        )
        band_measurement_rows.append(valley_to_valley_row)
        band_measurement_rows.append(opening_row)

canonical_band_count = len(canonical_band_records)

# Surface saturation as findings. A sub-ceiling knee means the sensor clipped below
# the container, which changes the ceiling every downstream count is judged against,
# so it is worth stating even when no cell ends up flagged. Saturated cells mean the
# affected integrated intensities are floors, not measures. Both are warnings with a
# suggested action, not hard stops; a clean image adds nothing here.
if sub_ceiling_knee_detected:
    BAND_DETECTION_FINDINGS.append(
        {
            "check_name": "effective_ceiling_below_container",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "a sub-ceiling saturation knee was detected at DN %d (container is "
            "%d); intensities are counted as saturated against the knee, not the "
            "container. Confirm against the saturation_derivation panel."
            % (effective_ceiling_value, CONTAINER_MAXIMUM_VALUE_16_BIT),
        }
    )
if saturated_cell_count > 0:
    BAND_DETECTION_FINDINGS.append(
        {
            "check_name": "no_cell_saturated",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "%d cell(s) have more than %.3f%% of their pixels at or above the "
            "effective ceiling %d; their integrated_intensity is a floor, not a "
            "measure. See saturation_status in band_measurements.csv and re-scan "
            "at a shorter exposure to quantify those bands."
            % (
                saturated_cell_count,
                saturation_warn_fraction * 100.0,
                effective_ceiling_value,
            ),
        }
    )

# Width-adequacy: if even the best-case lane's opened baseline still rides up into
# its tallest band, the cross-check width is too small everywhere and the opened
# nets are systematically low. Fired on the min-over-lanes climb so clusters, which
# legitimately ride up, do not raise a false alarm.
if minimum_lane_climb_into_peak > WIDTH_ADEQUACY_CLIMB_WARN_FRACTION:
    BAND_DETECTION_FINDINGS.append(
        {
            "check_name": "cross_check_width_clears_bands",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "the %s baseline rides %.1f%% into the tallest band of even the "
            "clearest lane; the width (%.2f mm) is likely too small and its "
            "nets are low. Increase --rolling-ball-width-millimetres and "
            "confirm against baseline_comparison.png."
            % (
                parsed_arguments.baseline_cross_check_method,
                minimum_lane_climb_into_peak * 100.0,
                parsed_arguments.rolling_ball_width_millimetres,
            ),
        }
    )

# Baseline-fragile cells: the two baselines disagree by more than the threshold on a
# signal-bearing, non-clustered cell, so which number to trust is not obvious there.
if baseline_fragile_cell_count > 0:
    BAND_DETECTION_FINDINGS.append(
        {
            "check_name": "baselines_agree_on_signal_cells",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "%d signal-bearing cell(s) disagree between valley-to-valley and "
            "the %s baseline by more than %.0f%% of the lane signal scale; treat "
            "their intensity as baseline-sensitive. See baseline_agreement_status "
            "in band_measurements.csv and baseline_comparison.png."
            % (
                baseline_fragile_cell_count,
                parsed_arguments.baseline_cross_check_method,
                BASELINE_FRAGILE_DISAGREEMENT_FRACTION * 100.0,
            ),
        }
    )

# Multiplicity: a window holding more than one maximum conflates sub-bands, so its
# single integrated number is not per-band there. This is where peak deconvolution
# would earn its cost; until then the cell carries sub_peak_count and the bracket.
if multi_maximum_window_count > 0:
    BAND_DETECTION_FINDINGS.append(
        {
            "check_name": "one_maximum_per_window",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "%d window(s) contain more than one maximum; their single "
            "integrated_intensity conflates sub-bands and is not a per-band "
            "number. See sub_peak_count in band_measurements.csv; peak "
            "deconvolution would be the resolving step." % multi_maximum_window_count,
        }
    )

# band_measurements.csv, long format per DESIGN 5.10: one row per lane x canonical
# band x baseline method (one baseline method in B2; B3 adds rolling-ball).
measurement_column_names = [
    "well_index",
    "lane_centre_stacking_pixels",
    "canonical_band_index",
    "canonical_position_migration_pixels",
    "canonical_position_migration_millimetres",
    "window_start_migration_pixels",
    "window_end_migration_pixels",
    "lane_strip_top_stacking_pixels",
    "lane_strip_bottom_stacking_pixels",
    "integration_method",
    "baseline_method",
    "integrated_intensity",
    "net_above_baseline_unclipped",
    "raw_window_sum",
    "baseline_mean_value",
    "occupancy",
    "canonical_member_well_count",
    "cross_lane_spread_pixels",
    "measurement_status",
    "saturated_pixel_count",
    "saturated_pixel_fraction",
    "saturation_status",
    "baseline_disagreement_fraction",
    "baseline_agreement_status",
    "sign_flip_vtv_clipped_opening_positive",
    "opening_baseline_elevated",
    "peak_apex_raw_value",
    "peak_apex_above_baseline",
    "peak_apex_migration_pixels",
    "sub_peak_count",
    "reported_intensity_bracket_width",
    "reported_value",
    "reported_value_basis",
]
measurement_csv_lines = [",".join(measurement_column_names)]
for measurement_row in band_measurement_rows:
    measurement_csv_lines.append(
        ",".join(str(measurement_row[name]) for name in measurement_column_names)
    )
band_measurements_csv_path = output_directory_path / BAND_MEASUREMENTS_CSV_FILENAME
band_measurements_csv_path.write_text("\n".join(measurement_csv_lines) + "\n")
emit_message("csv", "wrote " + str(band_measurements_csv_path))

# Consensus center map: per-lane detected centres by well, with canonical positions
# and their windows drawn, so the operator sees how the clustering grouped centres
# and where the flagged (wide-spread) bands are.
center_map_figure, center_map_axes = matplotlib.pyplot.subplots(figsize=(11.0, 4.5))
for lane_record in per_lane_records:
    for band_centre in lane_record["band_centre_migration_pixels"]:
        center_map_axes.plot(
            band_centre, lane_record["well_index"], "o", color="black", markersize=4
        )
for canonical_band in canonical_band_records:
    canonical_position = canonical_band["canonical_position_migration_pixels"]
    line_color = "orange" if canonical_band.get("spread_warns") else "red"
    center_map_axes.axvline(canonical_position, color=line_color, linewidth=0.9)
    center_map_axes.axvspan(
        canonical_band["window_start_migration_pixels"],
        canonical_band["window_end_migration_pixels"],
        color=line_color,
        alpha=0.08,
    )
center_map_axes.set_xlabel("migration position (pixels, crop-local)")
center_map_axes.set_ylabel("well index")
center_map_axes.set_title(
    "Consensus center map: black = per-lane detected centre; red line/band = "
    "canonical band and window; orange = spread reaches separation limit (%d canonical bands)"
    % canonical_band_count,
    fontsize="small",
)
center_map_figure.tight_layout()
consensus_center_map_path = output_directory_path / CONSENSUS_CENTER_MAP_FILENAME
center_map_figure.savefig(consensus_center_map_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(center_map_figure)
emit_message("map", "wrote " + str(consensus_center_map_path))

# Occupancy-coloured overlay: every canonical window on every lane, green where the
# lane locally detected a band, red where it is measured only because other lanes
# had a band there (the phantom cells, labelled not hidden).
occupancy_figure, occupancy_axes = matplotlib.pyplot.subplots(figsize=(9.0, 9.0))
occupancy_axes.imshow(
    display_normalized,
    cmap=DISPLAY_COLORMAP_NAME,
    interpolation="nearest",
    origin="upper",
)
for measurement_row in band_measurement_rows:
    if measurement_row["measurement_status"] == "empty_window":
        continue
    window_start = measurement_row["window_start_migration_pixels"]
    window_end = measurement_row["window_end_migration_pixels"]
    strip_top = measurement_row["lane_strip_top_stacking_pixels"]
    strip_bottom = measurement_row["lane_strip_bottom_stacking_pixels"]
    cell_color = (
        "green" if measurement_row["occupancy"] == "locally_detected" else "red"
    )
    if gel_migration_axis == "horizontal":
        cell_box = matplotlib.patches.Rectangle(
            (window_start, strip_top),
            window_end - window_start,
            strip_bottom - strip_top,
            fill=False,
            edgecolor=cell_color,
            linewidth=0.6,
        )
    else:
        cell_box = matplotlib.patches.Rectangle(
            (strip_top, window_start),
            strip_bottom - strip_top,
            window_end - window_start,
            fill=False,
            edgecolor=cell_color,
            linewidth=0.6,
        )
    occupancy_axes.add_patch(cell_box)
occupancy_axes.set_title(
    "Consensus windows on every lane (B2): "
    + input_tiff_absolute_path.name
    + "\ngreen = locally detected, red = consensus-only (measured, phantom-tagged); "
    "%d cells" % len(band_measurement_rows),
    fontsize="medium",
)
occupancy_axes.set_xlabel("crop column / 1")
occupancy_axes.set_ylabel("crop row / 1")
occupancy_figure.tight_layout()
consensus_windows_overlay_path = (
    output_directory_path / CONSENSUS_WINDOWS_OVERLAY_FILENAME
)
occupancy_figure.savefig(consensus_windows_overlay_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(occupancy_figure)
emit_message("overlay", "wrote " + str(consensus_windows_overlay_path))

# Overlay for the operator: display-scaled crop on gray_r so bands read dark and it
# matches Fiji. Each lane's strip is drawn as two boundary lines; each band window
# is a box spanning the strip, with the centre marked.
overlay_figure, overlay_axes = matplotlib.pyplot.subplots(figsize=(9.0, 9.0))
overlay_axes.imshow(
    display_normalized,
    cmap=DISPLAY_COLORMAP_NAME,
    interpolation="nearest",
    origin="upper",
)
for lane_record in per_lane_records:
    strip_top = lane_record["strip_top_stacking_pixels"]
    strip_bottom = lane_record["strip_bottom_stacking_pixels"]
    for band_record_position in range(len(lane_record["band_centre_migration_pixels"])):
        band_centre_migration = lane_record["band_centre_migration_pixels"][
            band_record_position
        ]
        merged = lane_record["merged_center_flags"][band_record_position]
        window_start = max(
            0, int(round(band_centre_migration - region_window_width_pixels / 2.0))
        )
        window_width = int(round(region_window_width_pixels))
        box_edge_color = "orange" if merged else "red"
        if gel_migration_axis == "horizontal":
            box = matplotlib.patches.Rectangle(
                (window_start, strip_top),
                window_width,
                strip_bottom - strip_top,
                fill=False,
                edgecolor=box_edge_color,
                linewidth=0.8,
            )
            overlay_axes.add_patch(box)
            overlay_axes.plot(
                [band_centre_migration, band_centre_migration],
                [strip_top, strip_bottom],
                color="blue",
                linewidth=0.5,
            )
        else:
            box = matplotlib.patches.Rectangle(
                (strip_top, window_start),
                strip_bottom - strip_top,
                window_width,
                fill=False,
                edgecolor=box_edge_color,
                linewidth=0.8,
            )
            overlay_axes.add_patch(box)
            overlay_axes.plot(
                [strip_top, strip_bottom],
                [band_centre_migration, band_centre_migration],
                color="blue",
                linewidth=0.5,
            )
overlay_axes.set_title(
    "Band windows per lane (B1, positions only): "
    + input_tiff_absolute_path.name
    + "\nmigration "
    + gel_migration_axis
    + "; %d lanes, %d bands; red=window, blue=centre, orange=merged"
    % (len(per_lane_records), total_band_count),
    fontsize="medium",
)
overlay_axes.set_xlabel("crop column / 1")
overlay_axes.set_ylabel("crop row / 1")
overlay_figure.tight_layout()
band_windows_overlay_path = output_directory_path / BAND_WINDOWS_OVERLAY_FILENAME
try:
    overlay_figure.savefig(band_windows_overlay_path, dpi=FIGURE_DOTS_PER_INCH)
except Exception as overlay_write_error:
    die(
        "overlay",
        "could not write "
        + str(band_windows_overlay_path)
        + ": "
        + str(overlay_write_error),
    )
matplotlib.pyplot.close(overlay_figure)
emit_message("overlay", "wrote " + str(band_windows_overlay_path))

# Per-lane profile panels, so the operator can see directly whether the peaks match
# the bands and whether the prominence gate is sane.
panel_count = len(per_lane_records)
profiles_figure, profiles_axes = matplotlib.pyplot.subplots(
    panel_count, 1, figsize=(9.0, max(2.0, 1.7 * panel_count)), squeeze=False
)
for panel_index in range(panel_count):
    lane_record = per_lane_records[panel_index]
    profile_values = lane_record["profile_smoothed"]
    panel_axes = profiles_axes[panel_index][0]
    panel_axes.plot(
        numpy.arange(profile_values.size), profile_values, color="black", linewidth=0.8
    )
    for band_record_position in range(len(lane_record["band_centre_migration_pixels"])):
        band_centre_migration = lane_record["band_centre_migration_pixels"][
            band_record_position
        ]
        panel_axes.axvline(band_centre_migration, color="blue", linewidth=0.6)
        panel_axes.axvspan(
            max(0, band_centre_migration - region_window_width_pixels / 2.0),
            band_centre_migration + region_window_width_pixels / 2.0,
            color="red",
            alpha=0.12,
        )
    panel_axes.set_title(
        "well index %d (stacking %d): %d bands"
        % (
            lane_record["well_index"],
            int(round(lane_record["lane_centre_stacking_pixels"])),
            lane_record["band_count"],
        ),
        fontsize="small",
    )
    panel_axes.set_xlabel("migration position (pixels, crop-local)")
profiles_figure.tight_layout()
lane_profiles_path = output_directory_path / LANE_MIGRATION_PROFILES_FILENAME
profiles_figure.savefig(lane_profiles_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(profiles_figure)
emit_message("profiles", "wrote " + str(lane_profiles_path))

# Band centres as data, long format, positions only.
csv_column_names = [
    "well_index",
    "lane_centre_stacking_pixels",
    "band_index_in_lane",
    "band_centre_migration_pixels",
    "band_centre_migration_millimetres",
    "window_start_migration_pixels",
    "window_end_migration_pixels",
    "lane_strip_top_stacking_pixels",
    "lane_strip_bottom_stacking_pixels",
    "lane_strip_measured_height_pixels",
    "lane_strip_fixed_height_pixels",
    "peak_profile_value",
    "merged_center_flag",
]
band_centers_csv_lines = [",".join(csv_column_names)]
for band_record in band_centers_csv_rows:
    band_centers_csv_lines.append(
        ",".join(str(band_record[name]) for name in csv_column_names)
    )
band_centers_csv_path = output_directory_path / BAND_CENTERS_CSV_FILENAME
band_centers_csv_path.write_text("\n".join(band_centers_csv_lines) + "\n")
emit_message("csv", "wrote " + str(band_centers_csv_path))

# Saturation derivation panel: the standing check that makes the effective-ceiling
# decision auditable at a glance. Left, the high-DN histogram of the full page and
# of the measured crop with the container and the effective ceiling marked, so a
# pile-up below the container (a knee) or its absence is visible. Right, the profile
# of the brightest crop band through its peak, so a resolved peak (unsaturated) is
# distinguishable from a flat top kissing the ceiling (clipped). Display-only arrays,
# built here from the raw crop, kept separate from the measured arrays above.
saturation_figure, saturation_axes = matplotlib.pyplot.subplots(1, 2, figsize=(13, 5))
high_tail_floor_value = knee_search_floor_value
histogram_bin_edges = numpy.arange(
    high_tail_floor_value, CONTAINER_MAXIMUM_VALUE_16_BIT + 200, 200
)
page_high_tail_values = raw_pixel_array[
    raw_pixel_array >= high_tail_floor_value
].ravel()
crop_high_tail_values = crop_pixels[crop_pixels >= high_tail_floor_value].ravel()
saturation_axes[0].hist(
    page_high_tail_values,
    bins=histogram_bin_edges,
    color="0.6",
    log=True,
    label="full page",
)
saturation_axes[0].hist(
    crop_high_tail_values,
    bins=histogram_bin_edges,
    color="C0",
    alpha=0.8,
    log=True,
    label="measured crop",
)
saturation_axes[0].axvline(
    CONTAINER_MAXIMUM_VALUE_16_BIT,
    color="r",
    linestyle="--",
    linewidth=1.5,
    label="container 65535",
)
saturation_axes[0].axvline(
    effective_ceiling_value,
    color="g",
    linestyle="-",
    linewidth=1.2,
    label="effective ceiling %d" % effective_ceiling_value,
)
saturation_axes[0].set_xlabel("stored DN")
saturation_axes[0].set_ylabel("pixel count (log)")
saturation_axes[0].set_title(
    "high-DN histogram; knee below container: %s"
    % ("yes" if sub_ceiling_knee_detected else "no")
)
saturation_axes[0].legend(fontsize=7, loc="upper center")

brightest_row_index, brightest_column_index = numpy.unravel_index(
    int(numpy.argmax(oriented_raw_stacking_by_migration)),
    oriented_raw_stacking_by_migration.shape,
)
peak_offsets = numpy.arange(-25, 26)
migration_lo = max(0, brightest_column_index - 25)
migration_hi = min(migration_extent_pixels, brightest_column_index + 26)
stacking_lo = max(0, brightest_row_index - 25)
stacking_hi = min(stacking_extent_pixels, brightest_row_index + 26)
migration_slice = oriented_raw_stacking_by_migration[
    brightest_row_index, migration_lo:migration_hi
]
stacking_slice = oriented_raw_stacking_by_migration[
    stacking_lo:stacking_hi, brightest_column_index
]
saturation_axes[1].plot(
    numpy.arange(migration_lo, migration_hi) - brightest_column_index,
    migration_slice,
    "-o",
    markersize=3,
    label="migration profile through peak",
)
saturation_axes[1].plot(
    numpy.arange(stacking_lo, stacking_hi) - brightest_row_index,
    stacking_slice,
    "-s",
    markersize=3,
    label="stacking profile through peak",
)
saturation_axes[1].axhline(
    CONTAINER_MAXIMUM_VALUE_16_BIT,
    color="r",
    linestyle="--",
    linewidth=1.5,
    label="container 65535",
)
saturation_axes[1].axhline(
    effective_ceiling_value,
    color="g",
    linestyle=":",
    linewidth=1.0,
    label="effective ceiling %d" % effective_ceiling_value,
)
saturation_axes[1].set_ylim(0, CONTAINER_MAXIMUM_VALUE_16_BIT * 1.03)
saturation_axes[1].set_xlabel("pixels from brightest peak")
saturation_axes[1].set_ylabel("raw stored DN")
saturation_axes[1].set_title(
    "brightest crop band: peak %d"
    % int(
        oriented_raw_stacking_by_migration[brightest_row_index, brightest_column_index]
    )
)
saturation_axes[1].legend(fontsize=7, loc="upper right")

saturation_figure.suptitle(
    "saturation derivation: dark offset %s, exposure mode %s, %d saturated cell(s)"
    % (str(dark_offset_value), str(exposure_mode_text), saturated_cell_count),
    fontsize=10,
)
saturation_figure.tight_layout()
saturation_derivation_path = output_directory_path / SATURATION_DERIVATION_FILENAME
saturation_figure.savefig(saturation_derivation_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(saturation_figure)
emit_message("saturation", "wrote " + str(saturation_derivation_path))

# Baseline comparison panel: the standing visualisation of the cross-check. One row
# per lane, each showing the raw integrated profile, the valley-to-valley reference
# points that set the local baseline, and the opened (or ball) baseline curve, with
# each canonical band centre marked and coloured by its agreement status. This is
# the operational form of the radius-probe overlay: it makes visible where the two
# baselines part company (the clusters) and where they agree (isolated bands), so a
# fragile flag can be eyeballed rather than trusted blind. Display-only.
comparison_lane_count = len(per_lane_records)
baseline_comparison_figure, baseline_comparison_axes = matplotlib.pyplot.subplots(
    max(1, comparison_lane_count),
    1,
    figsize=(11, 2.3 * max(1, comparison_lane_count)),
    sharex=True,
)
if comparison_lane_count <= 1:
    baseline_comparison_axes = [baseline_comparison_axes]
agreement_status_colour = {
    "agree": "tab:green",
    "fragile": "tab:red",
    "no_signal": "0.6",
}
for lane_axis, lane_record in zip(baseline_comparison_axes, per_lane_records):
    well_index = lane_record["well_index"]
    lane_migration_profile = lane_profile_by_well_index[well_index]
    lane_cross_check_baseline = cross_check_baseline_by_well_index[well_index]
    lane_axis.plot(
        lane_migration_profile,
        color="k",
        linewidth=0.8,
        label="well %d profile" % well_index,
    )
    # Valley-to-valley: the per-window sloped baseline actually integrated (the
    # adaptive, per-band baseline). Each segment is drawn as the full line between
    # its two flanking-minima (valley) reference points so it visibly spans valley
    # to valley, with the narrower sub-span that is actually integrated (the region
    # window) overdrawn thicker, and the valley anchors marked. This is display only.
    first_vtv_segment = True
    for (
        segment_start,
        segment_end,
        left_value,
        right_value,
        left_column,
        right_column,
    ) in valley_to_valley_display_segments_by_well[well_index]:
        if right_column != left_column:

            def baseline_at(columns_array):
                return left_value + (
                    (right_value - left_value)
                    * (columns_array - left_column)
                    / float(right_column - left_column)
                )
        else:

            def baseline_at(columns_array):
                return numpy.full(numpy.shape(columns_array), left_value)

        valley_span_columns = numpy.arange(left_column, right_column + 1)
        window_span_columns = numpy.arange(segment_start, segment_end)
        lane_axis.plot(
            valley_span_columns,
            baseline_at(valley_span_columns),
            color="tab:orange",
            linewidth=0.8,
            alpha=0.5,
            label="valley-to-valley (anchors)" if first_vtv_segment else None,
        )
        lane_axis.plot(
            window_span_columns,
            baseline_at(window_span_columns),
            color="tab:orange",
            linewidth=1.8,
            label="valley-to-valley (integrated window)" if first_vtv_segment else None,
        )
        lane_axis.plot(
            [left_column, right_column],
            [left_value, right_value],
            "o",
            color="tab:orange",
            markersize=2.5,
        )
        first_vtv_segment = False
    lane_axis.plot(
        lane_cross_check_baseline,
        color="tab:blue",
        linewidth=1.0,
        label="%s baseline" % parsed_arguments.baseline_cross_check_method,
    )
    for measurement_row in band_measurement_rows:
        if measurement_row["well_index"] != well_index:
            continue
        if measurement_row["baseline_method"] == "valley_to_valley":
            continue
        band_column = measurement_row["canonical_position_migration_pixels"]
        status_colour = agreement_status_colour.get(
            measurement_row["baseline_agreement_status"], "k"
        )
        lane_axis.axvline(
            band_column, color=status_colour, linewidth=0.8, alpha=0.4, linestyle=":"
        )
    lane_axis.set_ylabel("sum DN")
    lane_axis.legend(fontsize=6, loc="upper right")
baseline_comparison_axes[-1].set_xlabel("migration column (pixels)")
baseline_comparison_figure.suptitle(
    "baseline cross-check (width %.2f mm): orange = valley-to-valley (per-band, "
    "integrated), blue = %s (fixed-width cross-check); dotted ticks green=agree, "
    "red=fragile"
    % (
        parsed_arguments.rolling_ball_width_millimetres,
        parsed_arguments.baseline_cross_check_method,
    ),
    fontsize=9,
)
baseline_comparison_figure.tight_layout()
baseline_comparison_path = output_directory_path / BASELINE_COMPARISON_FILENAME
baseline_comparison_figure.savefig(baseline_comparison_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(baseline_comparison_figure)
emit_message("baseline", "wrote " + str(baseline_comparison_path))

# Shared helpers and per-cell lookups for the results figures below. status_colour
# keys the agreement verdict; cell_by_well_and_band pairs the two baseline rows so a
# metric from either is reachable per cell. All arrays here are display-only.
results_status_colour = {"agree": "tab:green", "fragile": "tab:red", "no_signal": "0.6"}
cell_rows_by_well_and_band = {}
for measurement_row in band_measurement_rows:
    cell_rows_by_well_and_band.setdefault(
        (measurement_row["well_index"], measurement_row["canonical_band_index"]), {}
    )[measurement_row["baseline_method"]] = measurement_row


def draw_lane_trace(target_axis, well_index_for_trace, transpose_axes):
    """Draw one lane's trace with both baselines, windows, and apex markers.

    transpose_axes False plots intensity on y against migration on x (trace read
    left to right). True swaps them so migration runs down the y-axis and intensity
    to the right, matching the gel-native vertical-lane orientation.
    """
    lane_trace_profile = lane_profile_by_well_index[well_index_for_trace]
    migration_axis_values = numpy.arange(lane_trace_profile.size)

    def plot_xy(migration_values, intensity_values, **plot_keywords):
        if transpose_axes:
            target_axis.plot(intensity_values, migration_values, **plot_keywords)
        else:
            target_axis.plot(migration_values, intensity_values, **plot_keywords)

    plot_xy(migration_axis_values, lane_trace_profile, color="k", linewidth=0.8)
    plot_xy(
        migration_axis_values,
        cross_check_baseline_by_well_index[well_index_for_trace],
        color="tab:blue",
        linewidth=1.0,
    )
    for (
        window_start_draw,
        window_end_draw,
        left_value_draw,
        right_value_draw,
        left_column_draw,
        right_column_draw,
    ) in valley_to_valley_display_segments_by_well[well_index_for_trace]:
        anchor_columns = numpy.arange(left_column_draw, right_column_draw + 1)
        window_columns_draw = numpy.arange(window_start_draw, window_end_draw)
        if right_column_draw != left_column_draw:
            slope_draw = (right_value_draw - left_value_draw) / float(
                right_column_draw - left_column_draw
            )
            anchor_line = left_value_draw + slope_draw * (
                anchor_columns - left_column_draw
            )
            window_line = left_value_draw + slope_draw * (
                window_columns_draw - left_column_draw
            )
        else:
            anchor_line = numpy.full(anchor_columns.shape, left_value_draw)
            window_line = numpy.full(window_columns_draw.shape, left_value_draw)
        plot_xy(
            anchor_columns, anchor_line, color="tab:orange", linewidth=0.8, alpha=0.5
        )
        plot_xy(window_columns_draw, window_line, color="tab:orange", linewidth=1.8)
    for opening_measurement_row in band_measurement_rows:
        if opening_measurement_row["well_index"] != well_index_for_trace:
            continue
        if opening_measurement_row["baseline_method"] == "valley_to_valley":
            continue
        if opening_measurement_row["occupancy"] != "locally_detected":
            continue
        apex_colour = results_status_colour.get(
            opening_measurement_row["baseline_agreement_status"], "k"
        )
        apex_marker = "<" if transpose_axes else "v"
        if transpose_axes:
            target_axis.plot(
                opening_measurement_row["peak_apex_raw_value"],
                opening_measurement_row["peak_apex_migration_pixels"],
                apex_marker,
                color=apex_colour,
                markersize=6,
            )
        else:
            target_axis.plot(
                opening_measurement_row["peak_apex_migration_pixels"],
                opening_measurement_row["peak_apex_raw_value"],
                apex_marker,
                color=apex_colour,
                markersize=6,
            )


# Figure: per-lane gel strip (image) over the strip-sum trace, migration aligned.
strip_trace_lane_count = len(per_lane_records)
strip_trace_figure = matplotlib.pyplot.figure(
    figsize=(12, 3.0 * strip_trace_lane_count)
)
strip_trace_grid = strip_trace_figure.add_gridspec(
    2 * strip_trace_lane_count,
    1,
    height_ratios=[1, 3] * strip_trace_lane_count,
    hspace=0.08,
)
for lane_row_index, lane_record in enumerate(per_lane_records):
    well_index = lane_record["well_index"]
    strip_top, strip_bottom = lane_strip_bounds_by_well_index[well_index]
    lane_strip_image = oriented_raw_stacking_by_migration[strip_top:strip_bottom, :]
    strip_axis = strip_trace_figure.add_subplot(strip_trace_grid[2 * lane_row_index, 0])
    strip_axis.imshow(
        lane_strip_image,
        cmap="gray_r",
        aspect="auto",
        vmin=float(numpy.percentile(lane_strip_image, 50)),
        vmax=float(numpy.percentile(lane_strip_image, 99.5)),
        extent=[0, lane_strip_image.shape[1], 0, lane_strip_image.shape[0]],
    )
    strip_axis.set_yticks([])
    strip_axis.set_xticks([])
    strip_axis.set_ylabel(
        "well %d\nstrip" % well_index, fontsize=7, rotation=0, ha="right", va="center"
    )
    trace_axis = strip_trace_figure.add_subplot(
        strip_trace_grid[2 * lane_row_index + 1, 0], sharex=strip_axis
    )
    draw_lane_trace(trace_axis, well_index, transpose_axes=False)
    trace_axis.set_ylabel("strip sum DN", fontsize=7)
trace_axis.set_xlabel("migration column (pixels)")
strip_trace_figure.suptitle(
    "Per-lane gel strip (top, image intensity) over strip-sum trace (bottom). black=profile, "
    "orange=valley-to-valley (thick=integrated window), blue=opening baseline, triangle=apex at "
    "detected bands (green agree / red fragile). Panels differ in units: strip is per-pixel image "
    "intensity, trace is a sum over the strip.",
    fontsize=7.5,
)
strip_trace_path = output_directory_path / STRIP_TRACE_FILENAME
strip_trace_figure.savefig(
    strip_trace_path, dpi=FIGURE_DOTS_PER_INCH, bbox_inches="tight"
)
matplotlib.pyplot.close(strip_trace_figure)
emit_message("figure", "wrote " + str(strip_trace_path))

# The gel-native rotated composite (figure_strip_trace_rotated) is intentionally
# not emitted: its trace did not rotate with the strip, so it was misleading. It is
# deferred to the figure-rework thread, where the rotation will be fixed before it
# is reinstated.

# Figure: per-lane metrics normalised to each lane's own max, vs band index, to show
# whether apex and region net rank and scale the bands alike (they diverge in clusters).
metric_figure, metric_axes = matplotlib.pyplot.subplots(
    strip_trace_lane_count, 1, figsize=(11, 2.2 * strip_trace_lane_count), sharex=True
)
if strip_trace_lane_count == 1:
    metric_axes = [metric_axes]
for metric_axis, lane_record in zip(metric_axes, per_lane_records):
    well_index = lane_record["well_index"]
    band_indices = sorted(
        {key[1] for key in cell_rows_by_well_and_band if key[0] == well_index}
    )
    region_net_values = numpy.array(
        [
            cell_rows_by_well_and_band[(well_index, band_index)]["valley_to_valley"][
                "integrated_intensity"
            ]
            for band_index in band_indices
        ],
        dtype=float,
    )
    apex_above_opening_values = numpy.array(
        [
            cell_rows_by_well_and_band[(well_index, band_index)][
                parsed_arguments.baseline_cross_check_method
            ]["peak_apex_above_baseline"]
            for band_index in band_indices
        ],
        dtype=float,
    )
    apex_raw_values = numpy.array(
        [
            cell_rows_by_well_and_band[(well_index, band_index)]["valley_to_valley"][
                "peak_apex_raw_value"
            ]
            for band_index in band_indices
        ],
        dtype=float,
    )

    def normalise_to_max(values_array):
        maximum_value = values_array.max()
        return values_array / maximum_value if maximum_value > 0 else values_array

    metric_axis.plot(
        band_indices,
        normalise_to_max(region_net_values),
        "-o",
        markersize=3,
        label="region net (vtv)",
    )
    metric_axis.plot(
        band_indices,
        normalise_to_max(apex_above_opening_values),
        "-s",
        markersize=3,
        label="apex above opening",
    )
    metric_axis.plot(
        band_indices,
        normalise_to_max(apex_raw_values),
        "-^",
        markersize=3,
        label="apex raw",
    )
    for band_index in band_indices:
        if (
            cell_rows_by_well_and_band[(well_index, band_index)]["valley_to_valley"][
                "baseline_agreement_status"
            ]
            == "fragile"
        ):
            metric_axis.axvline(band_index, color="tab:red", alpha=0.25)
    metric_axis.set_ylabel("well %d\n(norm)" % well_index, fontsize=7)
    metric_axis.legend(fontsize=6, loc="upper right", ncol=3)
metric_axes[-1].set_xlabel("canonical band index")
metric_figure.suptitle(
    "Per-lane metrics each normalised to that lane's max: do apex and region net rank/scale the "
    "bands alike? They overlay on isolated bands and diverge in clusters. red line = fragile cell.",
    fontsize=8,
)
metric_figure.tight_layout()
metric_compare_path = output_directory_path / METRIC_COMPARE_FILENAME
metric_figure.savefig(metric_compare_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(metric_figure)
emit_message("figure", "wrote " + str(metric_compare_path))

# Figure: per-cell scatter, region net vs apex above opening (log-log), a companion
# view of overall agreement. Cells where the region net clipped to zero cannot sit on
# a log axis, which is itself the point: the per-lane figure shows those explicitly.
scatter_figure, scatter_axis = matplotlib.pyplot.subplots(figsize=(6.5, 6))
for status_name, status_colour in results_status_colour.items():
    scatter_x_values = []
    scatter_y_values = []
    for (well_index, band_index), cell_rows in cell_rows_by_well_and_band.items():
        if cell_rows["valley_to_valley"]["baseline_agreement_status"] != status_name:
            continue
        region_net_value = cell_rows["valley_to_valley"]["integrated_intensity"]
        apex_above_opening_value = cell_rows[
            parsed_arguments.baseline_cross_check_method
        ]["peak_apex_above_baseline"]
        if region_net_value > 0 and apex_above_opening_value > 0:
            scatter_x_values.append(region_net_value)
            scatter_y_values.append(apex_above_opening_value)
    if scatter_x_values:
        scatter_axis.scatter(
            scatter_x_values,
            scatter_y_values,
            s=18,
            color=status_colour,
            alpha=0.8,
            label=status_name,
        )
scatter_axis.set_xscale("log")
scatter_axis.set_yscale("log")
scatter_axis.set_xlabel("region net, valley-to-valley (DN)")
scatter_axis.set_ylabel("apex above opening (DN)")
scatter_axis.set_title(
    "Region net vs apex above opening per cell (log-log); off-diagonal = divergence"
)
scatter_axis.legend(fontsize=7)
scatter_figure.tight_layout()
net_vs_apex_path = output_directory_path / NET_VS_APEX_FILENAME
scatter_figure.savefig(net_vs_apex_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(scatter_figure)
emit_message("figure", "wrote " + str(net_vs_apex_path))

# Overlay: the reported value drawn on the original oriented crop at each detected
# band, coloured by the basis on which it was chosen, so the published number can be
# eyeballed where it was measured.
reported_basis_colour = {
    "region_net": "tab:green",
    "opening_net_cluster_inclusive": "tab:red",
    "saturated_floor": "tab:purple",
}
overlay_reported_figure, overlay_reported_axis = matplotlib.pyplot.subplots(
    figsize=(
        12,
        12.0
        * oriented_raw_stacking_by_migration.shape[0]
        / max(1, oriented_raw_stacking_by_migration.shape[1]),
    )
)
overlay_reported_axis.imshow(
    oriented_raw_stacking_by_migration,
    cmap="gray_r",
    aspect="auto",
    vmin=float(numpy.percentile(oriented_raw_stacking_by_migration, 50)),
    vmax=float(numpy.percentile(oriented_raw_stacking_by_migration, 99.5)),
)
for measurement_row in band_measurement_rows:
    if measurement_row["baseline_method"] == "valley_to_valley":
        continue
    if measurement_row["occupancy"] != "locally_detected":
        continue
    annotation_colour = reported_basis_colour.get(
        measurement_row["reported_value_basis"], "k"
    )
    overlay_reported_axis.text(
        measurement_row["peak_apex_migration_pixels"],
        measurement_row["lane_centre_stacking_pixels"],
        "%.1e" % measurement_row["reported_value"],
        fontsize=5,
        color=annotation_colour,
        ha="center",
        va="center",
        rotation=90 if gel_migration_axis == "vertical" else 0,
    )
overlay_reported_axis.set_xlabel("migration column (pixels)")
overlay_reported_axis.set_ylabel("stacking row (pixels)")
overlay_reported_figure.suptitle(
    "Reported value per detected band on the crop (green=region net, red=opening net "
    "cluster-inclusive, purple=saturated floor). All values are areas (summed DN).",
    fontsize=8,
)
reported_values_overlay_path = output_directory_path / REPORTED_VALUES_OVERLAY_FILENAME
overlay_reported_figure.savefig(
    reported_values_overlay_path, dpi=FIGURE_DOTS_PER_INCH, bbox_inches="tight"
)
matplotlib.pyplot.close(overlay_reported_figure)
emit_message("figure", "wrote " + str(reported_values_overlay_path))


# Single exit point: assemble and write the one provenance report, then set the
# code. Hard stops and warnings are pooled across both halves of the stage.
failed_hard_stops = [
    finding["check_name"]
    for finding in (ANALYSIS_FINDINGS + BAND_DETECTION_FINDINGS)
    if finding["is_hard_stop"] and finding["status"] == "fail"
]
warnings = [
    finding["check_name"]
    for finding in (ANALYSIS_FINDINGS + BAND_DETECTION_FINDINGS)
    if finding["status"] == "warning"
]
overall_status = "fail" if failed_hard_stops else ("warning" if warnings else "pass")

# The crop migration origin in the full-image frame, so postprocessing can convert
# crop-local migration positions to any other frame, including distance-from-well
# once the well side is supplied.
migration_origin_full_image_pixels = (
    crop_x_pixels if gel_migration_axis == "horizontal" else crop_y_pixels
)

tilt_warning_fired = any(
    finding["check_name"] == "tilt_measurement" and finding["status"] == "warning"
    for finding in ANALYSIS_FINDINGS
)

gel_measurement_report = {
    "gel_measurement_report_schema_version": GEL_MEASUREMENT_REPORT_SCHEMA_VERSION,
    "analysis_report_schema_version": ANALYSIS_REPORT_SCHEMA_VERSION,
    "band_detection_report_schema_version": BAND_DETECTION_REPORT_SCHEMA_VERSION,
    "overall_status": overall_status,
    "failed_hard_stop_check_names": failed_hard_stops,
    "warning_check_names": warnings,
    "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "input_file": str(input_tiff_absolute_path),
    "read_from_stage1_report": str(validation_report_path),
    # Carried through untouched so no consumer can forget that the pixel values may
    # not be linear in signal.
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
        "migration_origin_full_image_pixels": migration_origin_full_image_pixels,
    },
    "tilt": {
        "warning_threshold_degrees": parsed_arguments.tilt_warning_threshold_degrees,
        "derived_tilt_angle_degrees": derived_tilt_angle_degrees,
        "warning_fired": tilt_warning_fired,
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
    "configuration_used": {
        "region_window_width_millimetres": parsed_arguments.region_window_width_millimetres,
        "region_window_width_pixels": round(region_window_width_pixels, 2),
        "minimum_band_separation_millimetres": parsed_arguments.minimum_band_separation_millimetres,
        "minimum_band_separation_pixels": round(minimum_band_separation_pixels, 2),
        "band_peak_prominence_fraction": parsed_arguments.band_peak_prominence_fraction,
        "lane_strip_cross_section_fraction": parsed_arguments.lane_strip_cross_section_fraction,
        "use_fixed_strip_height": parsed_arguments.use_fixed_strip_height,
        "lane_strip_fixed_half_height_millimetres": parsed_arguments.lane_strip_fixed_half_height_millimetres,
        "migration_profile_smoothing_pixels": migration_profile_smoothing_pixels,
        "band_cluster_tolerance_pixels": round(cluster_tolerance_pixels, 2),
        "band_cluster_tolerance_source": cluster_tolerance_source,
        "canonical_spread_warn_fraction": CANONICAL_SPREAD_WARN_FRACTION,
        "baseline_flank_search_pixels": baseline_flank_search_pixels,
        "saturation_warn_fraction": saturation_warn_fraction,
        "baseline_cross_check_method": parsed_arguments.baseline_cross_check_method,
        "rolling_ball_width_millimetres": parsed_arguments.rolling_ball_width_millimetres,
        "rolling_ball_width_pixels": round(rolling_ball_width_pixels, 2),
        "band_definition_source": band_definition_source,
    },
    "saturation": {
        "container_maximum_value": CONTAINER_MAXIMUM_VALUE_16_BIT,
        "effective_ceiling_value": effective_ceiling_value,
        "sub_ceiling_knee_detected": sub_ceiling_knee_detected,
        "page_pixels_at_effective_ceiling": page_pixels_at_effective_ceiling,
        "crop_pixels_at_effective_ceiling": crop_pixels_at_effective_ceiling,
        "saturated_cell_count": saturated_cell_count,
        "dark_offset_value": dark_offset_value,
        "exposure_mode_text": exposure_mode_text,
    },
    "baseline_cross_check": {
        "method": parsed_arguments.baseline_cross_check_method,
        "width_millimetres": parsed_arguments.rolling_ball_width_millimetres,
        "width_pixels": round(rolling_ball_width_pixels, 2),
        "fragile_disagreement_fraction": BASELINE_FRAGILE_DISAGREEMENT_FRACTION,
        "baseline_fragile_cell_count": baseline_fragile_cell_count,
        "opening_baseline_elevated_cell_count": sum(
            1
            for measurement_row in band_measurement_rows
            if measurement_row["baseline_method"] != "valley_to_valley"
            and measurement_row["opening_baseline_elevated"]
        ),
        "minimum_lane_climb_into_peak": round(minimum_lane_climb_into_peak, 4),
        "width_adequacy_climb_warn_fraction": WIDTH_ADEQUACY_CLIMB_WARN_FRACTION,
    },
    "peak_mode_and_multiplicity": {
        "sub_peak_prominence_fraction": parsed_arguments.band_peak_prominence_fraction,
        "multi_maximum_window_count": multi_maximum_window_count,
    },
    "total_band_count": total_band_count,
    "per_lane": [
        {
            "well_index": lane_record["well_index"],
            "lane_centre_stacking_pixels": lane_record["lane_centre_stacking_pixels"],
            "strip_top_stacking_pixels": lane_record["strip_top_stacking_pixels"],
            "strip_bottom_stacking_pixels": lane_record["strip_bottom_stacking_pixels"],
            "measured_strip_height_pixels": lane_record["measured_strip_height_pixels"],
            "fixed_strip_height_pixels": lane_record["fixed_strip_height_pixels"],
            "band_count": lane_record["band_count"],
            "band_centre_migration_pixels": lane_record["band_centre_migration_pixels"],
            "merged_center_flags": lane_record["merged_center_flags"],
        }
        for lane_record in per_lane_records
    ],
    "analysis_findings": ANALYSIS_FINDINGS,
    "band_detection_findings": BAND_DETECTION_FINDINGS,
    "consensus_bands": [
        {
            "canonical_band_index": canonical_position_index,
            "canonical_position_migration_pixels": int(
                round(
                    canonical_band_records[canonical_position_index][
                        "canonical_position_migration_pixels"
                    ]
                )
            ),
            "canonical_position_migration_millimetres": round(
                canonical_band_records[canonical_position_index][
                    "canonical_position_migration_pixels"
                ]
                * millimetres_per_pixel,
                4,
            ),
            "window_start_migration_pixels": canonical_band_records[
                canonical_position_index
            ]["window_start_migration_pixels"],
            "window_end_migration_pixels": canonical_band_records[
                canonical_position_index
            ]["window_end_migration_pixels"],
            "member_wells": canonical_band_records[canonical_position_index][
                "member_wells"
            ],
            "cross_lane_spread_pixels": canonical_band_records[
                canonical_position_index
            ]["cross_lane_spread_pixels"],
            "spread_warns": canonical_band_records[canonical_position_index].get(
                "spread_warns", False
            ),
        }
        for canonical_position_index in range(len(canonical_band_records))
    ],
    "integration_summary": {
        "canonical_band_count": canonical_band_count,
        "measured_cell_count": len(band_measurement_rows),
        "locally_detected_cell_count": locally_detected_cell_count,
        "consensus_only_cell_count": consensus_only_cell_count,
    },
    "outputs": {
        "lane_grid_overlay": str(overlay_output_path),
        "lane_profiles_plot": str(lane_profiles_plot_path),
        "lane_profiles_csv": str(lane_profiles_csv_path),
        "band_windows_overlay": str(band_windows_overlay_path),
        "lane_migration_profiles": str(lane_profiles_path),
        "band_centers_csv": str(band_centers_csv_path),
        "band_definitions_csv": str(band_definitions_csv_path),
        "band_measurements_csv": str(band_measurements_csv_path),
        "consensus_center_map": str(consensus_center_map_path),
        "consensus_windows_overlay": str(consensus_windows_overlay_path),
        "saturation_derivation": str(saturation_derivation_path),
        "baseline_comparison": str(baseline_comparison_path),
        "strip_trace": str(strip_trace_path),
        "metric_compare": str(metric_compare_path),
        "net_vs_apex": str(net_vs_apex_path),
        "reported_values_overlay": str(reported_values_overlay_path),
        "band_detection_report": str(
            output_directory_path / BAND_DETECTION_REPORT_FILENAME
        ),
    },
}
band_detection_report_path = output_directory_path / BAND_DETECTION_REPORT_FILENAME
band_detection_report_path.write_text(json.dumps(gel_measurement_report, indent=2))
emit_message(
    "report", "wrote " + str(band_detection_report_path) + "; overall " + overall_status
)

sys.exit(1 if failed_hard_stops else 0)
