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
Stage 3 of the gel densitometry pipeline, Slice B1: per-lane band-centre
detection only. Reads the Slice A geometry report (stage2_analysis_report.json)
for the crop, the migration axis, the pixel size, and the located lane centres,
then for each lane collapses the lane's own strip into a migration profile and
calls the band centres as prominent peaks of that profile. It places a region
window around each centre and draws everything for the operator to check.

It integrates nothing and computes no baseline or intensity; that is Slice B2/B3.
B1 exists to prove that the band centres and windows land where the bands are
before any number is measured on them.

Design decisions this slice commits to (see DESIGN.md sections 5.1, 5.2, 5.3,
5.4, and 14, and the s0002 measurement session):

1. Every lane is treated as an independent sample. Bands are detected per lane
   from that lane's own profile, NOT from a grid shared across lanes. On this gel
   the lanes carry two different protein-complex patterns, so a shared grid would
   manufacture phantom bands (DESIGN 5.2). Correspondence between lanes is a
   postprocessing concern, joined by band migration position and metadata, not by
   an ordinal index assigned here. This supersedes the vertical, shared-grid
   wording of DESIGN 5.1/5.2 for the axis-aware pipeline.
2. Band positions are reported in crop-local migration pixels and millimetres.
   The crop migration origin, the pixel size, and the axis are all carried in the
   report, so postprocessing can convert to distance-from-well or any other frame
   without re-running. Distance-from-well needs the well side, which is deferred.
3. Each lane's strip is bounded at the midpoints to its neighbouring lanes, so
   adjacent lanes whose doublets are wider than the comb pitch cannot bleed into
   each other or be double-counted. Within that bound the strip is measured from
   the lane's own stacking cross-section (its doublet width), which varies lane to
   lane. A command-line override forces a single fixed strip height instead, and
   both the measured and the fixed height are always written, so the adaptivity is
   auditable and pinnable.
4. The flip to top-left origin and lane/band numbering are NOT applied here.
   Integration is invariant to the flip, ordering does not matter for this slice,
   and applying the flip now would only relabel columns that postprocessing owns.

Shape follows CONVENTIONS.md: flat procedural, no helpers beyond the tagged stderr
emitter and the fail-fast exit, full descriptive names carrying units, ASCII only,
comments state why. Boilerplate is duplicated from stage 2 rather than shared,
because the two scripts are a pipeline for now and will be fused later; a shared
module before the second real reuse would be premature abstraction.

Usage:

    uv run measure_gel_bands.py '<scan>.tif'
    uv run measure_gel_bands.py --region-window-width-millimetres 2.5 '<scan>.tif'
    uv run measure_gel_bands.py --use-fixed-strip-height '<scan>.tif'
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


BAND_DETECTION_REPORT_SCHEMA_VERSION = 1

STAGE2_REPORT_FILENAME = "stage2_analysis_report.json"
BAND_DETECTION_REPORT_FILENAME = "band_detection_report.json"
BAND_WINDOWS_OVERLAY_FILENAME = "overlay_band_windows_labeled.png"
LANE_MIGRATION_PROFILES_FILENAME = "lane_migration_profiles.png"
BAND_CENTERS_CSV_FILENAME = "band_centers.csv"

CONTAINER_MAXIMUM_VALUE_16_BIT = 65535

# A migration column elevated down its whole length is a plate edge or scratch,
# not a lane; exclude columns whose median exceeds this multiple of the typical
# column median. Same rule and constant as the Slice A locator, kept identical so
# the two scripts see the same signal.
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

DISPLAY_COLORMAP_NAME = "gray_r"
FIGURE_DOTS_PER_INCH = 130


def emit_message(source_tag, message_text):
    sys.stderr.write("[" + source_tag + "] " + message_text + "\n")


def die(source_tag, message_text):
    emit_message(source_tag, message_text)
    sys.exit(2)


argument_parser = argparse.ArgumentParser(
    description="Stage 3 Slice B1: per-lane band-centre detection and windows, "
                "drawn and dumped. Reads the Slice A geometry report. No integration."
)
argument_parser.add_argument("input_tiff_path")
argument_parser.add_argument("--page-index", type=int, default=0)
argument_parser.add_argument(
    "--region-window-width-millimetres", type=float,
    default=DEFAULT_REGION_WINDOW_WIDTH_MILLIMETRES,
)
argument_parser.add_argument(
    "--minimum-band-separation-millimetres", type=float,
    default=DEFAULT_MINIMUM_BAND_SEPARATION_MILLIMETRES,
)
argument_parser.add_argument(
    "--band-peak-prominence-fraction", type=float,
    default=DEFAULT_BAND_PEAK_PROMINENCE_FRACTION,
)
argument_parser.add_argument(
    "--lane-strip-cross-section-fraction", type=float,
    default=DEFAULT_LANE_STRIP_CROSS_SECTION_FRACTION,
)
argument_parser.add_argument(
    "--use-fixed-strip-height", action="store_true",
    help="use one fixed strip height for every lane instead of the per-lane "
         "measured doublet width",
)
argument_parser.add_argument(
    "--lane-strip-fixed-half-height-millimetres", type=float,
    default=DEFAULT_LANE_STRIP_FIXED_HALF_HEIGHT_MILLIMETRES,
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
        + str(output_directory_path) + ". Run stage 1 and stage 2 first.")

stage2_report_path = output_directory_path / STAGE2_REPORT_FILENAME
if not stage2_report_path.is_file():
    die("stage2", "no Slice A report at " + str(stage2_report_path)
        + ". Run measure_gel_image.py first; this stage reads its lane geometry.")
try:
    stage2_report = json.loads(stage2_report_path.read_text())
except Exception as report_read_error:
    die("stage2", "Slice A report will not parse: " + str(report_read_error))

geometry_used = stage2_report.get("geometry_used")
lane_grid = stage2_report.get("lane_grid")
if geometry_used is None or lane_grid is None:
    die("stage2", "the Slice A report is missing geometry_used or lane_grid; it may "
        "be from an older schema. Re-run stage 2.")

gel_migration_axis = geometry_used["gel_migration_axis"]
crop_x_pixels = int(geometry_used["crop_x_pixels"])
crop_y_pixels = int(geometry_used["crop_y_pixels"])
crop_width_pixels = int(geometry_used["crop_width_pixels"])
crop_height_pixels = int(geometry_used["crop_height_pixels"])
micrometres_per_pixel = geometry_used["micrometres_per_pixel"]
if micrometres_per_pixel is None:
    die("stage2", "the Slice A report has no micrometres_per_pixel, so millimetre "
        "config values cannot be converted to pixels. A scan with no pixel size "
        "needs one supplied before band windows can be sized.")

lane_centres_stacking_pixels = lane_grid.get("lane_centres_stacking_pixels")
lane_well_indices = lane_grid.get("lane_well_indices")
lane_pitch_pixels = float(lane_grid.get("pitch_pixels"))
if not lane_centres_stacking_pixels:
    die("stage2", "the Slice A report found no populated lanes, so there is nothing "
        "to measure bands in. Check the Slice A overlay before proceeding.")
if lane_well_indices is None or len(lane_well_indices) != len(lane_centres_stacking_pixels):
    # Not fatal: well indices are provisional and only used as a label here.
    lane_well_indices = list(range(len(lane_centres_stacking_pixels)))

# Convert the millimetre config to pixels once, using the report's pixel size.
millimetres_per_pixel = micrometres_per_pixel / 1000.0
region_window_width_pixels = parsed_arguments.region_window_width_millimetres / millimetres_per_pixel
minimum_band_separation_pixels = max(
    1.0, parsed_arguments.minimum_band_separation_millimetres / millimetres_per_pixel
)
migration_profile_smoothing_pixels = max(
    1, int(round(DEFAULT_MIGRATION_PROFILE_SMOOTHING_MILLIMETRES / millimetres_per_pixel))
)
lane_strip_fixed_half_height_pixels = int(round(
    parsed_arguments.lane_strip_fixed_half_height_millimetres / millimetres_per_pixel
))

BAND_DETECTION_FINDINGS = []

# The read path, float64 for headroom, cropped exactly as Slice A cropped.
try:
    with tifffile.TiffFile(str(input_tiff_absolute_path)) as tiff_handle:
        raw_pixel_array = tiff_handle.pages[parsed_arguments.page_index].asarray()
except Exception as pixel_read_error:
    die("pixels", "page " + str(parsed_arguments.page_index) + " would not decode: "
        + str(pixel_read_error))
pixel_array = raw_pixel_array.astype(numpy.float64)
crop_pixels = pixel_array[
    crop_y_pixels:crop_y_pixels + crop_height_pixels,
    crop_x_pixels:crop_x_pixels + crop_width_pixels,
]

# Re-derive background and the edge-column mask exactly as Slice A did, rather than
# trust a passed array, so the two scripts cannot silently diverge on preprocessing.
plate_background_value = float(numpy.median(crop_pixels))
signal_above_plate = crop_pixels - plate_background_value
signal_above_plate[signal_above_plate < 0.0] = 0.0

# Orient so axis 0 is the stacking axis (lanes) and axis 1 is the migration axis
# (bands), whatever the scan orientation, so the per-lane logic is written once.
if gel_migration_axis == "horizontal":
    stacking_by_migration_signal = signal_above_plate
else:
    stacking_by_migration_signal = signal_above_plate.T
stacking_extent_pixels = stacking_by_migration_signal.shape[0]
migration_extent_pixels = stacking_by_migration_signal.shape[1]

per_migration_column_median = numpy.median(stacking_by_migration_signal, axis=0)
positive_column_medians = per_migration_column_median[per_migration_column_median > 0.0]
typical_column_median = (
    float(numpy.median(positive_column_medians)) if positive_column_medians.size else 1.0
)
per_migration_column_saturated_fraction = numpy.mean(
    (crop_pixels if gel_migration_axis == "horizontal" else crop_pixels.T)
    >= CONTAINER_MAXIMUM_VALUE_16_BIT, axis=0
)
usable_migration_column_mask = (
    (per_migration_column_median <= EDGE_COLUMN_MEDIAN_MULTIPLE * typical_column_median)
    & (per_migration_column_saturated_fraction < 0.05)
)
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
        lower_bound_stacking = 0.5 * (ordered_lane_centres[lane_position - 1] + lane_centre_stacking)
    else:
        lower_bound_stacking = lane_centre_stacking - lane_pitch_pixels
    if lane_position < len(ordered_lane_centres) - 1:
        upper_bound_stacking = 0.5 * (lane_centre_stacking + ordered_lane_centres[lane_position + 1])
    else:
        upper_bound_stacking = lane_centre_stacking + lane_pitch_pixels
    search_top = max(0, int(math.floor(lower_bound_stacking)))
    search_bottom = min(stacking_extent_pixels, int(math.ceil(upper_bound_stacking)))

    # The lane's stacking cross-section within its bounded territory, collapsed over
    # the usable migration columns. Its extent above a fraction of its own maximum
    # is the measured doublet strip.
    territory_signal = usable_stacking_by_migration_signal[search_top:search_bottom, :]
    lane_cross_section = territory_signal.sum(axis=1)
    cross_section_maximum = float(lane_cross_section.max()) if lane_cross_section.size else 0.0
    if cross_section_maximum > 0.0:
        cross_section_threshold = parsed_arguments.lane_strip_cross_section_fraction * cross_section_maximum
        rows_above_threshold = numpy.where(lane_cross_section > cross_section_threshold)[0]
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
    fixed_strip_top = max(search_top, int(round(lane_centre_stacking)) - lane_strip_fixed_half_height_pixels)
    fixed_strip_bottom = min(search_bottom, int(round(lane_centre_stacking)) + lane_strip_fixed_half_height_pixels)
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
    band_prominence_threshold = parsed_arguments.band_peak_prominence_fraction * lane_profile_maximum
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
            if other_position != band_position and abs(
                int(band_centre_columns[band_position]) - int(band_centre_columns[other_position])
            ) < region_window_width_pixels:
                merged_center_flags[band_position] = True

    lane_band_records = []
    for band_position in range(len(band_centre_columns)):
        band_centre_migration = int(band_centre_columns[band_position])
        window_start_migration = max(0, int(round(band_centre_migration - half_window_pixels)))
        window_end_migration = min(migration_extent_pixels, int(round(band_centre_migration + half_window_pixels)))
        band_record = {
            "well_index": well_index,
            "lane_centre_stacking_pixels": round(lane_centre_stacking, 2),
            "band_index_in_lane": band_position,
            "band_centre_migration_pixels": band_centre_migration,
            "band_centre_migration_millimetres": round(band_centre_migration * millimetres_per_pixel, 4),
            "window_start_migration_pixels": window_start_migration,
            "window_end_migration_pixels": window_end_migration,
            "lane_strip_top_stacking_pixels": strip_top,
            "lane_strip_bottom_stacking_pixels": strip_bottom,
            "lane_strip_measured_height_pixels": measured_strip_height_pixels,
            "lane_strip_fixed_height_pixels": fixed_strip_height_pixels,
            "peak_profile_value": round(float(lane_migration_profile_smoothed[band_centre_migration]), 2),
            "merged_center_flag": bool(merged_center_flags[band_position]),
        }
        lane_band_records.append(band_record)
        band_centers_csv_rows.append(band_record)

    if len(band_centre_columns) == 0:
        BAND_DETECTION_FINDINGS.append({
            "check_name": "lane_has_detectable_bands",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "lane at well index %d (stacking %d) yielded no band peaks above "
                      "the prominence floor; it may be empty or faint. Check the "
                      "profile panel." % (well_index, int(round(lane_centre_stacking))),
        })

    per_lane_records.append({
        "well_index": well_index,
        "lane_centre_stacking_pixels": round(lane_centre_stacking, 2),
        "strip_top_stacking_pixels": strip_top,
        "strip_bottom_stacking_pixels": strip_bottom,
        "measured_strip_height_pixels": measured_strip_height_pixels,
        "fixed_strip_height_pixels": fixed_strip_height_pixels,
        "band_count": len(band_centre_columns),
        "band_centre_migration_pixels": [int(value) for value in band_centre_columns],
        "merged_center_flags": merged_center_flags,
        "profile_smoothed": lane_migration_profile_smoothed,
    })

total_band_count = len(band_centers_csv_rows)

# Overlay for the operator: display-scaled crop on gray_r so bands read dark and it
# matches Fiji. Each lane's strip is drawn as two boundary lines; each band window
# is a box spanning the strip, with the centre marked. Display array is separate
# from the measured signal per DESIGN 5.9.
display_source = signal_above_plate
display_normalized = numpy.clip(
    display_source / numpy.percentile(display_source, 99.5), 0.0, 1.0
)
overlay_figure, overlay_axes = matplotlib.pyplot.subplots(figsize=(9.0, 9.0))
overlay_axes.imshow(display_normalized, cmap=DISPLAY_COLORMAP_NAME,
                    interpolation="nearest", origin="upper")
for lane_record in per_lane_records:
    strip_top = lane_record["strip_top_stacking_pixels"]
    strip_bottom = lane_record["strip_bottom_stacking_pixels"]
    for band_record_position in range(len(lane_record["band_centre_migration_pixels"])):
        band_centre_migration = lane_record["band_centre_migration_pixels"][band_record_position]
        merged = lane_record["merged_center_flags"][band_record_position]
        window_start = max(0, int(round(band_centre_migration - region_window_width_pixels / 2.0)))
        window_width = int(round(region_window_width_pixels))
        box_edge_color = "orange" if merged else "red"
        if gel_migration_axis == "horizontal":
            box = matplotlib.patches.Rectangle(
                (window_start, strip_top), window_width, strip_bottom - strip_top,
                fill=False, edgecolor=box_edge_color, linewidth=0.8,
            )
            overlay_axes.add_patch(box)
            overlay_axes.plot([band_centre_migration, band_centre_migration],
                              [strip_top, strip_bottom], color="blue", linewidth=0.5)
        else:
            box = matplotlib.patches.Rectangle(
                (strip_top, window_start), strip_bottom - strip_top, window_width,
                fill=False, edgecolor=box_edge_color, linewidth=0.8,
            )
            overlay_axes.add_patch(box)
            overlay_axes.plot([strip_top, strip_bottom],
                              [band_centre_migration, band_centre_migration],
                              color="blue", linewidth=0.5)
overlay_axes.set_title(
    "Band windows per lane (B1, positions only): " + input_tiff_absolute_path.name
    + "\nmigration " + gel_migration_axis
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
    die("overlay", "could not write " + str(band_windows_overlay_path) + ": " + str(overlay_write_error))
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
    panel_axes.plot(numpy.arange(profile_values.size), profile_values,
                    color="black", linewidth=0.8)
    for band_record_position in range(len(lane_record["band_centre_migration_pixels"])):
        band_centre_migration = lane_record["band_centre_migration_pixels"][band_record_position]
        panel_axes.axvline(band_centre_migration, color="blue", linewidth=0.6)
        panel_axes.axvspan(
            max(0, band_centre_migration - region_window_width_pixels / 2.0),
            band_centre_migration + region_window_width_pixels / 2.0,
            color="red", alpha=0.12,
        )
    panel_axes.set_title("well index %d (stacking %d): %d bands"
                         % (lane_record["well_index"],
                            int(round(lane_record["lane_centre_stacking_pixels"])),
                            lane_record["band_count"]),
                         fontsize="small")
    panel_axes.set_xlabel("migration position (pixels, crop-local)")
profiles_figure.tight_layout()
lane_profiles_path = output_directory_path / LANE_MIGRATION_PROFILES_FILENAME
profiles_figure.savefig(lane_profiles_path, dpi=FIGURE_DOTS_PER_INCH)
matplotlib.pyplot.close(profiles_figure)
emit_message("profiles", "wrote " + str(lane_profiles_path))

# Band centres as data, long format, positions only.
csv_column_names = [
    "well_index", "lane_centre_stacking_pixels", "band_index_in_lane",
    "band_centre_migration_pixels", "band_centre_migration_millimetres",
    "window_start_migration_pixels", "window_end_migration_pixels",
    "lane_strip_top_stacking_pixels", "lane_strip_bottom_stacking_pixels",
    "lane_strip_measured_height_pixels", "lane_strip_fixed_height_pixels",
    "peak_profile_value", "merged_center_flag",
]
band_centers_csv_lines = [",".join(csv_column_names)]
for band_record in band_centers_csv_rows:
    band_centers_csv_lines.append(",".join(str(band_record[name]) for name in csv_column_names))
band_centers_csv_path = output_directory_path / BAND_CENTERS_CSV_FILENAME
band_centers_csv_path.write_text("\n".join(band_centers_csv_lines) + "\n")
emit_message("csv", "wrote " + str(band_centers_csv_path))

# Single exit point: assemble and write the provenance report, then set the code.
failed_hard_stops = [
    finding["check_name"] for finding in BAND_DETECTION_FINDINGS
    if finding["is_hard_stop"] and finding["status"] == "fail"
]
warnings = [
    finding["check_name"] for finding in BAND_DETECTION_FINDINGS
    if finding["status"] == "warning"
]
overall_status = "fail" if failed_hard_stops else ("warning" if warnings else "pass")

# The crop migration origin in the full-image frame, so postprocessing can convert
# crop-local migration positions to any other frame, including distance-from-well
# once the well side is supplied.
migration_origin_full_image_pixels = (
    crop_x_pixels if gel_migration_axis == "horizontal" else crop_y_pixels
)

band_detection_report = {
    "band_detection_report_schema_version": BAND_DETECTION_REPORT_SCHEMA_VERSION,
    "overall_status": overall_status,
    "failed_hard_stop_check_names": failed_hard_stops,
    "warning_check_names": warnings,
    "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "input_file": str(input_tiff_absolute_path),
    "read_from_stage2_report": str(stage2_report_path),
    "geometry_used": {
        "gel_migration_axis": gel_migration_axis,
        "crop_x_pixels": crop_x_pixels,
        "crop_y_pixels": crop_y_pixels,
        "crop_width_pixels": crop_width_pixels,
        "crop_height_pixels": crop_height_pixels,
        "micrometres_per_pixel": micrometres_per_pixel,
        "migration_origin_full_image_pixels": migration_origin_full_image_pixels,
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
    "band_detection_findings": BAND_DETECTION_FINDINGS,
    "outputs": {
        "band_windows_overlay": str(band_windows_overlay_path),
        "lane_migration_profiles": str(lane_profiles_path),
        "band_centers_csv": str(band_centers_csv_path),
    },
}
band_detection_report_path = output_directory_path / BAND_DETECTION_REPORT_FILENAME
band_detection_report_path.write_text(json.dumps(band_detection_report, indent=2))
emit_message("report", "wrote " + str(band_detection_report_path) + "; overall " + overall_status)

sys.exit(1 if failed_hard_stops else 0)
