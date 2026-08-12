# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
#
# Date created: 2026-08-11
# Purpose: a SENSITIVITY harness for the gel densitometry pipeline. It re-runs the
# measurement chain on one gel under small, deliberate perturbations of the two
# operator-controlled inputs that the pipeline does not otherwise vary -- the crop
# rectangle (authored in Fiji, carried in the _preprocess.txt sidecar) and the
# baseline estimation knobs (stage-2 CLI) -- and reports how much each lane's
# normalized value moves as a RANGE. It answers "how much does the reported number
# depend on choices a human made", not "what is the number".
#
# This harness is OPTIONAL QC. It is not part of the path to a figure. It never
# modifies the repo, the uploaded inputs, or the real analysis directory: every run
# happens in a private scratch tree of per-cell directories, each holding its own
# copy of the image and its own (possibly perturbed) sidecar.
#
# It does NOT run R. The per-lane normalized value is recomputed here directly from
# the stage-2 band_measurements.csv, replicating the R normalization exactly:
#   per-lane total = sum over DISTINCT (well_index, canonical_band_index) of
#                    reported_value  (the two baseline rows carry a duplicated
#                    reported_value; the distinct() de-dups them before summing),
#   normalized     = 100 * per-lane total / reference-lane total,
# with the reference lane taken from the sample sheet role column (well_index =
# well_number - 1 under the default 'direct' flip). This duplication is deliberate;
# if gel_loading_normalization.R changes how it normalizes, this must be revisited.
#
# The gel measured here (s0002) has encoding_verified=FALSE, so every number is
# PROVISIONAL. The harness reports sensitivity of provisional numbers; it does not
# make them verified.
#
# Usage:
#   uv run sensitivity_harness.py <input_tiff> \
#       --pipeline-directory <dir with prepare/validate/measure scripts> \
#       --output-directory <scratch dir> \
#       [--sample-sheet <path>]
#
# The sidecar is located next to the image as <stem>_preprocess.txt unless the on-
# disk name and the sidecar's measured_against_input_filename disagree, in which
# case the harness writes a name-corrected base copy into scratch (scratch only).

import argparse
import csv
import pathlib
import shutil
import subprocess
import sys

# ==============================================================================
# Perturbation grid (edit these; one-factor-at-a-time around a nominal)
# ==============================================================================
# Crop nudges in pixels applied to the sidecar crop rectangle, as
# (label, delta_crop_x, delta_crop_y, delta_crop_width, delta_crop_height).
# The nominal (0,0,0,0) is always measured. Keep nudges within realistic hand-crop
# jitter so stage 1 still detects the expected lane count.
CROP_PERTURBATIONS_PIXELS = [
    ("crop_nominal", 0, 0, 0, 0),
    ("crop_shift_right_6", 6, 0, 0, 0),
    ("crop_shift_down_6", 0, 6, 0, 0),
    ("crop_grow_8", -4, -4, 8, 8),
    ("crop_shrink_8", 4, 4, -8, -8),
]
# Baseline knobs swept one at a time around the stage-2 defaults
# (flank 1.0 mm, rolling-ball width 4.0 mm).
BASELINE_FLANK_SEARCH_MILLIMETRES_SWEEP = [0.75, 1.0, 1.5]
ROLLING_BALL_WIDTH_MILLIMETRES_SWEEP = [3.0, 4.0, 5.0]
NOMINAL_FLANK_SEARCH_MILLIMETRES = 1.0
NOMINAL_ROLLING_BALL_WIDTH_MILLIMETRES = 4.0

ANALYSIS_DIRECTORY_SUFFIX = "_gel_analysis"
PREPROCESS_SIDECAR_SUFFIX = "_preprocess.txt"
SAMPLE_SHEET_SUFFIX = "_sample_sheet.csv"
BAND_MEASUREMENTS_FILENAME = "band_measurements.csv"
EMPTY_ROLE = "empty"


def log(tag, text):
    sys.stderr.write("[" + tag + "] " + text + "\n")


def fail(tag, text):
    sys.stderr.write("[" + tag + "] " + text + "\n")
    raise SystemExit(2)


# ==============================================================================
# Arguments and paths
# ==============================================================================
argument_parser = argparse.ArgumentParser(
    description="Sensitivity harness: perturb crop and baseline knobs on one gel "
    "and report per-lane normalized-value spread as a range."
)
argument_parser.add_argument("input_tiff_path")
argument_parser.add_argument("--pipeline-directory", required=True)
argument_parser.add_argument("--output-directory", required=True)
argument_parser.add_argument("--sample-sheet", default=None)
parsed_arguments = argument_parser.parse_args()

input_tiff_path = pathlib.Path(parsed_arguments.input_tiff_path).resolve()
if not input_tiff_path.is_file():
    fail("input", "input image not found: " + str(input_tiff_path))
pipeline_directory = pathlib.Path(parsed_arguments.pipeline_directory).resolve()
for stage_script_name in ("prepare_gel_image.py", "validate_gel_image.py", "measure_gel.py"):
    if not (pipeline_directory / stage_script_name).is_file():
        fail("pipeline", "missing stage script in pipeline directory: " + stage_script_name)
output_directory = pathlib.Path(parsed_arguments.output_directory).resolve()
output_directory.mkdir(parents=True, exist_ok=True)

image_stem = input_tiff_path.stem
image_filename = input_tiff_path.name

base_sidecar_path = input_tiff_path.with_name(image_stem + PREPROCESS_SIDECAR_SUFFIX)
if not base_sidecar_path.is_file():
    fail("sidecar", "preprocess sidecar not found next to image: " + str(base_sidecar_path))

if parsed_arguments.sample_sheet is not None:
    sample_sheet_path = pathlib.Path(parsed_arguments.sample_sheet).resolve()
else:
    sample_sheet_path = input_tiff_path.with_name(image_stem + SAMPLE_SHEET_SUFFIX)
if not sample_sheet_path.is_file():
    fail("sample_sheet", "sample sheet not found (needed for the reference lane): "
         + str(sample_sheet_path))

# ==============================================================================
# Reference lane from the sample sheet (R-free; well_index = well_number - 1)
# ==============================================================================
role_by_well_index = {}
reference_well_index = None
with open(sample_sheet_path, newline="") as sample_sheet_file:
    for sample_sheet_row in csv.DictReader(sample_sheet_file):
        physical_well_number = int(sample_sheet_row["well_number"])
        well_index = physical_well_number - 1
        role_text = (sample_sheet_row.get("role") or "").strip()
        role_by_well_index[well_index] = role_text
        if role_text == "reference":
            if reference_well_index is not None:
                fail("sample_sheet", "more than one reference lane in the sample sheet.")
            reference_well_index = well_index
if reference_well_index is None:
    fail("sample_sheet", "no reference lane (role=reference) in the sample sheet.")
log("reference", "reference well_index = " + str(reference_well_index))

# ==============================================================================
# Base sidecar text, with the measured_against_input_filename corrected to the on-
# disk name if they disagree (scratch only; the repo and uploads are untouched).
# ==============================================================================
base_sidecar_lines = base_sidecar_path.read_text().splitlines()
corrected_sidecar_lines = []
sidecar_crop_values = {}
for sidecar_line in base_sidecar_lines:
    if sidecar_line.startswith("measured_against_input_filename="):
        corrected_sidecar_lines.append("measured_against_input_filename=" + image_filename)
        continue
    for crop_field_name in ("crop_x", "crop_y", "crop_width", "crop_height"):
        if sidecar_line.startswith(crop_field_name + "="):
            sidecar_crop_values[crop_field_name] = int(sidecar_line.split("=", 1)[1])
    corrected_sidecar_lines.append(sidecar_line)
for crop_field_name in ("crop_x", "crop_y", "crop_width", "crop_height"):
    if crop_field_name not in sidecar_crop_values:
        fail("sidecar", "sidecar is missing crop field: " + crop_field_name)


def write_cell_sidecar(cell_directory, delta_crop_x, delta_crop_y, delta_crop_width, delta_crop_height):
    perturbed_lines = []
    for sidecar_line in corrected_sidecar_lines:
        replaced = False
        for crop_field_name, delta_value in (
            ("crop_x", delta_crop_x), ("crop_y", delta_crop_y),
            ("crop_width", delta_crop_width), ("crop_height", delta_crop_height),
        ):
            if sidecar_line.startswith(crop_field_name + "="):
                new_value = sidecar_crop_values[crop_field_name] + delta_value
                perturbed_lines.append(crop_field_name + "=" + str(new_value))
                replaced = True
                break
        if not replaced:
            perturbed_lines.append(sidecar_line)
    (cell_directory / (image_stem + PREPROCESS_SIDECAR_SUFFIX)).write_text(
        "\n".join(perturbed_lines) + "\n"
    )


# ==============================================================================
# Run one full chain into a crop-cell directory (stage 0 -> stage 1). Returns the
# analysis directory, or None if a stage failed.
# ==============================================================================
def build_crop_cell(cell_label, delta_crop_x, delta_crop_y, delta_crop_width, delta_crop_height):
    cell_directory = output_directory / cell_label
    if cell_directory.exists():
        shutil.rmtree(cell_directory)
    cell_directory.mkdir(parents=True)
    cell_image_path = cell_directory / image_filename
    shutil.copy2(input_tiff_path, cell_image_path)  # real copy: stage 2 resolves symlinks
    write_cell_sidecar(cell_directory, delta_crop_x, delta_crop_y, delta_crop_width, delta_crop_height)
    shutil.copy2(sample_sheet_path, cell_directory / (image_stem + SAMPLE_SHEET_SUFFIX))

    stage0 = subprocess.run(
        ["uv", "run", "prepare_gel_image.py", "--output-parent-directory",
         str(cell_directory), str(cell_image_path)],
        cwd=str(pipeline_directory), capture_output=True, text=True,
    )
    if stage0.returncode != 0:
        log("cell_failed", cell_label + ": stage 0 failed; skipping. tail: "
            + stage0.stderr.strip().splitlines()[-1] if stage0.stderr.strip() else cell_label)
        return None
    stage1 = subprocess.run(
        ["uv", "run", "validate_gel_image.py", "--output-parent-directory",
         str(cell_directory), str(cell_image_path)],
        cwd=str(pipeline_directory), capture_output=True, text=True,
    )
    if stage1.returncode != 0:
        log("cell_failed", cell_label + ": stage 1 failed; skipping.")
        return None
    return cell_directory, cell_image_path


# ==============================================================================
# Run stage 2 with given baseline knobs; return per-lane de-duped totals or None.
# ==============================================================================
def measure_and_total(cell_image_path, flank_search_millimetres, rolling_ball_width_millimetres):
    stage2 = subprocess.run(
        ["uv", "run", "measure_gel.py", str(cell_image_path),
         "--baseline-flank-search-millimetres", str(flank_search_millimetres),
         "--rolling-ball-width-millimetres", str(rolling_ball_width_millimetres)],
        cwd=str(pipeline_directory), capture_output=True, text=True,
    )
    analysis_directory = cell_image_path.parent / (image_stem + ANALYSIS_DIRECTORY_SUFFIX)
    band_measurements_path = analysis_directory / BAND_MEASUREMENTS_FILENAME
    if stage2.returncode != 0 or not band_measurements_path.is_file():
        return None, None
    seen_band_keys = set()
    per_lane_total = {}
    with open(band_measurements_path, newline="") as band_measurements_file:
        for band_row in csv.DictReader(band_measurements_file):
            well_index = int(band_row["well_index"])
            band_key = (well_index, band_row["canonical_band_index"])
            if band_key in seen_band_keys:
                continue
            seen_band_keys.add(band_key)
            per_lane_total[well_index] = per_lane_total.get(well_index, 0.0) + float(band_row["reported_value"])
    lane_grid = tuple(sorted(per_lane_total.keys()))
    return per_lane_total, lane_grid


# ==============================================================================
# Build the cell list (one-factor-at-a-time), run them, collect normalized values
# ==============================================================================
# Each collected record: (cell_label, factor, level, well_index, role, lane_total,
# normalized_value). normalized_value is normalized to the reference lane WITHIN
# that same cell (each cell has its own reference measurement).
collected_records = []
failed_cells = []
lane_grid_changed_cells = []
nominal_lane_grid = None


def collect_from_totals(cell_label, factor, level, measurement):
    global nominal_lane_grid
    per_lane_total, lane_grid = measurement
    if per_lane_total is None or reference_well_index not in per_lane_total:
        failed_cells.append((cell_label, factor, str(level)))
        return
    # A crop nudge can change how many lanes stage 1/2 detect; when it does, the
    # well_index labels reindex and are no longer comparable lane-to-lane across
    # cells. Such cells are EXCLUDED from the per-lane spread and reported as a
    # separate, coarser finding (this crop choice destabilizes lane detection).
    if nominal_lane_grid is not None and lane_grid != nominal_lane_grid:
        lane_grid_changed_cells.append(
            (cell_label, factor, str(level), len(nominal_lane_grid), len(lane_grid))
        )
        return
    reference_total = per_lane_total[reference_well_index]
    if reference_total <= 0:
        failed_cells.append((cell_label, factor, str(level) + " (nonpositive reference)"))
        return
    for well_index in sorted(per_lane_total):
        normalized_value = 100.0 * per_lane_total[well_index] / reference_total
        collected_records.append((
            cell_label, factor, str(level), well_index,
            role_by_well_index.get(well_index, "?"),
            per_lane_total[well_index], normalized_value,
        ))


# Nominal crop cell (reused for all baseline sweeps).
log("build", "building nominal crop cell and baseline sweeps")
nominal_build = build_crop_cell("crop_nominal", 0, 0, 0, 0)
if nominal_build is None:
    fail("nominal", "the nominal crop cell failed to build; cannot proceed.")
nominal_cell_directory, nominal_cell_image_path = nominal_build

# Nominal measurement establishes the reference lane grid all other cells match.
nominal_measurement = measure_and_total(
    nominal_cell_image_path,
    NOMINAL_FLANK_SEARCH_MILLIMETRES, NOMINAL_ROLLING_BALL_WIDTH_MILLIMETRES,
)
if nominal_measurement[0] is None or reference_well_index not in nominal_measurement[0]:
    fail("nominal", "nominal measurement failed; cannot establish the lane grid.")
nominal_lane_grid = nominal_measurement[1]
log("nominal", "nominal lane grid has " + str(len(nominal_lane_grid)) + " lanes")
collect_from_totals("crop_nominal", "nominal", "nominal", nominal_measurement)
# Flank sweep (crop nominal, rolling-ball nominal).
for flank_value in BASELINE_FLANK_SEARCH_MILLIMETRES_SWEEP:
    collect_from_totals(
        "crop_nominal", "flank_search_millimetres", flank_value,
        measure_and_total(nominal_cell_image_path, flank_value, NOMINAL_ROLLING_BALL_WIDTH_MILLIMETRES),
    )
# Rolling-ball width sweep (crop nominal, flank nominal).
for rolling_ball_value in ROLLING_BALL_WIDTH_MILLIMETRES_SWEEP:
    collect_from_totals(
        "crop_nominal", "rolling_ball_width_millimetres", rolling_ball_value,
        measure_and_total(nominal_cell_image_path, NOMINAL_FLANK_SEARCH_MILLIMETRES, rolling_ball_value),
    )

# Crop sweep (each nudge is its own full chain; baseline nominal).
for cell_label, dx, dy, dwidth, dheight in CROP_PERTURBATIONS_PIXELS:
    if cell_label == "crop_nominal":
        continue
    crop_build = build_crop_cell(cell_label, dx, dy, dwidth, dheight)
    if crop_build is None:
        failed_cells.append((cell_label, "crop", "build"))
        continue
    _, crop_cell_image_path = crop_build
    collect_from_totals(
        cell_label, "crop_pixels",
        "dx=%d,dy=%d,dw=%d,dh=%d" % (dx, dy, dwidth, dheight),
        measure_and_total(crop_cell_image_path,
                          NOMINAL_FLANK_SEARCH_MILLIMETRES, NOMINAL_ROLLING_BALL_WIDTH_MILLIMETRES),
    )

if not collected_records:
    fail("empty", "no cells produced measurements; nothing to summarise.")

# ==============================================================================
# Write the tidy per-cell values, the per-lane spread summary, and a short report
# ==============================================================================
per_cell_path = output_directory / "per_cell_lane_values.csv"
with open(per_cell_path, "w", newline="") as per_cell_file:
    writer = csv.writer(per_cell_file)
    writer.writerow(["cell_label", "factor", "level", "well_index", "role",
                     "lane_total_reported_value", "normalized_value_percent_of_reference"])
    for record in collected_records:
        writer.writerow(record)

# Per-lane spread across ALL cells.
normalized_values_by_well = {}
for record in collected_records:
    well_index = record[3]
    normalized_values_by_well.setdefault(well_index, []).append(record[6])

nominal_normalized_by_well = {}
for record in collected_records:
    if record[0] == "crop_nominal" and record[1] == "nominal":
        nominal_normalized_by_well[record[3]] = record[6]

summary_rows = []
worst_case_range_non_empty = 0.0
worst_case_well = None
for well_index in sorted(normalized_values_by_well):
    values = normalized_values_by_well[well_index]
    minimum_value = min(values)
    maximum_value = max(values)
    value_range = maximum_value - minimum_value
    nominal_value = nominal_normalized_by_well.get(well_index, float("nan"))
    fractional_range = (value_range / nominal_value) if nominal_value not in (0.0, float("nan")) else float("nan")
    role_text = role_by_well_index.get(well_index, "?")
    summary_rows.append([well_index, role_text, len(values), nominal_value,
                         minimum_value, maximum_value, value_range, fractional_range])
    if role_text != EMPTY_ROLE and value_range > worst_case_range_non_empty:
        worst_case_range_non_empty = value_range
        worst_case_well = well_index

summary_path = output_directory / "lane_sensitivity_summary.csv"
with open(summary_path, "w", newline="") as summary_file:
    writer = csv.writer(summary_file)
    writer.writerow(["well_index", "role", "cell_count", "nominal_normalized_value",
                     "minimum_normalized_value", "maximum_normalized_value",
                     "normalized_value_range", "fractional_range"])
    for summary_row in summary_rows:
        writer.writerow(summary_row)

report_path = output_directory / "sensitivity_report.md"
report_lines = []
report_lines.append("# Sensitivity report (PROVISIONAL)")
report_lines.append("")
report_lines.append("Gel: " + image_stem)
report_lines.append("")
report_lines.append("This gel's numbers are provisional (encoding_verified is FALSE upstream); "
                    "this harness reports how much they move, not whether they are correct.")
report_lines.append("")
report_lines.append("Perturbations applied, one factor at a time around a nominal:")
report_lines.append("  - crop rectangle nudges (pixels): "
                    + ", ".join(label for label, *_ in CROP_PERTURBATIONS_PIXELS))
report_lines.append("  - baseline flank search (mm): "
                    + ", ".join(str(v) for v in BASELINE_FLANK_SEARCH_MILLIMETRES_SWEEP))
report_lines.append("  - rolling-ball width (mm): "
                    + ", ".join(str(v) for v in ROLLING_BALL_WIDTH_MILLIMETRES_SWEEP))
report_lines.append("")
report_lines.append("Metric: per-lane normalized value (percent of the reference lane), recomputed "
                    "R-free from stage-2 band_measurements.csv. Spread is reported as the range "
                    "(max minus min) over all cells.")
report_lines.append("")
if worst_case_well is not None:
    report_lines.append("HEADLINE: worst-case normalized-value range across non-empty lanes is "
                        + ("%.3f" % worst_case_range_non_empty)
                        + " percent-of-reference (well_index " + str(worst_case_well) + ").")
report_lines.append("")
report_lines.append("Per-lane spread is in lane_sensitivity_summary.csv; every cell value is in "
                    "per_cell_lane_values.csv.")
if lane_grid_changed_cells:
    report_lines.append("")
    report_lines.append("Cells EXCLUDED from the per-lane spread because they changed the detected "
                        "lane count (the well_index labels reindex, so lane-to-lane comparison is "
                        "invalid). This is itself a sensitivity finding: these input choices "
                        "destabilize lane detection on this gel.")
    for changed_cell in lane_grid_changed_cells:
        report_lines.append("  - " + changed_cell[0] + " (" + changed_cell[1] + "=" + changed_cell[2]
                            + "): lanes " + str(changed_cell[3]) + " -> " + str(changed_cell[4]))
if failed_cells:
    report_lines.append("")
    report_lines.append("Cells that did not produce a measurement (excluded from the spread):")
    for failed_cell in failed_cells:
        report_lines.append("  - " + " / ".join(failed_cell))
report_path.write_text("\n".join(report_lines) + "\n")

log("done", "wrote " + per_cell_path.name + ", " + summary_path.name + ", " + report_path.name)
log("done", "worst-case non-empty normalized range: "
    + ("%.3f" % worst_case_range_non_empty) + " percent-of-reference")
