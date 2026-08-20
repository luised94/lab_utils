# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Pull one number per lane, across all samples, for a chosen band OR a chosen
migration region. This is the small, traceable step that turns the gel analysis
into the per-sample quantity an experiment needs. It does no scoring and picks no
"best" band: you choose the band or region (by eye, from the metrics and the
picker), and this writes exactly that, with full provenance.

Two mutually exclusive selection modes, matching the ways the pipeline can pull a
per-lane number:

  --band N [--quantity area|apex]
      Take consensus band N straight from band_measurements.csv for every measured
      lane. quantity=area is reported_area (region-net at the shared consensus
      window; smile-affected because the window is fixed in migration).
      quantity=apex is apex_height_above_baseline (the peak height at each lane's
      own apex; the smile-robust hedge).

  --region START_mm END_mm [--net-baseline none|straight]
      Integrate the width-summed profile (manual_lane_profiles.csv raw_value) over
      the migration window [START_mm, END_mm] for every lane. This is the same
      fixed-migration-window idea as a band, over a span you define, and it is
      smile-affected. net-baseline none sums raw_value as-is; straight also
      subtracts a trapezoidal baseline drawn between the two endpoint signal
      values (a local baseline that is only valid when the endpoints sit in
      valleys, i.e. a narrow single-peak span).

Addressing follows the rest of the pipeline: point at the gel analysis directory
(or any file inside it) and the standard filenames resolve. --profile-csv and
--band-measurements override when the convention does not hold.

Outputs, written into the gel directory so they travel with the analysis and are
checked every run:
    extract_band_<N>_<quantity>.csv         (band mode)   one row per lane
    extract_region_<A>-<B>mm_<baseline>.csv (region mode) one row per lane
    <same stem>_checks.json                 inputs, per-lane values, verdicts

Stdlib only. Flat and procedural; emit_message / die / record_check are the only
helpers, duplicated verbatim into each script per the house convention.
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys

# Standard filenames resolved inside the gel analysis directory.
PROFILE_CSV_FILENAME = "manual_lane_profiles.csv"
BAND_MEASUREMENTS_FILENAME = "band_measurements.csv"

# A region narrower than this many migration rows is almost certainly a typo (the
# two millimetre bounds swapped, or start == end) rather than a real window.
MINIMUM_REGION_ROW_COUNT = 2

# A summed profile is already background-subtracted and clipped at zero upstream,
# so a straight endpoint baseline that exceeds the summed signal means the
# endpoints were not valleys; the result is clipped here to avoid a negative area.
ZERO_AREA = 0.0

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
# terminal. The stem, selection, and gel id it needs are set before any check runs.
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
        hard_failure_checks_path = gel_analysis_directory / (
            output_stem + "_checks.json"
        )
        hard_failure_checks_path.write_text(
            json.dumps(
                {
                    "schema_version": "extract_lane_values_1",
                    "generated_at": datetime.datetime.now()
                    .astimezone()
                    .isoformat(timespec="seconds"),
                    "gel_id": gel_id,
                    "selection": selection_description,
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
    prog="extract_lane_values.py",
    description="Extract one per-lane value for a chosen band or migration region.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "gel_path",
    help="the gel analysis directory, or any file inside it",
)
selection_mode_group = argument_parser.add_mutually_exclusive_group(required=True)
selection_mode_group.add_argument(
    "--band",
    type=int,
    metavar="BAND_INDEX",
    help="consensus band index to take from band_measurements.csv",
)
selection_mode_group.add_argument(
    "--region",
    type=float,
    nargs=2,
    metavar=("START_MILLIMETRES", "END_MILLIMETRES"),
    help="migration window in millimetres to integrate from the profile",
)
argument_parser.add_argument(
    "--quantity",
    choices=["area", "apex"],
    default="area",
    help="band mode only: reported_area (default) or apex_height_above_baseline",
)
argument_parser.add_argument(
    "--net-baseline",
    dest="net_baseline",
    choices=["none", "straight"],
    default="none",
    help="region mode only: subtract a straight endpoint-to-endpoint baseline",
)
argument_parser.add_argument("--profile-csv", dest="profile_csv_override", default=None)
argument_parser.add_argument(
    "--band-measurements", dest="band_measurements_override", default=None
)
parsed_arguments = argument_parser.parse_args()

# Resolve the gel analysis directory from the directory itself or any file in it,
# matching the manifest and sheet-validator addressing rule.
given_path = pathlib.Path(parsed_arguments.gel_path)
if given_path.is_dir():
    gel_analysis_directory = given_path
elif given_path.is_file():
    gel_analysis_directory = given_path.parent
else:
    # gel_analysis_directory is referenced by record_check on early failure, so it
    # must exist before the first check; a bad path is fatal here, plainly.
    print("[input] ERROR: path does not exist: " + str(given_path), file=sys.stderr)
    sys.exit(1)

profile_csv_path = (
    pathlib.Path(parsed_arguments.profile_csv_override)
    if parsed_arguments.profile_csv_override
    else gel_analysis_directory / PROFILE_CSV_FILENAME
)
band_measurements_path = (
    pathlib.Path(parsed_arguments.band_measurements_override)
    if parsed_arguments.band_measurements_override
    else gel_analysis_directory / BAND_MEASUREMENTS_FILENAME
)

# =============================================================================
# Read band_measurements.csv. It carries identity (well_number, sample_label,
# analysis_role) per lane as well as the per-band measurements. utf-8-sig strips an
# Excel BOM and every cell is trimmed, because a hand-exported CSV carries a BOM
# and stray whitespace that would break the exact lane/band joins invisibly.
# =============================================================================

if not band_measurements_path.is_file():
    print("[input] ERROR: missing " + str(band_measurements_path), file=sys.stderr)
    sys.exit(1)
with band_measurements_path.open(newline="", encoding="utf-8-sig") as band_csv_handle:
    band_measurement_rows = [
        {
            column_name: (cell_value.strip() if isinstance(cell_value, str) else cell_value)
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in csv.DictReader(band_csv_handle)
    ]
if not band_measurement_rows:
    print("[input] ERROR: band_measurements.csv has no rows", file=sys.stderr)
    sys.exit(1)

gel_id = band_measurement_rows[0]["gel_id"]

identity_by_lane_index = {}
for measurement_row in band_measurement_rows:
    lane_index = int(measurement_row["lane_index"])
    if lane_index not in identity_by_lane_index:
        identity_by_lane_index[lane_index] = {
            "well_number": measurement_row["well_number"],
            "sample_label": measurement_row["sample_label"],
            "analysis_role": measurement_row["analysis_role"],
        }

# A lane present in the profile but absent from band_measurements (a ladder or an
# empty lane the analysis excluded) carries this sentinel rather than a blank, so a
# missing identity is visibly deliberate rather than a forgotten field.
IDENTITY_FOR_LANE_NOT_IN_MEASUREMENTS = {
    "well_number": "",
    "sample_label": "",
    "analysis_role": "not_in_measurements",
}

# =============================================================================
# BAND MODE
# =============================================================================

if parsed_arguments.band is not None:
    band_index = parsed_arguments.band
    quantity_name = parsed_arguments.quantity
    selection_description = {
        "mode": "band",
        "band_index": band_index,
        "quantity": quantity_name,
    }
    output_stem = "extract_band_%d_%s" % (band_index, quantity_name)

    selected_band_rows = [
        measurement_row
        for measurement_row in band_measurement_rows
        if int(measurement_row["canonical_band_index"]) == band_index
    ]
    record_check(
        "band index exists in measurements",
        "hard",
        len(selected_band_rows) > 0,
        "band %d present in %d lane rows" % (band_index, len(selected_band_rows)),
    )

    quantity_column_name = (
        "reported_area" if quantity_name == "area" else "apex_height_above_baseline"
    )
    canonical_position_millimetres = float(
        selected_band_rows[0]["canonical_position_millimetres"]
    )
    consensus_window_start_pixels = int(selected_band_rows[0]["window_start_pixels"])
    consensus_window_end_pixels = int(selected_band_rows[0]["window_end_pixels"])

    output_rows = []
    extracted_values = []
    fragile_lane_count = 0
    for measurement_row in sorted(
        selected_band_rows, key=lambda row: int(row["lane_index"])
    ):
        lane_index = int(measurement_row["lane_index"])
        identity = identity_by_lane_index[lane_index]
        extracted_value = float(measurement_row[quantity_column_name])
        extracted_values.append(extracted_value)
        if measurement_row["baseline_agreement_status"] == "fragile":
            fragile_lane_count += 1
        output_rows.append(
            {
                "gel_id": gel_id,
                "lane_index": lane_index,
                "well_number": identity["well_number"],
                "sample_label": identity["sample_label"],
                "analysis_role": identity["analysis_role"],
                "method": "band",
                "selector": "band_%d" % band_index,
                "quantity": quantity_name,
                "canonical_position_millimetres": measurement_row[
                    "canonical_position_millimetres"
                ],
                "window_start_pixels": measurement_row["window_start_pixels"],
                "window_end_pixels": measurement_row["window_end_pixels"],
                "value": round(extracted_value, 4),
                "reported_value_basis": measurement_row["reported_value_basis"],
                "region_net_area": measurement_row["region_net_area"],
                "opening_net_area": measurement_row["opening_net_area"],
                "apex_height_above_baseline": measurement_row[
                    "apex_height_above_baseline"
                ],
                "baseline_agreement_status": measurement_row[
                    "baseline_agreement_status"
                ],
                "saturation_status": measurement_row["saturation_status"],
                "band_occupancy": measurement_row["band_occupancy"],
                "band_is_well_supported": measurement_row["band_is_well_supported"],
                "cross_lane_spread_millimetres": measurement_row[
                    "cross_lane_spread_millimetres"
                ],
            }
        )

    record_check(
        "all selected values finite",
        "hard",
        all(
            extracted_value == extracted_value
            and abs(extracted_value) != float("inf")
            for extracted_value in extracted_values
        ),
        "%d values" % len(extracted_values),
    )
    record_check(
        "band is well supported",
        "soft",
        selected_band_rows[0]["band_is_well_supported"] == "yes",
        "band %d occupancy %s, well_supported=%s"
        % (
            band_index,
            selected_band_rows[0]["band_occupancy"],
            selected_band_rows[0]["band_is_well_supported"],
        ),
    )
    record_check(
        "baseline agreement across lanes",
        "soft",
        fragile_lane_count == 0,
        "%d of %d lanes fragile at band %d"
        % (fragile_lane_count, len(selected_band_rows), band_index),
    )
    selection_provenance = {
        "band_index": band_index,
        "quantity": quantity_name,
        "quantity_column": quantity_column_name,
        "canonical_position_millimetres": canonical_position_millimetres,
        "window_start_pixels": consensus_window_start_pixels,
        "window_end_pixels": consensus_window_end_pixels,
        "note": (
            "quantity=area is smile-affected (fixed migration window); "
            "quantity=apex is the smile-robust per-lane peak height"
        ),
    }
    profile_csv_used_path = None

# =============================================================================
# REGION MODE
# =============================================================================

else:
    region_start_millimetres, region_end_millimetres = parsed_arguments.region
    if region_end_millimetres < region_start_millimetres:
        region_start_millimetres, region_end_millimetres = (
            region_end_millimetres,
            region_start_millimetres,
        )
    net_baseline_name = parsed_arguments.net_baseline
    selection_description = {
        "mode": "region",
        "start_millimetres": region_start_millimetres,
        "end_millimetres": region_end_millimetres,
        "net_baseline": net_baseline_name,
    }
    output_stem = "extract_region_%g-%gmm_%s" % (
        region_start_millimetres,
        region_end_millimetres,
        net_baseline_name,
    )

    if not profile_csv_path.is_file():
        record_check(
            "profile csv present", "hard", False, "missing " + str(profile_csv_path)
        )
    with profile_csv_path.open(newline="", encoding="utf-8-sig") as profile_csv_handle:
        profile_rows = [
            {
                column_name: (
                    cell_value.strip() if isinstance(cell_value, str) else cell_value
                )
                for column_name, cell_value in raw_row.items()
            }
            for raw_row in csv.DictReader(profile_csv_handle)
        ]
    record_check(
        "profile csv has rows", "hard", len(profile_rows) > 0, str(profile_csv_path)
    )

    # Group the profile by lane, each an ordered list of (migration_pixels,
    # migration_millimetres, summed_signal). The background median is constant
    # across the file; keep it for provenance.
    profile_samples_by_lane_index = {}
    plate_background_median = None
    for profile_row in profile_rows:
        lane_index = int(profile_row["lane_index"])
        profile_samples_by_lane_index.setdefault(lane_index, []).append(
            (
                int(profile_row["migration_position_pixels"]),
                float(profile_row["migration_position_millimetres"]),
                float(profile_row["raw_value"]),
            )
        )
        if plate_background_median is None:
            plate_background_median = float(profile_row["plate_background_median"])
    for lane_index in profile_samples_by_lane_index:
        profile_samples_by_lane_index[lane_index].sort(
            key=lambda migration_sample: migration_sample[0]
        )

    # Which consensus bands fall inside the region, for context in the report.
    consensus_band_indices_in_region = sorted(
        {
            int(measurement_row["canonical_band_index"])
            for measurement_row in band_measurement_rows
            if region_start_millimetres
            <= float(measurement_row["canonical_position_millimetres"])
            <= region_end_millimetres
        }
    )

    output_rows = []
    extracted_values = []
    representative_window_row_count = None
    baseline_overshoot_lane_count = 0
    for lane_index in sorted(profile_samples_by_lane_index):
        migration_samples = profile_samples_by_lane_index[lane_index]
        samples_in_window = [
            (migration_pixels, migration_millimetres, summed_signal)
            for (migration_pixels, migration_millimetres, summed_signal) in migration_samples
            if region_start_millimetres
            <= migration_millimetres
            <= region_end_millimetres
        ]
        if not samples_in_window:
            extracted_value = None
            raw_window_sum = None
            baseline_subtracted_area = ZERO_AREA
            window_row_count = 0
            window_start_pixels = None
            window_end_pixels = None
        else:
            window_row_count = len(samples_in_window)
            window_start_pixels = samples_in_window[0][0]
            window_end_pixels = samples_in_window[-1][0]
            raw_window_sum = sum(
                summed_signal for (_, _, summed_signal) in samples_in_window
            )
            if net_baseline_name == "straight":
                left_endpoint_signal = samples_in_window[0][2]
                right_endpoint_signal = samples_in_window[-1][2]
                # trapezoidal area under the straight endpoint-to-endpoint line
                baseline_subtracted_area = (
                    0.5
                    * (left_endpoint_signal + right_endpoint_signal)
                    * window_row_count
                )
            else:
                baseline_subtracted_area = ZERO_AREA
            extracted_value = raw_window_sum - baseline_subtracted_area
            if extracted_value < ZERO_AREA:
                extracted_value = ZERO_AREA
                if raw_window_sum > ZERO_AREA:
                    baseline_overshoot_lane_count += 1
        if representative_window_row_count is None and window_row_count:
            representative_window_row_count = window_row_count
        if extracted_value is not None:
            extracted_values.append(extracted_value)
        identity = identity_by_lane_index.get(
            lane_index, IDENTITY_FOR_LANE_NOT_IN_MEASUREMENTS
        )
        output_rows.append(
            {
                "gel_id": gel_id,
                "lane_index": lane_index,
                "well_number": identity["well_number"],
                "sample_label": identity["sample_label"],
                "analysis_role": identity["analysis_role"],
                "method": "region",
                "selector": "region_%g-%gmm"
                % (region_start_millimetres, region_end_millimetres),
                "net_baseline": net_baseline_name,
                "region_start_millimetres": region_start_millimetres,
                "region_end_millimetres": region_end_millimetres,
                "window_start_pixels": "" if window_start_pixels is None else window_start_pixels,
                "window_end_pixels": "" if window_end_pixels is None else window_end_pixels,
                "row_count": window_row_count,
                "raw_sum": "" if raw_window_sum is None else round(raw_window_sum, 4),
                "baseline_subtracted": round(baseline_subtracted_area, 4),
                "value": "" if extracted_value is None else round(extracted_value, 4),
                "plate_background_median": plate_background_median,
            }
        )

    record_check(
        "region spans enough rows",
        "hard",
        (representative_window_row_count or 0) >= MINIMUM_REGION_ROW_COUNT,
        "region [%g, %g] mm covers ~%s rows"
        % (
            region_start_millimetres,
            region_end_millimetres,
            representative_window_row_count,
        ),
    )
    record_check(
        "all lanes have signal in the region",
        "soft",
        all(output_row["value"] != "" for output_row in output_rows),
        "%d of %d lanes empty in region"
        % (
            sum(1 for output_row in output_rows if output_row["value"] == ""),
            len(output_rows),
        ),
    )
    record_check(
        "straight baseline did not overshoot",
        "soft",
        baseline_overshoot_lane_count == 0,
        "%d lanes clipped to zero by the straight baseline (endpoints are not "
        "valleys; widen to valleys or use --net-baseline none)"
        % baseline_overshoot_lane_count,
    )
    selection_provenance = {
        "region_start_millimetres": region_start_millimetres,
        "region_end_millimetres": region_end_millimetres,
        "net_baseline": net_baseline_name,
        "consensus_bands_in_region": consensus_band_indices_in_region,
        "plate_background_median": plate_background_median,
        "note": (
            "raw_value is the width-summed, background-subtracted profile; a fixed "
            "migration window is smile-affected. net-baseline=straight removes a "
            "local endpoint-to-endpoint baseline and clips negatives at zero"
        ),
    }
    profile_csv_used_path = str(profile_csv_path)

# =============================================================================
# Summary statistics over the extracted values (the spread is what "is this
# consistent?" reduces to once a band or region is chosen).
# =============================================================================

present_values = [
    (float(output_row["value"]) if output_row["value"] != "" else None)
    if isinstance(output_row["value"], str)
    else output_row["value"]
    for output_row in output_rows
]
present_values = [value for value in present_values if value is not None]
if present_values:
    ordered_values = sorted(present_values)
    value_count = len(ordered_values)
    mean_value = sum(ordered_values) / value_count
    if value_count % 2:
        median_value = ordered_values[value_count // 2]
    else:
        median_value = 0.5 * (
            ordered_values[value_count // 2 - 1] + ordered_values[value_count // 2]
        )
    if value_count > 1 and mean_value != 0:
        sample_variance = sum(
            (value - mean_value) ** 2 for value in ordered_values
        ) / (value_count - 1)
        coefficient_of_variation = (sample_variance ** 0.5) / abs(mean_value)
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
else:
    value_summary = {"n": 0}

# =============================================================================
# Write the per-lane CSV and the checks JSON.
# =============================================================================

output_csv_path = gel_analysis_directory / (output_stem + ".csv")
output_column_names = list(output_rows[0].keys())
with output_csv_path.open("w", newline="", encoding="utf-8") as output_csv_handle:
    csv_writer = csv.DictWriter(output_csv_handle, fieldnames=output_column_names)
    csv_writer.writeheader()
    csv_writer.writerows(output_rows)
emit_message("output", "wrote " + str(output_csv_path))

checks_path = gel_analysis_directory / (output_stem + "_checks.json")
checks_path.write_text(
    json.dumps(
        {
            "schema_version": "extract_lane_values_1",
            "generated_at": datetime.datetime.now()
            .astimezone()
            .isoformat(timespec="seconds"),
            "gel_id": gel_id,
            "inputs": {
                "band_measurements": str(band_measurements_path),
                "profile_csv": profile_csv_used_path,
            },
            "selection": selection_description,
            "provenance": selection_provenance,
            "value_summary_over_all_output_lanes": value_summary,
            "checks": accumulated_check_records,
        },
        indent=2,
    )
    + "\n"
)
emit_message("output", "wrote " + str(checks_path))

# Console summary to stderr (stdout stays empty so the CSV path can be diffed):
# the per-lane values and their spread.
emit_message("summary", "gel " + gel_id)
emit_message("summary", "selection: " + json.dumps(selection_description))
for output_row in output_rows:
    emit_message(
        "value",
        "lane %2s well %2s role %-20s %-14s value %s"
        % (
            output_row["lane_index"],
            output_row["well_number"],
            output_row["analysis_role"],
            output_row["sample_label"],
            output_row["value"],
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
