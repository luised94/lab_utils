# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Pull one number per lane, across all samples, for a chosen band OR a chosen
migration region. This is the small, traceable step that turns the gel analysis
into the per-sample quantity an experiment needs. It does no scoring and picks no
"best" band: you choose the band or region (by eye, from the metrics and the
visualization), and this writes exactly that, with full provenance.

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
      fixed-migration-window idea as a band, but over a span you define. It is the
      "sum the peaks over the region" method, and it is smile-affected. baseline
      none sums raw_value as-is; straight also subtracts a trapezoidal baseline
      drawn between the two endpoint signal values (a local, transparent baseline).

Addressing follows the rest of the pipeline: point at the gel analysis directory
(or any file inside it) and the standard filenames resolve. --profile-csv and
--band-measurements override when the convention does not hold.

Outputs, written into the gel directory so they travel with the analysis and are
checked every run:
    extract_band_<N>_<quantity>.csv        (band mode)   one row per lane
    extract_region_<A>-<B>mm_<baseline>.csv (region mode) one row per lane
    <same stem>_checks.json                 inputs, per-lane values, verdicts

No non-standard-library dependency. Flat and procedural, ASCII only, matching the
other pipeline scripts (emit_message/die the only helpers).
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys

# =============================================================================
# Standard filenames resolved inside the gel analysis directory.
# =============================================================================

PROFILE_CSV_FILENAME = "manual_lane_profiles.csv"
BAND_MEASUREMENTS_FILENAME = "band_measurements.csv"

# A region narrower than this many migration rows is almost certainly a typo
# (mm swapped, or start==end) rather than a real integration window.
MINIMUM_REGION_ROW_COUNT = 2

RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


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
selection_group = argument_parser.add_mutually_exclusive_group(required=True)
selection_group.add_argument(
    "--band",
    type=int,
    metavar="N",
    help="consensus band index to take from band_measurements.csv",
)
selection_group.add_argument(
    "--region",
    type=float,
    nargs=2,
    metavar=("START_mm", "END_mm"),
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
argument_parser.add_argument("--profile-csv", dest="profile_csv", default=None)
argument_parser.add_argument(
    "--band-measurements", dest="band_measurements", default=None
)
parsed_arguments = argument_parser.parse_args()


# =============================================================================
# Resolve input files inside the directory (or from explicit overrides).
# =============================================================================


def resolve_directory(given_path):
    path_object = pathlib.Path(given_path)
    if path_object.is_dir():
        return path_object
    if path_object.is_file():
        return path_object.parent
    die("input", "path does not exist: " + str(given_path))


gel_directory = resolve_directory(parsed_arguments.gel_path)

profile_csv_path = (
    pathlib.Path(parsed_arguments.profile_csv)
    if parsed_arguments.profile_csv
    else gel_directory / PROFILE_CSV_FILENAME
)
band_measurements_path = (
    pathlib.Path(parsed_arguments.band_measurements)
    if parsed_arguments.band_measurements
    else gel_directory / BAND_MEASUREMENTS_FILENAME
)


def read_csv_rows(csv_path):
    if not csv_path.is_file():
        die("input", "missing required file: " + str(csv_path))
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


# =============================================================================
# Load identity (from band_measurements: it carries well/label/role per lane).
# =============================================================================

band_rows = read_csv_rows(band_measurements_path)
if not band_rows:
    die("input", "band_measurements.csv has no rows")

gel_id = band_rows[0]["gel_id"]

identity_by_lane = {}
for row in band_rows:
    lane_index_value = int(row["lane_index"])
    if lane_index_value not in identity_by_lane:
        identity_by_lane[lane_index_value] = {
            "well_number": row["well_number"],
            "sample_label": row["sample_label"],
            "analysis_role": row["analysis_role"],
        }


def identity_for(lane_index_value):
    return identity_by_lane.get(
        lane_index_value,
        {
            "well_number": "",
            "sample_label": "",
            "analysis_role": "not_in_measurements",
        },
    )


# =============================================================================
# Checks accumulation (same spirit as profile_checks.json: name, verdict, detail).
# =============================================================================

check_records = []


def record_check(check_name, severity, passed, detail_text):
    check_records.append(
        {
            "check": check_name,
            "severity": severity,
            "passed": bool(passed),
            "detail": detail_text,
        }
    )
    if severity == "hard" and not passed:
        write_checks_and_exit(check_name, detail_text)


def summarize_values(values):
    present = [v for v in values if v is not None]
    if not present:
        return {"n": 0}
    ordered = sorted(present)
    count = len(ordered)
    mean_value = sum(ordered) / count
    if count % 2:
        median_value = ordered[count // 2]
    else:
        median_value = 0.5 * (ordered[count // 2 - 1] + ordered[count // 2])
    if count > 1 and mean_value != 0:
        variance = sum((v - mean_value) ** 2 for v in ordered) / (count - 1)
        coefficient_of_variation = (variance ** 0.5) / abs(mean_value)
    else:
        coefficient_of_variation = None
    return {
        "n": count,
        "minimum": min(ordered),
        "median": median_value,
        "maximum": max(ordered),
        "coefficient_of_variation": (
            round(coefficient_of_variation, 4)
            if coefficient_of_variation is not None
            else None
        ),
    }


output_csv_path = None


def write_checks_and_exit(failed_check_name, failed_detail):
    # Write the checks file even on hard failure, so the run is diagnosable from
    # the artifact and not just the terminal (project convention).
    checks_path = gel_directory / (checks_stem + "_checks.json")
    checks_path.write_text(
        json.dumps(
            {
                "schema_version": "extract_lane_values_1",
                "generated_at": datetime.datetime.now()
                .astimezone()
                .isoformat(timespec="seconds"),
                "gel_id": gel_id,
                "selection": selection_description,
                "checks": check_records,
                "failed": {"check": failed_check_name, "detail": failed_detail},
            },
            indent=2,
        )
        + "\n"
    )
    die("checks", "hard check failed: " + failed_check_name + " (" + failed_detail + ")")


# =============================================================================
# BAND MODE
# =============================================================================

if parsed_arguments.band is not None:
    band_index_value = parsed_arguments.band
    quantity_name = parsed_arguments.quantity
    selection_description = {
        "mode": "band",
        "band_index": band_index_value,
        "quantity": quantity_name,
    }
    checks_stem = "extract_band_%d_%s" % (band_index_value, quantity_name)
    selected = [r for r in band_rows if int(r["canonical_band_index"]) == band_index_value]
    record_check(
        "band index exists in measurements",
        "hard",
        len(selected) > 0,
        "band %d present in %d lane rows" % (band_index_value, len(selected)),
    )

    quantity_column = (
        "reported_area" if quantity_name == "area" else "apex_height_above_baseline"
    )
    canonical_mm = float(selected[0]["canonical_position_millimetres"])
    window_start = int(selected[0]["window_start_pixels"])
    window_end = int(selected[0]["window_end_pixels"])

    output_rows = []
    values_for_summary = []
    fragile_lane_count = 0
    for row in sorted(selected, key=lambda r: int(r["lane_index"])):
        lane_index_value = int(row["lane_index"])
        identity = identity_for(lane_index_value)
        value = float(row[quantity_column])
        values_for_summary.append(value)
        if row["baseline_agreement_status"] == "fragile":
            fragile_lane_count += 1
        output_rows.append(
            {
                "gel_id": gel_id,
                "lane_index": lane_index_value,
                "well_number": identity["well_number"],
                "sample_label": identity["sample_label"],
                "analysis_role": identity["analysis_role"],
                "method": "band",
                "selector": "band_%d" % band_index_value,
                "quantity": quantity_name,
                "canonical_position_millimetres": row["canonical_position_millimetres"],
                "window_start_pixels": row["window_start_pixels"],
                "window_end_pixels": row["window_end_pixels"],
                "value": round(value, 4),
                "reported_value_basis": row["reported_value_basis"],
                "region_net_area": row["region_net_area"],
                "opening_net_area": row["opening_net_area"],
                "apex_height_above_baseline": row["apex_height_above_baseline"],
                "baseline_agreement_status": row["baseline_agreement_status"],
                "saturation_status": row["saturation_status"],
                "band_occupancy": row["band_occupancy"],
                "band_is_well_supported": row["band_is_well_supported"],
                "cross_lane_spread_millimetres": row["cross_lane_spread_millimetres"],
            }
        )

    record_check(
        "all selected values finite",
        "hard",
        all(v == v and abs(v) != float("inf") for v in values_for_summary),
        "%d values" % len(values_for_summary),
    )
    record_check(
        "band is well supported",
        "soft",
        selected[0]["band_is_well_supported"] == "yes",
        "band %d occupancy %s, well_supported=%s"
        % (band_index_value, selected[0]["band_occupancy"], selected[0]["band_is_well_supported"]),
    )
    record_check(
        "baseline agreement across lanes",
        "soft",
        fragile_lane_count == 0,
        "%d of %d lanes fragile at band %d" % (fragile_lane_count, len(selected), band_index_value),
    )
    provenance = {
        "band_index": band_index_value,
        "quantity": quantity_name,
        "quantity_column": quantity_column,
        "canonical_position_millimetres": canonical_mm,
        "window_start_pixels": window_start,
        "window_end_pixels": window_end,
        "note": (
            "quantity=area is smile-affected (fixed migration window); "
            "quantity=apex is the smile-robust per-lane peak height"
        ),
    }

# =============================================================================
# REGION MODE
# =============================================================================

else:
    region_start_mm, region_end_mm = parsed_arguments.region
    if region_end_mm < region_start_mm:
        region_start_mm, region_end_mm = region_end_mm, region_start_mm
    baseline_name = parsed_arguments.net_baseline
    selection_description = {
        "mode": "region",
        "start_millimetres": region_start_mm,
        "end_millimetres": region_end_mm,
        "net_baseline": baseline_name,
    }
    checks_stem = "extract_region_%g-%gmm_%s" % (
        region_start_mm,
        region_end_mm,
        baseline_name,
    )

    profile_rows = read_csv_rows(profile_csv_path)
    record_check(
        "profile csv has rows", "hard", len(profile_rows) > 0, str(profile_csv_path)
    )

    # Group the profile by lane, keeping (migration_pixels, mm, raw_value) ordered.
    profile_by_lane = {}
    plate_background_median = None
    for row in profile_rows:
        lane_index_value = int(row["lane_index"])
        profile_by_lane.setdefault(lane_index_value, []).append(
            (
                int(row["migration_position_pixels"]),
                float(row["migration_position_millimetres"]),
                float(row["raw_value"]),
            )
        )
        if plate_background_median is None:
            plate_background_median = float(row["plate_background_median"])
    for lane_index_value in profile_by_lane:
        profile_by_lane[lane_index_value].sort(key=lambda t: t[0])

    # Report which consensus bands fall inside the region, for context.
    bands_in_region = sorted(
        {
            int(r["canonical_band_index"])
            for r in band_rows
            if region_start_mm
            <= float(r["canonical_position_millimetres"])
            <= region_end_mm
        }
    )

    output_rows = []
    values_for_summary = []
    representative_row_count = None
    baseline_clipped_lane_count = 0
    for lane_index_value in sorted(profile_by_lane):
        samples = profile_by_lane[lane_index_value]
        in_window = [(px, mm, val) for (px, mm, val) in samples
                     if region_start_mm <= mm <= region_end_mm]
        if not in_window:
            value = None
            raw_sum = None
            baseline_sum = 0.0
            row_count = 0
            start_px = end_px = None
        else:
            row_count = len(in_window)
            start_px = in_window[0][0]
            end_px = in_window[-1][0]
            raw_sum = sum(val for (_, _, val) in in_window)
            if baseline_name == "straight":
                left_value = in_window[0][2]
                right_value = in_window[-1][2]
                # trapezoidal area under the straight endpoint-to-endpoint line
                baseline_sum = 0.5 * (left_value + right_value) * row_count
            else:
                baseline_sum = 0.0
            value = raw_sum - baseline_sum
            if value < 0:
                value = 0.0
                if raw_sum > 0:
                    baseline_clipped_lane_count += 1
        if representative_row_count is None and row_count:
            representative_row_count = row_count
        if value is not None:
            values_for_summary.append(value)
        identity = identity_for(lane_index_value)
        output_rows.append(
            {
                "gel_id": gel_id,
                "lane_index": lane_index_value,
                "well_number": identity["well_number"],
                "sample_label": identity["sample_label"],
                "analysis_role": identity["analysis_role"],
                "method": "region",
                "selector": "region_%g-%gmm" % (region_start_mm, region_end_mm),
                "net_baseline": baseline_name,
                "region_start_millimetres": region_start_mm,
                "region_end_millimetres": region_end_mm,
                "window_start_pixels": "" if start_px is None else start_px,
                "window_end_pixels": "" if end_px is None else end_px,
                "row_count": row_count,
                "raw_sum": "" if raw_sum is None else round(raw_sum, 4),
                "baseline_subtracted": round(baseline_sum, 4),
                "value": "" if value is None else round(value, 4),
                "plate_background_median": plate_background_median,
            }
        )

    record_check(
        "region spans enough rows",
        "hard",
        (representative_row_count or 0) >= MINIMUM_REGION_ROW_COUNT,
        "region [%g, %g] mm covers ~%s rows"
        % (region_start_mm, region_end_mm, representative_row_count),
    )
    record_check(
        "all lanes have signal in the region",
        "soft",
        all(r["value"] != "" for r in output_rows),
        "%d of %d lanes empty in region"
        % (sum(1 for r in output_rows if r["value"] == ""), len(output_rows)),
    )
    record_check(
        "straight baseline did not overshoot",
        "soft",
        baseline_clipped_lane_count == 0,
        "%d lanes clipped to zero by the straight baseline (endpoints are not "
        "valleys; widen the valleys or use --net-baseline none)"
        % baseline_clipped_lane_count,
    )
    provenance = {
        "region_start_millimetres": region_start_mm,
        "region_end_millimetres": region_end_mm,
        "net_baseline": baseline_name,
        "consensus_bands_in_region": bands_in_region,
        "plate_background_median": plate_background_median,
        "note": (
            "raw_value is the width-summed, background-subtracted profile; a fixed "
            "migration window is smile-affected. baseline=straight removes a local "
            "endpoint-to-endpoint baseline, clipped at zero"
        ),
    }

# =============================================================================
# Write outputs.
# =============================================================================

output_csv_path = gel_directory / (checks_stem + ".csv")
field_names = list(output_rows[0].keys())
with output_csv_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=field_names)
    writer.writeheader()
    writer.writerows(output_rows)
emit_message("output", "wrote " + str(output_csv_path))

value_summary = summarize_values(
    [
        (float(r["value"]) if r["value"] != "" else None)
        if isinstance(r["value"], str)
        else r["value"]
        for r in output_rows
    ]
)

checks_path = gel_directory / (checks_stem + "_checks.json")
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
                "profile_csv": str(profile_csv_path)
                if parsed_arguments.region is not None
                else None,
            },
            "selection": selection_description,
            "provenance": provenance,
            "value_summary_over_all_output_lanes": value_summary,
            "checks": check_records,
        },
        indent=2,
    )
    + "\n"
)
emit_message("output", "wrote " + str(checks_path))

# Console summary so the operator sees the shape of the result without opening a
# file: the per-lane values and the spread, which is what "is this consistent?"
# comes down to once a band or region is chosen.
emit_message("summary", "gel %s" % gel_id)
emit_message("summary", "selection: %s" % json.dumps(selection_description))
for row in output_rows:
    emit_message(
        "value",
        "lane %2s well %2s role %-16s %-14s value %s"
        % (
            row["lane_index"],
            row["well_number"],
            row["analysis_role"],
            row["sample_label"],
            row["value"],
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
for record in check_records:
    if record["severity"] == "soft" and not record["passed"]:
        emit_message("check", "WARNING " + record["check"] + ": " + record["detail"])
