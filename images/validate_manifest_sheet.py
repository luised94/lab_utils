# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Validate a replicate manifest: the file that ties replicate gels together for the
aggregation stage. Hard-gates structure and that every referenced analysis path
exists; warns when a referenced gel directory has not been analysed yet. Standard
library only; reads no image data.

Columns (see manifest_template.csv):
  experiment_id   groups replicate gels of the same experiment
  analysis_path   the gel_analysis directory, or any file inside it; may be
                  relative to the manifest's own location
  gel_id          stable id for one gel; unique across the manifest
  replicate       positive integer replicate number within the experiment
  notes           free text

The manifest lives in the experiment directory alongside the TIFFs, so relative
analysis_path values resolve against the manifest's directory.

Usage:
    uv run validate_manifest.py '<manifest.csv>'
"""

import argparse
import csv
import datetime
import json
import os
import pathlib
import sys

REQUIRED_MANIFEST_COLUMNS = (
    "experiment_id",
    "analysis_path",
    "gel_id",
    "replicate",
    "notes",
)
EXAMPLE_ROW_EXPERIMENT_ID = "example-remove-this-row"
# A finished gel directory contains at least the profile export; the measurements
# CSV appears only after the measure stage, so its absence is a warning, not a fail.
EXPECTED_PROFILE_FILENAME = "manual_lane_profiles.csv"
CHECKS_JSON_FILENAME = "manifest_checks.json"

ACCUMULATED_RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    ACCUMULATED_RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


argument_parser = argparse.ArgumentParser(
    prog="validate_manifest.py",
    description="Validate a replicate manifest and its referenced analysis paths.",
    allow_abbrev=False,
)
argument_parser.add_argument("manifest_csv_path", help="Path to the manifest CSV")
parsed_arguments = argument_parser.parse_args()

manifest_csv_path = pathlib.Path(os.path.abspath(parsed_arguments.manifest_csv_path))
if not manifest_csv_path.is_file():
    die("input", "manifest not found: " + str(manifest_csv_path))
manifest_directory = manifest_csv_path.parent
output_directory_path = manifest_directory

with open(manifest_csv_path, newline="", encoding="utf-8-sig") as manifest_file_handle:
    manifest_reader = csv.DictReader(manifest_file_handle)
    manifest_field_names = manifest_reader.fieldnames or []
    manifest_rows_raw = list(manifest_reader)
manifest_rows = []
for raw_row in manifest_rows_raw:
    trimmed_row = {}
    for column_name, cell_value in raw_row.items():
        trimmed_row[column_name] = (cell_value or "").strip()
    manifest_rows.append(trimmed_row)

missing_columns = [
    name for name in REQUIRED_MANIFEST_COLUMNS if name not in manifest_field_names
]
if missing_columns:
    die("schema", "manifest is missing required columns: " + ", ".join(missing_columns))

check_records = []


def record_check(check_name, severity, passed, detail_text):
    check_records.append(
        {
            "check_name": check_name,
            "severity": severity,
            "passed": bool(passed),
            "detail": detail_text,
        }
    )


# --- the shipped example row must be gone ---
example_rows_present = [
    row
    for row in manifest_rows
    if row.get("experiment_id", "") == EXAMPLE_ROW_EXPERIMENT_ID
]
record_check(
    "example_row_removed",
    "hard",
    len(example_rows_present) == 0,
    "the template's example row is still present; remove the row with experiment_id="
    + EXAMPLE_ROW_EXPERIMENT_ID,
)

# --- at least one real data row ---
real_rows = [
    row
    for row in manifest_rows
    if row.get("experiment_id", "") != EXAMPLE_ROW_EXPERIMENT_ID
]
record_check(
    "has_data_rows",
    "hard",
    len(real_rows) > 0,
    "manifest has %d data row(s) after excluding the example" % len(real_rows),
)

# --- experiment_id and gel_id non-blank; gel_id unique ---
blank_experiment_ids = [
    str(i + 1) for i, row in enumerate(real_rows) if row.get("experiment_id", "") == ""
]
record_check(
    "experiment_id_present",
    "hard",
    len(blank_experiment_ids) == 0,
    "rows with blank experiment_id: "
    + ("; ".join(blank_experiment_ids) if blank_experiment_ids else "none"),
)
gel_ids = [row.get("gel_id", "") for row in real_rows]
blank_gel_ids = [str(i + 1) for i, gid in enumerate(gel_ids) if gid == ""]
record_check(
    "gel_id_present",
    "hard",
    len(blank_gel_ids) == 0,
    "rows with blank gel_id: "
    + ("; ".join(blank_gel_ids) if blank_gel_ids else "none"),
)
duplicate_gel_ids = sorted(
    {gid for gid in gel_ids if gid != "" and gel_ids.count(gid) > 1}
)
record_check(
    "gel_id_unique",
    "hard",
    len(duplicate_gel_ids) == 0,
    "gel_id must be unique; duplicates: "
    + ("; ".join(duplicate_gel_ids) if duplicate_gel_ids else "none"),
)

# --- replicate is a positive integer ---
replicate_failures = []
for row in real_rows:
    replicate_text = row.get("replicate", "")
    try:
        replicate_value = int(replicate_text)
        if replicate_value < 1:
            replicate_failures.append(
                "gel_id %s: %r" % (row.get("gel_id"), replicate_text)
            )
    except ValueError:
        replicate_failures.append("gel_id %s: %r" % (row.get("gel_id"), replicate_text))
record_check(
    "replicate_positive_integer",
    "hard",
    len(replicate_failures) == 0,
    "non-positive-integer replicate at "
    + ("; ".join(replicate_failures) if replicate_failures else "none"),
)

# --- every analysis_path exists; resolved relative to the manifest's directory ---
missing_paths = []
incomplete_gel_directories = []
for row in real_rows:
    analysis_path_text = row.get("analysis_path", "")
    if analysis_path_text == "":
        missing_paths.append("gel_id %s: blank analysis_path" % row.get("gel_id"))
        continue
    candidate_path = pathlib.Path(analysis_path_text)
    if not candidate_path.is_absolute():
        candidate_path = manifest_directory / candidate_path
    if not candidate_path.exists():
        missing_paths.append("gel_id %s: %s" % (row.get("gel_id"), str(candidate_path)))
        continue
    # Resolve to the gel directory: the path may be the dir itself or a file in it.
    gel_directory = candidate_path if candidate_path.is_dir() else candidate_path.parent
    if not (gel_directory / EXPECTED_PROFILE_FILENAME).is_file():
        incomplete_gel_directories.append(
            "gel_id %s: %s (no %s)"
            % (row.get("gel_id"), str(gel_directory), EXPECTED_PROFILE_FILENAME)
        )
record_check(
    "analysis_paths_exist",
    "hard",
    len(missing_paths) == 0,
    "missing analysis_path(s): "
    + ("; ".join(missing_paths) if missing_paths else "none"),
)
record_check(
    "referenced_gels_have_profiles",
    "soft",
    len(incomplete_gel_directories) == 0,
    "gel dir(s) with no profile export yet: "
    + ("; ".join(incomplete_gel_directories) if incomplete_gel_directories else "none"),
)

# --- report replicate grouping (information, not a gate) ---
replicates_by_experiment = {}
for row in real_rows:
    replicates_by_experiment.setdefault(row.get("experiment_id", ""), []).append(
        row.get("replicate", "")
    )

hard_failures = [
    c for c in check_records if c["severity"] == "hard" and not c["passed"]
]
soft_warnings = [
    c for c in check_records if c["severity"] == "soft" and not c["passed"]
]

checks_report = {
    "schema_version": "manifest_checks_1",
    "generated_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "manifest_csv": str(manifest_csv_path),
    "data_row_count": len(real_rows),
    "replicates_by_experiment": replicates_by_experiment,
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
        "%d hard failure(s); manifest is not usable. See %s."
        % (len(hard_failures), str(checks_json_path)),
    )

emit_message(
    "done",
    "manifest valid; %d gel(s) across %d experiment(s); %d warning(s)."
    % (len(real_rows), len(replicates_by_experiment), len(soft_warnings)),
)
sys.exit(0)
