# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Validate a per-gel sample sheet against the lane profiles it describes, and report
the observed lane-to-well relationship. Hard-gates the structural invariants so a
malformed or self-contradictory sheet cannot reach analysis; warns on soft issues.

Reads only two files and no image data. Standard library only.

Column families (see the sample_sheet_template.csv the pipeline ships):
  identity        lane_index (image left-to-right, matches the profile CSV),
                  well_number (load order / sample-sheet order)
  disposition     lane_content   in {sample, ladder, empty}
                  analysis_role  in {reference, measured, excluded}
  experimental    condition_type in {input, experiment, positive_control,
                  negative_control, not_applicable}
  provenance      sample_label, prep_source, notes (and any factor columns, which
                  are permitted but not validated here)

The dangerous error this guards against is identity confusion: lane_index and
well_number are different integer spaces and may be inverted or permuted. This
tool never applies a flip. It reports the relationship it observes so the operator
sees it stated, and it hard-fails only on structural impossibility (a missing or
duplicated index), never on the mapping being non-identity.

Usage:
    uv run validate_sample_sheet.py '<sample_sheet.csv>' '<manual_lane_profiles.csv>'
"""

import argparse
import csv
import datetime
import json
import os
import pathlib
import sys

# =============================================================================
# Controlled vocabularies and required columns. Named so a reader sees the
# constraint, and so a value typo fails loudly instead of being accepted.
# =============================================================================

LANE_CONTENT_VOCABULARY = ("sample", "ladder", "empty")
ANALYSIS_ROLE_VOCABULARY = ("reference", "measured", "excluded")
CONDITION_TYPE_VOCABULARY = (
    "input",
    "experiment",
    "positive_control",
    "negative_control",
    "not_applicable",
)
# A ladder or empty lane holds no analyte, so it cannot be an experimental class
# and cannot be measured; these are the only legal dispositions for such lanes.
NON_ANALYTE_CONTENTS = ("ladder", "empty")
CONTROL_CONDITION_TYPES = ("positive_control", "negative_control")
REQUIRED_SHEET_COLUMNS = (
    "lane_index",
    "well_number",
    "lane_content",
    "analysis_role",
    "condition_type",
    "sample_label",
    "prep_source",
)
CHECKS_JSON_FILENAME = "sample_sheet_checks.json"

# =============================================================================
# The two permitted helpers, matching the other scripts' logging contract.
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
# Arguments
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="validate_sample_sheet.py",
    description="Validate a per-gel sample sheet against its lane profiles.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "sample_sheet_csv_path", help="Path to the sample sheet CSV"
)
argument_parser.add_argument(
    "profile_csv_path", help="Path to manual_lane_profiles.csv"
)
parsed_arguments = argument_parser.parse_args()

sample_sheet_csv_path = pathlib.Path(
    os.path.abspath(parsed_arguments.sample_sheet_csv_path)
)
profile_csv_path = pathlib.Path(os.path.abspath(parsed_arguments.profile_csv_path))
if not sample_sheet_csv_path.is_file():
    die("input", "sample sheet not found: " + str(sample_sheet_csv_path))
if not profile_csv_path.is_file():
    die("input", "profile CSV not found: " + str(profile_csv_path))
# The checks report lands next to the sample sheet, which lives in the gel_analysis
# directory, so it travels with the analysis it validates.
output_directory_path = sample_sheet_csv_path.parent

# =============================================================================
# Read the profile CSV for the cross-check: how many lanes were exported, and
# which lane_index values. utf-8-sig so an Excel BOM does not corrupt the header.
# =============================================================================

with open(profile_csv_path, newline="", encoding="utf-8-sig") as profile_file_handle:
    profile_reader = csv.DictReader(profile_file_handle)
    profile_field_names = profile_reader.fieldnames or []
    profile_rows = list(profile_reader)
if (
    "lane_index" not in profile_field_names
    or "comb_well_count" not in profile_field_names
):
    die("profile", "profile CSV is missing lane_index or comb_well_count; re-export.")
profile_lane_indices = sorted({int(row["lane_index"]) for row in profile_rows})
profile_comb_well_count = int(profile_rows[0]["comb_well_count"])
emit_message(
    "profile",
    "profile has %d lanes, comb_well_count=%d"
    % (len(profile_lane_indices), profile_comb_well_count),
)

# =============================================================================
# Read the sample sheet. utf-8-sig for the BOM; trim every cell, because a hand
# authored sheet exported from Excel carries stray whitespace that would break an
# exact-match join or a vocabulary check invisibly.
# =============================================================================

with open(sample_sheet_csv_path, newline="", encoding="utf-8-sig") as sheet_file_handle:
    sheet_reader = csv.DictReader(sheet_file_handle)
    sheet_field_names = sheet_reader.fieldnames or []
    sheet_rows_raw = list(sheet_reader)
sheet_rows = []
for raw_row in sheet_rows_raw:
    trimmed_row = {}
    for column_name, cell_value in raw_row.items():
        trimmed_row[column_name] = (cell_value or "").strip()
    sheet_rows.append(trimmed_row)

missing_columns = [
    name for name in REQUIRED_SHEET_COLUMNS if name not in sheet_field_names
]
if missing_columns:
    die(
        "schema",
        "sample sheet is missing required columns: " + ", ".join(missing_columns),
    )
if len(sheet_rows) == 0:
    die("schema", "sample sheet has a header but no rows.")

emit_message(
    "read",
    "sample sheet has %d rows and columns: %s"
    % (len(sheet_rows), ", ".join(sheet_field_names)),
)

# =============================================================================
# Checks. Each appends a verdict; the JSON is written even on hard failure so a
# failure is diagnosable from the file, then the process exits non-zero.
# =============================================================================

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


# --- row count vs comb ---
record_check(
    "row_count_matches_comb",
    "hard",
    len(sheet_rows) == profile_comb_well_count,
    "sheet rows=%d, profile comb_well_count=%d"
    % (len(sheet_rows), profile_comb_well_count),
)

# --- lane_index: integers, 1..comb, each once ---
sheet_lane_indices = []
lane_index_parse_failures = []
for row_number, row in enumerate(sheet_rows, start=1):
    try:
        sheet_lane_indices.append(int(row["lane_index"]))
    except (ValueError, KeyError):
        lane_index_parse_failures.append(
            "row %d: %r" % (row_number, row.get("lane_index"))
        )
record_check(
    "lane_index_integers",
    "hard",
    len(lane_index_parse_failures) == 0,
    "non-integer lane_index at " + "; ".join(lane_index_parse_failures)
    if lane_index_parse_failures
    else "all integer",
)
expected_index_set = set(range(1, profile_comb_well_count + 1))
record_check(
    "lane_index_is_full_set",
    "hard",
    sorted(sheet_lane_indices) == sorted(expected_index_set),
    "lane_index set = %s, expected 1..%d each once"
    % (sorted(sheet_lane_indices), profile_comb_well_count),
)

# --- well_number: integers, permutation of 1..comb, each once ---
sheet_well_numbers = []
well_number_parse_failures = []
for row_number, row in enumerate(sheet_rows, start=1):
    try:
        sheet_well_numbers.append(int(row["well_number"]))
    except (ValueError, KeyError):
        well_number_parse_failures.append(
            "row %d: %r" % (row_number, row.get("well_number"))
        )
record_check(
    "well_number_integers",
    "hard",
    len(well_number_parse_failures) == 0,
    "non-integer well_number at " + "; ".join(well_number_parse_failures)
    if well_number_parse_failures
    else "all integer",
)
record_check(
    "well_number_is_permutation",
    "hard",
    sorted(sheet_well_numbers) == sorted(expected_index_set),
    "well_number set = %s, expected a permutation of 1..%d"
    % (sorted(sheet_well_numbers), profile_comb_well_count),
)

# --- lane_index set matches the profile's lanes ---
record_check(
    "lane_index_matches_profile",
    "hard",
    sorted(sheet_lane_indices) == profile_lane_indices,
    "sheet lanes %s vs profile lanes %s"
    % (sorted(sheet_lane_indices), profile_lane_indices),
)

# --- controlled vocabularies, no blanks ---
for column_name, vocabulary in (
    ("lane_content", LANE_CONTENT_VOCABULARY),
    ("analysis_role", ANALYSIS_ROLE_VOCABULARY),
    ("condition_type", CONDITION_TYPE_VOCABULARY),
):
    offending = []
    for row in sheet_rows:
        value = row.get(column_name, "")
        if value == "" or value not in vocabulary:
            offending.append("lane %s: %r" % (row.get("lane_index"), value))
    record_check(
        column_name + "_in_vocabulary",
        "hard",
        len(offending) == 0,
        "allowed %s; offending: %s"
        % (list(vocabulary), "; ".join(offending) if offending else "none"),
    )

# --- cross-field legality: non-analyte lanes vs sample lanes ---
non_analyte_violations = []
sample_disposition_violations = []
for row in sheet_rows:
    lane_content_value = row.get("lane_content", "")
    condition_type_value = row.get("condition_type", "")
    analysis_role_value = row.get("analysis_role", "")
    lane_label = "lane %s" % row.get("lane_index")
    if lane_content_value in NON_ANALYTE_CONTENTS:
        if (
            condition_type_value != "not_applicable"
            or analysis_role_value != "excluded"
        ):
            non_analyte_violations.append(
                "%s (%s): condition_type=%s, analysis_role=%s"
                % (
                    lane_label,
                    lane_content_value,
                    condition_type_value,
                    analysis_role_value,
                )
            )
    elif lane_content_value == "sample":
        if condition_type_value == "not_applicable" or analysis_role_value not in (
            "reference",
            "measured",
        ):
            sample_disposition_violations.append(
                "%s: condition_type=%s, analysis_role=%s"
                % (lane_label, condition_type_value, analysis_role_value)
            )
record_check(
    "non_analyte_lanes_are_na_and_excluded",
    "hard",
    len(non_analyte_violations) == 0,
    "ladder/empty lanes must be condition_type=not_applicable and analysis_role=excluded; "
    + ("; ".join(non_analyte_violations) if non_analyte_violations else "ok"),
)
record_check(
    "sample_lanes_have_real_disposition",
    "hard",
    len(sample_disposition_violations) == 0,
    "sample lanes need a real condition_type and analysis_role in {reference, measured}; "
    + (
        "; ".join(sample_disposition_violations)
        if sample_disposition_violations
        else "ok"
    ),
)

# --- exactly one reference, and it is a sample ---
reference_rows = [
    row for row in sheet_rows if row.get("analysis_role", "") == "reference"
]
record_check(
    "exactly_one_reference",
    "hard",
    len(reference_rows) == 1,
    "found %d rows with analysis_role=reference (need exactly 1)" % len(reference_rows),
)
reference_is_sample = (
    len(reference_rows) == 1 and reference_rows[0].get("lane_content", "") == "sample"
)
record_check(
    "reference_is_a_sample",
    "hard",
    reference_is_sample,
    "the reference lane must have lane_content=sample"
    + (
        ""
        if reference_is_sample
        else "; reference lane_content=%s"
        % (reference_rows[0].get("lane_content") if reference_rows else "none")
    ),
)

# --- sample rows must carry a label ---
unlabelled_samples = [
    "lane %s" % row.get("lane_index")
    for row in sheet_rows
    if row.get("lane_content", "") == "sample" and row.get("sample_label", "") == ""
]
record_check(
    "sample_rows_have_label",
    "hard",
    len(unlabelled_samples) == 0,
    "sample lanes missing sample_label: "
    + ("; ".join(unlabelled_samples) if unlabelled_samples else "none"),
)

# --- soft: controls present, prep_source present on samples ---
present_condition_types = {row.get("condition_type", "") for row in sheet_rows}
record_check(
    "positive_control_present",
    "soft",
    "positive_control" in present_condition_types,
    "no lane marked condition_type=positive_control",
)
record_check(
    "negative_control_present",
    "soft",
    "negative_control" in present_condition_types,
    "no lane marked condition_type=negative_control",
)
samples_missing_prep_source = [
    "lane %s" % row.get("lane_index")
    for row in sheet_rows
    if row.get("lane_content", "") == "sample" and row.get("prep_source", "") == ""
]
record_check(
    "samples_have_prep_source",
    "soft",
    len(samples_missing_prep_source) == 0,
    "sample lanes missing prep_source: "
    + (
        "; ".join(samples_missing_prep_source)
        if samples_missing_prep_source
        else "none"
    ),
)

# =============================================================================
# Identity relationship: report, never enforce. Only meaningful if both index
# columns are well-formed sets, so guard on that.
# =============================================================================

identity_relationship = "undetermined"
lane_to_well_mapping = {}
if sorted(sheet_lane_indices) == sorted(expected_index_set) and sorted(
    sheet_well_numbers
) == sorted(expected_index_set):
    for row in sheet_rows:
        lane_to_well_mapping[int(row["lane_index"])] = int(row["well_number"])
    is_identity = all(lane_to_well_mapping[i] == i for i in expected_index_set)
    is_full_inversion = all(
        lane_to_well_mapping[i] == profile_comb_well_count + 1 - i
        for i in expected_index_set
    )
    if is_identity:
        identity_relationship = "identity (well_number == lane_index)"
    elif is_full_inversion:
        identity_relationship = "full_inversion (well_number == comb+1-lane_index)"
    else:
        identity_relationship = "permutation (neither identity nor a clean inversion)"
emit_message("identity", "lane-to-well relationship: " + identity_relationship)

# =============================================================================
# Report and exit
# =============================================================================

hard_failures = [
    c for c in check_records if c["severity"] == "hard" and not c["passed"]
]
soft_warnings = [
    c for c in check_records if c["severity"] == "soft" and not c["passed"]
]

checks_report = {
    "schema_version": "sample_sheet_checks_1",
    "generated_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "sample_sheet_csv": str(sample_sheet_csv_path),
    "profile_csv": str(profile_csv_path),
    "comb_well_count": profile_comb_well_count,
    "identity_relationship": identity_relationship,
    "lane_to_well_mapping": {
        str(k): v for k, v in sorted(lane_to_well_mapping.items())
    },
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
        "%d hard failure(s); sample sheet is not usable. See %s."
        % (len(hard_failures), str(checks_json_path)),
    )

emit_message(
    "done",
    "sample sheet valid; relationship is %s; %d warning(s)."
    % (identity_relationship, len(soft_warnings)),
)
sys.exit(0)
