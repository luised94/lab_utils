# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Hard-gate the replicate manifest before anything aggregates over it. The manifest
lives one level up from the gel analysis directories and lists, one row per gel,
which analyses are replicates of the same experiment. This validator refuses a
manifest that would silently corrupt an aggregate: a leftover template example row,
a duplicated gel_id, a non-integer replicate, or an analysis_path that does not
exist.

It also writes manifest_template.csv beside the manifest: the header plus one row
clearly marked is_example=true, so the "example row removed" gate always has
something concrete to catch and so an author has a correct shape to copy.

Hard gates (any failure stops the run after the checks file is written):
    manifest.csv present and has data rows
    required columns present: gel_id, replicate, analysis_path
    the template example row has been removed (no is_example=true row remains)
    gel_id is unique across rows
    replicate is a positive integer on every row
    every analysis_path exists, resolved relative to the manifest's directory

Soft gate (warns, does not stop):
    a referenced gel directory has no manual_lane_profiles.csv export yet

Addressing follows the rest of the pipeline: point at the manifest.csv or at the
directory that contains it. Stdlib only. Flat and procedural; emit_message / die /
record_check are the only helpers, duplicated verbatim per the house convention.
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys

# Standard filenames resolved inside the manifest's directory.
MANIFEST_FILENAME = "manifest.csv"
MANIFEST_TEMPLATE_FILENAME = "manifest_template.csv"

# The export the analyzer/extractor consume; a referenced gel that lacks it has not
# been through the ImageJ export yet, which is a warning rather than a failure.
PROFILE_EXPORT_FILENAME = "manual_lane_profiles.csv"

# Required manifest columns. Free condition columns beyond these are allowed and
# ignored here; the aggregator is what may group on them.
REQUIRED_MANIFEST_COLUMNS = ["gel_id", "replicate", "analysis_path"]

# The template example row carries this marker column set true, and this reserved
# gel_id, so the "example row removed" gate catches it whether the author cleared
# the column or only the value.
IS_EXAMPLE_COLUMN = "is_example"
EXAMPLE_GEL_ID = "EXAMPLE_replace_me"
IS_EXAMPLE_TRUE_TOKENS = {"true", "yes", "1"}

# The single example row written into the template. analysis_path is a relative
# placeholder the author overwrites with a real gel analysis directory.
MANIFEST_TEMPLATE_ROWS = [
    {
        "is_example": "true",
        "gel_id": EXAMPLE_GEL_ID,
        "replicate": "1",
        "analysis_path": "./replace_with_gel_analysis_directory",
    }
]

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
# terminal. manifest_directory is set before any check runs.
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
        hard_failure_checks_path = manifest_directory / "manifest_checks.json"
        hard_failure_checks_path.write_text(
            json.dumps(
                {
                    "schema_version": "validate_manifest_1",
                    "generated_at": datetime.datetime.now()
                    .astimezone()
                    .isoformat(timespec="seconds"),
                    "manifest": str(manifest_path),
                    "checks": accumulated_check_records,
                    "failed": {"check": check_name, "detail": detail_text},
                },
                indent=2,
            )
            + "\n"
        )
        die("checks", "hard check failed: " + check_name + " (" + detail_text + ")")


# =============================================================================
# Arguments and path resolution
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="validate_manifest.py",
    description="Hard-gate the replicate manifest and write a template beside it.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "manifest_path",
    help="the manifest.csv, or the directory that contains it",
)
parsed_arguments = argument_parser.parse_args()

given_path = pathlib.Path(parsed_arguments.manifest_path)
if given_path.is_file():
    manifest_path = given_path
    manifest_directory = given_path.parent
elif given_path.is_dir():
    manifest_directory = given_path
    manifest_path = given_path / MANIFEST_FILENAME
else:
    # manifest_directory is referenced by record_check on early failure, so a path
    # that does not exist at all is fatal here, plainly.
    print("[input] ERROR: path does not exist: " + str(given_path), file=sys.stderr)
    sys.exit(1)

# =============================================================================
# Always (re)write the template beside the manifest, so the scaffold with a
# catchable example row exists even on a first run where manifest.csv is absent.
# =============================================================================

manifest_template_path = manifest_directory / MANIFEST_TEMPLATE_FILENAME
template_column_names = [IS_EXAMPLE_COLUMN] + REQUIRED_MANIFEST_COLUMNS
with manifest_template_path.open("w", newline="", encoding="utf-8") as template_handle:
    template_writer = csv.DictWriter(template_handle, fieldnames=template_column_names)
    template_writer.writeheader()
    template_writer.writerows(MANIFEST_TEMPLATE_ROWS)
emit_message("output", "wrote " + str(manifest_template_path))

# =============================================================================
# Read the manifest. utf-8-sig strips an Excel BOM and every cell is trimmed,
# because a hand-exported CSV carries a BOM and stray whitespace that would break
# the exact uniqueness and vocabulary checks invisibly.
# =============================================================================

record_check(
    "manifest.csv present",
    "hard",
    manifest_path.is_file(),
    str(manifest_path),
)
with manifest_path.open(newline="", encoding="utf-8-sig") as manifest_handle:
    manifest_reader = csv.DictReader(manifest_handle)
    manifest_column_names = manifest_reader.fieldnames or []
    manifest_rows = [
        {
            column_name: (
                cell_value.strip() if isinstance(cell_value, str) else cell_value
            )
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in manifest_reader
    ]
record_check(
    "manifest has data rows",
    "hard",
    len(manifest_rows) > 0,
    "%d data row(s)" % len(manifest_rows),
)

missing_required_columns = [
    column_name
    for column_name in REQUIRED_MANIFEST_COLUMNS
    if column_name not in manifest_column_names
]
record_check(
    "required columns present",
    "hard",
    len(missing_required_columns) == 0,
    "missing: "
    + (", ".join(missing_required_columns) if missing_required_columns else "none")
    + "; present: "
    + ", ".join(manifest_column_names),
)

# =============================================================================
# Example row removed. A row counts as the template example if it is marked
# is_example=true OR still carries the reserved example gel_id, so the gate catches
# a half-cleared example either way. All later checks run over the non-example rows.
# =============================================================================

has_is_example_column = IS_EXAMPLE_COLUMN in manifest_column_names
example_row_numbers = []
non_example_rows = []
for row_ordinal, manifest_row in enumerate(manifest_rows, start=1):
    is_example_cell = manifest_row.get(IS_EXAMPLE_COLUMN, "")
    is_marked_example = (
        isinstance(is_example_cell, str)
        and is_example_cell.strip().lower() in IS_EXAMPLE_TRUE_TOKENS
    ) or manifest_row.get("gel_id", "") == EXAMPLE_GEL_ID
    if is_marked_example:
        example_row_numbers.append(row_ordinal)
    else:
        non_example_rows.append(manifest_row)
record_check(
    "template example row removed",
    "hard",
    len(example_row_numbers) == 0,
    (
        "example row(s) still present at data row(s) "
        + ", ".join(str(row_number) for row_number in example_row_numbers)
        if example_row_numbers
        else (
            "no example row"
            + (
                ""
                if has_is_example_column
                else " (no is_example column; only the reserved gel_id was checkable)"
            )
        )
    ),
)

# =============================================================================
# gel_id unique across the real (non-example) rows.
# =============================================================================

gel_id_counts = {}
for manifest_row in non_example_rows:
    gel_id_value = manifest_row.get("gel_id", "")
    gel_id_counts[gel_id_value] = gel_id_counts.get(gel_id_value, 0) + 1
duplicated_gel_ids = sorted(
    gel_id_value for gel_id_value, count in gel_id_counts.items() if count > 1
)
record_check(
    "gel_id is unique",
    "hard",
    len(duplicated_gel_ids) == 0,
    "duplicated gel_id(s): "
    + (", ".join(duplicated_gel_ids) if duplicated_gel_ids else "none"),
)

# =============================================================================
# replicate is a positive integer on every real row.
# =============================================================================

rows_with_bad_replicate = []
for row_ordinal, manifest_row in enumerate(non_example_rows, start=1):
    replicate_cell = manifest_row.get("replicate", "")
    replicate_is_positive_integer = (
        isinstance(replicate_cell, str)
        and replicate_cell.isdigit()
        and int(replicate_cell) >= 1
    )
    if not replicate_is_positive_integer:
        rows_with_bad_replicate.append(
            "row %d gel_id=%s replicate=%r"
            % (row_ordinal, manifest_row.get("gel_id", ""), replicate_cell)
        )
record_check(
    "replicate is a positive integer",
    "hard",
    len(rows_with_bad_replicate) == 0,
    "bad replicate(s): "
    + ("; ".join(rows_with_bad_replicate) if rows_with_bad_replicate else "none"),
)

# =============================================================================
# every analysis_path exists, resolved relative to the manifest's directory.
# =============================================================================

rows_with_missing_path = []
resolved_analysis_path_by_gel_id = {}
for manifest_row in non_example_rows:
    analysis_path_cell = manifest_row.get("analysis_path", "")
    candidate_path = pathlib.Path(analysis_path_cell)
    resolved_analysis_path = (
        candidate_path
        if candidate_path.is_absolute()
        else manifest_directory / candidate_path
    )
    resolved_analysis_path_by_gel_id[manifest_row.get("gel_id", "")] = (
        resolved_analysis_path
    )
    if not resolved_analysis_path.exists():
        rows_with_missing_path.append(
            "gel_id=%s analysis_path=%s -> %s"
            % (
                manifest_row.get("gel_id", ""),
                analysis_path_cell,
                str(resolved_analysis_path),
            )
        )
record_check(
    "every analysis_path exists",
    "hard",
    len(rows_with_missing_path) == 0,
    "missing path(s): "
    + ("; ".join(rows_with_missing_path) if rows_with_missing_path else "none"),
)

# =============================================================================
# Soft: a referenced gel directory has no profile export yet. Resolve the gel
# directory as the analysis_path itself when it is a directory, or its parent when
# it points at a file inside the directory.
# =============================================================================

gel_ids_without_profile_export = []
for gel_id_value, resolved_analysis_path in resolved_analysis_path_by_gel_id.items():
    gel_directory = (
        resolved_analysis_path
        if resolved_analysis_path.is_dir()
        else resolved_analysis_path.parent
    )
    if not (gel_directory / PROFILE_EXPORT_FILENAME).is_file():
        gel_ids_without_profile_export.append(gel_id_value)
record_check(
    "every referenced gel has a profile export",
    "soft",
    len(gel_ids_without_profile_export) == 0,
    "gel(s) without %s: %s"
    % (
        PROFILE_EXPORT_FILENAME,
        ", ".join(gel_ids_without_profile_export)
        if gel_ids_without_profile_export
        else "none",
    ),
)

# =============================================================================
# Write the checks JSON.
# =============================================================================

checks_path = manifest_directory / "manifest_checks.json"
checks_path.write_text(
    json.dumps(
        {
            "schema_version": "validate_manifest_1",
            "generated_at": datetime.datetime.now()
            .astimezone()
            .isoformat(timespec="seconds"),
            "manifest": str(manifest_path),
            "template": str(manifest_template_path),
            "summary": {
                "data_row_count": len(manifest_rows),
                "non_example_row_count": len(non_example_rows),
                "gel_ids": [
                    manifest_row.get("gel_id", "") for manifest_row in non_example_rows
                ],
                "extra_condition_columns": [
                    column_name
                    for column_name in manifest_column_names
                    if column_name not in REQUIRED_MANIFEST_COLUMNS
                    and column_name != IS_EXAMPLE_COLUMN
                ],
            },
            "checks": accumulated_check_records,
        },
        indent=2,
    )
    + "\n"
)
emit_message("output", "wrote " + str(checks_path))

emit_message(
    "summary",
    "manifest %s: %d row(s), %d after dropping example row(s)"
    % (str(manifest_path), len(manifest_rows), len(non_example_rows)),
)
for check_record in accumulated_check_records:
    if check_record["severity"] == "soft" and not check_record["passed"]:
        emit_message(
            "check", "WARNING " + check_record["check"] + ": " + check_record["detail"]
        )
emit_message("summary", "all hard gates passed")
