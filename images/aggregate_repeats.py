# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib>=3.7",
# ]
# ///
r"""
Stack replicate gels into per-sample means, spread, and a bar chart. This is the
aggregation step: it reads a manifest of replicate gels, loads each gel's
single_experiment_<selector>.csv, joins replicates by the stable biological name
(sample_label), and reports n / mean / sd / cv of the blank-corrected value across
replicates. There is no reference scoring anywhere in this pipeline; if you want a
ratio, name a denominator explicitly with --normalize-to.

The manifest is validated up front by the same hard gates validate_manifest.py
applies (duplicated here so this script is runnable on its own and refuses a bad
manifest rather than aggregating garbage). All replicates must share ONE selector,
so a "none" baseline extract is never silently averaged with a "straight" one; a
mismatch is a hard failure.

Options:
    --normalize-to SAMPLE_LABEL
        Within each gel, divide every value_corrected by that sample's
        value_corrected before aggregating (its own ratio becomes 1.0). Off by
        default. Fails if the sample is absent or zero in any gel.
    --group-by COLUMN
        Aggregate by a categorical column present in the single_experiment CSVs
        (e.g. condition_type) instead of by sample_label.
    --selector SELECTOR
        Use single_experiment_<SELECTOR>.csv in each gel explicitly, when a gel
        directory holds more than one single_experiment CSV.

Addressing follows the rest of the pipeline: point at the manifest.csv or at the
directory that contains it. Outputs are written beside the manifest, one level up
from the gels, because they are the cross-gel product:
    aggregate_<selector>.csv          one row per sample_label or group
    aggregate_<selector>.png          mean with SD error bars (NOT viewed)
    aggregate_<selector>_checks.json  inputs, settings, per-group values, verdicts

Flat and procedural; emit_message / die / record_check are the only helpers,
duplicated verbatim per the house convention.
"""

import argparse
import csv
import datetime
import json
import pathlib
import sys

# Never render interactively; this script only writes a PNG file and must never
# open a display. Set the non-interactive backend before importing pyplot.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot

# Standard filenames resolved inside the manifest's directory and the gel dirs.
MANIFEST_FILENAME = "manifest.csv"
SINGLE_EXPERIMENT_PREFIX = "single_experiment_"
SINGLE_EXPERIMENT_GLOB = "single_experiment_*.csv"

# Required manifest columns (duplicated from validate_manifest.py so this script
# refuses a bad manifest on its own). Free condition columns beyond these are
# allowed; --group-by can name one that also lives in the single_experiment CSVs.
REQUIRED_MANIFEST_COLUMNS = ["gel_id", "replicate", "analysis_path"]

# Columns the single_experiment CSV must carry for aggregation.
REQUIRED_SINGLE_EXPERIMENT_COLUMNS = ["sample_label", "condition_type", "value_corrected"]

# Template example markers (duplicated from validate_manifest.py) so the same
# leftover example row is caught here rather than aggregated.
IS_EXAMPLE_COLUMN = "is_example"
EXAMPLE_GEL_ID = "EXAMPLE_replace_me"
IS_EXAMPLE_TRUE_TOKENS = {"true", "yes", "1"}

# A group whose condition_type is not constant across its members is labelled this
# rather than silently taking one member's value.
MIXED_CONDITION_LABEL = "mixed"

# Distinct condition_type values are coloured from this fixed palette (assigned in
# sorted order so the same condition gets the same colour run to run). ASCII hex.
CONDITION_COLOUR_PALETTE = [
    "#4C72B0",
    "#DD8452",
    "#55A868",
    "#C44E52",
    "#8172B3",
    "#937860",
    "#DA8BC3",
    "#8C8C8C",
]
FALLBACK_CONDITION_COLOUR = "#333333"

PLOT_FIGURE_WIDTH_INCHES = 11.0
PLOT_FIGURE_HEIGHT_INCHES = 6.0
PLOT_DOTS_PER_INCH = 150
ERROR_BAR_CAP_SIZE_POINTS = 4

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
# terminal. output_stem starts selector-independent ("aggregate") so an early
# manifest failure is still writable; it becomes "aggregate_<selector>" once the
# selector is known, which is what the successful run and any late failure use.
accumulated_check_records = []
output_stem = "aggregate"


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
        hard_failure_checks_path = manifest_directory / (output_stem + "_checks.json")
        hard_failure_checks_path.write_text(
            json.dumps(
                {
                    "schema_version": "aggregate_repeats_1",
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
    prog="aggregate_repeats.py",
    description="Aggregate single-experiment values across replicate gels.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "manifest_path",
    help="the manifest.csv, or the directory that contains it",
)
argument_parser.add_argument(
    "--normalize-to",
    dest="normalize_to_sample_label",
    default=None,
    help="divide each gel's value_corrected by this sample_label before aggregating",
)
argument_parser.add_argument(
    "--group-by",
    dest="group_by_column",
    default=None,
    help="aggregate by this single_experiment column instead of sample_label",
)
argument_parser.add_argument(
    "--selector",
    dest="selector_override",
    default=None,
    help="use single_experiment_<SELECTOR>.csv when a gel dir has more than one",
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

normalize_to_sample_label = parsed_arguments.normalize_to_sample_label
group_by_column = parsed_arguments.group_by_column

# =============================================================================
# Validate the manifest up front with the same hard gates as validate_manifest.py
# (duplicated deliberately; this script must refuse a bad manifest on its own).
# utf-8-sig strips an Excel BOM and every cell is trimmed.
# =============================================================================

record_check(
    "manifest.csv present", "hard", manifest_path.is_file(), str(manifest_path)
)
with manifest_path.open(newline="", encoding="utf-8-sig") as manifest_handle:
    manifest_reader = csv.DictReader(manifest_handle)
    manifest_column_names = manifest_reader.fieldnames or []
    manifest_rows = [
        {
            column_name: (cell_value.strip() if isinstance(cell_value, str) else cell_value)
            for column_name, cell_value in raw_row.items()
        }
        for raw_row in manifest_reader
    ]
record_check(
    "manifest has data rows", "hard", len(manifest_rows) > 0, "%d data row(s)" % len(manifest_rows)
)
missing_required_columns = [
    column_name
    for column_name in REQUIRED_MANIFEST_COLUMNS
    if column_name not in manifest_column_names
]
record_check(
    "required manifest columns present",
    "hard",
    len(missing_required_columns) == 0,
    "missing: " + (", ".join(missing_required_columns) if missing_required_columns else "none"),
)

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
    "example row(s) at " + (", ".join(str(n) for n in example_row_numbers) if example_row_numbers else "none"),
)

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
    "duplicated: " + (", ".join(duplicated_gel_ids) if duplicated_gel_ids else "none"),
)

rows_with_bad_replicate = []
for row_ordinal, manifest_row in enumerate(non_example_rows, start=1):
    replicate_cell = manifest_row.get("replicate", "")
    if not (isinstance(replicate_cell, str) and replicate_cell.isdigit() and int(replicate_cell) >= 1):
        rows_with_bad_replicate.append(
            "row %d gel_id=%s replicate=%r"
            % (row_ordinal, manifest_row.get("gel_id", ""), replicate_cell)
        )
record_check(
    "replicate is a positive integer",
    "hard",
    len(rows_with_bad_replicate) == 0,
    "bad: " + ("; ".join(rows_with_bad_replicate) if rows_with_bad_replicate else "none"),
)

# Resolve each gel's analysis directory relative to the manifest and gate existence.
resolved_gel_directory_by_gel_id = {}
rows_with_missing_path = []
for manifest_row in non_example_rows:
    gel_id_value = manifest_row.get("gel_id", "")
    analysis_path_cell = manifest_row.get("analysis_path", "")
    candidate_path = pathlib.Path(analysis_path_cell)
    resolved_analysis_path = (
        candidate_path if candidate_path.is_absolute() else manifest_directory / candidate_path
    )
    if not resolved_analysis_path.exists():
        rows_with_missing_path.append(
            "gel_id=%s -> %s" % (gel_id_value, str(resolved_analysis_path))
        )
        continue
    resolved_gel_directory_by_gel_id[gel_id_value] = (
        resolved_analysis_path
        if resolved_analysis_path.is_dir()
        else resolved_analysis_path.parent
    )
record_check(
    "every analysis_path exists",
    "hard",
    len(rows_with_missing_path) == 0,
    "missing: " + ("; ".join(rows_with_missing_path) if rows_with_missing_path else "none"),
)

# =============================================================================
# Resolve ONE selector shared by every gel. With --selector, each gel must have
# that named file; otherwise auto-find the single single_experiment CSV per gel and
# require every gel to agree, so incompatible extractions cannot be averaged.
# =============================================================================

selector_by_gel_id = {}
single_experiment_path_by_gel_id = {}
gels_missing_single_experiment = []
gels_with_ambiguous_single_experiment = []
for gel_id_value, gel_directory in resolved_gel_directory_by_gel_id.items():
    if parsed_arguments.selector_override is not None:
        candidate_single_experiment = gel_directory / (
            SINGLE_EXPERIMENT_PREFIX + parsed_arguments.selector_override + ".csv"
        )
        if candidate_single_experiment.is_file():
            selector_by_gel_id[gel_id_value] = parsed_arguments.selector_override
            single_experiment_path_by_gel_id[gel_id_value] = candidate_single_experiment
        else:
            gels_missing_single_experiment.append(
                "%s (%s)" % (gel_id_value, str(candidate_single_experiment))
            )
    else:
        found_single_experiments = sorted(gel_directory.glob(SINGLE_EXPERIMENT_GLOB))
        if len(found_single_experiments) == 0:
            gels_missing_single_experiment.append(
                "%s (none in %s)" % (gel_id_value, str(gel_directory))
            )
        elif len(found_single_experiments) > 1:
            gels_with_ambiguous_single_experiment.append(
                "%s: %s"
                % (
                    gel_id_value,
                    ", ".join(path.name for path in found_single_experiments),
                )
            )
        else:
            single_experiment_path = found_single_experiments[0]
            selector_by_gel_id[gel_id_value] = single_experiment_path.stem[
                len(SINGLE_EXPERIMENT_PREFIX):
            ]
            single_experiment_path_by_gel_id[gel_id_value] = single_experiment_path

record_check(
    "every gel has its single_experiment CSV",
    "hard",
    len(gels_missing_single_experiment) == 0,
    "missing: " + ("; ".join(gels_missing_single_experiment) if gels_missing_single_experiment else "none"),
)
record_check(
    "single_experiment CSV is unambiguous per gel",
    "hard",
    len(gels_with_ambiguous_single_experiment) == 0,
    "ambiguous (pass --selector): "
    + ("; ".join(gels_with_ambiguous_single_experiment) if gels_with_ambiguous_single_experiment else "none"),
)

distinct_selectors = sorted(set(selector_by_gel_id.values()))
record_check(
    "all gels share one selector",
    "hard",
    len(distinct_selectors) == 1,
    "selectors seen: " + ", ".join(distinct_selectors),
)
selector = distinct_selectors[0]
output_stem = "aggregate_" + selector

# =============================================================================
# Load each gel's single_experiment rows, keyed by gel so normalization is per gel.
# =============================================================================

gels_missing_required_columns = []
gels_with_unparseable_values = []
gels_missing_group_by_column = []
rows_by_gel_id = {}
for gel_id_value in resolved_gel_directory_by_gel_id:
    single_experiment_path = single_experiment_path_by_gel_id[gel_id_value]
    with single_experiment_path.open(newline="", encoding="utf-8-sig") as single_experiment_handle:
        single_experiment_reader = csv.DictReader(single_experiment_handle)
        single_experiment_column_names = single_experiment_reader.fieldnames or []
        single_experiment_rows = [
            {
                column_name: (cell_value.strip() if isinstance(cell_value, str) else cell_value)
                for column_name, cell_value in raw_row.items()
            }
            for raw_row in single_experiment_reader
        ]
    missing_here = [
        column_name
        for column_name in REQUIRED_SINGLE_EXPERIMENT_COLUMNS
        if column_name not in single_experiment_column_names
    ]
    if missing_here:
        gels_missing_required_columns.append(
            "%s missing %s" % (gel_id_value, ", ".join(missing_here))
        )
        continue
    if group_by_column is not None and group_by_column not in single_experiment_column_names:
        gels_missing_group_by_column.append(gel_id_value)
        continue
    parsed_rows = []
    for single_experiment_row in single_experiment_rows:
        value_corrected_cell = single_experiment_row.get("value_corrected", "")
        try:
            value_corrected = float(value_corrected_cell)
        except (TypeError, ValueError):
            gels_with_unparseable_values.append(
                "%s sample_label=%s value_corrected=%r"
                % (gel_id_value, single_experiment_row.get("sample_label", ""), value_corrected_cell)
            )
            value_corrected = None
        parsed_rows.append(
            {
                "sample_label": single_experiment_row.get("sample_label", ""),
                "condition_type": single_experiment_row.get("condition_type", ""),
                "value_corrected": value_corrected,
                "group_by_value": (
                    single_experiment_row.get(group_by_column, "")
                    if group_by_column is not None
                    else None
                ),
            }
        )
    rows_by_gel_id[gel_id_value] = parsed_rows

record_check(
    "every single_experiment CSV has the required columns",
    "hard",
    len(gels_missing_required_columns) == 0,
    "; ".join(gels_missing_required_columns) if gels_missing_required_columns else "none",
)
record_check(
    "group-by column present in every single_experiment CSV",
    "hard",
    len(gels_missing_group_by_column) == 0,
    (
        "column %r missing in: %s" % (group_by_column, ", ".join(gels_missing_group_by_column))
        if gels_missing_group_by_column
        else ("group-by=%s" % group_by_column if group_by_column is not None else "grouping by sample_label")
    ),
)
record_check(
    "every value_corrected parses as a number",
    "hard",
    len(gels_with_unparseable_values) == 0,
    "; ".join(gels_with_unparseable_values) if gels_with_unparseable_values else "none",
)

# =============================================================================
# Optional per-gel normalization: divide each row's value_corrected by the named
# sample's value_corrected within the same gel. The denominator must be present
# exactly once and non-zero in every gel, or the ratio is undefined.
# =============================================================================

if normalize_to_sample_label is not None:
    gels_missing_normalizer = []
    gels_with_zero_normalizer = []
    normalizer_by_gel_id = {}
    for gel_id_value, parsed_rows in rows_by_gel_id.items():
        matching_normalizer_rows = [
            parsed_row
            for parsed_row in parsed_rows
            if parsed_row["sample_label"] == normalize_to_sample_label
        ]
        if len(matching_normalizer_rows) == 0:
            gels_missing_normalizer.append(gel_id_value)
        elif matching_normalizer_rows[0]["value_corrected"] in (None, 0.0):
            gels_with_zero_normalizer.append(gel_id_value)
        else:
            normalizer_by_gel_id[gel_id_value] = matching_normalizer_rows[0]["value_corrected"]
    record_check(
        "normalizer sample present in every gel",
        "hard",
        len(gels_missing_normalizer) == 0,
        "sample %r missing in: %s"
        % (normalize_to_sample_label, ", ".join(gels_missing_normalizer) if gels_missing_normalizer else "none"),
    )
    record_check(
        "normalizer value is non-zero in every gel",
        "hard",
        len(gels_with_zero_normalizer) == 0,
        "sample %r zero/blank in: %s"
        % (normalize_to_sample_label, ", ".join(gels_with_zero_normalizer) if gels_with_zero_normalizer else "none"),
    )
    for gel_id_value, parsed_rows in rows_by_gel_id.items():
        normalizer_value = normalizer_by_gel_id[gel_id_value]
        for parsed_row in parsed_rows:
            if parsed_row["value_corrected"] is not None:
                parsed_row["value_corrected"] = parsed_row["value_corrected"] / normalizer_value

# =============================================================================
# Group across all gels and compute n / mean / sd / cv. The group key is the
# --group-by column value, or sample_label by default. condition_type is carried
# through, collapsing to a single value or MIXED if it is not constant in a group.
# =============================================================================

values_by_group_key = {}
condition_types_by_group_key = {}
for gel_id_value, parsed_rows in rows_by_gel_id.items():
    for parsed_row in parsed_rows:
        if parsed_row["value_corrected"] is None:
            continue
        group_key = (
            parsed_row["group_by_value"]
            if group_by_column is not None
            else parsed_row["sample_label"]
        )
        values_by_group_key.setdefault(group_key, []).append(parsed_row["value_corrected"])
        condition_types_by_group_key.setdefault(group_key, set()).add(
            parsed_row["condition_type"]
        )

group_key_column_name = "group" if group_by_column is not None else "sample_label"
aggregate_rows = []
single_replicate_group_keys = []
mixed_condition_group_keys = []
for group_key in sorted(values_by_group_key):
    group_values = values_by_group_key[group_key]
    replicate_count = len(group_values)
    mean_value = sum(group_values) / replicate_count
    if replicate_count > 1:
        sample_variance = sum(
            (value - mean_value) ** 2 for value in group_values
        ) / (replicate_count - 1)
        standard_deviation = sample_variance ** 0.5
        coefficient_of_variation = (
            standard_deviation / abs(mean_value) if mean_value != 0 else None
        )
    else:
        # A single replicate has no spread to report; leave sd/cv empty rather than
        # writing a misleading zero.
        single_replicate_group_keys.append(str(group_key))
        standard_deviation = None
        coefficient_of_variation = None
    condition_types_here = sorted(condition_types_by_group_key[group_key])
    if len(condition_types_here) == 1:
        condition_type_value = condition_types_here[0]
    else:
        condition_type_value = MIXED_CONDITION_LABEL
        mixed_condition_group_keys.append(str(group_key))
    aggregate_rows.append(
        {
            group_key_column_name: group_key,
            "condition_type": condition_type_value,
            "n": replicate_count,
            "mean": round(mean_value, 4),
            "sd": ("" if standard_deviation is None else round(standard_deviation, 4)),
            "cv": ("" if coefficient_of_variation is None else round(coefficient_of_variation, 4)),
        }
    )

record_check(
    "every group has at least two replicates",
    "soft",
    len(single_replicate_group_keys) == 0,
    "single-replicate group(s) (sd/cv left blank): "
    + (", ".join(single_replicate_group_keys) if single_replicate_group_keys else "none"),
)
record_check(
    "condition_type is consistent within every group",
    "soft",
    len(mixed_condition_group_keys) == 0,
    "mixed-condition group(s): "
    + (", ".join(mixed_condition_group_keys) if mixed_condition_group_keys else "none"),
)
record_check(
    "at least one group aggregated",
    "hard",
    len(aggregate_rows) > 0,
    "%d group(s)" % len(aggregate_rows),
)

# =============================================================================
# Write the aggregate CSV.
# =============================================================================

output_csv_path = manifest_directory / (output_stem + ".csv")
output_column_names = list(aggregate_rows[0].keys())
with output_csv_path.open("w", newline="", encoding="utf-8") as output_csv_handle:
    csv_writer = csv.DictWriter(output_csv_handle, fieldnames=output_column_names)
    csv_writer.writeheader()
    csv_writer.writerows(aggregate_rows)
emit_message("output", "wrote " + str(output_csv_path))

# =============================================================================
# Write the bar chart: mean with SD error bars, coloured by condition_type. The
# PNG is written and never opened by this script.
# =============================================================================

condition_types_in_order = sorted(
    {aggregate_row["condition_type"] for aggregate_row in aggregate_rows}
)
colour_by_condition_type = {}
for condition_ordinal, condition_type_value in enumerate(condition_types_in_order):
    if condition_ordinal < len(CONDITION_COLOUR_PALETTE):
        colour_by_condition_type[condition_type_value] = CONDITION_COLOUR_PALETTE[condition_ordinal]
    else:
        colour_by_condition_type[condition_type_value] = FALLBACK_CONDITION_COLOUR

bar_positions = list(range(len(aggregate_rows)))
bar_heights = [aggregate_row["mean"] for aggregate_row in aggregate_rows]
# A blank sd (single replicate) draws no error bar.
bar_error_bars = [
    (0.0 if aggregate_row["sd"] == "" else aggregate_row["sd"])
    for aggregate_row in aggregate_rows
]
bar_colours = [
    colour_by_condition_type[aggregate_row["condition_type"]] for aggregate_row in aggregate_rows
]
bar_tick_labels = [
    str(aggregate_row[group_key_column_name]) + "\n(n=" + str(aggregate_row["n"]) + ")"
    for aggregate_row in aggregate_rows
]

figure_object, axes_object = matplotlib.pyplot.subplots(
    figsize=(PLOT_FIGURE_WIDTH_INCHES, PLOT_FIGURE_HEIGHT_INCHES)
)
axes_object.bar(
    bar_positions,
    bar_heights,
    yerr=bar_error_bars,
    capsize=ERROR_BAR_CAP_SIZE_POINTS,
    color=bar_colours,
)
axes_object.set_xticks(bar_positions)
axes_object.set_xticklabels(bar_tick_labels, rotation=45, ha="right", fontsize=8)
normalized_axis_note = (
    " (normalized to " + normalize_to_sample_label + ")"
    if normalize_to_sample_label is not None
    else ""
)
axes_object.set_ylabel("mean value_corrected" + normalized_axis_note)
axes_object.set_xlabel(group_key_column_name)
axes_object.set_title("aggregate across replicates\n" + selector)
legend_handles = [
    matplotlib.pyplot.Rectangle(
        (0, 0), 1, 1, color=colour_by_condition_type[condition_type_value]
    )
    for condition_type_value in condition_types_in_order
]
axes_object.legend(
    legend_handles, condition_types_in_order, title="condition_type", fontsize=8
)
figure_object.tight_layout()
output_png_path = manifest_directory / (output_stem + ".png")
figure_object.savefig(output_png_path, dpi=PLOT_DOTS_PER_INCH)
matplotlib.pyplot.close(figure_object)
emit_message("output", "wrote " + str(output_png_path))

# =============================================================================
# Write the checks JSON.
# =============================================================================

checks_path = manifest_directory / (output_stem + "_checks.json")
checks_path.write_text(
    json.dumps(
        {
            "schema_version": "aggregate_repeats_1",
            "generated_at": datetime.datetime.now()
            .astimezone()
            .isoformat(timespec="seconds"),
            "manifest": str(manifest_path),
            "selector": selector,
            "settings": {
                "normalize_to": normalize_to_sample_label,
                "group_by": group_by_column,
                "group_key_column": group_key_column_name,
            },
            "inputs": {
                gel_id_value: str(single_experiment_path_by_gel_id[gel_id_value])
                for gel_id_value in sorted(single_experiment_path_by_gel_id)
            },
            "summary": {
                "gel_count": len(rows_by_gel_id),
                "group_count": len(aggregate_rows),
            },
            "outputs": {
                "aggregate_csv": str(output_csv_path),
                "aggregate_png": str(output_png_path),
            },
            "checks": accumulated_check_records,
        },
        indent=2,
    )
    + "\n"
)
emit_message("output", "wrote " + str(checks_path))

# Console summary to stderr (stdout stays empty so paths can be diffed).
emit_message(
    "summary",
    "selector %s across %d gel(s); %d group(s)%s"
    % (
        selector,
        len(rows_by_gel_id),
        len(aggregate_rows),
        "" if normalize_to_sample_label is None else " normalized to " + normalize_to_sample_label,
    ),
)
for aggregate_row in aggregate_rows:
    emit_message(
        "group",
        "%-16s %-16s n=%s mean=%s sd=%s cv=%s"
        % (
            str(aggregate_row[group_key_column_name]),
            aggregate_row["condition_type"],
            aggregate_row["n"],
            aggregate_row["mean"],
            aggregate_row["sd"],
            aggregate_row["cv"],
        ),
    )
for check_record in accumulated_check_records:
    if check_record["severity"] == "soft" and not check_record["passed"]:
        emit_message(
            "check", "WARNING " + check_record["check"] + ": " + check_record["detail"]
        )
