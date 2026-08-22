# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
One-shot authoring tool for the n>1 aggregate fixture. NOT run by the harness and
NOT run at test time: it is the recorded procedure that PRODUCED the committed
second-replicate gel directory under tests/fixtures/, kept so the perturbation rule
is the source of truth rather than a pile of opaque numbers, and so the fixture can
be regenerated or re-perturbed by hand if it ever needs to change.

What it does: read the first (real) gel's manual_lane_profiles.csv, multiply every
raw_value by a fixed PER-LANE factor (see PER_LANE_RAW_VALUE_FACTOR below), and
write a second gel directory whose name -> gel_id differs from the first. The
sample_sheet.csv is copied verbatim, because the join that aggregate_repeats.py
performs is BY sample_label across gels: the two gels must share the same labels
for BK_blue (and every other sample) to land in the same group with n=2.

Why per-lane and not a single global scale: a single global factor makes gel two a
uniform multiple of gel one, so within any single-sample group the two replicate
values are perfectly proportional and the group standard deviation, while nonzero,
carries no independent information and a lane-swap bug between the two gels would
not move it. A DIFFERENT factor per lane makes each replicate value move by a
different amount, so the aggregated mean/sd/cv are non-degenerate and a join that
paired the wrong lanes across gels would change the numbers. The factors are held
deliberately close to 1.0 (roughly +/- 15 percent) so the perturbed gel still reads
as a plausible replicate of the same samples rather than a different experiment.

Why this is synthetic, stated plainly: there is no real second replicate of this
gel yet. This fixture exists to exercise the aggregator's join / mean / sd / cv
path, which has never run with n>1. The concrete numbers are not scientific; only
the machinery they flow through is under test. When a real replicate is available,
replace this generated directory with the real one and re-bless the goldens; the
harness block that drives the second gel does not care which it is.

Repo-weight contingency: manual_lane_profiles.csv is ~1.3 MB per gel (15 lanes x
1000 migration samples), so every committed replicate adds ~1.3 MB to the fixture
tree. --downsample K keeps every K-th migration row per lane (K=1 keeps all),
shrinking a fixture roughly K-fold if the tree ever needs to be trimmed. The
committed fixture is generated with K=1 (full resolution) so its goldens match a
full-resolution run; downsampling would move every golden and is offered only as
the lever to pull if fixture weight becomes a problem, not as the default.

Run (from images/), by hand, once:
    uv run tests/make_replicate_fixture.py
    uv run tests/make_replicate_fixture.py --downsample 10   (smaller, moves goldens)

Stdlib only. Flat and procedural; emit_message / die are the only helpers,
duplicated verbatim per the house convention.
"""

import argparse
import csv
import datetime
import pathlib
import shutil
import sys

# The first (real) gel this replicate is derived from, and the second (synthetic)
# gel directory this tool writes. Both names end in _gel_analysis because analyze_gel
# and serve derive gel_id and the expected TIFF stem from the directory name; the two
# names differ so the two gels carry different gel_ids (a hard requirement of the
# aggregator's "gel_id is unique" gate). The second name is the same sample imaged the
# next day at the next scan slot -- a plausible replicate label, not a real capture.
SOURCE_GEL_DIRECTORY_NAME = (
    "20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis"
)
REPLICATE_GEL_DIRECTORY_NAME = (
    "20260819_rotated_LM-0008_s0013_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis"
)

PROFILE_CSV_FILENAME = "manual_lane_profiles.csv"
SAMPLE_SHEET_FILENAME = "sample_sheet.csv"
RAW_VALUE_COLUMN_NAME = "raw_value"
LANE_INDEX_COLUMN_NAME = "lane_index"
MIGRATION_PIXELS_COLUMN_NAME = "migration_position_pixels"

# Per-lane multiplicative perturbation applied to raw_value. Keyed by lane_index
# (1..15 in this gel). Values sit within roughly +/- 15 percent of 1.0 and are
# intentionally non-uniform across lanes (see module docstring for why per-lane and
# not global). Lane 7 is BK_blue, the sample the harness goldens the aggregate mean
# on; its factor (1.08) is what makes the n=2 BK_blue mean differ from gel one's
# value, so the join is genuinely exercised rather than averaging two equal numbers.
# Any lane absent from this map would keep raw_value unchanged (factor 1.0), but all
# fifteen are listed explicitly so a lane is never silently left unperturbed.
PER_LANE_RAW_VALUE_FACTOR = {
    1: 1.00,
    2: 0.97,
    3: 1.05,
    4: 1.00,
    5: 0.91,
    6: 1.12,
    7: 1.08,
    8: 0.88,
    9: 1.03,
    10: 0.94,
    11: 1.10,
    12: 0.86,
    13: 1.06,
    14: 0.92,
    15: 1.00,
}

# raw_value is written by the ImageJ export with six decimal places; matching that
# precision keeps the replicate file byte-shaped like a real export and keeps the
# fixture hash stable across regenerations on any machine.
RAW_VALUE_OUTPUT_DECIMALS = 6


def emit_message(source_tag, message_text):
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    print("[" + source_tag + "] " + message_text, file=sys.stderr)


def die(source_tag, message_text):
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


argument_parser = argparse.ArgumentParser(
    prog="make_replicate_fixture.py",
    description="Author the synthetic second-replicate gel fixture for the n>1 aggregate test.",
    allow_abbrev=False,
)
argument_parser.add_argument(
    "--downsample",
    dest="downsample_stride",
    type=int,
    default=1,
    help="keep every K-th migration row per lane to shrink the fixture (K=1 keeps all; K>1 MOVES the goldens)",
)
parsed_arguments = argument_parser.parse_args()
downsample_stride = parsed_arguments.downsample_stride
if downsample_stride < 1:
    die("args", "--downsample must be a positive integer (got %d)" % downsample_stride)

# tests/ is this script's directory; fixtures/ sits beside it.
tests_directory = pathlib.Path(__file__).resolve().parent
fixtures_directory = tests_directory / "fixtures"
source_gel_directory = fixtures_directory / SOURCE_GEL_DIRECTORY_NAME
replicate_gel_directory = fixtures_directory / REPLICATE_GEL_DIRECTORY_NAME

source_profile_path = source_gel_directory / PROFILE_CSV_FILENAME
source_sample_sheet_path = source_gel_directory / SAMPLE_SHEET_FILENAME
if not source_profile_path.is_file():
    die("input", "source profile not found: " + str(source_profile_path))
if not source_sample_sheet_path.is_file():
    die("input", "source sample sheet not found: " + str(source_sample_sheet_path))

# Refuse to silently overwrite an existing replicate fixture; regeneration is a
# deliberate act (delete it first) so a stray run cannot move the committed goldens.
if replicate_gel_directory.exists():
    die(
        "output",
        "replicate gel directory already exists: "
        + str(replicate_gel_directory)
        + " (delete it first to regenerate)",
    )
replicate_gel_directory.mkdir(parents=True)

# Read the source profile, scale raw_value per lane, optionally downsample, and write
# the replicate profile. utf-8-sig strips an Excel BOM on read; the replicate is
# written as plain utf-8 with the same column order the source carried.
with source_profile_path.open(
    newline="", encoding="utf-8-sig"
) as source_profile_handle:
    source_profile_reader = csv.DictReader(source_profile_handle)
    profile_column_names = source_profile_reader.fieldnames or []
    source_profile_rows = [
        {
            column_name: (cell.strip() if isinstance(cell, str) else cell)
            for column_name, cell in raw_row.items()
        }
        for raw_row in source_profile_reader
    ]

for required_column in (
    LANE_INDEX_COLUMN_NAME,
    RAW_VALUE_COLUMN_NAME,
    MIGRATION_PIXELS_COLUMN_NAME,
):
    if required_column not in profile_column_names:
        die("input", "source profile missing required column: " + required_column)

# Downsampling keeps every K-th migration sample WITHIN each lane, so every lane is
# thinned identically and the region window still lands on the same migration range.
# Keeping is decided by each row's own migration_position_pixels modulo the stride
# rather than by a running counter, so lane boundaries in the file cannot shift which
# rows survive.
replicate_profile_rows = []
lanes_seen_in_order = []
for source_profile_row in source_profile_rows:
    lane_index = int(source_profile_row[LANE_INDEX_COLUMN_NAME])
    migration_pixels = int(source_profile_row[MIGRATION_PIXELS_COLUMN_NAME])
    if downsample_stride > 1 and (migration_pixels % downsample_stride) != 0:
        continue
    if lane_index not in PER_LANE_RAW_VALUE_FACTOR:
        die(
            "input",
            "lane %d has no perturbation factor; add it to PER_LANE_RAW_VALUE_FACTOR"
            % lane_index,
        )
    scaled_raw_value = (
        float(source_profile_row[RAW_VALUE_COLUMN_NAME])
        * PER_LANE_RAW_VALUE_FACTOR[lane_index]
    )
    replicate_profile_row = dict(source_profile_row)
    replicate_profile_row[RAW_VALUE_COLUMN_NAME] = "%.*f" % (
        RAW_VALUE_OUTPUT_DECIMALS,
        scaled_raw_value,
    )
    replicate_profile_rows.append(replicate_profile_row)
    if lane_index not in lanes_seen_in_order:
        lanes_seen_in_order.append(lane_index)

replicate_profile_path = replicate_gel_directory / PROFILE_CSV_FILENAME
with replicate_profile_path.open("w", newline="", encoding="utf-8") as replicate_handle:
    profile_writer = csv.DictWriter(replicate_handle, fieldnames=profile_column_names)
    profile_writer.writeheader()
    profile_writer.writerows(replicate_profile_rows)

# The sample sheet is copied verbatim: the aggregator joins by sample_label, so the
# replicate must carry the identical labels for its samples to co-group with gel one.
replicate_sample_sheet_path = replicate_gel_directory / SAMPLE_SHEET_FILENAME
shutil.copy2(source_sample_sheet_path, replicate_sample_sheet_path)

emit_message(
    "output",
    "wrote %s (%d rows across %d lanes, downsample stride %d)"
    % (
        str(replicate_profile_path),
        len(replicate_profile_rows),
        len(lanes_seen_in_order),
        downsample_stride,
    ),
)
emit_message("output", "copied sample sheet to " + str(replicate_sample_sheet_path))
emit_message(
    "summary",
    "replicate gel fixture authored; re-bless the harness goldens with --update",
)
