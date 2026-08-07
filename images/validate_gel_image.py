# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "tifffile",
#     "numpy",
#     "matplotlib",
# ]
# ///
r"""
Stage 1 of the gel densitometry pipeline: interrogate one TIFF, validate it, and
report. Emits input_file_validation_report.json, input_file_histogram.png, and
input_file_crop_preview.png when a preprocessing sidecar exists. Measures nothing,
finds no lanes, finds no bands, writes no modified image.

Scope is DESIGN.md section 7 stage 1 with the content of section 6, less the
.img-versus-.tif scatter, which was deferred by decision: it compares two
containers without a calibrated reference, so it can show that they differ but not
which is linear in signal. The linearity question is therefore recorded as
unverified rather than answered. See the linearity_evidence block in the report.

Shape is dictated by CONVENTIONS.md: flat procedural, no helpers beyond the tagged
stderr emitter and the fail-fast exit, full descriptive names carrying units,
ASCII only, comments stating why.

Two deliberate departures from the stage 0 scripts, both required by this script
having a report as its deliverable:

1. Failed checks accumulate in VALIDATION_FINDINGS and the report is always
   written, at a single exit point. die() is reserved for the unrecoverable, where
   there is nothing to report about: the path is unusable, the file will not open,
   it is not a TIFF, the output directory is missing. A validation failure that
   produced no validation report would be the worst available outcome for this
   script, and it is what calling die() at every hard stop would produce.
2. The output directory must already exist. Stage 0 creates it. Creating it here
   would let stage 1 run against a directory stage 0 never validated.

Usage. Single quotes, because real filenames here contain spaces, commas and
square brackets, and square brackets are glob metacharacters in bash:

    uv run validate_gel_image.py '/mnt/c/Users/liusm/.../20220303-...-[Phosphor].tif'
    uv run validate_gel_image.py --page-index 0 '<multi-page>.tif'
    uv run validate_gel_image.py -- '-filename-starting-with-a-dash.tif'
"""

import argparse
import datetime
import json
import math
import os
import pathlib
import re
import stat
import sys

import matplotlib

# Selected before pyplot is imported. There is no display in WSL by default, and
# the default interactive backend fails at import rather than at draw time, which
# would make this script unable to run for a reason unrelated to its job.
matplotlib.use("Agg")

import matplotlib.patches
import matplotlib.pyplot
import numpy
import tifffile

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIRECTORY_NAME_SUFFIX = "_gel_analysis"
VALIDATION_REPORT_FILENAME = "input_file_validation_report.json"
HISTOGRAM_IMAGE_FILENAME = "input_file_histogram.png"
CROP_PREVIEW_IMAGE_FILENAME = "input_file_crop_preview.png"

PREPROCESS_SIDECAR_FILENAME_SUFFIX = "_preprocess.txt"
INF_SIDECAR_SUFFIX = ".inf"

# Bumped whenever the report's shape changes, so stage 2 can refuse a report it
# does not understand instead of reading a missing key as absent.
VALIDATION_REPORT_SCHEMA_VERSION = 1

SUPPORTED_SIDECAR_SCHEMA_VERSION = 2

# Coordinates are meaningless without the frame they were measured in. This exists
# so that if the convention ever changes, old sidecars fail loudly instead of being
# silently reinterpreted against the new one.
REQUIRED_SIDECAR_FRAME_NAME = "raw_as_opened_unflipped"

# Duplicated from preprocess_sidecar_template.txt, which is the authority.
# CONVENTIONS.md section 1 forbids sharing between scripts. Drift is caught rather
# than tolerated: a missing key and an unknown key are both hard stops, so adding
# a key to the template without adding it here fails on the next run.
REQUIRED_SIDECAR_KEY_NAMES = (
    "schema_version",
    "measured_in_frame",
    "measured_against_input_filename",
    "measured_against_image_width_pixels",
    "measured_against_image_height_pixels",
    "gel_migration_axis",
    "coordinate_unit",
    "landmark_a_x",
    "landmark_a_y",
    "landmark_b_x",
    "landmark_b_y",
    "rotation_landmark_description",
    "crop_x",
    "crop_y",
    "crop_width",
    "crop_height",
    "expected_lane_count",
    "notes",
)

# The two free-text fields. Everything else must be non-empty.
SIDECAR_KEY_NAMES_PERMITTED_EMPTY = ("rotation_landmark_description", "notes")

# Always integers, and always pixels regardless of coordinate_unit: schema and the
# image dimensions read from the tags, and the lane count.
SIDECAR_INTEGER_KEY_NAMES = (
    "schema_version",
    "measured_against_image_width_pixels",
    "measured_against_image_height_pixels",
    "expected_lane_count",
)

# The geometry values. Numeric rather than integer because under
# coordinate_unit=centimetres they are fractional centimetres; they are converted
# to pixels once, below, using the pixel size read from the file's own tags.
SIDECAR_NUMERIC_KEY_NAMES = (
    "landmark_a_x",
    "landmark_a_y",
    "landmark_b_x",
    "landmark_b_y",
    "crop_x",
    "crop_y",
    "crop_width",
    "crop_height",
)

# Closed vocabularies. gel_migration_axis has no default (an old gel imaged in the
# other orientation must not be processed under a silent assumption); coordinate_unit
# distinguishes raw pixels from the calibrated centimetres Fiji reports.
SIDECAR_ENUM_KEY_ALLOWED_VALUES = {
    "gel_migration_axis": ("horizontal", "vertical"),
    "coordinate_unit": ("pixels", "centimetres"),
}

PRIVATE_TAG_CODE_FLOOR = 65000

# Generous, unlike the console inventory's limit: this goes into a JSON file
# nobody has to read on a terminal, and a truncated ImageDescription is exactly
# the provenance this script exists to preserve.
TAG_VALUE_CHARACTER_LIMIT = 20000

CONTAINER_MAXIMUM_VALUE_16_BIT = 65535
CONTAINER_VALUE_SLOT_COUNT_16_BIT = 65536

# 8 bits per sample means quantification is off the table, per DESIGN.md section 6.
MINIMUM_QUANTIFIABLE_BITS_PER_SAMPLE = 16

# A stride in the present values is the signature of fewer real bits left-shifted
# into a 16-bit container, which is the case DESIGN.md 5.12 warns about: a naive
# equality test against 65535 then reports zero saturation on a clipped image.
STRIDE_TO_EFFECTIVE_BITS_PER_SAMPLE = {1: 16, 2: 15, 4: 14, 8: 13, 16: 12}

# A pixel exceeding every one of its eight neighbours by this many robust standard
# deviations is an isolated spike: a cosmic ray or screen contamination on a
# phosphor plate, a cosmic ray or hot pixel on a cooled CCD. Reported, never
# filtered; DESIGN.md section 6 leaves the decision to the user.
ISOLATION_MARGIN_ROBUST_SIGMAS = 5.0

# Scale factor converting a median absolute deviation into a standard deviation
# for normally distributed data. Used instead of the standard deviation itself
# because a gel's band pixels inflate the latter by an amount that depends on how
# much sample was loaded.
MEDIAN_ABSOLUTE_DEVIATION_TO_SIGMA = 1.4826

# Clipping shows up as an anomalous spike at the maximum rather than as the
# maximum being any particular number. The comparison population is the present
# values immediately below it.
SATURATION_SPIKE_NEIGHBOURHOOD_VALUE_COUNT = 100
SATURATION_SPIKE_RATIO_WARNING_THRESHOLD = 10.0

# A large population at the floor means the image was already background
# subtracted or clipped at zero, which changes what a baseline means.
FLOOR_POPULATION_FRACTION_WARNING_THRESHOLD = 0.01

# Beyond this the landmarks were almost certainly not the level feature they were
# meant to be, or the two pairs were entered in the wrong order.
MAXIMUM_PLAUSIBLE_TILT_DEGREES = 10.0

# Precision of the derived angle is set by the span: a 1 px click error over 400 px
# is 0.14 degrees, and 0.5 degrees walks a band a full band height across a 1125 px
# image. Below the floor the angle is not worth deriving.
MINIMUM_LANDMARK_SPAN_PIXELS = 100
LANDMARK_SPAN_WARNING_PIXELS = 400

# Sanity bounds on the lane pitch implied by the crop width and the expected lane
# count, converted through the pixel size read from this file's own tags.
# DESIGN.md section 2 records 5 mm as the real pitch.
MINIMUM_PLAUSIBLE_LANE_PITCH_MILLIMETRES = 2.0
MAXIMUM_PLAUSIBLE_LANE_PITCH_MILLIMETRES = 20.0

CROP_AREA_FRACTION_WARNING_FLOOR = 0.02
CROP_AREA_FRACTION_WARNING_CEILING = 0.98

# The preview exists to be looked at, not measured, so it is downsampled by an
# integer stride to keep the file small.
PREVIEW_MAXIMUM_DIMENSION_PIXELS = 1400

HISTOGRAM_FIGURE_SIZE_INCHES = (11.0, 6.0)
PREVIEW_FIGURE_SIZE_INCHES = (9.0, 9.0)
FIGURE_DOTS_PER_INCH = 110

# PhotometricInterpretation on both instruments' files is MINISWHITE, so Fiji
# renders them with an inverting lookup table and bands appear dark. A default
# colormap would render the same array with bright bands, and anyone comparing the
# two concludes the pipeline inverted the image and starts flipping things.
DISPLAY_COLORMAP_NAME = "gray_r"

READABILITY_PROBE_BYTE_COUNT = 1

SURROUNDING_QUOTE_CHARACTERS = ('"', "'")

# Drive letter NOT followed by a separator: the residue of an unquoted Windows
# path whose backslashes the shell consumed as escapes.
WINDOWS_DRIVE_RELATIVE_PATH_PATTERN = re.compile(r"^([A-Za-z]):(?![/\\])")

NON_REGULAR_FILE_TYPE_PREDICATES = (
    ("directory", stat.S_ISDIR),
    ("named pipe", stat.S_ISFIFO),
    ("socket", stat.S_ISSOCK),
    ("character device", stat.S_ISCHR),
    ("block device", stat.S_ISBLK),
)

# Measuring one of these would be silent nonsense: they are 8-bit display
# renderings written alongside the 16-bit TIFF by the Amersham Imager 680.
DISPLAY_RENDERING_SUFFIXES = (".png", ".jpg", ".jpeg")

# Located by pattern rather than by line number. DESIGN.md section 12 established
# that section 2's positional line numbers are off by one from position 3 onward,
# and the omitted line turned out to be blank rather than the four characters
# section 12's arithmetic predicted, leaving 4 bytes still unexplained. Pattern
# matching is immune to all of that; positional indexing is not.
INF_HUMAN_TIMESTAMP_PATTERN = re.compile(
    r"^[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9]\d \d{2}:\d{2}:\d{2} \d{4}$"
)
INF_EPOCH_TIMESTAMP_PATTERN = re.compile(r"^\d{9,11}$")
INF_KEYED_LINE_PATTERN = re.compile(r"^([SH]),([^=]+)=(.*)$")

# =============================================================================
# The two permitted helpers (CONVENTIONS.md section 1)
# =============================================================================

ACCUMULATED_RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    # stderr, so stdout carries only the machine-readable summary and stays
    # pipeable. Every message also accumulates for the report, which
    # CONVENTIONS.md section 10 requires: the run log is part of the record.
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    ACCUMULATED_RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    # Reserved for the unrecoverable, where there is nothing to write a report
    # about. Every check that can fail with the file still readable appends to
    # VALIDATION_FINDINGS instead, so that the report exists either way.
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


# Findings accumulate here rather than exiting at the first failure. Each entry is
# a plain dict appended inline at the check site; no third helper, per
# CONVENTIONS.md section 1. is_hard_stop distinguishes a check whose failure
# invalidates the file from one that is merely reported, and only a failed hard
# stop changes the exit code.
VALIDATION_FINDINGS = []

# =============================================================================
# Arguments
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="validate_gel_image.py",
    description=(
        "Stage 1: interrogate and validate one TIFF, emit a validation report, a "
        "histogram, and a crop preview when a preprocessing sidecar exists."
    ),
    # An abbreviated flag that is unambiguous today becomes ambiguous the moment a
    # second flag shares its prefix, silently changing a command that used to work.
    allow_abbrev=False,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=(
        "Quote paths with single quotes. Real filenames contain spaces, commas\n"
        "and square brackets, and square brackets are glob metacharacters in\n"
        "bash. Use -- before a filename that starts with a dash.\n"
    ),
)
argument_parser.add_argument(
    "input_tiff_path",
    help="Path to the TIFF to validate. The .tif is the read path; never the .img or .gel.",
)
argument_parser.add_argument(
    "--page-index",
    type=int,
    default=None,
    help=(
        "TIFF page (IFD) to read. Required when the file has more than one page; "
        "DESIGN.md section 6 forbids defaulting to page zero."
    ),
)
argument_parser.add_argument(
    "--output-parent-directory",
    default=None,
    help=(
        "Directory holding <stem>" + OUTPUT_DIRECTORY_NAME_SUFFIX + ". Defaults to "
        "the directory holding the input, per DESIGN.md 5.9. Must already exist, "
        "along with the per-input directory that stage 0 created."
    ),
)
parsed_arguments = argument_parser.parse_args()

emit_message(
    "arguments", "input tiff path as given: " + repr(parsed_arguments.input_tiff_path)
)

# =============================================================================
# Path normalization
# =============================================================================

# Duplicated from prepare_gel_image.py rather than imported, per CONVENTIONS.md
# section 1. What was NOT duplicated, and why: the Windows drive-letter conversion
# and the UNC rejection. Stage 0 has already resolved and reported the path by the
# time this script runs, and it is documented as taking the /mnt form. The
# de-escaped-paste hard stop is kept, because that mistake is otherwise
# unintelligible.
working_path_text = parsed_arguments.input_tiff_path

whitespace_stripped_path_text = working_path_text.strip(" \t\r\n")
if whitespace_stripped_path_text != working_path_text:
    emit_message("path", "stripped surrounding whitespace from the supplied path")
    working_path_text = whitespace_stripped_path_text
if working_path_text == "":
    die("path", "path is empty after stripping surrounding whitespace")

for quote_character in SURROUNDING_QUOTE_CHARACTERS:
    if (
        len(working_path_text) >= 2
        and working_path_text.startswith(quote_character)
        and working_path_text.endswith(quote_character)
    ):
        emit_message(
            "path", "removed a surrounding pair of " + quote_character + " characters"
        )
        working_path_text = working_path_text[1:-1]
        break
if working_path_text == "":
    die("path", "path is empty after removing surrounding quotes")

if WINDOWS_DRIVE_RELATIVE_PATH_PATTERN.match(working_path_text) is not None:
    die(
        "path",
        "the path " + repr(working_path_text) + " has a drive letter but no separator "
        "after it. Backslash is the shell escape character, so an unquoted Windows path "
        "arrives with its separators already deleted and cannot be reconstructed. "
        "Re-run with the path inside single quotes, or supply the /mnt form.",
    )

if working_path_text.startswith("~"):
    try:
        tilde_expanded_path_text = os.path.expanduser(working_path_text)
    except RuntimeError as tilde_expansion_error:
        die("path", "cannot expand the leading tilde: " + str(tilde_expansion_error))
    if tilde_expanded_path_text == working_path_text:
        die(
            "path",
            "the leading tilde in " + repr(working_path_text) + " did not expand",
        )
    emit_message("path", "expanded leading tilde to " + tilde_expanded_path_text)
    working_path_text = tilde_expanded_path_text

input_tiff_absolute_path = pathlib.Path(os.path.abspath(working_path_text))
input_tiff_physical_path = pathlib.Path(working_path_text).resolve()
if input_tiff_physical_path != input_tiff_absolute_path:
    emit_message(
        "path",
        "path traverses a symlink: given form "
        + str(input_tiff_absolute_path)
        + " resolves to physical form "
        + str(input_tiff_physical_path),
    )
emit_message("path", "normalized to " + str(input_tiff_absolute_path))

# .suffix and .stem are safe on names containing dots, which the Amersham Imager
# 680 filenames do: 2026.07.10_20.06.32_LM-0008_s001_....tif. .suffixes is not,
# returning ['.07', '.10_20', '.06', '.32_LM-...', '.tif'] for that name, so it is
# never used here.
if input_tiff_absolute_path.suffix.lower() in DISPLAY_RENDERING_SUFFIXES:
    die(
        "arguments",
        "the input is a display rendering, not the 16-bit TIFF: "
        + str(input_tiff_absolute_path)
        + ". The Amersham Imager 680 writes a .png and "
        "a .jpg beside the .tif; both are 8-bit and quantifying either is silent "
        "nonsense. Pass the .tif.",
    )
if input_tiff_absolute_path.suffix.lower() == INF_SIDECAR_SUFFIX:
    die(
        "arguments",
        "the input is the .inf sidecar itself: "
        + str(input_tiff_absolute_path)
        + ". Pass the .tif; the sidecars are located from its stem.",
    )

# =============================================================================
# Existence, file type and readability
# =============================================================================

if not os.path.lexists(input_tiff_absolute_path):
    first_missing_path_component = input_tiff_absolute_path
    for candidate_ancestor_path in [input_tiff_absolute_path] + list(
        input_tiff_absolute_path.parents
    ):
        if os.path.lexists(candidate_ancestor_path):
            break
        first_missing_path_component = candidate_ancestor_path
    die(
        "existence",
        "input not found: "
        + str(input_tiff_absolute_path)
        + ". Highest path component that does not exist: "
        + str(first_missing_path_component),
    )

input_tiff_link_status = os.lstat(input_tiff_absolute_path)
if stat.S_ISLNK(input_tiff_link_status.st_mode) and not os.path.exists(
    input_tiff_absolute_path
):
    die(
        "existence",
        str(input_tiff_absolute_path)
        + " is a symlink whose target "
        + str(input_tiff_physical_path)
        + " does not exist.",
    )

input_tiff_file_status = os.stat(input_tiff_absolute_path)

# Rejecting non-regular files before opening anything: opening a named pipe blocks
# until a writer appears, so the readability probe below would hang rather than
# fail.
for (
    non_regular_type_name,
    non_regular_type_predicate,
) in NON_REGULAR_FILE_TYPE_PREDICATES:
    if non_regular_type_predicate(input_tiff_file_status.st_mode):
        directory_hint_text = ""
        if non_regular_type_name == "directory":
            directory_hint_text = (
                " Directory mode is stage 4 in DESIGN.md section 7; pass a single file."
            )
        die(
            "file type",
            str(input_tiff_absolute_path)
            + " is a "
            + non_regular_type_name
            + ", not a regular file."
            + directory_hint_text,
        )
if not stat.S_ISREG(input_tiff_file_status.st_mode):
    die(
        "file type",
        str(input_tiff_absolute_path)
        + " is not a regular file (st_mode "
        + oct(input_tiff_file_status.st_mode)
        + ").",
    )

# DESIGN.md section 3 names reading a partially synced Dropbox file as a real risk,
# and a zero byte file is exactly what an interrupted sync leaves behind.
if input_tiff_file_status.st_size == 0:
    die(
        "file size",
        str(input_tiff_absolute_path) + " is zero bytes. An interrupted or "
        "placeholder-only Dropbox sync produces exactly this.",
    )

# os.access is deliberately not used. DESIGN.md section 3 records that /mnt/c is
# DrvFs and reports 777 regardless of real permissions, and access() run as root
# returns true for nearly everything. Opening the file and reading a byte is the
# only test that answers the question actually being asked.
try:
    with open(input_tiff_absolute_path, "rb") as input_tiff_file_handle:
        readability_probe_bytes = input_tiff_file_handle.read(
            READABILITY_PROBE_BYTE_COUNT
        )
except PermissionError as permission_error:
    die(
        "readability",
        str(input_tiff_absolute_path)
        + " exists but cannot be opened for reading: "
        + str(permission_error),
    )
except OSError as open_error:
    die(
        "readability",
        str(input_tiff_absolute_path) + " failed to open: " + str(open_error),
    )
if len(readability_probe_bytes) != READABILITY_PROBE_BYTE_COUNT:
    die(
        "readability",
        str(input_tiff_absolute_path)
        + " opened but returned "
        + str(len(readability_probe_bytes))
        + " bytes instead of "
        + str(READABILITY_PROBE_BYTE_COUNT)
        + ".",
    )
emit_message(
    "readability",
    "opened and read 1 byte; reported size is "
    + str(input_tiff_file_status.st_size)
    + " bytes",
)

# =============================================================================
# Output directory, which stage 0 must already have created
# =============================================================================

if parsed_arguments.output_parent_directory is None:
    output_parent_directory_path = input_tiff_absolute_path.parent
else:
    output_parent_directory_path = pathlib.Path(
        os.path.abspath(
            os.path.expanduser(
                parsed_arguments.output_parent_directory.strip(" \t\r\n")
            )
        )
    )

output_directory_path = output_parent_directory_path / (
    input_tiff_absolute_path.stem + OUTPUT_DIRECTORY_NAME_SUFFIX
)

# Not created here. Stage 0 creates it, having verified the parent, and running
# stage 1 against a directory stage 0 never validated is the failure this prevents.
# Note that the .img and the .tif of one scan share a stem and therefore share this
# directory: one directory per scan, not per container, which is intended.
if not output_directory_path.is_dir():
    die(
        "output directory",
        str(output_directory_path) + " does not exist. Run prepare_gel_image.py on "
        "this file first; stage 1 does not create the output directory, so that it "
        "cannot run against a directory stage 0 never validated.",
    )

validation_report_output_path = output_directory_path / VALIDATION_REPORT_FILENAME
histogram_image_output_path = output_directory_path / HISTOGRAM_IMAGE_FILENAME
crop_preview_image_output_path = output_directory_path / CROP_PREVIEW_IMAGE_FILENAME
emit_message("output directory", "writing into " + str(output_directory_path))

# =============================================================================
# Open the TIFF and select the page
# =============================================================================

try:
    tiff_file_handle = tifffile.TiffFile(input_tiff_absolute_path)
except Exception as tiff_open_error:
    # Deliberately broad. Any failure to parse the file as TIFF is unrecoverable
    # for this script, and the exception type is itself the diagnostic.
    die(
        "tiff",
        "tifffile could not parse "
        + str(input_tiff_absolute_path)
        + " as TIFF: "
        + type(tiff_open_error).__name__
        + ": "
        + str(tiff_open_error),
    )

tiff_page_count = len(tiff_file_handle.pages)
emit_message(
    "tiff",
    "opened; byte order "
    + str(tiff_file_handle.byteorder)
    + ", bigtiff "
    + str(tiff_file_handle.is_bigtiff)
    + ", "
    + str(tiff_page_count)
    + " page(s)",
)

# DESIGN.md section 6 forbids defaulting to page zero. A single-page file has no
# ambiguity to resolve, so the requirement bites only above one.
if tiff_page_count == 0:
    die(
        "tiff", str(input_tiff_absolute_path) + " parsed as TIFF but contains no pages."
    )
if parsed_arguments.page_index is None:
    if tiff_page_count > 1:
        # The one hard stop that deliberately produces no report: without a page
        # there is nothing to inventory, and guessing is what section 6 forbids. The
        # message lists the pages so the choice can be made without a second tool.
        page_summary_texts = []
        for candidate_page_index, candidate_page in enumerate(tiff_file_handle.pages):
            page_summary_texts.append(
                str(candidate_page_index)
                + ": shape "
                + str(getattr(candidate_page, "shape", None))
                + ", dtype "
                + str(getattr(candidate_page, "dtype", None))
                + ", "
                + str(getattr(candidate_page, "bitspersample", None))
                + " bits"
            )
        die(
            "tiff",
            str(input_tiff_absolute_path)
            + " has "
            + str(tiff_page_count)
            + " pages, so "
            "--page-index is required: DESIGN.md section 6 forbids defaulting to page "
            "zero, because a reduced-resolution or differently-exposed page would be "
            "quantified without complaint. Pages are ["
            + "; ".join(page_summary_texts)
            + "].",
        )
    selected_page_index = 0
else:
    selected_page_index = parsed_arguments.page_index
    if selected_page_index < 0 or selected_page_index >= tiff_page_count:
        die(
            "tiff",
            "--page-index "
            + str(selected_page_index)
            + " is outside the range 0 to "
            + str(tiff_page_count - 1)
            + " for this file.",
        )
selected_tiff_page = tiff_file_handle.pages[selected_page_index]
emit_message("tiff", "reading page " + str(selected_page_index))

# =============================================================================
# Tag inventory, verbatim, into the report
# =============================================================================

# Values are stringified with repr for the same reason inventory_tiff_tags.py does
# it: a tag may hold arbitrary bytes, and CONVENTIONS.md section 4 requires ASCII.
# Enum tags additionally get their integer and their name recorded separately,
# because tifffile's enums are IntEnum and json.dumps serializes them silently as
# bare integers. PhotometricInterpretation reduced to "0" in the report would lose
# exactly the field whose name carries the most meaning.
standard_tag_records = []
private_tag_records = []
# Collected during this pass rather than looked up by code afterwards. Looking up
# by code returns the first tag with that code, so with ImageDescription present
# twice the first copy gets parsed twice and the second is never seen. Found by
# running: tifffile writes its own {"shape": ...} metadata as a second tag 270, and
# stage 0 testing saw duplicate 270s on a real file too.
image_description_tag_values = []
for tiff_tag in selected_tiff_page.tags:
    if tiff_tag.name == "ImageDescription":
        image_description_tag_values.append(tiff_tag.value)
    tag_value_text = repr(tiff_tag.value)
    tag_value_was_truncated = len(tag_value_text) > TAG_VALUE_CHARACTER_LIMIT
    if tag_value_was_truncated:
        tag_value_text = tag_value_text[:TAG_VALUE_CHARACTER_LIMIT]
    tag_record = {
        "code": int(tiff_tag.code),
        "name": str(tiff_tag.name),
        "dtype": str(tiff_tag.dtype),
        "count": int(tiff_tag.count),
        "value_repr": tag_value_text,
        "value_repr_was_truncated": tag_value_was_truncated,
        "enum_integer_value": None,
        "enum_name": None,
    }
    if isinstance(tiff_tag.value, int) and hasattr(tiff_tag.value, "name"):
        tag_record["enum_integer_value"] = int(tiff_tag.value)
        tag_record["enum_name"] = str(tiff_tag.value.name)
    if tiff_tag.code >= PRIVATE_TAG_CODE_FLOOR:
        private_tag_records.append(tag_record)
    else:
        standard_tag_records.append(tag_record)
standard_tag_records.sort(key=lambda record: record["code"])
private_tag_records.sort(key=lambda record: record["code"])
emit_message(
    "tiff",
    "inventoried "
    + str(len(standard_tag_records))
    + " standard tags and "
    + str(len(private_tag_records))
    + " private tags at or above "
    + str(PRIVATE_TAG_CODE_FLOOR),
)

# Duplicate tag codes are legal and do occur; tag 270 ImageDescription appeared
# twice during stage 0 testing. Anything reading a tag by name must therefore cope
# with several, which is why this is a list of values rather than a lookup.
tag_values_by_name = {}
for tag_record in standard_tag_records + private_tag_records:
    tag_values_by_name.setdefault(tag_record["name"], []).append(tag_record)
duplicated_tag_names = sorted(
    tag_name for tag_name, records in tag_values_by_name.items() if len(records) > 1
)
if len(duplicated_tag_names) > 0:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "tiff_duplicate_tag_codes",
            "status": "reported",
            "is_hard_stop": False,
            "detail": "tags present more than once, all copies recorded: "
            + ", ".join(duplicated_tag_names),
        }
    )

# =============================================================================
# TIFF hard stops (DESIGN.md section 6)
# =============================================================================

VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_opens_as_tiff",
        "status": "pass",
        "is_hard_stop": True,
        "detail": "tifffile parsed the file; byte order "
        + str(tiff_file_handle.byteorder)
        + ", bigtiff "
        + str(tiff_file_handle.is_bigtiff),
    }
)

VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_page_count_and_explicit_selection",
        "status": "pass",
        "is_hard_stop": True,
        "detail": "page count "
        + str(tiff_page_count)
        + ", page "
        + str(selected_page_index)
        + " selected "
        + (
            "by default because there is only one"
            if tiff_page_count == 1
            else "explicitly via --page-index"
        ),
    }
)

page_samples_per_pixel = selected_tiff_page.samplesperpixel
VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_samples_per_pixel_is_one",
        "status": "pass" if page_samples_per_pixel == 1 else "fail",
        "is_hard_stop": True,
        "detail": "SamplesPerPixel is "
        + str(page_samples_per_pixel)
        + (
            ". A value of 3 means RGB, which means something rendered this file "
            "rather than the instrument writing it."
            if page_samples_per_pixel != 1
            else ""
        ),
    }
)

page_bits_per_sample = selected_tiff_page.bitspersample
if isinstance(page_bits_per_sample, tuple):
    page_bits_per_sample = page_bits_per_sample[0]
VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_bits_per_sample_permits_quantification",
        "status": "pass"
        if page_bits_per_sample >= MINIMUM_QUANTIFIABLE_BITS_PER_SAMPLE
        else "fail",
        "is_hard_stop": True,
        "detail": "BitsPerSample is " + str(page_bits_per_sample) + "; the minimum for "
        "quantification is "
        + str(MINIMUM_QUANTIFIABLE_BITS_PER_SAMPLE)
        + (
            ". An 8-bit file is a display rendering, or a 16-bit file that was "
            "converted; either way quantification is off the table."
            if page_bits_per_sample < MINIMUM_QUANTIFIABLE_BITS_PER_SAMPLE
            else ""
        ),
    }
)

page_sample_format = selected_tiff_page.sampleformat
VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_sample_format",
        "status": "reported",
        "is_hard_stop": False,
        "detail": "SampleFormat is "
        + repr(page_sample_format)
        + " (1 unsigned integer, 2 signed integer, 3 floating point)",
    }
)

page_compression = selected_tiff_page.compression
page_compression_name = getattr(page_compression, "name", str(page_compression))
compression_is_lossless = int(page_compression) in (1, 5, 8, 32773, 32946, 34925)
VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_compression_is_not_lossy",
        "status": "pass" if compression_is_lossless else "fail",
        "is_hard_stop": True,
        "detail": "Compression is "
        + page_compression_name
        + " ("
        + str(int(page_compression))
        + ")"
        + (
            ""
            if compression_is_lossless
            else ". Lossy compression has already discarded pixel values, so no "
            "quantification from this file is defensible."
        ),
    }
)

page_photometric = selected_tiff_page.photometric
page_photometric_name = getattr(page_photometric, "name", str(page_photometric))
VALIDATION_FINDINGS.append(
    {
        "check_name": "tiff_photometric_interpretation",
        "status": "reported",
        "is_hard_stop": False,
        "detail": "PhotometricInterpretation is "
        + page_photometric_name
        + " ("
        + str(int(page_photometric))
        + "). Cross-checked against the "
        "histogram-inferred signal direction below.",
    }
)

page_image_width_pixels = int(selected_tiff_page.tags["ImageWidth"].value)
page_image_length_pixels = int(selected_tiff_page.tags["ImageLength"].value)

# =============================================================================
# Resolution, converted to micrometres per pixel from this file's own tags
# =============================================================================

# CONVENTIONS.md section 10 requires physical lengths in config to be converted
# using the pixel size read from the input's own header, and the two instruments
# make that non-negotiable rather than tidy: the Typhoon is 200 micrometres per
# pixel and the Amersham Imager 680 is 79.9, and repeat 2 of the Typhoon set has
# different dimensions from repeats 1 and 3.
resolution_unit_value = None
if "ResolutionUnit" in selected_tiff_page.tags:
    resolution_unit_value = int(selected_tiff_page.tags["ResolutionUnit"].value)
resolution_unit_to_micrometres = {2: 25400.0, 3: 10000.0}

micrometres_per_pixel_by_axis = {}
for resolution_tag_name, axis_name in (("XResolution", "x"), ("YResolution", "y")):
    if resolution_tag_name not in selected_tiff_page.tags:
        continue
    resolution_rational = selected_tiff_page.tags[resolution_tag_name].value
    if isinstance(resolution_rational, tuple) and len(resolution_rational) == 2:
        resolution_numerator, resolution_denominator = resolution_rational
    else:
        resolution_numerator, resolution_denominator = float(resolution_rational), 1.0
    if resolution_denominator == 0 or resolution_numerator == 0:
        continue
    pixels_per_unit = float(resolution_numerator) / float(resolution_denominator)
    if resolution_unit_value not in resolution_unit_to_micrometres:
        continue
    micrometres_per_pixel_by_axis[axis_name] = (
        resolution_unit_to_micrometres[resolution_unit_value] / pixels_per_unit
    )

if len(micrometres_per_pixel_by_axis) == 2:
    axis_difference = abs(
        micrometres_per_pixel_by_axis["x"] - micrometres_per_pixel_by_axis["y"]
    )
    VALIDATION_FINDINGS.append(
        {
            "check_name": "tiff_pixel_size_axes_agree",
            "status": "pass" if axis_difference < 1e-6 else "fail",
            "is_hard_stop": True,
            "detail": "x is %.4f and y is %.4f micrometres per pixel"
            % (micrometres_per_pixel_by_axis["x"], micrometres_per_pixel_by_axis["y"])
            + (
                ""
                if axis_difference < 1e-6
                else ". Anisotropic pixels are not handled: every millimetre-to-pixel "
                "conversion downstream assumes one scale factor."
            ),
        }
    )
    micrometres_per_pixel = micrometres_per_pixel_by_axis["x"]
elif len(micrometres_per_pixel_by_axis) == 1:
    micrometres_per_pixel = list(micrometres_per_pixel_by_axis.values())[0]
    VALIDATION_FINDINGS.append(
        {
            "check_name": "tiff_pixel_size_axes_agree",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "only one resolution axis is present; assuming square pixels at "
            "%.4f micrometres" % micrometres_per_pixel,
        }
    )
else:
    micrometres_per_pixel = None
    VALIDATION_FINDINGS.append(
        {
            "check_name": "tiff_pixel_size_is_readable",
            "status": "fail",
            "is_hard_stop": True,
            "detail": "no usable XResolution/YResolution with a supported ResolutionUnit "
            "(got ResolutionUnit " + repr(resolution_unit_value) + "). Every "
            "millimetre value in the config converts through this, so there is "
            "nothing to convert with.",
        }
    )

# =============================================================================
# ImageDescription, parsed as ordered pairs because keys repeat
# =============================================================================

# Split explicitly on CRLF and LF rather than with str.splitlines, which also
# splits on \x0b, \x0c, \x1c-\x1e and \x85. A tag holding arbitrary bytes decoded
# permissively would then be cut in places that are not line breaks.
image_description_pairs = []
image_description_raw_texts = []
for description_value in image_description_tag_values:
    if isinstance(description_value, bytes):
        description_text = description_value.decode("latin-1")
    else:
        description_text = str(description_value)
    image_description_raw_texts.append(description_text)
    for description_line in description_text.replace("\r\n", "\n").split("\n"):
        stripped_description_line = description_line.strip().strip("\r")
        if stripped_description_line == "" or "=" not in stripped_description_line:
            continue
        description_key, description_value_text = stripped_description_line.split(
            "=", 1
        )
        image_description_pairs.append(
            [description_key.strip(), description_value_text.strip()]
        )

# A dict would silently drop one of these. The Typhoon writes Shading twice, once
# as "On" and once as the shading file path, so this is not hypothetical.
image_description_key_occurrence_counts = {}
for pair_key, _ in image_description_pairs:
    image_description_key_occurrence_counts[pair_key] = (
        image_description_key_occurrence_counts.get(pair_key, 0) + 1
    )
image_description_duplicate_keys = sorted(
    pair_key
    for pair_key, occurrence_count in image_description_key_occurrence_counts.items()
    if occurrence_count > 1
)
if len(image_description_duplicate_keys) > 0:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "image_description_duplicate_keys",
            "status": "reported",
            "is_hard_stop": False,
            "detail": "keys appearing more than once, all occurrences recorded in order: "
            + ", ".join(image_description_duplicate_keys),
        }
    )

image_description_first_value_by_key = {}
for pair_key, pair_value in image_description_pairs:
    image_description_first_value_by_key.setdefault(pair_key, pair_value)

# =============================================================================
# Cross-checks between sources inside the TIFF
# =============================================================================

# Pixel size stated in ImageDescription against the resolution tags. The Typhoon
# writes "Pixel size=200 micrometer"; the Amersham Imager 680 does not write the
# field at all, so its absence is not a failure.
description_pixel_size_text = image_description_first_value_by_key.get("Pixel size")
if description_pixel_size_text is not None and micrometres_per_pixel is not None:
    description_pixel_size_match = re.match(
        r"^\s*([0-9]+(?:\.[0-9]+)?)", description_pixel_size_text
    )
    if description_pixel_size_match is None:
        VALIDATION_FINDINGS.append(
            {
                "check_name": "pixel_size_description_matches_resolution_tags",
                "status": "warning",
                "is_hard_stop": False,
                "detail": "ImageDescription Pixel size "
                + repr(description_pixel_size_text)
                + " does not begin with a number, so it could not be compared",
            }
        )
    else:
        description_pixel_size_micrometres = float(
            description_pixel_size_match.group(1)
        )
        pixel_size_difference = abs(
            description_pixel_size_micrometres - micrometres_per_pixel
        )
        VALIDATION_FINDINGS.append(
            {
                "check_name": "pixel_size_description_matches_resolution_tags",
                "status": "pass" if pixel_size_difference <= 0.05 else "fail",
                "is_hard_stop": True,
                "detail": "ImageDescription says %.4f micrometres per pixel and the "
                "resolution tags say %.4f"
                % (description_pixel_size_micrometres, micrometres_per_pixel)
                + (
                    ""
                    if pixel_size_difference <= 0.05
                    else ". Two sources inside one file disagree about the scale, which "
                    "DESIGN.md section 11 makes a hard stop rather than a preference."
                ),
            }
        )

# Orientation. Tag 274 is absent from both instruments' files, so the TIFF's own
# statement lives in ImageDescription. DESIGN.md section 11 required orientation to
# come from the TIFF rather than the .inf and made disagreement a hard stop; this
# is the form that requirement actually takes on these files.
orientation_sources = {}
if "Orientation" in tag_values_by_name:
    orientation_sources["tiff_tag_274"] = tag_values_by_name["Orientation"][0][
        "value_repr"
    ]
if "Orientation" in image_description_first_value_by_key:
    orientation_sources["tiff_image_description"] = (
        image_description_first_value_by_key["Orientation"]
    )

# =============================================================================
# Pixel data and the single-pass statistics
# =============================================================================

# Every failure below is recorded and the run continues to the report, rather than
# calling die(). The tags have already parsed at this point, so there is always
# something worth writing: an 8-bit or RGB or lossily compressed file produces a
# report naming exactly why it is unquantifiable, which is more use than an exit
# code. Found by running the failure fixtures, which exited silently under the
# first version of this section.
pixel_array = None
pixel_statistics_were_computed = False

try:
    pixel_array = selected_tiff_page.asarray()
except Exception as pixel_read_error:
    # Deliberately broad: lossy or unsupported compression, a truncated strip and a
    # missing codec all land here, and the exception type is itself the diagnostic.
    VALIDATION_FINDINGS.append(
        {
            "check_name": "pixel_data_decodes",
            "status": "fail",
            "is_hard_stop": True,
            "detail": "page "
            + str(selected_page_index)
            + " would not decode: "
            + type(pixel_read_error).__name__
            + ": "
            + str(pixel_read_error),
        }
    )

if pixel_array is not None and pixel_array.ndim != 2:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "pixel_data_is_single_channel_two_dimensional",
            "status": "fail",
            "is_hard_stop": True,
            "detail": "decoded array has shape "
            + str(pixel_array.shape)
            + "; stage 1 handles "
            "a single-channel two-dimensional page only. SamplesPerPixel above 1 is "
            "already a hard stop above, and a 3 there means something rendered this "
            "file rather than the instrument writing it.",
        }
    )
    pixel_array = None

if pixel_array is not None and pixel_array.dtype != numpy.uint16:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "pixel_dtype_is_unsigned_16_bit",
            "status": "fail",
            "is_hard_stop": True,
            "detail": "decoded dtype is "
            + str(pixel_array.dtype)
            + " rather than uint16. The "
            "value-count statistics assume a 16-bit unsigned container, so they are "
            "omitted from this report rather than computed on a container that "
            "cannot hold the range they describe.",
        }
    )
    pixel_array = None

if pixel_array is not None:
    pixel_array_height_pixels, pixel_array_width_pixels = pixel_array.shape
    VALIDATION_FINDINGS.append(
        {
            "check_name": "tiff_tag_dimensions_match_decoded_array",
            "status": "pass"
            if (
                pixel_array_width_pixels == page_image_width_pixels
                and pixel_array_height_pixels == page_image_length_pixels
            )
            else "fail",
            "is_hard_stop": True,
            "detail": "tags say ImageWidth "
            + str(page_image_width_pixels)
            + " by ImageLength "
            + str(page_image_length_pixels)
            + "; the decoded array is "
            + str(pixel_array_width_pixels)
            + " by "
            + str(pixel_array_height_pixels),
        }
    )
else:
    # Sidecar geometry can still be checked against the tags, which parsed, so a
    # bad crop box is still caught on a file whose pixels are unreadable.
    pixel_array_width_pixels = page_image_width_pixels
    pixel_array_height_pixels = page_image_length_pixels
    emit_message(
        "pixels",
        "pixel statistics omitted; the report records why and the tag-based checks "
        "below still run against ImageWidth "
        + str(page_image_width_pixels)
        + " by ImageLength "
        + str(page_image_length_pixels),
    )

if pixel_array is not None:
    pixel_statistics_were_computed = True

# Bound to None first so that the report can record their absence explicitly rather
# than omitting keys, which would make a consumer unable to distinguish "not
# computed" from "this version of the report did not have the field".
pixel_value_counts = None
present_pixel_values = None
total_pixel_count = None
minimum_pixel_value = None
maximum_pixel_value = None
mean_pixel_value = None
median_pixel_value = None
standard_deviation_pixel_value = None
median_absolute_deviation_counts = None
robust_standard_deviation_counts = None
histogram_mode_value = None
pixels_above_mode = None
pixels_below_mode = None
inferred_signal_direction = None
distinct_pixel_value_count = None
pixel_value_stride = None
effective_bits_per_sample = None
inferred_ceiling_value = None
at_inferred_ceiling_pixel_count = None
at_container_maximum_pixel_count = None
at_maximum_observed_pixel_count = None
saturation_spike_ratio = None
at_floor_pixel_count = None
floor_population_fraction = None
at_minimum_observed_pixel_count = None
isolated_extreme_pixel_count = None
isolation_margin_counts = None

# Guarded rather than assumed: an 8-bit, RGB, lossily compressed or otherwise
# undecodable page produces a report naming why, with this whole block omitted,
# instead of exiting without one.
if pixel_statistics_were_computed:
    # One bincount pass yields every statistic DESIGN.md section 6 asks for: the
    # distinct-value count, the floor and ceiling populations, the stride that reveals
    # fewer real bits left-shifted into 16, and the histogram itself. Measured at 0.07 s
    # for 1125x875 and 0.30 s for 2048x2816, so there is no reason to make five passes,
    # and no reason to sort or take percentiles over millions of values.
    pixel_value_counts = numpy.bincount(
        pixel_array.ravel(), minlength=CONTAINER_VALUE_SLOT_COUNT_16_BIT
    ).astype(numpy.int64)
    total_pixel_count = int(pixel_value_counts.sum())
    present_pixel_values = numpy.nonzero(pixel_value_counts)[0]
    minimum_pixel_value = int(present_pixel_values[0])
    maximum_pixel_value = int(present_pixel_values[-1])
    distinct_pixel_value_count = int(present_pixel_values.size)

    # Mean, variance, median and MAD computed from the value counts rather than from
    # the array, so nothing sorts and nothing allocates a float64 copy of the image.
    present_values_float = present_pixel_values.astype(numpy.float64)
    present_counts_float = pixel_value_counts[present_pixel_values].astype(
        numpy.float64
    )
    mean_pixel_value = float(
        (present_values_float * present_counts_float).sum() / total_pixel_count
    )
    variance_pixel_value = float(
        (((present_values_float - mean_pixel_value) ** 2) * present_counts_float).sum()
        / total_pixel_count
    )
    standard_deviation_pixel_value = math.sqrt(variance_pixel_value)

    cumulative_pixel_counts = numpy.cumsum(pixel_value_counts)
    median_pixel_value = int(
        numpy.searchsorted(
            cumulative_pixel_counts, (total_pixel_count + 1) // 2, side="left"
        )
    )

    # Absolute deviation d from the median is contributed by the two values median-d
    # and median+d, so the MAD comes out of the same count array without touching the
    # image again.
    maximum_absolute_deviation = max(
        median_pixel_value - minimum_pixel_value,
        maximum_pixel_value - median_pixel_value,
    )
    deviation_magnitudes = numpy.arange(0, maximum_absolute_deviation + 1)
    lower_values = median_pixel_value - deviation_magnitudes
    upper_values = median_pixel_value + deviation_magnitudes
    lower_in_range = lower_values >= 0
    upper_in_range = upper_values < CONTAINER_VALUE_SLOT_COUNT_16_BIT
    absolute_deviation_counts = numpy.zeros(
        maximum_absolute_deviation + 1, dtype=numpy.int64
    )
    absolute_deviation_counts[lower_in_range] += pixel_value_counts[
        lower_values[lower_in_range]
    ]
    absolute_deviation_counts[upper_in_range] += pixel_value_counts[
        upper_values[upper_in_range]
    ]
    # Deviation zero was counted once from each side.
    absolute_deviation_counts[0] = pixel_value_counts[median_pixel_value]
    median_absolute_deviation_counts = int(
        numpy.searchsorted(
            numpy.cumsum(absolute_deviation_counts),
            (total_pixel_count + 1) // 2,
            side="left",
        )
    )
    robust_standard_deviation_counts = (
        MEDIAN_ABSOLUTE_DEVIATION_TO_SIGMA * median_absolute_deviation_counts
    )

    # Stride over the present values. A stride of 16 means 12 real bits left-shifted
    # into a 16-bit container, in which case the ceiling is 65520 rather than 65535 and
    # an equality test against 65535 reports zero saturation on a clipped image.
    pixel_value_stride = (
        int(numpy.gcd.reduce(present_pixel_values))
        if distinct_pixel_value_count > 1
        else 1
    )
    if pixel_value_stride < 1:
        pixel_value_stride = 1
    effective_bits_per_sample = STRIDE_TO_EFFECTIVE_BITS_PER_SAMPLE.get(
        pixel_value_stride
    )
    if effective_bits_per_sample is None:
        inferred_ceiling_value = CONTAINER_MAXIMUM_VALUE_16_BIT
        effective_bit_depth_detail = (
            "stride "
            + str(pixel_value_stride)
            + " is not a power of two consistent with "
            "left-shifted data; treating the container maximum as the ceiling"
        )
    else:
        inferred_ceiling_value = (
            CONTAINER_MAXIMUM_VALUE_16_BIT // pixel_value_stride
        ) * pixel_value_stride
        effective_bit_depth_detail = (
            "stride "
            + str(pixel_value_stride)
            + " implies "
            + str(effective_bits_per_sample)
            + " effective bits, so the container ceiling is "
            + str(inferred_ceiling_value)
        )

    at_container_maximum_pixel_count = int(
        pixel_value_counts[CONTAINER_MAXIMUM_VALUE_16_BIT]
    )
    at_inferred_ceiling_pixel_count = int(pixel_value_counts[inferred_ceiling_value])
    at_maximum_observed_pixel_count = int(pixel_value_counts[maximum_pixel_value])
    at_floor_pixel_count = int(pixel_value_counts[0])
    at_minimum_observed_pixel_count = int(pixel_value_counts[minimum_pixel_value])

    # Clipping presents as an anomalous spike at the maximum, not as the maximum being
    # any particular number, and the instrument can saturate below the numeric ceiling.
    # The comparison population is the present values immediately below the maximum.
    comparison_values = present_pixel_values[
        max(
            0,
            distinct_pixel_value_count - 1 - SATURATION_SPIKE_NEIGHBOURHOOD_VALUE_COUNT,
        ) : distinct_pixel_value_count - 1
    ]
    if comparison_values.size > 0:
        comparison_mean_count = float(pixel_value_counts[comparison_values].mean())
    else:
        comparison_mean_count = 0.0
    if comparison_mean_count > 0.0:
        saturation_spike_ratio = at_maximum_observed_pixel_count / comparison_mean_count
    else:
        saturation_spike_ratio = None

    VALIDATION_FINDINGS.append(
        {
            "check_name": "effective_bit_depth_and_ceiling",
            "status": "reported",
            "is_hard_stop": False,
            "detail": effective_bit_depth_detail,
        }
    )

    saturation_status = "pass"
    saturation_detail = (
        "maximum observed value "
        + str(maximum_pixel_value)
        + " held by "
        + str(at_maximum_observed_pixel_count)
        + " pixels; "
        + str(at_container_maximum_pixel_count)
        + " pixels at the container maximum "
        + str(CONTAINER_MAXIMUM_VALUE_16_BIT)
    )
    if saturation_spike_ratio is not None:
        saturation_detail += (
            "; spike ratio at the maximum is %.2f against the %d present values below it"
            % (saturation_spike_ratio, comparison_values.size)
        )
    if at_container_maximum_pixel_count > 0 or (
        saturation_spike_ratio is not None
        and saturation_spike_ratio > SATURATION_SPIKE_RATIO_WARNING_THRESHOLD
    ):
        saturation_status = "warning"
        saturation_detail += (
            ". This is evidence of clipping. Stage 1 reports it and does not adjudicate: "
            "per-band at-ceiling counts and the plateau statistic decide whether any "
            "particular measurement is affected, and the plateau statistic is the primary "
            "detector because the detector can saturate below the numeric ceiling."
        )
    VALIDATION_FINDINGS.append(
        {
            "check_name": "whole_image_saturation",
            "status": saturation_status,
            "is_hard_stop": False,
            "detail": saturation_detail,
        }
    )

    floor_population_fraction = at_floor_pixel_count / total_pixel_count
    VALIDATION_FINDINGS.append(
        {
            "check_name": "floor_population",
            "status": (
                "warning"
                if floor_population_fraction
                > FLOOR_POPULATION_FRACTION_WARNING_THRESHOLD
                else "pass"
            ),
            "is_hard_stop": False,
            "detail": str(at_floor_pixel_count)
            + " pixels at zero (%.6f of the image), minimum "
            "observed value %d held by %d pixels"
            % (
                floor_population_fraction,
                minimum_pixel_value,
                at_minimum_observed_pixel_count,
            )
            + (
                ". A large floor population means the image was already background "
                "subtracted or clipped at zero, which changes what a baseline means."
                if floor_population_fraction
                > FLOOR_POPULATION_FRACTION_WARNING_THRESHOLD
                else ""
            ),
        }
    )

    # Polarity, inferred from the data rather than from the tag, then cross-checked
    # against the tag. The mode is the background, because background is almost all of
    # a gel; the long tail runs toward signal.
    histogram_mode_value = int(numpy.argmax(pixel_value_counts))
    span_above_mode = maximum_pixel_value - histogram_mode_value
    span_below_mode = histogram_mode_value - minimum_pixel_value
    pixels_above_mode = int(pixel_value_counts[histogram_mode_value + 1 :].sum())
    pixels_below_mode = int(pixel_value_counts[:histogram_mode_value].sum())
    if span_above_mode > span_below_mode:
        inferred_signal_direction = "higher_value_is_more_signal"
    elif span_below_mode > span_above_mode:
        inferred_signal_direction = "lower_value_is_more_signal"
    else:
        inferred_signal_direction = "indeterminate"

    # MINISWHITE renders zero as white and the maximum as black, so under MINISWHITE a
    # signal direction of higher-is-more renders bands dark on a light field, which is
    # the conventional autoradiograph appearance and is what Fiji shows.
    photometric_is_miniswhite = int(page_photometric) == 0
    polarity_is_consistent = (
        photometric_is_miniswhite
        and inferred_signal_direction == "higher_value_is_more_signal"
    ) or (
        not photometric_is_miniswhite
        and inferred_signal_direction == "lower_value_is_more_signal"
    )
    VALIDATION_FINDINGS.append(
        {
            "check_name": "polarity_histogram_against_photometric_interpretation",
            "status": "pass" if polarity_is_consistent else "warning",
            "is_hard_stop": False,
            "detail": "histogram mode is "
            + str(histogram_mode_value)
            + " with "
            + str(pixels_above_mode)
            + " pixels above and "
            + str(pixels_below_mode)
            + " below; the tail runs "
            + str(span_above_mode)
            + " counts up and "
            + str(span_below_mode)
            + " counts down, so the inference is "
            + inferred_signal_direction
            + ", against PhotometricInterpretation "
            + page_photometric_name
            + (
                ". Consistent: bands render dark on a light field, which is what Fiji "
                "shows."
                if polarity_is_consistent
                else ". Inconsistent. Integrating with the wrong sign inverts every ratio, "
                "so resolve this before measuring."
            ),
        }
    )

    # Isolated extreme pixels: cosmic rays and screen contamination on a phosphor
    # plate, cosmic rays and hot pixels on a CCD cooled to -25 C through a 77 second
    # exposure. Reported, never filtered; DESIGN.md section 6 leaves the decision to
    # the user. The margin is in robust standard deviations rather than counts so that
    # it means the same thing on both instruments.
    isolation_margin_counts = (
        ISOLATION_MARGIN_ROBUST_SIGMAS * robust_standard_deviation_counts
    )
    interior_pixel_values = pixel_array[1:-1, 1:-1]
    neighbour_maximum_values = numpy.zeros_like(interior_pixel_values)
    for row_offset in (-1, 0, 1):
        for column_offset in (-1, 0, 1):
            if row_offset == 0 and column_offset == 0:
                continue
            shifted_view = pixel_array[
                1 + row_offset : pixel_array_height_pixels - 1 + row_offset,
                1 + column_offset : pixel_array_width_pixels - 1 + column_offset,
            ]
            numpy.maximum(
                neighbour_maximum_values, shifted_view, out=neighbour_maximum_values
            )
    isolated_extreme_pixel_count = int(
        (
            interior_pixel_values.astype(numpy.int32)
            - neighbour_maximum_values.astype(numpy.int32)
            > isolation_margin_counts
        ).sum()
    )
    VALIDATION_FINDINGS.append(
        {
            "check_name": "isolated_extreme_pixels",
            "status": "reported",
            "is_hard_stop": False,
            "detail": str(isolated_extreme_pixel_count)
            + " interior pixels exceed all eight "
            "neighbours by more than %.1f counts (%.1f robust standard deviations, "
            "MAD %d). Reported, not filtered."
            % (
                isolation_margin_counts,
                ISOLATION_MARGIN_ROBUST_SIGMAS,
                median_absolute_deviation_counts,
            ),
        }
    )

    emit_message(
        "statistics",
        "min "
        + str(minimum_pixel_value)
        + ", max "
        + str(maximum_pixel_value)
        + ", mean %.2f, median %d, sd %.2f, distinct %d, stride %d, isolated extremes %d"
        % (
            mean_pixel_value,
            median_pixel_value,
            standard_deviation_pixel_value,
            distinct_pixel_value_count,
            pixel_value_stride,
            isolated_extreme_pixel_count,
        ),
    )

# =============================================================================
# The .inf sidecar, as provenance only
# =============================================================================

# Demoted from a validation target to pure provenance, by decision. The .tif
# carries dimensions, bits per sample, pixel size, timestamp and orientation
# authoritatively, so nothing here is read for a decision. That retires three
# defects at once rather than fixing them: section 2's positional line numbers are
# off by one, section 12's byte arithmetic did not close (the omitted line is blank,
# not four characters, leaving 4 bytes unexplained), and section 6's "line 14"
# count check points at the wrong line on every real file.
inf_sidecar_absolute_path = input_tiff_absolute_path.with_suffix(INF_SIDECAR_SUFFIX)
inf_sidecar_record = {
    "path": str(inf_sidecar_absolute_path),
    "is_present": False,
    "size_bytes": None,
    "line_count": None,
    "lines_verbatim": None,
    "keyed_pairs": None,
    "declared_keyed_line_count_text": None,
    "non_ascii_byte_count": None,
}
if os.path.isfile(inf_sidecar_absolute_path):
    try:
        inf_sidecar_bytes = inf_sidecar_absolute_path.read_bytes()
    except OSError as inf_read_error:
        inf_sidecar_bytes = None
        VALIDATION_FINDINGS.append(
            {
                "check_name": "inf_sidecar_readable",
                "status": "warning",
                "is_hard_stop": False,
                "detail": "sidecar exists but could not be read: "
                + str(inf_read_error),
            }
        )
    if inf_sidecar_bytes is not None:
        inf_sidecar_record["is_present"] = True
        inf_sidecar_record["size_bytes"] = len(inf_sidecar_bytes)
        inf_sidecar_record["non_ascii_byte_count"] = int(
            sum(1 for single_byte in inf_sidecar_bytes if single_byte > 127)
        )
        # latin-1 never raises and maps bytes one to one, so the record is faithful
        # even if the file is not what it claims to be. json.dumps escapes anything
        # outside ASCII on the way out, satisfying CONVENTIONS.md section 4.
        inf_sidecar_text = inf_sidecar_bytes.decode("latin-1")
        # Blank lines are preserved deliberately: there is at least one, at real
        # position 3, and the whole point of a verbatim record is that it is
        # verbatim.
        inf_sidecar_lines = inf_sidecar_text.replace("\r\n", "\n").split("\n")
        inf_sidecar_lines = [
            single_line.rstrip("\r") for single_line in inf_sidecar_lines
        ]
        inf_sidecar_record["line_count"] = len(inf_sidecar_lines)
        inf_sidecar_record["lines_verbatim"] = inf_sidecar_lines

        inf_keyed_pairs = []
        first_keyed_line_index = None
        for line_index, single_line in enumerate(inf_sidecar_lines):
            keyed_match = INF_KEYED_LINE_PATTERN.match(single_line.strip())
            if keyed_match is None:
                continue
            if first_keyed_line_index is None:
                first_keyed_line_index = line_index
            inf_keyed_pairs.append(
                [
                    keyed_match.group(1),
                    keyed_match.group(2).strip(),
                    keyed_match.group(3).strip(),
                ]
            )
        inf_sidecar_record["keyed_pairs"] = inf_keyed_pairs

        # Located relative to the first keyed line rather than by absolute index,
        # because every positional index in section 2's table is off by one.
        if first_keyed_line_index is not None and first_keyed_line_index > 0:
            declared_count_text = inf_sidecar_lines[first_keyed_line_index - 1].strip()
            inf_sidecar_record["declared_keyed_line_count_text"] = declared_count_text
            declared_count_matches = declared_count_text.isdigit() and int(
                declared_count_text
            ) == len(inf_keyed_pairs)
            VALIDATION_FINDINGS.append(
                {
                    "check_name": "inf_declared_keyed_line_count",
                    "status": "reported" if declared_count_matches else "warning",
                    "is_hard_stop": False,
                    "detail": "the line before the first keyed line reads "
                    + repr(declared_count_text)
                    + " and there are "
                    + str(len(inf_keyed_pairs))
                    + " keyed lines"
                    + (
                        ""
                        if declared_count_matches
                        else ". Recorded as provenance only; nothing reads the .inf for "
                        "a decision."
                    ),
                }
            )

        first_line_text = (
            inf_sidecar_lines[0].strip() if len(inf_sidecar_lines) > 0 else ""
        )
        VALIDATION_FINDINGS.append(
            {
                "check_name": "inf_magic_string",
                "status": "reported"
                if first_line_text == "FLA_IMAGE_FILE"
                else "warning",
                "is_hard_stop": False,
                "detail": "first line is "
                + repr(first_line_text)
                + ". Compared after stripping, because the file is CRLF and an "
                "unstripped comparison fails on every real file.",
            }
        )

        # Both timestamps live in the positional block, so section 6's cross-check
        # would have been lost with the positional parser. Located by pattern
        # instead, which is immune to the off-by-one that made positional indexing
        # unsafe, so the check survives the demotion.
        inf_human_timestamp_text = None
        inf_epoch_timestamp_text = None
        for single_line in inf_sidecar_lines:
            stripped_line = single_line.strip()
            if inf_human_timestamp_text is None and INF_HUMAN_TIMESTAMP_PATTERN.match(
                stripped_line
            ):
                inf_human_timestamp_text = stripped_line
            if inf_epoch_timestamp_text is None and INF_EPOCH_TIMESTAMP_PATTERN.match(
                stripped_line
            ):
                inf_epoch_timestamp_text = stripped_line
        if (
            inf_human_timestamp_text is not None
            and inf_epoch_timestamp_text is not None
        ):
            try:
                human_timestamp = datetime.datetime.strptime(
                    inf_human_timestamp_text, "%a %b %d %H:%M:%S %Y"
                )
            except ValueError:
                human_timestamp = None
            if human_timestamp is not None:
                epoch_timestamp = datetime.datetime.fromtimestamp(
                    int(inf_epoch_timestamp_text), datetime.timezone.utc
                ).replace(tzinfo=None)
                # Minutes and seconds must match exactly and the hour difference must
                # be a whole number of hours. A fixed UTC offset must NOT be asserted:
                # this scan is UTC-5 but any scan between the second Sunday in March
                # and the first Sunday in November is UTC-4, so a hardcoded offset
                # makes this fire on correct files. The whole-hours form still catches
                # a corrupted header, which is all it was ever for.
                timestamp_difference_seconds = (
                    epoch_timestamp - human_timestamp
                ).total_seconds()
                whole_hours_apart = (
                    abs(timestamp_difference_seconds % 3600.0) < 1e-6
                    and human_timestamp.minute == epoch_timestamp.minute
                    and human_timestamp.second == epoch_timestamp.second
                )
                VALIDATION_FINDINGS.append(
                    {
                        "check_name": "inf_timestamps_agree_modulo_whole_hours",
                        "status": "reported" if whole_hours_apart else "warning",
                        "is_hard_stop": False,
                        "detail": inf_human_timestamp_text
                        + " against epoch "
                        + inf_epoch_timestamp_text
                        + " ("
                        + epoch_timestamp.strftime("%a %b %d %H:%M:%S %Y")
                        + " UTC), "
                        + "%.0f seconds apart, %.4f hours"
                        % (
                            timestamp_difference_seconds,
                            timestamp_difference_seconds / 3600.0,
                        )
                        + (
                            ""
                            if whole_hours_apart
                            else ". Not a whole number of hours with matching minutes and "
                            "seconds, which is what a corrupted header looks like."
                        ),
                    }
                )

        # Orientation from the .inf, recorded and cross-checked against the TIFF's
        # own statement. DESIGN.md section 11 makes disagreement a hard stop.
        for keyed_prefix, keyed_name, keyed_value in inf_keyed_pairs:
            if keyed_name == "Orientation":
                orientation_sources["inf_keyed"] = keyed_value
                break
else:
    # Not a failure. The read path is the .tif, which is self-describing, and the
    # fluorescence and loading scans from the Amersham Imager 680 never have a .inf.
    VALIDATION_FINDINGS.append(
        {
            "check_name": "inf_sidecar_present",
            "status": "reported",
            "is_hard_stop": False,
            "detail": "no .inf beside this file, which is expected for Amersham Imager 680 "
            "output. Nothing is lost: the .tif carries dimensions, bits per "
            "sample, pixel size, timestamp and orientation.",
        }
    )

distinct_orientation_statements = sorted(
    {
        orientation_text.strip().lower()
        for orientation_text in orientation_sources.values()
    }
)
if len(orientation_sources) == 0:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "orientation_sources_agree",
            "status": "warning",
            "is_hard_stop": False,
            "detail": "no orientation statement in tag 274, in ImageDescription, or in a "
            ".inf. Row order is then unverified and must be confirmed by eye on "
            "the first overlay before any band label is trusted.",
        }
    )
else:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "orientation_sources_agree",
            "status": "pass" if len(distinct_orientation_statements) == 1 else "fail",
            "is_hard_stop": True,
            "detail": "; ".join(
                source_name + " says " + repr(orientation_text)
                for source_name, orientation_text in sorted(orientation_sources.items())
            )
            + (
                ". Agreed."
                if len(distinct_orientation_statements) == 1
                else ". Sources disagree, which DESIGN.md section 11 makes a hard stop rather "
                "than a preference: getting row order wrong inverts band labels while "
                "leaving lane order correct, producing entirely plausible wrong numbers."
            ),
        }
    )

# =============================================================================
# The preprocessing sidecar, when present
# =============================================================================

preprocess_sidecar_absolute_path = input_tiff_absolute_path.parent / (
    input_tiff_absolute_path.stem + PREPROCESS_SIDECAR_FILENAME_SUFFIX
)
preprocess_sidecar_record = {
    "path": str(preprocess_sidecar_absolute_path),
    "is_present": False,
    "values": None,
    "derived_tilt_angle_degrees": None,
    "landmark_span_pixels": None,
    "angle_uncertainty_degrees_per_pixel_of_click_error": None,
    "implied_lane_pitch_millimetres": None,
}
parsed_sidecar_values = None
if os.path.isfile(preprocess_sidecar_absolute_path):
    preprocess_sidecar_record["is_present"] = True
    try:
        preprocess_sidecar_text = preprocess_sidecar_absolute_path.read_bytes().decode(
            "latin-1"
        )
    except OSError as sidecar_read_error:
        preprocess_sidecar_text = None
        VALIDATION_FINDINGS.append(
            {
                "check_name": "preprocess_sidecar_readable",
                "status": "fail",
                "is_hard_stop": True,
                "detail": "sidecar exists but could not be read: "
                + str(sidecar_read_error),
            }
        )
    if preprocess_sidecar_text is not None:
        raw_sidecar_values = {}
        sidecar_parse_problems = []
        for line_number, single_line in enumerate(
            preprocess_sidecar_text.replace("\r\n", "\n").split("\n"), start=1
        ):
            stripped_line = single_line.strip().strip("\r")
            if stripped_line == "" or stripped_line.startswith("#"):
                continue
            if "=" not in stripped_line:
                sidecar_parse_problems.append(
                    "line " + str(line_number) + " has no = : " + repr(stripped_line)
                )
                continue
            sidecar_key, sidecar_value = stripped_line.split("=", 1)
            sidecar_key = sidecar_key.strip()
            sidecar_value = sidecar_value.strip()
            if sidecar_key in raw_sidecar_values:
                sidecar_parse_problems.append(
                    "key " + repr(sidecar_key) + " appears more than once"
                )
            raw_sidecar_values[sidecar_key] = sidecar_value

        # An unknown key is a hard stop rather than a warning, and this is the
        # mechanism that keeps the template and this script from drifting apart
        # under the no-sharing rule: a key added to one and not the other fails on
        # the next run instead of being silently ignored. It also catches a
        # misspelling, which would otherwise leave a required value at its default.
        unknown_sidecar_keys = sorted(
            set(raw_sidecar_values) - set(REQUIRED_SIDECAR_KEY_NAMES)
        )
        missing_sidecar_keys = sorted(
            set(REQUIRED_SIDECAR_KEY_NAMES) - set(raw_sidecar_values)
        )
        for sidecar_key in REQUIRED_SIDECAR_KEY_NAMES:
            if (
                sidecar_key in raw_sidecar_values
                and raw_sidecar_values[sidecar_key] == ""
                and sidecar_key not in SIDECAR_KEY_NAMES_PERMITTED_EMPTY
            ):
                sidecar_parse_problems.append("key " + repr(sidecar_key) + " is empty")
        for sidecar_key in SIDECAR_INTEGER_KEY_NAMES:
            if sidecar_key not in raw_sidecar_values:
                continue
            sidecar_value = raw_sidecar_values[sidecar_key]
            if not re.match(r"^-?\d+$", sidecar_value):
                sidecar_parse_problems.append(
                    "key "
                    + repr(sidecar_key)
                    + " is not an integer: "
                    + repr(sidecar_value)
                )
        # Geometry values may carry a decimal point because centimetres are
        # fractional; a bare integer is also valid (pixel coordinates).
        for sidecar_key in SIDECAR_NUMERIC_KEY_NAMES:
            if sidecar_key not in raw_sidecar_values:
                continue
            sidecar_value = raw_sidecar_values[sidecar_key]
            if not re.match(r"^-?\d+(?:\.\d+)?$", sidecar_value):
                sidecar_parse_problems.append(
                    "key "
                    + repr(sidecar_key)
                    + " is not a number: "
                    + repr(sidecar_value)
                )
        # Enum values must be one of a closed set, so a typo like coordinate_unit=px
        # is a hard stop rather than a silently unrecognised value.
        for sidecar_key, allowed_values in SIDECAR_ENUM_KEY_ALLOWED_VALUES.items():
            if sidecar_key not in raw_sidecar_values:
                continue
            sidecar_value = raw_sidecar_values[sidecar_key]
            if sidecar_value not in allowed_values:
                sidecar_parse_problems.append(
                    "key "
                    + repr(sidecar_key)
                    + " is "
                    + repr(sidecar_value)
                    + "; allowed: "
                    + repr(allowed_values)
                )

        sidecar_structure_is_valid = (
            len(unknown_sidecar_keys) == 0
            and len(missing_sidecar_keys) == 0
            and len(sidecar_parse_problems) == 0
        )
        VALIDATION_FINDINGS.append(
            {
                "check_name": "preprocess_sidecar_structure",
                "status": "pass" if sidecar_structure_is_valid else "fail",
                "is_hard_stop": True,
                "detail": (
                    "all " + str(len(REQUIRED_SIDECAR_KEY_NAMES)) + " keys present and "
                    "well formed"
                    if sidecar_structure_is_valid
                    else "unknown keys: "
                    + repr(unknown_sidecar_keys)
                    + "; missing keys: "
                    + repr(missing_sidecar_keys)
                    + "; other problems: "
                    + repr(sidecar_parse_problems)
                    + ". preprocess_sidecar_template.txt is the authority for the key "
                    "list."
                ),
            }
        )

        if sidecar_structure_is_valid:
            parsed_sidecar_values = dict(raw_sidecar_values)
            for sidecar_key in SIDECAR_INTEGER_KEY_NAMES:
                parsed_sidecar_values[sidecar_key] = int(
                    parsed_sidecar_values[sidecar_key]
                )
            for sidecar_key in SIDECAR_NUMERIC_KEY_NAMES:
                parsed_sidecar_values[sidecar_key] = float(
                    parsed_sidecar_values[sidecar_key]
                )
            preprocess_sidecar_record["values"] = parsed_sidecar_values

            gel_migration_axis = parsed_sidecar_values["gel_migration_axis"]
            coordinate_unit = parsed_sidecar_values["coordinate_unit"]
            # Convert every geometry value to pixels once, so every check below works
            # in pixel space regardless of the unit the operator recorded (F13). ImageJ
            # reports calibrated centimetres when the image carries a spatial scale, and
            # a centimetre value is a perfectly valid-looking pixel value, so the unit is
            # declared explicitly rather than guessed.
            if coordinate_unit == "pixels":
                pixels_per_coordinate_unit = 1.0
            elif micrometres_per_pixel is not None:
                # 1 cm = 10000 micrometres, so pixels per cm = 10000 / (micrometres/px).
                pixels_per_coordinate_unit = 10000.0 / micrometres_per_pixel
            else:
                pixels_per_coordinate_unit = None
            coordinate_conversion_ok = pixels_per_coordinate_unit is not None
            VALIDATION_FINDINGS.append(
                {
                    "check_name": "preprocess_sidecar_coordinate_unit_convertible",
                    "status": "pass" if coordinate_conversion_ok else "fail",
                    "is_hard_stop": True,
                    "detail": (
                        (
                            "coordinate_unit 'pixels' needs no conversion"
                            if coordinate_unit == "pixels"
                            else "coordinate_unit 'centimetres' converts at %.6f pixels "
                            "per centimetre from the file's pixel size"
                            % pixels_per_coordinate_unit
                        )
                        if coordinate_conversion_ok
                        else "coordinate_unit is 'centimetres' but the file carries no usable "
                        "pixel size, so the coordinates cannot be converted to pixels"
                    ),
                }
            )

            if coordinate_conversion_ok:
                landmark_a_x_pixels = (
                    parsed_sidecar_values["landmark_a_x"] * pixels_per_coordinate_unit
                )
                landmark_a_y_pixels = (
                    parsed_sidecar_values["landmark_a_y"] * pixels_per_coordinate_unit
                )
                landmark_b_x_pixels = (
                    parsed_sidecar_values["landmark_b_x"] * pixels_per_coordinate_unit
                )
                landmark_b_y_pixels = (
                    parsed_sidecar_values["landmark_b_y"] * pixels_per_coordinate_unit
                )
                # Crop bounds must be whole pixels for slicing and the preview rectangle.
                crop_x_pixels = int(
                    round(parsed_sidecar_values["crop_x"] * pixels_per_coordinate_unit)
                )
                crop_y_pixels = int(
                    round(parsed_sidecar_values["crop_y"] * pixels_per_coordinate_unit)
                )
                crop_width_pixels = int(
                    round(
                        parsed_sidecar_values["crop_width"] * pixels_per_coordinate_unit
                    )
                )
                crop_height_pixels = int(
                    round(
                        parsed_sidecar_values["crop_height"]
                        * pixels_per_coordinate_unit
                    )
                )
                preprocess_sidecar_record["geometry_pixels"] = {
                    "gel_migration_axis": gel_migration_axis,
                    "coordinate_unit": coordinate_unit,
                    "pixels_per_coordinate_unit": pixels_per_coordinate_unit,
                    "landmark_a_x_pixels": landmark_a_x_pixels,
                    "landmark_a_y_pixels": landmark_a_y_pixels,
                    "landmark_b_x_pixels": landmark_b_x_pixels,
                    "landmark_b_y_pixels": landmark_b_y_pixels,
                    "crop_x_pixels": crop_x_pixels,
                    "crop_y_pixels": crop_y_pixels,
                    "crop_width_pixels": crop_width_pixels,
                    "crop_height_pixels": crop_height_pixels,
                }

            VALIDATION_FINDINGS.append(
                {
                    "check_name": "preprocess_sidecar_schema_version",
                    "status": (
                        "pass"
                        if parsed_sidecar_values["schema_version"]
                        == SUPPORTED_SIDECAR_SCHEMA_VERSION
                        else "fail"
                    ),
                    "is_hard_stop": True,
                    "detail": "sidecar declares schema_version "
                    + str(parsed_sidecar_values["schema_version"])
                    + "; this script "
                    "understands " + str(SUPPORTED_SIDECAR_SCHEMA_VERSION),
                }
            )

            frame_name_matches = (
                parsed_sidecar_values["measured_in_frame"]
                == REQUIRED_SIDECAR_FRAME_NAME
            )
            VALIDATION_FINDINGS.append(
                {
                    "check_name": "preprocess_sidecar_measurement_frame",
                    "status": "pass" if frame_name_matches else "fail",
                    "is_hard_stop": True,
                    "detail": "measured_in_frame is "
                    + repr(parsed_sidecar_values["measured_in_frame"])
                    + "; required "
                    + repr(REQUIRED_SIDECAR_FRAME_NAME)
                    + (
                        ""
                        if frame_name_matches
                        else ". Coordinates are meaningless without the frame they were "
                        "measured in, and a crop measured on a flipped view lands on "
                        "the wrong bands while every number stays plausible."
                    ),
                }
            )

            # The filename check exists because the dimension check alone does not
            # catch the likeliest hand-editing error: copying one scan's sidecar to
            # another, updating the dimensions because they visibly differ, and
            # leaving the crop coordinates. Repeat 2 is 1250x750 where repeats 1 and
            # 3 are 1125x875, so that error is live rather than hypothetical.
            filename_matches = (
                parsed_sidecar_values["measured_against_input_filename"]
                == input_tiff_absolute_path.name
            )
            VALIDATION_FINDINGS.append(
                {
                    "check_name": "preprocess_sidecar_filename_matches_input",
                    "status": "pass" if filename_matches else "fail",
                    "is_hard_stop": True,
                    "detail": "sidecar was measured against "
                    + repr(parsed_sidecar_values["measured_against_input_filename"])
                    + " and the input is "
                    + repr(input_tiff_absolute_path.name),
                }
            )

            dimensions_match = (
                parsed_sidecar_values["measured_against_image_width_pixels"]
                == pixel_array_width_pixels
                and parsed_sidecar_values["measured_against_image_height_pixels"]
                == pixel_array_height_pixels
            )
            VALIDATION_FINDINGS.append(
                {
                    "check_name": "preprocess_sidecar_dimensions_match_input",
                    "status": "pass" if dimensions_match else "fail",
                    "is_hard_stop": True,
                    "detail": "sidecar was measured against "
                    + str(parsed_sidecar_values["measured_against_image_width_pixels"])
                    + " by "
                    + str(parsed_sidecar_values["measured_against_image_height_pixels"])
                    + " and the input is "
                    + str(pixel_array_width_pixels)
                    + " by "
                    + str(pixel_array_height_pixels),
                }
            )

            # Everything below works in pixel space, which only exists once the
            # coordinate conversion above succeeded. If it did not, the convertible
            # check already recorded a hard stop and there is nothing to check here.
            if coordinate_conversion_ok:
                landmarks_inside_image = all(
                    0 <= landmark_x < pixel_array_width_pixels
                    and 0 <= landmark_y < pixel_array_height_pixels
                    for landmark_x, landmark_y in (
                        (landmark_a_x_pixels, landmark_a_y_pixels),
                        (landmark_b_x_pixels, landmark_b_y_pixels),
                    )
                )
                VALIDATION_FINDINGS.append(
                    {
                        "check_name": "preprocess_sidecar_landmarks_inside_image",
                        "status": "pass" if landmarks_inside_image else "fail",
                        "is_hard_stop": True,
                        "detail": "landmarks a (%.1f, %.1f) and b (%.1f, %.1f) against an image "
                        "of %d by %d"
                        % (
                            landmark_a_x_pixels,
                            landmark_a_y_pixels,
                            landmark_b_x_pixels,
                            landmark_b_y_pixels,
                            pixel_array_width_pixels,
                            pixel_array_height_pixels,
                        ),
                    }
                )

                # Delta runs a -> b. The landmark line is the level feature; how its
                # angle becomes a tilt depends on gel_migration_axis, derived below.
                landmark_delta_x = landmark_b_x_pixels - landmark_a_x_pixels
                landmark_delta_y = landmark_b_y_pixels - landmark_a_y_pixels
                landmark_span_pixels = math.hypot(landmark_delta_x, landmark_delta_y)
                preprocess_sidecar_record["landmark_span_pixels"] = landmark_span_pixels

                # Sign convention, stated once here so that nothing downstream has to
                # re-derive it: the landmarks are in Fiji's frame, where y increases
                # downward. Stage 1 does not rotate; it reports. Whatever applies the
                # rotation must transform these two points by the same rotation and
                # assert the re-derived angle is near zero, because both Fiji's
                # downward y and the rotation library's convention invert the sign,
                # and the failure mode is a doubled tilt that still looks roughly
                # straight. ImageJ's Angle is the independent cross-check on the sign.
                if landmark_span_pixels >= MINIMUM_LANDMARK_SPAN_PIXELS:
                    # The landmark line is meant to be level and perpendicular to
                    # migration. For a vertical gel that level line runs left-to-right,
                    # so its target angle from the x axis is 0; for a horizontal gel
                    # (migration left-to-right) the level line runs top-to-bottom, so
                    # its target is 90. The tilt is the signed deviation of the actual
                    # line from that target, folded onto the interval nearest zero so
                    # that reversed landmarks or a near-vertical line do not read as a
                    # 180 degree tilt. A perfectly level gel therefore reads about zero
                    # in either orientation, and the old atan2(dy, dx) that hard-stopped
                    # on a level horizontal gel is gone.
                    landmark_line_angle_degrees = math.degrees(
                        math.atan2(landmark_delta_y, landmark_delta_x)
                    )
                    target_line_angle_degrees = (
                        0.0 if gel_migration_axis == "vertical" else 90.0
                    )
                    derived_tilt_angle_degrees = (
                        landmark_line_angle_degrees - target_line_angle_degrees
                    )
                    while derived_tilt_angle_degrees > 90.0:
                        derived_tilt_angle_degrees -= 180.0
                    while derived_tilt_angle_degrees <= -90.0:
                        derived_tilt_angle_degrees += 180.0
                    preprocess_sidecar_record["derived_tilt_angle_degrees"] = (
                        derived_tilt_angle_degrees
                    )
                    angle_uncertainty_degrees = math.degrees(
                        math.atan(1.0 / landmark_span_pixels)
                    )
                    preprocess_sidecar_record[
                        "angle_uncertainty_degrees_per_pixel_of_click_error"
                    ] = angle_uncertainty_degrees
                    VALIDATION_FINDINGS.append(
                        {
                            "check_name": "preprocess_sidecar_tilt_angle_is_plausible",
                            "status": (
                                "pass"
                                if abs(derived_tilt_angle_degrees)
                                <= MAXIMUM_PLAUSIBLE_TILT_DEGREES
                                else "fail"
                            ),
                            "is_hard_stop": True,
                            "detail": "derived tilt %.4f degrees over a %.1f pixel span, "
                            "uncertainty %.4f degrees per pixel of click error"
                            % (
                                derived_tilt_angle_degrees,
                                landmark_span_pixels,
                                angle_uncertainty_degrees,
                            )
                            + (
                                ""
                                if abs(derived_tilt_angle_degrees)
                                <= MAXIMUM_PLAUSIBLE_TILT_DEGREES
                                else ". Beyond %.1f degrees the landmarks were almost "
                                "certainly not the level feature intended, or the two "
                                "pairs were entered in the wrong order."
                                % MAXIMUM_PLAUSIBLE_TILT_DEGREES
                            ),
                        }
                    )
                    if landmark_span_pixels < LANDMARK_SPAN_WARNING_PIXELS:
                        VALIDATION_FINDINGS.append(
                            {
                                "check_name": "preprocess_sidecar_landmark_span_is_generous",
                                "status": "warning",
                                "is_hard_stop": False,
                                "detail": "span is only %.1f pixels, giving %.3f degrees of "
                                "uncertainty per pixel of click error. Re-measure with "
                                "the landmarks further apart: 0.5 degrees of residual "
                                "tilt walks a band a full band height across a 1125 "
                                "pixel image."
                                % (landmark_span_pixels, angle_uncertainty_degrees),
                            }
                        )
                else:
                    VALIDATION_FINDINGS.append(
                        {
                            "check_name": "preprocess_sidecar_tilt_angle_is_plausible",
                            "status": "fail",
                            "is_hard_stop": True,
                            "detail": "landmark span is only %.1f pixels, below the floor of %d, "
                            "so no angle worth having can be derived from it"
                            % (landmark_span_pixels, MINIMUM_LANDMARK_SPAN_PIXELS),
                        }
                    )

                crop_x = crop_x_pixels
                crop_y = crop_y_pixels
                crop_width = crop_width_pixels
                crop_height = crop_height_pixels
                crop_is_inside_image = (
                    crop_width > 0
                    and crop_height > 0
                    and 0 <= crop_x
                    and 0 <= crop_y
                    and crop_x + crop_width <= pixel_array_width_pixels
                    and crop_y + crop_height <= pixel_array_height_pixels
                )
                VALIDATION_FINDINGS.append(
                    {
                        "check_name": "preprocess_sidecar_crop_inside_image",
                        "status": "pass" if crop_is_inside_image else "fail",
                        "is_hard_stop": True,
                        "detail": "crop x %d y %d width %d height %d against an image of %d by %d"
                        % (
                            crop_x,
                            crop_y,
                            crop_width,
                            crop_height,
                            pixel_array_width_pixels,
                            pixel_array_height_pixels,
                        ),
                    }
                )

                if crop_is_inside_image:
                    crop_area_fraction = (crop_width * crop_height) / float(
                        pixel_array_width_pixels * pixel_array_height_pixels
                    )
                    crop_fraction_is_plausible = (
                        CROP_AREA_FRACTION_WARNING_FLOOR
                        <= crop_area_fraction
                        <= CROP_AREA_FRACTION_WARNING_CEILING
                    )
                    VALIDATION_FINDINGS.append(
                        {
                            "check_name": "preprocess_sidecar_crop_area_fraction",
                            "status": "pass"
                            if crop_fraction_is_plausible
                            else "warning",
                            "is_hard_stop": False,
                            "detail": "the crop covers %.4f of the image"
                            % crop_area_fraction
                            + (
                                ""
                                if crop_fraction_is_plausible
                                else ". Outside the plausible range. A crop that is too tight "
                                "biases every baseline in the image, and one that covers "
                                "everything has excluded neither the wells nor the plate "
                                "background."
                            ),
                        }
                    )

                    # Cross-check the geometry against the pixel size read from this
                    # file's own tags. DESIGN.md section 2 records a 5 mm lane pitch,
                    # which is 25 px at 200 micrometres and 63 px at 79.9, so the same
                    # millimetre figure has to be reached through the file. Lanes are
                    # stacked across the dimension perpendicular to migration, so the
                    # pitch is taken along the crop height for a horizontally migrating
                    # gel and the crop width for a vertical one.
                    expected_lane_count = parsed_sidecar_values["expected_lane_count"]
                    if gel_migration_axis == "horizontal":
                        lane_stacking_extent_pixels = crop_height
                    else:
                        lane_stacking_extent_pixels = crop_width
                    if expected_lane_count <= 0:
                        VALIDATION_FINDINGS.append(
                            {
                                "check_name": "preprocess_sidecar_expected_lane_count",
                                "status": "fail",
                                "is_hard_stop": True,
                                "detail": "expected_lane_count is "
                                + str(expected_lane_count)
                                + ", which must be a positive integer",
                            }
                        )
                    elif micrometres_per_pixel is not None:
                        implied_lane_pitch_millimetres = (
                            (lane_stacking_extent_pixels / float(expected_lane_count))
                            * micrometres_per_pixel
                            / 1000.0
                        )
                        preprocess_sidecar_record["implied_lane_pitch_millimetres"] = (
                            implied_lane_pitch_millimetres
                        )
                        pitch_is_plausible = (
                            MINIMUM_PLAUSIBLE_LANE_PITCH_MILLIMETRES
                            <= implied_lane_pitch_millimetres
                            <= MAXIMUM_PLAUSIBLE_LANE_PITCH_MILLIMETRES
                        )
                        VALIDATION_FINDINGS.append(
                            {
                                "check_name": "preprocess_sidecar_implied_lane_pitch",
                                "status": "pass" if pitch_is_plausible else "warning",
                                "is_hard_stop": False,
                                "detail": "%d lanes across %d crop pixels (the %s dimension) at "
                                "%.4f micrometres per pixel implies a lane pitch of "
                                "%.2f mm"
                                % (
                                    expected_lane_count,
                                    lane_stacking_extent_pixels,
                                    "height"
                                    if gel_migration_axis == "horizontal"
                                    else "width",
                                    micrometres_per_pixel,
                                    implied_lane_pitch_millimetres,
                                )
                                + (
                                    ""
                                    if pitch_is_plausible
                                    else ". Outside the plausible range of %.1f to %.1f mm, so "
                                    "either the lane count or the crop is wrong."
                                    % (
                                        MINIMUM_PLAUSIBLE_LANE_PITCH_MILLIMETRES,
                                        MAXIMUM_PLAUSIBLE_LANE_PITCH_MILLIMETRES,
                                    )
                                ),
                            }
                        )

                    # In-crop statistics, computed the same single-pass way as the
                    # whole-image block above but over the crop only. The whole-image
                    # mode is the bare plate, which is the majority population, so a
                    # baseline or an isolation margin derived from it is off by hundreds
                    # of counts (F2, F16). This reports the gel's own numbers. Note that
                    # a generous crop still includes a strip of plate at the well end, so
                    # these describe the gel only where the crop is gel; the baseline
                    # proper is a per-band local estimate in stage 2, not a crop-wide one.
                    if pixel_statistics_were_computed:
                        crop_pixels = pixel_array[
                            crop_y : crop_y + crop_height, crop_x : crop_x + crop_width
                        ]
                        crop_value_counts = numpy.bincount(
                            crop_pixels.ravel(),
                            minlength=CONTAINER_VALUE_SLOT_COUNT_16_BIT,
                        ).astype(numpy.int64)
                        crop_total_pixel_count = int(crop_value_counts.sum())
                        crop_present_values = numpy.nonzero(crop_value_counts)[0]
                        crop_minimum_value = int(crop_present_values[0])
                        crop_maximum_value = int(crop_present_values[-1])
                        crop_mode_value = int(crop_value_counts.argmax())
                        crop_present_float = crop_present_values.astype(numpy.float64)
                        crop_counts_float = crop_value_counts[
                            crop_present_values
                        ].astype(numpy.float64)
                        crop_mean_value = float(
                            (crop_present_float * crop_counts_float).sum()
                            / crop_total_pixel_count
                        )
                        crop_cumulative_counts = numpy.cumsum(crop_value_counts)
                        crop_median_value = int(
                            numpy.searchsorted(
                                crop_cumulative_counts,
                                (crop_total_pixel_count + 1) // 2,
                                side="left",
                            )
                        )
                        # MAD from the value counts, same trick as the whole-image block:
                        # deviation d is contributed by median-d and median+d.
                        crop_maximum_absolute_deviation = max(
                            crop_median_value - crop_minimum_value,
                            crop_maximum_value - crop_median_value,
                        )
                        crop_deviation_magnitudes = numpy.arange(
                            0, crop_maximum_absolute_deviation + 1
                        )
                        crop_lower_values = (
                            crop_median_value - crop_deviation_magnitudes
                        )
                        crop_upper_values = (
                            crop_median_value + crop_deviation_magnitudes
                        )
                        crop_lower_in_range = crop_lower_values >= 0
                        crop_upper_in_range = (
                            crop_upper_values < CONTAINER_VALUE_SLOT_COUNT_16_BIT
                        )
                        crop_absolute_deviation_counts = numpy.zeros(
                            crop_maximum_absolute_deviation + 1, dtype=numpy.int64
                        )
                        crop_absolute_deviation_counts[crop_lower_in_range] += (
                            crop_value_counts[crop_lower_values[crop_lower_in_range]]
                        )
                        crop_absolute_deviation_counts[crop_upper_in_range] += (
                            crop_value_counts[crop_upper_values[crop_upper_in_range]]
                        )
                        crop_absolute_deviation_counts[0] = crop_value_counts[
                            crop_median_value
                        ]
                        crop_median_absolute_deviation = int(
                            numpy.searchsorted(
                                numpy.cumsum(crop_absolute_deviation_counts),
                                (crop_total_pixel_count + 1) // 2,
                                side="left",
                            )
                        )
                        crop_robust_standard_deviation = (
                            MEDIAN_ABSOLUTE_DEVIATION_TO_SIGMA
                            * crop_median_absolute_deviation
                        )
                        crop_isolation_margin_counts = (
                            ISOLATION_MARGIN_ROBUST_SIGMAS
                            * crop_robust_standard_deviation
                        )
                        whole_image_mode_value = int(pixel_value_counts.argmax())
                        preprocess_sidecar_record["in_crop_statistics"] = {
                            "crop_pixel_count": crop_total_pixel_count,
                            "minimum_pixel_value": crop_minimum_value,
                            "maximum_pixel_value": crop_maximum_value,
                            "mode_value": crop_mode_value,
                            "median_pixel_value": crop_median_value,
                            "mean_pixel_value": crop_mean_value,
                            "median_absolute_deviation_counts": crop_median_absolute_deviation,
                            "robust_standard_deviation_counts": crop_robust_standard_deviation,
                            "isolation_margin_counts": crop_isolation_margin_counts,
                            "whole_image_mode_value_for_comparison": whole_image_mode_value,
                            "whole_image_median_value_for_comparison": median_pixel_value,
                        }
                        VALIDATION_FINDINGS.append(
                            {
                                "check_name": "preprocess_sidecar_in_crop_statistics",
                                "status": "reported",
                                "is_hard_stop": False,
                                "detail": "in-crop mode %d, median %d, MAD %d, isolation margin "
                                "%.1f counts over %d pixels; whole-image mode %d and "
                                "median %d for comparison. The whole-image mode is the "
                                "bare plate, so the baseline and the isolation margin "
                                "belong on these in-crop numbers, not the whole-image "
                                "ones."
                                % (
                                    crop_mode_value,
                                    crop_median_value,
                                    crop_median_absolute_deviation,
                                    crop_isolation_margin_counts,
                                    crop_total_pixel_count,
                                    whole_image_mode_value,
                                    median_pixel_value,
                                ),
                            }
                        )
else:
    VALIDATION_FINDINGS.append(
        {
            "check_name": "preprocess_sidecar_present",
            "status": "reported",
            "is_hard_stop": False,
            "detail": "no sidecar at "
            + str(preprocess_sidecar_absolute_path)
            + ". Statistics and the tag inventory are unaffected; there is no crop "
            "preview and no tilt angle. See PROTOCOL.md section 4.",
        }
    )

# =============================================================================
# Histogram
# =============================================================================

histogram_was_written = False
if not pixel_statistics_were_computed:
    emit_message(
        "histogram",
        "skipped: no pixel statistics were computed, so there is no distribution to plot. "
        "The report records which check made the pixels unquantifiable.",
    )
else:
    # Log-count y axis, because band pixels are a tiny minority of a gel and a
    # linear y axis hides them completely.
    histogram_figure, histogram_axes = matplotlib.pyplot.subplots(
        figsize=HISTOGRAM_FIGURE_SIZE_INCHES
    )
    histogram_axes.fill_between(
        present_pixel_values,
        1,
        pixel_value_counts[present_pixel_values],
        step="mid",
        linewidth=0,
    )
    histogram_axes.set_yscale("log")
    histogram_axes.set_xlabel("pixel value, counts")
    histogram_axes.set_ylabel("number of pixels, log scale")
    histogram_axes.set_title(
        "Value distribution: "
        + input_tiff_absolute_path.name
        + "\n"
        + str(pixel_array_width_pixels)
        + " by "
        + str(pixel_array_height_pixels)
        + " pixels, "
        + (
            "%.4f micrometres per pixel" % micrometres_per_pixel
            if micrometres_per_pixel is not None
            else "pixel size unavailable"
        )
    )
    for marker_value, marker_label, marker_style in (
        (minimum_pixel_value, "minimum " + str(minimum_pixel_value), "dotted"),
        (median_pixel_value, "median " + str(median_pixel_value), "dashed"),
        (histogram_mode_value, "mode " + str(histogram_mode_value), "dashdot"),
        (maximum_pixel_value, "maximum " + str(maximum_pixel_value), "solid"),
    ):
        histogram_axes.axvline(
            marker_value, linestyle=marker_style, linewidth=1.0, label=marker_label
        )
    # Limited to the observed range rather than the container range. Background sits
    # near 2,800 counts of a possible 65,535, so an axis spanning the container
    # crushes the entire distribution into a sliver a few pixels wide and the plot
    # stops doing the one job it has. Found by looking at the output.
    observed_value_margin = max(1.0, 0.02 * (maximum_pixel_value - minimum_pixel_value))
    histogram_axes.set_xlim(
        minimum_pixel_value - observed_value_margin,
        maximum_pixel_value + observed_value_margin,
    )
    if minimum_pixel_value <= inferred_ceiling_value <= maximum_pixel_value:
        histogram_axes.axvline(
            inferred_ceiling_value,
            linestyle=(0, (1, 3)),
            linewidth=1.0,
            label="inferred ceiling " + str(inferred_ceiling_value),
        )
    else:
        # Stated rather than drawn, because drawing it is what stretched the axis.
        histogram_axes.plot(
            [],
            [],
            linestyle=(0, (1, 3)),
            linewidth=1.0,
            label="inferred ceiling " + str(inferred_ceiling_value) + ", off scale",
        )
    histogram_axes.legend(loc="upper right", fontsize="small")
    histogram_axes.grid(True, which="both", linewidth=0.3, alpha=0.5)
    histogram_figure.tight_layout()
    try:
        histogram_figure.savefig(histogram_image_output_path, dpi=FIGURE_DOTS_PER_INCH)
    except OSError as histogram_write_error:
        die(
            "histogram",
            "could not write "
            + str(histogram_image_output_path)
            + ": "
            + str(histogram_write_error),
        )
    matplotlib.pyplot.close(histogram_figure)
    histogram_was_written = True
    emit_message("histogram", "wrote " + str(histogram_image_output_path))

# =============================================================================
# Crop preview, only when the sidecar parsed
# =============================================================================

crop_preview_was_written = False
if (
    parsed_sidecar_values is not None
    and pixel_statistics_were_computed
    and "geometry_pixels" in preprocess_sidecar_record
):
    preview_geometry = preprocess_sidecar_record["geometry_pixels"]
    # Downsampled by an integer stride: this exists to be looked at, not measured.
    preview_stride = max(
        1,
        int(
            math.ceil(
                max(pixel_array_height_pixels, pixel_array_width_pixels)
                / float(PREVIEW_MAXIMUM_DIMENSION_PIXELS)
            )
        ),
    )
    preview_array = pixel_array[::preview_stride, ::preview_stride]
    preview_figure, preview_axes = matplotlib.pyplot.subplots(
        figsize=PREVIEW_FIGURE_SIZE_INCHES
    )
    # gray_r, not gray. PhotometricInterpretation is MINISWHITE, so Fiji renders
    # these files with an inverting lookup table and bands appear dark. A default
    # colormap would show bright bands, and anyone comparing this against Fiji would
    # conclude the pipeline inverted the image and start flipping things to fix it.
    preview_axes.imshow(
        preview_array,
        cmap=DISPLAY_COLORMAP_NAME,
        vmin=numpy.percentile(preview_array, 1.0),
        vmax=numpy.percentile(preview_array, 99.9),
        interpolation="nearest",
        origin="upper",
    )
    preview_axes.add_patch(
        matplotlib.patches.Rectangle(
            (
                preview_geometry["crop_x_pixels"] / preview_stride - 0.5,
                preview_geometry["crop_y_pixels"] / preview_stride - 0.5,
            ),
            preview_geometry["crop_width_pixels"] / preview_stride,
            preview_geometry["crop_height_pixels"] / preview_stride,
            fill=False,
            linewidth=1.6,
            linestyle="solid",
            edgecolor="red",
        )
    )
    preview_axes.plot(
        [
            preview_geometry["landmark_a_x_pixels"] / preview_stride,
            preview_geometry["landmark_b_x_pixels"] / preview_stride,
        ],
        [
            preview_geometry["landmark_a_y_pixels"] / preview_stride,
            preview_geometry["landmark_b_y_pixels"] / preview_stride,
        ],
        marker="+",
        markersize=11,
        linewidth=1.2,
        color="blue",
    )
    preview_title_text = (
        "Crop and landmarks as measured, unflipped: "
        + input_tiff_absolute_path.name
        + "\nrendered with "
        + DISPLAY_COLORMAP_NAME
        + " to match Fiji's inverting LUT; migration axis "
        + preview_geometry["gel_migration_axis"]
    )
    if preprocess_sidecar_record["derived_tilt_angle_degrees"] is not None:
        preview_title_text += "\nderived tilt %.4f degrees over a %.0f pixel span" % (
            preprocess_sidecar_record["derived_tilt_angle_degrees"],
            preprocess_sidecar_record["landmark_span_pixels"],
        )
    preview_axes.set_title(preview_title_text, fontsize="medium")
    preview_axes.set_xlabel("column / " + str(preview_stride) + " (Fiji x)")
    preview_axes.set_ylabel(
        "row / " + str(preview_stride) + " (Fiji y, increasing downward)"
    )
    preview_figure.tight_layout()
    try:
        preview_figure.savefig(crop_preview_image_output_path, dpi=FIGURE_DOTS_PER_INCH)
        crop_preview_was_written = True
    except OSError as preview_write_error:
        die(
            "crop preview",
            "could not write "
            + str(crop_preview_image_output_path)
            + ": "
            + str(preview_write_error),
        )
    matplotlib.pyplot.close(preview_figure)
    emit_message(
        "crop preview",
        "wrote "
        + str(crop_preview_image_output_path)
        + " at stride "
        + str(preview_stride),
    )

# =============================================================================
# Report
# =============================================================================

failed_hard_stop_check_names = [
    finding["check_name"]
    for finding in VALIDATION_FINDINGS
    if finding["is_hard_stop"] and finding["status"] == "fail"
]
warning_check_names = [
    finding["check_name"]
    for finding in VALIDATION_FINDINGS
    if finding["status"] == "warning"
]
overall_status = "fail" if len(failed_hard_stop_check_names) > 0 else "pass"

validation_report = {
    "validation_report_schema_version": VALIDATION_REPORT_SCHEMA_VERSION,
    "overall_status": overall_status,
    "failed_hard_stop_check_names": failed_hard_stop_check_names,
    "warning_check_names": warning_check_names,
    "generated_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "input_file": {
        "given_path": parsed_arguments.input_tiff_path,
        "lexically_absolute_path": str(input_tiff_absolute_path),
        "symlink_resolved_path": str(input_tiff_physical_path),
        "filename": input_tiff_absolute_path.name,
        "size_bytes": int(input_tiff_file_status.st_size),
        "modified_at": datetime.datetime.fromtimestamp(
            input_tiff_file_status.st_mtime, datetime.timezone.utc
        ).isoformat(timespec="seconds"),
    },
    "tiff_container": {
        "byte_order": str(tiff_file_handle.byteorder),
        "is_bigtiff": bool(tiff_file_handle.is_bigtiff),
        "page_count": tiff_page_count,
        "selected_page_index": selected_page_index,
        "page_selection_was_explicit": parsed_arguments.page_index is not None,
    },
    "tiff_page": {
        "image_width_pixels": page_image_width_pixels,
        "image_length_pixels": page_image_length_pixels,
        "bits_per_sample": int(page_bits_per_sample),
        "samples_per_pixel": int(page_samples_per_pixel),
        "sample_format_repr": repr(page_sample_format),
        "compression_integer": int(page_compression),
        "compression_name": page_compression_name,
        "photometric_integer": int(page_photometric),
        "photometric_name": page_photometric_name,
        "decoded_dtype": str(pixel_array.dtype) if pixel_array is not None else None,
        "decoded_width_pixels": pixel_array_width_pixels
        if pixel_statistics_were_computed
        else None,
        "decoded_height_pixels": pixel_array_height_pixels
        if pixel_statistics_were_computed
        else None,
        "estimated_float64_footprint_bytes": (
            int(page_image_width_pixels * page_image_length_pixels * 8)
        ),
    },
    "geometry": {
        "micrometres_per_pixel": micrometres_per_pixel,
        "micrometres_per_pixel_by_axis": micrometres_per_pixel_by_axis,
        "resolution_unit_integer": resolution_unit_value,
        "image_extent_millimetres": (
            [
                pixel_array_width_pixels * micrometres_per_pixel / 1000.0,
                pixel_array_height_pixels * micrometres_per_pixel / 1000.0,
            ]
            if micrometres_per_pixel is not None
            else None
        ),
        "orientation_sources": orientation_sources,
    },
    "pixel_statistics_were_computed": pixel_statistics_were_computed,
    "pixel_statistics": None
    if not pixel_statistics_were_computed
    else {
        "total_pixel_count": total_pixel_count,
        "minimum_value": minimum_pixel_value,
        "maximum_value": maximum_pixel_value,
        "mean_value": mean_pixel_value,
        "median_value": median_pixel_value,
        "standard_deviation_value": standard_deviation_pixel_value,
        "median_absolute_deviation_counts": median_absolute_deviation_counts,
        "robust_standard_deviation_counts": robust_standard_deviation_counts,
        "histogram_mode_value": histogram_mode_value,
        "pixels_above_mode": pixels_above_mode,
        "pixels_below_mode": pixels_below_mode,
        "inferred_signal_direction": inferred_signal_direction,
        "distinct_value_count": distinct_pixel_value_count,
        "value_stride": pixel_value_stride,
        "effective_bits_per_sample": effective_bits_per_sample,
        "inferred_ceiling_value": inferred_ceiling_value,
        "at_inferred_ceiling_pixel_count": at_inferred_ceiling_pixel_count,
        "at_inferred_ceiling_pixel_fraction": at_inferred_ceiling_pixel_count
        / total_pixel_count,
        "at_container_maximum_pixel_count": at_container_maximum_pixel_count,
        "at_container_maximum_pixel_fraction": at_container_maximum_pixel_count
        / total_pixel_count,
        "at_maximum_observed_pixel_count": at_maximum_observed_pixel_count,
        "saturation_spike_ratio": saturation_spike_ratio,
        "at_floor_pixel_count": at_floor_pixel_count,
        "at_floor_pixel_fraction": floor_population_fraction,
        "at_minimum_observed_pixel_count": at_minimum_observed_pixel_count,
        "isolated_extreme_pixel_count": isolated_extreme_pixel_count,
        "isolation_margin_counts": isolation_margin_counts,
        "isolation_margin_robust_sigmas": ISOLATION_MARGIN_ROBUST_SIGMAS,
        "estimated_lane_tilt_angle_degrees": None,
    },
    "linearity_evidence": {
        "encoding_verified": False,
        "read_path_container": "tif",
        "vendor_claim_for_tif": "linear 16-bit grayscale, per Cytiva/GE documentation as "
        "reported in DESIGN.md section 10",
        "inf_scale_type": None,
        "inf_signal_process_2": None,
        "exposure_time_text": image_description_first_value_by_key.get("Exposure time"),
        "test_that_would_settle_it": (
            "A dilution series or a known-dose standard imaged on the same instrument "
            "and settings. The .img-versus-.tif scatter was deferred by decision "
            "because it compares two containers without a calibrated reference: it can "
            "show that they differ but not which is linear in signal, so DESIGN.md "
            "section 10's claim that it decides the question is an overclaim."
        ),
        "consequence_if_wrong": (
            "Every ratio this pipeline produces is compressed toward 1 and nothing in "
            "the pipeline notices. For Amersham Imager 680 fluorescence files the risk "
            "is much lower: a CCD cooled to -25 C is linear in photons by construction, "
            "but two of its images are comparable only after dividing by exposure time."
        ),
    },
    "image_description": {
        "raw_texts": image_description_raw_texts,
        "ordered_pairs": image_description_pairs,
        "duplicate_keys": image_description_duplicate_keys,
    },
    "tiff_tags": {
        "standard": standard_tag_records,
        "private_at_or_above_" + str(PRIVATE_TAG_CODE_FLOOR): private_tag_records,
        "duplicated_tag_names": duplicated_tag_names,
    },
    "inf_sidecar": inf_sidecar_record,
    "preprocess_sidecar": preprocess_sidecar_record,
    "outputs": {
        "validation_report": str(validation_report_output_path),
        "histogram_image": str(histogram_image_output_path)
        if histogram_was_written
        else None,
        "crop_preview_image": (
            str(crop_preview_image_output_path) if crop_preview_was_written else None
        ),
    },
    "validation_findings": VALIDATION_FINDINGS,
    "run_log_lines": ACCUMULATED_RUN_LOG_LINES,
}

for keyed_prefix, keyed_name, keyed_value in inf_sidecar_record["keyed_pairs"] or []:
    if keyed_name == "ScaleType":
        validation_report["linearity_evidence"]["inf_scale_type"] = keyed_value
    elif keyed_name == "SignalProcess2":
        validation_report["linearity_evidence"]["inf_signal_process_2"] = keyed_value

# ensure_ascii is the default and is relied on: CONVENTIONS.md section 4 requires
# ASCII output, and the .inf and ImageDescription were decoded permissively as
# latin-1 so that a byte above 127 is recorded rather than crashing the run.
try:
    validation_report_text = json.dumps(validation_report, indent=2, sort_keys=False)
except TypeError as serialization_error:
    # numpy scalars, numpy dtypes and bytes all raise here rather than serializing
    # wrongly, which is the good case. If this fires, a value reached the report
    # without being converted, and writing a partial report would be worse than
    # saying so.
    die(
        "report",
        "the report contains a value that will not serialize: "
        + str(serialization_error),
    )
try:
    validation_report_output_path.write_text(
        validation_report_text + "\n", encoding="ascii"
    )
except OSError as report_write_error:
    die(
        "report",
        "could not write "
        + str(validation_report_output_path)
        + ": "
        + str(report_write_error),
    )
emit_message("report", "wrote " + str(validation_report_output_path))

tiff_file_handle.close()

# stdout carries only this, tab separated, fixed field order, so it stays pipeable.
print("overall_status\t" + overall_status)
print("failed_hard_stop_count\t" + str(len(failed_hard_stop_check_names)))
print("warning_count\t" + str(len(warning_check_names)))
print("validation_report_path\t" + str(validation_report_output_path))
print("histogram_image_path\t" + str(histogram_image_output_path))
print(
    "crop_preview_image_path\t"
    + (
        str(crop_preview_image_output_path)
        if crop_preview_was_written
        else "not written"
    )
)

for finding in VALIDATION_FINDINGS:
    if finding["status"] == "fail":
        emit_message(
            "finding", "FAIL " + finding["check_name"] + ": " + finding["detail"]
        )
    elif finding["status"] == "warning":
        emit_message(
            "finding", "WARNING " + finding["check_name"] + ": " + finding["detail"]
        )

emit_message(
    "stage 1",
    "complete; overall status "
    + overall_status
    + ", "
    + str(len(failed_hard_stop_check_names))
    + " failed hard stop(s), "
    + str(len(warning_check_names))
    + " warning(s). The report was written either way, "
    "so a failure is diagnosable from the file rather than only from this scrollback.",
)
sys.exit(1 if overall_status == "fail" else 0)
