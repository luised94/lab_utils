# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "tifffile",
# ]
# ///
r"""
Dump every TIFF tag in a file verbatim, including private tags 65000 and above.
Reports. Interprets nothing. Validates nothing.

This exists because DESIGN.md section 6 requires a full tag inventory "so nothing
is silently ignored", and because the validation rules in section 6 cannot be
written until it is known which of those tags these files actually carry. Run this
first; the output decides what stage 1 checks.

Self-contained on purpose. Boilerplate is duplicated from prepare_gel_image.py
rather than shared, so this file stays a reference tool that can be run against any
TIFF years from now without the rest of the pipeline existing.

What was NOT duplicated, and why: the Windows drive-letter conversion. This tool is
run against paths that already exist on this machine, in the /mnt/c form. The
de-escaped-paste hard stop is kept, because that mistake produces an
unintelligible error otherwise.

Pixel data is read only to report shape and dtype, which tifffile takes from the
tags without decoding the image.

Usage. Single quotes, because real filenames here contain spaces, commas and
square brackets, and square brackets are glob metacharacters in bash:

    uv run inventory_tiff_tags.py '/mnt/c/Users/liusm/.../20220303-...-[Phosphor].tif'
    uv run inventory_tiff_tags.py -- '-filename-starting-with-a-dash.tif'
"""

import argparse
import datetime
import os
import pathlib
import re
import stat
import sys

import tifffile

# =============================================================================
# Configuration
# =============================================================================

# Anything at or above this code is vendor private. DESIGN.md section 6 names
# MD_FileTag and MD_ScalePixel in the 65000-65003 range and requires them dumped
# verbatim, so they are reported separately rather than buried in the main list.
PRIVATE_TAG_CODE_FLOOR = 65000

# Long tag values are truncated in the report. A full ColorMap or a JPEG table
# would otherwise bury everything else.
TAG_VALUE_CHARACTER_LIMIT = 300

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

# =============================================================================
# The two permitted helpers (CONVENTIONS.md section 1)
# =============================================================================

ACCUMULATED_RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    # stderr, so stdout carries only the inventory and can be redirected to a file
    # or diffed between two scans.
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
    prog="inventory_tiff_tags.py",
    description="Dump every TIFF tag verbatim. Reports only; interprets nothing.",
    allow_abbrev=False,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="Quote paths with single quotes. Use -- before a filename starting with a dash.\n",
)
argument_parser.add_argument("input_tiff_path", help="Path to the TIFF to inventory.")
parsed_arguments = argument_parser.parse_args()

# =============================================================================
# Path normalization
# =============================================================================

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
        emit_message("path", "removed a surrounding pair of " + quote_character + " characters")
        working_path_text = working_path_text[1:-1]
        break

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
        die("path", "the leading tilde in " + repr(working_path_text) + " did not expand")
    emit_message("path", "expanded leading tilde to " + tilde_expanded_path_text)
    working_path_text = tilde_expanded_path_text

input_tiff_absolute_path = pathlib.Path(os.path.abspath(working_path_text))
input_tiff_physical_path = pathlib.Path(working_path_text).resolve()
if input_tiff_physical_path != input_tiff_absolute_path:
    emit_message(
        "path",
        "path traverses a symlink: given form " + str(input_tiff_absolute_path)
        + " resolves to physical form " + str(input_tiff_physical_path),
    )
emit_message("path", "normalized to " + str(input_tiff_absolute_path))

# =============================================================================
# Existence, file type and readability
# =============================================================================

if not os.path.lexists(input_tiff_absolute_path):
    first_missing_path_component = input_tiff_absolute_path
    for candidate_ancestor_path in [input_tiff_absolute_path] + list(input_tiff_absolute_path.parents):
        if os.path.lexists(candidate_ancestor_path):
            break
        first_missing_path_component = candidate_ancestor_path
    die(
        "existence",
        "input not found: " + str(input_tiff_absolute_path)
        + ". Highest path component that does not exist: " + str(first_missing_path_component),
    )

input_tiff_link_status = os.lstat(input_tiff_absolute_path)
if stat.S_ISLNK(input_tiff_link_status.st_mode) and not os.path.exists(input_tiff_absolute_path):
    die(
        "existence",
        str(input_tiff_absolute_path) + " is a symlink whose target "
        + str(input_tiff_physical_path) + " does not exist.",
    )

input_tiff_file_status = os.stat(input_tiff_absolute_path)

# Rejecting non-regular files before opening anything: opening a named pipe blocks
# until a writer appears, so the probe below would hang rather than fail.
for non_regular_type_name, non_regular_type_predicate in NON_REGULAR_FILE_TYPE_PREDICATES:
    if non_regular_type_predicate(input_tiff_file_status.st_mode):
        die(
            "file type",
            str(input_tiff_absolute_path) + " is a " + non_regular_type_name
            + ", not a regular file.",
        )
if not stat.S_ISREG(input_tiff_file_status.st_mode):
    die("file type", str(input_tiff_absolute_path) + " is not a regular file.")

if input_tiff_file_status.st_size == 0:
    die("file size", str(input_tiff_absolute_path) + " is zero bytes.")

# os.access is deliberately not used: /mnt/c is DrvFs and reports 777 regardless of
# real permissions, and access() as root returns true for nearly everything.
try:
    with open(input_tiff_absolute_path, "rb") as input_tiff_file_handle:
        readability_probe_bytes = input_tiff_file_handle.read(READABILITY_PROBE_BYTE_COUNT)
except OSError as open_error:
    die("readability", str(input_tiff_absolute_path) + " failed to open: " + str(open_error))
if len(readability_probe_bytes) != READABILITY_PROBE_BYTE_COUNT:
    die("readability", str(input_tiff_absolute_path) + " opened but returned no bytes.")

emit_message(
    "readability",
    "opened and read 1 byte; reported size is " + str(input_tiff_file_status.st_size) + " bytes",
)

# =============================================================================
# Inventory
# =============================================================================

try:
    tiff_file_handle = tifffile.TiffFile(input_tiff_absolute_path)
except Exception as tiff_open_error:
    # Deliberately broad. The point of this tool is to find out what the file is,
    # so any failure to parse it as TIFF is a result to report, not a crash.
    die(
        "tiff",
        "tifffile could not parse " + str(input_tiff_absolute_path) + " as TIFF: "
        + type(tiff_open_error).__name__ + ": " + str(tiff_open_error),
    )

print("file_path\t" + str(input_tiff_absolute_path))
print("file_size_bytes\t" + str(input_tiff_file_status.st_size))
print("byte_order\t" + str(tiff_file_handle.byteorder))
print("is_bigtiff\t" + str(tiff_file_handle.is_bigtiff))
print("page_count\t" + str(len(tiff_file_handle.pages)))

# DESIGN.md section 6 makes more than one page require explicit selection and
# forbids defaulting to page zero, so the count is reported before anything else
# and every page is inventoried rather than just the first.
for page_index, tiff_page in enumerate(tiff_file_handle.pages):
    print()
    print("=== page " + str(page_index) + " ===")
    for page_attribute_name in ("shape", "dtype", "axes", "photometric", "compression",
                                "planarconfig", "predictor", "sampleformat",
                                "bitspersample", "samplesperpixel"):
        page_attribute_value = getattr(tiff_page, page_attribute_name, None)
        print("page_" + page_attribute_name + "\t" + repr(page_attribute_value))

    ordinary_tag_lines = []
    private_tag_lines = []
    for tiff_tag in tiff_page.tags:
        tag_value_text = repr(tiff_tag.value)
        if len(tag_value_text) > TAG_VALUE_CHARACTER_LIMIT:
            tag_value_text = (
                tag_value_text[:TAG_VALUE_CHARACTER_LIMIT]
                + "... [truncated, " + str(len(tag_value_text)) + " characters total]"
            )
        # repr keeps the output ASCII even when a tag holds arbitrary bytes, which
        # CONVENTIONS.md section 4 requires of every emitted message.
        tag_line = "\t".join([
            str(tiff_tag.code),
            str(tiff_tag.name),
            str(tiff_tag.dtype),
            str(tiff_tag.count),
            tag_value_text,
        ])
        if tiff_tag.code >= PRIVATE_TAG_CODE_FLOOR:
            private_tag_lines.append(tag_line)
        else:
            ordinary_tag_lines.append(tag_line)

    print()
    print("-- standard tags, page " + str(page_index) + " (code, name, dtype, count, value) --")
    for tag_line in sorted(ordinary_tag_lines, key=lambda line: int(line.split("\t")[0])):
        print(tag_line)

    print()
    print("-- private tags " + str(PRIVATE_TAG_CODE_FLOOR) + " and above, page "
          + str(page_index) + " --")
    if len(private_tag_lines) == 0:
        print("(none)")
    else:
        for tag_line in sorted(private_tag_lines, key=lambda line: int(line.split("\t")[0])):
            print(tag_line)

tiff_file_handle.close()

emit_message(
    "inventory",
    "complete; " + str(len(tiff_file_handle.pages)) + " page(s) inventoried, nothing interpreted",
)
sys.exit(0)
