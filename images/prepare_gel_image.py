# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
r"""
Stage 0 of the gel densitometry pipeline: argument handling, path normalization,
existence and readability checks, output directory creation. Prints resolved
absolute paths and exits. Reads no image data.

Scope is DESIGN.md section 7 stage 0. The shape of this file is dictated by
CONVENTIONS.md: flat procedural, no helpers other than the tagged stderr emitter
and the fail-fast exit, lengths and units in names.

No third-party dependencies on purpose. Stage 0 is what you run to find out
whether the environment and the paths are sane, so it must not be able to fail
for an environment reason of its own.

Messages go to stderr, tagged with source. Only resolved absolute paths go to
stdout, tab separated, so stdout stays pipeable.

Usage. Single quotes throughout, because real filenames here contain spaces,
commas and square brackets, and square brackets are glob metacharacters in bash:

    uv run prepare_gel_image.py '/home/luised94/gel_data/20220303-wt,4r,5,6-sofa-repeat-1-1000-[Phosphor].img'
    uv run prepare_gel_image.py 'C:\Users\liusm\MIT Dropbox\scan [Phosphor].img'
    uv run prepare_gel_image.py --output-parent-directory '~/gel_scratch' 'scan.img'
    uv run prepare_gel_image.py -- '-filename-starting-with-a-dash.img'

An unquoted Windows path cannot be recovered: backslash is the shell escape
character, so the separators are deleted before this script sees anything. The
script detects the residue and says so rather than guessing.
"""

import argparse
import datetime
import sys

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIRECTORY_NAME_SUFFIX = "_gel_analysis"

# =============================================================================
# The two permitted helpers (CONVENTIONS.md section 1)
# =============================================================================

# Every emitted message accumulates here as well as going to stderr.
# CONVENTIONS.md section 10 requires the run log to be serialized into the
# provenance JSON, which is stage 2. Nothing consumes this list yet. It exists
# now because adding it later means revisiting every call site instead of one
# line.
ACCUMULATED_RUN_LOG_LINES = []


def emit_message(source_tag, message_text):
    # stderr, not stdout: stdout carries only the resolved paths so the script
    # can be piped. The repo's existing Python scripts print informational
    # output to stdout instead; CONVENTIONS.md section 10 overrides that.
    timestamp_text = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    tagged_line = "[" + source_tag + "] " + message_text
    ACCUMULATED_RUN_LOG_LINES.append(timestamp_text + " " + tagged_line)
    print(tagged_line, file=sys.stderr)


def die(source_tag, message_text):
    # Single exit path for every failure so that the reason is always in the run
    # log. sys.exit("[ERROR] ...") as used elsewhere in this repo does reach
    # stderr and does exit 1, but bypasses the accumulator above, which is the
    # entire reason this function exists rather than reusing that idiom.
    emit_message(source_tag, "ERROR: " + message_text)
    sys.exit(1)


# =============================================================================
# Arguments
# =============================================================================

argument_parser = argparse.ArgumentParser(
    prog="prepare_gel_image.py",
    description=(
        "Stage 0: resolve the input path, verify it is a readable regular file, "
        "create the per-input output directory, print resolved absolute paths."
    ),
    # An abbreviated flag that happens to be unambiguous today becomes ambiguous
    # the moment a second flag shares its prefix, silently changing the meaning
    # of a command that used to work.
    allow_abbrev=False,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=(
        "Quote paths with single quotes. Real filenames contain spaces, commas\n"
        "and square brackets, and square brackets are glob metacharacters in\n"
        "bash. Use -- before a filename that starts with a dash.\n"
    ),
)
argument_parser.add_argument(
    "input_image_path",
    help=(
        "Path to the scan to prepare. WSL path, Windows path with a drive "
        "letter, ~ prefixed, relative or absolute."
    ),
)
argument_parser.add_argument(
    "--output-parent-directory",
    default=None,
    help=(
        "Directory in which to create <stem>" + OUTPUT_DIRECTORY_NAME_SUFFIX + ". "
        "Defaults to the directory holding the input, per DESIGN.md 5.9."
    ),
)
parsed_arguments = argument_parser.parse_args()

emit_message("arguments", "input image path as given: " + repr(parsed_arguments.input_image_path))
if parsed_arguments.output_parent_directory is not None:
    emit_message(
        "arguments",
        "output parent directory as given: " + repr(parsed_arguments.output_parent_directory),
    )

emit_message("stage 0", "argument handling complete")
sys.exit(0)
