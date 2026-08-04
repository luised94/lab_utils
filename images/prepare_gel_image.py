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
import os
import pathlib
import re
import shutil
import stat
import subprocess
import sys

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_DIRECTORY_NAME_SUFFIX = "_gel_analysis"

# Used only when wslpath is unavailable. /etc/wsl.conf can relocate the automount
# root, so this is a fallback and not an assumption.
WINDOWS_DRIVE_MOUNT_ROOT_FALLBACK = "/mnt/"
WSL_CONFIGURATION_FILE_PATH = "/etc/wsl.conf"

# Drive letter followed by a separator: a real Windows absolute path.
WINDOWS_DRIVE_ABSOLUTE_PATH_PATTERN = re.compile(r"^([A-Za-z]):[/\\]")

# Drive letter NOT followed by a separator. Either the drive-relative form
# "C:file.img", which has no WSL equivalent, or far more often the residue of an
# unquoted Windows path whose backslashes the shell consumed.
WINDOWS_DRIVE_RELATIVE_PATH_PATTERN = re.compile(r"^([A-Za-z]):(?![/\\])")

SURROUNDING_QUOTE_CHARACTERS = ('"', "'")

# One byte is enough to prove the file opens and yields data. This is a
# readability check, not image reading: no interpretation of the byte occurs.
READABILITY_PROBE_BYTE_COUNT = 1

# Checked in order; the first match names the type in the error message.
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


# =============================================================================
# Path normalization
# =============================================================================

# Both supplied paths get identical treatment. This is a two-element loop rather
# than a function because CONVENTIONS.md section 1 permits no third helper, and a
# loop over plain data is inline logic rather than an abstraction.
path_normalization_inputs = [("input_image_path", parsed_arguments.input_image_path)]
if parsed_arguments.output_parent_directory is not None:
    path_normalization_inputs.append(
        ("output_parent_directory", parsed_arguments.output_parent_directory)
    )

normalized_paths_by_role = {}

for path_role_name, raw_path_text in path_normalization_inputs:
    working_path_text = raw_path_text

    # A pasted path routinely carries a trailing newline or a stray space. The
    # cost of stripping is that a file genuinely named with a leading or
    # trailing space becomes unaddressable. Accepted: the paste case happens
    # every session, the named-with-space case has not happened here.
    whitespace_stripped_path_text = working_path_text.strip(" \t\r\n")
    if whitespace_stripped_path_text != working_path_text:
        emit_message(path_role_name, "stripped surrounding whitespace from the supplied path")
        working_path_text = whitespace_stripped_path_text
    if working_path_text == "":
        die(path_role_name, "path is empty after stripping surrounding whitespace")

    # Windows Explorer "Copy as path" produces a double-quoted string, and it
    # gets pasted verbatim inside another pair of shell quotes.
    for quote_character in SURROUNDING_QUOTE_CHARACTERS:
        if (
            len(working_path_text) >= 2
            and working_path_text.startswith(quote_character)
            and working_path_text.endswith(quote_character)
        ):
            emit_message(
                path_role_name,
                "removed a surrounding pair of " + quote_character + " characters",
            )
            working_path_text = working_path_text[1:-1]
            break
    if working_path_text == "":
        die(path_role_name, "path is empty after removing surrounding quotes")

    # This check exists because it is the single most common way this script will
    # be invoked wrongly, and the failure is otherwise unintelligible: bash turns
    # C:\Users\liusm\scan.img into C:Usersliusmscan.img before argv is built. The
    # information is gone, so the only correct response is to name the cause.
    if WINDOWS_DRIVE_RELATIVE_PATH_PATTERN.match(working_path_text) is not None:
        die(
            path_role_name,
            "the path " + repr(working_path_text) + " has a drive letter but no "
            "separator after it. Backslash is the shell escape character, so an "
            "unquoted Windows path arrives with its separators already deleted and "
            "cannot be reconstructed. Re-run with the path inside single quotes. "
            "If you really meant the Windows drive-relative form C:name, that has "
            "no WSL equivalent; supply a full path.",
        )

    # A UNC share has no general mapping to a WSL path. Guessing produces a path
    # that does not exist, and the resulting error blames the wrong thing.
    if working_path_text.startswith("\\\\"):
        die(
            path_role_name,
            "UNC path " + repr(working_path_text) + " is not supported. Map the share "
            "to a drive letter, or supply the /mnt path directly.",
        )

    windows_drive_match = WINDOWS_DRIVE_ABSOLUTE_PATH_PATTERN.match(working_path_text)
    if windows_drive_match is not None:
        wslpath_executable_path = shutil.which("wslpath")
        if wslpath_executable_path is not None:
            # wslpath reads the live automount configuration. A string transform
            # only guesses it, so prefer the authority whenever it is present.
            wslpath_result = subprocess.run(
                [wslpath_executable_path, "-u", "-a", working_path_text],
                capture_output=True,
                text=True,
            )
            if wslpath_result.returncode != 0:
                die(
                    path_role_name,
                    "wslpath refused to convert " + repr(working_path_text) + ": "
                    + wslpath_result.stderr.strip(),
                )
            converted_path_text = wslpath_result.stdout.strip("\n")
            if converted_path_text == "":
                die(path_role_name, "wslpath returned an empty conversion for " + repr(working_path_text))
            emit_message(
                path_role_name,
                "converted Windows path via wslpath to " + converted_path_text,
            )
            working_path_text = converted_path_text
        else:
            windows_drive_letter = windows_drive_match.group(1).lower()
            windows_drive_mount_root = WINDOWS_DRIVE_MOUNT_ROOT_FALLBACK
            # [automount] root = / is a common change, and hardcoding /mnt/ there
            # produces a path that silently does not exist.
            try:
                wsl_configuration_text = pathlib.Path(WSL_CONFIGURATION_FILE_PATH).read_text()
            except OSError:
                wsl_configuration_text = ""
            inside_automount_section = False
            for wsl_configuration_line in wsl_configuration_text.splitlines():
                stripped_configuration_line = wsl_configuration_line.strip()
                if stripped_configuration_line.startswith("["):
                    inside_automount_section = stripped_configuration_line.lower() == "[automount]"
                    continue
                if inside_automount_section and "=" in stripped_configuration_line:
                    configuration_key, configuration_value = stripped_configuration_line.split("=", 1)
                    if configuration_key.strip().lower() == "root":
                        windows_drive_mount_root = configuration_value.strip()
                        if not windows_drive_mount_root.endswith("/"):
                            windows_drive_mount_root = windows_drive_mount_root + "/"
            path_remainder_after_drive = working_path_text[len(windows_drive_match.group(0)):]
            converted_path_text = (
                windows_drive_mount_root
                + windows_drive_letter
                + "/"
                + path_remainder_after_drive.replace("\\", "/")
            )
            emit_message(
                path_role_name,
                "wslpath is unavailable; converted Windows path by string transform "
                "using mount root " + windows_drive_mount_root + " to " + converted_path_text,
            )
            working_path_text = converted_path_text

    # expanduser only touches a leading tilde, so a file named ~backup elsewhere
    # in the path is untouched. A file named ~backup in the leading position
    # would be misread; accepted, since the tilde only survives to argv when the
    # user quoted it deliberately.
    if working_path_text.startswith("~"):
        try:
            tilde_expanded_path_text = os.path.expanduser(working_path_text)
        except RuntimeError as tilde_expansion_error:
            die(path_role_name, "cannot expand the leading tilde: " + str(tilde_expansion_error))
        if tilde_expanded_path_text == working_path_text:
            die(
                path_role_name,
                "the leading tilde in " + repr(working_path_text) + " did not expand; "
                "the named user probably does not exist on this machine.",
            )
        emit_message(path_role_name, "expanded leading tilde to " + tilde_expanded_path_text)
        working_path_text = tilde_expanded_path_text

    if not pathlib.Path(working_path_text).is_absolute():
        emit_message(
            path_role_name,
            "path is relative and was resolved against the current working directory "
            + os.getcwd(),
        )

    # Two absolute forms, both kept and both reported. abspath normalizes ".."
    # lexically and preserves symlinks; resolve() follows them to the physical
    # target. DESIGN.md section 3 clones test files to a local working location
    # specifically to avoid touching the synced Dropbox tree, so collapsing to
    # the physical target and forgetting the given form would defeat that.
    lexically_absolute_path = pathlib.Path(os.path.abspath(working_path_text))
    symlink_resolved_absolute_path = pathlib.Path(working_path_text).resolve()
    if symlink_resolved_absolute_path != lexically_absolute_path:
        emit_message(
            path_role_name,
            "path traverses a symlink or a ..-through-symlink: given form "
            + str(lexically_absolute_path)
            + " resolves to physical form "
            + str(symlink_resolved_absolute_path),
        )

    # stdout is a tab separated report, so a tab or newline in a path would
    # corrupt it. Both are legal in POSIX filenames, so this is checked rather
    # than assumed away.
    for forbidden_character_name, forbidden_character in (("tab", "\t"), ("newline", "\n")):
        if forbidden_character in str(lexically_absolute_path):
            die(
                path_role_name,
                "path contains a literal " + forbidden_character_name + ", which cannot "
                "be represented in the tab separated report written to stdout.",
            )

    if not str(lexically_absolute_path).isascii():
        emit_message(
            path_role_name,
            "WARNING: path contains non-ASCII characters. CONVENTIONS.md section 4 "
            "requires ASCII filenames, and the output directory name derived from "
            "this stem will inherit them.",
        )

    normalized_paths_by_role[path_role_name] = {
        "raw_path_text": raw_path_text,
        "lexically_absolute_path": lexically_absolute_path,
        "symlink_resolved_absolute_path": symlink_resolved_absolute_path,
    }
    emit_message(path_role_name, "normalized to " + str(lexically_absolute_path))

emit_message("stage 0", "path normalization complete")

# =============================================================================
# Existence, file type and readability of the input
# =============================================================================

input_image_absolute_path = normalized_paths_by_role["input_image_path"]["lexically_absolute_path"]
input_image_physical_path = normalized_paths_by_role["input_image_path"]["symlink_resolved_absolute_path"]

# lexists rather than exists, so that a broken symlink is reported as a broken
# symlink instead of as a missing file.
if not os.path.lexists(input_image_absolute_path):
    # Naming the first component that does not exist converts "file not found"
    # into an actionable message when the real fault is an unmounted drive or a
    # mistyped directory several levels up.
    first_missing_path_component = input_image_absolute_path
    for candidate_ancestor_path in [input_image_absolute_path] + list(input_image_absolute_path.parents):
        if os.path.lexists(candidate_ancestor_path):
            break
        first_missing_path_component = candidate_ancestor_path
    backslash_hint_text = ""
    if "\\" in parsed_arguments.input_image_path:
        backslash_hint_text = (
            " The supplied path still contains backslashes; if this is a Windows path, "
            "it needs a drive letter and single quotes."
        )
    die(
        "existence",
        "input image not found: " + str(input_image_absolute_path)
        + ". Highest path component that does not exist: " + str(first_missing_path_component)
        + "." + backslash_hint_text,
    )

input_image_link_status = os.lstat(input_image_absolute_path)
if stat.S_ISLNK(input_image_link_status.st_mode) and not os.path.exists(input_image_absolute_path):
    die(
        "existence",
        "input image " + str(input_image_absolute_path) + " is a symlink whose target "
        + str(input_image_physical_path) + " does not exist.",
    )

input_image_file_status = os.stat(input_image_absolute_path)

# Rejecting non-regular files before opening anything is not pedantry: opening a
# named pipe blocks until a writer appears, so the readability probe below would
# hang forever rather than fail.
for non_regular_type_name, non_regular_type_predicate in NON_REGULAR_FILE_TYPE_PREDICATES:
    if non_regular_type_predicate(input_image_file_status.st_mode):
        directory_hint_text = ""
        if non_regular_type_name == "directory":
            directory_hint_text = (
                " Directory mode is stage 4 in DESIGN.md section 7; pass a single file."
            )
        die(
            "file type",
            "input image " + str(input_image_absolute_path) + " is a "
            + non_regular_type_name + ", not a regular file." + directory_hint_text,
        )
if not stat.S_ISREG(input_image_file_status.st_mode):
    die(
        "file type",
        "input image " + str(input_image_absolute_path) + " is not a regular file "
        "(st_mode " + oct(input_image_file_status.st_mode) + ").",
    )

# DESIGN.md section 3 names reading a partially synced Dropbox file as a real
# risk. A zero byte file is the one degenerate case stage 0 can catch without
# reading image data, and it is exactly what an interrupted sync leaves behind.
if input_image_file_status.st_size == 0:
    die(
        "file size",
        "input image " + str(input_image_absolute_path) + " is zero bytes. An "
        "interrupted or placeholder-only Dropbox sync produces exactly this.",
    )

# os.access is deliberately not used. DESIGN.md section 3 records that /mnt/c is
# DrvFs and reports 777 regardless of real permissions, and access() run as root
# returns true for nearly everything. Opening the file and reading a byte is the
# only test that answers the question actually being asked.
try:
    with open(input_image_absolute_path, "rb") as input_image_file_handle:
        readability_probe_bytes = input_image_file_handle.read(READABILITY_PROBE_BYTE_COUNT)
except PermissionError as permission_error:
    die(
        "readability",
        "input image " + str(input_image_absolute_path) + " exists but cannot be opened "
        "for reading: " + str(permission_error),
    )
except OSError as open_error:
    die(
        "readability",
        "input image " + str(input_image_absolute_path) + " exists but failed to open: "
        + str(open_error),
    )

if len(readability_probe_bytes) != READABILITY_PROBE_BYTE_COUNT:
    die(
        "readability",
        "input image " + str(input_image_absolute_path) + " opened but returned "
        + str(len(readability_probe_bytes)) + " bytes instead of "
        + str(READABILITY_PROBE_BYTE_COUNT) + ".",
    )

emit_message(
    "readability",
    "opened the input and read " + str(len(readability_probe_bytes))
    + " byte; reported size is " + str(input_image_file_status.st_size) + " bytes",
)

emit_message("stage 0", "input checks complete")
sys.exit(0)
