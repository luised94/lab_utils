#!/usr/bin/env python3
"""Tests for lw -- bare python, temp dirs, no test framework (SS9).

One test per command: a happy path, an error path, and a no-effect
assertion (a rejected command leaves state untouched). The no-effect
assertion is the test that directly exercises the effects boundary -- a
rejected command must never reach write_store_to_database, the counter
bump, or the filesystem.

For `link`, the no-effect assertion is mandatory across all FOUR rejection
paths (bad id, missing experiment, unknown relationship, duplicate edge),
each a distinct no-effect case.

Run: python3 test_lw.py
The harness sets a fresh temp LW_LAB_ROOT and re-imports lw so its
module-level path constants resolve into the sandbox.
"""

import datetime
import importlib.util
import os
import sys
import tempfile


# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

LW_MODULE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lw.py")

# A frozen clock so folder names and dates are deterministic across runs.
FROZEN_DATE = datetime.date(2026, 6, 24)
FROZEN_DATETIME = datetime.datetime(2026, 6, 24, 9, 30, 0)


def load_lw_into_sandbox(sandbox_lab_root):
    """Import a FRESH copy of lw with LW_LAB_ROOT pointed at the sandbox.

    lw resolves its paths at import time, so each test gets its own module
    instance bound to its own temp directory. The frozen clock is injected
    after import so date-derived output is deterministic.
    """
    os.environ["LW_LAB_ROOT"] = sandbox_lab_root
    module_specification = importlib.util.spec_from_file_location(
        "lw_under_test", LW_MODULE_PATH
    )
    lw = importlib.util.module_from_spec(module_specification)
    module_specification.loader.exec_module(lw)
    lw.startup_clock = {"today": FROZEN_DATE, "now": FROZEN_DATETIME}
    return lw


class CapturedRun:
    """Result of running a command: exit code plus captured stdout/stderr."""

    def __init__(self, exit_code, stdout_text, stderr_text):
        self.exit_code = exit_code
        self.stdout_text = stdout_text
        self.stderr_text = stderr_text


def run_command(lw, command_arguments):
    """Invoke lw.main with stdout/stderr captured, return a CapturedRun."""
    import io

    original_stdout = sys.stdout
    original_stderr = sys.stderr
    captured_stdout = io.StringIO()
    captured_stderr = io.StringIO()
    sys.stdout = captured_stdout
    sys.stderr = captured_stderr
    try:
        exit_code = lw.main(command_arguments)
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
    return CapturedRun(
        exit_code,
        captured_stdout.getvalue(),
        captured_stderr.getvalue(),
    )


def snapshot_store_state(lw):
    """Capture the full persistent state that a no-effect test must find
    unchanged: the two live tables plus the counter value plus the set of
    experiment folder names on disk.
    """
    store = lw.load_store_from_database(lw.DATABASE_PATH)
    counter_value = lw.read_next_experiment_number()
    if os.path.isdir(lw.EXPERIMENTS_DIRECTORY):
        folder_names = sorted(os.listdir(lw.EXPERIMENTS_DIRECTORY))
    else:
        folder_names = []
    return {
        "experiments": store["experiments"],
        "experiment_links": store["experiment_links"],
        "counter": counter_value,
        "folders": folder_names,
    }


# Minimal assertion helpers, each raising AssertionError with context.


def assert_equal(actual, expected, label):
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def assert_true(condition, label):
    if not condition:
        raise AssertionError(f"{label}: expected truthy, got {condition!r}")


def assert_contains(haystack, needle, label):
    if needle not in haystack:
        raise AssertionError(f"{label}: expected {needle!r} in {haystack!r}")


def assert_no_effect(before_state, after_state, label):
    """The effects-boundary assertion: nothing persistent changed."""
    for state_key in ("experiments", "experiment_links", "counter", "folders"):
        assert_equal(
            after_state[state_key],
            before_state[state_key],
            f"{label} / {state_key} unchanged",
        )


# ---------------------------------------------------------------------------
# Per-test fixtures
# ---------------------------------------------------------------------------


def fresh_initialized_lw(sandbox_lab_root):
    """A sandbox with `init` already run: dirs, empty DB, counter seeded."""
    lw = load_lw_into_sandbox(sandbox_lab_root)
    run_command(lw, ["init"])
    return lw


def seed_two_experiments(lw):
    """Create LM-0001 (purification) and LM-0002 (cloning) for link/show
    tests. Returns nothing; both rows and folders now exist."""
    run_command(lw, ["new", "purification", "prep-a", "--title", "Prep A"])
    run_command(lw, ["new", "cloning", "clone-b", "--title", "Clone B"])


# ---------------------------------------------------------------------------
# init
# ---------------------------------------------------------------------------


def test_init_happy_path(sandbox_lab_root):
    lw = load_lw_into_sandbox(sandbox_lab_root)
    result = run_command(lw, ["init"])
    assert_equal(result.exit_code, 0, "init exit code")
    assert_true(os.path.isdir(lw.LAB_ROOT), "init created LAB_ROOT")
    assert_true(
        os.path.isdir(lw.EXPERIMENTS_DIRECTORY),
        "init created experiments dir",
    )
    assert_true(os.path.exists(lw.DATABASE_PATH), "init created database")
    assert_equal(lw.read_next_experiment_number(), 1, "init seeded counter")


def test_init_idempotent_does_not_reset_counter(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    lw.write_next_experiment_number(42)
    result = run_command(lw, ["init"])
    assert_equal(result.exit_code, 0, "second init exit code")
    assert_equal(
        lw.read_next_experiment_number(),
        42,
        "re-init must NOT reset the counter",
    )


def test_command_before_init_fails_precondition(sandbox_lab_root):
    # error path for the store boundary: no DB yet.
    lw = load_lw_into_sandbox(sandbox_lab_root)
    result = run_command(lw, ["list"])
    assert_equal(result.exit_code, 1, "pre-init command exit code")
    assert_contains(
        result.stderr_text,
        "No lab database found",
        "pre-init precondition message",
    )


# ---------------------------------------------------------------------------
# new
# ---------------------------------------------------------------------------


def test_new_happy_path(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(
        lw, ["new", "purification", "orc4r-prep", "--title", "ORC4R prep"]
    )
    assert_equal(result.exit_code, 0, "new exit code")

    store = lw.load_store_from_database(lw.DATABASE_PATH)
    assert_equal(len(store["experiments"]), 1, "new inserted one row")
    experiment_row = store["experiments"][0]
    assert_equal(experiment_row["id"], "LM-0001", "new id")
    assert_equal(
        experiment_row["folder_name"],
        "20260624_LM-0001_pur_orc4r-prep",
        "new folder name (frozen clock)",
    )
    assert_equal(experiment_row["title"], "ORC4R prep", "new title")
    assert_equal(experiment_row["date_started"], "2026-06-24", "new date_started")
    # folder created on disk
    assert_true(
        os.path.isdir(
            os.path.join(lw.EXPERIMENTS_DIRECTORY, experiment_row["folder_name"])
        ),
        "new created the experiment folder",
    )
    # counter bumped LAST -> now 2
    assert_equal(lw.read_next_experiment_number(), 2, "new bumped counter")


def test_new_invalid_type_error_and_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["new", "nonsense", "x"])
    assert_equal(result.exit_code, 1, "new invalid-type exit code")
    assert_contains(result.stderr_text, "Unknown type", "new invalid-type message")
    after_state = snapshot_store_state(lw)
    # no row, no folder, and crucially NO counter bump on a rejected new.
    assert_no_effect(before_state, after_state, "new invalid type")


def test_new_assigns_sequential_ids(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    run_command(lw, ["new", "purification", "a"])
    run_command(lw, ["new", "cloning", "b"])
    store = lw.load_store_from_database(lw.DATABASE_PATH)
    assigned_ids = sorted(row["id"] for row in store["experiments"])
    assert_equal(assigned_ids, ["LM-0001", "LM-0002"], "sequential ids")
    assert_equal(lw.read_next_experiment_number(), 3, "counter after two new")


# ---------------------------------------------------------------------------
# link -- happy path + FOUR mandatory distinct no-effect rejection paths
# ---------------------------------------------------------------------------


def test_link_happy_path_normalizes_ids(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    # bare numbers normalize to LM-0001 / LM-0002
    result = run_command(lw, ["link", "1", "2", "follow_up_to"])
    assert_equal(result.exit_code, 0, "link exit code")
    store = lw.load_store_from_database(lw.DATABASE_PATH)
    assert_equal(len(store["experiment_links"]), 1, "link created one edge")
    edge = store["experiment_links"][0]
    assert_equal(edge["from_experiment"], "LM-0001", "link from normalized")
    assert_equal(edge["to_experiment"], "LM-0002", "link to normalized")
    assert_equal(edge["relationship"], "follow_up_to", "link relationship")


def test_link_rejection_bad_id_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["link", "banana", "2", "follow_up_to"])
    assert_equal(result.exit_code, 1, "link bad-id exit code")
    assert_contains(result.stderr_text, "Invalid experiment id", "link bad-id message")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "link bad id")


def test_link_rejection_missing_experiment_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["link", "1", "99", "follow_up_to"])
    assert_equal(result.exit_code, 1, "link missing-experiment exit code")
    assert_contains(result.stderr_text, "No such experiment", "link missing message")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "link missing experiment")


def test_link_rejection_unknown_relationship_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["link", "1", "2", "inspired_by"])
    assert_equal(result.exit_code, 1, "link unknown-relationship exit code")
    assert_contains(
        result.stderr_text,
        "Unknown relationship",
        "link unknown-relationship message",
    )
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "link unknown relationship")


def test_link_rejection_duplicate_edge_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    # establish the edge, THEN snapshot, THEN attempt the duplicate.
    run_command(lw, ["link", "1", "2", "follow_up_to"])
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["link", "1", "2", "follow_up_to"])
    assert_equal(result.exit_code, 1, "link duplicate exit code")
    assert_contains(result.stderr_text, "Link already exists", "link duplicate message")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "link duplicate edge")


def test_link_reverse_direction_is_not_a_duplicate(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    run_command(lw, ["link", "1", "2", "follow_up_to"])
    # (2 -> 1) is a distinct directed edge and must be allowed.
    result = run_command(lw, ["link", "2", "1", "replicate_of"])
    assert_equal(result.exit_code, 0, "reverse-direction link exit code")
    store = lw.load_store_from_database(lw.DATABASE_PATH)
    assert_equal(len(store["experiment_links"]), 2, "reverse direction added an edge")


# ---------------------------------------------------------------------------
# show
# ---------------------------------------------------------------------------


def test_show_happy_path_summary_and_full(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    run_command(lw, ["link", "1", "2", "follow_up_to"])

    store = lw.load_store_from_database(lw.DATABASE_PATH)
    folder_name = next(
        row["folder_name"] for row in store["experiments"] if row["id"] == "LM-0001"
    )
    experiment_folder = os.path.join(lw.EXPERIMENTS_DIRECTORY, folder_name)
    # date-front-loaded names: lexical sort == chronological
    for file_name in (
        "20260623_LM-0001_raw.csv",
        "20260625_LM-0001_notes.md",
        "20260624_LM-0001_gel.tiff",
    ):
        open(os.path.join(experiment_folder, file_name), "w").close()

    summary_result = run_command(lw, ["show", "1"])
    assert_equal(summary_result.exit_code, 0, "show summary exit code")
    assert_contains(
        summary_result.stdout_text,
        "Files: 3 (last: 20260625_LM-0001_notes.md)",
        "show summary newest file",
    )
    assert_contains(
        summary_result.stdout_text,
        "-> LM-0002 (follow_up_to)",
        "show renders outgoing link",
    )

    full_result = run_command(lw, ["show", "1", "--files"])
    assert_equal(full_result.exit_code, 0, "show --files exit code")
    assert_contains(full_result.stdout_text, "Files: 3", "show full count")
    assert_contains(
        full_result.stdout_text,
        "20260623_LM-0001_raw.csv",
        "show full lists each file",
    )


def test_show_missing_folder_is_rendered_not_raised(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    store = lw.load_store_from_database(lw.DATABASE_PATH)
    folder_name = next(
        row["folder_name"] for row in store["experiments"] if row["id"] == "LM-0001"
    )
    # user hand-deletes the folder: a NORMAL state, not an error.
    os.rmdir(os.path.join(lw.EXPERIMENTS_DIRECTORY, folder_name))
    result = run_command(lw, ["show", "1"])
    assert_equal(result.exit_code, 0, "missing folder is not an error exit")
    assert_contains(
        result.stdout_text,
        "Files: folder not found",
        "missing folder rendered line",
    )


def test_show_missing_experiment_error_and_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["show", "99"])
    assert_equal(result.exit_code, 1, "show missing exit code")
    assert_contains(result.stderr_text, "No such experiment", "show missing message")
    after_state = snapshot_store_state(lw)
    # show is read-only; a failed show must change nothing.
    assert_no_effect(before_state, after_state, "show missing experiment")


# ---------------------------------------------------------------------------
# list
# ---------------------------------------------------------------------------


def test_list_happy_path_and_filters(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)

    no_filter_result = run_command(lw, ["list"])
    assert_equal(no_filter_result.exit_code, 0, "list exit code")
    assert_contains(no_filter_result.stdout_text, "LM-0001", "list shows 0001")
    assert_contains(no_filter_result.stdout_text, "LM-0002", "list shows 0002")

    type_filter_result = run_command(lw, ["list", "--type", "cloning"])
    assert_contains(
        type_filter_result.stdout_text, "LM-0002", "type filter keeps cloning"
    )
    assert_true(
        "LM-0001" not in type_filter_result.stdout_text,
        "type filter drops purification",
    )

    status_filter_result = run_command(lw, ["list", "--status", "active"])
    assert_contains(status_filter_result.stdout_text, "LM-0001", "status filter active")


def test_list_no_match_is_normal_state_and_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["list", "--type", "purification", "--status", "done"])
    assert_equal(result.exit_code, 0, "list no-match exit code")
    assert_contains(result.stdout_text, "No experiments match", "list no-match message")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "list no-match")


# ---------------------------------------------------------------------------
# Hardening -- counter file corruption + counter/DB disagreement
# ---------------------------------------------------------------------------


def test_new_with_corrupt_counter_fails_cleanly_and_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    # hand-corrupt the counter as a user editing the filesystem might.
    with open(lw.COUNTER_FILE_PATH, "w") as counter_file:
        counter_file.write("not-a-number")
    # snapshot the DB + folders (cannot snapshot counter -- read raises).
    store_before = lw.load_store_from_database(lw.DATABASE_PATH)
    folders_before = sorted(os.listdir(lw.EXPERIMENTS_DIRECTORY))

    result = run_command(lw, ["new", "purification", "x"])
    assert_equal(result.exit_code, 1, "corrupt-counter exit code")
    assert_contains(result.stderr_text, "corrupt", "corrupt-counter message")
    assert_contains(result.stderr_text, "Fix it by hand", "corrupt-counter remedy hint")
    store_after = lw.load_store_from_database(lw.DATABASE_PATH)
    folders_after = sorted(os.listdir(lw.EXPERIMENTS_DIRECTORY))
    assert_equal(
        store_after["experiments"],
        store_before["experiments"],
        "corrupt counter / no row written",
    )
    assert_equal(folders_after, folders_before, "corrupt counter / no folder created")


def test_new_with_empty_counter_fails_cleanly(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    with open(lw.COUNTER_FILE_PATH, "w") as counter_file:
        counter_file.write("")
    result = run_command(lw, ["new", "purification", "x"])
    assert_equal(result.exit_code, 1, "empty-counter exit code")
    assert_contains(result.stderr_text, "empty", "empty-counter message")


def test_read_counter_missing_file_raises_counter_file_error(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    os.remove(lw.COUNTER_FILE_PATH)
    raised = False
    try:
        lw.read_next_experiment_number()
    except lw.CounterFileError as counter_error:
        raised = True
        assert_contains(str(counter_error), "not found", "missing-counter message")
    assert_true(raised, "missing counter file must raise CounterFileError")


def test_new_counter_db_disagreement_blocks_collision_no_effect(sandbox_lab_root):
    # The SS6.6 failure mode: counter points at an id the DB already has.
    lw = fresh_initialized_lw(sandbox_lab_root)
    run_command(lw, ["new", "purification", "first"])  # LM-0001, counter -> 2
    # hand-rewind the counter so it now points back at the existing LM-0001.
    lw.write_next_experiment_number(1)
    before_state = snapshot_store_state(lw)

    result = run_command(lw, ["new", "cloning", "second"])
    assert_equal(result.exit_code, 1, "disagreement exit code")
    assert_contains(result.stderr_text, "disagree", "disagreement message")
    assert_contains(result.stderr_text, "above the highest", "disagreement remedy hint")
    after_state = snapshot_store_state(lw)
    # crucial: the existing LM-0001 row is NOT clobbered by the rewrite, and
    # nothing new is written.
    assert_no_effect(before_state, after_state, "counter/DB disagreement")


def test_experiment_id_constants_are_coupled(sandbox_lab_root):
    # Guards the CRITICAL INVARIANT: build and normalize agree because both
    # read the same two constants. If a future edit hardcodes a width or
    # prefix at one site, this round-trip breaks.
    lw = load_lw_into_sandbox(sandbox_lab_root)
    built_id = f"{lw.EXPERIMENT_ID_PREFIX}{7:0{lw.EXPERIMENT_ID_NUMBER_WIDTH}d}"
    assert_equal(
        lw.normalize_experiment_id("7"),
        built_id,
        "normalize agrees with the built id format",
    )
    assert_equal(
        lw.normalize_experiment_id(built_id),
        built_id,
        "normalize is idempotent on a built id",
    )


# ---------------------------------------------------------------------------
# Dispatch usage messages -- condition stated, usage shown, next step given
# ---------------------------------------------------------------------------


def test_no_command_shows_full_overview(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, [])
    assert_equal(result.exit_code, 2, "no-command exit code")
    assert_contains(result.stderr_text, "No command given", "no-command condition")
    # full overview lists every command with its purpose line.
    for command_name in ("init", "new", "link", "show", "list"):
        assert_contains(
            result.stderr_text, command_name, f"overview lists {command_name}"
        )


def test_unknown_command_names_it_and_points_to_overview(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, ["frobnicate"])
    assert_equal(result.exit_code, 2, "unknown-command exit code")
    assert_contains(result.stderr_text, "frobnicate", "unknown-command names the input")


def test_new_too_few_args_lists_valid_types(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, ["new", "purification"])
    assert_equal(result.exit_code, 2, "new-too-few exit code")
    assert_contains(result.stderr_text, "1 given", "new-too-few states count")
    assert_contains(result.stderr_text, "purification", "new-too-few lists valid types")


def test_link_wrong_arg_count_lists_valid_relationships(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, ["link", "1", "2"])
    assert_equal(result.exit_code, 2, "link-arity exit code")
    assert_contains(result.stderr_text, "2 given", "link-arity states count")
    assert_contains(
        result.stderr_text, "follow_up_to", "link-arity lists relationships"
    )


def test_show_missing_id_points_to_list(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, ["show"])
    assert_equal(result.exit_code, 2, "show-missing-id exit code")
    assert_contains(
        result.stderr_text, "needs an experiment", "show-missing-id condition"
    )
    assert_contains(result.stderr_text, "lw list", "show-missing-id next step")


def test_show_extra_positional_suggests_files_flag(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, ["show", "1", "2"])
    assert_equal(result.exit_code, 2, "show-extra exit code")
    assert_contains(result.stderr_text, "--files", "show-extra suggests flag")


def test_list_stray_positional_is_rejected_with_flag_hint(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    result = run_command(lw, ["list", "active"])
    assert_equal(result.exit_code, 2, "list-stray exit code")
    assert_contains(result.stderr_text, "no positional", "list-stray condition")
    assert_contains(result.stderr_text, "--status", "list-stray flag hint")


def test_init_rejects_extra_arguments_and_suggests_new(sandbox_lab_root):
    # The real bug: `lw init cloning thing --title X` silently ran bare init.
    lw = load_lw_into_sandbox(sandbox_lab_root)
    result = run_command(lw, ["init", "cloning", "thing", "--title", "Some Title"])
    assert_equal(result.exit_code, 2, "init-extra-args exit code")
    assert_contains(
        result.stderr_text, "takes no arguments", "init-extra-args condition"
    )
    # the hint reconstructs the likely-intended `new` command.
    assert_contains(
        result.stderr_text, "lw new cloning thing", "init-extra-args suggests new"
    )


def test_new_rejects_extra_positionals(sandbox_lab_root):
    # multi-word title typed without --title -> caught, not silently dropped.
    lw = fresh_initialized_lw(sandbox_lab_root)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["new", "cloning", "thing", "Some", "Long", "Title"])
    assert_equal(result.exit_code, 2, "new-extra-positionals exit code")
    assert_contains(result.stderr_text, "positional arguments", "new-extra condition")
    assert_contains(result.stderr_text, "--title", "new-extra suggests --title")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "new extra positionals")


def test_new_rejects_unknown_flag(sandbox_lab_root):
    # a typo'd flag like --titel must not be silently swallowed.
    lw = fresh_initialized_lw(sandbox_lab_root)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["new", "cloning", "thing", "--titel", "Oops"])
    assert_equal(result.exit_code, 2, "new-unknown-flag exit code")
    assert_contains(result.stderr_text, "unexpected flag", "new-unknown-flag condition")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "new unknown flag")


def test_init_bare_still_works(sandbox_lab_root):
    # regression guard: the strictness must not break the normal init.
    lw = load_lw_into_sandbox(sandbox_lab_root)
    result = run_command(lw, ["init"])
    assert_equal(result.exit_code, 0, "bare init still works")
    assert_true(os.path.exists(lw.DATABASE_PATH), "bare init created database")


# ---------------------------------------------------------------------------
# complete -- mark done, built on the status-transition engine
# ---------------------------------------------------------------------------


def test_complete_happy_path_sets_status_and_date(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    result = run_command(lw, ["complete", "1"])
    assert_equal(result.exit_code, 0, "complete exit code")
    assert_contains(result.stdout_text, "complete", "complete confirmation")

    store = lw.load_store_from_database(lw.DATABASE_PATH)
    completed_row = next(row for row in store["experiments"] if row["id"] == "LM-0001")
    assert_equal(
        completed_row["status"],
        lw.COMPLETE_EXPERIMENT_STATUS,
        "complete set status",
    )
    assert_true(completed_row["date_completed"], "complete stamped date_completed")
    # the other experiment is untouched.
    other_row = next(row for row in store["experiments"] if row["id"] == "LM-0002")
    assert_equal(other_row["status"], "active", "complete left others active")


def test_complete_then_list_filter_reflects_status(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    run_command(lw, ["complete", "1"])
    completed_listing = run_command(lw, ["list", "--status", "complete"])
    assert_contains(
        completed_listing.stdout_text, "LM-0001", "completed shows in filter"
    )
    assert_true(
        "LM-0002" not in completed_listing.stdout_text,
        "active experiment absent from complete filter",
    )


def test_complete_already_complete_is_no_effect_error(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    run_command(lw, ["complete", "1"])
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["complete", "1"])
    assert_equal(result.exit_code, 1, "already-complete exit code")
    assert_contains(result.stderr_text, "already", "already-complete message")
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "complete already complete")


def test_complete_missing_experiment_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["complete", "99"])
    assert_equal(result.exit_code, 1, "complete-missing exit code")
    assert_contains(
        result.stderr_text, "No such experiment", "complete-missing message"
    )
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "complete missing experiment")


def test_complete_bad_id_no_effect(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    before_state = snapshot_store_state(lw)
    result = run_command(lw, ["complete", "banana"])
    assert_equal(result.exit_code, 1, "complete-bad-id exit code")
    assert_contains(
        result.stderr_text, "Invalid experiment id", "complete-bad-id message"
    )
    after_state = snapshot_store_state(lw)
    assert_no_effect(before_state, after_state, "complete bad id")


def test_complete_missing_id_argument_points_to_list(sandbox_lab_root):
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    result = run_command(lw, ["complete"])
    assert_equal(result.exit_code, 2, "complete-no-id exit code")
    assert_contains(
        result.stderr_text, "needs an experiment", "complete-no-id condition"
    )


def test_status_transition_engine_supports_reopen(sandbox_lab_root):
    # The engine is general by design: a future `reopen` is a wrapper that
    # targets the active status and clears the completion date. Exercise
    # that path directly to prove the core already supports it -- so adding
    # the verb later is a pure addition, not a redesign.
    lw = fresh_initialized_lw(sandbox_lab_root)
    seed_two_experiments(lw)
    run_command(lw, ["complete", "1"])

    store = lw.load_store_from_database(lw.DATABASE_PATH)
    reopen_effects = lw.build_status_transition_effects(
        store,
        {
            "id": "1",
            "target_status": lw.DEFAULT_EXPERIMENT_STATUS,
            "set_completion_date": "clear",
        },
        lw.startup_clock,
    )
    assert_true(reopen_effects["ok"], "engine reopen succeeds")
    reopened_row = next(
        row for row in reopen_effects["store"]["experiments"] if row["id"] == "LM-0001"
    )
    assert_equal(
        reopened_row["status"],
        lw.DEFAULT_EXPERIMENT_STATUS,
        "engine reopen restored active status",
    )
    assert_equal(
        reopened_row["date_completed"],
        None,
        "engine reopen cleared completion date",
    )


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

ALL_TESTS = [
    test_init_happy_path,
    test_init_idempotent_does_not_reset_counter,
    test_command_before_init_fails_precondition,
    test_new_happy_path,
    test_new_invalid_type_error_and_no_effect,
    test_new_assigns_sequential_ids,
    test_link_happy_path_normalizes_ids,
    test_link_rejection_bad_id_no_effect,
    test_link_rejection_missing_experiment_no_effect,
    test_link_rejection_unknown_relationship_no_effect,
    test_link_rejection_duplicate_edge_no_effect,
    test_link_reverse_direction_is_not_a_duplicate,
    test_show_happy_path_summary_and_full,
    test_show_missing_folder_is_rendered_not_raised,
    test_show_missing_experiment_error_and_no_effect,
    test_list_happy_path_and_filters,
    test_list_no_match_is_normal_state_and_no_effect,
    test_new_with_corrupt_counter_fails_cleanly_and_no_effect,
    test_new_with_empty_counter_fails_cleanly,
    test_read_counter_missing_file_raises_counter_file_error,
    test_new_counter_db_disagreement_blocks_collision_no_effect,
    test_experiment_id_constants_are_coupled,
    test_no_command_shows_full_overview,
    test_unknown_command_names_it_and_points_to_overview,
    test_new_too_few_args_lists_valid_types,
    test_link_wrong_arg_count_lists_valid_relationships,
    test_show_missing_id_points_to_list,
    test_show_extra_positional_suggests_files_flag,
    test_list_stray_positional_is_rejected_with_flag_hint,
    test_init_rejects_extra_arguments_and_suggests_new,
    test_new_rejects_extra_positionals,
    test_new_rejects_unknown_flag,
    test_init_bare_still_works,
    test_complete_happy_path_sets_status_and_date,
    test_complete_then_list_filter_reflects_status,
    test_complete_already_complete_is_no_effect_error,
    test_complete_missing_experiment_no_effect,
    test_complete_bad_id_no_effect,
    test_complete_missing_id_argument_points_to_list,
    test_status_transition_engine_supports_reopen,
]


def run_all_tests():
    passed_count = 0
    failed_tests = []
    for test_function in ALL_TESTS:
        with tempfile.TemporaryDirectory() as sandbox_lab_root:
            try:
                test_function(sandbox_lab_root)
                passed_count += 1
                print(f"  pass  {test_function.__name__}")
            except AssertionError as assertion_error:
                failed_tests.append((test_function.__name__, str(assertion_error)))
                print(f"  FAIL  {test_function.__name__}: {assertion_error}")
            except Exception as unexpected_error:  # noqa: BLE001
                failed_tests.append(
                    (test_function.__name__, f"ERROR: {unexpected_error!r}")
                )
                print(f"  ERROR {test_function.__name__}: {unexpected_error!r}")

    print()
    print(f"{passed_count}/{len(ALL_TESTS)} passed")
    if failed_tests:
        print("FAILURES:")
        for test_name, failure_detail in failed_tests:
            print(f"  {test_name}: {failure_detail}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(run_all_tests())
