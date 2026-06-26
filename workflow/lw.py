#!/usr/bin/env python3
"""lw -- personal lab experiment tracker.

Owns experiment identity and the link graph; nothing else. Experiment
details (yields, assay conditions, sample layouts, strains) live in Excel
sheets the user maintains by hand. Two tables, five commands (plus init).
"""

import datetime
import os
import sqlite3
import sys


# ---------------------------------------------------------------------------
# C1 -- data definitions + config + paths
# ---------------------------------------------------------------------------

# --- Explicit declarations: every value that was previously a magic literal
# --- is surfaced here as a named constant, so the coupled invariants below
# --- are single-sourced and a future reader can change them in one place.

# Experiment id format. CRITICAL INVARIANT: EXPERIMENT_ID_PREFIX and
# EXPERIMENT_ID_NUMBER_WIDTH are coupled. The prefix is matched (case
# -insensitively) when normalizing a raw id, and the width is the zero-pad
# applied when building an id from the counter. Both the build site
# (build_new_experiment_effects) and the normalize site
# (normalize_experiment_id) read these constants -- do NOT hardcode "LM-" or
# 4 anywhere else, or the two sites can silently drift apart.
EXPERIMENT_ID_PREFIX = "LM-"
EXPERIMENT_ID_NUMBER_WIDTH = 4

# Date formats. CRITICAL INVARIANT: FOLDER_DATE_FORMAT must remain
# date-front-loaded (YYYYMMDD, most-significant-first) because `show` relies
# on lexical sort equalling chronological order to find the newest file
# WITHOUT stat calls (SS3). Changing this format to anything not
# lexically-sortable-as-chronological silently breaks the "newest file"
# logic. ISO_DATE_FORMAT is the stored date_started value.
FOLDER_DATE_FORMAT = "%Y%m%d"
ISO_DATE_FORMAT = "%Y-%m-%d"

# Experiment row defaults and status values. DEFAULT_EXPERIMENT_STATUS is
# what `new` writes; COMPLETE_EXPERIMENT_STATUS is what `complete` writes.
# The status-transition engine reads these so the two ends of every
# transition (active <-> complete) are single-sourced.
DEFAULT_EXPERIMENT_STATUS = "active"
COMPLETE_EXPERIMENT_STATUS = "complete"

# Counter seed: the first experiment number a fresh `init` writes.
INITIAL_COUNTER_VALUE = 1

# `list` output column widths (status, type) for aligned columns.
LIST_STATUS_COLUMN_WIDTH = 8
LIST_TYPE_COLUMN_WIDTH = 13

# Filenames within LAB_ROOT (joined to LAB_ROOT to derive the paths below).
DATABASE_FILENAME = "lab.db"
COUNTER_FILENAME = "counter.txt"
COMMAND_LOG_FILENAME = "command_log.txt"

# Environment variables consulted for LAB_ROOT resolution, in precedence
# order. The second of each pair is a legacy alias kept because it is free.
LAB_ROOT_ENV_VARS = ("LW_LAB_ROOT", "LAB_ROOT")
WINDOWS_USER_ENV_VARS = ("LW_WINDOWS_USER", "MC_WINDOWS_USER")

# Folder-naming category codes. NOTE: `type` is NOT a data contract. The
# typed data it once dispatched to (purification_details, assay_details,
# samples) now lives in Excel. `type` is retained for exactly two uses:
# (a) it supplies the short code in the folder name (here), and (b) it
# filters `list`. Do not mistake `type` for a contract pointing at
# DB-resident data.
EXPERIMENT_TYPE_TO_CATEGORY_CODE = {
    "purification": "pur",
    "loading": "load",
    "gelshift": "gs",
    "atpase": "atp",
    "genetics": "gen",
    "computational": "comp",
    "cloning": "clon",
    "sequencing": "seq",
}
VALID_EXPERIMENT_TYPES = list(EXPERIMENT_TYPE_TO_CATEGORY_CODE.keys())

# Closed set of edge labels. The name encodes block-on-absence semantics:
# `link` REJECTS any relationship not in this list (blocks, does not warn).
# Extending is a deliberate one-line edit. All three are now plain edge
# labels -- `uses_prep_from` no longer has special behavior.
ALLOWED_LINK_RELATIONSHIPS = ["uses_prep_from", "replicate_of", "follow_up_to"]

# Captured once at startup, passed as the "clock" to transforms so the
# transforms never read the clock themselves.
startup_clock = {
    "today": datetime.date.today(),
    "now": datetime.datetime.now(),
}

# --- LAB_ROOT resolution ---
# Resolve in precedence order: an explicit LAB_ROOT env var wins; otherwise
# a Windows username env var derives the WSL desktop convention path;
# otherwise we cannot proceed and exit with a hint naming both options.


def first_environment_value(environment_variable_names):
    """Return the first set, non-empty environment variable value among the
    given names, or None if none is set."""
    for environment_variable_name in environment_variable_names:
        value = os.environ.get(environment_variable_name)
        if value:
            return value
    return None


lab_root_from_environment = first_environment_value(LAB_ROOT_ENV_VARS)
if lab_root_from_environment:
    LAB_ROOT = lab_root_from_environment
else:
    windows_username_from_environment = first_environment_value(WINDOWS_USER_ENV_VARS)
    if windows_username_from_environment:
        LAB_ROOT = f"/mnt/c/Users/{windows_username_from_environment}/Desktop/lab"
    else:
        print(
            "Error: lab root is not configured.\n"
            f"  Set one of {', '.join(LAB_ROOT_ENV_VARS)} to an absolute "
            "lab directory path,\n"
            f"  or set one of {', '.join(WINDOWS_USER_ENV_VARS)} to your "
            "Windows username\n"
            "  to use the /mnt/c/Users/<user>/Desktop/lab WSL convention.",
            file=sys.stderr,
        )
        sys.exit(2)

# --- derived paths ---
DATABASE_PATH = os.path.join(LAB_ROOT, DATABASE_FILENAME)
EXPERIMENTS_DIRECTORY = os.path.join(LAB_ROOT, "experiments")
COUNTER_FILE_PATH = os.path.join(LAB_ROOT, COUNTER_FILENAME)
COMMAND_LOG_PATH = os.path.join(LAB_ROOT, COMMAND_LOG_FILENAME)


# ---------------------------------------------------------------------------
# C2 -- effects-as-data primitives (the core contract)
# ---------------------------------------------------------------------------


def make_success_effects(standard_output_lines=None, store=None, exit_code=0):
    return {
        "ok": True,
        "stdout": standard_output_lines or [],
        "stderr": [],
        "store": store,
        "exit_code": exit_code,
    }


def make_failure_effects(error_message, exit_code=1):
    if isinstance(error_message, str):
        error_message = [error_message]
    return {
        "ok": False,
        "stdout": [],
        "stderr": error_message,
        "store": None,
        "exit_code": exit_code,
    }


def execute_effects(effects, database_path):
    """The ONLY function that performs IO for transformed commands.

    The stdout writes are guarded against BrokenPipeError so that piping a
    command into a reader that closes early (e.g. `lw show 3 | head`) exits
    cleanly instead of dumping a traceback. The DB commit still runs: the
    consumer closing the pipe does not mean the effect should be abandoned.
    """
    try:
        for output_line in effects.get("stdout", []):
            print(output_line)
    except BrokenPipeError:
        # Downstream reader closed the pipe. Redirect remaining stdout to
        # devnull so interpreter shutdown does not re-raise, then continue.
        devnull_file_descriptor = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull_file_descriptor, sys.stdout.fileno())
    for error_line in effects.get("stderr", []):
        print(error_line, file=sys.stderr)
    if effects.get("store") is not None:
        write_store_to_database(effects["store"], database_path)
    return effects.get("exit_code", 0)


# ---------------------------------------------------------------------------
# C3 -- id normalization
# ---------------------------------------------------------------------------


def normalize_experiment_id(raw_experiment_id):
    """'LM-0003' or '3' -> 'LM-0003'; None if invalid.

    Reads EXPERIMENT_ID_PREFIX and EXPERIMENT_ID_NUMBER_WIDTH so this site
    cannot drift from the id-building site in build_new_experiment_effects.
    """
    if raw_experiment_id.upper().startswith(EXPERIMENT_ID_PREFIX):
        return raw_experiment_id.upper()
    try:
        return (
            f"{EXPERIMENT_ID_PREFIX}"
            f"{int(raw_experiment_id):0{EXPERIMENT_ID_NUMBER_WIDTH}d}"
        )
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Store -- the ONLY place that touches SQL. One module owns the schema, the
# connection, load, and commit. Transforms receive plain dicts and return
# effects; they never touch this layer.
# ---------------------------------------------------------------------------

# C4 -- schema + load
# ---------------------------------------------------------------------------

# Live schema is exactly the two tables defined below.
#
# schema_notes: this DB file also contains eight ORPHANED tables --
#   strains, experiment_strains, files, purification_details,
#   assay_details, figure_panels, samples, stage_expectations
# -- from retired subsystems now migrated to manual Excel tracking. They are
# INTENTIONALLY RETIRED: no longer read or written by this code, and
# deliberately NOT dropped (leaving them untouched is zero-risk and
# reversible; a fresh DB is neither). A future reader seeing eight extra
# tables in the .db is seeing intended state, not a bug. The live schema is
# the two tables defined here.

DATABASE_SCHEMA = """
CREATE TABLE IF NOT EXISTS experiments (
    id TEXT PRIMARY KEY,
    type TEXT NOT NULL,
    title TEXT NOT NULL,
    shortname TEXT NOT NULL,
    date_started TEXT NOT NULL,
    date_completed TEXT,
    status TEXT DEFAULT 'active',
    folder_name TEXT NOT NULL,
    notes TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS experiment_links (
    from_experiment TEXT REFERENCES experiments(id),
    to_experiment   TEXT REFERENCES experiments(id),
    relationship    TEXT NOT NULL,
    PRIMARY KEY (from_experiment, to_experiment)
);
"""

# The two live tables, in the order load/commit walk them.
LIVE_TABLE_NAMES = ["experiments", "experiment_links"]


def load_store_from_database(database_path):
    """Read every live table into a dict-of-lists of row dicts.

    Reads ONLY the two live tables. The orphaned tables sit inert.
    """
    database_connection = sqlite3.connect(database_path)
    database_connection.row_factory = sqlite3.Row
    store = {}
    try:
        for table_name in LIVE_TABLE_NAMES:
            table_cursor = database_connection.execute(f"SELECT * FROM {table_name}")
            store[table_name] = [
                dict(table_row) for table_row in table_cursor.fetchall()
            ]
    finally:
        database_connection.close()
    return store


# C5 -- commit (whole-store rewrite + union-of-keys fix)
# ---------------------------------------------------------------------------


def write_store_to_database(store, database_path):
    """Rewrite the whole store in one transaction.

    This whole-store-rewrite is DELIBERATELY RETAINED. It works for a
    single-user, few-hundred-row DB. Do NOT replace with per-statement
    effects unless a real problem is reported (row count, a second writer).

    The column list per table is built from the UNION of all rows' keys,
    and values are read with experiment_row.get(column_name) rather than
    experiment_row[column_name] -- this is the fix for the KeyError the
    naive version raised on a missing dict key. Touches only the live
    tables that load_store_from_database loaded; orphaned tables are never
    read or written here.
    """
    database_connection = sqlite3.connect(database_path)
    try:
        database_connection.executescript(DATABASE_SCHEMA)
        database_cursor = database_connection.cursor()
        for table_name in LIVE_TABLE_NAMES:
            table_rows = store.get(table_name, [])
            database_cursor.execute(f"DELETE FROM {table_name}")
            if not table_rows:
                continue

            # Build the column list from the UNION of all rows' keys, in
            # first-seen order, so a ragged row missing a key cannot break
            # the INSERT.
            column_names = []
            seen_column_names = set()
            for table_row in table_rows:
                for column_name in table_row:
                    if column_name not in seen_column_names:
                        seen_column_names.add(column_name)
                        column_names.append(column_name)

            placeholder_list = ", ".join(["?"] * len(column_names))
            insert_statement = (
                f"INSERT INTO {table_name} "
                f"({', '.join(column_names)}) VALUES ({placeholder_list})"
            )
            for table_row in table_rows:
                row_values = [
                    table_row.get(column_name) for column_name in column_names
                ]
                database_cursor.execute(insert_statement, row_values)
        database_connection.commit()
    finally:
        database_connection.close()


# C6 -- missing-DB precondition
# ---------------------------------------------------------------------------


def require_existing_database(database_path):
    """Store-boundary precondition for every command except init.

    Returns None if the DB exists; otherwise a failure-effects dict
    describing the clean error. This runs BEFORE any transform --
    transforms assume a loaded store. Distinct from a missing experiment
    FOLDER on show, which is a normal rendered state, not a precondition
    error.
    """
    if not os.path.exists(database_path):
        return make_failure_effects("No lab database found. Run 'lw init' first.")
    return None


# C7 -- init (the bootstrap outsider)
# ---------------------------------------------------------------------------


def run_init_command():
    """Create LAB_ROOT/experiments, the DB, and seed counter.txt if absent.

    The outsider: runs before the store exists, so it cannot use the store
    machinery. Idempotent -- CREATE TABLE IF NOT EXISTS and the counter
    seed both no-op on re-run. Returns an exit code directly (it is its own
    edge, not a transform routed through execute_effects).
    """
    os.makedirs(LAB_ROOT, exist_ok=True)
    os.makedirs(EXPERIMENTS_DIRECTORY, exist_ok=True)

    database_connection = sqlite3.connect(DATABASE_PATH)
    try:
        database_connection.executescript(DATABASE_SCHEMA)
        database_connection.commit()
    finally:
        database_connection.close()

    if not os.path.exists(COUNTER_FILE_PATH):
        with open(COUNTER_FILE_PATH, "w") as counter_file:
            counter_file.write(str(INITIAL_COUNTER_VALUE))

    print(f"Initialized lab at {LAB_ROOT}")
    return 0


# ---------------------------------------------------------------------------
# new -- the #1 feature, decomposed: counter edge (C8) + pure transform (C9)
# + edge/wiring (C10). The one place transform purity, the counter-as-effect
# boundary, fixed bump-last ordering, and filesystem IO all intersect.
# ---------------------------------------------------------------------------

# C8 -- counter edge (read + bump), ordering-sensitive
# ---------------------------------------------------------------------------


class CounterFileError(Exception):
    """Raised when counter.txt is missing, empty, or not a plain integer.

    The counter lives on a filesystem the user hand-edits (SS3), so a
    corrupt counter is a foreseeable state, not an impossible one. The edge
    raises this with an actionable message rather than letting a bare
    ValueError/FileNotFoundError surface as an opaque traceback.
    """


def read_next_experiment_number():
    """Read the current 'next experiment number' from counter.txt.

    Filesystem edge -- never called from inside a transform. The value is
    passed INTO the new transform via args; the transform stays pure.

    Presupposition being guarded: counter.txt exists and holds a single
    plain integer. Both can be violated by hand-editing, so a missing,
    empty, or non-integer file raises CounterFileError with a remedy.
    """
    try:
        with open(COUNTER_FILE_PATH) as counter_file:
            raw_counter_contents = counter_file.read().strip()
    except FileNotFoundError:
        raise CounterFileError(
            f"Counter file not found at {COUNTER_FILE_PATH}. "
            "Run 'lw init' to create it."
        )

    if not raw_counter_contents:
        raise CounterFileError(
            f"Counter file {COUNTER_FILE_PATH} is empty; expected a plain "
            f"integer. Set it to the next experiment number (e.g. "
            f"{INITIAL_COUNTER_VALUE}) by hand and retry."
        )
    try:
        return int(raw_counter_contents)
    except ValueError:
        raise CounterFileError(
            f"Counter file {COUNTER_FILE_PATH} is corrupt: expected a plain "
            f"integer, got {raw_counter_contents!r}. Fix it by hand "
            "(set it above the highest existing LM- id) and retry."
        )


def write_next_experiment_number(next_experiment_number):
    """Write the bumped counter value to counter.txt.

    Filesystem edge. Called LAST in the new sequence (after row insert and
    mkdir), so a run that dies mid-way wastes a number (harmless) rather
    than reusing one (a PK collision -- a real failure). Bumping last fails
    toward the harmless direction. Honor the ordering.
    """
    with open(COUNTER_FILE_PATH, "w") as counter_file:
        counter_file.write(str(next_experiment_number))


# C9 -- new pure transform
# ---------------------------------------------------------------------------


def build_new_experiment_effects(store, arguments, clock):
    """(store, arguments, clock) -> effects. Pure: no IO, no clock-read,
    no print.

    arguments carries: type, shortname, title (optional), and the
    experiment number read at the edge (C8). Builds experiment_id /
    folder_name purely, produces the experiment row, and returns the next
    experiment number as an extra 'next_experiment_number' key on the
    effects dict for the new edge (C10) to write. The fixed success-effects
    shape is unchanged; the extra key rides alongside it.
    """
    experiment_type = arguments["type"]
    experiment_shortname = arguments["shortname"]
    experiment_title = arguments.get("title") or experiment_shortname
    current_experiment_number = arguments["experiment_number"]

    if experiment_type not in VALID_EXPERIMENT_TYPES:
        return make_failure_effects(
            f"Unknown type '{experiment_type}'. "
            f"Valid types: {', '.join(VALID_EXPERIMENT_TYPES)}"
        )

    experiment_id = (
        f"{EXPERIMENT_ID_PREFIX}"
        f"{current_experiment_number:0{EXPERIMENT_ID_NUMBER_WIDTH}d}"
    )

    # Collision precondition (the SS6.6 counter/DB-disagreement failure
    # mode). CRITICAL: the whole-store rewrite (write_store_to_database)
    # does DELETE-then-INSERT, so a duplicate id would NOT raise a PK
    # IntegrityError -- it would silently clobber the table and re-insert,
    # losing the older row. The collision must therefore be caught HERE, as
    # an explicit store-based precondition, not as a caught DB exception.
    existing_experiment_ids = {
        existing_row["id"] for existing_row in store["experiments"]
    }
    if experiment_id in existing_experiment_ids:
        return make_failure_effects(
            f"Experiment {experiment_id} already exists, but the counter "
            f"points at it -- counter.txt and the database disagree. Set "
            f"{COUNTER_FILE_PATH} above the highest existing "
            f"{EXPERIMENT_ID_PREFIX} id and retry."
        )

    compact_date_string = clock["today"].strftime(FOLDER_DATE_FORMAT)
    category_code = EXPERIMENT_TYPE_TO_CATEGORY_CODE[experiment_type]
    folder_name = (
        f"{compact_date_string}_{experiment_id}_{category_code}_{experiment_shortname}"
    )
    iso_date_string = clock["today"].strftime(ISO_DATE_FORMAT)

    experiment_row = {
        "id": experiment_id,
        "type": experiment_type,
        "title": experiment_title,
        "shortname": experiment_shortname,
        "date_started": iso_date_string,
        "date_completed": None,
        "status": DEFAULT_EXPERIMENT_STATUS,
        "folder_name": folder_name,
        "notes": None,
    }

    updated_store = dict(store)
    updated_store["experiments"] = list(store["experiments"]) + [experiment_row]

    effects = make_success_effects(
        standard_output_lines=[f"Created {experiment_id}: {folder_name}"],
        store=updated_store,
    )
    # Extra keys on top of the fixed success-effects shape: the next
    # experiment number the C10 edge writes LAST, and the folder name the
    # C10 edge needs for mkdir. Not part of execute_effects' duties.
    effects["next_experiment_number"] = current_experiment_number + 1
    effects["folder_name"] = folder_name
    return effects


# C10 -- new edge + wiring
# ---------------------------------------------------------------------------


def run_new_command(arguments):
    """Edge for `new`: read counter, run the pure transform, sequence the
    effects in the fixed fail-safe order, route through execute_effects.

    Order is deliberate (SS6.6): row insert (via execute_effects ->
    write_store_to_database) and mkdir happen FIRST; the counter bump
    happens LAST. A run that dies before the bump wastes a number
    (harmless); never reuses one.
    """
    try:
        current_experiment_number = read_next_experiment_number()
    except CounterFileError as counter_error:
        return execute_effects(make_failure_effects(str(counter_error)), DATABASE_PATH)
    arguments = dict(arguments)
    arguments["experiment_number"] = current_experiment_number

    effects = build_new_experiment_effects(
        load_store_from_database(DATABASE_PATH), arguments, startup_clock
    )

    if not effects["ok"]:
        # Rejected: no row, no mkdir, no bump. Just emit the failure.
        return execute_effects(effects, DATABASE_PATH)

    # FIRST: mkdir the experiment folder.
    experiment_folder_path = os.path.join(EXPERIMENTS_DIRECTORY, effects["folder_name"])
    os.makedirs(experiment_folder_path, exist_ok=True)

    # FIRST (cont.): commit the row + print stdout via execute_effects.
    exit_code = execute_effects(effects, DATABASE_PATH)

    # LAST: bump the counter, only after row + folder are durable.
    write_next_experiment_number(effects["next_experiment_number"])

    return exit_code


# ---------------------------------------------------------------------------
# link -- append a bare provenance edge. C11 transform + C12 wiring.
# ---------------------------------------------------------------------------

# C11 -- link transform (bare edge; the ~20-line version)
# ---------------------------------------------------------------------------


def build_link_effects(store, arguments, clock):
    """(store, arguments, clock) -> effects. Pure. Appends one bare edge.

    Four rejection paths, each a no-effect failure:
      1. bad id            -- either id fails normalize_experiment_id
      2. missing experiment-- either id not present in the store
      3. unknown relationship -- not in ALLOWED_LINK_RELATIONSHIPS
                                 (block, not warn)
      4. duplicate edge    -- (from, to) already present
    No assay_details, no type-specific behavior -- a plain edge label.
    """
    from_experiment_id = normalize_experiment_id(arguments["from_id"])
    to_experiment_id = normalize_experiment_id(arguments["to_id"])
    relationship = arguments["relationship"]

    # 1. bad id
    if from_experiment_id is None or to_experiment_id is None:
        if from_experiment_id is None:
            invalid_raw_id = arguments["from_id"]
        else:
            invalid_raw_id = arguments["to_id"]
        return make_failure_effects(f"Invalid experiment id: '{invalid_raw_id}'")

    # 2. missing experiment
    known_experiment_ids = {
        experiment_row["id"] for experiment_row in store["experiments"]
    }
    if from_experiment_id not in known_experiment_ids:
        return make_failure_effects(
            f"No such experiment: {from_experiment_id}. "
            "Run 'lw list' to see existing experiments."
        )
    if to_experiment_id not in known_experiment_ids:
        return make_failure_effects(
            f"No such experiment: {to_experiment_id}. "
            "Run 'lw list' to see existing experiments."
        )

    # 3. unknown relationship -- blocked, not warned
    if relationship not in ALLOWED_LINK_RELATIONSHIPS:
        return make_failure_effects(
            f"Unknown relationship '{relationship}'. "
            f"Allowed: {', '.join(ALLOWED_LINK_RELATIONSHIPS)}"
        )

    # 4. duplicate edge
    for existing_link in store["experiment_links"]:
        existing_pair_matches = (
            existing_link["from_experiment"] == from_experiment_id
            and existing_link["to_experiment"] == to_experiment_id
        )
        if existing_pair_matches:
            return make_failure_effects(
                f"Link already exists: {from_experiment_id} -> {to_experiment_id}"
            )

    new_link_edge = {
        "from_experiment": from_experiment_id,
        "to_experiment": to_experiment_id,
        "relationship": relationship,
    }
    updated_store = dict(store)
    updated_store["experiment_links"] = list(store["experiment_links"]) + [
        new_link_edge
    ]

    return make_success_effects(
        standard_output_lines=[
            f"Linked {from_experiment_id} -> {to_experiment_id} ({relationship})"
        ],
        store=updated_store,
    )


# C12 -- link wiring
# ---------------------------------------------------------------------------


def run_link_command(arguments):
    """Edge for `link`: run the pure transform, route through
    execute_effects."""
    effects = build_link_effects(
        load_store_from_database(DATABASE_PATH), arguments, startup_clock
    )
    return execute_effects(effects, DATABASE_PATH)


# ---------------------------------------------------------------------------
# show -- experiment row + links + folder file summary. C13.
# ---------------------------------------------------------------------------

# C13 -- show transform + listdir edge
# ---------------------------------------------------------------------------

# Sentinel the edge passes to the transform when the experiment folder is
# absent. Missing folder is a NORMAL state (the user hand-edits the
# filesystem), rendered not raised -- distinct from the missing-DB hard
# precondition (SS2).
EXPERIMENT_FOLDER_NOT_FOUND = object()


def build_show_effects(store, arguments, clock):
    """(store, arguments, clock) -> effects. Pure: receives the file
    listing (or the EXPERIMENT_FOLDER_NOT_FOUND sentinel) from the edge;
    never touches the FS.

    arguments carries: id, folder_file_names (list[str] |
    EXPERIMENT_FOLDER_NOT_FOUND), experiment_folder_path, show_full_listing
    (bool -- the --files flag).
    """
    experiment_id = normalize_experiment_id(arguments["id"])
    if experiment_id is None:
        return make_failure_effects(f"Invalid experiment id: '{arguments['id']}'")

    matching_experiment_rows = [
        experiment_row
        for experiment_row in store["experiments"]
        if experiment_row["id"] == experiment_id
    ]
    if not matching_experiment_rows:
        return make_failure_effects(
            f"No such experiment: {experiment_id}. "
            "Run 'lw list' to see existing experiments."
        )
    experiment_row = matching_experiment_rows[0]

    output_lines = []
    output_lines.append(f"{experiment_row['id']}: {experiment_row['title']}")
    output_lines.append(f"  type:      {experiment_row['type']}")
    output_lines.append(f"  shortname: {experiment_row['shortname']}")
    output_lines.append(f"  status:    {experiment_row['status']}")
    output_lines.append(f"  started:   {experiment_row['date_started']}")
    if experiment_row.get("date_completed"):
        output_lines.append(f"  completed: {experiment_row['date_completed']}")
    output_lines.append(f"  folder:    {experiment_row['folder_name']}")

    # Links: outgoing and incoming.
    outgoing_links = [
        link
        for link in store["experiment_links"]
        if link["from_experiment"] == experiment_id
    ]
    incoming_links = [
        link
        for link in store["experiment_links"]
        if link["to_experiment"] == experiment_id
    ]
    if outgoing_links:
        output_lines.append("  links out:")
        for link in outgoing_links:
            output_lines.append(
                f"    -> {link['to_experiment']} ({link['relationship']})"
            )
    if incoming_links:
        output_lines.append("  links in:")
        for link in incoming_links:
            output_lines.append(
                f"    <- {link['from_experiment']} ({link['relationship']})"
            )

    # Files: folder IS the record. Date-front-loaded names => lexical sort
    # equals chronological, so no stat calls.
    folder_file_names = arguments["folder_file_names"]
    if folder_file_names is EXPERIMENT_FOLDER_NOT_FOUND:
        output_lines.append(
            f"  Files: folder not found ({arguments['experiment_folder_path']})"
        )
    elif arguments.get("show_full_listing"):
        output_lines.append(f"  Files: {len(folder_file_names)}")
        for file_name in sorted(folder_file_names):
            output_lines.append(f"    {file_name}")
    else:
        if folder_file_names:
            newest_file_name = sorted(folder_file_names)[-1]
            output_lines.append(
                f"  Files: {len(folder_file_names)} (last: {newest_file_name})"
            )
        else:
            output_lines.append("  Files: 0")

    return make_success_effects(standard_output_lines=output_lines)


def run_show_command(arguments):
    """Edge for `show`: unconditional listdir at the FS edge, missing
    folder caught into a sentinel, then the pure transform formats the
    output.
    """
    store = load_store_from_database(DATABASE_PATH)
    experiment_id = normalize_experiment_id(arguments["id"])

    experiment_folder_path = ""
    folder_file_names = EXPERIMENT_FOLDER_NOT_FOUND
    if experiment_id is not None:
        matching_experiment_rows = [
            experiment_row
            for experiment_row in store["experiments"]
            if experiment_row["id"] == experiment_id
        ]
        if matching_experiment_rows:
            experiment_folder_path = os.path.join(
                EXPERIMENTS_DIRECTORY,
                matching_experiment_rows[0]["folder_name"],
            )
            try:
                folder_file_names = os.listdir(experiment_folder_path)
            except OSError:
                folder_file_names = EXPERIMENT_FOLDER_NOT_FOUND

    arguments = dict(arguments)
    arguments["folder_file_names"] = folder_file_names
    arguments["experiment_folder_path"] = experiment_folder_path
    effects = build_show_effects(store, arguments, startup_clock)
    return execute_effects(effects, DATABASE_PATH)


# ---------------------------------------------------------------------------
# list -- filter experiments by --type / --status. C14.
# ---------------------------------------------------------------------------

# C14 -- list transform + flag parse
# ---------------------------------------------------------------------------


def build_list_effects(store, arguments, clock):
    """(store, arguments, clock) -> effects. Pure. Absent flag = no filter
    on that field (an absent predicate is meaningful: it matches
    everything).
    """
    type_filter = arguments.get("type")
    status_filter = arguments.get("status")

    matching_experiment_rows = store["experiments"]
    if type_filter is not None:
        matching_experiment_rows = [
            experiment_row
            for experiment_row in matching_experiment_rows
            if experiment_row["type"] == type_filter
        ]
    if status_filter is not None:
        matching_experiment_rows = [
            experiment_row
            for experiment_row in matching_experiment_rows
            if experiment_row["status"] == status_filter
        ]

    if not matching_experiment_rows:
        return make_success_effects(standard_output_lines=["No experiments match."])

    output_lines = []
    sorted_experiment_rows = sorted(
        matching_experiment_rows,
        key=lambda experiment_row: experiment_row["id"],
    )
    for experiment_row in sorted_experiment_rows:
        output_lines.append(
            f"{experiment_row['id']}  "
            f"{experiment_row['status']:{LIST_STATUS_COLUMN_WIDTH}}  "
            f"{experiment_row['type']:{LIST_TYPE_COLUMN_WIDTH}}  "
            f"{experiment_row['title']}"
        )
    return make_success_effects(standard_output_lines=output_lines)


def run_list_command(arguments):
    """Edge for `list`: run the pure transform, route through
    execute_effects."""
    effects = build_list_effects(
        load_store_from_database(DATABASE_PATH), arguments, startup_clock
    )
    return execute_effects(effects, DATABASE_PATH)


# ---------------------------------------------------------------------------
# complete -- mark an experiment done. Built on a general status-transition
# engine so a future `reopen` (or a collapse into one `status` command) is a
# pure addition over an unchanged, already-tested core (see BUILD_SUMMARY).
# ---------------------------------------------------------------------------


def build_status_transition_effects(store, arguments, clock):
    """(store, arguments, clock) -> effects. Pure. The shared core of every
    status change: set one experiment's status, and set or clear its
    completion date.

    arguments carries:
      id                -- the experiment to transition (raw or normalized)
      target_status     -- the status to move the row INTO
      set_completion_date -- one of:
          "today"  : stamp date_completed with clock["today"]
          "clear"  : set date_completed back to None
          (a future "--on" flag would add an explicit-date policy here,
           which is why this is a policy argument rather than a bare bool)

    Rejection paths, each a no-effect failure:
      1. bad id           -- id fails normalize_experiment_id
      2. missing experiment
      3. already in target -- the row is already at target_status; rejected
                              (consistent with link's duplicate handling --
                              "you asked for a state that is already true"
                              is surfaced, not silently no-op'd)
    """
    experiment_id = normalize_experiment_id(arguments["id"])
    target_status = arguments["target_status"]
    completion_date_policy = arguments["set_completion_date"]

    # 1. bad id
    if experiment_id is None:
        return make_failure_effects(f"Invalid experiment id: '{arguments['id']}'")

    # 2. missing experiment
    matching_rows = [
        experiment_row
        for experiment_row in store["experiments"]
        if experiment_row["id"] == experiment_id
    ]
    if not matching_rows:
        return make_failure_effects(
            f"No such experiment: {experiment_id}. "
            "Run 'lw list' to see existing experiments."
        )
    current_row = matching_rows[0]

    # 3. already in the target state
    if current_row["status"] == target_status:
        return make_failure_effects(
            f"{experiment_id} is already '{target_status}'; nothing to do."
        )

    # Resolve the completion-date policy into a concrete value.
    if completion_date_policy == "today":
        new_completion_date = clock["today"].strftime(ISO_DATE_FORMAT)
    elif completion_date_policy == "clear":
        new_completion_date = None
    else:
        # Defensive: an unknown policy is a programming error, not user
        # input -- surface it loudly rather than writing a wrong row.
        return make_failure_effects(
            f"Internal error: unknown completion-date policy "
            f"'{completion_date_policy}'."
        )

    # Build the updated store with exactly this one row changed. The
    # transform stays pure: a fresh row dict, a fresh experiments list, a
    # fresh store dict -- the input store is not mutated.
    updated_row = dict(current_row)
    updated_row["status"] = target_status
    updated_row["date_completed"] = new_completion_date

    updated_experiments = [
        updated_row if experiment_row["id"] == experiment_id else experiment_row
        for experiment_row in store["experiments"]
    ]
    updated_store = dict(store)
    updated_store["experiments"] = updated_experiments

    if new_completion_date is not None:
        confirmation_line = (
            f"Marked {experiment_id} '{target_status}' "
            f"(completed {new_completion_date})."
        )
    else:
        confirmation_line = f"Marked {experiment_id} '{target_status}'."

    return make_success_effects(
        standard_output_lines=[confirmation_line],
        store=updated_store,
    )


def build_complete_effects(store, arguments, clock):
    """(store, arguments, clock) -> effects. Pure. The narrow `complete`
    verb: a thin wrapper that drives the status-transition engine toward
    the complete status, stamping today's completion date.

    This wrapper is the entire `complete`-specific surface. A future
    `reopen` would be the symmetric wrapper (target DEFAULT_EXPERIMENT_STATUS,
    policy "clear") over the same engine.
    """
    transition_arguments = dict(arguments)
    transition_arguments["target_status"] = COMPLETE_EXPERIMENT_STATUS
    transition_arguments["set_completion_date"] = "today"
    return build_status_transition_effects(store, transition_arguments, clock)


def run_complete_command(arguments):
    """Edge for `complete`: run the pure transform, route through
    execute_effects. No filesystem edge -- completion touches only the DB.
    """
    effects = build_complete_effects(
        load_store_from_database(DATABASE_PATH), arguments, startup_clock
    )
    return execute_effects(effects, DATABASE_PATH)


# ---------------------------------------------------------------------------
# Dispatch -- plain if/elif. C15.
# ---------------------------------------------------------------------------
# The interface is unified at the CONTRACT (every command parses to a
# transform whose effects flow through execute_effects), deliberately NOT at
# the parse: each branch parses its own argument shape because the shapes are
# genuinely different KINDS (required-ordered / optional-modifier /
# optional-free-text-tail / optional-predicate). A binder would just relocate
# this if/elif into a parse_shape wearing a costume. See SS5.


def separate_flags_from_positionals(
    command_arguments, value_flag_names=(), boolean_flag_names=()
):
    """Tiny shared parser: pull named flags out of command_arguments,
    return (positional_arguments, parsed_flags). Honest helper for a real
    shared need, not a binder -- it does NOT know command shapes, it just
    separates flags from positionals. value_flag_names take the next token;
    boolean_flag_names are presence-only.
    """
    positional_arguments = []
    parsed_flags = {}
    argument_index = 0
    while argument_index < len(command_arguments):
        current_token = command_arguments[argument_index]
        if current_token.startswith("--"):
            flag_name = current_token[2:]
        else:
            flag_name = None

        if flag_name in boolean_flag_names:
            parsed_flags[flag_name] = True
            argument_index += 1
        elif flag_name in value_flag_names:
            if argument_index + 1 < len(command_arguments):
                parsed_flags[flag_name] = command_arguments[argument_index + 1]
            else:
                parsed_flags[flag_name] = None
            argument_index += 2
        else:
            positional_arguments.append(current_token)
            argument_index += 1
    return positional_arguments, parsed_flags


def append_command_log_line(command_arguments):
    """Append one audit line to command_log.txt. An EFFECT performed at the
    edge after a successful run -- never from inside a transform. Optional
    and harmless; reintroducing it as an inline side-channel would be the
    entanglement this design removes (SS6.6).
    """
    log_line = f"{startup_clock['now'].isoformat()}  {' '.join(command_arguments)}\n"
    try:
        with open(COMMAND_LOG_PATH, "a") as log_file:
            log_file.write(log_line)
    except OSError:
        # Logging must never break a command.
        pass


# --- Usage text: one declarative source of truth, rendered by helpers so
# --- every error site shows a consistent (condition -> usage -> next step)
# --- shape instead of a bare one-line stub.

PROGRAM_NAME = "lw"

# command name -> (argument signature, one-line purpose)
COMMAND_USAGE = {
    "init": (
        "",
        "Create the lab directory, database, and counter (safe to re-run).",
    ),
    "new": (
        "<type> <shortname> [--title TITLE...]",
        "Register an experiment: make its folder, row, and id.",
    ),
    "link": (
        "<from> <to> <relationship>",
        "Record a directed provenance edge between two experiments.",
    ),
    "show": (
        "<id> [--files]",
        "Print an experiment's row, links, and folder file summary.",
    ),
    "list": (
        "[--type TYPE] [--status STATUS]",
        "List experiments, optionally filtered by type and/or status.",
    ),
    "complete": (
        "<id>",
        "Mark an experiment complete and stamp its completion date.",
    ),
}


def format_command_usage(command_name):
    """Return a single 'usage: lw <command> <signature>' line."""
    argument_signature, _purpose = COMMAND_USAGE[command_name]
    if argument_signature:
        return f"usage: {PROGRAM_NAME} {command_name} {argument_signature}"
    return f"usage: {PROGRAM_NAME} {command_name}"


def format_all_commands_usage():
    """Return the multi-line overview: every command, signature, purpose."""
    overview_lines = [f"usage: {PROGRAM_NAME} <command> [arguments]", "", "commands:"]
    for command_name in COMMAND_USAGE:
        argument_signature, command_purpose = COMMAND_USAGE[command_name]
        invocation = f"{command_name} {argument_signature}".rstrip()
        overview_lines.append(f"  {invocation}")
        overview_lines.append(f"      {command_purpose}")
    return "\n".join(overview_lines)


def print_usage_error(condition_line, command_name=None, next_step_line=None):
    """Emit a consistent (condition -> usage -> next step) error to stderr.

    condition_line states what triggered the error. If command_name is
    given, the specific command's usage line follows; otherwise the full
    command overview is shown. next_step_line, if given, is appended as an
    actionable hint.
    """
    print(condition_line, file=sys.stderr)
    if command_name is not None:
        print(format_command_usage(command_name), file=sys.stderr)
    else:
        print(format_all_commands_usage(), file=sys.stderr)
    if next_step_line is not None:
        print(next_step_line, file=sys.stderr)


def main(command_arguments=None):
    if command_arguments is None:
        command_arguments = list(sys.argv[1:])
    else:
        command_arguments = list(command_arguments)

    if not command_arguments:
        print_usage_error(
            "No command given.",
            next_step_line=(
                f"Run '{PROGRAM_NAME} init' first if you haven't, then "
                f"pick a command above."
            ),
        )
        return 2

    command_name = command_arguments[0]
    remaining_arguments = command_arguments[1:]

    # init is the outsider: runs before the store exists.
    if command_name == "init":
        if remaining_arguments:
            # init takes no arguments. Extra words here are almost always a
            # `new` invocation typed with the wrong verb.
            print_usage_error(
                f"'init' takes no arguments, got {', '.join(remaining_arguments)}.",
                command_name="init",
                next_step_line=(
                    f"Did you mean '{PROGRAM_NAME} new "
                    f"{' '.join(remaining_arguments)}'?"
                ),
            )
            return 2
        exit_code = run_init_command()
        if exit_code == 0:
            append_command_log_line(command_arguments)
        return exit_code

    # Every other command requires the DB. Precondition at the store
    # boundary, before any transform.
    missing_database_effects = require_existing_database(DATABASE_PATH)
    if missing_database_effects is not None:
        return execute_effects(missing_database_effects, DATABASE_PATH)

    if command_name == "new":
        # new <type> <shortname> [--title T...]: --title joins the tail,
        # so it must come last on the line. After the tail is stripped, the
        # remaining tokens must be exactly the two positionals -- no extra
        # words, no stray flags (a stray --flag here is usually a typo'd
        # --title, e.g. --titel, that would otherwise be silently dropped).
        experiment_title = None
        if "--title" in remaining_arguments:
            title_flag_index = remaining_arguments.index("--title")
            experiment_title = " ".join(remaining_arguments[title_flag_index + 1 :])
            remaining_arguments = remaining_arguments[:title_flag_index]

        unexpected_flags = [
            token for token in remaining_arguments if token.startswith("--")
        ]
        if unexpected_flags:
            print_usage_error(
                f"'new' got an unexpected flag: {', '.join(unexpected_flags)}.",
                command_name="new",
                next_step_line=("The only flag is '--title', and it must come last."),
            )
            return 2
        if len(remaining_arguments) < 2:
            print_usage_error(
                f"'new' needs <type> and <shortname> "
                f"({len(remaining_arguments)} given).",
                command_name="new",
                next_step_line=(f"Valid types: {', '.join(VALID_EXPERIMENT_TYPES)}."),
            )
            return 2
        if len(remaining_arguments) > 2:
            print_usage_error(
                f"'new' takes <type> and <shortname>, got "
                f"{len(remaining_arguments)} positional arguments: "
                f"{', '.join(remaining_arguments)}.",
                command_name="new",
                next_step_line=(
                    "A multi-word title goes after '--title', e.g. "
                    f"'{PROGRAM_NAME} new {remaining_arguments[0]} "
                    f"{remaining_arguments[1]} --title Some Long Title'."
                ),
            )
            return 2
        exit_code = run_new_command(
            {
                "type": remaining_arguments[0],
                "shortname": remaining_arguments[1],
                "title": experiment_title,
            }
        )

    elif command_name == "link":
        # link <from> <to> <relationship>: three required, ordered args.
        if len(remaining_arguments) != 3:
            print_usage_error(
                f"'link' needs exactly <from> <to> <relationship> "
                f"({len(remaining_arguments)} given).",
                command_name="link",
                next_step_line=(
                    f"Valid relationships: {', '.join(ALLOWED_LINK_RELATIONSHIPS)}."
                ),
            )
            return 2
        exit_code = run_link_command(
            {
                "from_id": remaining_arguments[0],
                "to_id": remaining_arguments[1],
                "relationship": remaining_arguments[2],
            }
        )

    elif command_name == "show":
        # show <id> [--files]: positional + optional bool modifier.
        positional_arguments, parsed_flags = separate_flags_from_positionals(
            remaining_arguments, boolean_flag_names=("files",)
        )
        if len(positional_arguments) == 0:
            print_usage_error(
                "'show' needs an experiment <id>.",
                command_name="show",
                next_step_line=(f"Run '{PROGRAM_NAME} list' to see existing ids."),
            )
            return 2
        if len(positional_arguments) > 1:
            print_usage_error(
                f"'show' takes one <id>, got "
                f"{len(positional_arguments)}: "
                f"{', '.join(positional_arguments)}.",
                command_name="show",
                next_step_line=(
                    "Did you mean to add '--files' rather than a second id?"
                ),
            )
            return 2
        exit_code = run_show_command(
            {
                "id": positional_arguments[0],
                "show_full_listing": parsed_flags.get("files", False),
            }
        )

    elif command_name == "list":
        # list [--type T] [--status S]: optional predicates.
        positional_arguments, parsed_flags = separate_flags_from_positionals(
            remaining_arguments, value_flag_names=("type", "status")
        )
        if positional_arguments:
            # `list` takes no positionals; a bare word here is almost
            # certainly a filter value typed without its flag.
            print_usage_error(
                f"'list' takes no positional arguments, got "
                f"{', '.join(positional_arguments)}.",
                command_name="list",
                next_step_line=(
                    "Filters use flags, e.g. "
                    f"'{PROGRAM_NAME} list --type purification "
                    f"--status active'."
                ),
            )
            return 2
        exit_code = run_list_command(
            {
                "type": parsed_flags.get("type"),
                "status": parsed_flags.get("status"),
            }
        )

    elif command_name == "complete":
        # complete <id>: one positional, like show. Reject missing or extra.
        positional_arguments, _parsed_flags = separate_flags_from_positionals(
            remaining_arguments
        )
        if len(positional_arguments) == 0:
            print_usage_error(
                "'complete' needs an experiment <id>.",
                command_name="complete",
                next_step_line=(f"Run '{PROGRAM_NAME} list' to see existing ids."),
            )
            return 2
        if len(positional_arguments) > 1:
            print_usage_error(
                f"'complete' takes one <id>, got "
                f"{len(positional_arguments)}: "
                f"{', '.join(positional_arguments)}.",
                command_name="complete",
                next_step_line="Complete one experiment at a time.",
            )
            return 2
        exit_code = run_complete_command({"id": positional_arguments[0]})

    else:
        print_usage_error(
            f"Unknown command: '{command_name}'.",
            next_step_line=(
                f"Run '{PROGRAM_NAME}' with no arguments to see all commands."
            ),
        )
        return 2

    if exit_code == 0:
        append_command_log_line(command_arguments)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
