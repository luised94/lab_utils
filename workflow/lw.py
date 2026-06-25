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

# Folder-naming category codes. NOTE: `type` is NOT a data contract. The
# typed data it once dispatched to (purification_details, assay_details,
# samples) now lives in Excel. `type` is retained for exactly two uses:
# (a) it supplies the short code in the folder name (here), and (b) it
# filters `list`. Do not mistake `type` for a contract pointing at
# DB-resident data.
CATEGORY_CODES = {
    "purification": "pur",
    "loading": "load",
    "gelshift": "gs",
    "atpase": "atp",
    "genetics": "gen",
    "computational": "comp",
    "cloning": "clon",
    "sequencing": "seq",
}
VALID_TYPES = list(CATEGORY_CODES.keys())

# Closed set of edge labels. The name encodes block-on-absence semantics:
# `link` REJECTS any relationship not in this list (blocks, does not warn).
# Extending is a deliberate one-line edit. All three are now plain edge
# labels -- `uses_prep_from` no longer has special behavior.
ALLOWED_RELATIONSHIPS = ["uses_prep_from", "replicate_of", "follow_up_to"]

# Captured once at startup, passed as the "clock" to transforms so the
# transforms never read the clock themselves.
clock = {"today": datetime.date.today(), "now": datetime.datetime.now()}

# --- LAB_ROOT resolution ---
_lab_root_env = os.environ.get("LW_LAB_ROOT") or os.environ.get("LAB_ROOT")
if _lab_root_env:
    LAB_ROOT = _lab_root_env
else:
    _win_user = os.environ.get("LW_WINDOWS_USER") or os.environ.get("MC_WINDOWS_USER")
    if _win_user:
        LAB_ROOT = f"/mnt/c/Users/{_win_user}/Desktop/lab"
    else:
        print("Error: Set LW_LAB_ROOT or LW_WINDOWS_USER environment variable.")
        sys.exit(1)

# --- derived paths ---
DB_PATH = os.path.join(LAB_ROOT, "lab.db")
EXPERIMENTS_DIR = os.path.join(LAB_ROOT, "experiments")
COUNTER_FILE = os.path.join(LAB_ROOT, "counter.txt")


# ---------------------------------------------------------------------------
# C2 -- effects-as-data primitives (the core contract)
# ---------------------------------------------------------------------------


def effects_ok(stdout=None, store=None, exit_code=0):
    return {
        "ok": True,
        "stdout": stdout or [],
        "stderr": [],
        "store": store,
        "exit_code": exit_code,
    }


def effects_fail(msg, exit_code=1):
    if isinstance(msg, str):
        msg = [msg]
    return {
        "ok": False,
        "stdout": [],
        "stderr": msg,
        "store": None,
        "exit_code": exit_code,
    }


def execute_effects(effects, db_path):
    """The ONLY function that performs IO for transformed commands."""
    for line in effects.get("stdout", []):
        print(line)
    for line in effects.get("stderr", []):
        print(line, file=sys.stderr)
    if effects.get("store") is not None:
        commit(effects["store"], db_path)
    return effects.get("exit_code", 0)


# ---------------------------------------------------------------------------
# C3 -- id normalization
# ---------------------------------------------------------------------------


def normalize_exp_id(raw_id):
    """'LM-0003' or '3' -> 'LM-0003'; None if invalid."""
    if raw_id.upper().startswith("LM-"):
        return raw_id.upper()
    try:
        return f"LM-{int(raw_id):04d}"
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

SCHEMA = """
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
LIVE_TABLES = ["experiments", "experiment_links"]


def load_store(db_path):
    """Read every live table into a dict-of-lists of row dicts.

    Reads ONLY the two live tables. The orphaned tables sit inert.
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    store = {}
    try:
        for table in LIVE_TABLES:
            cur = conn.execute(f"SELECT * FROM {table}")
            store[table] = [dict(row) for row in cur.fetchall()]
    finally:
        conn.close()
    return store


# C5 -- commit (whole-store rewrite + union-of-keys fix)
# ---------------------------------------------------------------------------


def commit(store, db_path):
    """Rewrite the whole store in one transaction.

    This whole-store-rewrite is DELIBERATELY RETAINED. It works for a
    single-user, few-hundred-row DB. Do NOT replace with per-statement
    effects unless a real problem is reported (row count, a second writer).

    The column list per table is built from the UNION of all rows' keys,
    and values are read with row.get(c) rather than row[c] -- this is the
    fix for the KeyError the naive version raised on a missing dict key.
    Touches only the live tables that load_store loaded; orphaned tables are
    never read or written here.
    """
    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA)
        cur = conn.cursor()
        for table in LIVE_TABLES:
            rows = store.get(table, [])
            cur.execute(f"DELETE FROM {table}")
            if not rows:
                continue
            columns = []
            seen = set()
            for r in rows:
                for k in r:
                    if k not in seen:
                        seen.add(k)
                        columns.append(k)
            placeholders = ", ".join(["?"] * len(columns))
            sql = f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({placeholders})"
            for row in rows:
                cur.execute(sql, [row.get(c) for c in columns])
        conn.commit()
    finally:
        conn.close()


# C6 -- missing-DB precondition
# ---------------------------------------------------------------------------


def require_db(db_path):
    """Store-boundary precondition for every command except init.

    Returns None if the DB exists; otherwise an effects_fail describing the
    clean error. This runs BEFORE any transform -- transforms assume a
    loaded store. Distinct from a missing experiment FOLDER on show, which
    is a normal rendered state, not a precondition error.
    """
    if not os.path.exists(db_path):
        return effects_fail("No lab database found. Run 'lw init' first.")
    return None


# C7 -- init (the bootstrap outsider)
# ---------------------------------------------------------------------------


def cmd_init():
    """Create LAB_ROOT/experiments, the DB, and seed counter.txt if absent.

    The outsider: runs before the store exists, so it cannot use the store
    machinery. Idempotent -- CREATE TABLE IF NOT EXISTS and the counter
    seed both no-op on re-run. Returns an exit code directly (it is its own
    edge, not a transform routed through execute_effects).
    """
    os.makedirs(LAB_ROOT, exist_ok=True)
    os.makedirs(EXPERIMENTS_DIR, exist_ok=True)

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.executescript(SCHEMA)
        conn.commit()
    finally:
        conn.close()

    if not os.path.exists(COUNTER_FILE):
        with open(COUNTER_FILE, "w") as f:
            f.write("1")

    print(f"Initialized lab at {LAB_ROOT}")
    return 0


# ---------------------------------------------------------------------------
# new -- the #1 feature, decomposed: counter edge (C8) + pure transform (C9)
# + edge/wiring (C10). The one place transform purity, the counter-as-effect
# boundary, fixed bump-last ordering, and filesystem IO all intersect.
# ---------------------------------------------------------------------------

# C8 -- counter edge (read + bump), ordering-sensitive
# ---------------------------------------------------------------------------


def read_counter():
    """Read the current 'next experiment number' from counter.txt.

    Filesystem edge -- never called from inside a transform. The value is
    passed INTO the new transform via args; the transform stays pure.
    """
    with open(COUNTER_FILE) as f:
        return int(f.read().strip())


def write_counter(value):
    """Write the bumped counter value to counter.txt.

    Filesystem edge. Called LAST in the new sequence (after row insert and
    mkdir), so a run that dies mid-way wastes a number (harmless) rather
    than reusing one (a PK collision -- a real failure). Bumping last fails
    toward the harmless direction. Honor the ordering.
    """
    with open(COUNTER_FILE, "w") as f:
        f.write(str(value))


# C9 -- new pure transform
# ---------------------------------------------------------------------------


def transform_new(store, args, clock):
    """(store, args, clock) -> effects. Pure: no IO, no clock-read, no print.

    args carries: type, shortname, title (optional), and the counter value
    read at the edge (C8). Builds exp_id / folder_name purely, produces the
    experiment row, and returns the new counter value as an extra
    'counter' key on the effects dict for the new edge (C10) to write. The
    fixed effects_ok shape is unchanged; the extra key rides alongside it.
    """
    exp_type = args["type"]
    shortname = args["shortname"]
    title = args.get("title") or shortname
    counter = args["counter"]

    if exp_type not in VALID_TYPES:
        return effects_fail(
            f"Unknown type '{exp_type}'. Valid types: {', '.join(VALID_TYPES)}"
        )

    exp_id = f"LM-{counter:04d}"
    today_compact = clock["today"].strftime("%Y%m%d")
    cat_code = CATEGORY_CODES[exp_type]
    folder_name = f"{today_compact}_{exp_id}_{cat_code}_{shortname}"
    today_iso = clock["today"].strftime("%Y-%m-%d")

    row = {
        "id": exp_id,
        "type": exp_type,
        "title": title,
        "shortname": shortname,
        "date_started": today_iso,
        "date_completed": None,
        "status": "active",
        "folder_name": folder_name,
        "notes": None,
    }

    new_store = dict(store)
    new_store["experiments"] = list(store["experiments"]) + [row]

    effects = effects_ok(
        stdout=[f"Created {exp_id}: {folder_name}"],
        store=new_store,
    )
    # Extra key on top of the fixed effects_ok shape: the new counter value
    # the C10 edge writes LAST. Not part of execute_effects' duties.
    effects["counter"] = counter + 1
    effects["folder_name"] = folder_name
    return effects


# C10 -- new edge + wiring
# ---------------------------------------------------------------------------


def cmd_new(args):
    """Edge for `new`: read counter, run the pure transform, sequence the
    effects in the fixed fail-safe order, route through execute_effects.

    Order is deliberate (SS6.6): row insert (via execute_effects -> commit)
    and mkdir happen FIRST; the counter bump happens LAST. A run that dies
    before the bump wastes a number (harmless); never reuses one.
    """
    counter = read_counter()
    args = dict(args)
    args["counter"] = counter

    effects = transform_new(load_store(DB_PATH), args, clock)

    if not effects["ok"]:
        # Rejected: no row, no mkdir, no bump. Just emit the failure.
        return execute_effects(effects, DB_PATH)

    # FIRST: mkdir the experiment folder.
    folder_path = os.path.join(EXPERIMENTS_DIR, effects["folder_name"])
    os.makedirs(folder_path, exist_ok=True)

    # FIRST (cont.): commit the row + print stdout via execute_effects.
    exit_code = execute_effects(effects, DB_PATH)

    # LAST: bump the counter, only after row + folder are durable.
    write_counter(effects["counter"])

    return exit_code


# ---------------------------------------------------------------------------
# link -- append a bare provenance edge. C11 transform + C12 wiring.
# ---------------------------------------------------------------------------

# C11 -- link transform (bare edge; the ~20-line version)
# ---------------------------------------------------------------------------


def transform_link(store, args, clock):
    """(store, args, clock) -> effects. Pure. Appends one bare edge.

    Four rejection paths, each a no-effect failure:
      1. bad id            -- either id fails normalize_exp_id
      2. missing experiment-- either id not present in the store
      3. unknown relationship -- not in ALLOWED_RELATIONSHIPS (block, not warn)
      4. duplicate edge    -- (from, to) already present
    No assay_details, no type-specific behavior -- a plain edge label.
    """
    from_id = normalize_exp_id(args["from_id"])
    to_id = normalize_exp_id(args["to_id"])
    relationship = args["relationship"]

    # 1. bad id
    if from_id is None or to_id is None:
        bad = args["from_id"] if from_id is None else args["to_id"]
        return effects_fail(f"Invalid experiment id: '{bad}'")

    # 2. missing experiment
    known = {r["id"] for r in store["experiments"]}
    if from_id not in known:
        return effects_fail(f"No such experiment: {from_id}")
    if to_id not in known:
        return effects_fail(f"No such experiment: {to_id}")

    # 3. unknown relationship -- blocked, not warned
    if relationship not in ALLOWED_RELATIONSHIPS:
        return effects_fail(
            f"Unknown relationship '{relationship}'. "
            f"Allowed: {', '.join(ALLOWED_RELATIONSHIPS)}"
        )

    # 4. duplicate edge
    for link in store["experiment_links"]:
        if link["from_experiment"] == from_id and link["to_experiment"] == to_id:
            return effects_fail(f"Link already exists: {from_id} -> {to_id}")

    edge = {
        "from_experiment": from_id,
        "to_experiment": to_id,
        "relationship": relationship,
    }
    new_store = dict(store)
    new_store["experiment_links"] = list(store["experiment_links"]) + [edge]

    return effects_ok(
        stdout=[f"Linked {from_id} -> {to_id} ({relationship})"],
        store=new_store,
    )


# C12 -- link wiring
# ---------------------------------------------------------------------------


def cmd_link(args):
    """Edge for `link`: run the pure transform, route through execute_effects."""
    effects = transform_link(load_store(DB_PATH), args, clock)
    return execute_effects(effects, DB_PATH)


# ---------------------------------------------------------------------------
# show -- experiment row + links + folder file summary. C13.
# ---------------------------------------------------------------------------

# C13 -- show transform + listdir edge
# ---------------------------------------------------------------------------

# Sentinel the edge passes to the transform when the experiment folder is
# absent. Missing folder is a NORMAL state (the user hand-edits the
# filesystem), rendered not raised -- distinct from the missing-DB hard
# precondition (SS2).
FOLDER_NOT_FOUND = object()


def transform_show(store, args, clock):
    """(store, args, clock) -> effects. Pure: receives the file listing (or
    the FOLDER_NOT_FOUND sentinel) from the edge; never touches the FS.

    args carries: id, files (list[str] | FOLDER_NOT_FOUND), folder_path,
    full (bool -- the --files flag).
    """
    exp_id = normalize_exp_id(args["id"])
    if exp_id is None:
        return effects_fail(f"Invalid experiment id: '{args['id']}'")

    rows = [r for r in store["experiments"] if r["id"] == exp_id]
    if not rows:
        return effects_fail(f"No such experiment: {exp_id}")
    row = rows[0]

    out = []
    out.append(f"{row['id']}: {row['title']}")
    out.append(f"  type:      {row['type']}")
    out.append(f"  shortname: {row['shortname']}")
    out.append(f"  status:    {row['status']}")
    out.append(f"  started:   {row['date_started']}")
    if row.get("date_completed"):
        out.append(f"  completed: {row['date_completed']}")
    out.append(f"  folder:    {row['folder_name']}")

    # Links: outgoing and incoming.
    outgoing = [l for l in store["experiment_links"] if l["from_experiment"] == exp_id]
    incoming = [l for l in store["experiment_links"] if l["to_experiment"] == exp_id]
    if outgoing:
        out.append("  links out:")
        for l in outgoing:
            out.append(f"    -> {l['to_experiment']} ({l['relationship']})")
    if incoming:
        out.append("  links in:")
        for l in incoming:
            out.append(f"    <- {l['from_experiment']} ({l['relationship']})")

    # Files: folder IS the record. Date-front-loaded names => lexical sort
    # equals chronological, so no stat calls.
    files = args["files"]
    if files is FOLDER_NOT_FOUND:
        out.append(f"  Files: folder not found ({args['folder_path']})")
    elif args.get("full"):
        out.append(f"  Files: {len(files)}")
        for name in sorted(files):
            out.append(f"    {name}")
    else:
        if files:
            newest = sorted(files)[-1]
            out.append(f"  Files: {len(files)} (last: {newest})")
        else:
            out.append("  Files: 0")

    return effects_ok(stdout=out)


def cmd_show(args):
    """Edge for `show`: unconditional listdir at the FS edge, missing folder
    caught into a sentinel, then the pure transform formats the output.
    """
    store = load_store(DB_PATH)
    exp_id = normalize_exp_id(args["id"])

    folder_path = ""
    files = FOLDER_NOT_FOUND
    if exp_id is not None:
        rows = [r for r in store["experiments"] if r["id"] == exp_id]
        if rows:
            folder_path = os.path.join(EXPERIMENTS_DIR, rows[0]["folder_name"])
            try:
                files = os.listdir(folder_path)
            except OSError:
                files = FOLDER_NOT_FOUND

    args = dict(args)
    args["files"] = files
    args["folder_path"] = folder_path
    effects = transform_show(store, args, clock)
    return execute_effects(effects, DB_PATH)


# ---------------------------------------------------------------------------
# list -- filter experiments by --type / --status. C14.
# ---------------------------------------------------------------------------

# C14 -- list transform + flag parse
# ---------------------------------------------------------------------------


def transform_list(store, args, clock):
    """(store, args, clock) -> effects. Pure. Absent flag = no filter on
    that field (an absent predicate is meaningful: it matches everything).
    """
    type_filter = args.get("type")
    status_filter = args.get("status")

    rows = store["experiments"]
    if type_filter is not None:
        rows = [r for r in rows if r["type"] == type_filter]
    if status_filter is not None:
        rows = [r for r in rows if r["status"] == status_filter]

    if not rows:
        return effects_ok(stdout=["No experiments match."])

    out = []
    for r in sorted(rows, key=lambda r: r["id"]):
        out.append(f"{r['id']}  {r['status']:8}  {r['type']:13}  {r['title']}")
    return effects_ok(stdout=out)


def cmd_list(args):
    """Edge for `list`: run the pure transform, route through execute_effects."""
    effects = transform_list(load_store(DB_PATH), args, clock)
    return execute_effects(effects, DB_PATH)
