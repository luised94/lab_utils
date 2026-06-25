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
