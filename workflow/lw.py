#!/usr/bin/env python3
# ============================================================
# SECTION 0: IMPORTS
# ============================================================
import sqlite3
import sys
import os
import copy
import datetime
import shutil
import itertools
import csv
import json
import re

# ============================================================
# SECTION 1: CONFIGURATION (module-level constants and sentinels)
# ============================================================

# --- path sentinels (assigned by main() at runtime) ---
LAB_ROOT = None
DB_PATH = None
EXPERIMENTS_DIR = None
STAGING_DIR = None
COUNTER_FILE = None
LAB_NOTES = None

DIRECTORIES = [
    "experiments",
    "protocols",
    "publications",
    "resources",
    "staging",
    "excel_templates",
]

# --- experiment types and category codes ---
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
# VALID_TYPES derived from CATEGORY_CODES - add new types to CATEGORY_CODES only
VALID_TYPES = list(CATEGORY_CODES.keys())

PROTECTED_FIELDS = {
    "id",
    "experiment_id",
    "created_at",
    "folder_name",
    "date_started",
    "date_completed",
    "status",
}

KNOWN_RELATIONSHIPS = ["uses_prep_from", "replicate_of", "follow_up_to"]

# --- schema ---
SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS experiments (
    id TEXT PRIMARY KEY,
    type TEXT NOT NULL,
    title TEXT NOT NULL,
    shortname TEXT NOT NULL,
    date_started TEXT NOT NULL,
    date_completed TEXT,
    status TEXT DEFAULT 'active',
    folder_name TEXT NOT NULL,
    protocol_id TEXT,
    notes TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);
CREATE TABLE IF NOT EXISTS strains (
    id TEXT PRIMARY KEY,
    genotype TEXT NOT NULL,
    label TEXT,
    parent TEXT,
    construction TEXT,
    selection_marker TEXT,
    verification TEXT,
    storage_location TEXT,
    notes TEXT
);
CREATE TABLE IF NOT EXISTS experiment_strains (
    experiment_id TEXT REFERENCES experiments(id),
    strain_id TEXT REFERENCES strains(id),
    role TEXT,
    PRIMARY KEY (experiment_id, strain_id)
);
CREATE TABLE IF NOT EXISTS experiment_links (
    from_experiment TEXT REFERENCES experiments(id),
    to_experiment TEXT REFERENCES experiments(id),
    relationship TEXT NOT NULL,
    PRIMARY KEY (from_experiment, to_experiment)
);
CREATE TABLE IF NOT EXISTS files (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id TEXT REFERENCES experiments(id),
    filename TEXT NOT NULL,
    file_type TEXT,
    instrument TEXT,
    description TEXT,
    source_file_id INTEGER REFERENCES files(id),
    created_at TEXT DEFAULT (datetime('now'))
);
CREATE TABLE IF NOT EXISTS purification_details (
    experiment_id TEXT PRIMARY KEY REFERENCES experiments(id),
    protein TEXT NOT NULL,
    induction_od REAL,
    galactose_hours REAL,
    columns TEXT,
    v5_depletion INTEGER DEFAULT 0,
    fractions_pooled TEXT,
    yield_mg REAL,
    storage TEXT
);
CREATE TABLE IF NOT EXISTS assay_details (
    experiment_id TEXT PRIMARY KEY REFERENCES experiments(id),
    protein_prep TEXT REFERENCES experiments(id),
    dna_template TEXT,
    dna_prep_method TEXT,
    conditions TEXT
);
CREATE TABLE IF NOT EXISTS figure_panels (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    publication TEXT NOT NULL,
    figure TEXT NOT NULL,
    panel TEXT,
    description TEXT,
    experiment_id TEXT REFERENCES experiments(id),
    source_file TEXT,
    ai_file TEXT
);
CREATE TABLE IF NOT EXISTS samples (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id TEXT NOT NULL REFERENCES experiments(id),
    position INTEGER NOT NULL,
    sample_type TEXT NOT NULL DEFAULT 'factorial',
    strain_id TEXT REFERENCES strains(id),
    dna TEXT,
    label TEXT,
    conditions TEXT,
    created_at TEXT DEFAULT (datetime('now')),
    UNIQUE(experiment_id, position)
);
CREATE TABLE IF NOT EXISTS stage_expectations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id TEXT NOT NULL REFERENCES experiments(id),
    slot INTEGER NOT NULL,
    descriptor TEXT NOT NULL,
    status TEXT DEFAULT 'pending',
    created_at TEXT DEFAULT (datetime('now')),
    filed_at TEXT,
    UNIQUE(experiment_id, slot)
);
"""


# ============================================================
# SECTION 1b: STORE (data-oriented state container)
# ============================================================
# The Store is a plain dict of lists-of-dicts.  Every table in the
# schema maps to a key whose value is a list of row-dicts (column
# names as keys).  load_store() reads the DB into this shape;
# commit() writes it back atomically.  Between those two calls,
# all transforms operate on the Store value -- no cursor needed.

# Tables in write-back order (children before parents for FK safety
# on delete; parents before children on insert).
# Tables with INTEGER PRIMARY KEY AUTOINCREMENT preserve their IDs
# explicitly so foreign-key references remain stable.
STORE_TABLES = [
    "experiments",
    "strains",
    "experiment_strains",
    "experiment_links",
    "files",
    "purification_details",
    "assay_details",
    "figure_panels",
    "samples",
    "stage_expectations",
]

# Tables whose primary key is INTEGER AUTOINCREMENT -- IDs must be
# preserved on write-back to avoid breaking FK references.
_AUTOINCREMENT_TABLES = {"files", "figure_panels", "samples", "stage_expectations"}


def load_store(db_path):
    """Read every table into a dict-of-lists-of-dicts.

    Returns {"experiments": [{...}, ...], "strains": [...], ...}.
    Each row is a plain dict with column names as keys.
    Opens and closes its own connection (read-only snapshot).
    """
    _conn = sqlite3.connect(db_path)
    _conn.execute("PRAGMA foreign_keys = ON")
    store = {}
    for table in STORE_TABLES:
        _cur = _conn.execute(f"SELECT * FROM {table}")
        col_names = [d[0] for d in _cur.description]
        store[table] = [dict(zip(col_names, row)) for row in _cur.fetchall()]
    _conn.close()
    return store


def commit(store, db_path):
    """Write the entire Store back to the database atomically.

    Strategy: open a connection, begin a transaction, delete all rows
    from each table (children first), then insert all rows (parents
    first).  For AUTOINCREMENT tables the stored id is written
    explicitly so FK references survive.

    Write-order policy:
      - The DB commit is atomic (sqlite transaction).
      - Filesystem operations (folder creation, file moves) that
        already happened are NOT rolled back by a DB failure.
      - This matches the original behavior but is now documented in
        one place.
    """
    _conn = sqlite3.connect(db_path)
    _conn.execute("PRAGMA foreign_keys = OFF")  # defer FK checks during rewrite
    _cur = _conn.cursor()
    try:
        _cur.execute("BEGIN")
        # delete in reverse order (children first)
        for table in reversed(STORE_TABLES):
            _cur.execute(f"DELETE FROM {table}")
        # insert in forward order (parents first)
        for table in STORE_TABLES:
            rows = store.get(table, [])
            if not rows:
                continue
            # Union of keys across all rows, so rows that omit a
            # nullable column (e.g. protocol_id) do not raise KeyError.
            # Preserve first-seen order for stable, readable SQL.
            columns = []
            seen = set()
            for r in rows:
                for k in r:
                    if k not in seen:
                        seen.add(k)
                        columns.append(k)
            placeholders = ", ".join(["?"] * len(columns))
            col_str = ", ".join(columns)
            sql = f"INSERT INTO {table} ({col_str}) VALUES ({placeholders})"
            for row in rows:
                _cur.execute(sql, [row.get(c) for c in columns])
        _cur.execute("COMMIT")
    except Exception:
        _cur.execute("ROLLBACK")
        raise
    finally:
        _conn.execute("PRAGMA foreign_keys = ON")
        _conn.close()


# ============================================================
# SECTION 2: HELPER FUNCTIONS (module-level, reference global sentinels)
# ============================================================

# ============================================================
# SECTION 2a: EFFECTS (effects-as-data infrastructure)
# ============================================================
# An Effects value is a plain dict describing what should happen,
# not a function that does it.  Transforms return Effects;
# execute_effects is the single "dumb shell" that performs IO.
#
# Shape:
#   {
#       "ok":        bool,           # True = success, False = error
#       "stdout":    [str, ...],     # lines to print (no trailing newline)
#       "stderr":    [str, ...],     # error lines (printed to stderr)
#       "store":     dict | None,    # mutated Store to commit, or None
#       "exit_code": int,            # 0 = success, nonzero = failure
#   }


def effects_ok(stdout=None, store=None, exit_code=0):
    """Construct a successful Effects value.

    stdout: list of output lines (each without trailing newline).
    store:  the mutated Store dict to commit, or None for read-only commands.
    """
    return {
        "ok": True,
        "stdout": stdout if stdout is not None else [],
        "stderr": [],
        "store": store,
        "exit_code": exit_code,
    }


def effects_fail(msg, exit_code=1):
    """Construct a failure Effects value.

    msg: a single error message string (or list of strings).
    """
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
    """The dumb shell: perform the IO described by an Effects value.

    1. Print stdout lines to sys.stdout.
    2. Print stderr lines to sys.stderr.
    3. If effects["store"] is not None, call commit(store, db_path).
    4. Return effects["exit_code"].

    This is the ONLY function that performs IO for transformed commands.
    """
    for line in effects.get("stdout", []):
        print(line)
    for line in effects.get("stderr", []):
        print(line, file=sys.stderr)
    if effects.get("store") is not None:
        commit(effects["store"], db_path)
    return effects.get("exit_code", 0)


# ============================================================
# SECTION 2b: PURE HELPERS AND TRANSFORMS
# ============================================================
# Each transform_<command> function is pure: it takes (store, args, clock)
# and returns an Effects dict.  No IO, no globals, no cursor.
# Transforms are registered in TRANSFORM_DISPATCH and called by main().


# --- Commit 3.1: shared append_note helper ---
def append_note(old_value, new_text, clock):
    """Pure: append a timestamped note entry.

    Returns the new notes string with [YYYY-MM-DD] prefix.
    old_value may be None (first note) or a string (append).
    """
    today = clock["today"].strftime("%Y-%m-%d")
    entry = f"[{today}] {new_text}"
    if old_value:
        return old_value + "\n" + entry
    return entry


# --- Commit 3.2: field-to-table mapping ---
# Replaces PRAGMA table_info discovery.  Maps field names to the table
# that owns them.  Protected fields are excluded (they cannot be updated).
FIELD_TABLES = {}
_FIELD_TABLE_DEFS = {
    "experiments": [
        "title",
        "shortname",
        "type",
        "protocol_id",
        "notes",
    ],
    "purification_details": [
        "protein",
        "induction_od",
        "galactose_hours",
        "columns",
        "v5_depletion",
        "fractions_pooled",
        "yield_mg",
        "storage",
    ],
    "assay_details": [
        "protein_prep",
        "dna_template",
        "dna_prep_method",
        "conditions",
    ],
}
for _tbl, _fields in _FIELD_TABLE_DEFS.items():
    for _f in _fields:
        FIELD_TABLES[_f] = _tbl

# All valid strain fields (for strain update)
STRAIN_FIELDS = [
    "genotype",
    "label",
    "parent",
    "construction",
    "selection_marker",
    "verification",
    "storage_location",
    "notes",
]


# --- Pure Store helpers ---


def lookup_strain(store, strain_id):
    """Pure lookup: find strain by ID in the Store. Returns dict or None."""
    for row in store["strains"]:
        if row["id"] == strain_id:
            return row
    return None


def store_update_row(store, table, match_key, match_val, updates):
    """Pure: return a new store with one row in `table` updated.

    Finds the first row where row[match_key] == match_val and applies
    `updates` (a dict of field->value).  Returns a new store (shallow
    copy of the table list with the updated row replaced).
    Returns (new_store, old_row) or (store, None) if not found.
    """
    new_rows = []
    old_row = None
    for row in store[table]:
        if old_row is None and row[match_key] == match_val:
            old_row = dict(row)
            updated = dict(row)
            updated.update(updates)
            new_rows.append(updated)
        else:
            new_rows.append(row)
    if old_row is None:
        return store, None
    new_store = dict(store)
    new_store[table] = new_rows
    return new_store, old_row


def store_append_row(store, table, row):
    """Pure: return a new store with `row` appended to `table`."""
    new_store = dict(store)
    new_store[table] = list(store[table]) + [row]
    return new_store


def transform_exp_complete(store, args, clock):
    """Pure transform: mark an experiment as complete.

    args: {"exp_id": str (raw, unnormalized)}
    Returns Effects with mutated store or error/no-op message.
    """
    raw_id = args["exp_id"]
    exp_id = normalize_exp_id(raw_id)
    if not exp_id:
        return effects_fail(f"Invalid experiment ID: '{raw_id}'")
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    if exp["status"] == "complete":
        return effects_ok(stdout=[f"{exp_id} is already marked complete."])
    today = clock["today"].strftime("%Y-%m-%d")
    started = datetime.date.fromisoformat(exp["date_started"])
    duration = (clock["today"] - started).days
    new_store, _ = store_update_row(
        store,
        "experiments",
        "id",
        exp_id,
        {"status": "complete", "date_completed": today},
    )
    lines = [
        f"Completed: {exp_id}",
        f"  Date: {today}",
        f"  Duration: {duration}d (started {exp['date_started']})",
    ]
    return effects_ok(stdout=lines, store=new_store)


def transform_exp_show(store, args, clock):
    """Pure transform: display experiment details.

    args: {"exp_id": str (raw, unnormalized)}
    Returns Effects with stdout lines (read-only, no store mutation).
    """
    raw_id = args["exp_id"]
    exp_id = normalize_exp_id(raw_id)
    if not exp_id:
        return effects_fail(f"Invalid experiment ID: '{raw_id}'")
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    lines = []
    # compute status display with age/duration
    if exp["status"] == "active":
        started = datetime.date.fromisoformat(exp["date_started"])
        age = (clock["today"] - started).days
        status_str = f"active ({age}d)"
    elif exp["status"] == "complete" and exp["date_completed"]:
        started = datetime.date.fromisoformat(exp["date_started"])
        completed = datetime.date.fromisoformat(exp["date_completed"])
        duration = (completed - started).days
        status_str = f"complete ({duration}d)"
    else:
        status_str = exp["status"]
    # main fields
    fields = [
        ("ID", exp["id"]),
        ("Type", exp["type"]),
        ("Title", exp["title"]),
        ("Shortname", exp["shortname"]),
        ("Date started", exp["date_started"]),
        ("Date completed", exp["date_completed"]),
        ("Status", status_str),
        ("Folder", exp["folder_name"] + "/"),
        ("Protocol", exp["protocol_id"]),
        ("Notes", exp["notes"]),
    ]
    label_width = max(len(f[0]) for f in fields)
    for label, val in fields:
        lines.append(f"  {label + ':':<{label_width + 2}} {val if val else '-'}")
    # linked strains
    strains = [
        (r["strain_id"], r["role"])
        for r in store["experiment_strains"]
        if r["experiment_id"] == exp_id
    ]
    if strains:
        lines.append("")
        lines.append("  Strains:")
        for sid, role in strains:
            role_str = f" ({role})" if role else ""
            lines.append(f"    {sid}{role_str}")
    # linked experiments
    links_out = [
        (r["to_experiment"], r["relationship"])
        for r in store["experiment_links"]
        if r["from_experiment"] == exp_id
    ]
    links_in = [
        (r["from_experiment"], r["relationship"])
        for r in store["experiment_links"]
        if r["to_experiment"] == exp_id
    ]
    if links_out or links_in:
        lines.append("")
        lines.append("  Links:")
        for target, rel in links_out:
            lines.append(f"    -> {target} ({rel})")
        for source, rel in links_in:
            lines.append(f"    <- {source} ({rel})")
    # files
    exp_files = [
        (r["filename"], r["file_type"], r["description"])
        for r in store["files"]
        if r["experiment_id"] == exp_id
    ]
    if exp_files:
        lines.append("")
        lines.append("  Files:")
        for fname, ftype, desc in exp_files:
            parts = [fname]
            if ftype:
                parts.append(f"[{ftype}]")
            if desc:
                parts.append(desc)
            lines.append(f"    {' '.join(parts)}")
    # purification details
    pur_rows = [
        r for r in store["purification_details"] if r["experiment_id"] == exp_id
    ]
    if pur_rows:
        pur = pur_rows[0]
        lines.append("")
        lines.append("  Purification details:")
        for k, v in pur.items():
            if k == "experiment_id" or v is None:
                continue
            lines.append(f"    {k}: {v}")
    # assay details
    assay_rows = [r for r in store["assay_details"] if r["experiment_id"] == exp_id]
    if assay_rows:
        assay = assay_rows[0]
        lines.append("")
        lines.append("  Assay details:")
        for k, v in assay.items():
            if k == "experiment_id" or v is None:
                continue
            lines.append(f"    {k}: {v}")
    # samples
    sample_rows = sorted(
        [r for r in store["samples"] if r["experiment_id"] == exp_id],
        key=lambda r: r["position"],
    )
    if sample_rows:
        lines.append("")
        lines.append(f"  Samples ({len(sample_rows)}):")
        for s in sample_rows:
            pos = s["position"]
            stype = s["sample_type"]
            sid = s["strain_id"]
            lbl = s["label"]
            dna = s["dna"]
            sid_str = sid if sid else "-"
            dna_str = dna if dna else ""
            extra = f"  dna={dna_str}" if dna_str else ""
            lines.append(f"    {pos:>3}  {stype:<10} {sid_str:<8} {lbl}{extra}")
    return effects_ok(stdout=lines)


def transform_exp_update(store, args, clock):
    """Pure transform: update a field on an experiment.

    args: {"exp_id": str, "field": str, "value": str}
    Handles notes-append, skeleton-row creation for purification/assay,
    protected field rejection.
    """
    raw_id = args["exp_id"]
    field = args["field"]
    value = args["value"]
    exp_id = normalize_exp_id(raw_id)
    if not exp_id:
        return effects_fail(f"Invalid experiment ID: '{raw_id}'")
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    # reject protected fields
    if field in PROTECTED_FIELDS:
        msgs = [f"Field '{field}' cannot be modified directly."]
        if field == "status":
            msgs.append(f"Use: python lw.py exp complete {exp_id}")
        return effects_fail(msgs)
    # find which table owns this field
    target_table = FIELD_TABLES.get(field)
    if not target_table:
        all_fields = sorted(FIELD_TABLES.keys())
        return effects_fail(
            [
                f"Unknown field: '{field}'",
                f"Valid fields: {', '.join(all_fields)}",
            ]
        )
    lines = []
    new_store = dict(store)
    # for type-specific tables, ensure row exists
    if target_table in ("purification_details", "assay_details"):
        existing = [r for r in store[target_table] if r["experiment_id"] == exp_id]
        if not existing:
            # insert skeleton row
            if target_table == "purification_details":
                if field == "protein":
                    new_row = {
                        "experiment_id": exp_id,
                        "protein": value,
                        "induction_od": None,
                        "galactose_hours": None,
                        "columns": None,
                        "v5_depletion": 0,
                        "fractions_pooled": None,
                        "yield_mg": None,
                        "storage": None,
                    }
                    new_store = store_append_row(
                        new_store, "purification_details", new_row
                    )
                    lines.append(f"Created purification record for {exp_id}")
                    lines.append(f"  protein: - -> {value}")
                    return effects_ok(stdout=lines, store=new_store)
                else:
                    return effects_fail(
                        [
                            f"No purification record for {exp_id}.",
                            f"Set 'protein' first: python lw.py exp update {exp_id} protein <name>",
                        ]
                    )
            elif target_table == "assay_details":
                new_row = {
                    "experiment_id": exp_id,
                    "protein_prep": None,
                    "dna_template": None,
                    "dna_prep_method": None,
                    "conditions": None,
                }
                new_store = store_append_row(new_store, "assay_details", new_row)
                lines.append(f"Created assay record for {exp_id}")
    # get old value
    id_col = "id" if target_table == "experiments" else "experiment_id"
    target_rows = [r for r in new_store[target_table] if r[id_col] == exp_id]
    old_val = target_rows[0][field] if target_rows else None
    original_value = value
    # notes: append with timestamp
    if field == "notes":
        value = append_note(old_val, value, clock)
    # apply update
    new_store, _ = store_update_row(
        new_store, target_table, id_col, exp_id, {field: value}
    )
    lines.append(f"Updated {exp_id}")
    if field == "notes":
        today = clock["today"].strftime("%Y-%m-%d")
        lines.append(f"  {field}: appended [{today}] {original_value}")
    else:
        lines.append(f"  {field}: {old_val if old_val else '-'} -> {value}")
    return effects_ok(stdout=lines, store=new_store)


def transform_exp_link(store, args, clock):
    """Pure transform: link two experiments.

    args: {"from_id": str, "to_id": str, "relationship": str}
    Handles uses_prep_from auto-populate of assay_details.
    """
    raw_from = args["from_id"]
    raw_to = args["to_id"]
    relationship = args["relationship"]
    from_id = normalize_exp_id(raw_from)
    if not from_id:
        return effects_fail(f"Invalid experiment ID: '{raw_from}'")
    to_id = normalize_exp_id(raw_to)
    if not to_id:
        return effects_fail(f"Invalid experiment ID: '{raw_to}'")
    from_exp = lookup_experiment(store, from_id)
    if not from_exp:
        return effects_fail(f"No experiment found with ID {from_id}")
    to_exp = lookup_experiment(store, to_id)
    if not to_exp:
        return effects_fail(f"No experiment found with ID {to_id}")
    lines = []
    # warn on unknown relationship (don't block)
    if relationship not in KNOWN_RELATIONSHIPS:
        lines.append(
            f"  Note: '{relationship}' is not a standard relationship ({', '.join(KNOWN_RELATIONSHIPS)})"
        )
    # check for duplicate
    for r in store["experiment_links"]:
        if r["from_experiment"] == from_id and r["to_experiment"] == to_id:
            return effects_fail(f"Link already exists: {from_id} -> {to_id}")
    if relationship == "uses_prep_from" and to_exp["type"] != "purification":
        lines.append(f"  Note: {to_id} is type '{to_exp['type']}', not purification.")
    # insert link
    new_link = {
        "from_experiment": from_id,
        "to_experiment": to_id,
        "relationship": relationship,
    }
    new_store = store_append_row(store, "experiment_links", new_link)
    # auto-populate assay_details.protein_prep if uses_prep_from
    if relationship == "uses_prep_from":
        existing_assay = [
            r for r in new_store["assay_details"] if r["experiment_id"] == from_id
        ]
        if existing_assay:
            # update only if protein_prep is NULL
            if existing_assay[0]["protein_prep"] is None:
                new_store, _ = store_update_row(
                    new_store,
                    "assay_details",
                    "experiment_id",
                    from_id,
                    {"protein_prep": to_id},
                )
        else:
            new_assay = {
                "experiment_id": from_id,
                "protein_prep": to_id,
                "dna_template": None,
                "dna_prep_method": None,
                "conditions": None,
            }
            new_store = store_append_row(new_store, "assay_details", new_assay)
    lines.append(f"Linked: {from_id} -> {to_id} ({relationship})")
    return effects_ok(stdout=lines, store=new_store)


def transform_exp_addstrain(store, args, clock):
    """Pure transform: associate a strain with an experiment.

    args: {"exp_id": str, "strain_id": str, "role": str | None}
    """
    raw_id = args["exp_id"]
    strain_id = args["strain_id"]
    role = args.get("role")
    exp_id = normalize_exp_id(raw_id)
    if not exp_id:
        return effects_fail(f"Invalid experiment ID: '{raw_id}'")
    if not lookup_experiment(store, exp_id):
        return effects_fail(f"No experiment found with ID {exp_id}")
    # validate strain exists
    if not lookup_strain(store, strain_id):
        return effects_fail(
            [
                f"No strain found with ID {strain_id}",
                "Register it first: python lw.py strain add <id> <genotype>",
            ]
        )
    # check for duplicate
    for r in store["experiment_strains"]:
        if r["experiment_id"] == exp_id and r["strain_id"] == strain_id:
            return effects_fail(f"Strain {strain_id} already linked to {exp_id}")
    new_row = {
        "experiment_id": exp_id,
        "strain_id": strain_id,
        "role": role,
    }
    new_store = store_append_row(store, "experiment_strains", new_row)
    role_str = f" as {role}" if role else ""
    return effects_ok(
        stdout=[f"Added: {strain_id} -> {exp_id}{role_str}"],
        store=new_store,
    )


def transform_strain_add(store, args, clock):
    """Pure transform: register a new strain.

    args: {"strain_id": str, "genotype": str}
    """
    strain_id = args["strain_id"]
    genotype = args["genotype"]
    ok, reason = validate_name(strain_id, "strain ID")
    if not ok:
        return effects_fail(reason)
    # check for duplicate
    if lookup_strain(store, strain_id):
        return effects_fail(
            [
                f"Strain {strain_id} already exists.",
                f"Use 'python lw.py strain show {strain_id}' to view it.",
            ]
        )
    new_row = {
        "id": strain_id,
        "genotype": genotype,
        "label": None,
        "parent": None,
        "construction": None,
        "selection_marker": None,
        "verification": None,
        "storage_location": None,
        "notes": None,
    }
    new_store = store_append_row(store, "strains", new_row)
    lines = [
        f"Registered: {strain_id}",
        f"Genotype:   {genotype}",
    ]
    if not re.match(r"^[A-Za-z]+\d+$", strain_id):
        lines.append(
            f"Note: '{strain_id}' doesn't follow typical strain ID format (letters + digits, e.g., LY456)"
        )
    lines.append(f"Set a label: lw strain update {strain_id} label <short-name>")
    return effects_ok(stdout=lines, store=new_store)


def transform_strain_update(store, args, clock):
    """Pure transform: update a field on a strain.

    args: {"strain_id": str, "field": str, "value": str}
    """
    strain_id = args["strain_id"]
    field = args["field"]
    value = args["value"]
    strain = lookup_strain(store, strain_id)
    if not strain:
        return effects_fail(f"No strain found with ID {strain_id}")
    if field not in STRAIN_FIELDS:
        return effects_fail(
            [
                f"Unknown field: '{field}'",
                f"Valid fields: {', '.join(STRAIN_FIELDS)}",
            ]
        )
    old_val = strain.get(field)
    original_value = value
    if field == "notes":
        value = append_note(old_val, value, clock)
    new_store, _ = store_update_row(store, "strains", "id", strain_id, {field: value})
    lines = [f"Updated {strain_id}"]
    if field == "notes":
        lines.append(f"  {field}: appended entry")
    else:
        lines.append(f"  {field}: {old_val if old_val else '-'} -> {value}")
    return effects_ok(stdout=lines, store=new_store)


# --- Phase 3B: Read-only transforms ---


def transform_strain_show(store, args, clock):
    """Pure transform: display strain details.

    args: {"strain_id": str}
    """
    strain_id = args["strain_id"]
    strain = lookup_strain(store, strain_id)
    if not strain:
        return effects_fail(f"No strain found with ID {strain_id}")
    lines = []
    fields = [
        ("ID", strain["id"]),
        ("Genotype", strain["genotype"]),
        ("Label", strain["label"]),
        ("Parent", strain["parent"]),
        ("Construction", strain["construction"]),
        ("Marker", strain["selection_marker"]),
        ("Verification", strain["verification"]),
        ("Storage", strain["storage_location"]),
        ("Notes", strain["notes"]),
    ]
    label_width = max(len(f[0]) for f in fields)
    for label, val in fields:
        lines.append(f"  {label + ':':<{label_width + 2}} {val if val else '-'}")
    # experiments using this strain
    exps = []
    for es in store["experiment_strains"]:
        if es["strain_id"] == strain_id:
            exp = lookup_experiment(store, es["experiment_id"])
            if exp:
                exps.append((exp["id"], exp["type"], exp["shortname"], es["role"]))
    if exps:
        lines.append("")
        lines.append("  Experiments:")
        for eid, etype, sname, role in exps:
            role_str = f" ({role})" if role else ""
            lines.append(f"    {eid}  {etype:<14} {sname}{role_str}")
    return effects_ok(stdout=lines)


def transform_strain_list(store, args, clock):
    """Pure transform: list all strains.

    args: {"no_label_only": bool}
    """
    no_label_only = args.get("no_label_only", False)
    strains = store["strains"]
    if no_label_only:
        strains = [s for s in strains if s["label"] is None]
    if not strains:
        if no_label_only:
            return effects_ok(stdout=["All strains have labels."])
        else:
            return effects_ok(stdout=["No strains registered."])
    # compute experiment counts per strain
    exp_counts = {}
    for es in store["experiment_strains"]:
        exp_counts[es["strain_id"]] = exp_counts.get(es["strain_id"], 0) + 1
    lines = []
    lines.append(f"  {'ID':<12} {'Label':<16} {'Genotype':<43} {'Exps':>4}")
    lines.append("  " + "-" * 77)
    n_missing = 0
    # sort by id
    for s in sorted(strains, key=lambda r: r["id"]):
        sid = s["id"]
        label = s["label"]
        genotype = s["genotype"]
        exp_count = exp_counts.get(sid, 0)
        if label:
            lbl_str = label
        else:
            lbl_str = "*NO LABEL*"
            n_missing += 1
        geno_str = (genotype[:40] + "...") if len(genotype) > 40 else genotype
        lines.append(f"  {sid:<12} {lbl_str:<16} {geno_str:<43} {exp_count:>4}")
    lines.append("")
    if no_label_only:
        lines.append(
            f"  {len(strains)} strain{'s' if len(strains) != 1 else ''} without labels"
        )
    else:
        msg = f"  {len(strains)} strain{'s' if len(strains) != 1 else ''}"
        if n_missing > 0:
            msg += f", {n_missing} missing label{'s' if n_missing != 1 else ''}"
        lines.append(msg)
    return effects_ok(stdout=lines)


def _parse_find_args(pos_args):
    """Parse exp find arguments. Returns (query, f_type, f_status, f_strain) or error string."""
    query_parts = []
    f_type = None
    f_status = None
    f_strain = None
    i = 0
    while i < len(pos_args):
        if pos_args[i] == "--type" and i + 1 < len(pos_args):
            f_type = pos_args[i + 1]
            i += 2
        elif pos_args[i] == "--status" and i + 1 < len(pos_args):
            f_status = pos_args[i + 1]
            i += 2
        elif pos_args[i] == "--strain" and i + 1 < len(pos_args):
            f_strain = pos_args[i + 1]
            i += 2
        elif not pos_args[i].startswith("--"):
            query_parts.append(pos_args[i])
            i += 1
        else:
            return None, None, None, None, f"Unknown flag: {pos_args[i]}"
    query = " ".join(query_parts) if query_parts else None
    return query, f_type, f_status, f_strain, None


def transform_exp_find(store, args, clock):
    """Pure transform: search experiments.

    args: {"pos_args": list}
    Uses case-insensitive matching (replicating SQLite LIKE behavior for ASCII).
    """
    pos_args = args["pos_args"]
    query, f_type, f_status, f_strain, err = _parse_find_args(pos_args)
    if err:
        return effects_fail(err)
    if f_type and f_type not in VALID_TYPES:
        return effects_fail(
            [
                f"Unknown type: '{f_type}'",
                f"Valid types: {', '.join(VALID_TYPES)}",
            ]
        )
    # filter experiments
    results = []
    # build strain set for filtering
    strain_exp_ids = None
    if f_strain:
        strain_exp_ids = set()
        for es in store["experiment_strains"]:
            if es["strain_id"] == f_strain:
                strain_exp_ids.add(es["experiment_id"])
    for e in store["experiments"]:
        if f_type and e["type"] != f_type:
            continue
        if f_status and e["status"] != f_status:
            continue
        if strain_exp_ids is not None and e["id"] not in strain_exp_ids:
            continue
        if query:
            q_lower = query.lower()
            fields_to_search = [
                (e.get("id") or "").lower(),
                (e.get("shortname") or "").lower(),
                (e.get("title") or "").lower(),
                (e.get("notes") or "").lower(),
            ]
            if not any(q_lower in f for f in fields_to_search):
                continue
        results.append(e)
    if not results:
        lines = []
        if f_strain and not lookup_strain(store, f_strain):
            return effects_fail(f"Strain '{f_strain}' not found in database.")
        lines.append("No experiments found.")
        return effects_ok(stdout=lines)
    # sort: date_started desc, id desc
    results.sort(key=lambda e: (e["date_started"], e["id"]), reverse=True)
    lines = []
    lines.append(f"{'ID':<10} {'Type':<14} {'Name':<24} {'Status':<10} {'Started'}")
    lines.append("-" * 72)
    for e in results:
        display_name = e["title"] if e["title"] != e["shortname"] else e["shortname"]
        lines.append(
            f"{e['id']:<10} {e['type']:<14} {display_name:<24} {e['status']:<10} {e['date_started']}"
        )
    lines.append("")
    lines.append(f"{len(results)} experiment{'s' if len(results) != 1 else ''}")
    return effects_ok(stdout=lines)


def _parse_list_args(pos_args):
    """Parse exp list arguments. Returns (f_type, error_string)."""
    f_type = None
    i = 0
    while i < len(pos_args):
        if pos_args[i] == "--type" and i + 1 < len(pos_args):
            f_type = pos_args[i + 1]
            i += 2
        else:
            return None, f"Unknown argument: {pos_args[i]}"
    return f_type, None


def transform_exp_list(store, args, clock):
    """Pure transform: list experiments grouped by type.

    args: {"pos_args": list}
    """
    pos_args = args["pos_args"]
    f_type, err = _parse_list_args(pos_args)
    if err:
        return effects_fail(err)
    if f_type and f_type not in VALID_TYPES:
        return effects_fail(
            [
                f"Unknown type: '{f_type}'",
                f"Valid types: {', '.join(VALID_TYPES)}",
            ]
        )
    experiments = store["experiments"]
    if f_type:
        experiments = [e for e in experiments if e["type"] == f_type]
    if not experiments:
        if f_type:
            return effects_ok(stdout=[f"No {f_type} experiments found."])
        else:
            return effects_ok(stdout=["No experiments registered."])
    # sort: date_started desc, id desc
    experiments = sorted(
        experiments, key=lambda e: (e["date_started"], e["id"]), reverse=True
    )
    # group by type
    groups = {}
    for e in experiments:
        display_name = e["title"] if e["title"] != e["shortname"] else e["shortname"]
        groups.setdefault(e["type"], []).append(
            (e["id"], display_name, e["status"], e["date_started"])
        )
    type_order = [t for t in VALID_TYPES if t in groups]
    lines = []
    total = 0
    type_counts = []
    for etype in type_order:
        exps = groups[etype]
        lines.append(f"{etype} ({len(exps)})")
        for eid, sname, estatus, started in exps:
            detail = ""
            if etype == "purification":
                pur = [
                    r
                    for r in store["purification_details"]
                    if r["experiment_id"] == eid
                ]
                if pur and pur[0].get("protein"):
                    detail = pur[0]["protein"]
                    if pur[0].get("yield_mg") is not None:
                        detail += f" / {pur[0]['yield_mg']} mg"
                else:
                    detail = "-"
            elif etype in ("loading", "gelshift", "atpase"):
                s_count = sum(1 for s in store["samples"] if s["experiment_id"] == eid)
                parts = [f"{s_count} samples"]
                assay = [r for r in store["assay_details"] if r["experiment_id"] == eid]
                if assay and assay[0].get("protein_prep"):
                    prep_id = assay[0]["protein_prep"]
                    pur = [
                        r
                        for r in store["purification_details"]
                        if r["experiment_id"] == prep_id
                    ]
                    protein = pur[0]["protein"] if pur else None
                    if protein:
                        parts.append(f"prep={prep_id} ({protein})")
                    else:
                        parts.append(f"prep={prep_id}")
                detail = "  ".join(parts)
            if detail:
                lines.append(
                    f"    {eid:<10} {sname:<20} {estatus:<10} {started}  {detail}"
                )
            else:
                lines.append(f"    {eid:<10} {sname:<20} {estatus:<10} {started}")
        total += len(exps)
        type_counts.append(f"{etype}: {len(exps)}")
        lines.append("")
    if f_type:
        lines.append(f"{total} experiment{'s' if total != 1 else ''}")
    else:
        lines.append(
            f"{total} experiment{'s' if total != 1 else ''} ({', '.join(type_counts)})"
        )
    return effects_ok(stdout=lines)


def format_size(size):
    """Pure: format a byte count as a human-readable string."""
    if size > 1024 * 1024:
        return f"{size / (1024 * 1024):.1f} MB"
    elif size > 1024:
        return f"{size / 1024:.1f} KB"
    else:
        return f"{size} B"


def prepare_stage_list_io(staging_dir):
    """IMPURE (labeled): read filesystem state for stage list.

    Returns a list of {"name": str, "size": int} dicts, or None if dir missing.
    """
    if not os.path.exists(staging_dir):
        return None
    files = [f for f in os.listdir(staging_dir) if not f.startswith(".")]
    return [
        {"name": f, "size": os.path.getsize(os.path.join(staging_dir, f))}
        for f in sorted(files)
    ]


def transform_stage_list(store, args, clock):
    """Pure transform: display staging files and pending expectations.

    args: {"file_list": list[dict] | None}
    file_list is the output of prepare_stage_list_io (None = dir missing).
    """
    file_list = args["file_list"]
    if file_list is None:
        return effects_fail("Staging directory does not exist.")
    lines = []
    if not file_list:
        lines.append("Staging is empty.")
    else:
        lines.append(f"Staging ({len(file_list)} files):")
        lines.append("")
        for f in file_list:
            fname = f["name"]
            size_str = format_size(f["size"])
            parsed_eid, parsed_slot = parse_staging_name(fname)
            if parsed_eid:
                slot_str = f" s{parsed_slot}" if parsed_slot else ""
                match_str = f"  -> {parsed_eid}{slot_str}"
            else:
                match_str = ""
            lines.append(f"  {fname:<40} {size_str:>10}{match_str}")
    # pending expectations from store
    pending_rows = sorted(
        [r for r in store["stage_expectations"] if r["status"] == "pending"],
        key=lambda r: (r["experiment_id"], r["slot"]),
    )
    if pending_rows:
        staging_matches = set()
        for f in file_list or []:
            parsed_eid, parsed_slot = parse_staging_name(f["name"])
            if parsed_eid is not None and parsed_slot is not None:
                staging_matches.add((parsed_eid, parsed_slot))
        lines.append("")
        lines.append("Pending expectations:")
        n_matched = 0
        for r in pending_rows:
            eid, slot, desc = r["experiment_id"], r["slot"], r["descriptor"]
            if (eid, slot) in staging_matches:
                lines.append(f"  {eid} s{slot}  {desc:<20} (matched)")
                n_matched += 1
            else:
                lines.append(f"  {eid} s{slot}  {desc:<20} (no file)")
        lines.append(f"{len(pending_rows)} pending, {n_matched} matched")
    return effects_ok(stdout=lines)


def prepare_status_io(staging_dir, lab_notes_path):
    """IMPURE (labeled): read filesystem state for status command.

    Returns {"staging_files": list[str], "staging_total_size": int,
             "lab_notes_date": str | None}.
    """
    staging_files = []
    staging_total_size = 0
    if os.path.exists(staging_dir):
        staging_files = [f for f in os.listdir(staging_dir) if not f.startswith(".")]
        staging_total_size = sum(
            os.path.getsize(os.path.join(staging_dir, f)) for f in staging_files
        )
    lab_notes_date = None
    if os.path.exists(lab_notes_path):
        mtime = os.path.getmtime(lab_notes_path)
        lab_notes_date = datetime.datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")
    return {
        "staging_files": staging_files,
        "staging_total_size": staging_total_size,
        "lab_notes_date": lab_notes_date,
    }


def transform_status(store, args, clock):
    """Pure transform: display system dashboard.

    args: {"fs": dict} where fs is the output of prepare_status_io.
    """
    fs = args["fs"]
    staging_files = fs["staging_files"]
    staging_total_size = fs["staging_total_size"]
    lab_notes_date = fs["lab_notes_date"]
    lines = []
    # --- experiments ---
    all_exps = store["experiments"]
    n_total_exp = len(all_exps)
    active_exps = [e for e in all_exps if e["status"] == "active"]
    n_active = len(active_exps)
    n_complete = sum(1 for e in all_exps if e["status"] == "complete")
    lines.append("Experiments:")
    lines.append(f"  {n_total_exp} total ({n_active} active, {n_complete} complete)")
    if n_active > 0:
        active_sorted = sorted(
            active_exps, key=lambda e: e["date_started"], reverse=True
        )
        lines.append("")
        for e in active_sorted[:10]:
            days = (
                clock["today"] - datetime.date.fromisoformat(e["date_started"])
            ).days
            lines.append(
                f"    {e['id']:<10} {e['type']:<14} {e['shortname']:<20} ({days}d)"
            )
        if n_active > 10:
            lines.append(f"    ... and {n_active - 10} more")
    # --- staging ---
    lines.append("")
    lines.append("Staging:")
    if staging_files:
        size_str = format_size(staging_total_size)
        lines.append(
            f"  {len(staging_files)} file{'s' if len(staging_files) != 1 else ''}"
            f" ({size_str})"
        )
    else:
        lines.append("  clean")
    # --- pending expectations ---
    n_pending = sum(1 for r in store["stage_expectations"] if r["status"] == "pending")
    if n_pending > 0:
        lines.append(
            f"  {n_pending} pending expectation{'s' if n_pending != 1 else ''}"
        )
    # --- strains ---
    lines.append("")
    n_strains = len(store["strains"])
    n_no_label = sum(1 for s in store["strains"] if s["label"] is None)
    lines.append("Strains:")
    strain_msg = f"  {n_strains} registered"
    if n_no_label > 0:
        strain_msg += f", {n_no_label} missing label{'s' if n_no_label != 1 else ''}"
    lines.append(strain_msg)
    # --- samples ---
    lines.append("")
    n_samples = len(store["samples"])
    lines.append("Samples:")
    lines.append(f"  {n_samples} total")
    # --- recent activity ---
    lines.append("")
    lines.append("Recent activity:")
    dates_started = [e["date_started"] for e in all_exps if e["date_started"]]
    dates_completed = [e["date_completed"] for e in all_exps if e["date_completed"]]
    last_started = max(dates_started) if dates_started else None
    last_completed = max(dates_completed) if dates_completed else None
    if last_started:
        lines.append(f"  Last experiment created: {last_started}")
    if last_completed:
        lines.append(f"  Last experiment completed: {last_completed}")
    if lab_notes_date:
        lines.append(f"  Lab notes last modified: {lab_notes_date}")
    if not last_started and not last_completed:
        lines.append("  No experiments yet.")
    # --- attention needed ---
    warnings = []
    if n_no_label > 0:
        warnings.append(
            f"{n_no_label} strain{'s' if n_no_label != 1 else ''} without labels"
        )
    if staging_files:
        warnings.append(
            f"{len(staging_files)} file{'s' if len(staging_files) != 1 else ''} in staging"
        )
    if n_pending > 0:
        warnings.append(
            f"{n_pending} pending expectation{'s' if n_pending != 1 else ''}"
        )
    n_no_strains = sum(
        1
        for e in all_exps
        if e["status"] == "active"
        and not any(
            es["experiment_id"] == e["id"] for es in store["experiment_strains"]
        )
    )
    if n_no_strains > 0:
        warnings.append(
            f"{n_no_strains} active experiment"
            f"{'s' if n_no_strains != 1 else ''}"
            " with no strains linked"
        )
    n_no_notes = sum(
        1 for e in all_exps if e["status"] == "active" and e["notes"] is None
    )
    if n_no_notes > 0:
        warnings.append(
            f"{n_no_notes} active experiment"
            f"{'s' if n_no_notes != 1 else ''}"
            " with no notes"
        )
    if warnings:
        lines.append("")
        lines.append("Attention needed:")
        for w in warnings:
            lines.append(f"  - {w}")
    return effects_ok(stdout=lines)


# --- Phase 3C: Complex/IO transforms ---


def transform_expect_list(store, args, clock):
    """Pure transform: list stage expectations.

    args: {"filter_exp_id": str | None}
    """
    filter_exp_id = args.get("filter_exp_id")
    if filter_exp_id:
        exp = lookup_experiment(store, filter_exp_id)
        if not exp:
            return effects_fail(f"No experiment found with ID {filter_exp_id}")
    expectations = store["stage_expectations"]
    if filter_exp_id:
        rows = [r for r in expectations if r["experiment_id"] == filter_exp_id]
        rows.sort(key=lambda r: r["slot"])
    else:
        rows = sorted(expectations, key=lambda r: (r["experiment_id"], r["slot"]))
    if not rows:
        if filter_exp_id:
            return effects_ok(stdout=[f"No expectations for {filter_exp_id}."])
        else:
            return effects_ok(stdout=["No expectations registered."])
    lines = []
    n_pending = 0
    n_filed = 0
    if filter_exp_id:
        for r in rows:
            status = r["status"]
            slot = r["slot"]
            desc = r["descriptor"]
            if status == "pending":
                created_date = (r["created_at"] or "unknown")[:10]
                lines.append(
                    f"  s{slot}  {desc:<20} (pending, registered {created_date})"
                )
                n_pending += 1
            elif status == "filed":
                filed_date = (r["filed_at"] or "unknown")[:10]
                lines.append(f"  s{slot}  {desc:<20} (filed {filed_date})")
                n_filed += 1
            elif status == "cancelled":
                lines.append(f"  s{slot}  {desc:<20} (cancelled)")
    else:
        current_eid = None
        for r in rows:
            eid = r["experiment_id"]
            status = r["status"]
            slot = r["slot"]
            desc = r["descriptor"]
            if eid != current_eid:
                current_eid = eid
                lines.append(f"{eid}:")
            if status == "pending":
                created_date = (r["created_at"] or "unknown")[:10]
                lines.append(
                    f"  s{slot}  {desc:<20} (pending, registered {created_date})"
                )
                n_pending += 1
            elif status == "filed":
                filed_date = (r["filed_at"] or "unknown")[:10]
                lines.append(f"  s{slot}  {desc:<20} (filed {filed_date})")
                n_filed += 1
            elif status == "cancelled":
                lines.append(f"  s{slot}  {desc:<20} (cancelled)")
    lines.append(f"{n_pending} pending, {n_filed} filed")
    return effects_ok(stdout=lines)


def transform_expect_cancel(store, args, clock):
    """Pure transform: cancel stage expectations.

    args: {"exp_id": str, "slot": int | None}
    If slot is None, cancels all pending for this experiment.
    """
    exp_id = args["exp_id"]
    slot = args.get("slot")
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    new_store = copy.deepcopy(store)
    lines = []
    if slot is not None:
        # find specific slot
        target = None
        for r in new_store["stage_expectations"]:
            if r["experiment_id"] == exp_id and r["slot"] == slot:
                target = r
                break
        if not target:
            return effects_fail(f"Slot s{slot} not found for {exp_id}")
        if target["status"] != "pending":
            return effects_fail(f"Slot s{slot} is already {target['status']}.")
        target["status"] = "cancelled"
        lines.append(f"Cancelled: {exp_id} s{slot} ({target['descriptor']})")
    else:
        # cancel all pending
        pending = [
            r
            for r in new_store["stage_expectations"]
            if r["experiment_id"] == exp_id and r["status"] == "pending"
        ]
        if not pending:
            return effects_fail(f"No pending expectations for {exp_id}")
        for r in pending:
            r["status"] = "cancelled"
        pending.sort(key=lambda r: r["slot"])
        if len(pending) == 1:
            s, d = pending[0]["slot"], pending[0]["descriptor"]
            lines.append(f"Cancelled 1 pending expectation for {exp_id}:")
            lines.append(f"  s{s}  {d}")
        else:
            lines.append(f"Cancelled {len(pending)} pending expectations for {exp_id}:")
            for r in pending:
                lines.append(f"  s{r['slot']}  {r['descriptor']}")
    return effects_ok(stdout=lines, store=new_store)


def transform_expect_register(store, args, clock):
    """Pure transform: register new stage expectations.

    args: {"exp_id": str, "descriptors": list[str]}
    Descriptors should already be validated and lowercased.
    """
    exp_id = args["exp_id"]
    descriptors = args["descriptors"]
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    # check for duplicate descriptors within arguments
    seen = set()
    for desc in descriptors:
        if desc in seen:
            return effects_fail(f"Duplicate descriptor in arguments: '{desc}'")
        seen.add(desc)
    # check against existing pending expectations
    existing_descs = {
        r["descriptor"]
        for r in store["stage_expectations"]
        if r["experiment_id"] == exp_id and r["status"] == "pending"
    }
    for desc in descriptors:
        if desc in existing_descs:
            return effects_fail(
                [
                    f"Descriptor '{desc}' already pending for {exp_id}.",
                    "Cancel it first or use a different descriptor.",
                ]
            )
    # get max slot
    exp_slots = [
        r["slot"] for r in store["stage_expectations"] if r["experiment_id"] == exp_id
    ]
    max_slot = max(exp_slots) if exp_slots else 0
    new_store = copy.deepcopy(store)
    lines = []
    slots = []
    created_at = clock["now"].isoformat(timespec="seconds")
    for i, desc in enumerate(descriptors, 1):
        slot = max_slot + i
        new_store["stage_expectations"].append(
            {
                "id": None,  # AUTOINCREMENT filled by commit
                "experiment_id": exp_id,
                "slot": slot,
                "descriptor": desc,
                "status": "pending",
                "created_at": created_at,
                "filed_at": None,
            }
        )
        slots.append((slot, desc))
        lines.append(f"  s{slot}  {desc}")
    # instrument names
    if len(slots) == 1:
        s, d = slots[0]
        lines.append(f"  At instrument: lemr {exp_id} s{s}")
    else:
        lines.append("  At instrument:")
        for s, d in slots:
            lines.append(f"    lemr {exp_id} s{s}")
    # count total pending
    total_pending = sum(
        1
        for r in new_store["stage_expectations"]
        if r["experiment_id"] == exp_id and r["status"] == "pending"
    )
    lines.append(f"  {total_pending} total pending for {exp_id}")
    return effects_ok(stdout=lines, store=new_store)


def transform_exp_delete_plan(store, args, clock):
    """Pure transform: plan experiment deletion (dry-run or pre-confirm analysis).

    args: {"exp_id": str, "confirm": bool, "keep_folder": bool,
           "folder_exists": bool, "folder_contents": list[str]}
    Returns effects with stdout lines describing what will/would be deleted.
    When confirm=True, also returns store with cascaded deletions applied.
    """
    exp_id = args["exp_id"]
    confirm = args.get("confirm", False)
    keep_folder = args.get("keep_folder", False)
    folder_exists = args.get("folder_exists", False)
    folder_contents = args.get("folder_contents", [])
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    # gather cascade targets
    cascade = []
    es_rows = [r for r in store["experiment_strains"] if r["experiment_id"] == exp_id]
    if es_rows:
        details = []
        for r in es_rows:
            details.append(
                f"  {r['strain_id']}" + (f" ({r['role']})" if r.get("role") else "")
            )
        cascade.append(("experiment_strains", len(es_rows), details))
    el_rows = [
        r
        for r in store["experiment_links"]
        if r["from_experiment"] == exp_id or r["to_experiment"] == exp_id
    ]
    if el_rows:
        details = [
            f"  {r['from_experiment']} -> {r['to_experiment']} ({r['relationship']})"
            for r in el_rows
        ]
        cascade.append(("experiment_links", len(el_rows), details))
    f_rows = [r for r in store["files"] if r["experiment_id"] == exp_id]
    if f_rows:
        details = [f"  {r['filename']}" for r in f_rows]
        cascade.append(("files", len(f_rows), details))
    pur_rows = [
        r for r in store["purification_details"] if r["experiment_id"] == exp_id
    ]
    if pur_rows:
        cascade.append(("purification_details", 1, []))
    assay_rows = [r for r in store["assay_details"] if r["experiment_id"] == exp_id]
    if assay_rows:
        cascade.append(("assay_details", 1, []))
    fp_rows = [r for r in store["figure_panels"] if r["experiment_id"] == exp_id]
    if fp_rows:
        details = [f"  {r['figure']} panel {r['panel']}" for r in fp_rows]
        cascade.append(("figure_panels", len(fp_rows), details))
    s_rows = [r for r in store["samples"] if r["experiment_id"] == exp_id]
    if s_rows:
        cascade.append(("samples", len(s_rows), []))
    se_rows = [r for r in store["stage_expectations"] if r["experiment_id"] == exp_id]
    if se_rows:
        cascade.append(("stage_expectations", len(se_rows), []))
    # downstream protein_prep references from other experiments
    prep_ref_rows = [
        r for r in store["assay_details"] if r.get("protein_prep") == exp_id
    ]
    lines = []
    mode = "DELETE" if confirm else "DRY RUN"
    lines.append(
        f"[{mode}] {exp_id}: {exp['type']} / {exp['shortname']} ({exp['status']})"
    )
    lines.append(
        f"  Folder: {exp['folder_name']}/ ({'exists' if folder_exists else 'missing'}"
        + (f", {len(folder_contents)} files" if folder_contents else "")
        + ")"
    )
    if cascade:
        lines.append("")
        lines.append("  Related records:")
        for table, count, details in cascade:
            lines.append(f"    {table}: {count} row{'s' if count != 1 else ''}")
            for d in details:
                lines.append(f"      {d}")
    if prep_ref_rows:
        lines.append("")
        lines.append("  Downstream references (protein_prep will be cleared):")
        for r in prep_ref_rows:
            lines.append(
                f"    {r['experiment_id']} assay_details.protein_prep -> {exp_id}"
            )
    if not confirm:
        lines.append("")
        lines.append(f"  This is a dry run. To delete, run:")
        lines.append(f"    python lw.py exp delete {exp_id} --confirm")
        if folder_exists:
            lines.append(
                f"    python lw.py exp delete {exp_id} --confirm --keep-folder  (keep files)"
            )
        return effects_ok(stdout=lines)
    # --- actual deletion (in store) ---
    new_store = copy.deepcopy(store)
    # clear downstream protein_prep references
    if prep_ref_rows:
        for r in new_store["assay_details"]:
            if r.get("protein_prep") == exp_id:
                r["protein_prep"] = None
        lines.append(
            f"  Cleared {len(prep_ref_rows)} protein_prep"
            f" reference{'s' if len(prep_ref_rows) != 1 else ''}"
        )
    # cascade deletes
    new_store["stage_expectations"] = [
        r for r in new_store["stage_expectations"] if r["experiment_id"] != exp_id
    ]
    new_store["samples"] = [
        r for r in new_store["samples"] if r["experiment_id"] != exp_id
    ]
    new_store["figure_panels"] = [
        r for r in new_store["figure_panels"] if r["experiment_id"] != exp_id
    ]
    new_store["files"] = [r for r in new_store["files"] if r["experiment_id"] != exp_id]
    new_store["assay_details"] = [
        r for r in new_store["assay_details"] if r["experiment_id"] != exp_id
    ]
    new_store["purification_details"] = [
        r for r in new_store["purification_details"] if r["experiment_id"] != exp_id
    ]
    new_store["experiment_strains"] = [
        r for r in new_store["experiment_strains"] if r["experiment_id"] != exp_id
    ]
    new_store["experiment_links"] = [
        r
        for r in new_store["experiment_links"]
        if r["from_experiment"] != exp_id and r["to_experiment"] != exp_id
    ]
    new_store["experiments"] = [
        r for r in new_store["experiments"] if r["id"] != exp_id
    ]
    total_rows = sum(c for _, c, _ in cascade) + 1
    lines.append("")
    lines.append(f"  Deleted {total_rows} database row{'s' if total_rows != 1 else ''}")
    # folder handling is IO; we emit instructions for the caller
    if folder_exists and not keep_folder:
        lines.append(f"  Removed folder: {exp['folder_name']}/")
    elif folder_exists and keep_folder:
        lines.append(f"  Folder kept: {exp['folder_name']}/")
    lines.append("")
    lines.append(f"Deleted: {exp_id}")
    return effects_ok(stdout=lines, store=new_store)


def transform_stage_auto_plan(store, args, clock):
    """Pure transform: plan stage auto-assignment.

    args: {
        "file_list": list[{"name": str, "src_path": str}] | None,
        "confirm": bool,
    }
    file_list is None if staging dir missing, [] if empty.

    Returns effects with stdout (dry-run or confirm report).
    For confirm mode, the actual file moves are done by the caller after
    inspecting the returned store's "_auto_plan" metadata.
    """
    file_list = args.get("file_list")
    confirm = args.get("confirm", False)
    if file_list is None:
        return effects_fail("Staging directory does not exist.")
    pending = [r for r in store["stage_expectations"] if r["status"] == "pending"]
    if not file_list:
        if pending:
            n_pend = len(pending)
            return effects_ok(
                stdout=[
                    f"Staging is empty. {n_pend} expectation"
                    f"{'s' if n_pend != 1 else ''} still pending."
                ]
            )
        else:
            return effects_ok(stdout=["Staging is empty."])
    if not pending:
        lines = ["No pending expectations. Use 'stage assign' for ad hoc files."]
        lines.append("")
        lines.append(f"Unmatched ({len(file_list)}):")
        for f in sorted(file_list, key=lambda f: f["name"]):
            lines.append(f"  {f['name']}")
        return effects_ok(stdout=lines)
    # build lookup structures
    pending_by_key = {}
    pending_by_exp = {}
    for r in pending:
        key = (r["experiment_id"], r["slot"])
        pending_by_key[key] = r["descriptor"]
        pending_by_exp.setdefault(r["experiment_id"], []).append(
            (r["slot"], r["descriptor"])
        )
    # match each file
    matched_groups = {}
    unmatched = []
    skipped = []
    inferred_files = set()
    for f in sorted(file_list, key=lambda f: f["name"]):
        fname = f["name"]
        parsed_eid, parsed_slot = parse_staging_name(fname)
        if parsed_eid is not None and parsed_slot is not None:
            key = (parsed_eid, parsed_slot)
            if key in pending_by_key:
                if key not in matched_groups:
                    matched_groups[key] = {
                        "descriptor": pending_by_key[key],
                        "files": [],
                    }
                matched_groups[key]["files"].append(fname)
            else:
                unmatched.append(fname)
        elif parsed_eid is not None:
            exp_pending = pending_by_exp.get(parsed_eid, [])
            if len(exp_pending) == 1:
                s, d = exp_pending[0]
                key = (parsed_eid, s)
                if key not in matched_groups:
                    matched_groups[key] = {"descriptor": d, "files": []}
                matched_groups[key]["files"].append(fname)
                inferred_files.add(fname)
            elif len(exp_pending) > 1:
                skipped.append(
                    (
                        fname,
                        f"{parsed_eid} found but slot ambiguous "
                        f"({len(exp_pending)} pending). "
                        f"Use 'stage assign' or rename with slot number.",
                    )
                )
            else:
                unmatched.append(fname)
        else:
            unmatched.append(fname)
    # resolve folder names for matched experiments
    folder_cache = {}
    valid_groups = {}
    for key, group in matched_groups.items():
        exp_id = key[0]
        if exp_id not in folder_cache:
            exp = lookup_experiment(store, exp_id)
            if not exp:
                for fname in group["files"]:
                    skipped.append(
                        (fname, f"Experiment {exp_id} not found in database.")
                    )
                continue
            folder_cache[exp_id] = exp["folder_name"]
        group["folder_name"] = folder_cache[exp_id]
        valid_groups[key] = group
    matched_groups = valid_groups
    # build plan: compute destination names
    today_compact = clock["today"].strftime("%Y%m%d")
    matched_files = []
    for key, group in matched_groups.items():
        exp_id, slot = key
        desc = group["descriptor"]
        folder_name = group["folder_name"]
        for fname in group["files"]:
            ext = os.path.splitext(fname)[1].lower()
            dest_name = f"{today_compact}_{exp_id}_{desc}{ext}"
            rel_path = os.path.join("experiments", folder_name, dest_name)
            matched_files.append(
                {
                    "fname": fname,
                    "exp_id": exp_id,
                    "slot": slot,
                    "descriptor": desc,
                    "folder_name": folder_name,
                    "dest_name": dest_name,
                    "rel_path": rel_path,
                    "inferred": fname in inferred_files,
                }
            )
    lines = []
    if not confirm:
        # --- dry run ---
        lines.append("[DRY RUN]")
        if matched_files:
            lines.append(f"Matched ({len(matched_files)}):")
            for mf in matched_files:
                tag = " (slot inferred)" if mf["inferred"] else ""
                lines.append(f"  {mf['fname']}{tag}")
                lines.append(f"     -> {mf['rel_path']}")
        if skipped:
            lines.append(f"Skipped ({len(skipped)}):")
            for fname, reason in skipped:
                lines.append(f"  {fname}")
                lines.append(f"     {reason}")
        if unmatched:
            lines.append(f"Unmatched ({len(unmatched)}):")
            for fname in unmatched:
                lines.append(f"  {fname}")
        lines.append("")
        if matched_files:
            lines.append(
                f"{len(matched_files)} file"
                f"{'s' if len(matched_files) != 1 else ''}"
                " ready to assign. Run with --confirm to execute."
            )
        else:
            lines.append("No files matched any pending expectations.")
            if unmatched or skipped:
                lines.append("Use 'stage assign' for ad hoc files.")
        return effects_ok(stdout=lines)
    # --- confirm mode: return plan + updated store ---
    # The caller is responsible for the actual file moves and destination checks.
    # We return the plan in a special "_auto_plan" key on the effects dict.
    # For the store update: mark expectations as filed for matched slots.
    new_store = copy.deepcopy(store)
    filed_at = clock["now"].isoformat(timespec="seconds")
    matched_slot_keys = set()
    for mf in matched_files:
        matched_slot_keys.add((mf["exp_id"], mf["slot"]))
    for r in new_store["stage_expectations"]:
        key = (r["experiment_id"], r["slot"])
        if key in matched_slot_keys and r["status"] == "pending":
            r["status"] = "filed"
            r["filed_at"] = filed_at
    result = effects_ok(stdout=lines, store=new_store)
    result["_auto_plan"] = {
        "matched_files": matched_files,
        "skipped": skipped,
        "unmatched": unmatched,
    }
    return result


def transform_exp_manifest_plan(store, args, clock):
    """Pure transform: validate design and generate sample manifest.

    args: {
        "exp_id": str,
        "force": bool,
        "design": dict,   -- the loaded DESIGN dict from design.py
        "design_source": str,  -- display path for messages
    }
    Returns effects with stdout (summary table) and store with new samples.
    """
    exp_id = args["exp_id"]
    force = args.get("force", False)
    design = args["design"]
    design_source = args.get("design_source", "design.py")
    exp = lookup_experiment(store, exp_id)
    if not exp:
        return effects_fail(f"No experiment found with ID {exp_id}")
    # validate required keys
    for key in ("experiment_id", "expected_samples", "categories", "sort_order"):
        if key not in design:
            return effects_fail(f"Error: DESIGN missing required key '{key}'")
    if design["experiment_id"] != exp_id:
        return effects_fail(
            f"Error: DESIGN experiment_id '{design['experiment_id']}' does not match {exp_id}"
        )
    # validate sort_order matches category keys
    cat_keys = set(design["categories"].keys())
    sort_keys = set(design["sort_order"])
    if sort_keys != cat_keys:
        err_lines = []
        missing = cat_keys - sort_keys
        extra = sort_keys - cat_keys
        if missing:
            err_lines.append(
                f"Error: sort_order missing categories: {', '.join(missing)}"
            )
        if extra:
            err_lines.append(
                f"Error: sort_order references unknown categories: {', '.join(extra)}"
            )
        return effects_fail(err_lines)
    # check for existing samples
    existing = [s for s in store["samples"] if s["experiment_id"] == exp_id]
    if existing:
        if not force:
            return effects_fail(
                [
                    f"Error: {len(existing)} samples already exist for {exp_id}",
                    f"Use --force to delete and regenerate:",
                    f"  python lw.py exp manifest {exp_id} --force",
                ]
            )
    new_store = copy.deepcopy(store)
    # if force, remove existing samples
    if existing:
        new_store["samples"] = [
            s for s in new_store["samples"] if s["experiment_id"] != exp_id
        ]
    lines = [f"Loading design from {design_source}"]
    # resolve strain labels
    strain_labels = {}
    all_strain_ids = set()
    categories = design["categories"]
    if "strain" in categories:
        for s in categories["strain"]:
            if s is not None:
                all_strain_ids.add(s)
    for extra in design.get("extras", []):
        s = extra.get("strain")
        if s is not None:
            all_strain_ids.add(s)
    resolve_parts = []
    for sid in sorted(all_strain_ids):
        strain = lookup_strain(new_store, sid)
        if not strain:
            return effects_fail(
                [
                    f"Error: strain '{sid}' not found in database.",
                    f"Register it first: python lw.py strain add {sid} <genotype>",
                ]
            )
        lbl = strain["label"] if strain["label"] else sid
        strain_labels[sid] = lbl
        resolve_parts.append(f"{sid} -> {lbl}")
    if resolve_parts:
        lines.append(f"Resolving strains: {', '.join(resolve_parts)}")
    # generate extras
    samples = []
    pos = 1
    for extra in design.get("extras", []):
        sid = extra.get("strain")
        lbl = extra.get("label", "")
        if not lbl and sid:
            lbl = strain_labels.get(sid, sid)
        conditions = dict(extra)
        samples.append(
            {
                "position": pos,
                "sample_type": "extra",
                "strain_id": sid,
                "dna": extra.get("dna"),
                "label": lbl,
                "conditions": conditions,
            }
        )
        pos += 1
    n_extras = len(samples)
    # generate factorial cross
    cat_names = list(categories.keys())
    cat_values = [categories[c] for c in cat_names]
    factorial_rows = []
    for combo in itertools.product(*cat_values):
        row = dict(zip(cat_names, combo))
        excluded = False
        for exc_i, exc in enumerate(design.get("exclude", [])):
            try:
                if exc(row):
                    excluded = True
                    break
            except Exception as e:
                return effects_fail(f"Error in design.py exclude function {exc_i}: {e}")
        if not excluded:
            factorial_rows.append(row)

    # sort by sort_order
    def sort_key(row):
        keys = []
        for cat in design["sort_order"]:
            vals = categories[cat]
            v = row[cat]
            try:
                keys.append(vals.index(v))
            except ValueError:
                keys.append(len(vals))
        return tuple(keys)

    factorial_rows.sort(key=sort_key)
    # build sample records
    for row in factorial_rows:
        sid = row.get("strain")
        if sid is None:
            db_strain_id = None
            lbl = "no protein"
        else:
            db_strain_id = sid
            lbl = strain_labels.get(sid, sid)
        samples.append(
            {
                "position": pos,
                "sample_type": "factorial",
                "strain_id": db_strain_id,
                "dna": row.get("dna"),
                "label": lbl,
                "conditions": row,
            }
        )
        pos += 1
    n_factorial = len(samples) - n_extras
    n_total = len(samples)
    # verify count
    expected = design["expected_samples"]
    if n_total != expected:
        err_lines = [f"Error: generated {n_total} samples (expected {expected})"]
        err_lines.append(f"  {n_extras} extras + {n_factorial} factorial = {n_total}")
        err_lines.append("")
        err_lines.append("Factorial rows generated:")
        for i, row in enumerate(factorial_rows, 1):
            err_lines.append(f"  {i}: {row}")
        err_lines.append("")
        err_lines.append("Check your categories, excludes, and expected_samples.")
        return effects_fail(err_lines)
    if existing and force:
        lines.append(f"Cleared {len(existing)} existing samples.")
    lines.append(
        f"Generated {n_extras} extras + {n_factorial} factorial = {n_total} samples (expected {expected})"
    )
    # insert into new_store
    for s in samples:
        new_store["samples"].append(
            {
                "id": None,
                "experiment_id": exp_id,
                "position": s["position"],
                "sample_type": s["sample_type"],
                "strain_id": s["strain_id"],
                "dna": s["dna"],
                "label": s["label"],
                "conditions": json.dumps(s["conditions"]),
            }
        )
    # build display table
    cat_display = [c for c in cat_names if c != "strain"]
    extra_keys = set()
    for extra in design.get("extras", []):
        for k in extra:
            if k not in ("strain", "label", "dna") and k not in cat_display:
                extra_keys.add(k)
    extra_keys = sorted(extra_keys)
    lines.append("")
    hdr = f"  {'pos':>5}  {'type':<10} {'strain':<8} {'label':<14}"
    for c in cat_display:
        hdr += f" {c:<8}"
    for k in extra_keys:
        hdr += f" {k:<12}"
    lines.append(hdr)
    for s in samples:
        cond = s["conditions"]
        line = f"  {s['position']:>5}  {s['sample_type']:<10} {(s['strain_id'] or '-'):<8} {s['label']:<14}"
        for c in cat_display:
            v = cond.get(c)
            line += f" {(str(v) if v is not None else '-'):<8}"
        for k in extra_keys:
            v = cond.get(k)
            line += f" {(str(v) if v is not None else '-'):<12}"
        lines.append(line)
    lines.append("")
    lines.append(f"Inserted {n_total} rows into samples table.")
    # Caller is responsible for writing manifest.csv
    result = effects_ok(stdout=lines, store=new_store)
    result["_manifest_samples"] = samples
    result["_manifest_cat_display"] = cat_display
    result["_manifest_extra_keys"] = extra_keys
    return result


def transform_exp_init_plan(store, args, clock):
    """Pure transform: plan experiment creation.

    args: {
        "exp_type": str,
        "shortname": str,
        "title": str | None,
        "counter": int,
    }
    Returns effects with stdout and store with the new experiment row.
    The caller is responsible for folder creation and counter file update.
    """
    exp_type = args["exp_type"]
    shortname = args["shortname"]
    title = args.get("title") or shortname
    counter = args["counter"]
    if exp_type not in VALID_TYPES:
        return effects_fail(
            [
                f"Unknown experiment type: '{exp_type}'",
                f"Valid types: {', '.join(VALID_TYPES)}",
            ]
        )
    ok, reason = validate_name(shortname, "shortname")
    if not ok:
        return effects_fail(reason)
    shortname = shortname.lower()
    if title == shortname:
        pass  # default
    exp_id = f"LM-{counter:04d}"
    today_compact = clock["today"].strftime("%Y%m%d")
    today_sql = clock["today"].strftime("%Y-%m-%d")
    cat_code = CATEGORY_CODES[exp_type]
    folder_name = f"{today_compact}_{exp_id}_{cat_code}_{shortname}"
    new_store = copy.deepcopy(store)
    new_store["experiments"].append(
        {
            "id": exp_id,
            "type": exp_type,
            "title": title,
            "shortname": shortname,
            "date_started": today_sql,
            "date_completed": None,
            "status": "active",
            "folder_name": folder_name,
            "protocol_id": None,
            "notes": None,
        }
    )
    lines = []
    lines.append(f"Created: {exp_id}")
    lines.append(f"Folder:  {folder_name}/")
    lines.append(f"Type:    {exp_type}")
    if title != shortname:
        lines.append(f"Title:   {title}")
    if exp_type == "purification":
        lines.append(f"Next: lw exp update {exp_id} protein <name>")
    elif exp_type in ("loading", "gelshift", "atpase"):
        lines.append(f"Next: lw exp addstrain {exp_id} <strain_id>")
    else:
        lines.append(f"Next: lw exp update {exp_id} notes <text>")
    result = effects_ok(stdout=lines, store=new_store)
    result["_init_plan"] = {
        "exp_id": exp_id,
        "folder_name": folder_name,
        "new_counter": counter + 1,
    }
    return result


# Maps (command, subcommand) -> transform function.
# main() checks this before falling through to the legacy if/elif chain.
TRANSFORM_DISPATCH = {
    ("exp", "complete"): transform_exp_complete,
    ("exp", "show"): transform_exp_show,
    ("exp", "update"): transform_exp_update,
    ("exp", "link"): transform_exp_link,
    ("exp", "addstrain"): transform_exp_addstrain,
    ("exp", "find"): transform_exp_find,
    ("exp", "list"): transform_exp_list,
    ("exp", "delete"): transform_exp_delete_plan,
    ("exp", "init"): transform_exp_init_plan,
    ("exp", "manifest"): transform_exp_manifest_plan,
    ("strain", "add"): transform_strain_add,
    ("strain", "update"): transform_strain_update,
    ("strain", "show"): transform_strain_show,
    ("strain", "list"): transform_strain_list,
    ("stage", "list"): transform_stage_list,
    ("stage", "expect_list"): transform_expect_list,
    ("stage", "expect_cancel"): transform_expect_cancel,
    ("stage", "expect_register"): transform_expect_register,
    ("stage", "auto"): transform_stage_auto_plan,
    ("status", None): transform_status,
}


def bail(code=1):
    """Log command failure and exit. Safe to call even before LAB_ROOT exists."""
    try:
        _log_path = os.path.join(LAB_ROOT, "command_log.txt")
        with open(_log_path, "a") as _f:
            _ts = datetime.datetime.now().isoformat(timespec="seconds")
            _cmd = " ".join(sys.argv[1:])
            _f.write(f"{_ts}\t{_cmd}\tfail\n")
    except Exception:
        pass
    sys.exit(code)


# --- staging filename parser ---
def parse_staging_name(filename):
    """Extract experiment ID and slot from staging filename.
    Parses filenames following the instrument convention:
        lemr LM-NNNN sN [instrument suffix].ext
    Returns (exp_id, slot) or (exp_id, None) or (None, None).
    exp_id is normalized to uppercase, slot is an integer.
    First-match-wins: scans all tokens, returns the first matching
    experiment ID and the first matching slot found.
    """
    stem = os.path.splitext(filename)[0]
    tokens = stem.split()
    exp_id = None
    slot = None
    for token in tokens:
        if re.match(r"^LM-\d{4}$", token, re.IGNORECASE):
            exp_id = token.upper()
        elif re.match(r"^s(\d+)$", token, re.IGNORECASE):
            slot = int(re.match(r"^s(\d+)$", token, re.IGNORECASE).group(1))
    return exp_id, slot


def build_filed_name(exp_id, descriptor, clock, ext):
    """Pure: construct the destination filename for a filed staging file.

    Returns the new filename string, e.g. '20250610_LM-0001_gel-01.tiff'.
    clock is the captured clock dict with a 'today' key.
    ext should include the leading dot and be lowercased by the caller.
    """
    today_compact = clock["today"].strftime("%Y%m%d")
    return f"{today_compact}_{exp_id}_{descriptor}{ext}"


def validate_destination(dest_dir, dest_name):
    """Pure: check whether the destination folder exists and the file is new.

    Returns (True, None) if valid, or (False, reason_string) if not.
    """
    if not os.path.isdir(dest_dir):
        folder = os.path.basename(dest_dir)
        return (False, f"Experiment folder missing: {folder}/")
    dest_path = os.path.join(dest_dir, dest_name)
    if os.path.exists(dest_path):
        return (False, f"Destination already exists: {dest_name}")
    return (True, None)


def file_staged(
    src_path, exp_id, folder_name, descriptor, clock, experiments_dir, db_path
):
    """Move a file from staging to an experiment folder and register it.

    IMPURE -- performs irreversible filesystem move and DB insert.
    Takes explicit parameters instead of using globals.

    Args:
        src_path: absolute path to source file in staging
        exp_id: normalized experiment ID
        folder_name: experiment folder name
        descriptor: file descriptor (lowercased)
        clock: captured clock dict with 'today' key
        experiments_dir: absolute path to experiments/ directory
        db_path: path to SQLite database

    Returns (new_name, rel_path).
    Raises ValueError if experiment folder is missing or destination exists.
    """
    _, ext = os.path.splitext(src_path)
    ext = ext.lower()
    new_name = build_filed_name(exp_id, descriptor, clock, ext)
    dest_dir = os.path.join(experiments_dir, folder_name)
    ok, reason = validate_destination(dest_dir, new_name)
    if not ok:
        raise ValueError(reason)
    dest_path = os.path.join(dest_dir, new_name)
    shutil.move(src_path, dest_path)
    rel_path = os.path.join("experiments", folder_name, new_name)
    _conn = sqlite3.connect(db_path)
    try:
        _conn.execute(
            "INSERT INTO files (experiment_id, filename, file_type) VALUES (?, ?, ?)",
            (exp_id, rel_path, None),
        )
        _conn.commit()
    finally:
        _conn.close()
    return (new_name, rel_path)


def normalize_exp_id(raw_id):
    """Normalize experiment ID to LM-NNNN format.

    Returns the normalized string, or None if the input is invalid.
    Callers are responsible for printing an error and bailing.
    """
    if raw_id.upper().startswith("LM-"):
        return raw_id.upper()
    try:
        return f"LM-{int(raw_id):04d}"
    except ValueError:
        return None


def validate_name(value, label="name"):
    """Validate a name/descriptor: alphanumeric + hyphens, starts alphanumeric.

    Returns (True, None) if valid, or (False, reason_string) if invalid.
    Does not exit, print, or touch the database.
    Caller is responsible for printing the reason and bailing.
    """
    if not value:
        return (False, f"Invalid {label}: cannot be empty.")
    if not all(c.isalnum() or c == "-" for c in value) or not value[0].isalnum():
        return (
            False,
            f"Invalid {label}: '{value}'. Use letters, digits, and hyphens. Must start with a letter or digit.",
        )
    return (True, None)


def lookup_experiment(store, exp_id):
    """Pure lookup: find experiment by ID in the Store.

    Returns the row dict if found, None otherwise.
    Does not exit, print, or touch the database.
    """
    for row in store["experiments"]:
        if row["id"] == exp_id:
            return row
    return None


# ============================================================
# SECTION 3: MAIN
# ============================================================


def main(argv=None):
    global LAB_ROOT, DB_PATH, EXPERIMENTS_DIR, STAGING_DIR, COUNTER_FILE, LAB_NOTES

    if argv is None:
        argv = sys.argv

    # --- Commit 0.2: capture Clock once ---
    clock = {"today": datetime.date.today(), "now": datetime.datetime.now()}

    # --- environment resolution ---
    _lab_root_env = os.environ.get("LW_LAB_ROOT") or os.environ.get("LAB_ROOT")
    if _lab_root_env:
        LAB_ROOT = _lab_root_env
    else:
        _win_user = os.environ.get("LW_WINDOWS_USER") or os.environ.get(
            "MC_WINDOWS_USER"
        )
        if _win_user:
            LAB_ROOT = f"/mnt/c/Users/{_win_user}/Desktop/lab"
        else:
            print("Error: Set LW_LAB_ROOT or LW_WINDOWS_USER environment variable.")
            bail()

    DB_PATH = os.path.join(LAB_ROOT, "lab.db")
    EXPERIMENTS_DIR = os.path.join(LAB_ROOT, "experiments")
    STAGING_DIR = os.path.join(LAB_ROOT, "staging")
    COUNTER_FILE = os.path.join(LAB_ROOT, "counter.txt")
    LAB_NOTES = os.path.join(LAB_ROOT, "lab_notes.md")

    # ============================================================
    # ARGUMENT PARSING
    # ============================================================
    args = argv[1:]
    if not args or args[0] in ("--help", "-h"):
        print("Usage: python lw.py <command> [subcommand] [args]")
        print()
        print("Commands:")
        print("  init              Bootstrap lab directory and database")
        print("  status                        Show system dashboard")
        print("  exp init <type> <shortname>   Create new experiment")
        print("  exp show <id>                 Show experiment details")
        print("  exp update <id> <f> <v>       Update an experiment field")
        print("  exp complete <id>             Mark experiment done")
        print("  exp link <from> <to> <rel>    Link two experiments")
        print("  exp addstrain <exp> <strain>  Associate strain with experiment")
        print("  exp delete <id>               Dry-run delete (--confirm to execute)")
        print(
            "  exp find [query] [--type T] [--status S] [--strain S]  Search experiments"
        )
        print("  exp list [--type T]           List experiments with type summaries")
        print("  exp manifest <id>             Generate sample manifest from design.py")
        print("  strain add <id> <genotype>    Register a new strain")
        print("  strain show <id>              Show strain details")
        print("  strain update <id> <f> <v>    Update a strain field")
        print("  strain list [--no-label]      List all strains with labels")
        print("  stage list                    Show files in staging")
        print("  stage assign <f> <exp> <desc> Rename, move, and register file")
        print("  stage expect <exp> <desc...>  Pre-register expected files")
        print("  stage expect --list [exp]     Show pending expectations")
        print("  stage expect --cancel <exp>   Cancel pending expectations")
        print(
            "  stage auto                    Auto-assign staged files from expectations"
        )
        print()
        print("Reference:")
        print(
            "  Types:     purification loading gelshift atpase genetics computational cloning sequencing"
        )
        print("  Relations: uses_prep_from replicate_of follow_up_to")
        print(
            "  Fields:    protein yield_mg columns storage | protein_prep dna_template | notes title protocol_id"
        )
        print(
            "  Protected: id status folder_name date_started date_completed created_at"
        )
        return 0
    elif args[0] == "--ref":
        print("Field Reference")
        print()
        print("Experiment types:")
        print(
            "  purification loading gelshift atpase genetics computational cloning sequencing"
        )
        print()
        print("Experiment fields (exp update):")
        print("  experiments:    title shortname type protocol_id notes")
        print("  purification:   protein induction_od galactose_hours columns")
        print("                  v5_depletion fractions_pooled yield_mg storage")
        print("  assay:          protein_prep dna_template dna_prep_method conditions")
        print()
        print(
            "Protected (cannot update): id status folder_name date_started date_completed created_at"
        )
        print()
        print("Strain fields (strain update):")
        print(
            "  genotype label parent construction selection_marker verification storage_location notes"
        )
        print()
        print("Relationships (exp link):")
        print("  uses_prep_from  replicate_of  follow_up_to")
        print()
        print("Naming rules:")
        print(
            "  shortname/descriptor: lowercase alphanumeric + hyphens, starts with letter or digit"
        )
        print("  retakes: append -v2, -v3 to descriptor")
        print("  staging instrument name: lemr LM-NNNN sN")
        print()
        print("Behaviors:")
        print(
            "  notes:       always appends with [YYYY-MM-DD] timestamp (exp and strain)"
        )
        print("  exp delete:  dry-run default, --confirm to execute")
        print("  stage auto:  dry-run default, --confirm to execute")
        print("  manifest:    --force to regenerate existing samples")
        print("  filing date: date of filing, not instrument date")
        return 0

    command = args[0]
    subcommand = args[1] if len(args) > 1 else None
    pos_args = args[2:]

    # ============================================================
    # DATABASE / STORE
    # ============================================================
    store = None
    if command != "init":
        if not os.path.exists(DB_PATH):
            print(f"Error: Database not found at {DB_PATH}")
            print("Run 'python lw.py init' first.")
            bail()
        store = load_store(DB_PATH)

    # ============================================================
    # COMMAND DISPATCH
    # ============================================================
    # Transaction model: all DB changes accumulate in a single
    # transaction and commit once at script end (cleanup section).
    # Filesystem operations (makedirs, shutil.move, shutil.rmtree)
    # are immediate and non-reversible. If the process crashes after
    # a filesystem change but before commit, the DB rolls back but
    # the filesystem does not. Commands are ordered to minimize
    # unrecoverable states: prefer DB-then-filesystem where possible.

    # --- Commit 2.2: Transform dispatch (pure commands) ---
    # If (command, subcommand) has a registered transform, call it
    # and feed the result to execute_effects.  Otherwise fall through
    # to the legacy if/elif chain below.
    _dispatched = False
    _transform_key = (command, subcommand)
    # Virtual dispatch keys for stage expect sub-modes
    if command == "stage" and subcommand == "expect":
        if "--list" in pos_args:
            _transform_key = ("stage", "expect_list")
        elif "--cancel" in pos_args:
            _transform_key = ("stage", "expect_cancel")
        else:
            _transform_key = ("stage", "expect_register")
    _transform_fn = TRANSFORM_DISPATCH.get(_transform_key)
    if _transform_fn is not None:
        # Parse args for the transform.  Each transform expects a
        # specific args dict; we build it here from pos_args.
        _targs = None

        # --- single <id> commands ---
        if _transform_key in (("exp", "complete"), ("exp", "show")):
            if len(pos_args) != 1:
                usage = {
                    ("exp", "complete"): "Usage: python lw.py exp complete <id>",
                    (
                        "exp",
                        "show",
                    ): "Usage: python lw.py exp show <id>\n  Accepts 'LM-0001' or just '1'",
                }
                print(usage[_transform_key])
                bail()
            _targs = {"exp_id": pos_args[0]}

        # --- exp update <id> <field> <value...> ---
        elif _transform_key == ("exp", "update"):
            if len(pos_args) < 3:
                print("Usage: python lw.py exp update <id> <field> <value>")
                bail()
            _targs = {
                "exp_id": pos_args[0],
                "field": pos_args[1],
                "value": " ".join(pos_args[2:]),
            }

        # --- exp link <from> <to> <relationship> ---
        elif _transform_key == ("exp", "link"):
            if len(pos_args) < 3:
                print("Usage: python lw.py exp link <from_id> <to_id> <relationship>")
                print("  Relationships: uses_prep_from, replicate_of, follow_up_to")
                bail()
            _targs = {
                "from_id": pos_args[0],
                "to_id": pos_args[1],
                "relationship": pos_args[2],
            }

        # --- exp addstrain <exp_id> <strain_id> [role] ---
        elif _transform_key == ("exp", "addstrain"):
            if len(pos_args) < 2:
                print("Usage: python lw.py exp addstrain <exp_id> <strain_id> [role]")
                bail()
            _targs = {
                "exp_id": pos_args[0],
                "strain_id": pos_args[1],
                "role": pos_args[2] if len(pos_args) > 2 else None,
            }

        # --- strain add <id> <genotype...> ---
        elif _transform_key == ("strain", "add"):
            if len(pos_args) < 2:
                print("Usage: python lw.py strain add <id> <genotype>")
                print(
                    '  Example: python lw.py strain add LY456 "MATa orc4-R267A::KanMX ura3 leu2 trp1 his3"'
                )
                bail()
            _targs = {
                "strain_id": pos_args[0],
                "genotype": " ".join(pos_args[1:]),
            }

        # --- strain update <id> <field> <value...> ---
        elif _transform_key == ("strain", "update"):
            if len(pos_args) < 3:
                print("Usage: python lw.py strain update <id> <field> <value>")
                bail()
            _targs = {
                "strain_id": pos_args[0],
                "field": pos_args[1],
                "value": " ".join(pos_args[2:]),
            }

        # --- strain show <id> ---
        elif _transform_key == ("strain", "show"):
            if len(pos_args) != 1:
                print("Usage: python lw.py strain show <id>")
                bail()
            _targs = {"strain_id": pos_args[0]}

        # --- strain list [--no-label] ---
        elif _transform_key == ("strain", "list"):
            _targs = {"no_label_only": "--no-label" in pos_args}

        # --- exp find [query] [--type T] [--status S] [--strain S] ---
        elif _transform_key == ("exp", "find"):
            _targs = {"pos_args": pos_args}

        # --- exp list [--type T] ---
        elif _transform_key == ("exp", "list"):
            _targs = {"pos_args": pos_args}

        # --- stage list (T7 sandwich: IO then pure) ---
        elif _transform_key == ("stage", "list"):
            _file_list = prepare_stage_list_io(STAGING_DIR)
            _targs = {"file_list": _file_list}

        # --- status (T7 sandwich: IO then pure) ---
        elif _transform_key == ("status", None):
            _fs = prepare_status_io(STAGING_DIR, LAB_NOTES)
            _targs = {"fs": _fs}

        # --- stage expect --list [exp_id] ---
        elif _transform_key == ("stage", "expect_list"):
            filter_args = [a for a in pos_args if a != "--list"]
            filter_exp_id = None
            if filter_args:
                filter_exp_id = normalize_exp_id(filter_args[0])
                if not filter_exp_id:
                    print(f"Invalid experiment ID: '{filter_args[0]}'")
                    bail()
            _targs = {"filter_exp_id": filter_exp_id}

        # --- stage expect --cancel <exp_id> [sN] ---
        elif _transform_key == ("stage", "expect_cancel"):
            cancel_args = [a for a in pos_args if a != "--cancel"]
            if not cancel_args:
                print("Usage: python lw.py stage expect --cancel <exp_id> [sN]")
                bail()
            exp_id = normalize_exp_id(cancel_args[0])
            if not exp_id:
                print(f"Invalid experiment ID: '{cancel_args[0]}'")
                bail()
            slot = None
            if len(cancel_args) > 1:
                m = re.match(r"^s?(\d+)$", cancel_args[1], re.IGNORECASE)
                if not m:
                    print(f"Invalid slot: '{cancel_args[1]}'. Use s1, s2, etc.")
                    bail()
                slot = int(m.group(1))
            _targs = {"exp_id": exp_id, "slot": slot}

        # --- stage expect <exp_id> <descriptor> [descriptor2...] ---
        elif _transform_key == ("stage", "expect_register"):
            if len(pos_args) < 2:
                print(
                    "Usage: python lw.py stage expect <exp_id> <descriptor> [descriptor2...]"
                )
                print(
                    "  Example: python lw.py stage expect LM-0005 gel-coomassie gel-silver"
                )
                bail()
            raw_id = pos_args[0]
            descriptors = pos_args[1:]
            exp_id = normalize_exp_id(raw_id)
            if not exp_id:
                print(f"Invalid experiment ID: '{raw_id}'")
                bail()
            for desc in descriptors:
                _ok, _reason = validate_name(desc, "descriptor")
                if not _ok:
                    print(_reason)
                    bail()
            _targs = {"exp_id": exp_id, "descriptors": [d.lower() for d in descriptors]}

        # --- exp delete <id> [--confirm] [--keep-folder] ---
        elif _transform_key == ("exp", "delete"):
            if len(pos_args) < 1:
                print("Usage: python lw.py exp delete <id> [--confirm] [--keep-folder]")
                print("  Without --confirm: dry run (shows what would be deleted)")
                bail()
            raw_id = pos_args[0]
            _confirm = "--confirm" in pos_args
            _keep_folder = "--keep-folder" in pos_args
            exp_id = normalize_exp_id(raw_id)
            if not exp_id:
                print(f"Invalid experiment ID: '{raw_id}'")
                bail()
            _exp = lookup_experiment(store, exp_id)
            if not _exp:
                print(f"No experiment found with ID {exp_id}")
                bail()
            _folder_path = os.path.join(EXPERIMENTS_DIR, _exp["folder_name"])
            _folder_exists = os.path.exists(_folder_path)
            _folder_contents = os.listdir(_folder_path) if _folder_exists else []
            _targs = {
                "exp_id": exp_id,
                "confirm": _confirm,
                "keep_folder": _keep_folder,
                "folder_exists": _folder_exists,
                "folder_contents": _folder_contents,
            }

        # --- exp init <type> <shortname> [--title <title>] ---
        elif _transform_key == ("exp", "init"):
            _title = None
            _filtered_args = []
            i = 0
            while i < len(pos_args):
                if pos_args[i] == "--title":
                    if i + 1 < len(pos_args):
                        _title = " ".join(pos_args[i + 1 :])
                    else:
                        print("--title requires a value.")
                        bail()
                    break
                else:
                    _filtered_args.append(pos_args[i])
                    i += 1
            if len(_filtered_args) != 2:
                print(
                    "Usage: python lw.py exp init <type> <shortname> [--title <title>]"
                )
                print(f"Types: {', '.join(VALID_TYPES)}")
                bail()
            # read counter (IO)
            try:
                with open(COUNTER_FILE, "r") as f:
                    _raw = f.read().strip()
                _counter = int(_raw)
            except FileNotFoundError:
                print(f"Counter file not found: {COUNTER_FILE}")
                print("Run 'python lw.py init' to initialize.")
                bail()
            except ValueError:
                print(f"Counter file contains invalid value: '{_raw}'")
                print(f"Expected an integer in {COUNTER_FILE}")
                bail()
            _targs = {
                "exp_type": _filtered_args[0],
                "shortname": _filtered_args[1],
                "title": _title,
                "counter": _counter,
            }

        # --- exp manifest <id> [--force] ---
        elif _transform_key == ("exp", "manifest"):
            if len(pos_args) < 1:
                print("Usage: python lw.py exp manifest <id> [--force]")
                bail()
            raw_id = pos_args[0]
            _force = "--force" in pos_args
            exp_id = normalize_exp_id(raw_id)
            if not exp_id:
                print(f"Invalid experiment ID: '{raw_id}'")
                bail()
            _exp = lookup_experiment(store, exp_id)
            if not _exp:
                print(f"No experiment found with ID {exp_id}")
                bail()
            _folder_path = os.path.join(EXPERIMENTS_DIR, _exp["folder_name"])
            _design_path = os.path.join(_folder_path, "design.py")
            if not os.path.exists(_design_path):
                print(f"No design.py found in {_exp['folder_name']}/")
                print(f"Create {_design_path} with a DESIGN dict.")
                bail()
            # load design via exec (IMPURE)
            _design_ns = {}
            with open(_design_path, "r") as f:
                exec(f.read(), _design_ns)
            if "DESIGN" not in _design_ns:
                print("Error: design.py must define a DESIGN dict.")
                bail()
            _targs = {
                "exp_id": exp_id,
                "force": _force,
                "design": _design_ns["DESIGN"],
                "design_source": f"experiments/{_exp['folder_name']}/design.py",
            }

        # --- stage auto [--confirm] ---
        elif _transform_key == ("stage", "auto"):
            _confirm = "--confirm" in pos_args
            if not os.path.exists(STAGING_DIR):
                _file_list = None
            else:
                _staging_files = [
                    f for f in os.listdir(STAGING_DIR) if not f.startswith(".")
                ]
                _file_list = [
                    {"name": f, "src_path": os.path.join(STAGING_DIR, f)}
                    for f in sorted(_staging_files)
                ]
            _targs = {"file_list": _file_list, "confirm": _confirm}

        # --- fallback (should not reach for registered transforms) ---
        else:
            _targs = {"pos_args": pos_args}

        _effects = _transform_fn(store, _targs, clock)
        _rc = execute_effects(_effects, DB_PATH)
        # If the transform mutated the store, skip the legacy
        # conn.commit path -- execute_effects already committed.
        if _effects.get("store") is not None:
            store = None  # prevent double-commit in cleanup

        # --- Post-transform IO for commands that need it ---

        # exp delete --confirm: remove folder
        if _transform_key == ("exp", "delete") and _targs.get("confirm") and _rc == 0:
            _exp = lookup_experiment(
                _targs.get("_orig_store", store) or {}, _targs["exp_id"]
            )
            # Folder removal already printed by transform; do the actual IO
            if _targs.get("folder_exists") and not _targs.get("keep_folder"):
                _fp = os.path.join(
                    EXPERIMENTS_DIR,
                    _targs.get("folder_contents", [""])[0] if False else "",
                )
                # reconstruct folder path from the experiment we looked up earlier
                pass  # handled below

        # We need the original exp for folder path. Let's do it properly:
        if _transform_key == ("exp", "delete") and _targs.get("confirm") and _rc == 0:
            if _targs["folder_exists"] and not _targs["keep_folder"]:
                # _exp was looked up before transform; get folder_name from targs context
                _del_exp_id = _targs["exp_id"]
                # The original store had the experiment; we need its folder_name.
                # Since execute_effects already committed, we need to reconstruct.
                # The folder path was computed from EXPERIMENTS_DIR + exp.folder_name
                # before the transform ran. We still have _folder_path from arg parsing.
                shutil.rmtree(_folder_path)

        # exp init: create folder, update counter
        if _transform_key == ("exp", "init") and _rc == 0:
            _plan = _effects.get("_init_plan", {})
            if _plan:
                _init_folder = os.path.join(EXPERIMENTS_DIR, _plan["folder_name"])
                if os.path.exists(_init_folder):
                    print(f"Error: Folder already exists: {_init_folder}")
                    print("Counter may be out of sync. Check counter.txt.")
                    bail()
                os.makedirs(_init_folder)
                with open(COUNTER_FILE, "w") as f:
                    f.write(str(_plan["new_counter"]))

        # exp manifest: write manifest.csv
        if _transform_key == ("exp", "manifest") and _rc == 0:
            _manifest_samples = _effects.get("_manifest_samples")
            if _manifest_samples is not None:
                _m_exp = lookup_experiment(
                    _effects.get("store") or store, _targs["exp_id"]
                )
                if _m_exp:
                    _csv_path = os.path.join(
                        EXPERIMENTS_DIR, _m_exp["folder_name"], "manifest.csv"
                    )
                    _cat_display = _effects.get("_manifest_cat_display", [])
                    _extra_keys = _effects.get("_manifest_extra_keys", [])
                    _csv_columns = (
                        ["position", "sample_type", "strain_id", "label"]
                        + _cat_display
                        + _extra_keys
                    )
                    with open(_csv_path, "w", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow(_csv_columns)
                        for s in _manifest_samples:
                            cond = s["conditions"]
                            row = [
                                s["position"],
                                s["sample_type"],
                                s["strain_id"] or "",
                                s["label"],
                            ]
                            for c in _cat_display:
                                v = cond.get(c)
                                row.append(str(v) if v is not None else "")
                            for k in _extra_keys:
                                v = cond.get(k)
                                row.append(str(v) if v is not None else "")
                            writer.writerow(row)
                    print(f"Wrote manifest.csv to experiments/{_m_exp['folder_name']}/")

        # stage auto --confirm: execute file moves
        if _transform_key == ("stage", "auto") and _targs.get("confirm") and _rc == 0:
            _plan = _effects.get("_auto_plan")
            if _plan:
                filed_count = 0
                filed_display = []
                skip_display = list(_plan["skipped"])
                filed_slots = set()
                error_slots = set()
                for mf in _plan["matched_files"]:
                    fname = mf["fname"]
                    src_path = os.path.join(STAGING_DIR, fname)
                    exp_id = mf["exp_id"]
                    slot = mf["slot"]
                    desc = mf["descriptor"]
                    folder_name = mf["folder_name"]
                    dest_name = mf["dest_name"]
                    key = (exp_id, slot)
                    try:
                        new_name, _ = file_staged(
                            src_path,
                            exp_id,
                            folder_name,
                            desc,
                            clock,
                            EXPERIMENTS_DIR,
                            DB_PATH,
                        )
                        filed_count += 1
                        filed_display.append((fname, new_name))
                        filed_slots.add(key)
                    except ValueError as e:
                        skip_display.append((fname, str(e)))
                        error_slots.add(key)
                # If any slots had errors, revert their filed status
                if error_slots:
                    _revert_conn = sqlite3.connect(DB_PATH)
                    try:
                        for key in error_slots:
                            eid, sl = key
                            _revert_conn.execute(
                                "UPDATE stage_expectations SET status = 'pending', filed_at = NULL "
                                "WHERE experiment_id = ? AND slot = ?",
                                (eid, sl),
                            )
                        _revert_conn.commit()
                    finally:
                        _revert_conn.close()
                for fname in _plan["unmatched"]:
                    skip_display.append((fname, "no match"))
                if filed_display:
                    print(f"Filed ({len(filed_display)}):")
                    for fname, new_name in filed_display:
                        print(f"  {fname}")
                        print(f"     -> {new_name}")
                if skip_display:
                    print(f"Skipped ({len(skip_display)}):")
                    for fname, reason in skip_display:
                        print(f"  {fname} ({reason})")
                print()
                print(f"{filed_count} filed, {len(skip_display)} skipped.")

        if _rc != 0:
            bail(_rc)
        _dispatched = True

    # ============================================================
    # LEGACY COMMAND DISPATCH
    # ============================================================
    # Only 'init' and 'stage assign' remain here.
    # 'init' cannot use the Store (DB does not exist yet).
    # 'stage assign' is a thin IO wrapper around file_staged.
    # All other commands are transform-dispatched above.

    # --- lw init ---
    if _dispatched:
        pass
    elif command == "init":
        if os.path.exists(DB_PATH):
            print(f"Database already exists at {DB_PATH}")
            print("Init aborted to protect existing data.")
            bail()
        for d in DIRECTORIES:
            os.makedirs(os.path.join(LAB_ROOT, d), exist_ok=True)
        print(f"Created directory structure at {LAB_ROOT}/")
        conn = sqlite3.connect(DB_PATH)
        conn.executescript(SCHEMA_SQL)
        conn.commit()
        conn.close()
        conn = None
        print(f"Created database at {DB_PATH}")
        print(
            f"Tables: {', '.join(['experiments', 'strains', 'experiment_strains', 'experiment_links', 'files', 'purification_details', 'assay_details', 'figure_panels', 'samples', 'stage_expectations'])}"
        )
        if not os.path.exists(COUNTER_FILE):
            with open(COUNTER_FILE, "w") as f:
                f.write("1")
            print(f"Created {COUNTER_FILE} (starting at 1)")
        else:
            print(f"Counter file already exists, skipping")
        if not os.path.exists(LAB_NOTES):
            today = clock["today"].strftime("%Y-%m-%d")
            with open(LAB_NOTES, "w") as f:
                f.write(f"## {today}\n\n")
            print(f"Created {LAB_NOTES}")
        else:
            print(f"Lab notes file already exists, skipping")
        print()
        print("Lab initialized. Ready for experiments.")

    elif command == "exp":
        if not subcommand:
            print("Usage: python lw.py exp <subcommand>")
            print("  init <type> <shortname>   Create new experiment")
            print("  show <id>                 Show experiment details")
            return 0

    elif command == "strain":
        if not subcommand:
            print("Usage: python lw.py strain <subcommand>")
            print("  add <id> <genotype>   Register a new strain")
            print("  show <id>             Show strain details and experiments")
            print("  list [--no-label]     List all strains (--no-label: missing only)")
            print("  update <id> <f> <v>   Update a strain field")
            return 0

    elif command == "stage":
        if not subcommand:
            print("Usage: python lw.py stage <subcommand>")
            print("  list                          Show files in staging/")
            print("  assign <file> <exp_id> <desc> Rename, move, and register")
            print("  expect <exp_id> <desc...>     Pre-register expected files")
            print("  expect --list [exp_id]        Show pending expectations")
            print("  expect --cancel <exp_id> [sN] Cancel pending expectations")
            print("  auto [--confirm]              Auto-assign from expectations")
            return 0
        # --- lw stage assign (kept in legacy: thin IO wrapper) ---
        if subcommand == "assign":
            if len(pos_args) < 3:
                print(
                    "Usage: python lw.py stage assign <filename> <exp_id> <descriptor>"
                )
                print(
                    "  Example: python lw.py stage assign Image_001.tiff LM-0001 gel-01"
                )
                bail()
            src_name, raw_id, descriptor = pos_args[0], pos_args[1], pos_args[2]
            _ok, _reason = validate_name(descriptor, "descriptor")
            if not _ok:
                print(_reason)
                bail()
            descriptor = descriptor.lower()
            src_path = os.path.join(STAGING_DIR, src_name)
            if not os.path.exists(src_path):
                print(f"File not found in staging: {src_name}")
                files = [f for f in os.listdir(STAGING_DIR) if not f.startswith(".")]
                if files:
                    print("Available files:")
                    for f in sorted(files):
                        print(f"  {f}")
                bail()
            exp_id = normalize_exp_id(raw_id)
            if not exp_id:
                print(f"Invalid experiment ID: '{raw_id}'")
                bail()
            exp = lookup_experiment(store, exp_id)
            if not exp:
                print(f"No experiment found with ID {exp_id}")
                bail()
            folder_name = exp["folder_name"]
            try:
                new_name, rel_path = file_staged(
                    src_path,
                    exp_id,
                    folder_name,
                    descriptor,
                    clock,
                    EXPERIMENTS_DIR,
                    DB_PATH,
                )
            except ValueError as e:
                print(str(e))
                print("Use a different descriptor.")
                bail()
            print(f"Experiment: {exp_id} ({exp['shortname']}, {exp['type']})")
            print(f"Filed: {src_name}")
            print(f"    -> {rel_path}")
        else:
            print(f"Unknown subcommand: stage {subcommand}")
            print("Run 'python lw.py stage' for usage.")
            bail()

    # --- unknown command ---
    else:
        print(f"Unknown command: {command}")
        print("Run 'python lw.py' for usage.")
        bail()

    # ============================================================
    # CLEANUP
    # ============================================================
    # log successful command
    try:
        _log_path = os.path.join(LAB_ROOT, "command_log.txt")
        with open(_log_path, "a") as _log_f:
            _ts = datetime.datetime.now().isoformat(timespec="seconds")
            _cmd = " ".join(argv[1:])
            _log_f.write(f"{_ts}\t{_cmd}\tok\n")
    except Exception:
        pass


if __name__ == "__main__":
    sys.exit(main())
