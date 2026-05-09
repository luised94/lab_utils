#!/usr/bin/env python3
# ============================================================
# SECTION 0: IMPORTS
# ============================================================
import sqlite3
import sys
import os
import datetime
import shutil
import itertools
import csv
import json
# ============================================================
# SECTION 1: CONFIGURATION
# ============================================================
# --- paths ---
# LAB_ROOT from env var, otherwise construct from MC_WINDOWS_USER
_lab_root_env = os.environ.get("LAB_ROOT")
if _lab_root_env:
    LAB_ROOT = _lab_root_env
else:
    _win_user = os.environ.get("MC_WINDOWS_USER")
    if _win_user:
        LAB_ROOT = f"/mnt/c/Users/{_win_user}/Desktop/lab"
    else:
        print("Error: Set LAB_ROOT or MC_WINDOWS_USER environment variable.")
        sys.exit(1)
DB_PATH = os.path.join(LAB_ROOT, "lab.db")
EXPERIMENTS_DIR = os.path.join(LAB_ROOT, "experiments")
STAGING_DIR = os.path.join(LAB_ROOT, "staging")
COUNTER_FILE = os.path.join(LAB_ROOT, "counter.txt")
LAB_NOTES = os.path.join(LAB_ROOT, "lab_notes.md")
DIRECTORIES = [
    "experiments",
    "protocols",
    "publications",
    "resources",
    "staging",
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
VALID_TYPES = list(CATEGORY_CODES.keys())
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
"""
# ============================================================
# SECTION 2: ARGUMENT PARSING
# ============================================================
args = sys.argv[1:]
if not args:
    print("Usage: python lw.py <command> [subcommand] [args]")
    print()
    print("Commands:")
    print("  init              Bootstrap lab directory and database")
    print("  exp init <type> <shortname>   Create new experiment")
    print("  exp show <id>                 Show experiment details")
    print("  exp update <id> <f> <v>       Update an experiment field")
    print("  exp complete <id>             Mark experiment done")
    print("  exp link <from> <to> <rel>    Link two experiments")
    print("  exp addstrain <exp> <strain>  Associate strain with experiment")
    print("  exp delete <id>               Dry-run delete (--confirm to execute)")
    print("  exp find [query] [--type T] [--status S] [--strain S]  Search experiments")
    print("  exp manifest <id>             Generate sample manifest from design.py")
    print("  strain add <id> <genotype>    Register a new strain")
    print("  strain show <id>              Show strain details")
    print("  strain update <id> <f> <v>    Update a strain field")
    print("  stage list                    Show files in staging")
    print("  stage assign <f> <exp> <desc> Rename, move, and register file")
    sys.exit(0)
command = args[0]
subcommand = args[1] if len(args) > 1 else None
pos_args = args[2:]
# ============================================================
# SECTION 3: DATABASE CONNECTION
# ============================================================
conn = None
cursor = None
if command != "init":
    if not os.path.exists(DB_PATH):
        print(f"Error: Database not found at {DB_PATH}")
        print("Run 'python lw.py init' first.")
        sys.exit(1)
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    cursor = conn.cursor()
# ============================================================
# SECTION 4: COMMAND DISPATCH
# ============================================================
# --- lw init ---
if command == "init":
    if os.path.exists(DB_PATH):
        print(f"Database already exists at {DB_PATH}")
        print("Init aborted to protect existing data.")
        sys.exit(1)
    # create directories
    for d in DIRECTORIES:
        os.makedirs(os.path.join(LAB_ROOT, d), exist_ok=True)
    print(f"Created directory structure at {LAB_ROOT}/")
    # create database
    conn = sqlite3.connect(DB_PATH)
    conn.executescript(SCHEMA_SQL)
    conn.commit()
    conn.close()
    conn = None
    print(f"Created database at {DB_PATH}")
    print(f"Tables: {', '.join(['experiments', 'strains', 'experiment_strains', 'experiment_links', 'files', 'purification_details', 'assay_details', 'figure_panels', 'samples'])}")
    # seed counter
    if not os.path.exists(COUNTER_FILE):
        with open(COUNTER_FILE, "w") as f:
            f.write("1")
        print(f"Created {COUNTER_FILE} (starting at 1)")
    else:
        print(f"Counter file already exists, skipping")
    # seed lab notes
    if not os.path.exists(LAB_NOTES):
        today = datetime.date.today().strftime("%Y-%m-%d")
        with open(LAB_NOTES, "w") as f:
            f.write(f"## {today}\n\n")
        print(f"Created {LAB_NOTES}")
    else:
        print(f"Lab notes file already exists, skipping")
    print()
    print("Lab initialized. Ready for experiments.")
# --- lw exp ---
elif command == "exp":
    if not subcommand:
        print("Usage: python lw.py exp <subcommand>")
        print("  init <type> <shortname>   Create new experiment")
        print("  show <id>                 Show experiment details")
        sys.exit(0)
    # --- lw exp init <type> <shortname> ---
    if subcommand == "init":
        if len(pos_args) != 2:
            print("Usage: python lw.py exp init <type> <shortname>")
            print(f"Types: {', '.join(VALID_TYPES)}")
            sys.exit(1)
        exp_type, shortname = pos_args[0], pos_args[1]
        # validate type
        if exp_type not in VALID_TYPES:
            print(f"Unknown experiment type: '{exp_type}'")
            print(f"Valid types: {', '.join(VALID_TYPES)}")
            sys.exit(1)
        # validate shortname: lowercase alphanumeric and hyphens only
        if not all(c.isalnum() or c == "-" for c in shortname) or not shortname[0].isalnum():
            print("Invalid shortname. Use lowercase letters, digits, and hyphens only.")
            print("Must start with a letter or digit.")
            sys.exit(1)
        shortname = shortname.lower()
        # read and increment counter
        with open(COUNTER_FILE, "r") as f:
            counter = int(f.read().strip())
        exp_id = f"LM-{counter:04d}"
        # build folder name
        today_compact = datetime.date.today().strftime("%Y%m%d")
        today_sql = datetime.date.today().strftime("%Y-%m-%d")
        cat_code = CATEGORY_CODES[exp_type]
        folder_name = f"{today_compact}_{exp_id}_{cat_code}_{shortname}"
        folder_path = os.path.join(EXPERIMENTS_DIR, folder_name)
        # guard: folder shouldn't exist
        if os.path.exists(folder_path):
            print(f"Error: Folder already exists: {folder_path}")
            print("Counter may be out of sync. Check counter.txt.")
            sys.exit(1)
        # create folder
        os.makedirs(folder_path)
        # insert record
        cursor.execute(
            "INSERT INTO experiments (id, type, title, shortname, date_started, status, folder_name) VALUES (?, ?, ?, ?, ?, 'active', ?)",
            (exp_id, exp_type, shortname, shortname, today_sql, folder_name),
        )
        # increment counter
        with open(COUNTER_FILE, "w") as f:
            f.write(str(counter + 1))
        print(f"Created: {exp_id}")
        print(f"Folder:  {folder_name}/")
        print(f"Type:    {exp_type}")
    # --- lw exp update <id> <field> <value> ---
    elif subcommand == "update":
        if len(pos_args) < 3:
            print("Usage: python lw.py exp update <id> <field> <value>")
            sys.exit(1)
        raw_id, field, value = pos_args[0], pos_args[1], pos_args[2]
        # normalize ID
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        # verify experiment exists
        cursor.execute("SELECT id, type FROM experiments WHERE id = ?", (exp_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        exp_type = row[1]
        # find which table owns this field
        target_table = None
        for table in ["experiments", "purification_details", "assay_details"]:
            cursor.execute(f"PRAGMA table_info({table})")
            columns = [r[1] for r in cursor.fetchall()]
            if field in columns:
                target_table = table
                break
        if not target_table:
            # collect all valid fields across tables
            all_fields = []
            for table in ["experiments", "purification_details", "assay_details"]:
                cursor.execute(f"PRAGMA table_info({table})")
                all_fields.extend(r[1] for r in cursor.fetchall() if r[1] not in ("id", "experiment_id", "created_at"))
            print(f"Unknown field: '{field}'")
            print(f"Valid fields: {', '.join(sorted(set(all_fields)))}")
            sys.exit(1)
        # for type-specific tables, ensure row exists
        if target_table in ("purification_details", "assay_details"):
            cursor.execute(f"SELECT experiment_id FROM {target_table} WHERE experiment_id = ?", (exp_id,))
            if not cursor.fetchone():
                # insert skeleton row
                if target_table == "purification_details":
                    if field == "protein":
                        cursor.execute(
                            "INSERT INTO purification_details (experiment_id, protein) VALUES (?, ?)",
                            (exp_id, value),
                        )
                        print(f"Created purification record for {exp_id}")
                        print(f"  protein: -  {value}")
                    else:
                        print(f"No purification record for {exp_id}.")
                        print(f"Set 'protein' first: python lw.py exp update {exp_id} protein <name>")
                        sys.exit(1)
                    # already inserted, skip the UPDATE below
                    target_table = None
                elif target_table == "assay_details":
                    cursor.execute(
                        "INSERT INTO assay_details (experiment_id) VALUES (?)",
                        (exp_id,),
                    )
                    print(f"Created assay record for {exp_id}")
        if target_table:
            # get old value
            id_col = "id" if target_table == "experiments" else "experiment_id"
            cursor.execute(f"SELECT {field} FROM {target_table} WHERE {id_col} = ?", (exp_id,))
            old_row = cursor.fetchone()
            old_val = old_row[0] if old_row else None
            # notes fields: always prefix with timestamp, append if existing
            if field == "notes":
                today = datetime.date.today().strftime("%Y-%m-%d")
                if old_val:
                    value = old_val + f"\n[{today}] " + value
                else:
                    value = f"[{today}] " + value
            cursor.execute(
                f"UPDATE {target_table} SET {field} = ? WHERE {id_col} = ?",
                (value, exp_id),
            )
            print(f"Updated {exp_id}")
            if field == "notes":
                print(f"  {field}: appended entry")
            else:
                print(f"  {field}: {old_val if old_val else '-'}  {value}")
    # --- lw exp complete <id> ---
    elif subcommand == "complete":
        if len(pos_args) != 1:
            print("Usage: python lw.py exp complete <id>")
            sys.exit(1)
        raw_id = pos_args[0]
        # normalize ID
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        cursor.execute("SELECT status FROM experiments WHERE id = ?", (exp_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        if row[0] == "complete":
            print(f"{exp_id} is already marked complete.")
            sys.exit(0)
        today = datetime.date.today().strftime("%Y-%m-%d")
        cursor.execute(
            "UPDATE experiments SET status = 'complete', date_completed = ? WHERE id = ?",
            (today, exp_id),
        )
        print(f"Completed: {exp_id}")
        print(f"  Date: {today}")
    # --- lw exp link <from_id> <to_id> <relationship> ---
    elif subcommand == "link":
        if len(pos_args) < 3:
            print("Usage: python lw.py exp link <from_id> <to_id> <relationship>")
            print("  Relationships: uses_prep_from, replicate_of, follow_up_to")
            sys.exit(1)
        raw_from, raw_to, relationship = pos_args[0], pos_args[1], pos_args[2]
        # normalize IDs
        for label, raw in [("from", raw_from), ("to", raw_to)]:
            if raw.upper().startswith("LM-"):
                pass
            else:
                try:
                    int(raw)
                except ValueError:
                    print(f"Invalid {label} experiment ID: '{raw}'")
                    sys.exit(1)
        from_id = raw_from.upper() if raw_from.upper().startswith("LM-") else f"LM-{int(raw_from):04d}"
        to_id = raw_to.upper() if raw_to.upper().startswith("LM-") else f"LM-{int(raw_to):04d}"
        # validate both experiments exist
        for eid in (from_id, to_id):
            cursor.execute("SELECT id FROM experiments WHERE id = ?", (eid,))
            if not cursor.fetchone():
                print(f"No experiment found with ID {eid}")
                sys.exit(1)
        # warn on unknown relationship (don't block)
        known_rels = ["uses_prep_from", "replicate_of", "follow_up_to"]
        if relationship not in known_rels:
            print(f"  Note: '{relationship}' is not a standard relationship ({', '.join(known_rels)})")
        # check for duplicate
        cursor.execute(
            "SELECT 1 FROM experiment_links WHERE from_experiment = ? AND to_experiment = ?",
            (from_id, to_id),
        )
        if cursor.fetchone():
            print(f"Link already exists: {from_id}  {to_id}")
            sys.exit(1)
        cursor.execute(
            "INSERT INTO experiment_links (from_experiment, to_experiment, relationship) VALUES (?, ?, ?)",
            (from_id, to_id, relationship),
        )
        # auto-populate assay_details.protein_prep if linking to a purification
        if relationship == "uses_prep_from":
            cursor.execute(
                "SELECT experiment_id FROM assay_details WHERE experiment_id = ?",
                (from_id,),
            )
            if cursor.fetchone():
                cursor.execute(
                    "UPDATE assay_details SET protein_prep = ? WHERE experiment_id = ? AND protein_prep IS NULL",
                    (to_id, from_id),
                )
            else:
                cursor.execute(
                    "INSERT INTO assay_details (experiment_id, protein_prep) VALUES (?, ?)",
                    (from_id, to_id),
                )
        print(f"Linked: {from_id}  {to_id} ({relationship})")
    # --- lw exp addstrain <exp_id> <strain_id> [role] ---
    elif subcommand == "addstrain":
        if len(pos_args) < 2:
            print("Usage: python lw.py exp addstrain <exp_id> <strain_id> [role]")
            sys.exit(1)
        raw_id, strain_id = pos_args[0], pos_args[1]
        role = pos_args[2] if len(pos_args) > 2 else None
        # normalize experiment ID
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        # validate experiment exists
        cursor.execute("SELECT id FROM experiments WHERE id = ?", (exp_id,))
        if not cursor.fetchone():
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        # validate strain exists
        cursor.execute("SELECT id FROM strains WHERE id = ?", (strain_id,))
        if not cursor.fetchone():
            print(f"No strain found with ID {strain_id}")
            print("Register it first: python lw.py strain add <id> <genotype>")
            sys.exit(1)
        # check for duplicate
        cursor.execute(
            "SELECT 1 FROM experiment_strains WHERE experiment_id = ? AND strain_id = ?",
            (exp_id, strain_id),
        )
        if cursor.fetchone():
            print(f"Strain {strain_id} already linked to {exp_id}")
            sys.exit(1)
        cursor.execute(
            "INSERT INTO experiment_strains (experiment_id, strain_id, role) VALUES (?, ?, ?)",
            (exp_id, strain_id, role),
        )
        role_str = f" as {role}" if role else ""
        print(f"Added: {strain_id}  {exp_id}{role_str}")
    # --- lw exp show <id> ---
    elif subcommand == "show":
        if len(pos_args) != 1:
            print("Usage: python lw.py exp show <id>")
            print("  Accepts 'LM-0001' or just '1'")
            sys.exit(1)
        raw_id = pos_args[0]
        # normalize: accept "1", "0001", "LM-0001"
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        # query experiment
        cursor.execute("SELECT * FROM experiments WHERE id = ?", (exp_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        col_names = [d[0] for d in cursor.description]
        exp = dict(zip(col_names, row))
        # display
        fields = [
            ("ID", exp["id"]),
            ("Type", exp["type"]),
            ("Title", exp["title"]),
            ("Shortname", exp["shortname"]),
            ("Date started", exp["date_started"]),
            ("Date completed", exp["date_completed"]),
            ("Status", exp["status"]),
            ("Folder", exp["folder_name"] + "/"),
            ("Protocol", exp["protocol_id"]),
            ("Notes", exp["notes"]),
        ]
        label_width = max(len(f[0]) for f in fields)
        for label, val in fields:
            print(f"  {label + ':':<{label_width + 2}} {val if val else '-'}")
        # linked strains
        cursor.execute(
            "SELECT strain_id, role FROM experiment_strains WHERE experiment_id = ?",
            (exp_id,),
        )
        strains = cursor.fetchall()
        if strains:
            print()
            print("  Strains:")
            for sid, role in strains:
                role_str = f" ({role})" if role else ""
                print(f"    {sid}{role_str}")
        # linked experiments
        cursor.execute(
            "SELECT to_experiment, relationship FROM experiment_links WHERE from_experiment = ?",
            (exp_id,),
        )
        links_out = cursor.fetchall()
        cursor.execute(
            "SELECT from_experiment, relationship FROM experiment_links WHERE to_experiment = ?",
            (exp_id,),
        )
        links_in = cursor.fetchall()
        if links_out or links_in:
            print()
            print("  Links:")
            for target, rel in links_out:
                print(f"     {target} ({rel})")
            for source, rel in links_in:
                print(f"     {source} ({rel})")
        # files
        cursor.execute(
            "SELECT filename, file_type, description FROM files WHERE experiment_id = ?",
            (exp_id,),
        )
        exp_files = cursor.fetchall()
        if exp_files:
            print()
            print("  Files:")
            for fname, ftype, desc in exp_files:
                parts = [fname]
                if ftype:
                    parts.append(f"[{ftype}]")
                if desc:
                    parts.append(desc)
                print(f"    {' '.join(parts)}")
        # purification details
        cursor.execute(
            "SELECT * FROM purification_details WHERE experiment_id = ?",
            (exp_id,),
        )
        pur_row = cursor.fetchone()
        if pur_row:
            pur_cols = [d[0] for d in cursor.description]
            pur = dict(zip(pur_cols, pur_row))
            print()
            print("  Purification details:")
            for k, v in pur.items():
                if k == "experiment_id" or v is None:
                    continue
                print(f"    {k}: {v}")
        # assay details
        cursor.execute(
            "SELECT * FROM assay_details WHERE experiment_id = ?",
            (exp_id,),
        )
        assay_row = cursor.fetchone()
        if assay_row:
            assay_cols = [d[0] for d in cursor.description]
            assay = dict(zip(assay_cols, assay_row))
            print()
            print("  Assay details:")
            for k, v in assay.items():
                if k == "experiment_id" or v is None:
                    continue
                print(f"    {k}: {v}")
        # samples
        cursor.execute(
            "SELECT position, sample_type, strain_id, label, dna, conditions "
            "FROM samples WHERE experiment_id = ? ORDER BY position",
            (exp_id,),
        )
        sample_rows = cursor.fetchall()
        if sample_rows:
            print()
            print(f"  Samples ({len(sample_rows)}):")
            for pos, stype, sid, lbl, dna, cond in sample_rows:
                sid_str = sid if sid else "-"
                dna_str = dna if dna else ""
                extra = f"  dna={dna_str}" if dna_str else ""
                print(f"    {pos:>3}  {stype:<10} {sid_str:<8} {lbl}{extra}")
    # --- lw exp delete <id> [--confirm] [--keep-folder] ---
    elif subcommand == "delete":
        if len(pos_args) < 1:
            print("Usage: python lw.py exp delete <id> [--confirm] [--keep-folder]")
            print("  Without --confirm: dry run (shows what would be deleted)")
            sys.exit(1)
        raw_id = pos_args[0]
        confirm = "--confirm" in pos_args
        keep_folder = "--keep-folder" in pos_args
        # normalize ID
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        # verify experiment exists
        cursor.execute("SELECT * FROM experiments WHERE id = ?", (exp_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        col_names = [d[0] for d in cursor.description]
        exp = dict(zip(col_names, row))
        # gather cascade targets
        cascade = []
        cursor.execute(
            "SELECT strain_id, role FROM experiment_strains WHERE experiment_id = ?",
            (exp_id,),
        )
        es_rows = cursor.fetchall()
        if es_rows:
            cascade.append(("experiment_strains", len(es_rows),
                [f"  {sid}" + (f" ({role})" if role else "") for sid, role in es_rows]))
        cursor.execute(
            "SELECT from_experiment, to_experiment, relationship FROM experiment_links "
            "WHERE from_experiment = ? OR to_experiment = ?",
            (exp_id, exp_id),
        )
        el_rows = cursor.fetchall()
        if el_rows:
            cascade.append(("experiment_links", len(el_rows),
                [f"  {fr} -> {to} ({rel})" for fr, to, rel in el_rows]))
        cursor.execute(
            "SELECT filename FROM files WHERE experiment_id = ?",
            (exp_id,),
        )
        f_rows = cursor.fetchall()
        if f_rows:
            cascade.append(("files", len(f_rows),
                [f"  {fn}" for (fn,) in f_rows]))
        cursor.execute(
            "SELECT experiment_id FROM purification_details WHERE experiment_id = ?",
            (exp_id,),
        )
        if cursor.fetchone():
            cascade.append(("purification_details", 1, []))
        cursor.execute(
            "SELECT experiment_id FROM assay_details WHERE experiment_id = ?",
            (exp_id,),
        )
        if cursor.fetchone():
            cascade.append(("assay_details", 1, []))
        cursor.execute(
            "SELECT id, figure, panel FROM figure_panels WHERE experiment_id = ?",
            (exp_id,),
        )
        fp_rows = cursor.fetchall()
        if fp_rows:
            cascade.append(("figure_panels", len(fp_rows),
                [f"  {fig} panel {pan}" for _, fig, pan in fp_rows]))
        cursor.execute(
            "SELECT COUNT(*) FROM samples WHERE experiment_id = ?",
            (exp_id,),
        )
        s_count = cursor.fetchone()[0]
        if s_count:
            cascade.append(("samples", s_count, []))
        # folder
        folder_path = os.path.join(EXPERIMENTS_DIR, exp["folder_name"])
        folder_exists = os.path.exists(folder_path)
        folder_contents = os.listdir(folder_path) if folder_exists else []
        # print summary
        mode = "DELETE" if confirm else "DRY RUN"
        print(f"[{mode}] {exp_id}: {exp['type']} / {exp['shortname']} ({exp['status']})")
        print(f"  Folder: {exp['folder_name']}/ ({'exists' if folder_exists else 'missing'}"
              + (f", {len(folder_contents)} files" if folder_contents else "") + ")")
        if cascade:
            print()
            print("  Related records:")
            for table, count, details in cascade:
                print(f"    {table}: {count} row{'s' if count != 1 else ''}")
                for d in details:
                    print(f"      {d}")
        if not confirm:
            print()
            print(f"  This is a dry run. To delete, run:")
            print(f"    python lw.py exp delete {exp_id} --confirm")
            if folder_exists:
                print(f"    python lw.py exp delete {exp_id} --confirm --keep-folder  (keep files)")
            sys.exit(0)
        # --- actual deletion ---
        cursor.execute("DELETE FROM samples WHERE experiment_id = ?", (exp_id,))
        cursor.execute("DELETE FROM figure_panels WHERE experiment_id = ?", (exp_id,))
        cursor.execute("DELETE FROM files WHERE experiment_id = ?", (exp_id,))
        cursor.execute("DELETE FROM assay_details WHERE experiment_id = ?", (exp_id,))
        cursor.execute("DELETE FROM purification_details WHERE experiment_id = ?", (exp_id,))
        cursor.execute("DELETE FROM experiment_strains WHERE experiment_id = ?", (exp_id,))
        cursor.execute(
            "DELETE FROM experiment_links WHERE from_experiment = ? OR to_experiment = ?",
            (exp_id, exp_id),
        )
        cursor.execute("DELETE FROM experiments WHERE id = ?", (exp_id,))
        total_rows = sum(c for _, c, _ in cascade) + 1  # +1 for experiments row
        print()
        print(f"  Deleted {total_rows} database row{'s' if total_rows != 1 else ''}")
        if folder_exists and not keep_folder:
            shutil.rmtree(folder_path)
            print(f"  Removed folder: {exp['folder_name']}/")
        elif folder_exists and keep_folder:
            print(f"  Folder kept: {exp['folder_name']}/")
        print()
        print(f"Deleted: {exp_id}")
    # --- lw exp find [query] [--type T] [--status S] [--strain S] ---
    elif subcommand == "find":
        # parse flags and free-text query
        query = None
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
                query = pos_args[i]
                i += 1
            else:
                print(f"Unknown flag: {pos_args[i]}")
                sys.exit(1)
        # validate type if given
        if f_type and f_type not in VALID_TYPES:
            print(f"Unknown type: '{f_type}'")
            print(f"Valid types: {', '.join(VALID_TYPES)}")
            sys.exit(1)
        # build query
        sql = "SELECT DISTINCT e.id, e.type, e.shortname, e.status, e.date_started FROM experiments e"
        conditions = []
        params = []
        if f_strain:
            sql += " JOIN experiment_strains es ON e.id = es.experiment_id"
            conditions.append("es.strain_id = ?")
            params.append(f_strain)
        if query:
            like = f"%{query}%"
            conditions.append(
                "(e.id LIKE ? OR e.shortname LIKE ? OR e.title LIKE ? OR e.notes LIKE ?)"
            )
            params.extend([like, like, like, like])
        if f_type:
            conditions.append("e.type = ?")
            params.append(f_type)
        if f_status:
            conditions.append("e.status = ?")
            params.append(f_status)
        if conditions:
            sql += " WHERE " + " AND ".join(conditions)
        sql += " ORDER BY e.date_started DESC, e.id DESC"
        cursor.execute(sql, params)
        rows = cursor.fetchall()
        if not rows:
            print("No experiments found.")
            sys.exit(0)
        # display
        print(f"{'ID':<10} {'Type':<14} {'Shortname':<24} {'Status':<10} {'Started'}")
        print("-" * 72)
        for eid, etype, sname, status, started in rows:
            print(f"{eid:<10} {etype:<14} {sname:<24} {status:<10} {started}")
        print()
        print(f"{len(rows)} experiment{'s' if len(rows) != 1 else ''}")
    # --- lw exp manifest <id> [--force] ---
    elif subcommand == "manifest":
        if len(pos_args) < 1:
            print("Usage: python lw.py exp manifest <id> [--force]")
            sys.exit(1)
        raw_id = pos_args[0]
        force = "--force" in pos_args
        # normalize ID
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        # verify experiment exists
        cursor.execute("SELECT folder_name FROM experiments WHERE id = ?", (exp_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        folder_name = row[0]
        folder_path = os.path.join(EXPERIMENTS_DIR, folder_name)
        # locate design.py
        design_path = os.path.join(folder_path, "design.py")
        if not os.path.exists(design_path):
            print(f"No design.py found in {folder_name}/")
            print(f"Create {design_path} with a DESIGN dict.")
            sys.exit(1)
        print(f"Loading design from experiments/{folder_name}/design.py")
        # load config via exec
        design_ns = {}
        with open(design_path, "r") as f:
            exec(f.read(), design_ns)
        if "DESIGN" not in design_ns:
            print("Error: design.py must define a DESIGN dict.")
            sys.exit(1)
        design = design_ns["DESIGN"]
        # validate required keys
        for key in ("experiment_id", "expected_samples", "categories", "sort_order"):
            if key not in design:
                print(f"Error: DESIGN missing required key '{key}'")
                sys.exit(1)
        if design["experiment_id"] != exp_id:
            print(f"Error: DESIGN experiment_id '{design['experiment_id']}' does not match {exp_id}")
            sys.exit(1)
        # validate sort_order matches category keys
        cat_keys = set(design["categories"].keys())
        sort_keys = set(design["sort_order"])
        if sort_keys != cat_keys:
            missing = cat_keys - sort_keys
            extra = sort_keys - cat_keys
            if missing:
                print(f"Error: sort_order missing categories: {', '.join(missing)}")
            if extra:
                print(f"Error: sort_order references unknown categories: {', '.join(extra)}")
            sys.exit(1)
        # check for existing samples
        cursor.execute("SELECT COUNT(*) FROM samples WHERE experiment_id = ?", (exp_id,))
        existing = cursor.fetchone()[0]
        if existing > 0:
            if not force:
                print(f"Error: {existing} samples already exist for {exp_id}")
                print(f"Use --force to delete and regenerate:")
                print(f"  python lw.py exp manifest {exp_id} --force")
                sys.exit(1)
            cursor.execute("DELETE FROM samples WHERE experiment_id = ?", (exp_id,))
            print(f"Cleared {existing} existing samples.")
        # --- resolve strain labels ---
        strain_labels = {}
        all_strain_ids = set()
        # collect from categories
        if "strain" in design["categories"]:
            for s in design["categories"]["strain"]:
                if s is not None:
                    all_strain_ids.add(s)
        # collect from extras
        for extra in design.get("extras", []):
            s = extra.get("strain")
            if s is not None:
                all_strain_ids.add(s)
        # look up each strain
        resolve_parts = []
        for sid in sorted(all_strain_ids):
            cursor.execute("SELECT id, label FROM strains WHERE id = ?", (sid,))
            row = cursor.fetchone()
            if not row:
                print(f"Error: strain '{sid}' not found in database.")
                print(f"Register it first: python lw.py strain add {sid} <genotype>")
                sys.exit(1)
            lbl = row[1] if row[1] else sid  # fall back to ID if no label
            strain_labels[sid] = lbl
            resolve_parts.append(f"{sid}\u2192{lbl}")
        if resolve_parts:
            print(f"Resolving strains: {', '.join(resolve_parts)}")
        # --- generate extras ---
        samples = []
        pos = 1
        for extra in design.get("extras", []):
            sid = extra.get("strain")
            lbl = extra.get("label", "")
            if not lbl and sid:
                lbl = strain_labels.get(sid, sid)
            conditions = dict(extra)  # copy all fields into conditions
            samples.append({
                "position": pos,
                "sample_type": "extra",
                "strain_id": sid,
                "dna": extra.get("dna"),
                "label": lbl,
                "conditions": conditions,
            })
            pos += 1
        n_extras = len(samples)
        # --- generate factorial cross ---
        categories = design["categories"]
        cat_names = list(categories.keys())
        cat_values = [categories[c] for c in cat_names]
        factorial_rows = []
        for combo in itertools.product(*cat_values):
            row = dict(zip(cat_names, combo))
            excluded = False
            for exc in design.get("exclude", []):
                if exc(row):
                    excluded = True
                    break
            if not excluded:
                factorial_rows.append(row)
        # sort by sort_order using category-level ordering
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
            samples.append({
                "position": pos,
                "sample_type": "factorial",
                "strain_id": db_strain_id,
                "dna": row.get("dna"),
                "label": lbl,
                "conditions": row,
            })
            pos += 1
        n_factorial = len(samples) - n_extras
        n_total = len(samples)
        # verify count
        expected = design["expected_samples"]
        if n_total != expected:
            print(f"Error: generated {n_total} samples (expected {expected})")
            print(f"  {n_extras} extras + {n_factorial} factorial = {n_total}")
            # print the factorial rows for debugging
            print()
            print("Factorial rows generated:")
            for i, row in enumerate(factorial_rows, 1):
                print(f"  {i}: {row}")
            print()
            print("Check your categories, excludes, and expected_samples.")
            sys.exit(1)
        print(f"Generated {n_extras} extras + {n_factorial} factorial = {n_total} samples (expected {expected})")
        # --- insert into database ---
        for s in samples:
            cursor.execute(
                "INSERT INTO samples (experiment_id, position, sample_type, strain_id, dna, label, conditions) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (exp_id, s["position"], s["sample_type"], s["strain_id"],
                 s["dna"], s["label"], json.dumps(s["conditions"])),
            )
        # --- determine display columns ---
        # categories except "strain" (shown as strain_id column already)
        cat_display = [c for c in cat_names if c != "strain"]
        # extra-specific keys not already covered
        extra_keys = set()
        for extra in design.get("extras", []):
            for k in extra:
                if k not in ("strain", "label", "dna") and k not in cat_display:
                    extra_keys.add(k)
        extra_keys = sorted(extra_keys)
        # --- print summary table ---
        print()
        hdr = f"  {'pos':>5}  {'type':<10} {'strain':<8} {'label':<14}"
        for c in cat_display:
            hdr += f" {c:<8}"
        for k in extra_keys:
            hdr += f" {k:<12}"
        print(hdr)
        for s in samples:
            cond = s["conditions"]
            line = f"  {s['position']:>5}  {s['sample_type']:<10} {(s['strain_id'] or '-'):<8} {s['label']:<14}"
            for c in cat_display:
                v = cond.get(c)
                line += f" {(str(v) if v is not None else '-'):<8}"
            for k in extra_keys:
                v = cond.get(k)
                line += f" {(str(v) if v is not None else '-'):<12}"
            print(line)
        # --- write manifest.csv ---
        csv_path = os.path.join(folder_path, "manifest.csv")
        csv_columns = ["position", "sample_type", "strain_id", "label"] + cat_display + extra_keys
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(csv_columns)
            for s in samples:
                cond = s["conditions"]
                row = [
                    s["position"],
                    s["sample_type"],
                    s["strain_id"] or "",
                    s["label"],
                ]
                for c in cat_display:
                    v = cond.get(c)
                    row.append(str(v) if v is not None else "")
                for k in extra_keys:
                    v = cond.get(k)
                    row.append(str(v) if v is not None else "")
                writer.writerow(row)
        print()
        print(f"Inserted {n_total} rows into samples table.")
        print(f"Wrote manifest.csv to experiments/{folder_name}/")
    else:
        print(f"Unknown subcommand: exp {subcommand}")
        print("Run 'python lw.py exp' for usage.")
        sys.exit(1)
# --- lw strain ---
elif command == "strain":
    if not subcommand:
        print("Usage: python lw.py strain <subcommand>")
        print("  add <id> <genotype>   Register a new strain")
        print("  show <id>             Show strain details and experiments")
        sys.exit(0)
    # --- lw strain update <id> <field> <value> ---
    if subcommand == "update":
        if len(pos_args) < 3:
            print("Usage: python lw.py strain update <id> <field> <value>")
            sys.exit(1)
        strain_id, field, value = pos_args[0], pos_args[1], pos_args[2]
        # verify strain exists
        cursor.execute("SELECT id FROM strains WHERE id = ?", (strain_id,))
        if not cursor.fetchone():
            print(f"No strain found with ID {strain_id}")
            sys.exit(1)
        # validate field
        cursor.execute("PRAGMA table_info(strains)")
        valid_fields = [r[1] for r in cursor.fetchall() if r[1] != "id"]
        if field not in valid_fields:
            print(f"Unknown field: '{field}'")
            print(f"Valid fields: {', '.join(valid_fields)}")
            sys.exit(1)
        # get old value and update
        cursor.execute(f"SELECT {field} FROM strains WHERE id = ?", (strain_id,))
        old_val = cursor.fetchone()[0]
        # notes: always prefix with timestamp, append if existing
        if field == "notes":
            today = datetime.date.today().strftime("%Y-%m-%d")
            if old_val:
                value = old_val + f"\n[{today}] " + value
            else:
                value = f"[{today}] " + value
        cursor.execute(f"UPDATE strains SET {field} = ? WHERE id = ?", (value, strain_id))
        print(f"Updated {strain_id}")
        if field == "notes":
            print(f"  {field}: appended entry")
        else:
            print(f"  {field}: {old_val if old_val else '-'}  {value}")
    # --- lw strain add <id> <genotype> ---
    if subcommand == "add":
        if len(pos_args) < 2:
            print("Usage: python lw.py strain add <id> <genotype>")
            print('  Example: python lw.py strain add LY456 "MATa orc4-R267A::KanMX ura3 leu2 trp1 his3"')
            sys.exit(1)
        strain_id = pos_args[0]
        genotype = pos_args[1]
        # check for duplicate
        cursor.execute("SELECT id FROM strains WHERE id = ?", (strain_id,))
        if cursor.fetchone():
            print(f"Strain {strain_id} already exists.")
            print(f"Use 'python lw.py strain show {strain_id}' to view it.")
            sys.exit(1)
        cursor.execute(
            "INSERT INTO strains (id, genotype) VALUES (?, ?)",
            (strain_id, genotype),
        )
        print(f"Registered: {strain_id}")
        print(f"Genotype:   {genotype}")
    # --- lw strain show <id> ---
    elif subcommand == "show":
        if len(pos_args) != 1:
            print("Usage: python lw.py strain show <id>")
            sys.exit(1)
        strain_id = pos_args[0]
        cursor.execute("SELECT * FROM strains WHERE id = ?", (strain_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No strain found with ID {strain_id}")
            sys.exit(1)
        col_names = [d[0] for d in cursor.description]
        strain = dict(zip(col_names, row))
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
            print(f"  {label + ':':<{label_width + 2}} {val if val else '-'}")
        # experiments using this strain
        cursor.execute(
            "SELECT e.id, e.type, e.shortname, es.role FROM experiment_strains es JOIN experiments e ON es.experiment_id = e.id WHERE es.strain_id = ?",
            (strain_id,),
        )
        exps = cursor.fetchall()
        if exps:
            print()
            print("  Experiments:")
            for eid, etype, sname, role in exps:
                role_str = f" ({role})" if role else ""
                print(f"    {eid}  {etype:<14} {sname}{role_str}")
    else:
        print(f"Unknown subcommand: strain {subcommand}")
        print("Run 'python lw.py strain' for usage.")
        sys.exit(1)
# --- lw stage ---
elif command == "stage":
    if not subcommand:
        print("Usage: python lw.py stage <subcommand>")
        print("  list                          Show files in staging/")
        print("  assign <file> <exp_id> <desc> Rename, move, and register")
        sys.exit(0)
    # --- lw stage list ---
    if subcommand == "list":
        if not os.path.exists(STAGING_DIR):
            print("Staging directory does not exist.")
            sys.exit(1)
        files = [f for f in os.listdir(STAGING_DIR) if not f.startswith(".")]
        if not files:
            print("Staging is empty.")
            sys.exit(0)
        # load known experiment IDs for matching
        cursor.execute("SELECT id, folder_name FROM experiments")
        experiments = cursor.fetchall()
        print(f"Staging ({len(files)} files):")
        print()
        for fname in sorted(files):
            size = os.path.getsize(os.path.join(STAGING_DIR, fname))
            if size > 1024 * 1024:
                size_str = f"{size / (1024 * 1024):.1f} MB"
            elif size > 1024:
                size_str = f"{size / 1024:.1f} KB"
            else:
                size_str = f"{size} B"
            # check if filename contains any known experiment ID
            matches = [eid for eid, _ in experiments if eid in fname]
            match_str = f"   {', '.join(matches)}" if matches else ""
            print(f"  {fname:<40} {size_str:>10}{match_str}")
    # --- lw stage assign <filename> <exp_id> <descriptor> ---
    elif subcommand == "assign":
        if len(pos_args) < 3:
            print("Usage: python lw.py stage assign <filename> <exp_id> <descriptor>")
            print("  Example: python lw.py stage assign Image_001.tiff LM-0001 gel-01")
            sys.exit(1)
        src_name, raw_id, descriptor = pos_args[0], pos_args[1], pos_args[2]
        # validate source file exists
        src_path = os.path.join(STAGING_DIR, src_name)
        if not os.path.exists(src_path):
            print(f"File not found in staging: {src_name}")
            # list what's there
            files = [f for f in os.listdir(STAGING_DIR) if not f.startswith(".")]
            if files:
                print("Available files:")
                for f in sorted(files):
                    print(f"  {f}")
            sys.exit(1)
        # normalize experiment ID
        if raw_id.upper().startswith("LM-"):
            exp_id = raw_id.upper()
        else:
            try:
                exp_id = f"LM-{int(raw_id):04d}"
            except ValueError:
                print(f"Invalid experiment ID: '{raw_id}'")
                sys.exit(1)
        # validate experiment exists and get folder
        cursor.execute("SELECT folder_name FROM experiments WHERE id = ?", (exp_id,))
        row = cursor.fetchone()
        if not row:
            print(f"No experiment found with ID {exp_id}")
            sys.exit(1)
        folder_name = row[0]
        # build new filename
        today_compact = datetime.date.today().strftime("%Y%m%d")
        _, ext = os.path.splitext(src_name)
        ext = ext.lower()
        new_name = f"{today_compact}_{exp_id}_{descriptor}{ext}"
        # build destination path
        dest_dir = os.path.join(EXPERIMENTS_DIR, folder_name)
        dest_path = os.path.join(dest_dir, new_name)
        if os.path.exists(dest_path):
            print(f"Destination already exists: {new_name}")
            print("Use a different descriptor.")
            sys.exit(1)
        # move file
        shutil.move(src_path, dest_path)
        # register in files table (path relative to LAB_ROOT)
        rel_path = os.path.join("experiments", folder_name, new_name)
        cursor.execute(
            "INSERT INTO files (experiment_id, filename, file_type) VALUES (?, ?, ?)",
            (exp_id, rel_path, None),
        )
        print(f"Filed: {src_name}")
        print(f"     {rel_path}")
    else:
        print(f"Unknown subcommand: stage {subcommand}")
        print("Run 'python lw.py stage' for usage.")
        sys.exit(1)
# --- unknown command ---
else:
    print(f"Unknown command: {command}")
    print("Run 'python lw.py' for usage.")
    sys.exit(1)
# ============================================================
# SECTION 5: CLEANUP
# ============================================================
if conn:
    conn.commit()
    conn.close()
