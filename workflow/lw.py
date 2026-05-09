# ============================================================
# SECTION 0: IMPORTS
# ============================================================
import sqlite3
import sys
import os
import datetime

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
    print(f"Tables: {', '.join(t for t in ['experiments', 'strains', 'experiment_strains', 'experiment_links', 'files', 'purification_details', 'assay_details', 'figure_panels'])}")

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
                print(f"     {target} ({rel})")
            for source, rel in links_in:
                print(f"     {source} ({rel})")

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

    else:
        print(f"Unknown subcommand: exp {subcommand}")
        print("Run 'python lw.py exp' for usage.")
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
