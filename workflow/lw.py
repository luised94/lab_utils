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
    sys.exit(0)

command = args[0]

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
