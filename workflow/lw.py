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
    "purification": "pur", "loading": "load", "gelshift": "gs",
    "atpase": "atp", "genetics": "gen", "computational": "comp",
    "cloning": "clon", "sequencing": "seq",
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
DB_PATH         = os.path.join(LAB_ROOT, "lab.db")
EXPERIMENTS_DIR = os.path.join(LAB_ROOT, "experiments")
COUNTER_FILE    = os.path.join(LAB_ROOT, "counter.txt")


# ---------------------------------------------------------------------------
# C2 -- effects-as-data primitives (the core contract)
# ---------------------------------------------------------------------------

def effects_ok(stdout=None, store=None, exit_code=0):
    return {"ok": True, "stdout": stdout or [], "stderr": [],
            "store": store, "exit_code": exit_code}


def effects_fail(msg, exit_code=1):
    if isinstance(msg, str):
        msg = [msg]
    return {"ok": False, "stdout": [], "stderr": msg,
            "store": None, "exit_code": exit_code}


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
