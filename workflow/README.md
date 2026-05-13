# lw -- Lab Workflow

A personal laboratory information management system built from SQLite,
plain text, and CLI tools. Tracks experiments, strains, files, and
figure provenance for a single-user wet-lab workflow.

Built for one person. Designed to be portable across labs.


## What This Is

Three systems, each doing what it does best:

  SQLite (lab.db)    Structured metadata, cross-references, queries.
  Filesystem         Binary data (TIFFs, gels), scripts, plots.
  Markdown files     Notes, protocols, observations.

Given any figure, trace back to raw data, experiment, strain, and
verification. Given any strain, find all experiments that used it.
Given any experiment, find all downstream outputs. Navigation is by
CLI query, not folder browsing.


## Environment Setup

### Prerequisites

- WSL on Windows (or native Linux/macOS)
- Python 3.8+ managed with uv
- SQLite3
- Neovim >= 0.7 with telescope.nvim (for editor integration)
- rsync (for backup)

### Shell Configuration (lw.sh)

Source lw.sh from your shell rc file. It requires one environment
variable:

  export MC_WINDOWS_USER="yourname"        # or LW_WINDOWS_USER

lw.sh derives all paths and provides these aliases:

  lw             Run lw.py via uv (main entry point for all commands)
  lwnotes        Open lab_notes.md in $EDITOR
  lwdb           Open lab.db in sqlite3 shell
  lwcheck        Run database integrity check
  lwbackup       Run backup script (see Backup section)
  lwcd [id]      CD to lab root, or to experiment folder matching *id*

Tab completion is configured for all lw subcommands and nested
subcommands (exp init, stage auto, etc.).

Examples:

  lw status
  lwcd 0005                   # cd to experiments/*LM-0005*/
  lwdb "SELECT * FROM strains"
  lwnotes                     # quick note in nvim

### Neovim Integration (lw.lua)

Requires env vars LW_LAB_ROOT, LW_DB, LW_NOTES (set by lw.sh).
All keymaps use <leader>l prefix.

  <leader>ln     Open lab_notes.md
  <leader>lh     Prepend today's date header (## YYYY-MM-DD) if missing
  <leader>le     Pick experiment from DB, insert ID into buffer
  <leader>ls     Pick strain from DB, insert ID into buffer
  <leader>lf     Telescope file finder scoped to lab root
                 (excludes .git, __pycache__, .venv, staging)
  <leader>lp     Telescope file finder scoped to protocols/
  <leader>lx     Find files in current experiment folder
                 (auto-detects from buffer path, or prompts to pick)
  <leader>lt     Scan lab_notes.md for TODO lines, populate quickfix

Commands (callable from command mode):
  :LwPickExperiment    Same as <leader>le
  :LwPickStrain        Same as <leader>ls
  :LwTodos             Same as <leader>lt


## Quick Start

Initialize the lab directory:

  lw init

This creates the directory structure, database (10 tables), counter
file (starting at 1), and lab_notes.md.

Register strains:

  lw strain add LY101 "MATa ura3 leu2 trp1 his3"
  lw strain update LY101 label wild-type

Create an experiment:

  lw exp init purification orc-wt --title "Wild-type ORC purification"

This prints the experiment ID (LM-0001), creates the folder, and
shows a next-step hint.

Add metadata:

  lw exp update LM-0001 protein ORC
  lw exp addstrain LM-0001 LY101 expression-strain
  lw exp update LM-0001 notes Started 6L culture at OD 0.8

Check system state:

  lw status


## Command Reference

### System Commands

  lw init
    Bootstrap lab directory and database. Refuses to run if lab.db
    already exists.

  lw status
    System dashboard: experiment counts, staging contents, pending
    expectations, strain summary, recent activity, attention warnings.
    Warnings appear unconditionally and go away only when addressed.

  lw --help
    Show usage summary for all commands.

### Experiment Commands

  lw exp init <type> <shortname> [--title "Human-readable name"]
    Create a new experiment. Assigns next LM-NNNN ID, creates folder.
    Title defaults to shortname if not specified.
    Types: purification, loading, gelshift, atpase, genetics,
           computational, cloning, sequencing
    Prints type-sensitive next-step hint.

  lw exp show <id>
    Full detail view. Accepts "LM-0001" or bare "1".
    Shows age (active experiments) or duration (complete).
    Includes strains, links, files, type-specific details, samples.

  lw exp update <id> <field> <value...>
    Update a field. Multi-word values accepted without quoting.
    Scans experiments, purification_details, and assay_details tables.
    Notes field: always appends with [YYYY-MM-DD] timestamp prefix.
    All other fields: overwrite.
    Protected fields that cannot be modified:
      id, experiment_id, created_at, folder_name,
      date_started, date_completed, status
    For status changes, use "lw exp complete".
    For purification records, set "protein" first (creates the row).

  lw exp complete <id>
    Mark experiment done. Sets date_completed, shows duration.
    Already-complete experiments are skipped.

  lw exp link <from> <to> <relationship>
    Connect two experiments.
    Standard relationships: uses_prep_from, replicate_of, follow_up_to
    Non-standard relationships accepted with a warning.
    uses_prep_from auto-populates assay_details.protein_prep.
    Warns if target of uses_prep_from is not a purification.

  lw exp addstrain <exp_id> <strain_id> [role]
    Associate a registered strain with an experiment.
    Strain must exist in database (register first with strain add).

  lw exp delete <id> [--confirm] [--keep-folder]
    Dry-run by default: shows what would be deleted.
    --confirm executes the deletion.
    --keep-folder preserves the experiment directory.
    Cleans downstream protein_prep references before cascade.
    Cascade covers: experiment_strains, experiment_links, files,
    purification_details, assay_details, figure_panels, samples,
    stage_expectations.

  lw exp find [query...] [--type T] [--status S] [--strain S]
    Search experiments by text (matches id, shortname, title, notes),
    type, status, and/or strain. Multi-word queries accepted.
    Displays title when it differs from shortname.
    If --strain filter finds nothing, checks whether strain exists.

  lw exp list [--type T]
    Type-grouped listing with type-specific summaries.
    Purification: protein name + yield.
    Assays: sample count + protein_prep with protein name lookup.

  lw exp manifest <id> [--force]
    Generate sample manifest from design.py in experiment folder.
    Creates rows in samples table and writes manifest.csv.
    All strains referenced must be registered in database first.
    --force required to regenerate (deletes existing samples first).
    Extras are positioned first (positions 1..N), factorial samples follow.
    See docs/design_template.py for design.py format.

### Strain Commands

  lw strain add <id> <genotype>
    Register a new strain. Validates ID format (alphanumeric + hyphens).
    Warns on non-standard patterns (typical: letters + digits, e.g. LY456).
    Prompts to set a label afterward.

  lw strain show <id>
    Show all strain fields and linked experiments.

  lw strain update <id> <field> <value...>
    Update a strain field. Multi-word values accepted.
    Notes field: appends with timestamp (same as exp update).
    Fields: genotype, label, parent, construction, selection_marker,
            verification, storage_location, notes.

  lw strain list [--no-label]
    List all strains with labels, genotypes, and experiment counts.
    --no-label: show only strains missing labels (pre-manifest QC).

### Staging Commands

  lw stage list
    Show files in staging/ with sizes and parsed experiment/slot matches.
    Also shows pending expectations with match status.

  lw stage assign <filename> <exp_id> <descriptor>
    Manually rename, move, and register a file from staging.
    Descriptor: alphanumeric + hyphens (validated).
    Confirms experiment identity before filing.
    File renamed to YYYYMMDD_LM-NNNN_descriptor.ext (filing date).

  lw stage expect <exp_id> <descriptor> [descriptor2...]
    Pre-register expected files. Assigns slot numbers (s1, s2...).
    Prints instrument naming convention for each slot.
    Checks for duplicate descriptors (within args and existing pending).

  lw stage expect --list [exp_id]
    Show all expectations (or filter by experiment).
    Displays status: pending, filed, cancelled.

  lw stage expect --cancel <exp_id> [sN]
    Cancel pending expectations. Without slot: cancels all pending.
    With slot (s1, S1, or bare 1): cancels single expectation.

  lw stage auto [--confirm]
    Match staged files to pending expectations by experiment ID + slot.
    Dry-run by default: shows matched, skipped, and unmatched files.
    Annotates inferred slot matches (single-pending-per-experiment).
    --confirm executes: renames, moves, registers, marks expectations filed.


## Directory Structure

  lab/
  |-- experiments/              One folder per experiment (flat)
  |   +-- YYYYMMDD_LM-NNNN_CAT_shortname/
  |       +-- design.py         Factorial design (for manifest)
  |-- protocols/                Markdown, one per protocol
  |-- publications/             One subfolder per publication
  |-- resources/                Shared computational assets
  |-- staging/                  Instrument dump area
  |-- docs/                     Documentation
  |   |-- workflow_recipes.md
  |   +-- design_template.py
  |-- lab_notes.md              Daily notes, date headers
  |-- lab.db                    SQLite database
  |-- counter.txt               Next experiment number
  |-- command_log.txt           Append-only success/failure log
  +-- README.md                 This file


## Database Schema

10 tables. Full DDL in lw.py Section 1 (SCHEMA_SQL).

  experiments          Core record: id, type, title, shortname, dates,
                       status, folder_name, protocol_id, notes.
  strains              Strain registry: id, genotype, label, parent,
                       construction, marker, verification, storage, notes.
  experiment_strains   Many-to-many: experiment <-> strain with role.
  experiment_links     Experiment relationships (uses_prep_from, etc.).
  files                Registered files with paths relative to lab root.
  purification_details Type-specific fields for purification experiments.
  assay_details        Type-specific fields for assay experiments.
  figure_panels        Links publication figures to experiments and files.
  samples              Factorial manifest rows (position, strain, conditions).
  stage_expectations   Pre-registered file expectations with slot numbers.


## Naming Conventions

### Experiment IDs: LM-NNNN
  Initials + zero-padded 4-digit sequential. Globally unique, portable,
  never reused. Counter stored in counter.txt.

### Folder Names: YYYYMMDD_LM-NNNN_CAT_shortname
  Category codes:
    pur   purification       load  loading
    gs    gelshift           atp   atpase
    gen   genetics           comp  computational
    clon  cloning            seq   sequencing

### Shortnames
  Filesystem-safe: lowercase alphanumeric + hyphens. Must start with
  a letter or digit. Set at creation, used in folder name.

### Titles
  Human-readable, for recall. Can contain spaces and mixed case.
  Defaults to shortname. Set via --title flag on exp init, or
  updated via "lw exp update LM-NNNN title descriptive name here".

### File Names: YYYYMMDD_LM-NNNN_descriptor.ext
  Date is the filing date (when the file was assigned to the experiment),
  not the instrument date. Descriptor is lowercase alphanumeric + hyphens.

### Derived Files: YYYYMMDD_prefix_LM-NNNN_descriptor.ext
  Prefixes: bc-enhanced, quantified, cropped, plot.
  Example: 20260513_bc-enhanced_LM-0005_gel-coomassie.tif

### Time Disambiguation
  When multiple files share the same date, experiment, and descriptor,
  append _HHMM before the extension:
  20260513_LM-0005_gel-coomassie_1430.tif

### Retakes
  Append -v2, -v3 to descriptor. Both versions kept. Not system-enforced.
  Example: gel-coomassie-v2

### Instrument Naming Convention
  At all instruments, name files:  lemr LM-NNNN sN
    lemr     Personal search prefix (findable on instrument software)
    LM-NNNN  Experiment ID (parsed by staging system)
    sN       Slot number from stage expect

  No underscores (Imager 680 restriction). Spaces as delimiters.
  Instrument software may append its own suffix (timestamps, channels).


## Key Behaviors

### Notes Append
  "lw exp update <id> notes <text>" always prepends [YYYY-MM-DD] and
  appends to existing content. All other fields overwrite.
  Same behavior for "lw strain update <id> notes <text>".

### Protected Fields
  Cannot be modified via exp update:
    id, experiment_id, created_at, folder_name,
    date_started, date_completed, status
  Status changes go through "lw exp complete".

### Dry-Run Defaults
  "lw exp delete" and "lw stage auto" are dry-run by default.
  Both require --confirm to execute changes.
  Dry runs show exactly what would happen, including cascade targets.

### Filing Dates
  stage assign and stage auto use today's date for the renamed file,
  not the date embedded in the original instrument filename.
  Rationale: filing date is consistent and predictable.

### Command Log
  Successful commands: logged to command_log.txt at script exit.
  Failed commands: logged via bail() before exit.
  Format: ISO timestamp <tab> command args <tab> ok|fail
  Commands that exit early with nothing to do (dry-run display,
  empty results, already-complete) are not logged.

### Transaction Model
  All DB changes accumulate during command execution and commit once
  at script end. Filesystem operations (file moves, directory creation)
  are immediate and non-reversible. Operations ordered to minimize
  unrecoverable states: prefer DB-then-filesystem.

### Sub-tasks
  Bradfords, OD readings, agarose gels for verification -- these are
  files within the parent experiment folder, not separate experiments.
  A purification is one experiment; the bradford measuring its fractions
  is a file in that experiment's folder.


## Backup

Backup is handled by lw_backup.sh, invoked via the lwbackup alias.

Requires: LW_LAB_ROOT, LW_DB, LW_BACKUP_DIR env vars, plus sqlite3
and rsync.

Process:
  1. Creates a consistent DB snapshot using "sqlite3 .backup".
  2. rsync -av (additive, never deletes) syncs lab/ to backup dir,
     excluding lab.db and the temp snapshot.
  3. Copies snapshot to backup dir as lab.db. Removes temp file.
  4. Weekly: copies lab.db to snapshots/lab_YYYY-WNN.db (if not
     already present). Prunes snapshots older than 84 days.

Dry-run mode: lwbackup --dry-run (simulates all actions, no changes).

The backup destination is typically a Dropbox-synced folder, providing
off-machine redundancy. The weekly snapshots provide point-in-time
recovery independent of Dropbox versioning.


## Design Principles

Single-file tool. lw.py is one Python file, flat procedural, standard
library only. No OOP, no classes, no argparse. sys.argv parsing,
sqlite3 for database, if/elif dispatch. bail() for all error exits.

Pull features. Build when friction demands, not in anticipation.
The system grows by discovering what's missing during real use.

Invisible to collaborators. The PI does not use git, LaTeX, or
programmatic tools. Shared artifacts are Word-compatible. The system
never leaks into shared workflows.

Forward-only migration. New experiments use the new system. Old data
stays where it is. No bulk import required.


## Known Limitations

- protocol_id is a soft foreign key (no validation that file exists).
- JSON fields (purification_details.columns, assay_details.conditions)
  written as raw strings via exp update. No append/parse logic.
- No strain find command. Use strain list or raw SQL.
- Schema migration for existing databases requires manual ALTER TABLE
  or reinitialization. Not urgent while DB has minimal production data.
- Status warnings (no strains, no notes) fire on freshly created and
  computational experiments. Accepted trade-off.
- Fig commands (fig link, fig trace) not yet implemented.
  Manual tracking via lab_notes.md until then.
