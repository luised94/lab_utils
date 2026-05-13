# Lab Workflow System

A personal laboratory information management system (LIMS) built on SQLite,
plain text, and a single CLI tool. Tracks experiments, strains, instrument
files, and sample manifests for biochemistry research.

Built for one person. Designed to be portable across labs.


## What It Does

Given any figure, trace back to raw data, experiment, strain, and
verification. Given any strain, find all experiments that used it. Given any
experiment, find all downstream outputs. Navigation is by CLI query, not
folder browsing.

Three systems, each doing what it does best:

    SQLite (lab.db)     Structured metadata, cross-references, queries
    Filesystem          Binary data (TIFFs, gels), scripts, plots
    Markdown files      Notes, protocols, observations


## Quick Start

### First-time setup

    python lw.py init

This creates the directory structure, SQLite database, counter file, and
lab notes file. Run it once.

### Create your first experiment

    python lw.py exp init purification orc4-wt-prep
    # Created: LM-0001
    # Folder:  20260513_LM-0001_pur_orc4-wt-prep/
    # Type:    purification
    # Next: lw exp update LM-0001 protein <name>

    python lw.py exp update LM-0001 protein ORC-WT

### Register a strain

    python lw.py strain add LY100 "MATa ura3 leu2 trp1 his3"
    # Registered: LY100
    # Set a label: lw strain update LY100 label <short-name>

    python lw.py strain update LY100 label "wild-type"

### Link a strain to an experiment

    python lw.py exp addstrain LM-0001 LY100

### Check system state

    python lw.py status


## Commands

All commands go through a single entry point: `python lw.py <command>`.
For brevity, examples below use `lw` as an alias.

### General

    lw --help                   Show usage and all commands
    lw init                     Bootstrap lab directory and database (run once)
    lw status                   System dashboard: counts, staging, warnings

### Experiments

    lw exp init <type> <shortname> [--title "Human-readable name"]

Create a new experiment. Assigns the next sequential ID (LM-0001, LM-0002,
...), creates a folder, and inserts a database record. The shortname becomes
part of the folder name and must be lowercase alphanumeric with hyphens.
The title defaults to the shortname if not specified.

Valid types: purification, loading, gelshift, atpase, genetics,
computational, cloning, sequencing.

After creation, the system prints a context-aware next step:
- purification: set the protein name
- loading/gelshift/atpase: add a strain
- other types: add notes

    lw exp show <id>

Display full experiment details including type-specific fields, linked
strains, linked experiments, registered files, and samples. Active
experiments show age in days; completed experiments show duration.

Accepts full ID (LM-0001) or bare number (1).

    lw exp update <id> <field> <value>

Update a single field. The system automatically finds which table owns
the field (experiments, purification_details, or assay_details). If a
purification or assay record does not exist yet, it creates one.

Multi-word values do not need quoting:

    lw exp update LM-0001 notes gel looks clean today

Some fields cannot be changed this way. See "Protected Fields" below.

    lw exp complete <id>

Mark an experiment as done. Sets the completion date and prints the
duration. An already-complete experiment is left unchanged.

    lw exp link <from_id> <to_id> <relationship>

Connect two experiments. Standard relationships:
- uses_prep_from: downstream assay uses protein from a purification
- replicate_of: same experiment repeated
- follow_up_to: continuation of earlier work

When linking with uses_prep_from, the system auto-populates the
protein_prep field in assay_details and warns if the target is not
a purification experiment.

    lw exp addstrain <exp_id> <strain_id> [role]

Associate a registered strain with an experiment. The strain must already
exist in the database (register it first with `strain add`).

    lw exp delete <id> [--confirm] [--keep-folder]

Delete an experiment and all related records. By default this is a dry
run -- it shows what would be deleted without doing anything. Pass
--confirm to execute. Pass --keep-folder to delete the database records
but preserve the experiment directory and its files.

The dry run shows:
- Related records in all tables (strains, links, files, samples, etc.)
- Downstream protein_prep references that would be cleared
- Folder contents

    lw exp find [query] [--type T] [--status S] [--strain S]

Search experiments. The query matches against ID, shortname, title, and
notes. Filters can be combined. Multi-word queries do not need quoting.
Results are flat, sorted by date (newest first). Shows title when it
differs from shortname.

If --strain returns no results, the system checks whether the strain
exists in the database and reports if it does not.

    lw exp list [--type T]

Browse experiments grouped by type. Each group shows type-specific
summaries:
- purification: protein name and yield
- loading/gelshift/atpase: sample count and linked protein prep
- other types: basic info only

    lw exp manifest <id> [--force]

Generate a sample manifest from a design.py file in the experiment
folder. See "Sample Manifests" below for the design.py format.

If samples already exist for this experiment, the command aborts unless
--force is passed, which deletes existing samples first.

### Strains

    lw strain add <id> <genotype>

Register a new strain. The ID must be alphanumeric with optional hyphens
(e.g., LY100, LY456, orc4-R267A-1). The system warns if the ID does
not follow the typical format of letters followed by digits (e.g., LY456)
but registers it anyway.

After registration, the system prompts you to set a label.

    lw strain show <id>

Display strain details and all experiments using this strain.

    lw strain update <id> <field> <value>

Update a strain field. Multi-word values do not need quoting.
Valid fields: genotype, label, parent, construction, selection_marker,
verification, storage_location, notes.

    lw strain list [--no-label]

List all strains with labels, genotypes, and experiment counts. Pass
--no-label to show only strains missing labels (useful before generating
manifests).

### Staging

The staging system handles instrument file management. Files from
instruments go into the staging/ directory, then get renamed and moved
to the correct experiment folder.

Two workflows: manual assignment and pre-registered expectations.

#### Manual assignment

    lw stage list

Show files in staging/ with sizes. Also shows pending expectations and
which files match them.

    lw stage assign <filename> <exp_id> <descriptor>

Rename a file from staging/ and move it to the experiment folder. The
descriptor becomes part of the new filename. Before filing, the system
confirms the experiment identity (ID, shortname, type).

The descriptor must be alphanumeric with hyphens (e.g., gel-coomassie,
scan-01, blot-anti-orc4).

Example:

    lw stage assign "lemr LM-0005 s1.tif" LM-0005 gel-coomassie
    # Experiment: LM-0005 (orc4-loading, loading)
    # Filed: lemr LM-0005 s1.tif
    #     -> experiments/20260513_LM-0005_load_orc4-loading/20260513_LM-0005_gel-coomassie.tif

#### Pre-registered expectations

Register what you plan to image before going to the instrument.

    lw stage expect <exp_id> <descriptor> [descriptor2 ...]

Pre-register one or more expected files. Each gets a slot number (s1, s2,
...). The system prints the name to use at the instrument.

    lw stage expect LM-0005 gel-coomassie gel-silver
    #   s1  gel-coomassie
    #   s2  gel-silver
    #   At instrument:
    #     lemr LM-0005 s1
    #     lemr LM-0005 s2
    #   2 total pending for LM-0005

At the instrument, name each file using the printed convention:
`lemr LM-NNNN sN` (no underscores -- some instruments forbid them).

After imaging, copy files to staging/ and run:

    lw stage auto

Dry run: shows which files match which expectations. Then:

    lw stage auto --confirm

Executes the moves. Files without a slot number but with only one pending
expectation for that experiment are matched automatically (annotated as
"slot inferred" in the dry run).

    lw stage expect --list [exp_id]

Show all expectations (pending and filed). Filter by experiment if given.

    lw stage expect --cancel <exp_id> [sN]

Cancel pending expectations. Cancel a single slot or all pending for
an experiment.


## Directory Structure

    lab/
    |-- experiments/              One folder per experiment (flat)
    |   +-- YYYYMMDD_LM-NNNN_CAT_shortname/
    |       +-- design.py         For factorial assays (used by exp manifest)
    |-- protocols/                Markdown files, one per protocol
    |-- publications/             One subfolder per publication
    |-- resources/                Shared computational assets
    |-- staging/                  Instrument file dump area
    |-- docs/                     Documentation
    |-- lab_notes.md              Single file, date headers
    |-- lab.db                    SQLite database
    |-- counter.txt               Next experiment number
    |-- command_log.txt           Success/failure log
    +-- README.md                 This file


## Database Schema

Ten tables in lab.db:

    experiments          Core experiment records (id, type, title, status, dates)
    strains              Strain registry (id, genotype, label, construction, etc.)
    experiment_strains   Many-to-many: which strains are in which experiments
    experiment_links     Directed links between experiments (uses_prep_from, etc.)
    files                Registered files with paths relative to lab root
    purification_details Per-experiment: protein, yield, columns, fractions
    assay_details        Per-experiment: protein_prep link, DNA template, conditions
    figure_panels        Publication figure provenance (panel -> experiment + file)
    samples              Per-experiment sample manifests (position, strain, conditions)
    stage_expectations   Pre-registered staging slots (pending/filed/cancelled)

The experiments table is the hub. Type-specific metadata lives in
purification_details (for purification experiments) and assay_details
(for loading, gelshift, atpase experiments). The `exp update` command
searches across all three tables automatically.


## Naming Conventions

### Experiment IDs: LM-NNNN

Initials + zero-padded 4-digit sequential number. Globally unique, never
reused. The counter lives in counter.txt and increments on each `exp init`.

### Folder names: YYYYMMDD_LM-NNNN_CAT_shortname

Category codes map to experiment types:

    purification -> pur       genetics     -> gen
    loading      -> load      computational -> comp
    gelshift     -> gs        cloning      -> clon
    atpase       -> atp       sequencing   -> seq

### File names: YYYYMMDD_LM-NNNN_descriptor.ext

The filing date is the date the file was assigned (today), not the date
the instrument captured it. This keeps file naming consistent regardless
of when files are moved out of staging.

Derived files use a prefix: YYYYMMDD_prefix_LM-NNNN_descriptor.ext

    bc-enhanced    Brightness/contrast adjusted
    quantified     Quantification output
    cropped        Cropped region
    plot           Plot or graph

If a time disambiguation is needed, append _HHMM before the extension.

### Shortnames

Lowercase alphanumeric characters and hyphens only. Must start with a
letter or digit. Used in folder names and as a filesystem-safe identifier.
Examples: orc4-wt-prep, mcm-loading-01, ars1-gelshift.

### Titles

Free-form human-readable names for recall. Set with `--title` on
`exp init` or via `exp update <id> title "descriptive name"`. Default
to the shortname if not specified. Displayed in search results and
listings when different from the shortname.

### Retakes

Not system-enforced. Append -v2, -v3 to the descriptor when re-imaging.
Both versions are kept in the experiment folder.

### Instrument naming

At all instruments, name files: `lemr LM-NNNN sN`

- lemr: personal search prefix (findable in instrument software)
- LM-NNNN: experiment ID (parsed by the staging system)
- sN: slot number from pre-registered expectation

No underscores (Imager 680 restriction). Spaces as delimiters.


## Sample Manifests

For factorial assays (loading, gelshift, ATPase), define the experimental
design in a `design.py` file inside the experiment folder, then generate
the manifest:

    lw exp manifest LM-0050
    lw exp manifest LM-0050 --force    # regenerate, replacing existing

### design.py format

The file must define a DESIGN dict:

    DESIGN = {
        "experiment_id": "LM-0050",       # must match the experiment
        "expected_samples": 28,           # total: extras + factorial
        "extras": [                       # optional, positioned first
            {"label": "input", "strain": "LY100", "amount": "5% total"},
            {"label": "input", "strain": "LY456", "amount": "5% total"},
        ],
        "categories": {                   # factorial dimensions
            "strain": ["LY100", "LY456", "LY789", None],
            "dna": ["ARS1", "A-B2-"],
            "kglut": ["250", "300"],
        },
        "sort_order": ["strain", "dna", "kglut"],  # must match category keys
        "exclude": [                      # optional, lambdas returning True to drop
            lambda r: r["strain"] is None and r["kglut"] != "250",
        ],
    }

Required keys: experiment_id, expected_samples, categories, sort_order.
Optional keys: extras, exclude.

- Categories define the factorial dimensions. Use None for no-protein
  controls (rendered as "no protein" in the manifest).
- sort_order must list exactly the same keys as categories.
- exclude lambdas receive a dict of {category: value} and return True
  to drop that combination.
- extras are non-factorial samples (inputs, markers). Each must have a
  "label" key. Optional: strain, dna, amount, or any other key (stored
  in conditions JSON).
- expected_samples is the total count (extras + factorial after excludes).
  The command aborts if the generated count does not match.

All strains referenced in categories or extras must be registered in the
database before running the manifest command.

The manifest command:
1. Loads design.py
2. Resolves strain labels from the database
3. Generates extras (positions 1..N)
4. Generates factorial cross product, applies excludes, sorts
5. Verifies total matches expected_samples
6. Inserts into samples table
7. Writes manifest.csv to the experiment folder

See docs/design_template.py for an annotated template.


## Key Behaviors

### Notes always append

The `notes` field on experiments and strains always appends with a
date prefix. It never overwrites:

    lw exp update LM-0001 notes gel looks clean
    # notes: appended [2026-05-13] gel looks clean

    lw exp update LM-0001 notes repeated with fresh buffer
    # notes: appended [2026-05-13] repeated with fresh buffer

All other fields overwrite the previous value.

### Protected fields

These fields cannot be changed via `exp update`:

    id, experiment_id, created_at, folder_name,
    date_started, date_completed, status

To change status, use `exp complete`. The other fields are immutable
after creation.

### Dry-run defaults

Two commands default to dry run (show what would happen, do nothing):
- `exp delete` -- pass --confirm to execute
- `stage auto` -- pass --confirm to execute

### Filing date convention

When a file is assigned from staging, the date in the new filename is
today's date, not the date shown in the instrument's original filename.
This keeps naming predictable regardless of when files leave staging.

### Command log

Every successful command is logged to lab/command_log.txt:

    2026-05-13T14:30:00  exp init purification orc4-prep  ok

Failed commands (those that exit via bail()) are also logged:

    2026-05-13T14:31:00  exp update LM-9999 protein ORC  fail

The log is append-only. Commands that exit early with nothing to do
(dry-run display, empty results, already-complete experiments) are
not logged.


## Design Principles

**SQLite from day 1.** No flat files for metadata. The database schema
is the format. CLI insertion is lower friction than file editing.

**Flat experiment directories.** No hierarchy by assay type. The
category code in the folder name provides soft filtering. Linked
experiments (e.g., a purification and its downstream assays) live in
separate folders and are connected via experiment_links in the database.

**Single Python file.** lw.py uses only the standard library. No OOP,
no classes, no argparse, no external dependencies. Data-oriented
procedural style.

**Forward-only migration.** Old experiments stay in the old system. New
experiments use the new system. No bulk import.

**Pull features.** Commands are built when friction demands them, not
speculatively. If something is not implemented, it was not yet needed.

**Invisible to collaborators.** The system is personal. Shared artifacts
(manuscripts, figures) are Word-compatible. Collaborators never need to
interact with the CLI, database, or directory structure.

**Transaction model.** All database changes accumulate during command
execution and commit once at the end. Filesystem operations (file moves,
directory creation) are immediate. A crash after a filesystem change but
before commit leaves the filesystem changed and the database unchanged.
For a single-user CLI tool, this is acceptable.


## Environment

### Required environment variables

Set one of:
- LAB_ROOT: full path to the lab directory
- MC_WINDOWS_USER: Windows username (constructs /mnt/c/Users/$USER/Desktop/lab)

### Alias setup (optional)

Add to your shell profile:

    alias lw='python /path/to/lw.py'


## Not Yet Implemented

These features are planned but do not exist yet:

- fig link / fig trace: figure provenance commands (figure_panels table
  exists but has no CLI commands)
- Strain import from old system
- Old experiment registration/migration
- Task tracking module
- Inventory management (currently in Excel)
- Document compilation pipeline


## Known Limitations

- protocol_id is a soft foreign key (no validation that the file exists)
- JSON fields (purification columns, assay conditions) are written as raw
  strings via exp update; no structured editing
- No strain find command (strain list covers basic browsing)
- Schema migration for existing databases requires manual intervention
  or reinitialization
- Status warnings (no strains, no notes) fire on freshly created and
  computational experiments; this is by design
