# Workflow Recipes

Recipes for the current `lw` tool: six commands (`init`, `new`, `link`,
`show`, `list`, `complete`) over two tables (`experiments`,
`experiment_links`). Experiment details beyond identity now live in Excel,
not the database -- `lw` owns the folder, the id, the status, and the
provenance edges; everything else is yours to track outside it.

Source of truth for command behavior is `lw.py` (`COMMAND_USAGE`). The
pre-migration workflow (strains, staging, manifests, figure panels) is
retired; its old recipes and the code that drove them are parked in
`archive/` if you ever reconstitute a subsystem.

Aliases referenced below come from `lw.sh`: `lw` (the tool), `lwdb` (sqlite
shell on the database), `lwcheck` (integrity + FK check on live tables),
`lwnotes` (open lab notes), `lwbackup` (run the backup), `lwcd` (cd into an
experiment folder by fuzzy match).


## Recipe 1: Register and close an experiment

Goal: the full lifecycle of one experiment through `lw` -- create it, look
at it, finish it.

### One-time setup

  lw init

Creates the lab directory, database, and id counter. Safe to re-run; it
never resets the counter.

### Create the experiment

  lw new purification orc-r267a --title "Orc4-R267A ORC purification"

Output:

  Created LM-0005: 20260513_LM-0005_pur_orc-r267a

`new <type> <shortname> [--title ...]` makes the folder, the database row,
and the id in one step. Valid types (with their folder codes): purification
(pur), loading (load), gelshift (gs), atpase (atp), genetics (gen),
computational (comp), cloning (clon), sequencing (seq). The folder name is
`YYYYMMDD_LM-NNNN_<code>_<shortname>` -- date-front-loaded so lexical sort is
chronological.

### Do the work

The experiment's folder now exists under `experiments/`. Put gels, data,
design files, and notes there directly. Detailed metadata (strains, yields,
conditions, sample manifests) goes in Excel or your notes -- `lw` does not
track it anymore. Jump to the folder anytime:

  lwcd orc-r267a        # fuzzy-matches the folder name

### Review

  lw show LM-0005

Output:

  LM-0005: Orc4-R267A ORC purification
    type:      purification
    shortname: orc-r267a
    status:    active
    started:   2026-05-13
    folder:    20260513_LM-0005_pur_orc-r267a
    Files: 6 (last: 20260515_LM-0005_gel-coomassie.tif)

Add `--files` for the full file listing instead of just the count and
newest:

  lw show LM-0005 --files

If the folder was moved or deleted, the Files line reads "folder not found"
rather than erroring -- `show` still prints the row.

### Close it

  lw complete LM-0005

Sets status to complete and stamps the completion date. Running it on an
already-complete experiment is rejected (no-effect), so it is safe to retry.

  Nvim shortcut: <leader>lx while editing any file inside the experiment
  folder telescopes through all files in that experiment.


## Recipe 2: Linking experiments (provenance)

Goal: record directed edges between experiments to build a provenance chain.

  lw link LM-0008 LM-0005 uses_prep_from

Reads as: LM-0008 uses the prep from LM-0005. `link <from> <to>
<relationship>` records a directed edge. Output:

  Linked LM-0008 -> LM-0005 (uses_prep_from)

### The three relationships

The relationship must be one of exactly three:

  uses_prep_from     downstream assay uses an upstream prep
  replicate_of       a repeat of an earlier experiment
  follow_up_to       a follow-up to an earlier experiment

These are plain edge labels with no side effects (`uses_prep_from` no longer
auto-populates anything -- that behavior retired with the details tables).
Anything outside the three is REJECTED, not warned -- the link is not
recorded. To add a fourth relationship, extend `ALLOWED_LINK_RELATIONSHIPS`
in `lw.py` (a deliberate one-line edit).

`link` rejects four cases with no effect: a malformed id, an experiment that
does not exist, an unknown relationship, and a duplicate of an existing edge.
A reverse-direction edge is not a duplicate.

### Viewing links

  lw show LM-0008

The links appear split by direction:

    links out:
      -> LM-0005 (uses_prep_from)
    links in:
      <- LM-0010 (replicate_of)
      <- LM-0012 (follow_up_to)


## Recipe 3: Listing and filtering

Goal: see your experiments, optionally narrowed.

  lw list

Output is a flat list sorted by id (id, status, type, title):

  LM-0003  complete  purification  ORC wild-type purification
  LM-0005  active    purification  Orc4-R267A ORC purification
  LM-0008  active    loading       MCM loading: WT vs R267A ORC

Filter by type, status, or both (absent filter means "match everything"):

  lw list --type purification
  lw list --status active
  lw list --type loading --status active

This pairs with `complete`: after closing experiments, `lw list --status
active` shows only what is still open.

### Raw queries when you need more

The two live tables are open to SQL via `lwdb`. Experiments linked to a
purification:

  lwdb "SELECT from_experiment, relationship
        FROM experiment_links
        WHERE to_experiment = 'LM-0005'"

Active experiments started this month:

  lwdb "SELECT id, title FROM experiments
        WHERE status = 'active' AND date_started >= '2026-06-01'
        ORDER BY date_started DESC"

Only `experiments` and `experiment_links` are live. The retired tables
(strains, files, samples, purification_details, and the rest) may still
exist in older databases but are not maintained -- do not query them as
current.


## Recipe 4: Daily workflow (Neovim + CLI)

Goal: the typical pattern combining nvim keymaps and the CLI. All keymaps
are under `<leader>l`.

### Review and plan

  lw list --status active     # what is still open
  lwnotes                     # open lab notes
  <leader>lh                  # insert today's date header

Write the day's plan in lab notes. Reference experiments inline with
`<leader>le` to pick an id from the database and insert it.

### While working

  <leader>le     Pick an experiment, insert its id
  <leader>lf     Find files anywhere in the lab root (telescope)
  <leader>lp     Find and open a protocol (telescope)
  <leader>lx     Telescope through files in the current experiment folder
  <leader>ln     Open lab notes
  <leader>lt     List TODOs from lab notes in the quickfix list

`<leader>lx` detects the experiment from the buffer's path (it must sit
inside an `experiments/.../` folder). Outside one, it prompts you to pick an
active experiment first.

### End of day

  lw complete LM-0005         # if an experiment finished
  <leader>lt                  # review outstanding TODOs
  lwbackup                    # back everything up


## Recipe 5: Backup

Goal: maintain consistent, recoverable backups of the lab directory.

  lwbackup

This runs `lw_backup.sh`, which:

  1. Snapshots the database safely (sqlite3 .backup, not a raw file copy)
  2. Syncs all files via rsync, additively (never deletes from the backup)
  3. Creates weekly database snapshots (lab_YYYY-WNN.db) automatically
  4. Prunes weekly snapshots older than 84 days

Dry run first if you want to see what would happen:

  lwbackup --dry-run

### Recovery

  1. Copy the backup directory contents to a new lab location.
  2. The database in the backup is a consistent snapshot, ready to use.
  3. For point-in-time recovery, weekly snapshots are under snapshots/.

Because rsync is additive, a file deleted from the working directory is
still preserved in the backup. Back up after any session that creates or
modifies data; the `lwbackup` alias makes it a one-command habit.


## Recipe 6: Health check

Goal: confirm the database is structurally sound.

  lwcheck

Runs `PRAGMA integrity_check` (global) plus `foreign_key_check` scoped to the
two live tables. It is scoped deliberately: an unscoped FK check would walk
the retired tables, which still carry foreign keys into `experiments` and
would report violations you no longer care about. A clean run prints nothing
alarming; investigate any output.
