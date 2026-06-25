# lw -- personal lab experiment tracker

A small command-line tool that owns two things and nothing else:

1. **Experiment identity** -- creates experiment folders in a fixed naming
   convention and records each experiment as a row in a SQLite database.
2. **The link graph** -- records directed provenance relationships between
   experiments so you can trace what came from what.

Everything else about an experiment (protein yields, assay conditions,
sample layouts, strains) lives in Excel sheets you maintain by hand. This
tool deliberately does not store that data. It is two tables and five
commands.

---

## Install / setup

`lw` is a single Python 3 file (`lw.py`) with no third-party dependencies
(standard library only: `sqlite3`, `os`, `sys`, `datetime`). Copy it
somewhere on your `PATH`, or run it directly with `python3 lw.py ...`.

You must tell `lw` where your lab directory lives via one environment
variable. It checks them in this precedence order:

| Variable | Meaning |
|---|---|
| `LW_LAB_ROOT` | Absolute path to your lab directory. Wins if set. |
| `LAB_ROOT` | Legacy alias for `LW_LAB_ROOT`. |
| `LW_WINDOWS_USER` | Your Windows username; derives `/mnt/c/Users/<user>/Desktop/lab` (the WSL convention). |
| `MC_WINDOWS_USER` | Legacy alias for `LW_WINDOWS_USER`. |

If none is set, `lw` exits with a message naming both options. For example:

```sh
export LW_LAB_ROOT="$HOME/lab"
python3 lw.py init
```

---

## Quick start

```sh
lw init                                  # create dirs, database, counter
lw new purification orc4r-prep --title ORC4R prep round two
lw new cloning pET28-orc
lw link 1 2 follow_up_to                 # LM-0002 follows up on LM-0001
lw show 1                                # row + links + file summary
lw show 1 --files                        # full folder listing
lw list                                  # all experiments
lw list --type cloning --status active   # filtered
```

---

## Commands

### `init`

Creates the lab directory, the `experiments/` subdirectory, the SQLite
database, and seeds the counter file. Safe to re-run: it uses
`CREATE TABLE IF NOT EXISTS` and will **not** reset an existing counter.

### `new <type> <shortname> [--title TITLE...]`

The core feature. Validates the type, allocates the next experiment id
(`LM-0001`, `LM-0002`, ...), builds the folder name, creates the folder,
and registers the database row.

- `<type>` must be one of the valid types (see below); a typo is rejected
  rather than silently producing a wrong folder.
- `<shortname>` is a short slug used in the folder name.
- `--title` is optional free text; **everything after `--title` is joined
  into the title**, so it must come last on the line. If omitted, the
  title defaults to the shortname.

**Valid types** (each maps to a short folder code):

| Type | Code | | Type | Code |
|---|---|---|---|---|
| `purification` | `pur` | | `genetics` | `gen` |
| `loading` | `load` | | `computational` | `comp` |
| `gelshift` | `gs` | | `cloning` | `clon` |
| `atpase` | `atp` | | `sequencing` | `seq` |

### `link <from> <to> <relationship>`

Records a directed provenance edge. All three arguments are required and
ordered. `<from>` and `<to>` accept either form of id (`LM-0003` or just
`3`). The edge is directed: `link 1 2` and `link 2 1` are distinct.

`<relationship>` must be one of a **closed set**, and an unknown
relationship is **blocked, not warned** -- because in a provenance tool a
typo'd relationship silently fragments the graph:

- `uses_prep_from`
- `replicate_of`
- `follow_up_to`

To add a new relationship, edit `ALLOWED_LINK_RELATIONSHIPS` in `lw.py`
(a deliberate one-line change).

### `show <id> [--files]`

Prints the experiment's row, its incoming and outgoing links, and a
summary of the files in its folder. The folder **is** the file record;
there is no files table.

- Default: a one-line file summary, e.g.
  `Files: 7 (last: 20260624_LM-0003_gel-coomassie.tiff)`.
- `--files`: the full sorted listing.

Because folder and file names front-load the date (`YYYYMMDD_...`), a
plain lexical sort is already chronological -- so "newest file" needs no
filesystem stat calls.

If the experiment **folder** is missing (you moved or deleted it by hand),
that is a normal state: `show` renders `Files: folder not found (<path>)`
rather than failing.

### `list [--type TYPE] [--status STATUS]`

Lists experiments, one per line, sorted by id. Each filter is optional;
an absent filter matches everything. The filters are flags -- a bare word
(`lw list active`) is rejected with a hint, not silently ignored.

---

## Directory layout

```
$LAB_ROOT/
|-- lab.db                                  the database (two live tables)
|-- counter.txt                             next experiment number (plain int)
|-- command_log.txt                         one audit line per successful command
|-- experiments/
|   `-- 20260624_LM-0003_pur_orc4r-prep/    one folder per experiment
|       `-- ... your files, named by hand ...
```

A folder name is `YYYYMMDD_<id>_<typecode>_<shortname>`. The date is the
**creation** date (back-dating is not currently supported -- see Known
sharp edges).

---

## Database

Two live tables:

- **`experiments`** -- one row per experiment (id, type, title, shortname,
  dates, status, folder name, notes).
- **`experiment_links`** -- the directed edges
  (`from_experiment`, `to_experiment`, `relationship`), primary key on the
  `(from, to)` pair.

> **Note on existing databases with extra tables.** If your `.db` predates
> this tool, it may contain eight additional tables from retired
> subsystems (strains, samples, assay details, staging, and so on). `lw`
> **does not read or write them, and deliberately does not drop them** --
> leaving them untouched is zero-risk and reversible. Seeing eight extra
> tables in the file is expected, not a bug. The live schema is the two
> tables above.

---

## Design notes (why it is shaped this way)

These are the load-bearing decisions; the source comments mark them too.

- **Pure transforms, effects at the edges.** Each command is a pure
  `(store, args, clock) -> effects` function that computes a description of
  what should happen and touches no IO. A single `execute_effects`
  performs the database write and stdout; filesystem effects (folder
  creation, directory listing, counter read/bump) happen at named edges.
  This is what keeps the tool testable and prevents the IO entanglement an
  earlier version suffered from.

- **The counter is bumped last.** When creating an experiment, the row and
  folder are written first and the counter is incremented last. If a run
  dies midway, this wastes an experiment number (harmless) rather than
  reusing one (which would collide on the primary key). If the counter and
  database ever disagree, `new` detects it and tells you how to fix it by
  hand.

- **Whole-store rewrite on commit.** The database is small and
  single-user, so writes rewrite the whole store in one transaction rather
  than tracking per-row statements. This is a deliberate scale bet, correct
  for a few hundred rows and one writer; it would need revisiting for a
  larger or multi-writer database.

- **`type` is not a data contract.** It exists only to pick the folder
  code and to filter `list`. The typed data it once pointed at now lives in
  Excel.

---

## Known sharp edges

Deliberately out of scope for a single-user personal tool. Listed so they
don't surprise you.

- **No input sanitizing on `shortname` / `title`.** Whatever you pass goes
  into the folder name and the database. A slash or other
  filesystem-special character in a shortname can produce a surprising
  folder path. Keep shortnames to simple slugs.

- **Folder/row/counter are not transactionally coupled.** The fixed
  "row and folder first, counter last" ordering is the safety mechanism,
  not a rollback. A process killed at exactly the wrong moment can leave a
  row without its folder, or a wasted counter number. Both are
  notice-and-hand-fix events, not data loss; `new` will tell you if the
  counter and database disagree.

- **No concurrent-writer protection.** Two `lw` processes writing at once
  can clobber each other (last write wins on the whole-store rewrite).
  This tool assumes one person running one command at a time.

- **Back-dating is not supported.** The folder date and `date_started` are
  always today's date. If you need to record an experiment that started
  earlier, edit the row by hand.

---

## Running the tests

```sh
python3 test_lw.py
```

29 tests, no framework, each in a fresh temporary lab directory. Every
command has a happy path, an error path, and a no-effect assertion (a
rejected command must leave the database, counter, and folders untouched).
`link`'s four rejection paths and the counter-corruption / counter-DB
-disagreement cases are covered explicitly.
