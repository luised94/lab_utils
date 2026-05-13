# Workflow Recipes


## Recipe 1: Purification Experiment End-to-End

Goal: Set up, execute, and close a protein purification from culture
to frozen aliquots.

### Register strains (if not already done)

  lw strain add LY101 "MATa ura3 leu2 trp1 his3"
  lw strain update LY101 label wild-type
  lw strain add LY456 "MATa ORC1-TEV-3xFLAG::URA3 orc4-R267A::KanMX"
  lw strain update LY456 label orc4-R267A

Check existing strains:

  lw strain list

### Create the experiment

  lw exp init purification orc-r267a --title "Orc4-R267A ORC purification"

Output:
  Created: LM-0005
  Folder:  20260513_LM-0005_pur_orc-r267a/
  Type:    purification
  Title:   Orc4-R267A ORC purification
  Next: lw exp update LM-0005 protein <name>

### Add metadata

  lw exp update LM-0005 protein ORC-R267A
  lw exp addstrain LM-0005 LY456 expression-strain
  lw exp update LM-0005 notes Started 6L culture in YPR at OD 0.4

Later, as the purification progresses:

  lw exp update LM-0005 induction_od 1.2
  lw exp update LM-0005 galactose_hours 4.5
  lw exp update LM-0005 columns anti-FLAG -> HiTrap Heparin
  lw exp update LM-0005 notes Induced at OD 1.2, 2% galactose
  lw exp update LM-0005 notes Lysed in buffer H-150, 3 passes French press

### Pre-register expected images

Before going to the instrument:

  lw stage expect LM-0005 gel-coomassie gel-silver

Output:
  s1  gel-coomassie
  s2  gel-silver
  At instrument:
    lemr LM-0005 s1
    lemr LM-0005 s2
  2 total pending for LM-0005

### Image at instrument

At the Imager 680 or Typhoon, name each acquisition using the
convention printed by stage expect:

  lemr LM-0005 s1     (the coomassie gel)
  lemr LM-0005 s2     (the silver stain)

The instrument appends its own suffix (timestamp, channel).

### Transfer and file images

Copy files from instrument to staging/ (drag-and-drop or USB).
For Imager 680 ZIPs, unzip first.

Check what's in staging:

  lw stage list

Dry-run the auto-assignment:

  lw stage auto

Review the output. If correct:

  lw stage auto --confirm

### Record results and close

  lw exp update LM-0005 fractions_pooled 12-18
  lw exp update LM-0005 yield_mg 2.4
  lw exp update LM-0005 storage -80C box 3 positions A1-A6
  lw exp update LM-0005 notes Good yield, clean gel, minor contaminant at 45 kDa
  lw exp complete LM-0005

Output:
  Completed: LM-0005
    Date: 2026-05-15
    Duration: 2d (started 2026-05-13)

### Review

  lw exp show LM-0005

Nvim shortcut: <leader>lx while editing a file in the experiment
folder to telescope through all files in that experiment.


## Recipe 2: Loading Assay with Manifest

Goal: Set up a factorial loading assay using a design.py file to
generate a sample manifest.

### Prerequisites

- Upstream purification exists (e.g. LM-0005 ORC-R267A, LM-0003 ORC-WT)
- All strains registered with labels (check: lw strain list --no-label)

### Create the experiment

  lw exp init loading mcm-load-orc-comparison --title "MCM loading: WT vs R267A ORC"

### Link to upstream purification

  lw exp link LM-0008 LM-0005 uses_prep_from
  lw exp link LM-0008 LM-0003 uses_prep_from

This auto-populates assay_details.protein_prep for the first link.

### Write design.py

Copy the template:

  cp docs/design_template.py experiments/20260515_LM-0008_load_mcm-load-orc-comparison/design.py

Edit the design file (see docs/design_template.py for full annotation):

  DESIGN = {
      "experiment_id": "LM-0008",
      "expected_samples": 14,
      "categories": {
          "strain": ["LY101", "LY456", None],
          "dna": ["pARS1", "pARS1-B2mut"],
          "salt": ["100mM", "300mM"],
      },
      "exclude": [
          lambda row: row["strain"] is None and row["salt"] == "300mM",
      ],
      "extras": [
          {"strain": "LY101", "label": "input", "dna": "pARS1"},
      ],
      "sort_order": ["strain", "dna", "salt"],
  }

### Generate the manifest

  lw exp manifest LM-0008

Output shows: strain resolution, sample count verification, summary
table, and confirms manifest.csv was written.

If you need to regenerate after editing design.py:

  lw exp manifest LM-0008 --force

### Add strains to experiment

  lw exp addstrain LM-0008 LY101
  lw exp addstrain LM-0008 LY456

### Proceed with staging

Pre-register expected images as in Recipe 1:

  lw stage expect LM-0008 gel-load gel-input

Image, transfer, and file using the same staging workflow.

### Nvim shortcuts during setup

  <leader>le    Insert an experiment ID while writing notes
  <leader>ls    Insert a strain ID while writing design.py
  <leader>lp    Find a protocol file for reference


## Recipe 3: Staging -- Pre-Registered Workflow

Goal: Use the three-phase staging model (plan, image, file) for
systematic file management.

### Phase 1: Plan

Before going to the instrument, register what you expect to image:

  lw stage expect LM-0005 gel-coomassie gel-silver bradford-raw

Output:
  s1  gel-coomassie
  s2  gel-silver
  s3  bradford-raw
  At instrument:
    lemr LM-0005 s1
    lemr LM-0005 s2
    lemr LM-0005 s3
  3 total pending for LM-0005

Check all pending expectations:

  lw stage expect --list
  lw stage expect --list LM-0005

### Phase 2: Image

At every instrument, use the same convention:

  lemr LM-0005 s1

The instrument appends its own metadata. The important part is that
LM-0005 and s1 appear as tokens in the filename.

### Phase 3: File

Copy all files to staging/. Then:

  lw stage list          # verify files are detected and parsed
  lw stage auto          # dry run -- review matches
  lw stage auto --confirm

If a file has the experiment ID but no slot, and that experiment has
only one pending expectation, the slot is inferred (shown as
"slot inferred" in dry-run output).

### Canceling expectations

Changed plans and won't image something:

  lw stage expect --cancel LM-0005 s3     # cancel one slot
  lw stage expect --cancel LM-0005        # cancel all pending

### Retakes

Not system-enforced. Append -v2 to descriptor:

  lw stage expect LM-0005 gel-coomassie-v2

Both versions are kept in the experiment folder.


## Recipe 4: Staging -- Manual (Ad Hoc)

Goal: File a one-off image without pre-registration.

For files that weren't planned ahead of time (a quick gel photo,
a screenshot, exported data):

  lw stage assign random_photo.tif LM-0005 gel-quick-check

Output:
  Experiment: LM-0005 (orc-r267a, purification)
  Filed: random_photo.tif
      -> experiments/20260513_LM-0005_pur_orc-r267a/20260515_LM-0005_gel-quick-check.tif

The file is renamed with today's date, moved to the experiment folder,
and registered in the files table.

Manual assign coexists fully with the pre-registration workflow. Use
whichever fits.


## Recipe 5: Linking Experiments

Goal: Connect related experiments to build a provenance chain.

### Purification -> Downstream assay

  lw exp link LM-0008 LM-0005 uses_prep_from

This means: LM-0008 (loading assay) uses protein from LM-0005
(purification). The assay_details.protein_prep field is auto-populated.

### Replicates

  lw exp link LM-0010 LM-0008 replicate_of

### Follow-ups

  lw exp link LM-0012 LM-0008 follow_up_to

### Viewing links

  lw exp show LM-0008

Links section shows:
  -> LM-0005 (uses_prep_from)
  <- LM-0010 (replicate_of)
  <- LM-0012 (follow_up_to)

### Custom relationships

Non-standard relationships are accepted with a note:

  lw exp link LM-0015 LM-0005 troubleshoots

Output includes:
  Note: 'troubleshoots' is not a standard relationship
        (uses_prep_from, replicate_of, follow_up_to)


## Recipe 6: Figure Provenance (Manual)

Goal: Track which experiments contribute to which publication figures
until fig commands are implemented.

The figure_panels table exists in the schema but has no CLI commands
yet. Track provenance manually in lab_notes.md:

  ## 2026-05-20

  Figure 3A: LM-0005 gel-coomassie (Coomassie stain of purification)
  Figure 3B: LM-0005 gel-silver (Silver stain, same gel)
  Figure 4A: LM-0008 gel-load (MCM loading, WT vs R267A)
  Source AI file: publications/plos-one-2026/figures/figure-3.ai

Or insert directly via SQL if you want the database record now:

  lwdb
  INSERT INTO figure_panels (publication, figure, panel, experiment_id,
    source_file, ai_file)
  VALUES ('plos-one-2026', 'Figure 3', 'A', 'LM-0005',
    'experiments/20260513_LM-0005_pur_orc-r267a/20260515_LM-0005_gel-coomassie.tif',
    'publications/plos-one-2026/figures/figure-3.ai');

Query provenance:

  lwdb "SELECT figure, panel, experiment_id, source_file
        FROM figure_panels
        WHERE publication = 'plos-one-2026'
        ORDER BY figure, panel"


## Recipe 7: Common Queries

### CLI queries

Find all purification experiments:

  lw exp find --type purification

Find active experiments using a specific strain:

  lw exp find --strain LY456 --status active

Search by keyword (searches id, shortname, title, notes):

  lw exp find mcm loading

List all experiments grouped by type:

  lw exp list

List only loading experiments:

  lw exp list --type loading

Strains without labels (pre-manifest QC):

  lw strain list --no-label

### Raw SQL queries (via lwdb)

Files registered to an experiment:

  lwdb "SELECT filename FROM files WHERE experiment_id = 'LM-0005'"

All experiments linked to a purification:

  lwdb "SELECT from_experiment, relationship
        FROM experiment_links
        WHERE to_experiment = 'LM-0005'"

Experiments with no files:

  lwdb "SELECT id, shortname FROM experiments
        WHERE status = 'active'
        AND id NOT IN (SELECT DISTINCT experiment_id FROM files)"

Purification yields:

  lwdb "SELECT e.id, p.protein, p.yield_mg
        FROM experiments e
        JOIN purification_details p ON e.id = p.experiment_id
        ORDER BY e.date_started DESC"

Sample manifest for an experiment:

  lwdb "SELECT position, strain_id, label, dna, conditions
        FROM samples
        WHERE experiment_id = 'LM-0008'
        ORDER BY position"

Strain usage across all experiments:

  lwdb "SELECT s.id, s.label, COUNT(es.experiment_id) as n_exps
        FROM strains s
        LEFT JOIN experiment_strains es ON s.id = es.strain_id
        GROUP BY s.id
        ORDER BY n_exps DESC"


## Recipe 8: Daily Workflow (Neovim + CLI)

Goal: Typical daily pattern combining nvim keymaps and CLI.

### Morning: review and plan

  lw status                    # check what needs attention
  lwnotes                     # open notes, review yesterday
  <leader>lh                  # add today's date header

Write today's plan in lab_notes.md. Reference experiments inline
using <leader>le to insert IDs directly from the database.

### During experiments

While editing lab_notes.md or any file:

  <leader>le     Insert experiment ID (pick from sorted list)
  <leader>ls     Insert strain ID
  <leader>lp     Find and open a protocol

While working in an experiment folder:

  <leader>lx     Telescope through files in this experiment

### After imaging

  lw stage list              # see what landed in staging
  lw stage auto              # dry run
  lw stage auto --confirm    # file everything

### End of day

  lw exp update LM-0005 notes Completed column chromatography, fractions on ice
  <leader>lt                  # review TODOs in lab_notes.md
  lwbackup                    # back up everything


## Recipe 9: Backup

Goal: Maintain consistent, recoverable backups of the lab directory.

### Regular backup

  lwbackup

This runs lw_backup.sh which:
  1. Snapshots lab.db safely (sqlite3 .backup, not raw file copy)
  2. Syncs all files via rsync (additive, never deletes from backup)
  3. Creates weekly DB snapshots (lab_YYYY-WNN.db) automatically
  4. Prunes weekly snapshots older than 84 days

### Dry-run (see what would happen)

  lwbackup --dry-run

### What gets backed up

Everything in lab/ except lab.db itself (replaced by consistent
snapshot). The rsync is additive: files deleted from lab/ are
preserved in the backup.

### Recovery

To restore from backup:
  1. Copy the backup directory contents to a new lab/ location
  2. lab.db in the backup is a consistent snapshot, ready to use
  3. For point-in-time recovery, weekly snapshots are in
     snapshots/lab_YYYY-WNN.db

### When to back up

After any session that creates or modifies data. The lwbackup alias
makes it a one-command habit. Consider adding it to your end-of-day
routine (see Recipe 8).
