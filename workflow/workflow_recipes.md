# Workflow Recipes

Step-by-step guides for common workflows. Each recipe lists the goal,
prerequisites, commands, and expected output.

All examples use realistic experiment IDs and strain names from ORC/MCM
biochemistry. Substitute your own values.


## Recipe 1: Set Up a Purification Experiment

Goal: Create a purification experiment, register strains, and record
initial metadata.

Prerequisites: System initialized (`lw init` has been run).

### Steps

1. Create the experiment:

        lw exp init purification orc4-r267a-prep --title "Orc4-R267A ORC purification"

   Output:

        Created: LM-0001
        Folder:  20260513_LM-0001_pur_orc4-r267a-prep/
        Type:    purification
        Title:   Orc4-R267A ORC purification
        Next: lw exp update LM-0001 protein <name>

2. Set the protein name (required for purification records):

        lw exp update LM-0001 protein ORC-R267A

3. Register the strain if it does not exist yet:

        lw strain add LY456 "MATa orc4-R267A::KanMX ura3 leu2 trp1 his3"
        lw strain update LY456 label "Orc4-R267A"

4. Link the strain to the experiment:

        lw exp addstrain LM-0001 LY456

5. Add notes as you go:

        lw exp update LM-0001 notes 1L culture, induced at OD 1.2

   Notes always append with a date stamp. Add more entries anytime:

        lw exp update LM-0001 notes calmodulin elution looks clean

6. Record purification-specific details:

        lw exp update LM-0001 columns "calmodulin, MonoQ"
        lw exp update LM-0001 yield_mg 0.45
        lw exp update LM-0001 storage "-80C box 3 slot 12"

7. When finished:

        lw exp complete LM-0001
        # Completed: LM-0001
        # Duration: 3d (started 2026-05-13)


## Recipe 2: Set Up a Loading Assay with Manifest

Goal: Create a loading assay, define a factorial design, generate a
sample manifest, and link to an upstream purification.

Prerequisites: Upstream purification exists (e.g., LM-0001). Strains
registered and labeled.

### Steps

1. Create the experiment:

        lw exp init loading mcm-loading-ars1 --title "MCM loading with ARS1, testing R267A"

2. Link to the upstream purification:

        lw exp link LM-0002 LM-0001 uses_prep_from

   This auto-populates assay_details.protein_prep.

3. Verify strains are registered and labeled:

        lw strain list --no-label

   If any strains lack labels, set them now. Labels appear in the
   manifest and are required for readability:

        lw strain update LY100 label "WT ORC"

4. Create the design file. In the experiment folder, create design.py:

        experiments/20260513_LM-0002_load_mcm-loading-ars1/design.py

   Contents (see docs/design_template.py for a full annotated template):

        DESIGN = {
            "experiment_id": "LM-0002",
            "expected_samples": 14,
            "extras": [
                {"label": "input", "strain": "LY100", "amount": "5% total"},
                {"label": "input", "strain": "LY456", "amount": "5% total"},
            ],
            "categories": {
                "strain": ["LY100", "LY456", None],
                "dna": ["ARS1", "A-B2-"],
            },
            "sort_order": ["strain", "dna"],
            "exclude": [
                lambda r: r["strain"] is None and r["dna"] == "A-B2-",
            ],
        }

   Expected count: 2 extras + (3 strains x 2 DNAs - 1 excluded) = 2 + 5 = 7.
   Wait -- that does not add up. Let me recalculate:
   3 strains x 2 DNAs = 6. Minus 1 excluded (None + A-B2-) = 5. Plus
   2 extras = 7. But expected_samples says 14. Adjust your categories
   and expected count to match your actual design.

5. Generate the manifest:

        lw exp manifest LM-0002

   Output shows strain label resolution, sample table, and writes
   manifest.csv to the experiment folder.

   If you need to regenerate after fixing the design:

        lw exp manifest LM-0002 --force

6. Add strains to the experiment:

        lw exp addstrain LM-0002 LY100
        lw exp addstrain LM-0002 LY456

7. Verify everything:

        lw exp show LM-0002


## Recipe 3: Process Instrument Files Through Staging

Goal: Get files from an instrument into the correct experiment folder
with proper names, using the pre-registration workflow.

Prerequisites: Experiment exists.

### Steps

1. Before going to the instrument, register what you plan to image:

        lw stage expect LM-0002 gel-coomassie gel-silver membrane-anti-mcm

   Output:

        s1  gel-coomassie
        s2  gel-silver
        s3  membrane-anti-mcm
        At instrument:
          lemr LM-0002 s1
          lemr LM-0002 s2
          lemr LM-0002 s3
        3 total pending for LM-0002

2. At the instrument, name each acquisition using the printed names.
   For the Imager 680, type `lemr LM-0002 s1` as the file prefix --
   the instrument appends its own timestamp and channel suffix.

3. Transfer files to your computer. For the Imager 680, unzip the
   downloaded archive first.

4. Copy (or move) all files into the staging/ directory.

5. Preview what will happen:

        lw stage auto

   Output shows matched files, their destinations, any skipped or
   unmatched files. Files matched by inference (experiment ID found
   but no slot number, with only one pending expectation) are annotated
   with "(slot inferred)".

6. Execute:

        lw stage auto --confirm

   Output:

        Filed (3):
          lemr LM-0002 s1 2026.05.13_14.30.00_Cy5.tif
             -> 20260513_LM-0002_gel-coomassie.tif
          lemr LM-0002 s2 2026.05.13_14.35.00_Cy5.tif
             -> 20260513_LM-0002_gel-silver.tif
          lemr LM-0002 s3 2026.05.13_14.40.00_Cy5.tif
             -> 20260513_LM-0002_membrane-anti-mcm.tif

### Manual assignment (ad hoc files)

If you have a file that was not pre-registered:

    lw stage assign "some_scan.tif" LM-0002 supplemental-scan

The system confirms the experiment identity before filing:

    Experiment: LM-0002 (mcm-loading-ars1, loading)
    Filed: some_scan.tif
        -> experiments/20260513_LM-0002_load_mcm-loading-ars1/20260513_LM-0002_supplemental-scan.tif

### Managing expectations

Check pending expectations:

    lw stage expect --list
    lw stage expect --list LM-0002

Cancel expectations you no longer need:

    lw stage expect --cancel LM-0002 s3    # cancel one slot
    lw stage expect --cancel LM-0002       # cancel all pending


## Recipe 4: Link a Downstream Assay to an Upstream Purification

Goal: Connect an assay experiment to the purification that produced
the protein it uses.

Prerequisites: Both experiments exist.

### Steps

1. Link them:

        lw exp link LM-0005 LM-0001 uses_prep_from

   This does two things:
   - Creates a record in experiment_links
   - Auto-populates assay_details.protein_prep for LM-0005

   If LM-0001 is not a purification experiment, the system warns but
   creates the link anyway.

2. Verify:

        lw exp show LM-0005

   Under "Links:" you will see:

        -> LM-0001 (uses_prep_from)

   Under "Assay details:" you will see:

        protein_prep: LM-0001

3. Other relationship types:

        lw exp link LM-0006 LM-0005 replicate_of
        lw exp link LM-0007 LM-0005 follow_up_to

   Non-standard relationships are allowed with a warning.


## Recipe 5: Track Figure Provenance (Manual Workflow)

Goal: Record which experiment data appears in which publication figure.
The fig link and fig trace commands are not yet implemented, so this
uses manual tracking in lab_notes.md and SQL.

Prerequisites: Experiments with processed images. Publication folder
exists.

### Steps

1. In lab_notes.md, record figure panel sources:

        ## 2026-05-13

        Figure 3A: LM-0001, gel-coomassie (Coomassie stain of purified ORC)
        Figure 3B: LM-0002, membrane-anti-mcm (MCM loading Western)
        Figure 3C: LM-0005, plot-quantification (loading quantification)

2. For formal tracking, insert directly into the database:

        sqlite3 lab.db "INSERT INTO figure_panels
            (publication, figure, panel, description, experiment_id, source_file)
            VALUES (
                'orc4-r267a-paper',
                'Figure 3',
                'A',
                'Coomassie stain of purified ORC',
                'LM-0001',
                'experiments/20260513_LM-0001_pur_orc4-r267a-prep/20260513_LM-0001_gel-coomassie.tif'
            );"

3. Query provenance:

        sqlite3 lab.db "SELECT figure, panel, experiment_id, source_file
            FROM figure_panels
            WHERE publication = 'orc4-r267a-paper'
            ORDER BY figure, panel;"

When fig link and fig trace commands are implemented, they will replace
this manual workflow.


## Recipe 6: Common Queries

### Find all experiments using a specific strain

    lw exp find --strain LY456

### Find all active purification experiments

    lw exp find --type purification --status active

### Free-text search across experiments

    lw exp find orc4

Searches ID, shortname, title, and notes.

### List all purifications with yields

    lw exp list --type purification

Shows protein name and yield for each purification.

### SQL: Find all samples using a specific DNA template

    sqlite3 lab.db "SELECT s.experiment_id, s.position, s.label, s.dna
        FROM samples s
        JOIN experiments e ON s.experiment_id = e.id
        WHERE s.dna = 'ARS1'
        ORDER BY e.date_started, s.position;"

### SQL: Count samples per experiment

    sqlite3 lab.db "SELECT experiment_id, COUNT(*) as n
        FROM samples
        GROUP BY experiment_id
        ORDER BY n DESC;"

### SQL: Find experiments by a condition stored in JSON

    sqlite3 lab.db "SELECT experiment_id, position, label
        FROM samples
        WHERE json_extract(conditions, '$.kglut') = '300';"

### SQL: Trace a figure panel back to its experiment

    sqlite3 lab.db "SELECT fp.figure, fp.panel, fp.description,
            e.type, e.shortname, e.date_started
        FROM figure_panels fp
        JOIN experiments e ON fp.experiment_id = e.id
        WHERE fp.publication = 'orc4-r267a-paper'
        ORDER BY fp.figure, fp.panel;"

### Check the command log

    tail -20 lab/command_log.txt

Or filter for a specific experiment:

    grep LM-0005 lab/command_log.txt
