# How to measure a gel: the live pipeline, end to end

Onboarding operator guide. Assumes no prior knowledge of this repo. Tells you what
to click and what to type to go from a scanned gel image to per-sample numbers and,
across repeats, an aggregate. A little "why" is sprinkled in where it stops you
making a mistake; the full reasoning lives in `DESIGN.md`, the schema in
`sample_sheet_schema.md`, and the extract/picker detail in
`gel_extract_cli_README.md`.

This replaces the old `PROTOCOL.md` (now in `archive/`), which described the
retired automated pipeline (`prepare_gel_image` / `validate_gel_image` /
`measure_gel` / the R scripts). Those are archived. If a habit or note tells you to
run a "prepare" step or a "measure" step, that habit is from the old path; the live
path below has neither.

ASCII only, to match the repo.

===============================================================================
Ground rules (read once)
===============================================================================

- Run every command from `images/`:
      cd ~/personal_repos/lab_utils/images
- Quote every path in single or double quotes. Real gel filenames contain spaces,
  commas, and square brackets (e.g. `20260818 rotated, LM-0008 [s0012].tif`).
  Unquoted, bash splits them and the command fails or, worse, half-works.
- The scripts are run with `uv run`, which reads each script's own dependency block.
  Do not run them with a bare `python`; the deps will be missing.
- The unit of work is a gel ANALYSIS DIRECTORY named `<image_stem>_gel_analysis`.
  The ImageJ export macro creates it next to the image. Everything for one gel
  lives inside it; the scripts take the directory (or any file in it) as their
  argument and find the standard files by name.

===============================================================================
The whole path at a glance
===============================================================================

```
ImageJ (Windows):
  draw lane ROIs  ->  run export macro  ->  <stem>_gel_analysis/ created, with:
                                              manual_lane_profiles.csv  (the signal)
                                              lane_rois.zip             (the ROIs)
                                              export provenance
  author sample sheet in Excel -> export sample_sheet.csv INTO that folder

WSL (from images/), per gel:
  validate_sample_sheet.py   gate the sheet against the profiles
  analyze_gel.py             detrend, detect bands, measure -> band_measurements.csv
  extract_lane_values.py     (or the picker) -> one number per lane
  plot_single_experiment.py  blank-correct -> single_experiment_<selector>.csv + png

WSL (from images/), once across repeats:
  validate_manifest.py       gate the manifest that lists the repeat gels
  aggregate_repeats.py       stack the repeats -> means, sd, cv, plot
```

===============================================================================
STEP 1 - ImageJ: draw lanes and export  (Windows side)
===============================================================================

1. Open the gel image in ImageJ/Fiji. Save it to disk first if it is not already;
   the export macro writes next to the image and refuses to run on an unsaved image.
2. Draw one rectangular ROI per lane, left to right, and add each to the ROI
   Manager. For evenly-spaced lanes, the `Distribute Selection.ijm` macro tiles a
   copies of one ROI at a fixed spacing (see `distribute_selection.ijm` and its
   header).
3. Run `export_lane_profiles.ijm`. It asks for the comb well count, then:
   - creates `<image_stem>_gel_analysis/` next to the image (this is the folder
     everything else uses; there is no separate "prepare" step),
   - writes `manual_lane_profiles.csv` (the per-lane migration profile; `raw_value`
     is the width-summed signal down each lane),
   - saves `lane_rois.zip` (the exact rectangles, for replay),
   - writes an export provenance file.

Why lane_index is not biological: the macro numbers lanes by image position, left
to right. Which sample sits in which lane is recorded by YOU in the sample sheet
(next step), via `well_number`. Keep them straight; the validator reports the
relationship but never guesses a flip.

===============================================================================
STEP 2 - Author the sample sheet  (Excel, exported into the folder)
===============================================================================

Fill one row per lane and export as `sample_sheet.csv` INTO the
`<stem>_gel_analysis/` folder. The columns, the vocabulary, and the not_applicable
vs blank rule are in `sample_sheet_schema.md` -- read it once; it is the source of
truth. The essentials:

- `lane_index`   image position 1..comb, matches the profile CSV. You do not choose
                 it; it is the image coordinate.
- `well_number`  YOUR load order, a permutation of 1..comb. This is the stable
                 identity of a sample.
- `sample_label` the sample's name; this is what the aggregator joins repeats on, so
                 the SAME sample must carry the SAME label across gels.
- non-sample lanes (ladders, empties) use the `not_applicable` sentinel, never a
  blank. A blank is treated as a mistake; `not_applicable` means "deliberately
  nothing here." (The validator hard-fails a blank where a sentinel is required.)

===============================================================================
STEP 3 - Cross to WSL and run the per-gel stages
===============================================================================

The clicking happened in Windows; the scripts run in WSL. The friction is getting
the Windows paths across and looping the stages over a session's worth of gels. The
session driver does both.

Open the driver next to this doc. In nvim, with this file open:

    :vsplit session_template.sh      or put the cursor on the name below and press gf
                                     session_template.sh

(For `gf` to resolve the name, nvim's working directory should be `images/`, or add
`:set path+=.` and `:set suffixesadd=.sh`. Two panes: this doc on the left, the
driver on the right.)

The driver is a TEMPLATE. Copy it to an active, gitignored working copy so the
template stays clean and your pasted paths never get committed:

    cp session_template.sh session.sh      # session.sh is gitignored
    nvim session.sh

In `session.sh`, paste your gel paths into `GEL_PASTED_PATHS`, one per line. In
Windows File Explorer, right-click the gel image (or its `_gel_analysis` folder) and
"Copy as path"; paste verbatim (the driver strips the quotes and runs `wslpath` for
you). Set `REGION_START_MILLIMETRES` / `REGION_END_MILLIMETRES` to the migration
window you are quantifying. Then:

    bash session.sh                  # runs the readiness report + per-gel stages

What it does, and why it is shaped this way:
- Resolves every pasted Windows path to its `_gel_analysis` directory (image path or
  folder path both work).
- Prints a readiness report: for each gel, whether the profile CSV, the ROI zip, and
  the sample sheet exist. It proceeds on the gels that have all three and reports the
  rest rather than aborting -- a half-prepared gel does not block your ready ones.
- Runs validate_sample_sheet -> analyze_gel -> extract_lane_values -> plot for each
  ready gel. A real failure on one gel (a bad sheet, a failed analyze) skips that
  gel's downstream stages but continues the session.

Existence is all the driver checks; schema validation stays in the scripts, on
purpose, so there is one source of truth for what a valid input is.

Running a single gel by hand (no driver) is of course fine:

    uv run validate_sample_sheet.py '<stem>_gel_analysis'
    uv run analyze_gel.py           '<stem>_gel_analysis'
    uv run extract_lane_values.py   '<stem>_gel_analysis' --region 31.3 46.1
    uv run plot_single_experiment.py '<stem>_gel_analysis' \
        --extract-csv '<stem>_gel_analysis/extract_region_31.3-46.1mm_none.csv'

Choosing the number: `extract_lane_values.py` has two modes. `--region A B`
integrates a migration window you choose (millimetres); `--band N` takes consensus
band N from `band_measurements.csv`. If you would rather see the gel and drag the
window, use the picker instead of the CLI:

    uv run serve_gel_picker.py '<stem>_gel_analysis'
    # open the printed localhost URL, drag a span or pick a band, Write CSV, Close

The picker shells `extract_lane_values.py` on Write CSV, so the browser value is
byte-for-byte the CLI value; they cannot drift. Detail in
`gel_extract_cli_README.md`.

===============================================================================
STEP 4 - Aggregate across repeats
===============================================================================

Once each repeat gel has its `single_experiment_<selector>.csv` (the plot step's
output), stack them.

1. Author `manifest.csv` (one row per repeat gel; see `manifest_template.csv` for
   the columns). Each row's `analysis_path` points at a gel's `_gel_analysis`
   directory; each `gel_id` must be unique.
2. Validate and aggregate:

       uv run validate_manifest.py '<dir with manifest.csv>'
       uv run aggregate_repeats.py '<dir with manifest.csv>' \
           --selector 'region_31.3-46.1mm_none'

   `--selector` picks which `single_experiment_<selector>.csv` to read when a gel
   directory has more than one. `--normalize-to <sample_label>` divides each gel's
   values by a reference sample before averaging, if you want ratios rather than raw.

CAVEAT you will hit with real repeats: the aggregator currently requires every gel
to share the SAME selector, i.e. the same absolute migration window in millimetres.
Two repeats of the same gel imaged at different positions on the platen migrate to
slightly different absolute mm ranges, so their windows differ and the aggregator
refuses the join. This is a known limitation, deferred to the region-definition
design decision; the intended fix is to require equal WIDTH and warn (not fail) on a
different offset. Until then, if two repeats will not aggregate, check whether their
region windows differ -- and log it in the friction section below, because your real
repeats are exactly the data that will settle how this should behave.

===============================================================================
Friction / deviations log  (append as you go)
===============================================================================

The point of this section: as you quantify real gels, write down anything awkward,
wrong, or surprising, keyed to the step number above. This is the raw material for
the next round of refinements, and it is how deferred decisions (the region caveat,
measurement-selection ergonomics, plotter styling) get settled by real use instead
of guessing. One line each is fine.

- (step N) what happened / what was awkward / what you expected instead
-
