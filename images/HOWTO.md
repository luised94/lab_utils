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
  author sample sheet in Excel -> export sample_sheet.csv INTO the <stem>_gel_analysis/ folder

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
2. WRONG Draw one rectangular ROI per lane, left to right, and add each to the ROI
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
   - ADD other is also an option

Why lane_index is not biological: the macro numbers lanes by image position, left
to right. Which sample sits in which lane is recorded by YOU in the sample sheet
(next step), via `well_number`. Keep them straight; the validator reports the
relationship but never guesses a flip.

===============================================================================
STEP 2 - Author the sample sheet  (Excel, exported into the folder)
===============================================================================

Fill one row per lane in the experiments' respective excel file and export as `sample_sheet.csv` INTO the
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
All tif filepaths has to be absolute paths. Use printf in an accesbile tmux pane to be able to copy paste easily.

    uv run validate_sample_sheet.py '<stem>.tif'
    uv run analyze_gel.py           '<stem>.tif'
    uv run serve_gel_picker.py      '<stem>.tif'
    uv run extract_lane_values.py   '<stem>.tif' --region MM.M MM,M # copy from the bottom of the html file while picking the region.
    uv run calculate_gel_shift_ratio.py  '<stem>.tif'
    Rscript plot_gel_shift_ratio.R  '<path-to-csv>/manifest.csv' <path-to-csv> # printf '%s\n' "$(pwd)"/**/gel_shift_ratio_*.csv to display the files to add to manifest csv
    OR
    uv run calculate_loading_value.py '<stem>.tif'
    uv run plot_loading.R 'manifest.csv'
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

REPEATS AT DIFFERENT PLATEN POSITIONS (the common real case). Two repeats of the
same gel imaged at different positions migrate to slightly different absolute mm
ranges, so the SAME band sits at a different window on each -- e.g. 31.3-46.1mm on
one, 29.8-44.6mm on another. The aggregator handles this: it requires the same
window WIDTH (end - start) and the same baseline across gels -- the proxy for the
same molecular-size range under one imager calibration -- and only WARNS on a
different offset. So extract EACH gel over its own window (same width, positioned
over the band on that gel), then aggregate:

  - If every gel used the identical window, pin `--selector` as above.
  - If gels used different offsets (same width), do NOT pin `--selector`; each gel's
    selector string differs. Omit it and let the aggregator auto-discover one
    `single_experiment_*.csv` per gel. It joins them, warns about the offset
    difference in the checks JSON, and names the output by the shared width, e.g.
    `aggregate_region_14.8mm_none_multi.csv` (the `_multi` marks an offset-spanning
    aggregate). For auto-discovery to be unambiguous, keep exactly one
    `single_experiment_*.csv` per gel directory (extract only the window you intend
    to aggregate for that gel).

If two repeats will NOT aggregate, the checks JSON says why: a hard failure on
"region extractions share one width" means the windows are genuinely different sizes
(a different size range) -- re-extract them at matching widths. A width check assumes
one imager calibration across gels; if you ever change imager settings between
repeats, that assumption breaks and the width proxy is no longer faithful (it would
have to be verified from the TIFF). Log anything surprising in the friction section.
