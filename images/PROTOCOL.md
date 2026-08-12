# Measurement protocol - gel densitometry pipeline

Operator instructions. Companion to `DESIGN.md`, which records why; this records
what to type. Self-contained: everything here can be followed without reading the
design.

State as of 2026-08-11. This describes the committed pipeline (stage 0 through the
replicate aggregator). Record any issue, deviation, or friction you hit against the
section number here; the "Deviations log" stub at the end is the place for it.

Run every command from `~/personal_repos/lab_utils/images`.

**Quote every path with single quotes.** Real filenames contain spaces, commas,
square brackets, and characters like `& ; | ( ) $ ` `` ` `` `! * ? # ~`. Single
quotes neutralize all of them. Note that inside *double* quotes an interactive
bash still history-expands `!`, which is one more reason to use single quotes.

The whole path, at a glance:

```
Fiji (manual): landmarks + crop  -> <stem>_preprocess.txt
  stage 0  prepare_gel_image.py     -> <stem>_gel_analysis/ (created)
  stage 1  validate_gel_image.py    -> input_file_validation_report.json (+ preview)
  stage 2  measure_gel.py           -> band_measurements.csv, band_detection_report.json, overlays
  R (per gel)  gel_loading_normalization.R  -> loading_normalization_per_lane.csv (+ per-band, plot)
  R (across)   aggregate_loading_replicates.R -> summary, matrices, aggregate_bar.pdf
  Illustrator: open the per-gel plot and aggregate_bar.pdf (vector, editable text)
```

Two R stages are separate on purpose: the first normalizes one gel against its own
reference lane; the second combines replicates. The sensitivity harness and the
statistics reference are OPTIONAL and are not on this path (appendices A and B).

**Every number is PROVISIONAL until `encoding_verified` is true.** For the Typhoon
`.tif` files linearity is unresolved (section 6 and section 12). The pipeline
carries the flag forward into every report and every CSV; it never blocks you, and
it never silently upgrades. Until it is true, treat ratios as provisional.

---

## 0. The two things that will catch you out

The image is mirrored top-to-bottom on screen. The scan origin is bottom-left, so
file row 0 is the physical bottom of the gel, and Fiji draws row 0 at the top of
the window. You do not correct for this: measure the image exactly as it opens,
and the pipeline flips it later.

Identify the orientation by eye before measuring. Find the wells (the loading
points); migration runs away from them. If the wells run down the left or right
edge and bands march sideways, the migration axis is horizontal; if the wells run
across the top or bottom and bands march up or down, it is vertical. Record this
in `gel_migration_axis`. Do not assume: gels get imaged both ways, and the wells
are only faintly visible, so if the thick cast edge of the cassette is clearer,
use it to fix the axis and confirm with the wells.

Second, less dangerous: Fiji renders these files through an inverting lookup
table, because the TIFF declares `PhotometricInterpretation = MINISWHITE`. Bands
look dark on a light field. Every image this pipeline emits uses a matching
reversed colormap, so the two agree. Bright bands on a dark field from this
pipeline mean something is wrong with the pipeline, not your file.

---

## 1. Which files to process

**Process the `.tif` only.** Never the `.img`, `.gel`, `.png` or `.jpg`.

- The Typhoon (phosphor, radioactivity) writes four files per scan: `.gel`,
  `.img`, `.inf`, `.tif`. The `.tif` is the read path. The `.inf` is recorded as
  provenance if it is present and is not required.
- The Amersham Imager 680 (fluorescence, loading) writes a `.tif` plus a `.png`
  and a `.jpg`. Only the `.tif` is 16-bit. The `.png` and `.jpg` are 8-bit
  display renderings and measuring one would be silent nonsense.

Both instruments are handled by the same scripts with no flags, because
everything that differs is read from the file's own tags. Pixel size is 200
micrometres on the Typhoon and 79.9 on the Imager 680, and nothing anywhere
assumes either.

To list the `.tif` files in an experiment directory with full, space-safe paths:

```
printf '%s\n' "$(pwd)"/*.tif
```

If no `.tif` matches, bash prints the literal pattern `.../*.tif` rather than
nothing. The printed paths still need single-quoting when pasted back.

---

## 2. Inventory the tags, once per new instrument or new file type

Reports everything, interprets nothing, validates nothing. Run it when you meet a
file type for the first time, or when something looks wrong.

```
uv run inventory_tiff_tags.py '<experiment directory>/<scan>.tif' > tags.txt
```

stdout is the inventory and stderr is the commentary, so this redirects cleanly.

---

## 3. Stage 0 - create the output directory

```
uv run prepare_gel_image.py '<path>.tif'
```

Prints resolved paths and exits, reading no image data. It creates
`<stem>_gel_analysis` beside the input and locates the scanner `.inf` sidecar by
stem if present (`inf_sidecar_is_present` is reported either way; its absence is
fine).

The output directory is **per scan, not per file**: the `.img` and the `.tif` of
one scan share a stem and resolve to the same directory. If you already ran this
against the `.img`, the directory exists and you will see a not-empty warning. That
is expected. The parent directory must already exist; only the leaf is created, so
a typo in `--output-parent-directory` fails instead of building a tree somewhere
wrong.

---

## 4. Measure the sidecar in Fiji

Once per scan, the only manual step, and the one the pipeline cannot do for you:
the crop rectangle and the two landmarks cannot be derived from pixels. Copy the
template to `<stem>_preprocess.txt` next to the `.tif` first. Open the `.tif` in
Fiji (drag and drop). Do NOT adjust brightness/contrast and Apply, convert to
8-bit, rotate, or save; any of those destroys data or discards tags.

Do these in order. The line macro in step 4 only works while a line is the active
selection, so the rectangle must come after it.

1. **Preflight, once per session.** `Analyze > Set Measurements`: tick "Bounding
   rectangle", and confirm "Invert Y coordinates" is UNTICKED. If it is ticked,
   every y you record is silently mirrored. Confirm the origin is 0,0.

2. **Record identity.** Fill `measured_against_input_filename` (filename only,
   and it must match the on-disk `.tif` name EXACTLY - see the caution below),
   `measured_against_image_width_pixels` and `_height_pixels`,
   `gel_migration_axis` (from section 0), and `coordinate_unit`.

3. **Draw the landmark line.** Use the well line if you can see it, else the
   cassette cast edge. With the straight-line tool, draw from one far end to the
   other, as far apart as possible. Do NOT hold Shift: it snaps the line straight
   and reports a tilt of exactly zero that looks like success but has destroyed
   the measurement. Zoom to 100% at each end before placing the point.

4. **Measure and capture the line.** `Ctrl+M`, record the `Angle`. Then, with the
   line still selected, `Plugins > New > Macro`, paste and run (`Ctrl+R`):

   ```
   getLine(x1, y1, x2, y2, w); print(x1+", "+y1+", "+x2+", "+y2);
   ```

   Record the endpoints as `landmark_a` and `landmark_b`, and note in
   `rotation_landmark_description` what you clicked and which end is a.

5. **Draw the crop rectangle and `Ctrl+M`.** Contain every band plus 20 to 50 px
   of clean gel on all four sides; do not clip a band or its flanking gel. Glance
   at the reported Max: below 65535 means no ceiling-clipped pixel is inside the
   crop. Record `crop_x`, `crop_y`, `crop_width`, `crop_height`.

6. **Finish.** Fill `expected_lane_count` (empty lanes included), save.

> **Caution that will cost you a run.** `measured_against_input_filename` must be
> byte-for-byte the on-disk `.tif` name. A scan whose timestamp is dotted in one
> place and underscored in another (`2026.07.14` vs `2026_07_14`) trips stage 1's
> exact-name hard stop. Fix the sidecar to match the file on disk; do not rename
> the file.

---

## 5. Stage 1 - validate

```
uv run validate_gel_image.py '<path>.tif'
```

Writes into `<stem>_gel_analysis`:

- `input_file_validation_report.json` - every check, every tag, every statistic,
  and the `preprocess_sidecar.geometry_pixels` block the measurement stage reads.
- `input_file_histogram.png` - log-count y axis.
- `input_file_crop_preview.png` - only if the sidecar exists.

Exit code is 0 if every check passed, 1 if any hard stop failed. The report is
written either way.

**Look at `input_file_crop_preview.png` before going further.** It draws your crop
box and your two landmarks onto a display-scaled copy. This is the only check on
numbers you typed by hand, and a transposed digit is obvious in the picture and
invisible in the text.

---

## 6. What stage 1 will not tell you

**Whether the pixel values are linear in signal.** For the Typhoon files this is
unresolved. The report carries `encoding_verified: false` and names the test that
would change it. If the `.tif` is not linear, every ratio this pipeline produces is
compressed toward 1 and nothing downstream will notice. For the Imager 680
fluorescence files the position is much stronger (a cooled CCD is linear by
construction), but those files record an exposure time and two images are only
comparable after dividing by it; the report extracts the exposure.

**Whether the gel is tilted more than you think.** Stage 1 reports the angle from
your landmarks and does not find the tilt itself. The pipeline never rotates: tilt
is measured and warned above 0.25 degrees, never applied.

**Whether a band is saturated.** Stage 1 reports whole-image ceiling statistics
only. Per-band saturation and the plateau statistic come at stage 2.

---

## 7. Stage 2 - measure

```
uv run measure_gel.py '<path>.tif'
```

Reads the validated `.tif` plus the stage 1 report, crops once to the sidecar
geometry, detects lanes and band centres, integrates, and writes into
`<stem>_gel_analysis`:

- `band_measurements.csv` - one row per band per baseline method; the join key is
  `well_index` (0-based). `reported_value` is the per-band reported quantity.
- `band_detection_report.json` - provenance, including `gel_id` (the input stem)
  and `encoding_verified` carried forward from stage 1. Schema version 3.
- overlays: `overlay_reported_values.png`, the lane-grid overlay, lane profiles.

**Look at the overlays.** The reported-values overlay is the fastest way to catch
a mis-detected lane or a window that missed a band. Read the stderr summary for the
overall status and any warnings (saturation, baseline disagreement, tilt).

Defaults are chosen to be sensible; the baseline knobs
(`--baseline-flank-search-millimetres`, default 1.0; `--rolling-ball-width-millimetres`,
default 4.0; method `--baseline-cross-check-method`, default `opening`) rarely need
changing. If you do change them, record it: the numbers move (Appendix A quantifies
how much on one gel).

> `well_index` is the physical lane position (0..14), the join key everywhere
> downstream. `loaded_lane_index` also appears in the CSV; it is not the join key.
> Do not join on it.

---

## 8. Per-gel normalization (R)

Author the sample sheet, once per gel. Copy `sample_sheet_template.csv` to
`<stem>_sample_sheet.csv` next to the `.tif`. Columns:

```
well_number, sample_label, role, prep_source, mutation, suppressor,
salt_condition_mM, normalize_on_band, notes
```

- `well_number` is 1..15, each once (maps to `well_index = well_number - 1` under
  the default `direct` flip).
- `role` is one of `reference`, `sample`, `empty`. There must be **exactly one**
  `reference` lane; it is pinned to 100 and everything else is expressed as a
  percent of it. `empty` lanes are dropped from normalization.
- The factor columns (`mutation`, `suppressor`, `salt_condition_mM`, ...) carry
  within-gel identity. Fill the one(s) your experiment varies; leaving them blank
  is allowed but see the note in section 9 about the condition key.
- The repeat relationship does NOT live here. Repeats are recorded only in the
  aggregation manifest (section 9).

Run:

```
Rscript gel_loading_normalization.R '<stem>_gel_analysis' ['<stem>_sample_sheet.csv']
```

The sample sheet is found automatically by stem if you omit the second argument.
Outputs into `<stem>_gel_analysis/loading_normalization/`:

- `loading_normalization_per_lane.csv` - one row per lane, with `gel_id` first,
  `whole_lane_total`, `percent_of_reference_whole_lane`, and the control-band
  columns. This is the file the aggregator consumes.
- the per-band CSV and a per-gel plot (vector PDF).

Every output carries `gel_id` as its first column and `encoding_verified` so a
provisional gel stays labelled provisional.

---

## 9. Aggregate replicates (R)

Author a manifest: one row per gel, listing which gels are replicates of which
experiment. Copy `manifest_template.csv` and fill:

```
experiment_id, analysis_path, gel_id, replicate, notes
```

- `experiment_id` groups replicates into one experiment; `replicate` is the label
  within it (unique per experiment). These two are the only fields no file knows -
  you author them.
- `analysis_path` is the `<stem>_gel_analysis` directory, or any file inside it.
  Paths may be relative to the manifest's own location, so keep the manifest next
  to the data and use short paths.
- `gel_id` is checked against the per-lane CSV; a mismatch is a hard stop. This is
  the safety net that stops you pointing a row at the wrong gel.

Two constants at the top of `aggregate_loading_replicates.R` are the only per-run
knobs:

- `CONDITION_KEY_COLUMNS` - the sample-sheet factor column(s) that define "the same
  condition across replicates" (for example `c("suppressor")`). It is per-run.
- `NORMALIZATION_COLUMN` - `"whole_lane"` (default) or `"control_band"`.

Run:

```
Rscript aggregate_loading_replicates.R '<manifest>.csv'
```

Outputs into `<manifest_dir>/aggregate_analysis/`:

- `combined_raw_data.csv` - every gel's lanes stacked (the read-together state;
  empties kept; pre-normalized totals preserved).
- `processed_data.csv` - empties dropped, the reported value selected,
  `condition_label` added.
- per `experiment_id`: `summary_data.csv` (mean, sd, se per condition),
  `matrix_per_replicate_long.csv`, `matrix_summary_long.csv`,
  `fold_difference_matrix_mean_wide.csv`,
  `percent_difference_matrix_mean_wide.csv`, and `aggregate_bar.pdf` (vector).

The matrices are all-against-all, computed within each replicate then averaged.
Descriptive only: mean, sd, se, and the two difference matrices. No inferential
test is run here (see Appendix B).

> **By design, not a bug.** If `CONDITION_KEY_COLUMNS` is blank on every lane (a
> fast verification run with an unfilled sheet), the aggregator warns, still writes
> the combined and processed states plus a per-well fallback summary, and skips the
> matrices and plot. One condition key set applies to every experiment in a single
> manifest run; for experiments needing different keys, run the manifest once per
> experiment. Duplicate condition within a replicate stops with a message rather
> than silently averaging.

---

## 10. Open in Illustrator

Open the per-gel plot (section 8) and `aggregate_bar.pdf` (section 9) directly in
Illustrator. Both are written with the cairo PDF device: text is real text (not
outlines and not a rasterized image), so fonts are editable and swappable, and the
vector geometry scales without loss. If a preferred font (Arial) is not installed
on the machine that ran R, the plot falls back to a generic sans family; the text
is still editable, so restyle in Illustrator as needed.

That is the end of the path to a figure.

---

## Appendix A. Sensitivity harness (optional QC)

Not on the path to a figure. It re-measures one gel under small, deliberate
perturbations of the crop rectangle and the baseline knobs, and reports how much
each lane's normalized value moves as a range.

```
uv run sensitivity_harness.py '<path>.tif' \
    --pipeline-directory . \
    --output-directory '<scratch dir>' \
    [--sample-sheet '<stem>_sample_sheet.csv']
```

Reads `lane_sensitivity_summary.csv` (per lane: nominal, min, max, range) and
`sensitivity_report.md` (headline worst-case range plus any cell that changed the
detected lane count, which is excluded from the spread and flagged separately). It
never touches the repo, the inputs, or the real analysis directory. Provisional
numbers stay provisional; it measures their stability, not their correctness.

## Appendix B. Statistics reference (optional)

`loading_statistics_reference.R` is a self-contained, simulated-data reference for
the lab's house-style inferential machinery (paired t blocked on replicate, Holm
within one named family, exact p-values, effect sizes, an nlme consistency
diagnostic). It is kept for experiments that warrant inference; it is not wired
into the pipeline and reads no pipeline file.

```
Rscript loading_statistics_reference.R ['<output directory>']
```

## Appendix C. Quick tests in Fiji

Cheap checks that each catch one specific silent failure. None alter the file.

**Saturation-location test.** `Image > Adjust > Threshold` (`Ctrl+Shift+T`). Drag
the lower slider to 65535, leave the upper at 65535, so only ceiling pixels are
red. Do NOT press Apply. If any red is inside your crop, those lanes are saturated.

**Exact line endpoints.** With a line selected, run the `getLine` macro from
section 4. If it prints `-1, -1, -1, -1`, a line is not the active selection
(you probably measured a rectangle after the line); redraw and run again.

**Centimetre-versus-pixel check.** Hover any pixel and read the status bar. A value
in parentheses like `x=3.80 (476)` means the image is calibrated (centimetres, then
pixel); a bare integer means pixels. This decides `coordinate_unit`.

---

## Deviations log

Record issues, deviations, and friction here as you run the pipeline, one line
each: date, section number, what happened, what you did. This is the raw material
for the next refinement pass.

- (add entries here)
