# Gel densitometry pipeline - design

Why the pipeline is shaped the way it is: architecture, decisions with the
alternatives that were rejected, and the data contracts between stages. The
companion `PROTOCOL.md` records what to type; this records why.

State as of 2026-08-11. This rewrite folds the prior verification logs (former
sections 9-14) into the body as settled facts. The verbatim primary-source records
they contain - the `.inf` transcription and its arithmetic, the vendor and ImageJ
reconciliation, and the live s001/s0002 measurement sessions - are provenance and
are RETAINED from the prior DESIGN revision (see the Appendix directive); they are
not reproduced here, because paraphrasing measured records loses fidelity.

---

## 1. Purpose and scope

### Goal

Replace a manual ImageJ workflow (crop, rotate, set lanes, plot lanes, draw the
baseline by hand, wand-select peak areas) with a scripted, auditable pipeline. It
now runs from a validated image all the way to a vector figure: per-band
integrated intensities, per-gel normalization against a reference lane, and
replicate aggregation into descriptive summaries and difference matrices, with
annotated overlays and provenance at every step.

### Scope correction (this supersedes the prior "stops at the CSV" statement)

The prior design scoped normalization, percent/fold computation, and plots OUT,
and declared that the pipeline stops at the measurement CSV. That is no longer
true and was the largest single change since the first revision. The pipeline now
owns two R stages as well:

- `gel_loading_normalization.R` - per gel, normalizes lanes to the reference lane
  (reference = 100 percent) and emits per-lane and per-band CSVs and a per-gel plot.
- `aggregate_loading_replicates.R` - across replicates, emits descriptive summaries
  (mean, sd, se) and all-against-all percent- and fold-difference matrices, and a
  vector aggregate plot.

What remains out of scope: molecular-weight calibration, smiling correction,
overlapping-band deconvolution, and cross-gel absolute normalization. Inferential
statistics are also out of the pipeline by decision (section 5.16); a standalone
reference for the house-style machinery is kept separately and reads no pipeline
file.

### The measurement actually needed

Within-gel comparison. One reference lane normalized to 100 percent, other lanes
expressed relative to it. Then, across replicate gels, the descriptive spread and
the pairwise differences. No cross-gel absolute normalization, no molecular-weight
calibration. This simplification is deliberate and should be preserved.

---

## 2. Architecture as built

One measurement program, not a chain of slices. The former Slice A (crop, lane
geometry) and Slice B (band detection, integration) were merged into a single
`measure_gel.py`: the image is read once and cropped once, and the band stage
consumes the lane geometry from memory, so the two halves cannot disagree about the
crop. The old on-disk `stage2_analysis_report.json` handoff between them is gone.

The full path, and who owns each step:

```
Fiji (manual)  landmarks + crop rectangle -> <stem>_preprocess.txt
stage 0  prepare_gel_image.py       creates <stem>_gel_analysis; locates .inf if present
stage 1  validate_gel_image.py      writes input_file_validation_report.json,
                                     deriving preprocess_sidecar.geometry_pixels + tilt
stage 2  measure_gel.py             reads that report; crops; detects; integrates;
                                     writes band_measurements.csv + band_detection_report.json
R (per gel)  gel_loading_normalization.R    -> loading_normalization_per_lane.csv (+ per-band, plot)
R (across)   aggregate_loading_replicates.R -> summaries, matrices, aggregate_bar.pdf
Illustrator  open the per-gel plot and aggregate_bar.pdf (vector, editable text)
```

Three architectural facts that are easy to get wrong:

- **The manual Fiji step is essential and cannot be derived from pixels.** The two
  well-line landmarks and the crop rectangle are human input carried in the
  `_preprocess.txt` sidecar. The pipeline is honest that this step precedes it.
- **Tilt is measured but never applied.** Stage 1 derives the tilt angle from the
  two landmarks and warns above 0.25 degrees; nothing ever rotates. Removing the
  in-line resampler means the crop the locator sees is the crop the integrator
  sees, and a two-point landmark on an already-resampled image would read a tilt of
  zero, so a derived tilt of 0.0 on a re-imaged gel is expected, not reassuring.
- **`well_index` is the join key, not `loaded_lane_index`.** `well_index` (0..14)
  is the physical lane position and is the key every downstream consumer joins on.
  `loaded_lane_index` also appears and is not the join key.

---

## 3. Confirmed facts about the input files (folded to conclusions)

Established from real files and real tag dumps, then reconciled against vendor and
ImageJ sources over successive verification passes. Conclusions only; the verbatim
records are retained per the Appendix.

- **Two instruments, one code path.** Typhoon (phosphor/radioactivity) writes
  `.gel/.img/.inf/.tif`; the Amersham Imager 680 (fluorescence/loading) writes
  `.tif/.png/.jpg`. Only the 16-bit `.tif` is measured. No flags distinguish them:
  everything that differs is read from the file's own tags. Pixel size is 200
  micrometres (Typhoon) or 79.9 (Imager 680); nothing assumes either.
- **Read the `.tif`; the `.inf` is provenance, not a validation target.** Earlier
  drafts treated the `.inf` as authoritative and even mis-transcribed it by one
  positional field; the resolved position is that the `.inf` is optional (recorded
  if present, tolerated if absent) and the `.tif` tags are the source of truth.
- **Linearity is UNRESOLVED for the Typhoon `.tif` and is the load-bearing caveat.**
  The prior "data is linear" statement did not survive reconciliation. Vendor docs
  call the `.tif` linear 16-bit grayscale and no source claims it is transformed,
  but that is not verified against these files and no cheap test settles it. Reports
  carry `encoding_verified: false` and name the test that would change it. If the
  `.tif` is not linear, every ratio is compressed toward 1 and nothing downstream
  notices. For the cooled-CCD fluorescence files the case is stronger, but those
  record an exposure time and are only comparable after dividing by it.
- **Orientation is bottom-left origin; polarity is MINISWHITE.** The image is
  mirrored top-to-bottom on screen; the pipeline flips it, so the operator measures
  the image as it opens. Every emitted image uses a matching reversed colormap.
- **The internal filename differs from the disk filename.** The `.inf`/TIFF may
  carry a different internal name than the file on disk. This is the root of the
  filename-mangling hazard: the sidecar's `measured_against_input_filename` must
  match the on-disk `.tif` name exactly or stage 1 hard-stops.
- **The numeric ceiling is below the container maximum.** Saturation detection
  cannot assume 65535 is the only at-ceiling value; the plateau statistic is the
  primary detector because the sensor can saturate below the numeric ceiling.

---

## 4. Coding conventions

The user's stated working defaults, not negotiable within this project.

- **Flat procedural.** Data in, data out. No classes, no OOP, no helper functions
  or wrappers built for a single call site. Inline logic where it is used. The only
  permitted helpers are the tagged-message emitter and the fail-fast error exit.
- **No premature abstraction.** Add it when the second real call site exists.
- **No indirection or nesting that does not pay for itself.** One longer readable
  block beats three short ones that force jumping around.
- **Full descriptive names carrying domain and unit information.** No abbreviations
  except tokens that behave as nouns (usb, tiff). No single-letter names, loop
  variables included. `lane_pitch_millimetres` not `pitch`.
- **Comments state why**, not what: the constraint, the failure prevented, the
  alternative rejected, the non-obvious platform detail.
- **ASCII only.**
- **Verify by executing, not asserting.** Show the output.
- **Messages to stderr**, tagged with source, so stdout stays clean for piping.
  Messages also accumulate and are serialized into the provenance JSON.

R stages follow the same spirit (mostly flat, minimal helpers, descriptive names,
ASCII). The operator's own downstream analysis scripts are house-style references
and are left untouched.

---

## 5. Design decisions, with rejected alternatives

Decisions 5.1 through 5.13 are folded to their settled form (corrections from the
later verification passes are applied in place). Decisions 5.14 onward are new since
the first revision.

**5.1 Grid model, not per-lane band search.** Lanes and band rows are a grid the
operator seeds with an expected lane count; bands are found within the grid. A free
per-lane search was rejected as fragile on faint and empty lanes.

**5.2 Expected lane count supplied by the operator.** Cheaper and more robust than
inferring it, and the operator always knows it. Empty lanes are counted.

**5.3 Both integration modes, always, clearly labeled.** Region and peak modes are
both emitted per band and labeled, rather than forcing one choice at measure time.

**5.4 Baselines: dual, and LOCAL.** Both baseline methods are reported. The live
s001/s0002 sessions corrected an earlier global-baseline assumption: because a real
crop includes bare plate, the baseline must be local per lane. The two baseline
rows per band carry a duplicated `reported_value` (see the de-dup contract in 6).

**5.5 Instrument handled by tags, not a flag.** An earlier "modality is a mandatory
argument" decision was superseded: nothing that differs between instruments needs a
flag, because it is all in the tags.

**5.6 Preprocessing writes a sidecar, not a new image.** The Fiji step records
landmarks and a crop rectangle into `_preprocess.txt`; it never writes a modified
image. No data is destroyed and the crop is reproducible.

**5.7 / 5.8 Rotation and crop come from the operator as landmarks.** The operator
supplies two landmarks and a crop rectangle, not an angle. Tilt is derived and
warned, never applied (section 2). Crop is a rectangle in the measured frame.

**5.9 Output directory: one per scan.** `<stem>_gel_analysis`, shared by the `.img`
and `.tif` of one scan (same stem), created by stage 0.

**5.10 CSV shape: long, not wide.** One row per band per baseline method. Joins and
aggregation are done downstream on `well_index`.

**5.11 Read the `.tif`.** Settled against real files; the `.img` is big-endian and
not a TIFF, and the `.tif` is the comparable, tag-bearing container.

**5.12 Saturation diagnostics.** Whole-image ceiling stats at stage 1; per-band
saturation and the plateau statistic at stage 2, with the ceiling inferred rather
than assumed to be 65535.

**5.13 Sensitivity analysis as a substitute for operator variance.** Realized as
the optional sensitivity harness (5.17): rather than ask the operator to repeat the
manual crop many times, perturb the crop and baseline knobs programmatically and
report the spread.

**5.14 Stable `gel_id` and additive schema (new).** Every measurement and R output
carries `gel_id` (the input TIFF stem) as a first-class field: top-level in
`band_detection_report.json` and as the first CSV column throughout. The report
schema was bumped 2 -> 3 additively; the R consumer reads `gel_id` from the report
when schema >= 3 and falls back to the directory name with a warning otherwise.
Rejected: keying downstream joins on directory paths, which are not stable across
machines or moves. The repeat relationship was removed from the sample sheet: a
gel does not know it is a replicate, so that fact lives only in the aggregation
manifest.

**5.15 Per-gel normalization is a join-normalize-plot layer, and re-measures
nothing (new).** `gel_loading_normalization.R` consumes `band_measurements.csv` and
the sample sheet, pins the single reference lane to 100, and expresses the rest as
percent-of-reference. It does not touch pixels. Exactly one `reference` role is
required; `empty` roles are excluded from normalization.

**5.16 Replicate aggregation is DESCRIPTIVE; inference is separated (new).** The
aggregator emits mean/sd/se and two all-against-all matrices, computed within each
replicate then averaged. Formulas, i relative to j:
`fold(i,j) = value_i / value_j`, `percent_difference(i,j) = (value_i - value_j)/value_j*100`.
It writes every intermediate state (read-together, processed, right-before-plot).
Inferential testing (paired t blocked on replicate, Holm within one named family,
exact p-values, effect sizes, an nlme consistency diagnostic) was deliberately kept
OUT of the pipeline and lifted into a standalone, simulated-data reference so the
machinery survives without adding statistics to routine output. The condition key
is a per-run constant; a blank key degrades to a per-well fallback with a warning;
one key set applies per manifest run. Rejected: baking a fixed factor model into
the pipeline, and silently averaging within-gel duplicates.

**5.17 Sensitivity harness perturbs crop and baseline, not tilt (new).** Because the
pipeline never rotates, a rotation knob would test a transform the pipeline does not
perform; the harness perturbs the crop rectangle (sidecar copies, full chain
re-run) and the baseline knobs (stage-2 CLI, stage-2 re-run) one factor at a time,
and reports per-lane normalized spread as a range. It recomputes the simple
normalization itself rather than calling R (a documented, deliberate duplication).
A perturbation that changes the detected lane count reindexes `well_index` and is
excluded from the spread and reported as its own coarser finding.

---

## 6. Data contracts between stages

The interfaces a consumer must not break.

- **Sidecar `_preprocess.txt`** (human, Fiji): `measured_against_input_filename`
  (must equal the on-disk name), image dimensions, `gel_migration_axis`,
  `coordinate_unit`, two landmarks, `crop_x/crop_y/crop_width/crop_height`,
  `expected_lane_count`.
- **Stage 1 report** `input_file_validation_report.json`: carries
  `preprocess_sidecar.geometry_pixels` (`crop_x_pixels` etc., `gel_migration_axis`)
  that stage 2 reads, the derived tilt, echoed geometry, and `encoding_verified`.
- **Stage 2 `band_measurements.csv`**: long form, join key `well_index` (0..14).
  `reported_value` is the per-band reported quantity and is DUPLICATED across the
  two baseline rows; a per-lane total de-dups on distinct
  `(well_index, canonical_band_index)` before summing. `gel_id` is the first column.
  `band_detection_report.json` carries `gel_id` top-level and `encoding_verified`
  forward; schema version 3.
- **Sample sheet `<stem>_sample_sheet.csv`** (human): `well_number` (1..15, each
  once; `well_index = well_number - 1` under the default `direct` flip),
  `sample_label`, `role` in {reference, sample, empty} with exactly one reference,
  and the factor columns. No `replicate` column (removed by decision 5.14).
- **Per-lane `loading_normalization_per_lane.csv`** (R): `gel_id` first,
  `whole_lane_total`, `percent_of_reference_whole_lane`, control-band columns,
  `role`, factor columns, provenance. This is the aggregator's input.
- **Manifest** (human): `experiment_id`, `analysis_path` (dir or any file inside;
  relative to the manifest), `gel_id` (verified against the per-lane CSV; hard stop
  on mismatch), `replicate` (unique per experiment), optional `notes`.
- **Aggregate outputs**: `combined_raw_data.csv` (stacked, empties kept,
  pre-normalized preserved), `processed_data.csv` (empties dropped, reported value
  selected, `condition_label`), and per experiment `summary_data.csv`, the matrix
  files, and `aggregate_bar.pdf` (vector cairo, editable text).

---

## 7. Verification history (folded)

The prior logs (former sections 9-14) recorded, and this design now takes as
settled: the section-2 arithmetic held except for the one-off `.inf` transcription,
corrected; the `.inf` was demoted from validation target to provenance; the `.img`
is big-endian and not a TIFF; the output location and which-file-to-read questions
were resolved to the `.tif` and one directory per scan; the ceiling-below-container
inference was implemented and tested; orientation is recoverable from
`ImageDescription`; and the live horizontal-gel session (s001) forced the
local-baseline correction in 5.4. The three-way encoding scatter did NOT decide
linearity, which is why `encoding_verified` remains false.

Added since: 3.1 (stable `gel_id`, schema 2 -> 3) verified end to end on s0002; 3.2
(aggregator + standalone stats reference) verified on synthetic multi-gel fixtures
across keyed, keyless, mismatch, and incomplete-pairing paths, with a hand
recomputation of the within-replicate-then-average matrices matching the files
exactly, and the plot confirmed vector; 3.3 (sensitivity harness) verified on the
real s0002 gel, including the lane-grid-change exclusion.

---

## 8. Open questions and deferred work

- **Encoding verification.** `encoding_verified` is false for the Typhoon `.tif`
  until the named linearity test is run; until then every ratio is provisional.
- **EMSA / radioactivity readout is a different, not-yet-built R path.** The loading
  path measures percent-of-reference; an EMSA readout is fraction bound
  = shifted/(shifted+free). The aggregator consumes the loading per-lane CSV; a
  fraction-bound producer does not exist yet, and real EMSA numbers also wait on
  pristine s0003 sidecars (blocked on operator Fiji measurement; current s0003
  geometry is untrusted).
- **R-in-harness deferred.** The harness recomputes normalization itself; driving
  `gel_loading_normalization.R` per cell is a later non-breaking extension. If the R
  normalization formula changes, the harness duplication must be revisited.
- **Multi-key manifests.** One condition key set per manifest run; differing keys
  per experiment require one manifest per experiment.

---

## Appendix. Retained primary-source records (directive, not reproduced)

The following from the prior DESIGN revision are measured, primary-source records
and must be RETAINED verbatim rather than paraphrased: the verbatim `.inf`
transcription and its arithmetic confirmation (former section 2), the vendor and
ImageJ reconciliation (former section 10), the "settled by running against real
files" results (former section 11), the `.inf` off-by-one correction (former
section 12), the stage-1 verification log (former section 13), and the live s001
session (former section 14). Their CONCLUSIONS are folded into sections 3, 5, and 7
above; the records themselves stay as provenance. When merging this rewrite, append
those sections unchanged.
