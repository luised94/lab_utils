# Gel lane extraction: picker and CLI

The last friction step before plotting. Given one gel's analysis, produce one
traceable number per lane (per sample), either for a chosen consensus band or for
a migration span you draw. That per-lane CSV is what gets plotted for this gel, or
stacked across repeats by the (separate, later) aggregator into means and plots.

Position in the pipeline: this runs after analyze_gel.py has written
band_measurements.csv, and it also reads the ImageJ export manual_lane_profiles.csv
for the raw profile and the ROIs. It does no band detection of its own; it reads
the profile, the consensus bands, and the metrics, and lets you choose.

## The two ways to run

Command line, when you already know the band or span:

    uv run extract_lane_values.py '<gel_analysis_dir>' --band 3 --quantity apex
    uv run extract_lane_values.py '<gel_analysis_dir>' --region 34 45

Picker, when you want to see the gel and choose:

    uv run serve_gel_picker.py '<gel_analysis_dir>'
    # open http://localhost:8080, choose a band or drag a span, Write CSV, Close

The picker never computes the number itself. On "Write CSV" it runs
extract_lane_values.py as a subprocess, so the value the browser writes is
byte-for-byte the command-line value; the two cannot drift.

## What a selection means

- Band: consensus band N taken straight from band_measurements.csv. quantity=area
  is the region-net area at the shared, fixed migration window (smile-affected);
  quantity=apex is the per-lane apex height above baseline (smile-robust).
- Region (span): the raw sum of the profile (raw_value) over the migration window,
  per lane. raw_value is already the width-summed signal with the plate-background
  median subtracted and clipped at zero, so this uses the pipeline's background
  correction. A fixed migration window is smile-affected; drawing a span is you
  choosing to treat that stretch as one quantity, independent of how the detector
  split it. Preferred when a lane (a negative control, a defective prep) may not
  contain the detected bands at all.

## Files

    extract_lane_values.py   the extractor (stdlib only); writes the per-lane CSV
                             and a checks JSON into the gel directory
    serve_gel_picker.py      the local picker server (numpy, tifffile); no web
                             framework, PNG encoded with stdlib zlib
    gel_picker.html          the picker page; fetches /payload, POSTs the selection
    make_payload.py          builds an inlined static copy of the page (no server),
                             for a frozen, shareable picker

## Deferred, on purpose (with the trigger to revisit)

- Fuller per-lane peaks overlay. The picker flags each lane's per-band apex. The
  complete set of per-lane detected peaks lives in band_candidates_per_lane.csv;
  add it to the payload when that file is supplied and a richer peak overlay is
  wanted.
- Valley-to-valley net baseline for a span. The span default is the raw sum
  (plate-background already removed). The band measurement additionally removes a
  local valley-to-valley baseline; add that as a span option only if a real gel
  shows the raw span sum is insufficient, rather than building it on spec.
- Leaner path without analyze_gel. Pure integrate-span needs only the ImageJ
  export plus lane identity, not band detection. If band mode and the guiding
  metrics stop being needed, the picker could read the profile and the sample
  sheet directly and drop the band_measurements dependency.
- Aggregation across replicate gels (means, SD, plots) is a separate step, later.

## Conventions

ASCII only. Flat and procedural; emit_message / die (and record_check in the
extractor) are the only helpers, duplicated into each script rather than shared,
so each stays runnable on its own. uv inline dependency block at the top of each
script. Every stage writes a checkable artifact (the per-lane CSV and a checks
JSON) that is inspected every run.
