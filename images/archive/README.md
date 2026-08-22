# Archived scripts: why they are here, and what is worth extracting

These are the OLD automated TIFF->measure pipeline and its old harness, superseded
by: the ImageJ macro (export_lane_profiles.ijm) for profile export, analyze_gel.py
for the manual measure stage, and tests/run_pipeline_tests.py for regression. They
are kept, not deleted, per the house convention ("kept for reference and logic
reuse, not run") and because they contain diagnostic and QC logic the live pipeline
does not yet have. Nothing live imports or runs them.

This file is the map of what is worth reaching back for, and what is not, so a
future thread does not have to re-read ~7,600 lines to find out. Each item is
marked ALREADY LIVE (do not re-extract), GAP (genuinely absent; a real candidate),
or DROP (tied to the old automated path; not worth carrying). ASCII only.

===============================================================================
ALREADY LIVE - the load-bearing math was inlined into analyze_gel.py; do NOT
re-extract, you would be duplicating what is already there
===============================================================================

The core measurement contract from measure_gel.py is already present, copied
verbatim into analyze_gel.py per the no-sharing rule (analyze_gel's comments cite
archive/measure_gel.py at each site):

  - valley-to-valley baseline (flanking minima, sloped line)
  - rolling-ball / opening cross-check baseline and the fragility test (the two
    nets disagree by >= 0.25 of lane scale -> report the cluster-inclusive net)
  - consensus band clustering across lanes
  - the reported-value contract (net area vs apex height; saturation handling)
  - saturation via the at_ceiling_count column (analyze_gel consumes it directly)

If a future reorg wants these as a SHARED module rather than duplicated, that is
open design decision #1 (shared-vs-duplicate); resolve it there, do not pre-extract.

===============================================================================
GAP - genuinely absent from the live pipeline; the real extraction candidates,
in rough value order. Each is a diagnostic/QC layer, not core math.
===============================================================================

1. TIFF-quality input gate  (from archive/validate_gel_image.py)   HIGH VALUE
   The live pipeline starts at manual_lane_profiles.csv and trusts it. It never
   checks the TIFF the macro measured. validate_gel_image.py gates exactly the
   things that silently corrupt densitometry at the source, as named hard/soft
   checks with a checks JSON in the same idiom the live scripts use:
       - compression is not lossy (a JPEG-compressed gel is unmeasurable)
       - SamplesPerPixel == 1, BitsPerSample adequate, SampleFormat sane
       - effective bit depth and ceiling (real dynamic range vs container)
       - saturation-spike and floor-population fractions
       - resolution / mm-per-pixel present and consistent with ImageDescription
   Why it is a gap worth closing: these are pre-conditions the macro assumes but
   nothing verifies. Extraction shape: a stdlib+tifffile validate_gel_image step
   that runs BEFORE the macro (or beside it), writing input_file_validation_report
   .json, wired into the harness as a new early stage. Note it needs a TIFF
   fixture, which the harness has so far avoided (CSV-only); pin the TIFF by
   sha256 per the handoff's TIFF note if this is built.

2. Measurement diagnostic overlays  (from archive/measure_gel.py)   MEDIUM VALUE
   measure_gel.py emits a rich set of QC figures the live analyze_gel.py does not:
       overlay_reported_values.png        the reported number on each band, on the gel
       baseline_comparison.png            valley-to-valley vs rolling-ball, per lane
       saturation_derivation.png          how each saturation flag was derived
       figure_net_vs_apex.png             net-area vs apex-height per band (the
                                          smile-sensitivity picture)
       figure_metric_compare.png          metric agreement across lanes
       consensus_center_map.png           band centre alignment across lanes
   analyze_gel.py already emits band_detection_overlay.png and
   band_measurement_overlay.png, so the BASIC visual gate exists. The above are
   the DEEPER diagnostics for when a gel looks wrong and you need to see why.
   Extraction shape: port individual figures into analyze_gel.py as optional
   outputs (guarded by a --diagnostics flag) rather than lifting the whole plotting
   layer; each figure is largely self-contained. Do this only when a real gel
   forces the question ("why is this band's number off?"), not speculatively.

3. Baseline / saturation warning THRESHOLDS as tunables  (from measure_gel.py)  LOW
   measure_gel.py exposes --saturation-warn-fraction, tilt-warning threshold,
   width-adequacy climb fraction, canonical-spread warn fraction, etc. analyze_gel
   hard-codes its equivalents. If the live pipeline ever needs these tunable, the
   values and their meanings are here. Low value now: hard-coded is fine until a
   gel needs a different threshold.

===============================================================================
DROP - tied to the old automated path; not worth extracting
===============================================================================

  - prepare_gel_image.py (Stage 0 path normalisation / sidecar handling): the
    manual path does not use it; the macro handles its own I/O.
  - inventory_tiff_tags.py as a standalone step: its useful tag-reading is already
    folded into validate_gel_image.py's checks (item 1 above). Keep only if item 1
    is extracted and wants the raw-tag dump for provenance.
  - sensitivity_harness.py: the old harness. Fully superseded by
    tests/run_pipeline_tests.py. Reference only for how Stage 0->1->2 was wired.
  - validate_manifest_sheet.py: strict subset of the live validate_manifest.py.
    Nothing to extract; the live one is better in every respect.

===============================================================================
Recovery
===============================================================================

git mv preserved history:  git log --follow archive/<file>
Original SHA these were archived from is in the commit that created this file.
