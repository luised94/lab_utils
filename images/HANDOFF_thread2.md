# Handoff: stage-1 provenance/geometry + empty-lane prediction -> next thread

Self-contained. Pasteable into a fresh context. Records what this thread landed,
what is verified vs provisional, and what is deferred or blocked and why.

## What this thread delivered (7 commits over the received baseline)

All in `images/`. Baseline commit is `baseline (as received)`.

Stage 1 (`validate_gel_image.py`, final md5 089eb8c8):
- `input_file_provenance` (warn): flags an ImageJ resave stamp or an unrecognized
  acquisition signature against a named allowlist. Catches the wrong-file case that
  nearly measured an ImageJ derivative instead of the pristine scan.
- Sidecar schema v3: rotation-provenance fields (`rotation_applied_degrees`,
  `rotation_interpolation_method` default none, `rotation_enlarged_canvas`,
  `rotation_reference`), version-aware so v2 sidecars still validate. Cross-checks
  warn on a bilinear/bicubic resample and on sidecar-vs-file rotation mismatch.
- `floor_population_within_crop` (warn): splits the zero floor inside vs outside the
  crop, escalates only when zeros are in the measured region.
- Crop-preview true-axis reference gridlines (dashed teal), which is what made the
  s0003 smile visible.
- Template `preprocess_sidecar_template.txt` (md5 b1e2bac1) carries the v3 block.

`FINDINGS.md` (md5 03f9e1e8): the s0003 geometry/provenance probes, the ones that
worked (tag inspection, zero-location, robust ridge quadratic), and the four that
FAILED and must not be retried (straight-line fit to a curve, cross-correlation
across mismatched ends, deskew/Radon hijacked by blobs, fixed-window tracker that
clamps on real tilt). Section 4 is the deferred smile-metric brief.

Stage 2 (`measure_gel.py`, final md5 b59ca42f):
- Empty-lane / pitch prediction: loaded-but-defective comb positions (negative
  control, dead mutant) now get a real integrated value at the consensus positions
  instead of an absent row. The comb fitted to populated lanes places every missing
  position; integration measures it (consensus_only), so the value is a real
  near-zero, not an asserted 0.
- New band_measurements columns: `loaded_lane_index` (0..loaded-1, the flip-agnostic
  key R's plate map joins on), `lane_detection_status` (detected|predicted_from_comb),
  `prediction_span` (interpolated within the populated span, extrapolated beyond).
  `well_index` stays the comb slot.
- Predicted lanes use the median detected strip height (fair footprint), clamped to
  neighbour midpoints so they cannot bleed into an adjacent loaded lane.
- Guards (warn): fallback pitch; a comb slot at/beyond expected_lane_count; a
  predicted lane carrying >40% of the strongest band (missed lane or strong-neighbour
  bleed). Index-origin marker on the consensus overlay for flip-checking by eye.

## Verified vs provisional

- s0002 (fluorescence, REAL geometry): the trustworthy demonstration. Detected wells
  [0,11,12,13,14] byte-identical to pre-feature; wells 1-10 predicted, interpolated,
  real near-zero; 120 -> 360 rows; self-check quiet; band_centers stays detection-only.
- s0003 phosphor (radioactivity): machinery-generalization check ONLY. It runs
  cross-instrument (200 um Typhoon, low background) and emits the same contract, but
  its geometry is NOT trustworthy: the only sidecar was measured on the ImageJ-rotated
  derivative, and an auto-crop from pixels made detection worse (8 lanes -> 2). Its
  lane count, assignments, and any smile-split are provisional until a proper sidecar.

## Blocked on operator work (cannot be done from pixels)

- Phosphor needs a Fiji-measured sidecar on the PRISTINE original
  (`..._orc1,3,4_...`, not `..._rotated_...`), with `expected_lane_count=10` (comb is
  10; 9 loaded, all visible because radioactivity is sensitive). Crop the loaded lanes
  only; the comb predicts the rest. For an edge defective control, pad the crop ~one
  comb-pitch past the last visible lane on that side; interior defective lanes are
  automatic; truly unloaded edge wells can be ignored.

## Deferred by decision (do not block R's columns)

- Smile correction. The s0003 smile (loading line bows ~14 px, ~31% of a lane pitch,
  edge-worst) makes the consensus split one band into two canonical bands (top and
  bottom of the bow), which over/under-estimates those rows. THE PRIMARY FIX IS THE
  BENCH, not software: the distortion is a drying artifact (vacuum-heat dry, then
  paper-towel padding into the phosphor cassette puts uneven edge tension; edge-worst
  matches). Software options, ranked: (1) warn/flag smile-split canonical bands
  (member lanes partition by stacking, not overlap) - safe, inert on straight gels,
  unvalidated-active until a good smiling sidecar exists; (2) integrate detected cells
  at their own smile-following centre - drops consensus-window comparability, breaks
  byte-identical; (3) full ridge-model smile metric (FINDINGS section 4) - correct,
  larger. None is validatable without a correctly-cropped smiling gel, so none was
  built. Recommendation: fix drying; keep software to flag-only if built at all.
- Schema doc (one page: what a row means in each of the four CSVs, units, coordinate
  frame, the two gotchas that occupancy carries locally_detected/consensus_only and
  that 120 rows is 60 cells x 2 baseline methods, plus the four new lane columns).
- Metadata CSV (echo micrometres_per_pixel, migration origin, encoding_verified,
  overall_status + warnings, schema versions, gel id, so R reads CSV-only).
- R script linking to measure_gel.py output. House style captured from the two
  loading scripts (source()-driven, script-dir + env-var path resolution, tidyverse
  unqualified, fail-fast Cairo/Arial, PLOT_CONFIG + cairo_pdf, exact-p never stars).
  Different contract: read band_measurements.csv (not pre-normalized Excel), collapse
  to one row per (well, canonical_band) or filter one baseline_method, honor the trust
  flags, and do gel-shift normalization (fraction bound = shifted/(shifted+free), or
  normalize to a reference band). loaded_lane_index is the plate-map key; the placement
  flip lives in R's plate map by decision, not the sidecar.
- Also pre-existing/pipeline-internal, unchanged: background/mask de-dup (blocked on
  the ceiling choice), run-log-into-JSON, the rotated-vs-unrotated crop inconsistency
  note. encoding_verified stays false until a dilution series proves linearity.

## Open questions to carry

- Does drying improve enough to drop the smile entirely (then software smile handling
  is unnecessary), or is a flag-only detector still wanted for the archive gels?
- Phosphor constants: fragility 0.15, elevation, rolling-ball width are s0002-tuned;
  the one approximate phosphor run showed a much higher baseline-clip rate (confounded
  by geometry). Revisit once a proper phosphor sidecar exists.

## What to bring to the next thread

- measure_gel.py, validate_gel_image.py, prepare_gel_image.py, the template,
  DESIGN.md, CONVENTIONS.md, PROTOCOL.md, FINDINGS.md, this handoff.
- s0002 TIF + sidecar (real geometry). The phosphor TIFs (both), and a NEW
  Fiji-measured sidecar for the pristine original with expected_lane_count=10.
- The two reference R loading scripts (style reference; stay as-is).
