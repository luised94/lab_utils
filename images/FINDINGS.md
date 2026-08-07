# Findings - gel geometry and provenance probes

Self-contained lab notebook for the geometry and file-provenance work. Records
what the second test gel (the radioactivity / phosphor s0003) taught us, which
probes were reliable, and which probes failed and must not be re-tried blind. The
failed list is the load-bearing part: without it, the next hard gel invites the
same dead ends.

Companion to DESIGN.md (why) and PROTOCOL.md (what to type). This file is where
observations and rejected approaches live so they survive past one work thread.

---

## 1. Observations from s0003 (phosphor)

- **The gel is not tilted; it smiles.** The loading line is a curve, not a
  straight tilted line. A straight-line fit to the loading ridge leaves 4.3 px
  RMS residual; a quadratic leaves 1.6 px. The bow is ~14 px (2.7 mm at
  200 um/px, ~31 percent of a lane pitch) off the chord between the two ends,
  worst mid-span and at the edge lanes. Consistent with drying distortion
  (edges pulled and twisted as the gel dries). This is the first gel in the set
  where the geometry departs from straight-and-parallel; s0002's wells sit
  0.22 degrees off vertical, effectively straight, so nothing here surfaced
  before.

- **A DNA-binding lane shows up as a ridge spike.** One sample binds DNA and its
  band shifts wholesale, so that lane's ridge point jumps off the smile curve.
  Any curve fit to the loading ridge must be robust (reject outliers past N
  residual-sigma, refit) or one shifted band drags the whole curve. The rejected
  points are not noise to discard: a lane sitting off the smile by more than the
  fit residual is a candidate shift/bind, which in an EMSA context is the signal
  of interest. Outlier rejection and shift detection are the same code.

- **The 2-landmark model is structurally blind to smile.** Two landmark points
  define a chord, and a chord cannot represent a curve. The tilt the sidecar
  derives is the chord's angle; on a smiling gel the chord can read near-zero
  tilt while the middle bows out 14 px. This is not a threshold that is tuned
  wrong. There is no smile term in the model at all.

- **The delivered s0003 .tif was an ImageJ-resaved derivative.** The file named
  `..._rotated_...` carried `ImageDescription = ImageJ=1.54g`, had its Amersham
  Typhoon instrument tags stripped, and was resampled: ImageJ had rotated it
  ~6 degrees (measured from the data-rectangle edge), which tilted the frame and
  filled the exposed corners with zeros. The pristine original
  (`..._orc1,3,4_...`, no `rotated`) carries the full Typhoon tags, 200 um/px,
  and is un-resampled and straight. Rotating in ImageJ to "straighten" cannot fix
  a smile anyway; it only resamples the pixels and tilts the frame. The pristine
  scan is strictly the better input. The sidecar pointed at the derivative, which
  is how the wrong file nearly got measured.

- **The zero floor is a scan-pad border, not clipping.** ~21 percent of pixels
  are exactly zero in the pristine original too, so the floor is intrinsic to the
  Typhoon export (unscanned canvas, whole top rows at 0), not an ImageJ artifact.
  It sits entirely outside the crop; the gel interior has zero exact-zeros. The
  radioactivity background is genuinely low, not clipped. The existing
  whole-image floor_population warning fires on this border and is benign here,
  which motivated the floor-location split (report inside-crop vs outside-crop and
  escalate only when it is inside).

---

## 2. Probes that worked (reliable, reusable)

- **TIFF tag inspection.** Instrument tags (Model, Make, Software) plus the
  ImageDescription decide provenance cheaply and decisively: Typhoon vs Imager 680
  vs an ImageJ resave. This is the basis of the provenance guard added to stage 1.

- **Zero-population location.** Splitting the exact-zero count into inside-crop vs
  outside-crop tells a benign scan-pad border apart from real in-gel clipping.

- **Robust ridge quadratic.** Trace the loading-line ridge row by row, fit a
  quadratic with outlier rejection, and use the line-vs-quadratic residual ratio
  as the tilt-versus-curve discriminator (a curve fits far better than a line)
  and the sag off the endpoint chord as the smile magnitude. This is the honest
  smile metric and the basis for the deferred smile work (Option A below).

- **Eyeball against a true-axis reference.** Overlaying true horizontal/vertical
  reference lines on the image is what made the smile visible when automatic
  angle estimators disagreed. This motivated the crop-preview reference gridlines.

---

## 3. Probes that FAILED - do not re-integrate

Each of these was tried on s0003 and gave a wrong or meaningless answer. They are
recorded so nobody spends another afternoon rediscovering the failure.

- **Straight-line fit to the loading column.** Reported ~1 to 6 degrees of
  "tilt" that did not exist. The feature is a curve; a line through a curve
  returns the average of the bow, which is a tilt artifact, not a tilt.

- **Cross-correlation of lane profiles between the two ends of the gel.** The far
  ends hold different content (wells at one end, migration front at the other),
  so the lane patterns do not correspond and the correlation peak was ~0.03.
  Meaningless shift estimate.

- **Deskew / Radon sweep (rotate, collapse, maximize profile peakiness).** Pinned
  to the search-range boundary. The metric is dominated by the strong well column
  and the bound-DNA blob, not by the lane periodicity, so it optimizes alignment
  of the wrong feature.

- **Fixed-window lane-centre tracker.** Tracks each lane's centroid in a fixed
  y-window across migration strips. On a genuinely tilted or strongly smiling
  lane the centroid leaves the window and the estimate clamps toward the window
  centre, so a real tilt reads as ~0. It cannot distinguish flat from tilted and
  must not be used as a tilt detector.

General lesson: global automatic tilt/deskew estimators are unreliable on real
gels with dominant blobs and non-straight geometry. Prefer tag inspection,
zero-location, and a robust ridge fit over any whole-image angle search.

---

## 4. Deferred: the smile metric (Option A)

When built, this belongs in stage 2 (measure), because it needs the crop and the
lane geometry. Design brief:

- Take the two existing landmarks only as a search corridor. Do not add a third
  landmark: two landmarks stay as the operator's job; the curve is discovered
  from the pixels, so hundreds of traced rows outvote one shifted band. A 3-point
  parabola would break on exactly the DNA-binding lane above.
- Trace the loading-line ridge between the landmarks, fit a robust quadratic,
  report sag in px, mm, and fraction of lane pitch (pitch-fraction is the
  scale-free unit that travels across the two instruments), warn above a
  threshold, and flag outlier lanes as candidate shifts.
- Degrade to the current straight chord, and say so, when the loading line is too
  faint to trace. Do not guess.
- Consequence to model: handling smile means the lane grid stops being straight
  parallel lines and becomes per-lane curves, touching the grid fit, the
  integration windows, and every overlay. Larger than a tilt warning; scope it
  with the geometry/figure rework, not as a quick change.
