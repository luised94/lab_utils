# Gel densitometry pipeline - design and handoff

Status: design agreed, no code written yet.
Target location: `~/personal_repos/lab_utils/images/`
Date of this document: 2026-08-03

This document is the output of a design conversation. It is intended to be read
before any implementation begins, and to live in the repository as the record of
why the pipeline is shaped the way it is.

---

## 1. Purpose and scope

### Goal

Replace a manual ImageJ workflow (crop, rotate, set lanes, plot lanes, draw
baseline by hand, wand-select peak areas) with a scripted pipeline that emits a
CSV of per-band integrated intensities, plus annotated overlay images that make
the numbers auditable by eye.

### Downstream use, explicitly out of scope for this pipeline

The user has a separate script (likely R) that consumes the CSV and computes
percent-of-wildtype, fold differences, and plots. None of that belongs here.
This pipeline stops at the CSV.

### The measurement actually needed

Within-gel comparison only. Wildtype lane normalized to 100 percent, other lanes
expressed relative to it, with a negative control lane establishing background.
No cross-gel normalization, no molecular weight calibration, no absolute
quantification. This is a substantial simplification and should be preserved.

### In scope

- Single TIFF/IMG in, validated against its own header
- Grid-based lane and band-row detection
- Per-cell integrated intensity, region mode and peak mode, dual baselines
- Saturation and plateau diagnostics per cell
- Annotated overlay images and per-lane profile plots
- Tidy CSV
- Provenance JSON capturing config, instrument header, git state, checksums

### Out of scope, deliberately

- Normalization, percent/fold computation, statistics, bar and line plots
- Molecular weight calibration from a ladder (ladder often absent; sizes not needed)
- Smiling correction
- Overlapping-band deconvolution
- Cross-gel comparison
- Optical density transform for transmitted-light stains (structure allows it,
  we do not build it now)

### On the boundary, flag before folding in

- Multi-page TIFF handling beyond detecting and reporting page count
- Automatic gel-edge detection for cropping
- Any per-lane vertical offset tracking to handle smiling

---

## 2. Confirmed facts about the input files

All of the following was established from a real `ls -la` and a real `.inf` file,
not assumed.

### The Typhoon writes four files per scan

Example scan (repeat 1):

```
20220303-wt,4r,5,6-sofa-repeat-1-1000-[Phosphor].gel   1,971,357 bytes
20220303-wt,4r,5,6-sofa-repeat-1-1000-[Phosphor].img   1,968,750 bytes
20220303-wt,4r,5,6-sofa-repeat-1-1000-[Phosphor].inf       1,170 bytes
20220303-wt,4r,5,6-sofa-repeat-1-1000-[Phosphor].tif   1,970,177 bytes
```

### The `.inf` is a plain-text header, positional then key-value

Verbatim contents of repeat 1:

```
FLA_IMAGE_FILE
20220303-114817-[Phosphor].img
200
200
16
1125
875
1
5
Thu Mar 03 11:51:36 2022
1646326296
Amersham Typhoon
*** more info ***
42
S,RangeLow=2
S,RangeHigh=258
H,ScaleType=Linear
S,MethodName=-
S,LaserAndFilter=[Phosphor]
S,LaserName=635 nm
S,FilterName=IP BP390
S,PMTType=Bi-alkali
S,V=300
H,LaserPowerMode=Low
H,Filter=0/2
S,MainStartPosition=26[mm]
S,SubStartPosition=20[mm]
S,ScanMode=Phosphor
S,ScanNumber=1/1
S,StageAndAreaName=20x40
S,ScanSpeed=Normal
H,Correction1=1.0802
H,Correction2=Off
H,Correction3=Off
H,Correction3A=On
H,Correction4=On
H,Correction5=On
H,Correction6=On
H,Correction7=26531/26553
H,Correction8=Off
H,Correction8_DateTime=-
H,Correction9=Off
H,Correction9_DateTime=-
S,Shading=C:\ProgramData\GE Healthcare\Amersham Typhoon\Data\ShadingData\Phosphor_3_0_1_1.shading
H,Shading_DateTime=2020/02/17 21:57:42
H,CorrectionA=Off
H,Correction10=Off
S,Invert=Off
H,SubMotorMode=1
H,SignalProcess1=On
H,SignalProcess2=Linear
S,Orientation=bottom-left
S,SerialNumber=06230698
S,FirmwareVersion=303
S,FPGAVersion=10
S,Software=Amersham Typhoon Scanner Control Software 2.0.0.6
```

### Positional field interpretation

| Line | Value | Interpretation | Confidence |
|------|-------|----------------|------------|
| 1 | `FLA_IMAGE_FILE` | Magic string, format identifier | High |
| 2 | `20220303-114817-[Phosphor].img` | Internal acquisition filename | High |
| 3 | `200` | Pixel size, main axis, micrometres | High |
| 4 | `200` | Pixel size, sub axis, micrometres | High |
| 5 | `16` | Bits per pixel | High |
| 6 | `1125` | Width in pixels | High, confirmed by arithmetic |
| 7 | `875` | Height in pixels | High, confirmed by arithmetic |
| 8 | `1` | Unknown | None - do not interpret |
| 9 | `5` | Unknown | None - do not interpret |
| 10 | `Thu Mar 03 11:51:36 2022` | Acquisition time, human readable | High |
| 11 | `1646326296` | Acquisition time, Unix epoch | High, cross-checks line 10 |
| 12 | `Amersham Typhoon` | Instrument | High |
| 13 | `*** more info ***` | Section separator | High |
| 14 | `42` | Count of key-value lines following | High, verified by counting |

### Arithmetic confirmation of the raw format

```
1125 * 875 * 2 bytes = 1,968,750
.img filesize        = 1,968,750
```

Exact match. The `.img` is headerless, uncompressed, 16-bit, no offset, no
row padding.

Repeat 2 is 1,875,000 bytes = 937,500 pixels, a different scan area. Dimensions
must always come from the file's own `.inf`. Never hardcode.

### Container overhead

```
.img  1,968,750  (raw baseline)
.tif  1,970,177  (+1,427 bytes)
.gel  1,971,357  (+2,607 bytes, ~1,180 more than .tif)
```

Consistent with both being uncompressed 16-bit TIFF wrappers around the same
pixel data, with `.gel` carrying roughly 1.2 KB of additional private metadata,
presumably the MD_ tags (65000-65003).

### Data is linear

`H,ScaleType=Linear` and `H,SignalProcess2=Linear`.

This defuses the main risk with phosphorimager data. Fuji/GE/Amersham
instruments can write square-root-encoded data, where raw pixel values must be
squared and scaled to recover PSL (photostimulated luminescence). Integrating
square-root-encoded values as if linear compresses ratios toward 1 - a real 3x
reads smaller - and normalizing to wildtype does not fix it.

These files are linear. No inverse transform required.

This must still be verified empirically. A header field is a claim; the pixels
are the evidence. See stage 1 cross-check below.

### Orientation is bottom-left

`S,Orientation=bottom-left`.

The first row of bytes in the `.img` file is the **bottom** row of the gel.
Naively reshaping with numpy produces a vertically flipped image.

Why this matters, in descending severity:

1. Band row indices invert. Band 1 becomes band N. A vertical flip does not
   affect lane ordering, so the downstream percent-of-wildtype comparison would
   still pair lanes correctly and produce entirely plausible numbers with every
   band label wrong. This is the exact "plausible but wrong" failure mode the
   pipeline is meant to prevent.
2. The overlay will not match what the user sees in ImageJ, making visual
   verification confusing instead of reassuring.
3. The `.gel`/`.tif` may already have been flipped by the writing software,
   since TIFF carries its own orientation tag. `.img` and `.gel` may therefore
   disagree in row order while containing identical data.

Handling: read the field explicitly, flip to top-left convention immediately
after load, record in provenance that a flip was applied and why.

### Polarity

`S,Invert=Off`. Phosphor emission: bright bands on dark background. No optical
density transform. Must agree with histogram-mode polarity inference.

### Pixel size 200 micrometres per pixel

Image is 225 mm x 175 mm. A 5 mm lane is 25 pixels; a 2 mm band is 10 pixels.

**Consequence for config design:** every length parameter (smoothing window,
minimum peak separation, rolling-ball radius, region window height) is specified
in **millimetres** in the config and converted to pixels using the pixel size
from the `.inf`. Otherwise a scan at different resolution silently changes the
analysis.

### Size is a non-issue

About 1 megapixel. 2 MB as 16-bit, 8 MB as float64. Everything runs in well
under a second. No performance considerations anywhere in this design.

### Internal filename differs from disk filename

`.inf` line 2 says `20220303-114817-[Phosphor].img`; the disk name is
`20220303-wt,4r,5,6-sofa-repeat-1-1000-[Phosphor].img`. Files were renamed after
scanning.

Consequence: match sidecar files by **disk stem**, never by internal name.
Record both in provenance.

### Timestamps cross-check

`Thu Mar 03 11:51:36 2022` and epoch `1646326296` agree (11:51 EST = 16:51 UTC).
Two independent representations of the same instant, so the script can assert
they match and catch a corrupted header.

Disk mtime is 2022-03-10, a week later. **mtime is not acquisition time.** The
`.inf` is authoritative.

### Uninterpreted fields, to be recorded verbatim and flagged

- `S,RangeLow=2` / `S,RangeHigh=258` - possibly display windowing, possibly
  gradation/latitude. Unknown. Worth a web query.
- Positional lines 8 (`1`) and 9 (`5`) - unknown.

### PMT voltage is comparability metadata

`S,V=300`. If repeats were scanned at different PMT voltages, absolute
intensities are not comparable between gels. Within-gel percent-of-wildtype is
unaffected, so this does not harm the analysis, but the script should surface it
so a difference between repeats is noticed rather than missed.

---

## 3. Environment

- WSL Ubuntu on Windows, CPU only
- Python managed with `uv`; scripts self-contained via PEP 723 inline metadata
  (`# /// script`). Lock/requirements file written once dependencies are settled.
- Repo at `~/personal_repos/lab_utils/`, script goes in `images/`
- `.gitignore` reportedly already covers Python and uv artifacts - **verify**

### Data location and the working copy

Source data lives under `/mnt/c/Users/liusm/MIT Dropbox/...`. This is a Dropbox
folder synced from a work device; the home device can see it but is not where the
work is meant to happen.

**Decision:** clone test files to a local working location for development. Point
at the real paths only when running for real. This avoids syncing every
intermediate, avoids reading a partially-synced file mid-write, and avoids DrvFs
metadata slowness.

Note that `/mnt/c` is DrvFs: permissions show as 777 regardless of what is set.
Harmless at this scale.

### Filenames are hostile

Real filenames contain spaces, commas, and square brackets. Brackets are glob
metacharacters in bash. Path handling must be robust to all of this and any
documented shell examples must be quoted.

---

## 4. Coding conventions

These are the user's stated working defaults and are not negotiable within this
project.

- **Flat procedural.** Data in, data out. No classes, no OOP, no helper
  functions or wrappers built for a single call site. Inline logic where it is
  used. The only permitted helpers are the tagged-message emitter and the
  fail-fast error exit.
- **No premature abstraction.** Add it when the second real call site exists.
- **No indirection or nesting that does not pay for itself.** One longer readable
  block beats three short ones that force jumping around.
- **Full descriptive names carrying domain and unit information.** No
  abbreviations except tokens that behave as nouns (usb, tiff). No single-letter
  names, loop variables included. `WORKTREE_ROOT` not `ROOT`.
  `lane_pitch_millimetres` not `pitch`.
- **Comments state why**, not what: the constraint, the failure prevented, the
  alternative rejected, the non-obvious platform detail.
- **ASCII only.**
- **Verify by executing, not asserting.** Show the output.
- **Messages to stderr**, tagged with source, so stdout stays clean for piping.
  All messages also accumulate in memory and are serialized into the provenance
  JSON, so the run log is part of the record rather than lost scrollback.

---

## 5. Design decisions, with rejected alternatives

### 5.1 Grid model, not per-lane band search

**Chosen:** treat the gel as a grid. Lanes are columns found by collapsing the
image vertically into a 1D profile and fitting a *periodic* grid (pitch from
autocorrelation, then offset). Band rows are found by collapsing across all lanes
to find y positions where signal accumulates anywhere. Measure every
intersection, present or absent.

**Rejected:** find lanes, then find bands independently within each lane.

**Why:** the negative control lane is near-empty by design. Independent
per-lane band-finding returns fewer bands there, so "band 1" in the negative
control is not the same species as "band 1" in wildtype. The comparison silently
compares the wrong things, and the failure is invisible in the numbers.

**Consequence:** absent bands yield approximately zero, not a missing row. The
CSV is a complete rectangular table requiring no alignment logic downstream.

**Tradeoff named:** the grid model assumes bands run roughly horizontally and
that a species sits at comparable height in every lane. Smiling breaks this.
Mitigation (bounded per-lane row offset) deferred until we see whether real gels
need it.

### 5.2 Expected lane count supplied by the user

**The parameter is geometric, not experimental:** the number of comb positions
spanned by the cropped image, *counting skipped ones*. If the crop contains wells
1 through 13 and well 7 was left empty as a separator, the answer is 13.

**Chosen:** user supplies the count; the tool independently infers occupancy and
reports both. Disagreement is printed, not resolved.

**Rejected:** tool infers lane count entirely. A deliberate loading gap and a
genuinely empty negative control lane are indistinguishable to software, and
guessing risks silently dropping the negative control - the one lane whose
emptiness is the actual measurement.

**Practice note:** crop so the left edge sits in the gap before the first lane
and the right edge in the gap after the last. Makes the count unambiguous and
the grid offset fit better conditioned.

**Related decision:** gels containing two side-by-side experiments get **cropped
into two images and analyzed separately.** Their band rows may sit at different
heights, which would corrupt the global row-detection model into a union of both
row sets with phantom bands in each half. One extra crop removes the problem
entirely; block-aware row detection is not worth building.

### 5.3 Both integration modes, always, clearly labeled

**Region mode (default):** fixed row windows applied identically to every lane.
Integrate everything inside the box.

**Peak mode:** bounds are the valleys flanking a local maximum. Matches the
existing ImageJ wand workflow.

**Why region is primary:** the window is identical across lanes, so the same
vertical span is guaranteed to be compared in wildtype and mutant. No risk that
peak-finding chose slightly different bounds and manufactured a difference. Also
handles smeary signal (cleavage ladders, degradation products) correctly, and
handles the empty negative control naturally.

**Important asymmetry:** on an empty lane, region mode integrates its window and
returns approximately zero, which is correct and usable. Peak mode has no peak,
and the honest return is **NaN**, not zero - otherwise "measured, and empty" is
indistinguishable from "could not measure." Peak-mode cells with no detectable
peak emit NaN plus a `measurement_status` column stating why. The negative
control lane will likely be all NaN in peak mode and near-zero in region mode.
That is correct behaviour.

### 5.4 Dual baselines, both reported

**Ranked at design time:**

1. Valley-to-valley linear baseline. Score 8. Matches the existing manual
   workflow, so old and new numbers are comparable. Fails on badly overlapping
   peaks.
2. Rolling-ball / morphological opening on the 1D profile. Score 7. Robust to
   sloped background. Costs one magic number (radius) that changes results.
3. Asymmetric least squares. Score 6. Principled, two opaque parameters, hard to
   explain later.
4. Fitted Gaussian/pseudo-Voigt deconvolution. Score 4. Right answer for
   overlapping bands, substantial project. Not now.
5. Fixed constant from a blank region. Score 3. Wrong whenever illumination is
   uneven, which on gels it usually is.

**Chosen: 1 and 2 both, computed and written side by side.** If they agree within
a few percent the result is not a baseline artifact; if they diverge, the tool is
reporting that the measurement is fragile. Costs almost nothing and is the
strongest form of surfacing an assumption - not merely disclosing it but
quantifying what it is worth.

### 5.5 Modality is a mandatory argument, cross-checked

Three separable things are bundled in "modality," with different detectability:

- **Polarity** (bright-on-dark vs dark-on-light): **reliably detectable.** The
  histogram mode is the background level, since background dominates by pixel
  count. Mode near the bottom with a long upward tail means emission.
- **Encoding transform:** detectable from header/tags, not from pixels.
- **Stain chemistry** (Coomassie vs silver vs SYBR): **not detectable.** Matters
  because silver is non-linear with load in a way no transform fixes.

**Chosen:** mandatory command-line argument, no default. Script infers polarity
independently and cross-checks against the declaration. Agreement logged,
disagreement is a hard stop with both values printed. Catches the two realistic
mistakes: mislabeling a file, and a header field that lies.

Current work is phosphor and fluorescent - both emission, both linear, neither
needs the OD transform. That branch stays unwritten. Structure accommodates it
as one explicit branch; no framework.

### 5.6 Preprocessing writes a sidecar, not a new image

**Ranked:**

1. Preprocessing writes JSON sidecar with crop box and rotation angle; analysis
   reads the original file plus sidecar and applies both itself. Score 9.
2. Preprocessing writes a new cropped/rotated file. Score 6. Breaks the
   provenance chain, live risk of rotating twice.
3. Crop and rotation as command-line arguments each run. Score 5. Coordinates
   live in shell history rather than a file.

**Chosen: 1.** The measured array derives from the original in one traceable
step, crop and angle are recorded numbers in provenance, re-cropping means
editing four integers, and the input checksum stays the checksum of the
instrument's own output.

### 5.7 Rotation: tacit knowledge worth stating

Rotating by a non-multiple of 90 degrees requires **resampling** - the new pixel
grid does not align with the old, so every output pixel is interpolated from
several input pixels. This slightly blurs data and changes individual values. It
does not meaningfully change integrated band areas (interpolation approximately
conserves total signal), which is why it is acceptable, but it happens.

Consequences:

- **Rotate exactly once, from the original.** Never rotate an already-rotated
  file. Decision 5.6 enforces this structurally.
- **Bilinear (order 1) is the right default.** Nearest-neighbour preserves exact
  values but distorts geometry and jags band edges. Higher orders can overshoot
  and produce values outside the original range, breaking saturation checks.
- **Work in floating point** through measurement; convert to integer only for
  display images.

If the gel is within about half a degree of straight, skip rotation and record
that. Lane profiling tolerates small tilt (it just widens each lane's peak);
it only matters when tilt approaches a meaningful fraction of lane pitch. Stage 1
measures the tilt and reports whether correction is worth it.

### 5.8 Crop box acquisition

**Ranked:**

1. matplotlib TkAgg interactive selector. Score 8 on Windows 11 (WSLg gives
   working GUI with `sudo apt install python3-tk`), score 4 on Windows 10
   (needs VcXsrv plus DISPLAY fiddling).
2. Browser-based selector, HTML with drag-rectangle. Score 7. Sidesteps X
   entirely, costs hand-written JavaScript.
3. Read coordinates from ImageJ, pass on command line. Score 6. Zero
   infrastructure, works today, tedious, never breaks.
4. Headless auto-propose with preview and override. Score 5. The override path
   means building option 3 anyway.

**Chosen: build option 3's plumbing first** - crop box as four integers into the
sidecar JSON - and layer the selector on top later. The interactive picker is a
front end for the same four numbers. With the JSON contract in place first, the
GUI is a convenience that can be added or abandoned without touching anything
downstream, and the project is never blocked on WSL display issues.

Check `echo $WSL_DISTRO_NAME` and Windows version; on Win11 option 1 likely
works with no configuration.

### 5.9 Output directory: one per input file

**Chosen:** `<original_filename_stem>_gel_analysis/` created next to the source
file (or next to the working copy).

**Rejected:** outputs written directly into the input directory with a filename
prefix. Sets a trap for directory mode - a batch run sees its own previous
outputs as inputs. Prefix filtering solves it only until a prefix collides or
something gets renamed.

Inside the directory, fixed names. The directory already carries the identity, so
repeating the filename in every member is noise.

```
input_file_validation_report.json
input_file_histogram.png
provenance_and_configuration.json
run_log.txt
preprocessing_crop_and_rotation.json     (read, and copied in for the record)
overlay_lane_grid.png
overlay_band_boxes_labeled.png
lane_intensity_profiles.png
band_measurements.csv
```

Reruns overwrite directory contents, so there is never ambiguity about which
outputs belong to which run. Provenance JSON carries timestamp and input
checksum to prove it.

**Every image written for human eyes gets display scaling applied and is a
separate array from the measured data.** Contrast adjustment is display-only and
must never touch the analyzed array. The code should make this obvious.

### 5.10 CSV shape: long, not wide

**Chosen:** one row per cell **per method combination**, with
`integration_method` and `baseline_method` as columns.

**Rejected:** wide, with four intensity columns per cell row.

**Why:** downstream is R. Long is what R wants - `filter(integration_method ==
"region")` beats selecting the right column by name, and faceting a plot by
method is free. Adding a method later changes row counts, not column headers.

**Cost named:** a longer file, and a required filter step before anything else.
Annoying if opened in a spreadsheet.

Columns include at minimum: lane index, lane occupancy status, band row index,
integration method, baseline method, integrated intensity, background value used,
box coordinates, at-ceiling pixel count within the box, plateau statistic,
measurement status.

### 5.11 Which file to read

**Provisional ranking, to be settled empirically in stage 1:**

1. `.img` + `.inf` raw. Maximally direct: `numpy.fromfile`, reshape, flip.
   Plain-text header stating the parameters. Eliminates the TIFF library as a
   source of ambiguity. Now strongly supported by the exact filesize arithmetic.
   Risks: byte order, row order (documented as bottom-left), any header offset
   (arithmetic says none).
2. `.gel` via a TIFF library, with MD_ tags read and reported.
3. `.tif`. Convenient but its relationship to detector output is least
   specified. Do not build on it without evidence.

**Decisive cross-check, cheap, belongs in stage 1:** load the same scan from two
formats and plot pixel-for-pixel.

- Straight line through origin: identical data up to scale.
- Parabola: one is square-root encoded relative to the other, and `ScaleType`
  does not mean what we think.
- S-curve or flat ends: one is display-mapped and disqualified.

This also doubles as a **flip verification**. If two formats only correlate after
flipping one, the containers disagree in row order and this reveals which needed
it.

### 5.12 Saturation diagnostics

Whole-image saturation counting is insufficient - it says nothing about *where*.

**Per measured cell, emitted as CSV columns:**

- Count of at-ceiling pixels inside the box. Turns "clipping is only in regions
  I do not care about" from a background belief into a checkable per-row claim.
- A flatness/plateau statistic on the cell's 1D profile. A saturated band has a
  flat-topped profile rather than a peak, catching soft saturation that never
  reached the literal ceiling.

**The ceiling must be inferred, not assumed.** 16-bit TIFFs from scientific
cameras frequently carry 12-bit or 14-bit data in a 16-bit container. The ceiling
is then 4095 or 16383, not 65535, and a naive `== 65535` check reports zero
saturation on a badly clipped image. Infer from the data (maximum value, and
whether values are quantized to a stride - 12-bit left-aligned in 16 bits gives
multiples of 16) and **report what was concluded.**

### 5.13 Sensitivity analysis as a substitute for operator variance

Measuring how much two humans disagree is the correct bar but costs days the user
does not have.

**Cheap proxy:** run the same gel several times with deliberately perturbed
parameters - crop shifted a few pixels, rotation nudged half a degree, baseline
radius moved - and report the spread in final values.

This answers the real question: is this 2.3x, or is it "somewhere between 1.8x
and 2.9x depending on where the box was drawn." Small spread means the number is
real. Large spread means it was never there and no amount of manual care would
have found it.

Deferred to stage 4, but **config parameters should be trivially overridable from
the command line** so this is nearly free later.

---

## 6. Validation checklist

Nothing here is inferred from filename or file size. Every item is read and
reported.

### For `.inf` (hard stop on failure)

- Line 1 equals `FLA_IMAGE_FILE`
- Line 14 count matches the actual number of key-value lines following
- Bits per pixel is 16
- Width and height parse as positive integers
- `width * height * (bits/8) == .img filesize` **exactly**
- Human-readable timestamp and epoch timestamp agree
- `Orientation` field present and understood
- `ScaleType` present; anything other than `Linear` is a hard stop until the
  transform is implemented

### For TIFF containers (hard stop on failure)

- Opens as TIFF
- IFD (page) count; more than one requires explicit page selection, never
  default to page zero
- `SamplesPerPixel` must be 1; a 3 means RGB, which means something rendered it
- `BitsPerSample`; 8 means quantification is off the table
- `SampleFormat` (unsigned int vs float)
- `PhotometricInterpretation`, cross-checked against histogram-inferred polarity
- `Compression`; any lossy compression is a hard stop

### TIFF provenance (reported, never assumed)

- `ImageDescription`, `Software`, `Make`, `Model`, `DateTime`
- `XResolution`, `YResolution`, `ResolutionUnit`
- **All private tags 65000+ dumped verbatim** (MD_FileTag, MD_ScalePixel, etc.)
- Full tag inventory, so nothing is silently ignored

### Data statistics (computed, plotted)

- Dimensions and estimated memory footprint
- Minimum, maximum, mean, median, standard deviation
- **Histogram plotted on log-count y-axis.** Band pixels are a tiny minority;
  linear y hides them completely.
- Effective bit depth inference and the resulting ceiling
- Count and fraction of pixels at that ceiling
- Count and fraction at the floor. A large floor population means the image was
  already background-subtracted or clipped at zero, which changes what baselines
  mean.
- **Isolated extreme pixel count** - single-pixel spikes with no bright
  neighbours. On phosphor screens these come from cosmic rays and screen
  contamination; on any camera from hot pixels. **Reported, not silently
  filtered** - the user decides.
- Estimated lane tilt angle

---

## 7. Staging plan

Each stage ends at something independently verifiable by running it, not at a
code-organization boundary.

### Stage 0 - skeleton, environment, path handling

`prepare_gel_image` script plus the analysis script's argument handling. Path
normalization: WSL paths, Windows paths with backslashes and drive letters,
spaces, commas, square brackets, `~`, relative paths. Existence and readability
checks. Output directory creation.

Ends by printing resolved absolute paths and exiting. No image reading.

Verify against: a WSL path, a pasted Windows path, a path with spaces and
brackets, a nonexistent file, a directory.

### Stage 1 - file interrogation and validation

Everything in section 6. Emits `input_file_validation_report.json` and
`input_file_histogram.png`. Plus the `.img`-versus-`.gel` scatter cross-check
from 5.11.

No analysis, no lanes.

**Run this on all 4-6 real files before writing another line.** It is valuable
standalone: it answers "are my images actually what I think they are," currently
an unverified assumption carried across every gel.

### Stage 2 - config, provenance, geometry

Config block as flat data at the top of the file, lengths in millimetres.
Provenance serialization: git commit and dirty flag for `lab_utils`, script
hash, input file SHA-256, full `.inf` parsed verbatim into JSON, all config
values, package versions, timestamps, the accumulated run log.

Then polarity inference and cross-check, orientation flip, rotation and crop from
the sidecar, lane grid fit.

Ends by writing `overlay_lane_grid.png` on a display-scaled copy. Verify by eye.

### Stage 3 - band rows and measurement

Row detection, per-cell integration in region and peak modes, dual baselines,
saturation and plateau statistics per cell.

Emits `overlay_band_boxes_labeled.png`, `lane_intensity_profiles.png`,
`band_measurements.csv`.

### Stage 4 - directory mode and hardening

Glob a directory following a filename convention, skip already-generated output
directories, aggregate. Sensitivity analysis.

---

## 8. Open questions

### Blocking nothing, but worth resolving early

- Whether ImageJ silently applies a transform when opening `.gel` files. If it
  does, and these files are already linear, historical manual measurements may
  not be comparable to new pipeline output. Worth knowing **before** comparing
  old and new numbers and concluding something is broken.
- Meaning of `RangeLow` / `RangeHigh` and positional lines 8 and 9.
- Byte order of `.img`. Almost certainly little-endian; verifiable by checking
  whether the histogram looks sane or like noise.

### Requires user input

- Confirm `.gitignore` already covers Python and uv artifacts.
- Windows version, for the matplotlib GUI question.
- Whether outputs should live next to the source data, in the repo, or in a
  local scratch directory.

### Deferred by decision

- Per-lane vertical offset tracking for smiling
- OD transform for transmitted-light stains
- Multi-page TIFF handling
