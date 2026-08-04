# Measurement protocol - gel densitometry pipeline

Operator instructions. Companion to `DESIGN.md`, which records why; this records
what to type. Self-contained: everything here can be followed without reading the
design.

Run every command from `~/personal_repos/lab_utils/images`.

**Quote every path with single quotes.** Real filenames contain spaces, commas and
square brackets, and square brackets are glob metacharacters in bash.

---

## 0. The one thing that will catch you out

**The image opens upside down.** The scan origin is bottom-left, so row 0 of the
file is the physical bottom of the gel, and Fiji draws row 0 at the top of the
window. On screen the gel is mirrored vertically, which means:

> **The wells appear at the BOTTOM of the Fiji window, not the top.**

Everything you measure for the sidecar is measured in this as-opened frame. Do not
flip it, do not rotate it, do not save it. The pipeline handles the flip.

Use this as a free check on the metadata: whichever end of the window the wells
are at tells you the row order directly, without trusting any tag. If the wells
are at the top instead, the file is not what the pipeline believes it is. Stop and
say so.

Second thing, less dangerous but confusing: Fiji renders these files with an
inverting lookup table, because the TIFF declares `PhotometricInterpretation =
MINISWHITE`. Bands look dark on light in Fiji. Every image this pipeline emits
uses a matching reversed colormap, so the two agree. If you ever see an image from
this pipeline with bright bands on a dark field, something is wrong with the
pipeline, not with your file.

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

---

## 2. Inventory the tags, once per new instrument or new file type

Reports everything, interprets nothing, validates nothing. Run it when you meet a
file type for the first time, or when something looks wrong.

```
uv run inventory_tiff_tags.py '/mnt/c/Users/liusm/MIT Dropbox/Luis Martinez/lab-backup/experiments/<experiment>/<scan>.tif'
```

stdout is the inventory and stderr is the commentary, so this redirects cleanly:

```
uv run inventory_tiff_tags.py '<path>.tif' > tags.txt
```

---

## 3. Create the output directory

```
uv run prepare_gel_image.py '<path>.tif'
```

Prints resolved paths and exits, reading no image data. It creates
`<stem>_gel_analysis` beside the input.

The output directory is **per scan, not per file**: the `.img` and the `.tif` of
one scan share a stem, so both resolve to the same directory. That is intended.
If you already ran this against the `.img` at some point, the directory exists and
you will see a warning that it is not empty. That is expected.

The parent directory must already exist. Only the leaf is created, so a typo in
`--output-parent-directory` fails instead of building a tree somewhere wrong.

---

## 4. Measure the sidecar in Fiji

Once per scan. This is the only manual step.

1. Copy the template:

   ```
   cp preprocess_sidecar_template.txt '/mnt/c/.../<stem>_preprocess.txt'
   ```

2. Open the `.tif` in Fiji. Drag and drop is fine. **Do not** adjust
   brightness/contrast and press Apply, **do not** convert to 8-bit, **do not**
   rotate, **do not** save. Any of those either destroys the data or discards the
   tags. Fiji is being used here purely as a ruler.

3. Read the dimensions from `Image > Properties` and fill in
   `measured_against_image_width_pixels` and
   `measured_against_image_height_pixels`. Fill in
   `measured_against_input_filename` with the filename only, no directory.

4. **Landmarks for the tilt.** Find something that should be level on the
   physical gel: the line of wells is usually best, or a band shared by the
   outermost lanes. Hover over the leftmost instance and read `x=..., y=...` from
   the status bar; record it as the `left` pair. Hover over the rightmost
   instance and record it as the `right` pair.

   **Choose the two points as far apart as possible.** The angle is derived from
   them, so the span sets the precision:

   | Span between landmarks | Uncertainty from a 1 px click error |
   |---|---|
   | 400 px | 0.14 degrees |
   | 800 px | 0.07 degrees |
   | 1000 px | 0.06 degrees |

   For scale, on a 1125 px wide image a residual tilt of 0.5 degrees shifts a
   band by about 10 px end to end, which is a full band height at 200
   micrometres per pixel. 0.1 degrees shifts it by 2 px, which does not matter.

   Write down what you clicked in `rotation_landmark_description`. It is the only
   record.

5. **Crop box.** Drag a rectangle with the rectangle tool and read `x`, `y`, `w`,
   `h` from the status bar.

   The crop is **not** a measurement boundary and does not need to be precise. Be
   generous. Exclude:

   - the wells, which are a near-saturated step artifact that a rolling-ball
     baseline treats as a giant feature, bleeding the correction into your
     topmost real band. Remember they are at the **bottom** of the window.
   - everything outside the gel edge, which is plate background: a different
     pixel population that skews the histogram, the inferred floor and the
     baseline.

   Keep a margin of clean gel, roughly 20 to 50 px, outside the outermost lanes
   and beyond the topmost and bottommost bands. The baseline needs somewhere to
   be measured. A generous crop costs a slightly slower run; a tight one biases
   every baseline in the image.

   The pipeline crops first and rotates second, and rotation leaves small blank
   wedges in the corners, about `(width / 2) * tan(angle)`. On a 910 px crop at 1
   degree that is 8 px. Your margin absorbs it, which is one more reason for the
   margin.

6. Fill in `expected_lane_count`, counting empty lanes.

7. Save the file. Leave `notes` empty if you have nothing to say.

---

## 5. Validate

```
uv run validate_gel_image.py '<path>.tif'
```

Writes into `<stem>_gel_analysis`:

- `input_file_validation_report.json` - every check, every tag, every statistic
- `input_file_histogram.png` - log-count y axis
- `input_file_crop_preview.png` - only if the sidecar exists

Exit code is 0 if every check passed and 1 if any hard stop failed. The report is
written either way, so a failure is diagnosable from the file and not only from
the scrollback.

**Look at `input_file_crop_preview.png` before going further.** It draws your crop
box and your two landmarks back onto a display-scaled copy. This is the only check
on numbers you typed by hand, and a transposed digit is obvious in the picture and
invisible in the text.

The sidecar is optional at this stage. Without it you still get the full tag
inventory and every statistic; you just get no preview and no tilt angle.

---

## 6. What the report will not tell you

**Whether the pixel values are linear in signal.** For the Typhoon files this is
unresolved. Vendor documentation calls the `.tif` linear 16-bit grayscale and it
is the one container no source claims is transformed, but that has not been
verified against these files, and no cheap test settles it. The report carries
`encoding_verified: false` and names the test that would change it. If the `.tif`
is not linear, every ratio this pipeline produces is compressed toward 1 and
nothing in the pipeline will notice.

For the Imager 680 fluorescence files the position is much stronger: a CCD cooled
to -25 C is linear in photons by construction, and no transform is claimed.
**But** those files record an exposure time, so two fluorescence images are only
comparable after dividing by it. The report extracts and records the exposure.

**Whether the gel is tilted more than you think.** Stage 1 reports the angle
derived from your landmarks and does not attempt to find the tilt itself.

**Whether a band is saturated.** Stage 1 reports whole-image ceiling statistics
only. Per-band saturation and the plateau statistic come later, and the plateau
statistic is the primary detector because the instrument can saturate below the
numeric ceiling.
