# Measurement protocol - gel densitometry pipeline

Operator instructions. Companion to `DESIGN.md`, which records why; this records
what to type. Self-contained: everything here can be followed without reading the
design.

Run every command from `~/personal_repos/lab_utils/images`.

**Quote every path with single quotes.** Real filenames contain spaces, commas,
square brackets, and characters like `& ; | ( ) $ ` `` ` `` `! * ? # ~`. Single
quotes neutralize all of them. Note that inside *double* quotes an interactive
bash still history-expands `!`, which is one more reason to use single quotes.

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

Two gotchas. If no `.tif` matches, bash prints the literal pattern `.../*.tif`
rather than nothing, so an empty directory is obvious but looks odd. And the
printed paths still need single-quoting when you paste them back into a command.
To capture one path in a variable that survives a `cd`, quote it once at
assignment:

```
GEL="$(pwd)/<name>.tif"; ls -la "$GEL"
```

---

## 2. Inventory the tags, once per new instrument or new file type

Reports everything, interprets nothing, validates nothing. Run it when you meet a
file type for the first time, or when something looks wrong.

```
uv run inventory_tiff_tags.py '<experiment directory>/<scan>.tif'
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

Once per scan, the only manual step. Copy the template to `<stem>_preprocess.txt`
next to the `.tif` first. Open the `.tif` in Fiji (drag and drop). Do NOT adjust
brightness/contrast and Apply, convert to 8-bit, rotate, or save; any of those
destroys data or discards tags.

Do these in order. The order matters: the line macro in step 4 only works while a
line is the active selection, so the rectangle must come after it.

1. **Preflight, once per session.** `Analyze > Set Measurements`: tick "Bounding
   rectangle" (so Measure reports BX, BY, Width, Height), and confirm "Invert Y
   coordinates" is UNTICKED. If it is ticked, every y you record is silently
   mirrored. Confirm the origin is 0,0.

2. **Record identity.** Fill `measured_against_input_filename` (filename only),
   `measured_against_image_width_pixels` and `_height_pixels` (from
   `Image > Properties` or the title bar), `gel_migration_axis` (from section 0),
   and `coordinate_unit` (see the template's how-to-tell note).

3. **Draw the landmark line.** Use the well line if you can see it, else the
   cassette cast edge. With the straight-line tool, draw from one far end to the
   other, as far apart as possible. Do NOT hold Shift: it snaps the line
   perfectly straight, forcing the sideways drift to zero and reporting a tilt of
   exactly zero that looks like success but has destroyed the measurement. Zoom to
   100% at each end before placing the point.

4. **Measure and capture the line.** `Ctrl+M`, and record the `Angle` value (an
   independent cross-check on the tilt sign). Then, with the line still selected,
   `Plugins > New > Macro`, paste and run (`Ctrl+R`):

   ```
   getLine(x1, y1, x2, y2, w); print(x1+", "+y1+", "+x2+", "+y2);
   ```

   `getLine` returns the endpoints in raw pixels regardless of any centimetre
   calibration. Record them as `landmark_a` and `landmark_b`, and note in
   `rotation_landmark_description` what you clicked and which end is a.

5. **Draw the crop rectangle and `Ctrl+M`.** Follow the crop rule in the template:
   contain every band plus 20 to 50 px of clean gel on all four sides; do not clip
   a band or its flanking gel; include or exclude the pure well strip per that
   rule. Glance at the reported Max: below 65535 means no saturated
   (ceiling-clipped) pixel is inside the crop.

6. **Finish.** Fill `expected_lane_count` (empty lanes included), save, leave
   `notes` empty if you have nothing to add.

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

---

## Appendix. Quick tests in Fiji

Cheap checks that each catch one specific silent failure. None of them alter the
file.

**Saturation-location test.** `Image > Adjust > Threshold` (`Ctrl+Shift+T`). Drag
the lower slider up to 65535 and leave the upper at 65535, so only ceiling pixels
(the maximum a 16-bit pixel can hold) are highlighted in red. Do NOT press Apply.
Look at where the red sits: if any is inside your crop, those lanes are
saturated and their numbers are compromised; if it is all outside, the measured
lanes are clean. This is the direct test of whether at-ceiling pixels touch your
measurement region.

**Whole-image read.** `Ctrl+Shift+A` to deselect, `Ctrl+A` to select all, `Ctrl+M`.
The Min, Max and Mean now describe every pixel. On these gels the majority
population is bare plate, so the image-wide mode and median describe the plate,
not the gel; the gel's own statistics come from inside the crop.

**Exact line endpoints.** With a line selected, `Plugins > New > Macro`, run
`getLine(x1,y1,x2,y2,w); print(x1+", "+y1+", "+x2+", "+y2);`. Returns raw pixels
with no centimetre rounding. If it prints `-1, -1, -1, -1`, a line is not the
active selection (you probably measured a rectangle after the line); redraw the
line and run it again.

**Centimetre-versus-pixel check.** Hover any pixel and read the status bar. A
value in parentheses, like `x=3.80 (476)`, means the image is calibrated: the
leading number is centimetres, the parenthesised number is the pixel. A bare
integer means you are already in pixels. This decides `coordinate_unit` and
prevents recording centimetres into a field read as pixels.
