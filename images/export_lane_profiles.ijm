// export_lane_profiles.ijm
// ---------------------------------------------------------------------------
// Manual densitometry bridge for the gel pipeline.
//
// What it does, in one pass, on the CURRENTLY OPEN rotated 16-bit image:
//   1) Saves the lane rectangles you added to the ROI Manager as a RoiSet.zip
//      (replayable: roiManager("Open", path) restores the exact rectangles).
//   2) Exports a raw per-lane migration profile CSV whose raw_value matches
//      measure_gel.py: signal = pixel - median(whole image), clipped at 0,
//      then summed across the lane width at each migration position.
//
// The CSV is intentionally raw. detrended_value and band integration are done
// in Python so that math stays identical to measure_gel.py's own functions.
//
// PREP (do this before running):
//   - Open the rotated .tif. Do NOT invert. Do NOT apply "Invert Peaks".
//   - Draw each lane with the rectangle tool; press  t  to add it to the
//     ROI Manager. Order does not matter; the macro records both drawn order
//     and a left-to-right (or top-to-bottom) position order.
//   - Then run this macro (Plugins > Macros > Run...).
//
// ASSUMPTIONS:
//   - LANES_ARE_VERTICAL = true means migration runs down the image (Y axis)
//     and lanes are side by side across X. Set false for horizontal gels.
//   - Raw pixel values are bright-on-dark (bands = high values). Confirmed for
//     these MinIsWhite scans: the LUT inverts display only, not the data.
// ---------------------------------------------------------------------------

// -------------------- operator settings --------------------
LANES_ARE_VERTICAL = true;   // migration axis = Y when true, X when false
// -----------------------------------------------------------

if (nImages == 0) exit("No image open. Open the rotated .tif first.");
n_rois = roiManager("count");
if (n_rois == 0)
    exit("ROI Manager is empty. Draw each lane and press t to add it, then re-run.");

image_title = getTitle();

// Output directory: reuse the pipeline's <stem>_gel_analysis folder next to the
// image if it exists, otherwise create it, so files land where the rest of the
// pipeline expects them.
image_dir = getDirectory("image");
if (image_dir == "")
    exit("This image has no directory on disk. Save it first, then re-run.");
stem = image_title;
stem = replace(stem, "\\.tif$", "");
stem = replace(stem, "\\.tiff$", "");
out_dir = image_dir + stem + "_gel_analysis" + File.separator;
if (!File.exists(out_dir)) File.makeDirectory(out_dir);

// Plate background = whole-image median, matching measure_gel.py line 774.
run("Select None");
List.setMeasurements;
plate_bg = List.getValue("Median");

// Pixel size for the millimetre column. getVoxelSize returns size per pixel in
// the image's own unit; convert to mm. If the scale was zeroed to pixels, leave
// mm blank rather than emit a wrong number.
getVoxelSize(vx_w, vx_h, vx_d, vx_unit);
mm_per_pixel_migration = 0;   // 0 means "unknown / not set"
axis_pixel_size = vx_h;                       // migration = Y -> use pixel height
if (!LANES_ARE_VERTICAL) axis_pixel_size = vx_w;  // migration = X -> pixel width
u = toLowerCase(vx_unit);
if (u == "cm" || u == "centimeter" || u == "centimetre")
    mm_per_pixel_migration = axis_pixel_size * 10.0;
else if (u == "mm" || u == "millimeter" || u == "millimetre")
    mm_per_pixel_migration = axis_pixel_size;
else if (u == "um" || u == "micron" || u == "microns" || u == "micrometer" || u == "micrometre")
    mm_per_pixel_migration = axis_pixel_size / 1000.0;
// any other unit (e.g. "pixel") leaves mm_per_pixel_migration = 0

// Determine position order of lanes so a stable lane_index can be assigned.
// Vertical lanes: order by ROI centre X (left to right).
// Horizontal lanes: order by ROI centre Y (top to bottom).
centre_key = newArray(n_rois);
for (i = 0; i < n_rois; i++) {
    roiManager("select", i);
    Roi.getBounds(rx, ry, rw, rh);
    if (LANES_ARE_VERTICAL) centre_key[i] = rx + rw / 2.0;
    else                    centre_key[i] = ry + rh / 2.0;
}
// rank[i] = number of ROIs whose centre_key is strictly less (ties broken by index)
position_rank = newArray(n_rois);
for (i = 0; i < n_rois; i++) {
    rank = 0;
    for (j = 0; j < n_rois; j++) {
        if (centre_key[j] < centre_key[i]) rank++;
        else if (centre_key[j] == centre_key[i] && j < i) rank++;
    }
    position_rank[i] = rank;   // 0-based; lane_index below is 1-based
}

// Build the CSV. One row per (lane, migration position).
header = "lane_index,drawn_order,roi_name,lane_detection_status,prediction_span,"
       + "roi_x,roi_y,roi_w,roi_h,plate_background_median,"
       + "migration_position_pixels,migration_position_millimetres,raw_value";
csv = header + "\n";

setBatchMode(true);
for (i = 0; i < n_rois; i++) {
    roiManager("select", i);
    roi_name = Roi.getName;
    if (roi_name == "") roi_name = "lane_" + IJ.pad(i + 1, 3);
    Roi.getBounds(rx, ry, rw, rh);
    lane_index = position_rank[i] + 1;   // 1-based, position-sorted
    drawn_order = i + 1;                  // 1-based, order added to manager

    if (LANES_ARE_VERTICAL) {
        // migration = Y: one value per row, summed across the lane width (X)
        n_pos = rh;
        for (p = 0; p < n_pos; p++) {
            y = ry + p;
            row_sum = 0.0;
            for (x = rx; x < rx + rw; x++) {
                v = getPixel(x, y) - plate_bg;
                if (v < 0) v = 0;
                row_sum += v;
            }
            csv = csv + line_for(lane_index, drawn_order, roi_name, rx, ry, rw, rh,
                                 plate_bg, p, mm_per_pixel_migration, row_sum);
        }
    } else {
        // migration = X: one value per column, summed across the lane height (Y)
        n_pos = rw;
        for (p = 0; p < n_pos; p++) {
            x = rx + p;
            col_sum = 0.0;
            for (y = ry; y < ry + rh; y++) {
                v = getPixel(x, y) - plate_bg;
                if (v < 0) v = 0;
                col_sum += v;
            }
            csv = csv + line_for(lane_index, drawn_order, roi_name, rx, ry, rw, rh,
                                 plate_bg, p, mm_per_pixel_migration, col_sum);
        }
    }
    showProgress(i + 1, n_rois);
}
setBatchMode(false);

// Write the profile CSV.
profile_path = out_dir + "manual_lane_profiles.csv";
File.saveString(csv, profile_path);

// Save the exact rectangles for replay.
roi_zip_path = out_dir + "manual_lane_rois.zip";
roiManager("Deselect");
roiManager("Save", roi_zip_path);

// Small provenance sidecar so the export is self-describing.
prov = "image_title\t" + image_title + "\n"
     + "image_dir\t" + image_dir + "\n"
     + "lane_count\t" + n_rois + "\n"
     + "lanes_are_vertical\t" + LANES_ARE_VERTICAL + "\n"
     + "plate_background_median\t" + plate_bg + "\n"
     + "voxel_unit\t" + vx_unit + "\n"
     + "mm_per_pixel_migration\t" + mm_per_pixel_migration + "\n"
     + "profile_csv\t" + profile_path + "\n"
     + "roi_zip\t" + roi_zip_path + "\n";
File.saveString(prov, out_dir + "manual_export_provenance.txt");

print("[export] wrote " + profile_path);
print("[export] wrote " + roi_zip_path);
print("[export] plate background median = " + plate_bg
      + " ; mm/pixel = " + mm_per_pixel_migration
      + " ; lanes = " + n_rois);

// Helper: format one CSV row. mm column left blank when scale is unknown.
function line_for(lane_index, drawn_order, roi_name, rx, ry, rw, rh,
                  plate_bg, pos, mm_per_pixel, value) {
    mm_text = "";
    if (mm_per_pixel > 0) mm_text = d2s(pos * mm_per_pixel, 6);
    return lane_index + "," + drawn_order + "," + roi_name
         + ",manual,manual,"
         + rx + "," + ry + "," + rw + "," + rh + ","
         + d2s(plate_bg, 6) + ","
         + pos + "," + mm_text + "," + d2s(value, 6) + "\n";
}
