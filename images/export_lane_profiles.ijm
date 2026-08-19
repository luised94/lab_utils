// export_lane_profiles.ijm
// ---------------------------------------------------------------------------
// Manual densitometry bridge for the gel pipeline. Runs on the CURRENTLY OPEN
// rotated 16-bit image and, in one pass:
//   1) Saves the lane rectangles from the ROI Manager as a RoiSet zip that
//      restores the exact rectangles later (replay for the analysis).
//   2) Writes a raw per-lane migration profile CSV whose raw_value uses the
//      same background convention as measure_gel.py: signal = pixel minus the
//      whole-image median, negatives clipped to zero, summed across the lane
//      width at each migration position.
//
// detrended_value and band integration are intentionally NOT done here. They
// stay in Python so that math matches measure_gel.py's own functions exactly.
//
// PREP before running:
//   - Open the rotated .tif. Do not invert, do not apply Invert Peaks. On these
//     MinIsWhite scans the inverting LUT is display only: getPixel returns the
//     raw stored value, which is high at bands, which is what the sum needs.
//   - Draw each lane with the rectangle tool and press t to add it to the ROI
//     Manager. Draw order does not matter; the CSV records both the draw order
//     and a position-sorted lane index.
//   - Then run this macro (Plugins > Macros > Run...).
// ---------------------------------------------------------------------------

// A number-first "+" chain is evaluated as arithmetic by the macro interpreter
// and collapses to NaN, so every string built below starts from "" to force
// string concatenation, and rows accumulate in the String buffer for speed on
// images with tens of thousands of migration rows.

MIGRATION_AXIS_IS_VERTICAL = true;   // true: migration runs down (Y), lanes across X
OUTPUT_SUBDIRECTORY_SUFFIX = "_gel_analysis";
PROFILE_CSV_FILENAME = "manual_lane_profiles.csv";
ROI_ZIP_FILENAME = "manual_lane_rois.zip";
PROVENANCE_FILENAME = "manual_export_provenance.txt";
CENTIMETRES_TO_MILLIMETRES = 10.0;
MICROMETRES_PER_MILLIMETRE = 1000.0;
// Standard comb sizes offered in the up-front dialog, plus a custom escape hatch.
// The chosen value is stamped into every row so the export records the comb the
// operator intended, and is asserted against the ROI count below.
COMB_WELL_COUNT_CHOICES = newArray("15", "12", "10", "other");
CSV_HEADER = "lane_index,drawn_order,roi_name,lane_detection_status,prediction_span,"
           + "comb_well_count,"
           + "roi_x,roi_y,roi_w,roi_h,plate_background_median,"
           + "migration_position_pixels,migration_position_millimetres,raw_value";

if (nImages == 0)
    exit("No image open. Open the rotated .tif first.");
lane_roi_count = roiManager("count");
if (lane_roi_count == 0)
    exit("ROI Manager is empty. Draw each lane, press t to add it, then re-run.");

// Ask the comb size once, up front, before any pixels are read. Making it an
// explicit choice (rather than inferring it from the ROI count) is what lets the
// assertion below catch a missed or extra box instead of silently trusting it.
Dialog.create("Gel comb");
Dialog.addChoice("Wells in the comb:", COMB_WELL_COUNT_CHOICES, COMB_WELL_COUNT_CHOICES[0]);
Dialog.addNumber("If 'other', wells:", lane_roi_count);
Dialog.show();
comb_well_count_choice = Dialog.getChoice();
comb_well_count_custom = Dialog.getNumber();
if (comb_well_count_choice == "other")
    comb_well_count = comb_well_count_custom;
else
    comb_well_count = parseInt(comb_well_count_choice);

// Draw every comb position, loaded or not, so lane_index maps one-to-one onto
// well number across gels. A mismatch here means a box was missed or doubled.
if (lane_roi_count != comb_well_count)
    exit("Drew " + lane_roi_count + " lanes but the comb has " + comb_well_count
       + ". Add or remove ROIs so every well is drawn exactly once, then re-run.");

image_title = getTitle();

// Reuse the pipeline's <stem>_gel_analysis folder next to the image if it exists
// (prepare_gel_image.py makes it), otherwise create it, so files land where the
// rest of the pipeline expects them.
image_directory = getDirectory("image");
if (image_directory == "")
    exit("This image has no directory on disk. Save it first, then re-run.");
image_stem = image_title;
image_stem = replace(image_stem, "\\.tif$", "");
image_stem = replace(image_stem, "\\.tiff$", "");
output_directory = image_directory + image_stem + OUTPUT_SUBDIRECTORY_SUFFIX + File.separator;
if (!File.exists(output_directory))
    File.makeDirectory(output_directory);

// Plate background is the whole-image median (matches measure_gel.py). Select
// None first so the statistic covers the whole image, not a leftover selection.
run("Select None");
List.setMeasurements;
plate_background_median = List.getValue("Median");

// Pixel size for the millimetre column, taken from the image calibration. If the
// scale was zeroed to pixels the unit is not a length, so the mm column is left
// blank rather than filled with a wrong number.
getVoxelSize(voxel_width, voxel_height, voxel_depth, voxel_unit);
if (MIGRATION_AXIS_IS_VERTICAL)
    pixel_size_along_migration_axis = voxel_height;
else
    pixel_size_along_migration_axis = voxel_width;
millimetres_per_migration_pixel = 0;   // 0 means unknown / not calibrated
voxel_unit_lowercase = toLowerCase(voxel_unit);
if (voxel_unit_lowercase == "cm" || voxel_unit_lowercase == "centimeter" || voxel_unit_lowercase == "centimetre")
    millimetres_per_migration_pixel = pixel_size_along_migration_axis * CENTIMETRES_TO_MILLIMETRES;
else if (voxel_unit_lowercase == "mm" || voxel_unit_lowercase == "millimeter" || voxel_unit_lowercase == "millimetre")
    millimetres_per_migration_pixel = pixel_size_along_migration_axis;
else if (voxel_unit_lowercase == "um" || voxel_unit_lowercase == "micron" || voxel_unit_lowercase == "microns" || voxel_unit_lowercase == "micrometer" || voxel_unit_lowercase == "micrometre")
    millimetres_per_migration_pixel = pixel_size_along_migration_axis / MICROMETRES_PER_MILLIMETRE;

// A stable lane_index that does not depend on draw order: rank each ROI by the
// centre coordinate across the lanes (X for vertical lanes, Y for horizontal),
// so lane 1 is the physically first lane regardless of the order it was added.
lane_centre_coordinate = newArray(lane_roi_count);
for (ranked_roi_index = 0; ranked_roi_index < lane_roi_count; ranked_roi_index++) {
    roiManager("select", ranked_roi_index);
    Roi.getBounds(ranked_roi_x, ranked_roi_y, ranked_roi_width, ranked_roi_height);
    if (MIGRATION_AXIS_IS_VERTICAL)
        lane_centre_coordinate[ranked_roi_index] = ranked_roi_x + ranked_roi_width / 2.0;
    else
        lane_centre_coordinate[ranked_roi_index] = ranked_roi_y + ranked_roi_height / 2.0;
}
lane_index_by_roi = newArray(lane_roi_count);
for (outer_roi_index = 0; outer_roi_index < lane_roi_count; outer_roi_index++) {
    lanes_positioned_before = 0;
    for (inner_roi_index = 0; inner_roi_index < lane_roi_count; inner_roi_index++) {
        if (lane_centre_coordinate[inner_roi_index] < lane_centre_coordinate[outer_roi_index])
            lanes_positioned_before++;
        // Break centre-coordinate ties by draw order so two lanes never share an index.
        else if (lane_centre_coordinate[inner_roi_index] == lane_centre_coordinate[outer_roi_index] && inner_roi_index < outer_roi_index)
            lanes_positioned_before++;
    }
    lane_index_by_roi[outer_roi_index] = lanes_positioned_before + 1;   // 1-based
}

String.resetBuffer;
String.append(CSV_HEADER + "\n");

setBatchMode(true);
for (lane_roi_index = 0; lane_roi_index < lane_roi_count; lane_roi_index++) {
    roiManager("select", lane_roi_index);
    lane_roi_name = Roi.getName;
    if (lane_roi_name == "")
        lane_roi_name = "lane_" + IJ.pad(lane_roi_index + 1, 3);
    Roi.getBounds(roi_bounds_x, roi_bounds_y, roi_bounds_width, roi_bounds_height);
    lane_index_for_this_roi = lane_index_by_roi[lane_roi_index];
    drawn_order_number = lane_roi_index + 1;

    if (MIGRATION_AXIS_IS_VERTICAL) {
        migration_axis_length = roi_bounds_height;
        width_axis_length = roi_bounds_width;
    } else {
        migration_axis_length = roi_bounds_width;
        width_axis_length = roi_bounds_height;
    }

    for (migration_position_index = 0; migration_position_index < migration_axis_length; migration_position_index++) {
        summed_signal_across_width = 0.0;
        for (width_column_index = 0; width_column_index < width_axis_length; width_column_index++) {
            if (MIGRATION_AXIS_IS_VERTICAL) {
                current_pixel_x = roi_bounds_x + width_column_index;
                current_pixel_y = roi_bounds_y + migration_position_index;
            } else {
                current_pixel_x = roi_bounds_x + migration_position_index;
                current_pixel_y = roi_bounds_y + width_column_index;
            }
            background_subtracted_value = getPixel(current_pixel_x, current_pixel_y) - plate_background_median;
            // Clip negatives so plate-noise dips below the median cannot subtract
            // from real band signal; matches measure_gel.py signal_above_plate.
            if (background_subtracted_value < 0)
                background_subtracted_value = 0;
            summed_signal_across_width += background_subtracted_value;
        }

        migration_millimetres_text = "";
        if (millimetres_per_migration_pixel > 0)
            migration_millimetres_text = d2s(migration_position_index * millimetres_per_migration_pixel, 6);

        // Leading "" forces string concatenation over the whole chain; without it
        // the interpreter reads lane_index_for_this_roi + "," as arithmetic (NaN).
        output_row_text = "" + lane_index_for_this_roi
                        + "," + drawn_order_number
                        + "," + lane_roi_name
                        + ",manual,manual,"
                        + comb_well_count + ","
                        + roi_bounds_x + "," + roi_bounds_y
                        + "," + roi_bounds_width + "," + roi_bounds_height
                        + "," + d2s(plate_background_median, 6)
                        + "," + migration_position_index
                        + "," + migration_millimetres_text
                        + "," + d2s(summed_signal_across_width, 6);
        String.append(output_row_text + "\n");
    }
    showProgress(lane_roi_index + 1, lane_roi_count);
}
setBatchMode(false);

profile_csv_path = output_directory + PROFILE_CSV_FILENAME;
File.saveString(String.buffer, profile_csv_path);

// Save the exact rectangles for replay. Deselect first so the whole set is saved
// rather than only the single selected ROI.
roi_zip_path = output_directory + ROI_ZIP_FILENAME;
roiManager("Deselect");
roiManager("Save", roi_zip_path);

provenance_text = "" + "image_title\t" + image_title + "\n"
                + "image_directory\t" + image_directory + "\n"
                + "lane_count\t" + lane_roi_count + "\n"
                + "comb_well_count\t" + comb_well_count + "\n"
                + "migration_axis_is_vertical\t" + MIGRATION_AXIS_IS_VERTICAL + "\n"
                + "plate_background_median\t" + plate_background_median + "\n"
                + "voxel_unit\t" + voxel_unit + "\n"
                + "millimetres_per_migration_pixel\t" + millimetres_per_migration_pixel + "\n"
                + "profile_csv\t" + profile_csv_path + "\n"
                + "roi_zip\t" + roi_zip_path + "\n";
File.saveString(provenance_text, output_directory + PROVENANCE_FILENAME);

print("[export] wrote " + profile_csv_path);
print("[export] wrote " + roi_zip_path);
print("[export] plate background median = " + plate_background_median
      + " ; mm/pixel = " + millimetres_per_migration_pixel
      + " ; lanes = " + lane_roi_count);
