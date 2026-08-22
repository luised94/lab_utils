// distribute_selection.ijm
//
// Draw one rectangular ROI (a lane), then tile evenly-spaced copies of it. Built for
// laying out gel lanes: draw the first lane, get the rest at a fixed pitch, add them
// to the ROI Manager, and hand the set to export_lane_profiles.ijm.
//
// Layout is Row (lanes left-to-right, the gel case) or Column (top-to-bottom). Grid
// (2D tiling) was removed: gel lanes are a single row, so Grid had no real use here,
// and its perfect-square case placed a copy on top of another. Row and Column are the
// two real call sites; nothing tiles in 2D.
//
// Conventions follow the repo's Python style where the macro language allows: full
// descriptive names, no abbreviations, comments state why. ImageJ macro has no local
// scoping or typed arrays, so some globals and the newArray idiom are unavoidable.
// ASCII only.

requires("1.53f");

if (nImages == 0)
    exit("Open an image first.");

if (selectionType() < 0)
    waitForUser("Draw ROI", "Draw a rectangle around the first lane, then click OK.");

if (selectionType() < 0)
    exit("No selection. Draw a rectangle and run again.");

getSelectionBounds(source_x, source_y, source_width, source_height);
source_roi_name = Roi.getName;
if (source_roi_name == "")
    source_roi_name = "source";

getPixelSize(unit, pixel_width, pixel_height);
// A gap given in physical units only makes sense on a calibrated image; otherwise the
// gap is in pixels. "pixel"/"pixels" units mean uncalibrated.
image_is_calibrated = (pixel_width != 1 || pixel_height != 1) && unit != "pixel" && unit != "pixels";

number_of_copies = 5;
gap_x = 10;
gap_y = 10;
layout = "Row (X)";
gaps_in_physical_units = image_is_calibrated;
paste_pixel_content = true;
add_copies_to_roi_manager = true;
keep_original_selection = true;

layouts = newArray("Row (X)", "Column (Y)");

Dialog.create("Distribute copies");
Dialog.addMessage("Source: " + source_width + " x " + source_height + " px at (" + source_x + ", " + source_y + ")");
Dialog.addNumber("Number of copies (not including source)", number_of_copies, 0, 6, "");
Dialog.addChoice("Layout", layouts, layout);
Dialog.addNumber("Gap X (between bounding boxes)", gap_x, 2, 8, "");
Dialog.addNumber("Gap Y (between bounding boxes)", gap_y, 2, 8, "");
if (image_is_calibrated)
    Dialog.addCheckbox("Gaps are in " + unit + " (unchecked = pixels)", gaps_in_physical_units);
Dialog.addCheckbox("Paste pixel content", paste_pixel_content);
Dialog.addCheckbox("Add copies to ROI Manager", add_copies_to_roi_manager);
Dialog.addCheckbox("Keep original selection", keep_original_selection);
Dialog.addMessage("A preview (yellow source, cyan copies) is shown before applying.");
Dialog.show();

number_of_copies = Dialog.getNumber();
layout = Dialog.getChoice();
gap_x = Dialog.getNumber();
gap_y = Dialog.getNumber();
if (image_is_calibrated)
    gaps_in_physical_units = Dialog.getCheckbox();
paste_pixel_content = Dialog.getCheckbox();
add_copies_to_roi_manager = Dialog.getCheckbox();
keep_original_selection = Dialog.getCheckbox();

if (number_of_copies < 1)
    exit("Need at least 1 copy.");

// Convert a physical-unit gap to pixels for placement; placement math is all in pixels.
gap_x_pixels = gap_x;
gap_y_pixels = gap_y;
if (image_is_calibrated && gaps_in_physical_units) {
    gap_x_pixels = gap_x / pixel_width;
    gap_y_pixels = gap_y / pixel_height;
}

// Pitch between copies is one bounding box plus the gap.
step_x = source_width + gap_x_pixels;
step_y = source_height + gap_y_pixels;

copy_coordinates = planPositions(source_x, source_y, number_of_copies, layout, step_x, step_y);

// Preview: draw the source and the planned copies as an overlay, then a single
// blocking confirm. (This is not a live-updating preview; to change the numbers,
// cancel the confirm and re-run the macro.)
Overlay.remove;
drawPreview(source_x, source_y, source_width, source_height, copy_coordinates);
Overlay.show;

getDimensions(image_width, image_height, channel_count, slice_count, frame_count);
off_canvas_count = countOffCanvas(copy_coordinates, source_width, source_height, image_width, image_height);

preview_message = "Preview: " + number_of_copies + " copy(ies), step "
    + d2s(step_x, 1) + " x " + d2s(step_y, 1) + " px.";
if (off_canvas_count > 0)
    preview_message = preview_message + "\nWarning: " + off_canvas_count + " copy(ies) extend past the image.";
preview_message = preview_message + "\n\nOK applies this layout. To change values, close this and re-run the macro.";
showMessage("Preview", preview_message);

// Copy the source pixels once, before placement, so every Paste uses the same content.
makeRectangle(source_x, source_y, source_width, source_height);
if (selectionType() >= 0 && paste_pixel_content)
    run("Copy");

if (add_copies_to_roi_manager) {
    roiManager("reset");
    roiManager("Add");
    roiManager("Select", 0);
    roiManager("Rename", source_roi_name + "_0");
}

// Batch mode hides the per-paste redraws; much faster for many copies.
setBatchMode("hide");
placed_count = 0;
for (copy_index = 0; copy_index < copy_coordinates.length / 2; copy_index++) {
    copy_x = copy_coordinates[copy_index * 2];
    copy_y = copy_coordinates[copy_index * 2 + 1];
    makeRectangle(round(copy_x), round(copy_y), source_width, source_height);
    if (paste_pixel_content)
        run("Paste");
    if (add_copies_to_roi_manager) {
        roiManager("Add");
        roiManager("Select", roiManager("count") - 1);
        roiManager("Rename", source_roi_name + "_" + (copy_index + 1));
    }
    placed_count++;
}
setBatchMode("exit and display");

if (keep_original_selection)
    makeRectangle(source_x, source_y, source_width, source_height);
else
    run("Select None");

showStatus("Placed " + placed_count + " copies.");

// Return a flat [x0,y0, x1,y1, ...] array of copy top-left corners. Row places copies
// to the right of the source; Column below it. The source itself is not in the array
// (it is added to the ROI Manager separately as _0).
function planPositions(origin_x, origin_y, count, mode, pitch_x, pitch_y) {
    coordinates = newArray(count * 2);
    if (mode == "Row (X)") {
        for (index = 0; index < count; index++) {
            coordinates[index * 2] = origin_x + (index + 1) * pitch_x;
            coordinates[index * 2 + 1] = origin_y;
        }
    } else {
        // Column (Y)
        for (index = 0; index < count; index++) {
            coordinates[index * 2] = origin_x;
            coordinates[index * 2 + 1] = origin_y + (index + 1) * pitch_y;
        }
    }
    return coordinates;
}

function drawPreview(origin_x, origin_y, width, height, coordinates) {
    Overlay.remove;
    makeRectangle(origin_x, origin_y, width, height);
    Overlay.addSelection("yellow", 2);
    for (index = 0; index < coordinates.length / 2; index++) {
        makeRectangle(round(coordinates[index * 2]), round(coordinates[index * 2 + 1]), width, height);
        Overlay.addSelection("cyan", 1);
    }
    Overlay.show;
}

function countOffCanvas(coordinates, width, height, image_width, image_height) {
    off_count = 0;
    for (index = 0; index < coordinates.length / 2; index++) {
        corner_x = coordinates[index * 2];
        corner_y = coordinates[index * 2 + 1];
        if (corner_x < 0 || corner_y < 0 || corner_x + width > image_width || corner_y + height > image_height)
            off_count++;
    }
    return off_count;
}
