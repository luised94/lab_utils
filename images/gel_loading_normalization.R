# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., jsonlite::fromJSON).
# Date created: 2026-08-10
# Purpose: consume the gel densitometry pipeline output (measure_gel.py) and turn
# per-band intensities into per-lane loading-normalized results plus diagnostic
# plots. This is a join -> normalize -> plot layer. It does NOT re-measure pixels
# (measure_gel.py already produced reported_value per (well, canonical band) with a
# recorded basis) and does NOT run inferential statistics (a single gel has no
# replicates; testing belongs in a later multi-gel script over these per-lane CSVs).
#
# Usage: pass ONE path, relative or absolute. Everything else is derived from it.
#   Rscript gel_loading_normalization.R <analysis_path> [sample_sheet_csv]
# where <analysis_path> is the gel analysis directory (named "<stem>_gel_analysis")
# or any file inside it. The script normalizes the path, resolves the analysis
# directory, and from it derives band_measurements.csv, band_detection_report.json,
# the sample sheet ("<stem>_sample_sheet.csv", searched in the analysis directory
# then its parent), and the output directory. A second argument overrides the
# sample-sheet path. It can also be source()d after setting COMMAND_LINE_PATH.
#
# Analysis choices are constants below (edit and re-run): CONTROL_BAND_INDEX,
# FLIP_DIRECTION, BASELINE_METHOD, REQUIRE_FONT.
#
# Output (into <analysis_dir>/loading_normalization):
#   loading_normalization_per_lane.csv   one row per well: both normalizations,
#                                         trust flags, provenance carried per row.
#   loading_normalization_per_band.csv   one row per (well, band).
#   plot_reconciliation_flipcheck.pdf    banding pattern per lane; flip-direction check.
#   plot_normalized_results.pdf          both normalizations, reference at 100%.
#   plot_baseline_basis_transparency.pdf where the two baselines diverged.

# ==============================================================================
# GIT STATE REFERENCE (manual-fill at deposit time; no runtime git calls)
# ==============================================================================
# Commit hash:    ____________________________________________
# Branch:         ____________________________________________
# Tag / release:  ____________________________________________
# Snapshot date:  ____________________________________________
# Repository URL: ____________________________________________
# ==============================================================================

# ==============================================================================
# Analysis-choice constants (edit these; they are the only per-run knobs)
# ==============================================================================
CONTROL_BAND_INDEX <- 1L   # loading-control canonical_band_index for band-normalization
FLIP_DIRECTION <- "direct" # 'direct' (well_index = well_number - 1) or 'flipped'
BASELINE_METHOD <- "valley_to_valley" # which baseline row to keep ('valley_to_valley'|'opening')
REQUIRE_FONT <- FALSE      # TRUE stops if no preferred font is installed

# ==============================================================================
# Fixed conventions (rarely changed)
# ==============================================================================
BAND_MEASUREMENTS_FILENAME <- "band_measurements.csv"
BAND_DETECTION_REPORT_FILENAME <- "band_detection_report.json"
ANALYSIS_DIRECTORY_SUFFIX <- "_gel_analysis"
SAMPLE_SHEET_SUFFIX <- "_sample_sheet.csv"
OUTPUT_SUBDIRECTORY_NAME <- "loading_normalization"
EXPECTED_WELL_COUNT <- 15L
ROLE_VOCABULARY <- c("reference", "sample", "empty")
ROLES_EXCLUDED_FROM_NORMALIZATION <- c("empty")

stopifnot(
    "FLIP_DIRECTION must be 'direct' or 'flipped'." =
        FLIP_DIRECTION %in% c("direct", "flipped"),
    "BASELINE_METHOD must be 'valley_to_valley' or 'opening'." =
        BASELINE_METHOD %in% c("valley_to_valley", "opening")
)

# ==============================================================================
# Resolve the one input path, then derive everything from it
# ==============================================================================
# Accept the path from the command line, or from COMMAND_LINE_PATH when source()d.
command_line_arguments <- commandArgs(trailingOnly = TRUE)
if (length(command_line_arguments) >= 1) {
    input_path_argument <- command_line_arguments[1]
} else if (exists("COMMAND_LINE_PATH")) {
    input_path_argument <- COMMAND_LINE_PATH
} else {
    stop(
        "No path provided. Usage:\n",
        "  Rscript gel_loading_normalization.R <analysis_path> [sample_sheet_csv]"
    )
}
sample_sheet_override <- if (length(command_line_arguments) >= 2) {
    command_line_arguments[2]
} else if (exists("SAMPLE_SHEET_OVERRIDE_PATH")) {
    SAMPLE_SHEET_OVERRIDE_PATH
} else {
    NA_character_
}

# normalizePath collapses "..", resolves symlinks, and makes the path absolute so
# every derived path is unambiguous regardless of the working directory.
normalized_input_path <- normalizePath(input_path_argument, mustWork = TRUE)
analysis_directory <- if (dir.exists(normalized_input_path)) {
    normalized_input_path
} else {
    dirname(normalized_input_path)
}
gel_stem <- sub(paste0(ANALYSIS_DIRECTORY_SUFFIX, "$"), "", basename(analysis_directory))

band_measurements_filepath <- file.path(analysis_directory, BAND_MEASUREMENTS_FILENAME)
band_detection_report_filepath <- file.path(analysis_directory, BAND_DETECTION_REPORT_FILENAME)

# Sample sheet: explicit override, else derived name in the analysis directory,
# else the same name in the parent directory.
derive_sample_sheet_path <- function() {
    if (!is.na(sample_sheet_override)) {
        return(normalizePath(sample_sheet_override, mustWork = FALSE))
    }
    derived_name <- paste0(gel_stem, SAMPLE_SHEET_SUFFIX)
    candidate_paths <- c(
        file.path(analysis_directory, derived_name),
        file.path(dirname(analysis_directory), derived_name)
    )
    existing <- candidate_paths[file.exists(candidate_paths)]
    if (length(existing) >= 1) existing[1] else candidate_paths[1]
}
sample_sheet_filepath <- derive_sample_sheet_path()
output_directory <- file.path(analysis_directory, OUTPUT_SUBDIRECTORY_NAME)

# Fail fast if any required input is missing, naming every one that is absent.
required_inputs <- c(
    "band measurements" = band_measurements_filepath,
    "gel detection report" = band_detection_report_filepath,
    "sample sheet" = sample_sheet_filepath
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
    stop(
        "Missing required input(s) derived from path '", input_path_argument, "':\n",
        paste0("  ", names(missing_inputs), ": ", missing_inputs, collapse = "\n")
    )
}

# ==============================================================================
# Configuration
# ==============================================================================
library(tidyverse)

stopifnot(
    "Cairo graphics not available. Install libcairo2-dev and rebuild R." =
        capabilities("cairo")
)

installed_font_families <- function() {
    raw_lines <- tryCatch(
        system2("fc-list", args = c(":", "family"), stdout = TRUE, stderr = FALSE),
        error = function(e) character(0)
    )
    unique(trimws(unlist(strsplit(raw_lines, ","))))
}
PREFERRED_FONT_FAMILIES <- c("Arial", "Liberation Sans", "DejaVu Sans")
present_font_families <- installed_font_families()
resolved_font_family <- PREFERRED_FONT_FAMILIES[
    PREFERRED_FONT_FAMILIES %in% present_font_families
][1]
if (is.na(resolved_font_family)) {
    if (REQUIRE_FONT) {
        stop("No preferred font installed and REQUIRE_FONT=TRUE: ",
             paste(PREFERRED_FONT_FAMILIES, collapse = ", "))
    }
    resolved_font_family <- "sans"
    warning("No preferred font found; falling back to the generic 'sans' family.")
}
message("Using font family: ", resolved_font_family)

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

PLOT_CONFIG <- list(
    detected_fill = "#377EB8",
    predicted_fill = "#BDBDBD",
    reference_outline = "#E41A1C",
    signal_gradient = c("#FFFFFF", "#08519C"),
    divergence_gradient = c("#FFFFFF", "#D73027"),
    theme = list(base_family = resolved_font_family, base_size = 11),
    output = list(device = cairo_pdf, width = 7.2, height = 5.6)
)
message("Analysis directory: ", analysis_directory)
message("Gel stem: ", gel_stem)

if (!dir.exists(output_directory)) {
    dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# Gel-level provenance (band_detection_report.json, schema >= 2)
# ==============================================================================
gel_report <- jsonlite::fromJSON(band_detection_report_filepath, simplifyVector = TRUE)
band_detection_report_schema_version <- gel_report[["band_detection_report_schema_version"]]
if (is.null(band_detection_report_schema_version) ||
    band_detection_report_schema_version < 2) {
    warning(
        "band_detection_report schema version is ",
        if (is.null(band_detection_report_schema_version)) "absent"
        else band_detection_report_schema_version,
        "; expected >= 2 for the self-sufficient provenance block."
    )
}
gel_overall_status <- gel_report[["overall_status"]]
gel_warning_check_names <- gel_report[["warning_check_names"]]
gel_encoding_verified <- gel_report[["encoding_verified"]]
gel_provenance_block <- gel_report[["provenance"]]
gel_input_sha256 <- if (!is.null(gel_provenance_block)) gel_provenance_block[["input_file_sha256"]] else NA_character_
gel_exposure_time_text <- if (!is.null(gel_provenance_block)) gel_provenance_block[["exposure_time_text"]] else NA_character_

if (isFALSE(gel_encoding_verified)) {
    message("NOTE: encoding_verified is FALSE. Every normalized number is ",
            "PROVISIONAL until a dilution series confirms linearity.")
}
message("Gel status: ", gel_overall_status,
        " | warnings: ", paste(gel_warning_check_names, collapse = ", "))

# ==============================================================================
# Per-band measurements: read, validate, de-duplicate across baselines
# ==============================================================================
raw_band_measurements <- readr::read_csv(band_measurements_filepath, show_col_types = FALSE)

REQUIRED_MEASUREMENT_COLUMNS <- c(
    "well_index", "lane_detection_status", "canonical_band_index",
    "canonical_position_migration_millimetres", "baseline_method", "occupancy",
    "measurement_status", "saturation_status", "baseline_disagreement_fraction",
    "baseline_agreement_status", "reported_value", "reported_value_basis"
)
stopifnot(
    "band_measurements.csv is missing required columns." =
        all(REQUIRED_MEASUREMENT_COLUMNS %in% colnames(raw_band_measurements))
)

# reported_value is one chosen value per (well, band), duplicated across the two
# baseline rows; the baseline pair carries only divergence diagnostics. Assert that
# invariant so a future pipeline change that made them differ cannot slip through
# and be silently summed twice.
reported_value_disagreements <- raw_band_measurements %>%
    group_by(.data$well_index, .data$canonical_band_index) %>%
    summarise(distinct_reported_values = n_distinct(.data$reported_value), .groups = "drop") %>%
    filter(.data$distinct_reported_values > 1)
if (nrow(reported_value_disagreements) > 0) {
    stop(
        "reported_value differs across baseline rows for ",
        nrow(reported_value_disagreements),
        " (well, band) cells; the de-duplication assumption is violated."
    )
}

per_band <- raw_band_measurements %>%
    filter(.data$baseline_method == BASELINE_METHOD)
stopifnot(
    "Chosen BASELINE_METHOD produced no rows." = nrow(per_band) > 0,
    "De-duplication did not yield one row per (well, band)." =
        nrow(per_band) == raw_band_measurements %>%
            distinct(.data$well_index, .data$canonical_band_index) %>% nrow()
)
message("Per-band rows after de-duplication: ", nrow(per_band),
        " (", n_distinct(per_band$well_index), " wells x ",
        n_distinct(per_band$canonical_band_index), " bands)")

# ==============================================================================
# Sample sheet: read (CSV), defend against whitespace, validate
# ==============================================================================
raw_sample_sheet <- readr::read_csv(sample_sheet_filepath, show_col_types = FALSE)

# Trim leading/trailing whitespace on every character column (operator sheets carry
# cells like "BK Box 1 clear " whose trailing space would split factor levels).
sample_sheet <- raw_sample_sheet %>%
    mutate(across(where(is.character), ~ str_trim(.x)))

REQUIRED_SHEET_COLUMNS <- c("well_number", "sample_label", "role")
stopifnot(
    "sample sheet is missing required columns (well_number, sample_label, role)." =
        all(REQUIRED_SHEET_COLUMNS %in% colnames(sample_sheet)),
    "sample sheet must have exactly one row per well." =
        nrow(sample_sheet) == EXPECTED_WELL_COUNT,
    "well_number must be the integers 1..15, each once." =
        setequal(sample_sheet$well_number, seq_len(EXPECTED_WELL_COUNT)),
    "role values must be within the controlled vocabulary." =
        all(sample_sheet$role %in% ROLE_VOCABULARY)
)
has_normalize_on_band_column <- "normalize_on_band" %in% colnames(sample_sheet)

reference_rows <- sample_sheet %>% filter(.data$role == "reference")
if (nrow(reference_rows) != 1) {
    stop("Normalization needs exactly one row with role='reference'; found ",
         nrow(reference_rows), ".")
}
message("Sample sheet validated: ", nrow(sample_sheet), " wells, reference = well ",
        reference_rows$well_number, " (", reference_rows$sample_label, ").")

# ==============================================================================
# Join key: map physical well_number to pipeline well_index via the flip flag
# ==============================================================================
# well_index is the pipeline's flip-corrected 0..14 comb index; it cannot know the
# operator's physical numbering direction. The flip flag chooses the mapping; the
# reconciliation plot is how the operator confirms it.
sample_sheet <- sample_sheet %>%
    mutate(
        well_index = if (FLIP_DIRECTION == "direct") {
            .data$well_number - 1L
        } else {
            EXPECTED_WELL_COUNT - .data$well_number
        }
    )

# ==============================================================================
# Join and validate: every well maps; loaded-but-undetected is a loading failure
# ==============================================================================
lane_detection_by_well <- per_band %>%
    distinct(.data$well_index, .data$lane_detection_status)

joined_identity <- sample_sheet %>%
    left_join(lane_detection_by_well, by = "well_index")

unmatched_sheet_wells <- joined_identity %>% filter(is.na(.data$lane_detection_status))
if (nrow(unmatched_sheet_wells) > 0) {
    stop(
        "These sample-sheet wells have no matching well_index in the measurements ",
        "(check FLIP_DIRECTION): well_number ",
        paste(unmatched_sheet_wells$well_number, collapse = ", ")
    )
}
measurement_wells_missing_from_sheet <- setdiff(
    unique(per_band$well_index), sample_sheet$well_index
)
if (length(measurement_wells_missing_from_sheet) > 0) {
    warning("Measured well_index values absent from the sample sheet: ",
            paste(measurement_wells_missing_from_sheet, collapse = ", "))
}

joined_identity <- joined_identity %>%
    mutate(
        is_loading_failure = .data$role == "sample" &
            .data$lane_detection_status == "predicted_from_comb"
    )
loading_failures <- joined_identity %>% filter(.data$is_loading_failure)
if (nrow(loading_failures) > 0) {
    warning("Loading failure(s): well_number ",
            paste(loading_failures$well_number, collapse = ", "),
            " marked as loaded samples but the lane was predicted, not detected.")
}

# ==============================================================================
# Per-band table joined to identity (for output and the transparency plot)
# ==============================================================================
per_band_with_identity <- per_band %>%
    left_join(
        joined_identity %>% select(
            "well_index", "well_number", "sample_label", "role", "is_loading_failure"
        ),
        by = "well_index"
    ) %>%
    mutate(
        gel_overall_status = gel_overall_status,
        encoding_verified = gel_encoding_verified,
        input_file_sha256 = gel_input_sha256
    )

# ==============================================================================
# Normalization 1: whole-lane total, relative to the reference lane
# ==============================================================================
whole_lane_totals <- per_band %>%
    group_by(.data$well_index) %>%
    summarise(whole_lane_total = sum(.data$reported_value), .groups = "drop")

reference_well_index <- sample_sheet %>%
    filter(.data$role == "reference") %>%
    pull(.data$well_index)
reference_whole_lane_total <- whole_lane_totals %>%
    filter(.data$well_index == reference_well_index) %>%
    pull(.data$whole_lane_total)
stopifnot(
    "Reference lane whole-lane total is not positive; cannot normalize." =
        length(reference_whole_lane_total) == 1 && reference_whole_lane_total > 0
)

# ==============================================================================
# Normalization 2: designated loading-control band, relative to the reference
# ==============================================================================
control_band_by_well <- joined_identity %>%
    mutate(
        control_band_index = if (has_normalize_on_band_column) {
            coalesce(as.integer(.data$normalize_on_band), CONTROL_BAND_INDEX)
        } else {
            CONTROL_BAND_INDEX
        }
    ) %>%
    select("well_index", "control_band_index")

control_band_values <- per_band %>%
    inner_join(control_band_by_well, by = "well_index") %>%
    filter(.data$canonical_band_index == .data$control_band_index) %>%
    select("well_index", "control_band_index", control_band_value = "reported_value")

reference_control_band_value <- control_band_values %>%
    filter(.data$well_index == reference_well_index) %>%
    pull(.data$control_band_value)
if (length(reference_control_band_value) != 1 || reference_control_band_value <= 0) {
    warning("Reference lane has no positive signal at control band index ",
            CONTROL_BAND_INDEX, "; band-normalization will be NA.")
    reference_control_band_value <- NA_real_
}

# ==============================================================================
# Assemble per-lane results
# ==============================================================================
per_lane_results <- joined_identity %>%
    left_join(whole_lane_totals, by = "well_index") %>%
    left_join(control_band_values, by = "well_index") %>%
    mutate(
        percent_of_reference_whole_lane =
            100 * .data$whole_lane_total / reference_whole_lane_total,
        percent_of_reference_control_band =
            100 * .data$control_band_value / reference_control_band_value,
        is_excluded_from_normalization = .data$role %in% ROLES_EXCLUDED_FROM_NORMALIZATION,
        gel_overall_status = gel_overall_status,
        encoding_verified = gel_encoding_verified,
        input_file_sha256 = gel_input_sha256,
        exposure_time_text = gel_exposure_time_text,
        flip_direction = FLIP_DIRECTION
    ) %>%
    arrange(.data$well_number)

# ==============================================================================
# Write CSV outputs
# ==============================================================================
per_lane_csv_filepath <- file.path(output_directory, "loading_normalization_per_lane.csv")
per_band_csv_filepath <- file.path(output_directory, "loading_normalization_per_band.csv")
if (!file.exists(per_lane_csv_filepath) || OVERWRITE_CSVS) {
    write.csv(per_lane_results, per_lane_csv_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(per_lane_csv_filepath))
}
if (!file.exists(per_band_csv_filepath) || OVERWRITE_CSVS) {
    write.csv(per_band_with_identity, per_band_csv_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(per_band_csv_filepath))
}

# ==============================================================================
# Shared plot helpers
# ==============================================================================
# Captions are wrapped and the provisional-encoding note moved to a subtitle so
# long text is not clipped at the device's bottom edge.
wrapped_caption <- function(caption_text) {
    str_wrap(caption_text, width = 110)
}
provisional_subtitle <- function() {
    encoding_note <- if (isTRUE(gel_encoding_verified)) {
        "encoding_verified = TRUE"
    } else {
        "PROVISIONAL: encoding_verified = FALSE (values not confirmed linear)"
    }
    paste0(encoding_note, "  |  flip = ", FLIP_DIRECTION,
           "  |  gel status = ", gel_overall_status)
}
base_theme <- function() {
    theme_classic(
        base_size = PLOT_CONFIG$theme$base_size,
        base_family = PLOT_CONFIG$theme$base_family
    ) +
        theme(
            plot.subtitle = element_text(size = 8),
            plot.caption = element_text(size = 7, hjust = 0),
            plot.caption.position = "plot",
            plot.margin = margin(t = 6, r = 10, b = 12, l = 6)
        )
}
save_plot <- function(plot_object, filepath) {
    if (!file.exists(filepath) || OVERWRITE_PLOTS) {
        ggsave(filepath, plot_object,
               device = PLOT_CONFIG$output$device,
               width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
        message("Saved plot: ", basename(filepath))
    }
}

lane_axis_labels <- per_lane_results %>%
    transmute(
        .data$well_number,
        axis_label = paste0(
            "W", .data$well_number, "  ", .data$sample_label,
            ifelse(.data$lane_detection_status == "detected", "", "  (pred)")
        )
    ) %>%
    arrange(.data$well_number)

# ==============================================================================
# Plot 1: reconciliation / flip check
# ==============================================================================
reconciliation_data <- per_band_with_identity %>%
    left_join(lane_axis_labels, by = "well_number") %>%
    mutate(axis_label = factor(.data$axis_label, levels = rev(lane_axis_labels$axis_label)))
reconciliation_plot <- ggplot(
    reconciliation_data,
    aes(x = factor(.data$canonical_band_index), y = .data$axis_label,
        fill = .data$reported_value)
) +
    geom_tile(color = "grey85", linewidth = 0.3) +
    geom_point(
        data = reconciliation_data %>% filter(.data$occupancy == "locally_detected"),
        aes(x = factor(.data$canonical_band_index), y = .data$axis_label),
        inherit.aes = FALSE, shape = 3, size = 1.1, color = "grey20"
    ) +
    scale_fill_gradientn(colors = PLOT_CONFIG$signal_gradient, name = "reported_value",
                         trans = "sqrt") +
    labs(
        x = "canonical_band_index (migration order)",
        y = "lane: physical well + operator label",
        title = "Reconciliation / flip check: banding pattern per lane",
        subtitle = provisional_subtitle(),
        caption = wrapped_caption(paste0(
            "Crosses mark locally-detected bands; blank cells were measured as ~0. ",
            "Confirm the detected lanes sit on the wells you loaded; flip FLIP_DIRECTION if not."
        ))
    ) +
    base_theme() +
    theme(axis.text.y = element_text(size = 7))
save_plot(reconciliation_plot, file.path(output_directory, "plot_reconciliation_flipcheck.pdf"))

# ==============================================================================
# Plot 2: normalized results, both methods, reference at 100%
# ==============================================================================
results_long <- per_lane_results %>%
    select("well_number", "sample_label", "role", "lane_detection_status",
           "is_excluded_from_normalization",
           whole_lane = "percent_of_reference_whole_lane",
           control_band = "percent_of_reference_control_band") %>%
    pivot_longer(c("whole_lane", "control_band"),
                 names_to = "normalization_method", values_to = "percent_of_reference") %>%
    mutate(
        normalization_method = recode(
            .data$normalization_method,
            whole_lane = "Whole-lane total",
            control_band = paste0("Control band (index ", CONTROL_BAND_INDEX, ")")
        ),
        lane_label = factor(
            paste0("W", .data$well_number, " ", .data$sample_label),
            levels = paste0("W", per_lane_results$well_number, " ", per_lane_results$sample_label)
        ),
        lane_class = ifelse(
            .data$role == "reference", "reference",
            ifelse(.data$role == "empty", "empty",
                   ifelse(.data$lane_detection_status == "detected",
                          "loaded (detected)", "loaded (predicted)"))
        )
    )
results_plot <- ggplot(
    results_long,
    aes(x = .data$lane_label, y = .data$percent_of_reference, fill = .data$lane_class)
) +
    geom_col(color = "black", linewidth = 0.3, width = 0.75) +
    geom_hline(yintercept = 100, linetype = "dashed", color = PLOT_CONFIG$reference_outline) +
    facet_wrap(~ .data$normalization_method, ncol = 1, scales = "free_y") +
    scale_fill_manual(
        values = c("loaded (detected)" = PLOT_CONFIG$detected_fill,
                   "loaded (predicted)" = PLOT_CONFIG$predicted_fill,
                   "reference" = PLOT_CONFIG$reference_outline,
                   "empty" = "#E0E0E0"),
        name = "lane"
    ) +
    labs(
        x = "lane (physical well + label)",
        y = "percent of reference (%)",
        title = "Loading normalization (reference = 100%)",
        subtitle = provisional_subtitle(),
        caption = wrapped_caption(
            "Empty lanes shown as measured noise, not dropped. Dashed line = reference (100%)."
        )
    ) +
    base_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
save_plot(results_plot, file.path(output_directory, "plot_normalized_results.pdf"))

# ==============================================================================
# Plot 3: baseline / basis transparency
# ==============================================================================
transparency_data <- per_band_with_identity %>%
    left_join(lane_axis_labels, by = "well_number") %>%
    mutate(axis_label = factor(.data$axis_label, levels = rev(lane_axis_labels$axis_label)))
transparency_plot <- ggplot(
    transparency_data,
    aes(x = factor(.data$canonical_band_index), y = .data$axis_label,
        fill = .data$baseline_disagreement_fraction)
) +
    geom_tile(color = "grey85", linewidth = 0.3) +
    geom_point(
        data = transparency_data %>%
            filter(.data$reported_value_basis == "opening_net_cluster_inclusive"),
        aes(x = factor(.data$canonical_band_index), y = .data$axis_label),
        inherit.aes = FALSE, shape = 4, size = 1.3, color = "grey10"
    ) +
    scale_fill_gradientn(colors = PLOT_CONFIG$divergence_gradient, name = "baseline\ndisagreement") +
    labs(
        x = "canonical_band_index",
        y = "lane: physical well + operator label",
        title = "Baseline divergence and reported-value basis",
        subtitle = provisional_subtitle(),
        caption = wrapped_caption(paste0(
            "Fill = fraction disagreement between the two baselines. ",
            "X marks cells whose reported value used the opening-cluster basis rather than region-net."
        ))
    ) +
    base_theme() +
    theme(axis.text.y = element_text(size = 7))
save_plot(transparency_plot, file.path(output_directory, "plot_baseline_basis_transparency.pdf"))

# ==============================================================================
# Console summary
# ==============================================================================
message("\nPer-lane normalized results:")
print(as.data.frame(
    per_lane_results %>% select(
        "well_number", "sample_label", "role", "lane_detection_status",
        "whole_lane_total", "percent_of_reference_whole_lane",
        "percent_of_reference_control_band", "is_loading_failure"
    )
))
message("\nScript complete. Outputs in: ", output_directory)
