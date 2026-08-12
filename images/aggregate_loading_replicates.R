# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced.
# Date created: 2026-08-11
# Purpose: combine several single-gel loading results (each produced by
# gel_loading_normalization.R as loading_normalization_per_lane.csv) into
# replicate-level DESCRIPTIVE results: mean, standard deviation, standard error,
# and all-against-all percent-difference and fold-difference matrices. It does
# NOT run inferential statistics (see loading_statistics_reference.R for that
# machinery, kept separate on purpose). It writes every intermediate state to
# disk: the per-gel rows read together, the processed analysis frame, and the
# summary frame handed to the plot. Pre-normalized totals are preserved
# alongside the percent-of-reference values.
#
# The repeat relationship lives ONLY in the manifest (one row per gel). The
# per-gel sample sheet carries within-gel identity (conditions, roles); the
# manifest carries experiment_id and replicate. This script is single-purpose:
# it does not re-measure pixels and does not re-normalize (each gel's reference
# lane is already pinned to 100 by gel_loading_normalization.R).
#
# Usage:
#   Rscript aggregate_loading_replicates.R <manifest_csv>
# The manifest is the only required argument. Paths inside it may be relative to
# the manifest's own directory, and analysis_path may be an analysis directory
# ("<stem>_gel_analysis") or any file inside it. Outputs are written under
# <manifest_directory>/aggregate_analysis. It can also be source()d after
# setting COMMAND_LINE_MANIFEST_PATH.
#
# Manifest columns (see manifest_template.csv):
#   experiment_id  required   groups replicates into one experiment
#   analysis_path  required   the gel analysis dir, or any file inside it
#   gel_id         required   verified against the per-lane CSV; hard stop on mismatch
#   replicate      required   biological replicate label within the experiment
#   notes          optional   free text, carried through untouched

# ==============================================================================
# Analysis-choice constants (edit these; they are the only per-run knobs)
# ==============================================================================
# CONDITION_KEY_COLUMNS names the sample-sheet factor column(s) that define "the
# same condition across replicates". It is intentionally per-run: set it to the
# one column your experiment varies (e.g. "suppressor" or "salt_condition_mM").
# If these columns exist but are entirely blank (a fast verification run where
# the sheet was not filled), the script WARNS and still writes the read-together
# and processed intermediates plus a per-well fallback summary, but skips the
# condition-level products (matrices, condition summary, plot). One key set
# applies to every experiment in a single manifest run; for experiments needing
# different keys, run the manifest once per experiment.
CONDITION_KEY_COLUMNS <- c("suppressor")

# Which per-lane normalization is the reported value. "whole_lane" uses
# percent_of_reference_whole_lane; "control_band" uses
# percent_of_reference_control_band. The matching pre-normalized total is
# preserved in the outputs either way.
NORMALIZATION_COLUMN <- "whole_lane"

# ==============================================================================
# Fixed conventions (rarely changed)
# ==============================================================================
PER_LANE_SUBDIRECTORY_NAME <- "loading_normalization"
PER_LANE_FILENAME <- "loading_normalization_per_lane.csv"
ANALYSIS_DIRECTORY_SUFFIX <- "_gel_analysis"
OUTPUT_DIRECTORY_NAME <- "aggregate_analysis"
REFERENCE_ROLE <- "reference"
ROLES_INCLUDED_IN_ANALYSIS <- c("reference", "sample")  # role == "empty" is dropped
REFERENCE_PERCENT_TARGET <- 100
REFERENCE_PERCENT_TOLERANCE <- 1e-4
REQUIRE_FONT <- FALSE

stopifnot(
    "NORMALIZATION_COLUMN must be 'whole_lane' or 'control_band'." =
        NORMALIZATION_COLUMN %in% c("whole_lane", "control_band"),
    "CONDITION_KEY_COLUMNS must name at least one column." =
        length(CONDITION_KEY_COLUMNS) >= 1
)
reported_value_source_column <- if (NORMALIZATION_COLUMN == "whole_lane") {
    "percent_of_reference_whole_lane"
} else {
    "percent_of_reference_control_band"
}
prenormalized_source_column <- if (NORMALIZATION_COLUMN == "whole_lane") {
    "whole_lane_total"
} else {
    "control_band_value"
}

# ==============================================================================
# Resolve the manifest path, then derive everything from it
# ==============================================================================
command_line_arguments <- commandArgs(trailingOnly = TRUE)
if (length(command_line_arguments) >= 1) {
    manifest_path_argument <- command_line_arguments[1]
} else if (exists("COMMAND_LINE_MANIFEST_PATH")) {
    manifest_path_argument <- COMMAND_LINE_MANIFEST_PATH
} else {
    stop(
        "No manifest provided. Usage:\n",
        "  Rscript aggregate_loading_replicates.R <manifest_csv>"
    )
}
manifest_absolute_path <- normalizePath(manifest_path_argument, mustWork = TRUE)
manifest_directory <- dirname(manifest_absolute_path)
output_directory <- file.path(manifest_directory, OUTPUT_DIRECTORY_NAME)
if (!dir.exists(output_directory)) {
    dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
}

library(tidyverse)

# ==============================================================================
# Read and validate the manifest
# ==============================================================================
raw_manifest <- readr::read_csv(manifest_absolute_path, show_col_types = FALSE)
# Trim whitespace on character columns so a trailing space cannot split an
# experiment_id or a replicate label into two.
manifest <- raw_manifest %>% mutate(across(where(is.character), ~ str_trim(.x)))

REQUIRED_MANIFEST_COLUMNS <- c("experiment_id", "analysis_path", "gel_id", "replicate")
missing_manifest_columns <- setdiff(REQUIRED_MANIFEST_COLUMNS, colnames(manifest))
if (length(missing_manifest_columns) > 0) {
    stop("Manifest is missing required column(s): ",
         paste(missing_manifest_columns, collapse = ", "),
         ". Required: ", paste(REQUIRED_MANIFEST_COLUMNS, collapse = ", "), ".")
}
manifest <- manifest %>% mutate(
    experiment_id = as.character(.data$experiment_id),
    gel_id = as.character(.data$gel_id),
    replicate = as.character(.data$replicate)
)
# Row-level completeness: name the offending manifest row rather than failing
# with a vague message, so the fix is obvious.
for (manifest_row_index in seq_len(nrow(manifest))) {
    for (required_column in REQUIRED_MANIFEST_COLUMNS) {
        cell_value <- manifest[[required_column]][manifest_row_index]
        if (is.na(cell_value) || cell_value == "") {
            stop("Manifest row ", manifest_row_index, " has an empty required field '",
                 required_column, "'. Fill it or remove the row.")
        }
    }
}
# replicate must be unique within an experiment (two rows cannot both be
# 'replicate 1' of the same experiment).
duplicate_replicates <- manifest %>%
    count(.data$experiment_id, .data$replicate) %>%
    filter(.data$n > 1)
if (nrow(duplicate_replicates) > 0) {
    stop("Manifest has repeated replicate labels within an experiment:\n",
         paste(capture.output(print(as.data.frame(duplicate_replicates))), collapse = "\n"))
}
message("Manifest validated: ", nrow(manifest), " gel(s), ",
        n_distinct(manifest$experiment_id), " experiment(s).")

# ==============================================================================
# Resolve each per-lane CSV, read it, verify gel_id, stack
# ==============================================================================
stacked_row_frames <- list()
for (manifest_row_index in seq_len(nrow(manifest))) {
    analysis_path_as_written <- manifest$analysis_path[manifest_row_index]
    # Relative paths resolve against the manifest directory so the manifest can
    # sit next to the data and use short paths.
    analysis_path_candidate <- if (
        grepl("^(/|~|[A-Za-z]:)", analysis_path_as_written)
    ) {
        analysis_path_as_written
    } else {
        file.path(manifest_directory, analysis_path_as_written)
    }
    if (!file.exists(analysis_path_candidate) && !dir.exists(analysis_path_candidate)) {
        stop("Manifest row ", manifest_row_index, ": analysis_path does not exist: ",
             analysis_path_candidate, " (from '", analysis_path_as_written, "').")
    }
    analysis_path_resolved <- normalizePath(analysis_path_candidate, mustWork = TRUE)
    analysis_directory <- if (dir.exists(analysis_path_resolved)) {
        analysis_path_resolved
    } else {
        dirname(analysis_path_resolved)
    }
    per_lane_filepath <- file.path(analysis_directory, PER_LANE_SUBDIRECTORY_NAME, PER_LANE_FILENAME)
    if (!file.exists(per_lane_filepath)) {
        stop("Manifest row ", manifest_row_index, ": expected per-lane file not found: ",
             per_lane_filepath, ". Run gel_loading_normalization.R on that gel first.")
    }

    one_gel_rows <- readr::read_csv(per_lane_filepath, show_col_types = FALSE)
    if (!("gel_id" %in% colnames(one_gel_rows))) {
        stop("Manifest row ", manifest_row_index, ": per-lane file has no gel_id column: ",
             per_lane_filepath, ". Re-run gel_loading_normalization.R (schema >= 3).")
    }
    distinct_gel_ids_in_file <- unique(as.character(one_gel_rows$gel_id))
    if (length(distinct_gel_ids_in_file) != 1) {
        stop("Manifest row ", manifest_row_index, ": per-lane file carries more than one gel_id: ",
             paste(distinct_gel_ids_in_file, collapse = ", "), ".")
    }
    # Keyed on both (safety/traceability): the manifest gel_id must match the file.
    if (!identical(distinct_gel_ids_in_file[1], manifest$gel_id[manifest_row_index])) {
        stop("Manifest row ", manifest_row_index, ": gel_id mismatch.\n",
             "  manifest gel_id:  ", manifest$gel_id[manifest_row_index], "\n",
             "  per-lane gel_id:  ", distinct_gel_ids_in_file[1], "\n",
             "  file: ", per_lane_filepath)
    }
    # A stray replicate column would mean an old sample sheet; the manifest is
    # authoritative, so warn and drop it.
    if ("replicate" %in% colnames(one_gel_rows)) {
        warning("Manifest row ", manifest_row_index, ": per-lane file carries a 'replicate' ",
                "column (old sample sheet). Ignoring it; the manifest is authoritative.")
        one_gel_rows <- one_gel_rows %>% select(-"replicate")
    }

    one_gel_rows <- one_gel_rows %>%
        mutate(
            experiment_id = manifest$experiment_id[manifest_row_index],
            replicate = manifest$replicate[manifest_row_index],
            manifest_notes = if ("notes" %in% colnames(manifest)) manifest$notes[manifest_row_index] else NA_character_
        )
    stacked_row_frames[[length(stacked_row_frames) + 1]] <- one_gel_rows
}
combined_raw_data <- bind_rows(stacked_row_frames)
# Put the identity columns first for readability.
combined_raw_data <- combined_raw_data %>%
    relocate("experiment_id", "replicate", "gel_id")

combined_raw_csv_path <- file.path(output_directory, "combined_raw_data.csv")
write.csv(combined_raw_data, combined_raw_csv_path, row.names = FALSE)
message("Wrote read-together state: ", basename(combined_raw_csv_path),
        " (", nrow(combined_raw_data), " rows).")

# ==============================================================================
# Build the processed analysis frame: drop empty lanes, select the reported
# value, verify the per-replicate reference anchor.
# ==============================================================================
required_value_columns <- c("role", "well_number", "sample_label",
                            reported_value_source_column, prenormalized_source_column)
missing_value_columns <- setdiff(required_value_columns, colnames(combined_raw_data))
if (length(missing_value_columns) > 0) {
    stop("Per-lane data is missing expected column(s): ",
         paste(missing_value_columns, collapse = ", "), ".")
}

processed_data <- combined_raw_data %>%
    filter(.data$role %in% ROLES_INCLUDED_IN_ANALYSIS) %>%
    mutate(
        reported_value = .data[[reported_value_source_column]],
        prenormalized_value = .data[[prenormalized_source_column]],
        reported_value_source = reported_value_source_column,
        normalization_column = NORMALIZATION_COLUMN
    )

# Exactly one reference per (experiment, replicate), pinned to ~100.
reference_counts <- processed_data %>%
    filter(.data$role == REFERENCE_ROLE) %>%
    count(.data$experiment_id, .data$replicate)
bad_reference_counts <- reference_counts %>% filter(.data$n != 1)
if (nrow(bad_reference_counts) > 0) {
    stop("Every (experiment, replicate) needs exactly one reference lane. Offending:\n",
         paste(capture.output(print(as.data.frame(bad_reference_counts))), collapse = "\n"))
}
reference_value_check <- processed_data %>%
    filter(.data$role == REFERENCE_ROLE) %>%
    filter(abs(.data$reported_value - REFERENCE_PERCENT_TARGET) > REFERENCE_PERCENT_TOLERANCE)
if (nrow(reference_value_check) > 0) {
    stop("Reference lane reported value is not ~100 in some (experiment, replicate); ",
         "the per-gel normalization anchor is violated:\n",
         paste(capture.output(print(as.data.frame(
             reference_value_check %>% select("experiment_id", "replicate", "reported_value")
         ))), collapse = "\n"))
}

# ==============================================================================
# Condition key: present, present-but-blank, or absent-column
# ==============================================================================
missing_key_columns <- setdiff(CONDITION_KEY_COLUMNS, colnames(processed_data))
if (length(missing_key_columns) > 0) {
    stop("CONDITION_KEY_COLUMNS names column(s) not present in the per-lane data: ",
         paste(missing_key_columns, collapse = ", "),
         ". Check the constant against the sample-sheet columns.")
}
# Coerce key columns to character and mark blanks (NA or empty after trimming).
for (key_column in CONDITION_KEY_COLUMNS) {
    processed_data[[key_column]] <- as.character(processed_data[[key_column]])
}
row_key_is_blank <- processed_data %>%
    transmute(across(all_of(CONDITION_KEY_COLUMNS),
                     ~ is.na(.x) | str_trim(.x) == "")) %>%
    reduce(`&`)
condition_key_is_absent <- all(row_key_is_blank)
condition_key_is_mixed <- (!condition_key_is_absent) && any(row_key_is_blank)

if (condition_key_is_mixed) {
    offending <- processed_data %>%
        mutate(row_number_in_processed = row_number()) %>%
        filter(row_key_is_blank) %>%
        select("row_number_in_processed", "experiment_id", "replicate", "well_number", "sample_label")
    stop("The condition key (", paste(CONDITION_KEY_COLUMNS, collapse = ", "),
         ") is set on some analysis rows but blank on others. Fill it on every ",
         "reference/sample lane, or clear it entirely for a keyless verification run. ",
         "Offending rows:\n",
         paste(capture.output(print(as.data.frame(offending))), collapse = "\n"))
}

# ==============================================================================
# Font resolution for plotting (graceful fallback; plots are edited in Illustrator)
# ==============================================================================
installed_families_raw <- tryCatch(
    system2("fc-list", args = c(":", "family"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
)
present_font_families <- unique(trimws(unlist(strsplit(installed_families_raw, ","))))
PREFERRED_FONT_FAMILIES <- c("Arial", "Liberation Sans", "DejaVu Sans")
resolved_font_family <- PREFERRED_FONT_FAMILIES[PREFERRED_FONT_FAMILIES %in% present_font_families][1]
if (is.na(resolved_font_family)) {
    if (REQUIRE_FONT) {
        stop("No preferred font installed and REQUIRE_FONT=TRUE: ",
             paste(PREFERRED_FONT_FAMILIES, collapse = ", "))
    }
    resolved_font_family <- "sans"
    warning("No preferred font found; falling back to the generic 'sans' family. ",
            "Text stays editable in the vector PDF regardless.")
}

PLOT_CONFIG <- list(
    reference_line_color = "#E41A1C",
    bar_fill = "#377EB8",
    bar = list(width = 0.7, color = "black", linewidth = 0.4),
    errorbar = list(width = 0.25, linewidth = 0.6),
    point = list(size = 2, fill = "grey30", color = "black", stroke = 0.4,
                 jitter_width = 0.12, jitter_seed = 42),
    output = list(device = cairo_pdf, width = 6.0, height = 4.5),
    base_size = 12
)

# ==============================================================================
# Keyless mode: write per-well fallback summary, skip condition-level products
# ==============================================================================
if (condition_key_is_absent) {
    warning("Condition key (", paste(CONDITION_KEY_COLUMNS, collapse = ", "),
            ") is blank on every analysis row. This is a keyless run: writing the ",
            "read-together and processed states and a per-well fallback summary, ",
            "and SKIPPING the matrices, condition summary, and plot.")
    processed_data$condition_label <- NA_character_
    processed_csv_path <- file.path(output_directory, "processed_data.csv")
    write.csv(processed_data, processed_csv_path, row.names = FALSE)
    message("Wrote processed state (keyless): ", basename(processed_csv_path),
            " (", nrow(processed_data), " analysis rows; empties dropped).")
    per_well_summary <- processed_data %>%
        group_by(.data$experiment_id, .data$well_number, .data$sample_label, .data$role) %>%
        summarise(
            replicate_count = n(),
            mean_reported_value = mean(.data$reported_value, na.rm = TRUE),
            sd_reported_value = sd(.data$reported_value, na.rm = TRUE),
            se_reported_value = sd(.data$reported_value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
        )
    for (experiment_identifier in unique(per_well_summary$experiment_id)) {
        experiment_output_directory <- file.path(
            output_directory, gsub("[^A-Za-z0-9._-]", "_", experiment_identifier)
        )
        if (!dir.exists(experiment_output_directory)) {
            dir.create(experiment_output_directory, showWarnings = FALSE, recursive = TRUE)
        }
        one_experiment_summary <- per_well_summary %>%
            filter(.data$experiment_id == experiment_identifier)
        summary_csv_path <- file.path(experiment_output_directory, "summary_data_per_well.csv")
        write.csv(one_experiment_summary, summary_csv_path, row.names = FALSE)
        message("Keyless: wrote per-well summary for '", experiment_identifier, "': ",
                basename(summary_csv_path))
    }
    message("Keyless mode: no condition-level products. Outputs in: ", output_directory)
} else {

# ==============================================================================
# Keyed mode: condition label, one row per (experiment, replicate, condition)
# ==============================================================================
processed_data <- processed_data %>%
    mutate(condition_label = do.call(paste, c(across(all_of(CONDITION_KEY_COLUMNS)), sep = " | ")))

processed_csv_path <- file.path(output_directory, "processed_data.csv")
write.csv(processed_data, processed_csv_path, row.names = FALSE)
message("Wrote processed state: ", basename(processed_csv_path),
        " (", nrow(processed_data), " analysis rows; empties dropped; condition_label added).")

# Assert one analysis row per (experiment, replicate, condition). Duplicates are
# most likely intended technical replicates loaded in the same gel; the aggregator
# does not silently average them.
condition_row_counts <- processed_data %>%
    count(.data$experiment_id, .data$replicate, .data$condition_label)
duplicate_condition_rows <- condition_row_counts %>% filter(.data$n > 1)
if (nrow(duplicate_condition_rows) > 0) {
    stop("More than one analysis lane maps to the same (experiment, replicate, condition):\n",
         paste(capture.output(print(as.data.frame(duplicate_condition_rows))), collapse = "\n"),
         "\nIf these are technical replicates, average them in the sample sheet, ",
         "or tell me to average within-gel duplicates.")
}

experiment_identifiers <- unique(processed_data$experiment_id)
combined_summary_frames <- list()

for (experiment_identifier in experiment_identifiers) {
    experiment_output_directory <- file.path(
        output_directory, gsub("[^A-Za-z0-9._-]", "_", experiment_identifier)
    )
    if (!dir.exists(experiment_output_directory)) {
        dir.create(experiment_output_directory, showWarnings = FALSE, recursive = TRUE)
    }
    one_experiment <- processed_data %>% filter(.data$experiment_id == experiment_identifier)

    experiment_conditions <- sort(unique(one_experiment$condition_label))
    experiment_replicates <- sort(unique(one_experiment$replicate))

    # ---- summary: mean / sd / se per condition -----------------------------
    condition_summary <- one_experiment %>%
        group_by(.data$condition_label) %>%
        summarise(
            role = paste(sort(unique(.data$role)), collapse = "/"),
            replicate_count = n(),
            mean_reported_value = mean(.data$reported_value, na.rm = TRUE),
            sd_reported_value = sd(.data$reported_value, na.rm = TRUE),
            se_reported_value = sd(.data$reported_value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
        ) %>%
        mutate(experiment_id = experiment_identifier) %>%
        relocate("experiment_id")
    summary_csv_path <- file.path(experiment_output_directory, "summary_data.csv")
    write.csv(condition_summary, summary_csv_path, row.names = FALSE)
    message("Experiment '", experiment_identifier, "': wrote right-before-plot state: ",
            basename(summary_csv_path), " (", nrow(condition_summary), " conditions).")
    combined_summary_frames[[length(combined_summary_frames) + 1]] <- condition_summary

    # ---- pairing completeness for the matrices ------------------------------
    pairing_table <- one_experiment %>% count(.data$replicate, .data$condition_label)
    complete_pairing <- nrow(pairing_table) ==
        length(experiment_replicates) * length(experiment_conditions)
    if (!complete_pairing) {
        present_pairs <- one_experiment %>%
            distinct(.data$replicate, .data$condition_label)
        warning("Experiment '", experiment_identifier, "': incomplete pairing (not every ",
                "condition is present in every replicate). Wrote the summary; SKIPPING the ",
                "matrices and plot for this experiment. Present (replicate, condition) pairs:\n",
                paste(capture.output(print(as.data.frame(present_pairs))), collapse = "\n"))
        next
    }

    # ---- denominator guard: every reported value must be positive -----------
    nonpositive_values <- one_experiment %>%
        filter(!is.na(.data$reported_value)) %>%
        filter(.data$reported_value <= 0) %>%
        select("replicate", "condition_label", "reported_value")
    if (nrow(nonpositive_values) > 0) {
        stop("Experiment '", experiment_identifier, "': a reported value is <= 0, which ",
             "makes a fold/percent-difference denominator undefined:\n",
             paste(capture.output(print(as.data.frame(nonpositive_values))), collapse = "\n"))
    }

    # ---- all-against-all, within each replicate -----------------------------
    # Look up value by (replicate, condition), then form every ordered pair
    # (condition_i, condition_j) within a replicate. i relative to j:
    #   fold_difference(i,j)        = value_i / value_j
    #   percent_difference(i,j)     = (value_i - value_j) / value_j * 100
    value_lookup <- one_experiment %>%
        select("replicate", "condition_label", "reported_value")

    matrix_replicate <- character(0)
    matrix_condition_i <- character(0)
    matrix_condition_j <- character(0)
    matrix_fold_difference <- numeric(0)
    matrix_percent_difference <- numeric(0)
    for (replicate_label in experiment_replicates) {
        replicate_values <- value_lookup %>% filter(.data$replicate == replicate_label)
        for (condition_i_label in experiment_conditions) {
            value_i <- replicate_values$reported_value[replicate_values$condition_label == condition_i_label]
            for (condition_j_label in experiment_conditions) {
                value_j <- replicate_values$reported_value[replicate_values$condition_label == condition_j_label]
                matrix_replicate <- c(matrix_replicate, replicate_label)
                matrix_condition_i <- c(matrix_condition_i, condition_i_label)
                matrix_condition_j <- c(matrix_condition_j, condition_j_label)
                matrix_fold_difference <- c(matrix_fold_difference, value_i / value_j)
                matrix_percent_difference <- c(matrix_percent_difference,
                                                (value_i - value_j) / value_j * 100)
            }
        }
    }
    matrix_per_replicate_long <- data.frame(
        experiment_id = experiment_identifier,
        replicate = matrix_replicate,
        condition_i = matrix_condition_i,
        condition_j = matrix_condition_j,
        fold_difference = matrix_fold_difference,
        percent_difference = matrix_percent_difference,
        stringsAsFactors = FALSE
    )
    matrix_long_csv_path <- file.path(experiment_output_directory, "matrix_per_replicate_long.csv")
    write.csv(matrix_per_replicate_long, matrix_long_csv_path, row.names = FALSE)

    # Average across replicates (the reported matrix values), plus sd/se.
    matrix_summary_long <- matrix_per_replicate_long %>%
        group_by(.data$condition_i, .data$condition_j) %>%
        summarise(
            replicate_count = n(),
            mean_fold_difference = mean(.data$fold_difference),
            sd_fold_difference = sd(.data$fold_difference),
            se_fold_difference = sd(.data$fold_difference) / sqrt(n()),
            mean_percent_difference = mean(.data$percent_difference),
            sd_percent_difference = sd(.data$percent_difference),
            se_percent_difference = sd(.data$percent_difference) / sqrt(n()),
            .groups = "drop"
        ) %>%
        mutate(experiment_id = experiment_identifier) %>%
        relocate("experiment_id")
    matrix_summary_csv_path <- file.path(experiment_output_directory, "matrix_summary_long.csv")
    write.csv(matrix_summary_long, matrix_summary_csv_path, row.names = FALSE)

    # Wide (K x K) mean matrices, rows = condition_i, columns = condition_j.
    fold_matrix_wide <- matrix_summary_long %>%
        select("condition_i", "condition_j", "mean_fold_difference") %>%
        pivot_wider(names_from = "condition_j", values_from = "mean_fold_difference")
    percent_matrix_wide <- matrix_summary_long %>%
        select("condition_i", "condition_j", "mean_percent_difference") %>%
        pivot_wider(names_from = "condition_j", values_from = "mean_percent_difference")
    write.csv(fold_matrix_wide,
              file.path(experiment_output_directory, "fold_difference_matrix_mean_wide.csv"),
              row.names = FALSE)
    write.csv(percent_matrix_wide,
              file.path(experiment_output_directory, "percent_difference_matrix_mean_wide.csv"),
              row.names = FALSE)
    message("Experiment '", experiment_identifier, "': wrote fold and percent-difference ",
            "matrices (", length(experiment_conditions), " x ",
            length(experiment_conditions), ", averaged over ",
            length(experiment_replicates), " replicates).")

    # ---- plot: mean + sd bars, per-replicate points, reference at 100 -------
    plot_frame <- condition_summary %>%
        mutate(condition_label = factor(.data$condition_label, levels = experiment_conditions))
    point_frame <- one_experiment %>%
        mutate(condition_label = factor(.data$condition_label, levels = experiment_conditions))
    aggregate_bar_plot <- ggplot(
        plot_frame,
        aes(x = .data$condition_label, y = .data$mean_reported_value)
    ) +
        geom_col(fill = PLOT_CONFIG$bar_fill, width = PLOT_CONFIG$bar$width,
                 color = PLOT_CONFIG$bar$color, linewidth = PLOT_CONFIG$bar$linewidth) +
        geom_errorbar(
            aes(ymin = pmax(0, .data$mean_reported_value - .data$sd_reported_value),
                ymax = .data$mean_reported_value + .data$sd_reported_value),
            width = PLOT_CONFIG$errorbar$width, linewidth = PLOT_CONFIG$errorbar$linewidth
        ) +
        geom_point(
            data = point_frame,
            aes(x = .data$condition_label, y = .data$reported_value, shape = factor(.data$replicate)),
            position = position_jitter(width = PLOT_CONFIG$point$jitter_width,
                                       seed = PLOT_CONFIG$point$jitter_seed),
            size = PLOT_CONFIG$point$size, fill = PLOT_CONFIG$point$fill,
            color = PLOT_CONFIG$point$color, stroke = PLOT_CONFIG$point$stroke,
            inherit.aes = FALSE
        ) +
        geom_hline(yintercept = REFERENCE_PERCENT_TARGET, linetype = "dashed",
                   color = PLOT_CONFIG$reference_line_color) +
        scale_shape_discrete(name = "Replicate") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(
            x = paste(CONDITION_KEY_COLUMNS, collapse = " | "),
            y = paste0("percent of reference (", NORMALIZATION_COLUMN, ")"),
            title = paste0("Loading aggregate: ", experiment_identifier),
            caption = paste0("Bars = mean +/- sd; points = per-replicate values; ",
                             "dashed line = reference (100%). n = ",
                             length(experiment_replicates), " replicates. ",
                             "Descriptive only; no inferential test.")
        ) +
        theme_classic(base_size = PLOT_CONFIG$base_size, base_family = resolved_font_family) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.caption = element_text(size = 7, hjust = 0))

    plot_output_path <- file.path(experiment_output_directory, "aggregate_bar.pdf")
    plotting_outcome <- tryCatch({
        ggsave(plot_output_path, aggregate_bar_plot, device = PLOT_CONFIG$output$device,
               width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
        "ok"
    }, error = function(e) conditionMessage(e))
    if (identical(plotting_outcome, "ok")) {
        message("Experiment '", experiment_identifier, "': wrote plot ", basename(plot_output_path),
                " (vector PDF; text stays editable in Illustrator).")
    } else {
        warning("Experiment '", experiment_identifier, "': plotting failed (", plotting_outcome,
                "). All CSVs were written before plotting, so the data outputs are complete.")
    }
}

# Combined summary across experiments, for convenience.
if (length(combined_summary_frames) > 0) {
    combined_summary <- bind_rows(combined_summary_frames)
    write.csv(combined_summary, file.path(output_directory, "summary_data_all_experiments.csv"),
              row.names = FALSE)
}

}  # end keyed mode

message("Aggregate script complete. Outputs in: ", output_directory)
