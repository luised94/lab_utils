# ==============================================================================
# ATPase Percent Hydrolysis Analysis - Consolidated
# ==============================================================================
# Date created: 2026-04-11
# Consolidates analysis of ORC ATPase assays across Orc1, Orc3, Orc4, Orc5,
# and Orc6 suppressor mutants in the 4R background.
# Input: Excel files containing ImageJ band intensity measurements.
# Output: Combined raw data CSV, processed data CSV, summary CSV, and
#         percent hydrolysis timecourse plot (PDF).
# Usage: source("260411_orc4r-screen_analyze-sofa-atpase.R")
# ==============================================================================

# ==============================================================================
# GIT STATE REFERENCE (manual-fill at deposit time; no runtime git calls)
# ==============================================================================
# Commit hash:    ____________________________________________
# Branch:         ____________________________________________
# Tag / release:  ____________________________________________
# Snapshot date:  ____________________________________________
# Repository URL: ____________________________________________
# ==============================================================================

library(readxl)  # C3: migrated from library(xlsx)
library(tidyverse)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
message("=== CONFIGURATION ===")

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

TIMEPOINTS <- c(0, 15, 45, 90)

# ------------------------------------------------------------------------------
# Script location (C1): resolve the directory of THIS file under source().
# Only source() invocation is supported. Rscript and interactive paste do not
# set 'ofile' in any call frame, so we stop() with a clear message otherwise.
# ------------------------------------------------------------------------------
script_path_under_source <- NULL
for (frame_index in seq_len(sys.nframe())) {
    candidate_ofile <- sys.frame(frame_index)$ofile
    if (!is.null(candidate_ofile)) {
        script_path_under_source <- candidate_ofile
    }
}
if (is.null(script_path_under_source)) {
    stop(
        "This script must be run via source(\"260411_orc4r-screen_analyze-sofa-atpase.R\").\n",
        "Rscript and interactive invocation are not supported ",
        "(no script path is available to resolve data locations)."
    )
}
SCRIPT_DIRECTORY <- dirname(normalizePath(script_path_under_source))

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
# MC_DROPBOX_PATH is the original (Dropbox) data home. It may be unset in a
# Zenodo deposit, where data sits alongside this script instead. Resolution
# (below, after the registry is defined) picks the live location.
MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")

# ------------------------------------------------------------------------------
# Sample ordering and colors
# ------------------------------------------------------------------------------
# Factor order determines legend and axis ordering in all outputs.
# Biological logic: WT control first, then 4R mutant, then suppressors
# ordered by ORC subunit number.
SAMPLE_DISPLAY_ORDER <- c(
    "WT", "4R",
    "Orc1_E495K_4R", "Orc3_P481L_4R", "Orc4_P225S_4R",
    "Orc5_E104K_4R", "Orc6_E304K_4R"
)

# Hex values taken from RColorBrewer::brewer.pal(7, "Dark2"), assigned
# in SAMPLE_DISPLAY_ORDER so each sample keeps its color across figures.
SAMPLE_COLORS <- c(
    "WT"            = "#1B9E77",
    "4R"            = "#D95F02",
    "Orc1_E495K_4R" = "#7570B3",
    "Orc3_P481L_4R" = "#E7298A",
    "Orc4_P225S_4R" = "#66A61E",
    "Orc5_E104K_4R" = "#E6AB02",
    "Orc6_E304K_4R" = "#A6761D"
)

# Relative paths from BASE_EXPERIMENT_DIRECTORY to each Excel file.
# Some files contain multiple experiments on different sheets.
FILE_ORC3_PL_1 <- file.path(
    "2020_09_03 ATPase Analysis of 4R supps",
    "2020_07_13 Orc3 P481L 4R",
    "2020_07_13 Orc3 P481L 4R ATPase.xlsx"
)
FILE_ORC3_PL_2 <- file.path(
    "2020_09_03 ATPase Analysis of 4R supps",
    "2020_07_27 Orc3 P481L 4R #2",
    "2020_07_27 Orc3 P481L 4R ATPase #2.xlsx"
)
FILE_ORC3_PL_3_ORC4_PS_1 <- file.path(
    "2020_09_03 ATPase Analysis of 4R supps",
    "2020_08_13 Orc3 PL 4R #3 and Orc4 PS 4R",
    "2020_08_13 Orc3 PL 4R #3 and Orc4 PS 4R #1.xlsx"
)
FILE_ORC4_PS_2_3 <- file.path(
    "2020_09_03 ATPase Analysis of 4R supps",
    "2020_08_24 Orc4 PS 4R #2 and #3",
    "2020_08_24 Orc4 PS 4R #2 and #3.xlsx"
)
FILE_ORC1_EK <- file.path(
    "2020_09_03 ATPase Analysis of 4R supps",
    "2020_09_21 Orc1 EK 4R",
    "2020_09_21 Orc1 EK 4R.xlsx"
)
FILE_ORC5_ORC6_EK <- file.path(
    "2021_07_29 ORC4R O5 and O6 EK",
    "2022_03_01 ORC4R O5,O6 EK.xlsx"
)

# ------------------------------------------------------------------------------
# Experiment registry
# ------------------------------------------------------------------------------
# --- Physical context ---
# Each ATPase experiment resolves radiolabeled ATP and ADP on a TLC
# (thin-layer chromatography) plate. ADP migrates further from the
# origin than ATP. The plate is exposed to a phosphor screen, scanned,
# and quantified in ImageJ by selecting each band with the wand tool.
#
# ImageJ records one row per band selection. Bands alternate between
# ADP and ATP for each lane (one lane = one reaction timepoint). The
# order in which ADP vs ATP was selected with the wand varied between
# experiment batches, producing two conventions:
#
#   Layouts A/B (2020 experiments): odd ImageJ index = ADP, even = ATP
#   Layout C   (2021 experiments): even ImageJ index = ADP, odd = ATP
#
# The No_ORC control (no enzyme, measures spontaneous hydrolysis) was
# recorded differently across batches:
#
# Layout A: 16 raw rows. No_ORC occupies rows 1-4 (all timepoints).
#           Odd ImageJ index = ADP, next row = ATP.
#           Background = No_ORC at t=90 (processed row 4).
#           Used only for 2020_07_13 (first experiment).
#
# Layout B: 17 raw rows. Four samples in rows 1-16 (4 timepoints each).
#           Single No_ORC at row 17 (t=90 only).
#           Odd ImageJ index = ADP, next row = ATP.
#           Background = row 17.
#
# Layout C: 13 raw rows. Single No_ORC at row 1 (t=90 only).
#           Three samples in rows 2-13 (4 timepoints each).
#           Even ImageJ index = ADP, previous row = ATP.
#           Background = row 1.

EXPERIMENT_REGISTRY <- list(
    list(
        label = "2020_07_13",
        relative_file_path = FILE_ORC3_PL_1,
        sheet_index = 4,
        layout = "A",
        sample_names = c("No_ORC", "WT", "WT_DNA", "Orc3_P481L_4R"),
        expected_raw_row_count = 32
    ),
    list(
        label = "2020_07_27",
        relative_file_path = FILE_ORC3_PL_2,
        sheet_index = 4,
        layout = "B",
        sample_names = c("WT", "WT_DNA", "4R", "Orc3_P481L_4R"),
        expected_raw_row_count = 34
    ),
    list(
        label = "2020_08_13a",
        relative_file_path = FILE_ORC3_PL_3_ORC4_PS_1,
        sheet_index = 4,
        layout = "B",
        sample_names = c("WT", "WT_DNA", "4R", "Orc3_P481L_4R"),
        expected_raw_row_count = 34
    ),
    list(
        label = "2020_08_13b",
        relative_file_path = FILE_ORC3_PL_3_ORC4_PS_1,
        sheet_index = 5,
        layout = "B",
        sample_names = c("WT", "WT_DNA", "4R", "Orc4_P225S_4R"),
        expected_raw_row_count = 34
    ),
    list(
        label = "2020_08_24",
        relative_file_path = FILE_ORC4_PS_2_3,
        sheet_index = 4,
        layout = "B",
        sample_names = c("WT", "4R", "Orc4_P225S_4R", "Orc4_P225S_4R"),
        expected_raw_row_count = 34
    ),
    list(
        label = "2020_09_21a",
        relative_file_path = FILE_ORC1_EK,
        sheet_index = 4,
        layout = "B",
        sample_names = c("WT", "WT_DNA", "4R", "Orc1_E495K_4R"),
        expected_raw_row_count = 34
    ),
    list(
        label = "2020_09_21b",
        relative_file_path = FILE_ORC1_EK,
        sheet_index = 5,
        layout = "B",
        sample_names = c("WT", "WT_DNA", "4R", "Orc1_E495K_4R"),
        expected_raw_row_count = 34
    ),
    list(
        label = "2021_07_29a",
        relative_file_path = FILE_ORC5_ORC6_EK,
        sheet_index = 6,
        layout = "C",
        sample_names = c("WT", "4R", "Orc6_E304K_4R"),
        expected_raw_row_count = 26
    ),
    list(
        label = "2021_07_29b",
        relative_file_path = FILE_ORC5_ORC6_EK,
        sheet_index = 7,
        layout = "C",
        sample_names = c("WT", "4R", "Orc6_E304K_4R"),
        expected_raw_row_count = 26
    ),
    list(
        label = "2021_07_29c",
        relative_file_path = FILE_ORC5_ORC6_EK,
        sheet_index = 8,
        layout = "C",
        sample_names = c("WT", "4R", "Orc6_E304K_4R"),
        expected_raw_row_count = 26
    ),
    list(
        label = "2021_07_29d",
        relative_file_path = FILE_ORC5_ORC6_EK,
        sheet_index = 9,
        layout = "C",
        sample_names = c("WT", "4R", "Orc5_E104K_4R"),
        expected_raw_row_count = 26
    ),
    list(
        label = "2021_07_29e",
        relative_file_path = FILE_ORC5_ORC6_EK,
        sheet_index = 10,
        layout = "C",
        sample_names = c("WT", "4R", "Orc5_E104K_4R"),
        expected_raw_row_count = 26
    ),
    list(
        label = "2021_07_29f",
        relative_file_path = FILE_ORC5_ORC6_EK,
        sheet_index = 11,
        layout = "C",
        sample_names = c("WT", "4R", "Orc5_E104K_4R"),
        expected_raw_row_count = 26
    )
)

# ------------------------------------------------------------------------------
# Path resolution (C1): script-relative (Zenodo co-located) -> MC_DROPBOX_PATH
# -> stop() naming both. Probe with the first registry file's relative path.
# ------------------------------------------------------------------------------
probe_relative_file_path <- EXPERIMENT_REGISTRY[[1]]$relative_file_path
script_relative_base_directory <- SCRIPT_DIRECTORY
if (nchar(MC_DROPBOX_PATH) > 0) {
    dropbox_base_directory <- file.path(
        MC_DROPBOX_PATH, "Lab", "Experiments", "ATPase"
    )
} else {
    dropbox_base_directory <- NA_character_
}

if (file.exists(file.path(script_relative_base_directory, probe_relative_file_path))) {
    BASE_EXPERIMENT_DIRECTORY <- script_relative_base_directory
    OUTPUT_DIRECTORY <- SCRIPT_DIRECTORY
} else if (!is.na(dropbox_base_directory) &&
    file.exists(file.path(dropbox_base_directory, probe_relative_file_path))) {
    BASE_EXPERIMENT_DIRECTORY <- dropbox_base_directory
    OUTPUT_DIRECTORY <- file.path(
        dropbox_base_directory,
        "2020_09_03 ATPase Analysis of 4R supps",
        "consolidated_analysis"
    )
} else {
    stop(
        "ATPase input files not found in either supported location:\n",
        "  script-relative: ",
        file.path(script_relative_base_directory, probe_relative_file_path), "\n",
        "  MC_DROPBOX_PATH:  ",
        if (is.na(dropbox_base_directory)) {
            "<MC_DROPBOX_PATH not set>"
        } else {
            file.path(dropbox_base_directory, probe_relative_file_path)
        }
    )
}

if (!dir.exists(OUTPUT_DIRECTORY)) {
    dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)
    message("Created output directory: ", OUTPUT_DIRECTORY)
}

# ------------------------------------------------------------------------------
# Validate file paths
# ------------------------------------------------------------------------------
message("Validating file paths...")
all_experiment_labels <- vapply(
    EXPERIMENT_REGISTRY, function(entry) entry$label, character(1)
)
if (anyDuplicated(all_experiment_labels)) {
    stop(
        "Duplicate experiment labels in registry: ",
        paste(all_experiment_labels[duplicated(all_experiment_labels)], collapse = ", ")
    )
}

unique_relative_paths <- unique(vapply(
    EXPERIMENT_REGISTRY, function(entry) entry$relative_file_path, character(1)
))

for (relative_path in unique_relative_paths) {
    full_path <- file.path(BASE_EXPERIMENT_DIRECTORY, relative_path)
    if (!file.exists(full_path)) {
        stop("File not found: ", full_path)
    }
}

message("All ", length(unique_relative_paths), " input files found.")
message("Experiment registry contains ", length(EXPERIMENT_REGISTRY), " experiments.")
# Validate that every sample name in the registry is a recognized sample.
# Catches typos when adding new experiments.
all_registry_sample_names <- unique(unlist(
    lapply(EXPERIMENT_REGISTRY, function(entry) entry$sample_names)
))
unknown_samples <- setdiff(
    all_registry_sample_names,
    c(SAMPLE_DISPLAY_ORDER, "No_ORC", "WT_DNA")
)
if (length(unknown_samples) > 0) {
    stop(
        "Unknown sample names in registry: ",
        paste(unknown_samples, collapse = ", "),
        "\nAdd to SAMPLE_DISPLAY_ORDER or to the exclusion list."
    )
}
message("All registry sample names are recognized.")
message("Output directory: ", OUTPUT_DIRECTORY)
message("=== CONFIGURATION COMPLETE ===")

# ==============================================================================
# DATA LOADING - READ AND VALIDATE RAW SHEETS
# ==============================================================================
message("=== DATA LOADING ===")

EXPECTED_NUMBER_OF_COLUMNS <- 2
RAW_COLUMN_NAMES <- c("imagej_index", "intensity")

raw_data_list <- vector(mode = "list", length = length(EXPERIMENT_REGISTRY))
names(raw_data_list) <- vapply(
    EXPERIMENT_REGISTRY, function(entry) entry$label, character(1)
)

for (registry_index in seq_along(EXPERIMENT_REGISTRY)) {
    current_experiment <- EXPERIMENT_REGISTRY[[registry_index]]
    current_label <- current_experiment$label
    current_full_path <- file.path(
        BASE_EXPERIMENT_DIRECTORY, current_experiment$relative_file_path
    )

    message(
        "  Reading ", current_label,
        " (sheet ", current_experiment$sheet_index, ")..."
    )

    # C3: readxl replaces xlsx. read_excel(sheet=) is 1-based, matching the
    # former sheetIndex= numbering. col_names = FALSE mirrors header = FALSE.
    current_raw_data <- readxl::read_excel(
        current_full_path,
        col_names = FALSE,
        sheet = current_experiment$sheet_index
    )
    current_raw_data <- as.data.frame(current_raw_data)

    # -- Validate dimensions BEFORE assigning column names (C3) --
    # readxl returns a tibble; an unexpected column count would make the
    # length-2 RAW_COLUMN_NAMES assignment crash with an opaque length
    # error. Checking dimensions first lets the clear assertion fire.
    if (nrow(current_raw_data) != current_experiment$expected_raw_row_count) {
        stop(
            "Row count mismatch for ", current_label, ": ",
            "expected ", current_experiment$expected_raw_row_count,
            ", got ", nrow(current_raw_data)
        )
    }
    if (ncol(current_raw_data) != EXPECTED_NUMBER_OF_COLUMNS) {
        stop(
            "Column count mismatch for ", current_label, ": ",
            "expected ", EXPECTED_NUMBER_OF_COLUMNS,
            ", got ", ncol(current_raw_data)
        )
    }

    colnames(current_raw_data) <- RAW_COLUMN_NAMES

    # -- Validate content --
    if (anyNA(current_raw_data$intensity)) {
        stop("NA values found in intensity column for ", current_label)
    }
    if (any(current_raw_data$intensity < 0)) {
        stop("Negative intensity values found for ", current_label)
    }

    # ImageJ indices should be sequential 1:N. A gap or reordering
    # indicates the Excel sheet was modified or rows were deleted.
    # NOTE (C3): readxl parses imagej_index as double; as.integer() coerces
    # it back so identical(integer, integer) still fires for any gap/reorder.
    expected_indices <- seq_len(nrow(current_raw_data))
    if (!identical(as.integer(current_raw_data$imagej_index), expected_indices)) {
        stop(
            "Non-sequential imagej_index for ", current_label,
            ". Expected 1:", nrow(current_raw_data),
            ". Check Excel sheet for missing or reordered rows."
        )
    }

    # Raw row count must be even: every ADP band has a paired ATP band.
    if (nrow(current_raw_data) %% 2 != 0) {
        stop(
            "Odd number of rows (", nrow(current_raw_data), ") for ",
            current_label, ". ADP/ATP bands must come in pairs."
        )
    }

    raw_data_list[[current_label]] <- current_raw_data
}

message(
    "All ", length(raw_data_list), " sheets loaded and validated."
)
message("=== DATA LOADING COMPLETE ===")

# ==============================================================================
# DATA RESHAPING - PAIR ADP/ATP AND ASSIGN METADATA
# ==============================================================================
message("=== DATA RESHAPING ===")

reshaped_data_list <- vector(mode = "list", length = length(EXPERIMENT_REGISTRY))

for (registry_index in seq_along(EXPERIMENT_REGISTRY)) {
    current_experiment <- EXPERIMENT_REGISTRY[[registry_index]]
    current_label <- current_experiment$label
    current_layout <- current_experiment$layout
    current_sample_names <- current_experiment$sample_names
    current_raw_data <- raw_data_list[[current_label]]

    message("  Reshaping ", current_label, " (layout ", current_layout, ")...")

    # -- Pair ADP and ATP intensities based on layout --
    # The wand selection order in ImageJ determines which raw rows are
    # ADP vs ATP. See layout documentation in the configuration section.
    if (current_layout %in% c("A", "B")) {
        # Odd-indexed rows are ADP, the following even-indexed row is ATP.
        odd_row_positions <- seq(1, nrow(current_raw_data), by = 2)
        paired_data <- data.frame(
            adp_intensity = current_raw_data$intensity[odd_row_positions],
            atp_intensity = current_raw_data$intensity[odd_row_positions + 1]
        )
    } else if (current_layout == "C") {
        # Even-indexed rows are ADP, the preceding odd-indexed row is ATP.
        even_row_positions <- seq(2, nrow(current_raw_data), by = 2)
        paired_data <- data.frame(
            adp_intensity = current_raw_data$intensity[even_row_positions],
            atp_intensity = current_raw_data$intensity[even_row_positions - 1]
        )
    } else {
        stop("Unknown layout '", current_layout, "' for ", current_label)
    }

    # -- Assign timepoints --
    if (current_layout == "A") {
        # All samples have all 4 timepoints. No standalone No_ORC row.
        paired_data$timepoint <- rep(TIMEPOINTS, times = length(current_sample_names))
    } else if (current_layout == "B") {
        # 4 samples with all timepoints, then single No_ORC at t=90.
        paired_data$timepoint <- c(
            rep(TIMEPOINTS, times = length(current_sample_names)), 90
        )
    } else if (current_layout == "C") {
        # Single No_ORC at t=90 first, then 3 samples with all timepoints.
        paired_data$timepoint <- c(
            90, rep(TIMEPOINTS, times = length(current_sample_names))
        )
    }

    # -- Assign sample names --
    if (current_layout == "A") {
        # Each sample occupies 4 consecutive rows.
        paired_data$sample <- rep(current_sample_names, each = 4)
    } else if (current_layout == "B") {
        # 4 samples with 4 rows each, then one No_ORC row.
        paired_data$sample <- c(
            rep(current_sample_names, each = 4), "No_ORC"
        )
    } else if (current_layout == "C") {
        # One No_ORC row, then 3 samples with 4 rows each.
        paired_data$sample <- c(
            "No_ORC", rep(current_sample_names, each = 4)
        )
    }

    paired_data$experiment_label <- current_label

    # -- Validate paired data dimensions --
    if (current_layout == "A") {
        expected_paired_rows <- length(current_sample_names) * 4
    } else if (current_layout == "B") {
        expected_paired_rows <- length(current_sample_names) * 4 + 1
    } else if (current_layout == "C") {
        expected_paired_rows <- length(current_sample_names) * 4 + 1
    }

    if (nrow(paired_data) != expected_paired_rows) {
        stop(
            "Paired row count mismatch for ", current_label, ": ",
            "expected ", expected_paired_rows, ", got ", nrow(paired_data)
        )
    }

    # -- Validate ADP/ATP assignment --
    # If ADP and ATP columns were swapped, percent would be near 1.0
    # for No_ORC (almost all signal in the "ADP" slot). Correct assignment
    # gives No_ORC percent near 0.05-0.20 (low spontaneous hydrolysis).
    quick_percent_check <- paired_data$adp_intensity /
        (paired_data$adp_intensity + paired_data$atp_intensity)
    if (any(quick_percent_check < 0 | quick_percent_check > 1, na.rm = TRUE)) {
        stop(
            "percent_adp outside [0, 1] for ", current_label,
            ". ADP/ATP intensities may be swapped or negative."
        )
    }

    # Each experiment must have exactly one No_ORC measurement at t=90
    # for background subtraction downstream.
    no_orc_t90_count <- sum(
        paired_data$sample == "No_ORC" & paired_data$timepoint == 90
    )
    if (no_orc_t90_count != 1) {
        stop(
            "Expected 1 No_ORC at t=90 for ", current_label,
            ", found ", no_orc_t90_count,
            ". Check sample_names and layout in the registry."
        )
    }

    reshaped_data_list[[registry_index]] <- paired_data
}

combined_raw_data <- do.call(rbind, reshaped_data_list)
rownames(combined_raw_data) <- NULL

message("Combined raw data: ", nrow(combined_raw_data), " rows x ",
    ncol(combined_raw_data), " columns.")

# -- Per-experiment row counts --
message("Row counts per experiment:")
print(table(combined_raw_data$experiment_label))

# -- Write to CSV --
combined_raw_csv_path <- file.path(OUTPUT_DIRECTORY, "combined_raw_data.csv")
if (!file.exists(combined_raw_csv_path) || OVERWRITE_CSVS) {
    write.csv(combined_raw_data, combined_raw_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(combined_raw_csv_path))
} else {
    message("Skipped CSV (already exists): ", basename(combined_raw_csv_path))
}

message("=== DATA RESHAPING COMPLETE ===")

# ==============================================================================
# PROCESSING - PERCENT HYDROLYSIS AND BACKGROUND SUBTRACTION
# ==============================================================================
message("=== PROCESSING ===")

processed_data <- combined_raw_data
processed_data$percent_adp <- processed_data$adp_intensity /
    (processed_data$adp_intensity + processed_data$atp_intensity)

# Percent ADP must be in [0, 1]. Values outside this range indicate
# a problem in the ADP/ATP pairing or negative intensities.
if (any(processed_data$percent_adp < 0 | processed_data$percent_adp > 1)) {
    bad_rows <- which(
        processed_data$percent_adp < 0 | processed_data$percent_adp > 1
    )
    stop(
        "percent_adp outside [0, 1] at rows: ",
        paste(head(bad_rows, 5), collapse = ", "),
        ". Check ADP/ATP pairing for these experiments."
    )
}

# Background subtraction: subtract the No_ORC percent_adp at t=90 from all
# rows within the same experiment.
processed_data$percent_adp_corrected <- NA

for (current_label in unique(processed_data$experiment_label)) {
    current_rows <- processed_data$experiment_label == current_label
    background_row <- current_rows &
        processed_data$sample == "No_ORC" &
        processed_data$timepoint == 90

    background_count <- sum(background_row)
    if (background_count != 1) {
        stop(
            "Expected exactly 1 No_ORC t=90 row for ", current_label,
            ", found ", background_count
        )
    }

    background_value <- processed_data$percent_adp[background_row]
    processed_data$percent_adp_corrected[current_rows] <-
        processed_data$percent_adp[current_rows] - background_value
}

if (anyNA(processed_data$percent_adp_corrected)) {
    stop("NA values found in percent_adp_corrected after background subtraction.")
}

# Sanity check: WT at t=0 should be near zero after background subtraction.
# A large value (above 0.15) suggests the wrong row was used as background
# or the ADP/ATP assignment is incorrect for that experiment.
wt_t0_corrected <- processed_data$percent_adp_corrected[
    processed_data$sample == "WT" & processed_data$timepoint == 0
]
wt_t0_max_absolute <- max(abs(wt_t0_corrected))
if (wt_t0_max_absolute > 0.15) {
    warning(
        "Largest |WT t=0 corrected| is ", round(wt_t0_max_absolute, 4),
        " (threshold: 0.15). Inspect background subtraction."
    )
} else {
    message(
        "WT t=0 sanity check passed (max |corrected| = ",
        round(wt_t0_max_absolute, 4), ")."
    )
}

# Note: negative corrected values are possible when a sample has less
# apparent ADP than the No_ORC background. These are biologically
# meaningless but are retained for transparency. They typically occur
# at t=0 and are small in magnitude.

message("Processing complete. ", nrow(processed_data), " rows.")

# -- Write processed data CSV --
processed_data_csv_path <- file.path(OUTPUT_DIRECTORY, "processed_data.csv")
if (!file.exists(processed_data_csv_path) || OVERWRITE_CSVS) {
    write.csv(processed_data, processed_data_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(processed_data_csv_path))
} else {
    message("Skipped CSV (already exists): ", basename(processed_data_csv_path))
}

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
message("=== SUMMARY STATISTICS ===")

# Filter out control samples before summarizing.
samples_to_exclude <- c("No_ORC", "WT_DNA")
plotting_data <- processed_data[
    !(processed_data$sample %in% samples_to_exclude),
]

plotting_data$sample <- factor(
    plotting_data$sample,
    levels = SAMPLE_DISPLAY_ORDER,
    ordered = TRUE
)

message(
    "Excluded samples: ", paste(samples_to_exclude, collapse = ", "),
    ". Rows remaining: ", nrow(plotting_data)
)

summary_data <- plotting_data %>%
    group_by(timepoint, sample) %>%
    summarise(
        mean_percent_adp_corrected = mean(percent_adp_corrected, na.rm = TRUE),
        sd_percent_adp_corrected = sd(percent_adp_corrected, na.rm = TRUE),
        replicate_count = n(),
        .groups = "drop"
    )

message("Summary computed: ", nrow(summary_data), " rows.")
message("Replicate counts per sample:")
replicate_counts_at_t0 <- summary_data[
    summary_data$timepoint == 0, c("sample", "replicate_count")
]
print(replicate_counts_at_t0)

# Replicate counts should be consistent across timepoints for each sample.
# If they differ, rows were lost during processing.
replicate_count_check <- plotting_data %>%
    group_by(sample) %>%
    summarise(
        timepoints_observed = n_distinct(timepoint),
        .groups = "drop"
    )
samples_missing_timepoints <- replicate_count_check$sample[
    replicate_count_check$timepoints_observed != length(TIMEPOINTS)
]
if (length(samples_missing_timepoints) > 0) {
    warning(
        "Samples with fewer than ", length(TIMEPOINTS), " timepoints: ",
        paste(samples_missing_timepoints, collapse = ", "),
        ". Rows may have been lost during processing."
    )
}

# WT and 4R should have more replicates than individual mutants since
# they were measured in nearly every experiment.
wt_replicate_count <- replicate_counts_at_t0$replicate_count[
    replicate_counts_at_t0$sample == "WT"
]
if (wt_replicate_count < 5) {
    warning(
        "WT replicate count is ", wt_replicate_count,
        ", expected at least 5. Data may be incomplete."
    )
}

# -- Write summary CSV --
summary_csv_path <- file.path(OUTPUT_DIRECTORY, "summary_data.csv")
if (!file.exists(summary_csv_path) || OVERWRITE_CSVS) {
    write.csv(summary_data, summary_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(summary_csv_path))
} else {
    message("Skipped CSV (already exists): ", basename(summary_csv_path))
}

# ==============================================================================
# PLOTTING
# ==============================================================================
message("=== PLOTTING ===")

atpase_timecourse_plot <- ggplot(
    summary_data,
    aes(
        x = timepoint,
        y = mean_percent_adp_corrected,
        color = sample
    )
) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
        aes(
            ymin = mean_percent_adp_corrected - sd_percent_adp_corrected,
            ymax = mean_percent_adp_corrected + sd_percent_adp_corrected
        ),
        width = 3,
        linewidth = 0.5
    ) +
    scale_color_manual(values = SAMPLE_COLORS) +
    labs(
        title = "ORC ATPase Timecourse",
        x = "Time (min)",
        y = "Percent hydrolysis (background-corrected)",
        color = "Sample"
    ) +
    theme_classic(base_size = 13) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
    )

plot_pdf_path <- file.path(OUTPUT_DIRECTORY, "atpase_timecourse_plot.pdf")
if (!file.exists(plot_pdf_path) || OVERWRITE_PLOTS) {
    ggsave(plot_pdf_path, atpase_timecourse_plot, width = 8, height = 5)
    message("Saved plot: ", basename(plot_pdf_path))
} else {
    message("Skipped plot (already exists): ", basename(plot_pdf_path))
}

# ==============================================================================
# COMPLETE
# ==============================================================================
message("=== SCRIPT COMPLETE ===")
message("Output directory: ", OUTPUT_DIRECTORY)
message("Files written:")
for (output_file in list.files(OUTPUT_DIRECTORY)) {
    message("  ", output_file)
}
