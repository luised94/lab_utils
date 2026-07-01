# ==============================================================================
# ATPase Percent Hydrolysis Analysis - Consolidated
# ==============================================================================
# Date created: 2026-04-11
# Consolidates analysis of ORC ATPase assays across Orc1, Orc3, Orc4, Orc5,
# and Orc6 suppressor mutants in the 4R background.
# Input: Excel files containing ImageJ band intensity measurements.
# Output: Combined raw data CSV, processed data CSV, summary CSV, statistics
#         CSVs (Thread 3), and percent hydrolysis timecourse plot (PDF).
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
library(nlme)    # C17 (Thread 3): mixed models; ships with R, recorded via renv (C0d)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
message("=== CONFIGURATION ===")

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

TIMEPOINTS <- c(0, 15, 45, 90)

# Thread 3: subtraction-only negative-value tripwire (NOT a floor gate). A wild
# negative signals an upstream sign/load error and must STOP; small near-zero
# scatter passes silently. Bound mirrors the existing WT-t0 background threshold
# (0.15). See DECISION_LOG [Thread 3]: ATPase does NOT floor (differs from loading).
ATPASE_NEGATIVE_SANITY_BOUND <- -0.15

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
# (see full layout documentation retained below)
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
# Path resolution (C1)
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

    current_raw_data <- readxl::read_excel(
        current_full_path,
        col_names = FALSE,
        sheet = current_experiment$sheet_index
    )
    current_raw_data <- as.data.frame(current_raw_data)

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

    if (anyNA(current_raw_data$intensity)) {
        stop("NA values found in intensity column for ", current_label)
    }
    if (any(current_raw_data$intensity < 0)) {
        stop("Negative intensity values found for ", current_label)
    }

    expected_indices <- seq_len(nrow(current_raw_data))
    if (!identical(as.integer(current_raw_data$imagej_index), expected_indices)) {
        stop(
            "Non-sequential imagej_index for ", current_label,
            ". Expected 1:", nrow(current_raw_data),
            ". Check Excel sheet for missing or reordered rows."
        )
    }

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
        odd_row_positions <- seq(1, nrow(current_raw_data), by = 2)
        paired_data <- data.frame(
            adp_intensity = current_raw_data$intensity[odd_row_positions],
            atp_intensity = current_raw_data$intensity[odd_row_positions + 1]
        )
    } else if (current_layout == "C") {
        even_row_positions <- seq(2, nrow(current_raw_data), by = 2)
        paired_data <- data.frame(
            adp_intensity = current_raw_data$intensity[even_row_positions],
            atp_intensity = current_raw_data$intensity[even_row_positions - 1]
        )
    } else {
        stop("Unknown layout '", current_layout, "' for ", current_label)
    }

    if (current_layout == "A") {
        paired_data$timepoint <- rep(TIMEPOINTS, times = length(current_sample_names))
    } else if (current_layout == "B") {
        paired_data$timepoint <- c(
            rep(TIMEPOINTS, times = length(current_sample_names)), 90
        )
    } else if (current_layout == "C") {
        paired_data$timepoint <- c(
            90, rep(TIMEPOINTS, times = length(current_sample_names))
        )
    }

    if (current_layout == "A") {
        paired_data$sample <- rep(current_sample_names, each = 4)
    } else if (current_layout == "B") {
        paired_data$sample <- c(
            rep(current_sample_names, each = 4), "No_ORC"
        )
    } else if (current_layout == "C") {
        paired_data$sample <- c(
            "No_ORC", rep(current_sample_names, each = 4)
        )
    }

    paired_data$experiment_label <- current_label

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

    quick_percent_check <- paired_data$adp_intensity /
        (paired_data$adp_intensity + paired_data$atp_intensity)
    if (any(quick_percent_check < 0 | quick_percent_check > 1, na.rm = TRUE)) {
        stop(
            "percent_adp outside [0, 1] for ", current_label,
            ". ADP/ATP intensities may be swapped or negative."
        )
    }

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

message("Row counts per experiment:")
print(table(combined_raw_data$experiment_label))

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

# ------------------------------------------------------------------------------
# THREAD 3: ATPase NEGATIVE-VALUE HANDLING -- DETECT + REPORT + TRIPWIRE.
# DELIBERATELY DIFFERENT FROM LOADING (which FLOORED). Rationale (see
# DECISION_LOG [Thread 3] + STATISTICAL_METHODS.md): percent_adp_corrected is
# analysed by mixed model + per-timepoint DIFFERENCES -- both subtraction-based,
# negative-safe, NO division by a near-zero quantity. Small negatives are
# EXPECTED noise around a true zero in the near-zero 4R/suppressor group, and
# that noise is DATA the model must absorb. Flooring would erase exactly the
# variance the model needs AND would bias the near-zero group mean upward (not
# cleanly one-directional safe here). So: NO floor, NO floored column, NO pmax.
# (a) DETECT + REPORT; (b) KEEP RAW; (c) SANITY TRIPWIRE (membership-checked
# before value-checked, per the Thread-1 vacuous-guard lesson).
# (d) DIVISION-COLUMN CHECK: the only division is adp/(adp+atp), bounded [0,1]
#     BEFORE subtraction; there is NO fold-over-4R / percent-of-WT column. None
#     is introduced. (Verified against the received file.)
# ------------------------------------------------------------------------------
stopifnot("percent_adp_corrected" %in% names(processed_data))   # membership BEFORE value
negative_corrected_rows <- processed_data[
    processed_data$percent_adp_corrected < 0,
    c("sample", "experiment_label", "timepoint",
      "adp_intensity", "atp_intensity", "percent_adp", "percent_adp_corrected")
]
message(
    "Negative percent_adp_corrected rows detected: ",
    nrow(negative_corrected_rows),
    " (NOT floored; reported for transparency; flow raw into all models/tests)."
)
if (nrow(negative_corrected_rows) > 0) {
    print(negative_corrected_rows)
}
# (c) tripwire: a wild negative => upstream sign/load error => STOP. Small
#     near-zero scatter passes silently.
if (any(processed_data$percent_adp_corrected < ATPASE_NEGATIVE_SANITY_BOUND)) {
    stop(
        "Implausibly large negative percent_adp_corrected (< ",
        ATPASE_NEGATIVE_SANITY_BOUND,
        ") detected -- this signals an upstream sign/load error, not near-zero ",
        "scatter. Inspect the offending rows above."
    )
}
message(
    "Negative-value sanity tripwire passed (most-negative = ",
    round(min(processed_data$percent_adp_corrected), 5),
    " >= bound ", ATPASE_NEGATIVE_SANITY_BOUND, ")."
)

negative_values_csv_path <- file.path(OUTPUT_DIRECTORY, "atpase_negative_values.csv")
if (!file.exists(negative_values_csv_path) || OVERWRITE_CSVS) {
    write.csv(negative_corrected_rows, negative_values_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(negative_values_csv_path))
}

# Note: negative corrected values are biologically meaningless (a sample reading
# below the No_ORC background) but are RETAINED RAW for the reasons above.

message("Processing complete. ", nrow(processed_data), " rows.")

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

wt_replicate_count <- replicate_counts_at_t0$replicate_count[
    replicate_counts_at_t0$sample == "WT"
]
if (wt_replicate_count < 5) {
    warning(
        "WT replicate count is ", wt_replicate_count,
        ", expected at least 5. Data may be incomplete."
    )
}

summary_csv_path <- file.path(OUTPUT_DIRECTORY, "summary_data.csv")
if (!file.exists(summary_csv_path) || OVERWRITE_CSVS) {
    write.csv(summary_data, summary_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(summary_csv_path))
} else {
    message("Skipped CSV (already exists): ", basename(summary_csv_path))
}

# ==============================================================================
# STATISTICAL ANALYSIS (Thread 3: C15 diagnostics, C16 per-timepoint contrasts,
# C17 pooled mixed models, C18a results CSVs)
# ==============================================================================
message("=== STATISTICAL ANALYSIS ===")

# ------------------------------------------------------------------------------
# Build the analysis frame (exclude controls). Block = experiment_label.
# ------------------------------------------------------------------------------
atpase_analysis_data <- processed_data[
    !(processed_data$sample %in% c("No_ORC", "WT_DNA")),
]

# n=3 HONESTY caveat (I4), reused verbatim from the consolidated log:
N3_CAVEAT <- paste0(
    "Per-cell replicate counts are 2-3 for the mutants; normality CANNOT be ",
    "meaningfully tested (Shapiro-Wilk at these n has almost no power and is ",
    "NOT a license to assume normality). Justification rests on the ",
    "PAIRED/BLOCKED design plus the assay's established track record -- NOT on ",
    "a normality test passing. Shapiro is DECORATIVE."
)

# ------------------------------------------------------------------------------
# C15: ASSUMPTION DIAGNOSTICS
# Per-cell {n, mean, sd}; pooled + per-timepoint within-cell-residual Shapiro
# (DECORATIVE). Unequal-n and grossly unequal variance noted explicitly.
# ------------------------------------------------------------------------------
message("--- C15 diagnostics ---")
message(
    "NOTE (C15): per-cell replicate counts are UNEQUAL (WT n=13, 4R n=12, ",
    "Orc1 n=2, Orc3/Orc4/Orc5/Orc6 n=3 lanes) and within-cell variances are ",
    "GROSSLY unequal (WT cells large variance; near-floor mutant cells ~0, ",
    "several exactly 0). Bartlett/Levene are UNDEFINED on zero-variance cells ",
    "and are OMITTED (as in loading). The homoscedastic-residual assumption of ",
    "the C17 lme is therefore VIOLATED; the lme is reported WITH that caveat, ",
    "and varIdent(form = ~1 | sample) is noted as a documented, NOT-adopted, ",
    "heteroscedastic alternative."
)

atpase_analysis_data <- atpase_analysis_data %>%
    group_by(sample, timepoint) %>%
    mutate(
        cell_n = n(),
        cell_mean = mean(percent_adp_corrected),
        cell_sd = sd(percent_adp_corrected),
        within_cell_residual = percent_adp_corrected - cell_mean
    ) %>%
    ungroup() %>%
    as.data.frame()

atpase_cell_stats <- atpase_analysis_data %>%
    group_by(sample, timepoint) %>%
    summarise(
        n = dplyr::first(cell_n),
        mean_corrected = dplyr::first(cell_mean),
        sd_corrected = dplyr::first(cell_sd),
        .groups = "drop"
    ) %>%
    as.data.frame()

cellstats_csv_path <- file.path(OUTPUT_DIRECTORY, "atpase_diagnostics_cellstats.csv")
if (!file.exists(cellstats_csv_path) || OVERWRITE_CSVS) {
    write.csv(atpase_cell_stats, cellstats_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(cellstats_csv_path))
}

# Shapiro on pooled within-cell residuals + per-timepoint splits (DECORATIVE).
shapiro_rows <- data.frame()
pooled_shapiro <- tryCatch(
    shapiro.test(atpase_analysis_data$within_cell_residual),
    error = function(e) NULL
)
shapiro_rows <- rbind(shapiro_rows, data.frame(
    scope = "pooled_within_cell_residuals",
    shapiro_W = if (is.null(pooled_shapiro)) NA_real_ else unname(pooled_shapiro$statistic),
    shapiro_p = if (is.null(pooled_shapiro)) NA_real_ else pooled_shapiro$p.value,
    n = length(atpase_analysis_data$within_cell_residual),
    caveat = N3_CAVEAT,
    stringsAsFactors = FALSE
))
for (current_timepoint in TIMEPOINTS) {
    tp_residuals <- atpase_analysis_data$within_cell_residual[
        atpase_analysis_data$timepoint == current_timepoint
    ]
    tp_shapiro <- tryCatch(shapiro.test(tp_residuals), error = function(e) NULL)
    shapiro_rows <- rbind(shapiro_rows, data.frame(
        scope = paste0("within_cell_residuals_timepoint_", current_timepoint),
        shapiro_W = if (is.null(tp_shapiro)) NA_real_ else unname(tp_shapiro$statistic),
        shapiro_p = if (is.null(tp_shapiro)) NA_real_ else tp_shapiro$p.value,
        n = length(tp_residuals),
        caveat = N3_CAVEAT,
        stringsAsFactors = FALSE
    ))
}
shapiro_csv_path <- file.path(OUTPUT_DIRECTORY, "atpase_diagnostics_shapiro.csv")
if (!file.exists(shapiro_csv_path) || OVERWRITE_CSVS) {
    write.csv(shapiro_rows, shapiro_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(shapiro_csv_path))
}
print(shapiro_rows[, c("scope", "shapiro_W", "shapiro_p", "n")])

# ------------------------------------------------------------------------------
# C16: PER-TIMEPOINT WT-CONTRASTS (PRIMARY/plotted; the "not restored to WT"
# half of the delta-free claim -- rejectable, load-bearing).
# Block = experiment_label. WT is SUBSET to each suppressor's SHARED blocks
# (D1): a "WT vs OrcX" paired test uses ONLY the experiments where BOTH appear.
# WT's full replicate count does NOT strengthen these contrasts.
# Block-mean aggregation collapses the two Orc4 lanes in 2020_08_24, so the
# block is the paired unit. Holm WITHIN each timepoint (I2), NEVER across
# timepoints. Wilcoxon companions DECORATIVE (I4).
# ------------------------------------------------------------------------------
message("--- C16 per-timepoint WT-contrasts (WT subset to shared blocks) ---")

block_means <- atpase_analysis_data %>%
    group_by(sample, experiment_label, timepoint) %>%
    summarise(block_mean_corrected = mean(percent_adp_corrected), .groups = "drop") %>%
    as.data.frame()

CONTRAST_SAMPLES <- c(
    "4R", "Orc1_E495K_4R", "Orc3_P481L_4R",
    "Orc4_P225S_4R", "Orc5_E104K_4R", "Orc6_E304K_4R"
)

per_timepoint_contrast_rows <- data.frame()
for (current_timepoint in TIMEPOINTS) {
    timepoint_raw_p <- rep(NA_real_, length(CONTRAST_SAMPLES))
    timepoint_rows_list <- vector("list", length(CONTRAST_SAMPLES))

    for (k in seq_along(CONTRAST_SAMPLES)) {
        current_group <- CONTRAST_SAMPLES[k]
        wt_rows <- block_means[
            block_means$sample == "WT" & block_means$timepoint == current_timepoint,
        ]
        grp_rows <- block_means[
            block_means$sample == current_group & block_means$timepoint == current_timepoint,
        ]
        shared_experiments <- sort(intersect(
            wt_rows$experiment_label, grp_rows$experiment_label
        ))
        n_pairs <- length(shared_experiments)

        wt_vec <- wt_rows$block_mean_corrected[
            match(shared_experiments, wt_rows$experiment_label)
        ]
        grp_vec <- grp_rows$block_mean_corrected[
            match(shared_experiments, grp_rows$experiment_label)
        ]
        paired_diff <- wt_vec - grp_vec

        estimate_wt_minus_group <- if (n_pairs >= 1) mean(paired_diff) else NA_real_
        is_computable <- FALSE
        ci_lower <- NA_real_; ci_upper <- NA_real_
        t_statistic <- NA_real_; degrees_freedom <- NA_real_
        raw_p <- NA_real_; wilcoxon_p <- NA_real_; contrast_note <- ""

        if (n_pairs < 2) {
            contrast_note <- paste0(
                "only ", n_pairs, " shared block(s); paired t-test needs >= 2"
            )
        } else if (sd(paired_diff) == 0) {
            contrast_note <- paste0(
                "zero variance in paired differences (all values identical at ",
                "this timepoint, e.g. both groups at zero); paired t-test undefined"
            )
        } else {
            tt <- t.test(wt_vec, grp_vec, paired = TRUE)
            is_computable <- TRUE
            estimate_wt_minus_group <- unname(tt$estimate)
            ci_lower <- tt$conf.int[1]
            ci_upper <- tt$conf.int[2]
            t_statistic <- unname(tt$statistic)
            degrees_freedom <- unname(tt$parameter)
            raw_p <- tt$p.value
            wilcoxon_p <- tryCatch(
                suppressWarnings(wilcox.test(wt_vec, grp_vec, paired = TRUE)$p.value),
                error = function(e) NA_real_
            )
            contrast_note <- "ok"
        }

        timepoint_raw_p[k] <- if (is_computable) raw_p else NA_real_
        timepoint_rows_list[[k]] <- data.frame(
            timepoint = current_timepoint,
            comparison = paste0("WT_vs_", current_group),
            group_label = current_group,
            claim = "not_restored_to_WT (rejectable; load-bearing half of delta-free framing)",
            test = "paired t-test; block=experiment_label; WT subset to shared blocks; block-mean aggregated",
            holm_family = paste0("WT-contrasts within timepoint=", current_timepoint, " (NEVER across timepoints, I2)"),
            n_pairs = n_pairs,
            shared_experiments = paste(shared_experiments, collapse = ";"),
            estimate_wt_minus_group = estimate_wt_minus_group,
            ci_lower = ci_lower,
            ci_upper = ci_upper,
            t_statistic = t_statistic,
            df = degrees_freedom,
            raw_p_value = raw_p,
            holm_adjusted_p_value = NA_real_,
            holm_family_size = NA_integer_,
            wilcoxon_p_value_decorative = wilcoxon_p,
            computable = is_computable,
            note = contrast_note,
            stringsAsFactors = FALSE
        )
    }

    computable_idx <- which(!is.na(timepoint_raw_p))
    holm_p <- rep(NA_real_, length(CONTRAST_SAMPLES))
    if (length(computable_idx) > 0) {
        holm_p[computable_idx] <- p.adjust(
            timepoint_raw_p[computable_idx], method = "holm"
        )
    }
    for (k in seq_along(CONTRAST_SAMPLES)) {
        timepoint_rows_list[[k]]$holm_adjusted_p_value <- holm_p[k]
        timepoint_rows_list[[k]]$holm_family_size <- length(computable_idx)
    }
    per_timepoint_contrast_rows <- rbind(
        per_timepoint_contrast_rows, do.call(rbind, timepoint_rows_list)
    )
}

message("Per-timepoint WT-contrast pairing and p-values:")
print(per_timepoint_contrast_rows[, c(
    "timepoint", "comparison", "n_pairs", "computable",
    "estimate_wt_minus_group", "raw_p_value", "holm_adjusted_p_value", "note"
)])

per_timepoint_csv_path <- file.path(
    OUTPUT_DIRECTORY, "atpase_per_timepoint_wt_contrasts.csv"
)
if (!file.exists(per_timepoint_csv_path) || OVERWRITE_CSVS) {
    write.csv(per_timepoint_contrast_rows, per_timepoint_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(per_timepoint_csv_path))
}

# ------------------------------------------------------------------------------
# C17: POOLED MIXED MODELS (nlme). PRIMARY home for the 4R-floor-proximity
# claim (D2 -- NOT a stand-in diagnostic).
# Two complementary fits, each answering one question (procedural, no helpers):
#  (1) ADDITIVE model, reference = 4R -> gives the suppressor-vs-4R difference
#      POOLED across timepoints + 95% CI directly from intervals(). This is the
#      reported 4R-floor-proximity estimate. (A "pooled across timepoints"
#      contrast IS the additive sample effect under no-interaction.)
#  (2) INTERACTION model sample*timepoint -> tests whether the curves diverge
#      (sample:timepoint) and the sample main effect, via anova(). The
#      interaction test is the check that the additive pooling is reasonable
#      for the near-floor groups (expected: WT drives any interaction; the
#      suppressor curves stay flat near 4R).
# timepoint as FACTOR (accumulation saturates -> not linear); the numeric/linear
# trend alternative is documented and NOT adopted (only 4 levels; saturating
# kinetics; factor is the conservative choice).
# Caveat (C15): WT-vs-suppressor variance heterogeneity violates the
# homoscedastic-residual assumption; CIs for the near-floor contrasts are if
# anything CONSERVATIVELY WIDE (residual variance is inflated by WT spread).
# Both fits wrapped in tryCatch -> non-convergence is a REPORTED finding.
# ------------------------------------------------------------------------------
message("--- C17 pooled mixed models (nlme) ---")

atpase_analysis_data$sample_factor <- factor(
    atpase_analysis_data$sample, levels = SAMPLE_DISPLAY_ORDER
)
atpase_analysis_data$sample_ref4R <- relevel(
    factor(atpase_analysis_data$sample, levels = SAMPLE_DISPLAY_ORDER), ref = "4R"
)
atpase_analysis_data$timepoint_factor <- factor(
    atpase_analysis_data$timepoint, levels = TIMEPOINTS
)
atpase_analysis_data$experiment_block <- factor(atpase_analysis_data$experiment_label)

# (1) additive pooled model, reference = 4R
pooled_additive_fit <- tryCatch(
    nlme::lme(
        fixed = percent_adp_corrected ~ sample_ref4R + timepoint_factor,
        random = ~ 1 | experiment_block,
        data = atpase_analysis_data,
        method = "REML"
    ),
    error = function(e) e
)

pooled_contrast_rows <- data.frame()
if (inherits(pooled_additive_fit, "lme")) {
    pooled_ttable <- summary(pooled_additive_fit)$tTable
    pooled_intervals <- tryCatch(
        intervals(pooled_additive_fit, which = "fixed")$fixed,
        error = function(e) NULL
    )
    coef_names <- rownames(pooled_ttable)
    sample_coef_names <- coef_names[grepl("^sample_ref4R", coef_names)]
    for (cn in sample_coef_names) {
        lvl <- sub("^sample_ref4R", "", cn)
        ci_l <- if (!is.null(pooled_intervals)) pooled_intervals[cn, "lower"] else NA_real_
        ci_u <- if (!is.null(pooled_intervals)) pooled_intervals[cn, "upper"] else NA_real_
        pooled_contrast_rows <- rbind(pooled_contrast_rows, data.frame(
            contrast = paste0(lvl, "_minus_4R"),
            level = lvl,
            estimate_minus_4R = unname(pooled_ttable[cn, "Value"]),
            ci_lower = ci_l,
            ci_upper = ci_u,
            std_error = unname(pooled_ttable[cn, "Std.Error"]),
            df = unname(pooled_ttable[cn, "DF"]),
            t_value = unname(pooled_ttable[cn, "t-value"]),
            p_value = unname(pooled_ttable[cn, "p-value"]),
            model = "lme percent_adp_corrected ~ sample + timepoint, ~1|experiment (REML); ref=4R; POOLED across timepoints",
            interpretation = if (lvl == "WT") {
                "WT elevation above 4R (pooled); NOT the floor-proximity claim"
            } else {
                "4R-floor-proximity (suppressor-4R pooled); near 0 with CI bracketing 0 supports 'near the 4R floor'"
            },
            converged = TRUE,
            stringsAsFactors = FALSE
        ))
    }
    message("Pooled additive lme CONVERGED. Suppressor-vs-4R pooled contrasts:")
    print(pooled_contrast_rows[, c(
        "contrast", "estimate_minus_4R", "ci_lower", "ci_upper", "p_value"
    )])
} else {
    message("Pooled additive lme DID NOT CONVERGE: ",
        conditionMessage(pooled_additive_fit))
    pooled_contrast_rows <- data.frame(
        contrast = "CONVERGENCE_FAILURE", level = NA_character_,
        estimate_minus_4R = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        std_error = NA_real_, df = NA_real_, t_value = NA_real_, p_value = NA_real_,
        model = "lme percent_adp_corrected ~ sample + timepoint, ~1|experiment (REML); ref=4R",
        interpretation = conditionMessage(pooled_additive_fit),
        converged = FALSE, stringsAsFactors = FALSE
    )
}

pooled_contrasts_csv_path <- file.path(
    OUTPUT_DIRECTORY, "atpase_pooled_contrasts_vs_4R.csv"
)
if (!file.exists(pooled_contrasts_csv_path) || OVERWRITE_CSVS) {
    write.csv(pooled_contrast_rows, pooled_contrasts_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(pooled_contrasts_csv_path))
}

# (2) interaction model: divergence + sample main effect
interaction_fit <- tryCatch(
    nlme::lme(
        fixed = percent_adp_corrected ~ sample_factor * timepoint_factor,
        random = ~ 1 | experiment_block,
        data = atpase_analysis_data,
        method = "REML"
    ),
    error = function(e) e
)

interaction_anova_rows <- data.frame()
if (inherits(interaction_fit, "lme")) {
    interaction_anova <- anova(interaction_fit)
    interaction_anova_rows <- data.frame(
        term = rownames(interaction_anova),
        numDF = interaction_anova[["numDF"]],
        denDF = interaction_anova[["denDF"]],
        F_value = interaction_anova[["F-value"]],
        p_value = interaction_anova[["p-value"]],
        model = "lme percent_adp_corrected ~ sample*timepoint, ~1|experiment (REML); sequential anova",
        converged = TRUE,
        stringsAsFactors = FALSE
    )
    message("Interaction lme CONVERGED. Sequential anova (divergence = sample:timepoint):")
    print(interaction_anova_rows)
} else {
    message("Interaction lme DID NOT CONVERGE: ",
        conditionMessage(interaction_fit))
    interaction_anova_rows <- data.frame(
        term = "CONVERGENCE_FAILURE", numDF = NA_real_, denDF = NA_real_,
        F_value = NA_real_, p_value = NA_real_,
        model = conditionMessage(interaction_fit), converged = FALSE,
        stringsAsFactors = FALSE
    )
}

interaction_anova_csv_path <- file.path(
    OUTPUT_DIRECTORY, "atpase_interaction_anova.csv"
)
if (!file.exists(interaction_anova_csv_path) || OVERWRITE_CSVS) {
    write.csv(interaction_anova_rows, interaction_anova_csv_path, row.names = FALSE)
    message("Saved CSV: ", basename(interaction_anova_csv_path))
}

message("=== STATISTICAL ANALYSIS COMPLETE ===")

# ==============================================================================
# PLOTTING (legacy in-script timecourse; the ANNOTATED stats figure is produced
# by 260411_orc4r-screen_plot-sofa-atpase.R, which reads the C18a CSV)
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
