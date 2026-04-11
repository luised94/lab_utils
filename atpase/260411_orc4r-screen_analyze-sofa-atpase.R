# ==============================================================================
# ATPase Percent Hydrolysis Analysis - Consolidated
# ==============================================================================
# Date created: 2026-04-11
# Consolidates analysis of ORC ATPase assays across Orc1, Orc3, Orc4, Orc5,
# and Orc6 suppressor mutants in the 4R background.
# Input: Excel files containing ImageJ band intensity measurements.
# Output: Combined raw data CSV, processed data CSV, summary CSV, and
#         percent hydrolysis timecourse plot (PDF).
# Usage: Source in R REPL or run via Rscript.
# ==============================================================================

library(xlsx)
library(tidyverse)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
message("=== CONFIGURATION ===")

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

TIMEPOINTS <- c(0, 15, 45, 90)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")
if (nchar(MC_DROPBOX_PATH) == 0) {
    stop(
        "MC_DROPBOX_PATH not defined in environment.\n",
        "Set via Sys.setenv(MC_DROPBOX_PATH = '/your/path') or export in shell."
    )
}

BASE_EXPERIMENT_DIRECTORY <- file.path(
    MC_DROPBOX_PATH,
    "Lab", "Experiments",
    "ATPase"
)

OUTPUT_DIRECTORY <- file.path(MC_DROPBOX_PATH, "Lab", "Experiments",
    "ATPase", "2020_09_03 ATPase Analysis of 4R supps",
    "consolidated_analysis"
)

if (!dir.exists(OUTPUT_DIRECTORY)) {
    dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)
    message("Created output directory: ", OUTPUT_DIRECTORY)
}

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
# Layout A: 16 raw rows. No_ORC occupies rows 1-4 (all timepoints).
#           Odd ImageJ index = ADP, next row = ATP.
#           Background = No_ORC at t=90 (processed row 4).
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
# Validate file paths
# ------------------------------------------------------------------------------
message("Validating file paths...")

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
message("Output directory: ", OUTPUT_DIRECTORY)
message("=== CONFIGURATION COMPLETE ===")
