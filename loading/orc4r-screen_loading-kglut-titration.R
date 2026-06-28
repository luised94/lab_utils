# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., readxl::read_excel).
# Date updated: 2026-04-22
# Data of the xlsx was produced by analyzing tiff files using imagej.
# Usage: source("orc4r-screen_loading-kglut-titration.R")
# The output is the plot called faceted_by_kglut_plot.
# All other plots kept for reference.
# Prerequisites:
#   1. Either MC_DROPBOX_PATH is set (original data home) OR the input Excel
#      file sits in the same directory as this script (Zenodo co-located).
#   2. Font setup (one-time per machine):
#      a. Install system fonts (WSL/Ubuntu/Debian):
#           sudo apt update
#           sudo apt install ttf-mscorefonts-installer
#           sudo fc-cache -fv
#         Verify: fc-list | grep -i arial
#      b. Import fonts into R extrafont database:
#           renv::install(c("extrafont", "ragg"))
#           library(extrafont); font_import(prompt = FALSE)
#         Verify: any(grepl("Arial", extrafont::fonts()))
#      c. Per-session (handled automatically by runtime check below):
#           extrafont::loadfonts(device = "pdf")
#   3. renv lockfile: renv_loading-analysis.lock
#      Restore: renv::restore(lockfile = "renv_loading-analysis.lock")
#   4. R must be compiled with Cairo support: capabilities("cairo") == TRUE
#      If FALSE, install libcairo2-dev (Debian/Ubuntu) and rebuild R.
#   5. Input: single Excel file, sheets 2/3/4 (sheet 1 = metadata).
#      14 rows per sheet:
#        Row 1:    input lane (filtered by input == "no")
#        Row 2:    no-ORC negative control (filtered by orc != "None")
#                  Also used as background for intensity subtraction.
#        Rows 3-14: 4 conditions x 3 kglut concentrations
#      Three sheets = 3 biological replicates.
#      WT/ORC4R/+4sofa: n=3 across all sheets.
#      Rotating suppressor (+1/+3/+5sofa): n=1 per sheet.
#      14 rows -> filter input=="no" -> 13 -> filter orc!="None" -> 12 per sheet.

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
# Configuration
# ==============================================================================

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
        "This script must be run via source(\"orc4r-screen_loading-kglut-titration.R\").\n",
        "Rscript and interactive invocation are not supported ",
        "(no script path is available to resolve data locations)."
    )
}
SCRIPT_DIRECTORY <- dirname(normalizePath(script_path_under_source))

# MC_DROPBOX_PATH is the original (Dropbox) data home; may be unset in a
# Zenodo deposit where the input Excel sits alongside this script.
MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")

library(readxl)
library(tidyverse)

# Runtime checks: fail fast with actionable diagnostics.
stopifnot(
    "Cairo graphics not available. Install libcairo2-dev and rebuild R." =
        capabilities("cairo")
)
extrafont::loadfonts(device = "pdf", quiet = TRUE)

if (!("Arial" %in% extrafont::fonts())) {
    stop(
        "Arial font not found in extrafont database.\n",
        "See Prerequisites section at top of this script for full setup:\n",
        "  1. Install system fonts: sudo apt install ttf-mscorefonts-installer\n",
        "  2. Rebuild font cache: sudo fc-cache -fv\n",
        "  3. Import into R: library(extrafont); font_import(prompt = FALSE)\n",
        "  4. Restart R and re-source this script."
    )
}

message("Font and Cairo checks passed.")

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

INPUT_FILENAME <- "Analysis.xlsx"
OUTPUT_DIRECTORY <- "~/data/loading_analysis"
EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication"

# Define CSV file paths
summary_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_kglut-titration_per-kglut-summary.csv")
full_data_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_kglut-titration_full-data.csv")

# Path resolution (C1): script-relative (Zenodo co-located) -> MC_DROPBOX_PATH
# -> stop() naming both attempted absolute paths.
script_relative_input_filepath <- file.path(SCRIPT_DIRECTORY, INPUT_FILENAME)
if (nchar(MC_DROPBOX_PATH) > 0) {
    dropbox_input_filepath <- file.path(MC_DROPBOX_PATH, EXPERIMENT_DIRECTORY, INPUT_FILENAME)
} else {
    dropbox_input_filepath <- NA_character_
}
if (file.exists(script_relative_input_filepath)) {
    INPUT_FILEPATH <- script_relative_input_filepath
} else if (!is.na(dropbox_input_filepath) && file.exists(dropbox_input_filepath)) {
    INPUT_FILEPATH <- dropbox_input_filepath
} else {
    stop(
        "Input file not found in either supported location:\n",
        "  script-relative: ", script_relative_input_filepath, "\n",
        "  MC_DROPBOX_PATH:  ",
        if (is.na(dropbox_input_filepath)) "<MC_DROPBOX_PATH not set>" else dropbox_input_filepath
    )
}

# ==============================================================================
# File Validation
# ==============================================================================

if (!dir.exists(OUTPUT_DIRECTORY)) {
  dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)

}

if (!file.exists(INPUT_FILEPATH)) {
  stop("INPUT_FILEPATH does not exist: ", INPUT_FILEPATH)

}

sheet_indices <- c(2, 3, 4)
EXPECTED_NUMBER_OF_ROWS <- 14
EXPECTED_NUMBER_OF_COLUMNS <- 7

# Row index of the no-ORC negative control lane in each sheet.
# Used for background intensity subtraction. This assumption breaks if
# rows are reordered in the Excel source.
# TODO: long-term fix - add a named marker column (e.g., "background") to
# the Excel data and use programmatic lookup instead of row index.
BACKGROUND_ROW_INDEX <- 2

COLUMN_TO_REMOVE <- "...1" # Column has row numbers from imagej export/copy-paste

# @NOTE: readxl trims trailing whitespace from column names by default.
# If a column name has a trailing space in the Excel cell (e.g., "ORC "),
# readxl reads it as "ORC" (without the dot suffix that xlsx produced).
# Unnamed columns are read as "...1", "...2", etc.
REQUIRED_COLUMNS <- c(
  "...1", "Intensity", "Lane",
  "ORC", "Suppressor",
  "kGlut", "Input"
)

# Define custom orderings
factor_order <- list(
  "suppressor" = c("None", "1EK", "3PL", "4PS", "5EK"),
  "kglut" = c("250", "300", "350"),
  "label" = c("WT", "ORC4R", "+1sofa", "+3sofa", "+4sofa", "+5sofa")
)

# Centralized plot configuration. Consumed by ggplot calls in the plot section.
# Any parameter set to NULL means "use ggplot default."
# fill_colors covers all conditions across both scripts for cross-figure
# consistency. +1sofa and +4sofa colors are swapped relative to default Set1
# so +4sofa keeps its color in the companion single-panel figure.
PLOT_CONFIG <- list(
    fill_colors = c(
        "WT" = "#E41A1C", "ORC4R" = "#377EB8",
        "+1sofa" = "#FF7F00", "+3sofa" = "#984EA3",
        "+4sofa" = "#4DAF4A", "+5sofa" = "#FFFF33",
        "+6sofa" = "#A65628"
    ),
    bar = list(width = 0.7, color = "black", linewidth = 0.4),
    errorbar = list(width = 0.25, linewidth = 0.6),
    point = list(
        size = 2, fill = "grey30", color = "black", stroke = 0.5,
        jitter_width = 0.15, jitter_seed = 42
    ),
    replicate_shapes = c("1" = 21, "2" = 24, "3" = 22),
    theme = list(base_family = "Arial", base_size = 12, legend_position = "bottom"),
    output = list(
        device = cairo_pdf,
        width = 7.5,
        height = 4.5
    )
)


message("Paths for input and output set...")

# ==============================================================================
# Data Loading
# ==============================================================================
df_lst <- vector(mode = "list", length = length(sheet_indices))
temp_df <- data.frame(matrix(
  NA,
  nrow = EXPECTED_NUMBER_OF_ROWS,
  ncol = EXPECTED_NUMBER_OF_COLUMNS
))
loading_data <- data.frame(matrix(
  NA,
  nrow = EXPECTED_NUMBER_OF_ROWS * length(sheet_indices),
  ncol = EXPECTED_NUMBER_OF_COLUMNS
))

message("Reading in loading assay sheets...")

df_count <- 0
for (sheet_idx in sheet_indices){
  df_count <- df_count + 1

  temp_df <- readxl::read_excel(INPUT_FILEPATH, sheet = sheet_idx)

  stopifnot(
    "Number of rows in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      nrow(temp_df) == EXPECTED_NUMBER_OF_ROWS,
    "Number of columns in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      ncol(temp_df) == EXPECTED_NUMBER_OF_COLUMNS,
    "df_lst[[df_count]] colnames not identical to REQUIRED_COLUMNS." =
      identical(REQUIRED_COLUMNS, colnames(temp_df))
  )



  temp_df <- dplyr::rename(
      temp_df,
      intensity = Intensity,
      lane = Lane,
      orc = ORC,
      suppressor = Suppressor,
      kglut = kGlut,
      input = Input
  )

  temp_df <- temp_df %>%
    mutate(
      # Subtract background (Row 2) from all rows
      net_intensity = intensity - intensity[BACKGROUND_ROW_INDEX],
      replicate = df_count
    ) %>%
    # Correct for input volume (0.5 factor) using the net_intensity
    mutate(relative_to_input = (net_intensity / net_intensity[input == "yes"]) * 0.5) %>%
    filter(input == "no")
  df_lst[[df_count]] <- temp_df

} # end read-data for loop

loading_data <- do.call(rbind, df_lst)
loading_data <- loading_data[, names(loading_data) != COLUMN_TO_REMOVE]

message("loading_data preparation complete...")

# ==============================================================================
# Preprocessing
# ==============================================================================

# Remove negative control rows (orc == "None") before label assignment.
# These rows have no meaningful label and would cause indexing errors in
# downstream paired calculations that subset by label.
loading_data <- loading_data %>%
    filter(orc != "None")
message("Negative control rows removed (orc == 'None').")


loading_data <- loading_data %>%
  mutate(label = case_when(
    orc == "WT" & suppressor == "None" ~ "WT",
    orc == "RA" & suppressor == "None" ~ "ORC4R",
    orc == "RA" & suppressor == "4PS" ~ "+4sofa",
    orc == "RA" & suppressor == "1EK" ~ "+1sofa",
    orc == "RA" & suppressor == "3PL" ~ "+3sofa",
    orc == "RA" & suppressor == "5EK" ~ "+5sofa"
  ))
stopifnot(
    "NA labels found after case_when mapping. Check for unmatched orc/suppressor combinations." =
        sum(is.na(loading_data$label)) == 0
)

message("Adjust loading_data to factor order...")
for (col_name in names(factor_order)) {

  loading_data[[col_name]] <- factor(
    loading_data[[col_name]],
    levels = factor_order[[col_name]],
    ordered = FALSE
  )
}

# Row-count assertion: each label appears exactly once per replicate x kglut.
label_counts <- loading_data %>% count(replicate, kglut, label)
bad_labels <- label_counts %>% filter(n != 1)
if (nrow(bad_labels) > 0) {
    stop("Expected exactly 1 row per replicate x kglut x label. Offending groups:\n",
         paste(capture.output(print(as.data.frame(bad_labels))), collapse = "\n"))
}

# WT uniqueness: exactly one WT per replicate x kglut.
wt_counts <- loading_data %>% filter(label == "WT") %>% count(replicate, kglut)
bad_wt <- wt_counts %>% filter(n != 1)
if (nrow(bad_wt) > 0) {
    stop("Expected exactly 1 WT per replicate x kglut. Offending groups:\n",
         paste(capture.output(print(as.data.frame(bad_wt))), collapse = "\n"))
}

# ORC4R uniqueness: exactly one ORC4R per replicate x kglut.
orc4r_counts <- loading_data %>% filter(label == "ORC4R") %>% count(replicate, kglut)
bad_orc4r <- orc4r_counts %>% filter(n != 1)
if (nrow(bad_orc4r) > 0) {
    stop("Expected exactly 1 ORC4R per replicate x kglut. Offending groups:\n",
         paste(capture.output(print(as.data.frame(bad_orc4r))), collapse = "\n"))
}
message("Row-count and uniqueness assertions passed.")

loading_data <- loading_data %>%
  # Step 1: Normalize to the WT of the SAME salt concentration
  group_by(replicate, kglut) %>%
  mutate(
    normalized_to_wildtype_per_kglut = relative_to_input / relative_to_input[orc == "WT"][1]
  ) %>%
  # Step 2: Normalize to the WT of the 250 mM concentration (Global Baseline)
  group_by(replicate) %>%
  mutate(
    normalized_to_wildtype_250mm = relative_to_input / relative_to_input[orc == "WT" & kglut == "250"][1]
  ) %>%
  ungroup()

# Convert per-kglut normalization to percentage scale (WT = 100 within each kglut).
# This is the plotting value and the basis for all derived calculations.
loading_data <- loading_data %>%
    mutate(percent_wildtype = normalized_to_wildtype_per_kglut * 100)

message("Percent wildtype column computed.")

# Global baseline: all values as percentage of WT at 250 mM (within same replicate).
# Unlike per-kglut normalization, WT is NOT 100% at 300/350 mM - this shows
# how WT loading itself declines with increasing salt.
loading_data <- loading_data %>%
    mutate(percent_wildtype_global = normalized_to_wildtype_250mm * 100)
message("Percent wildtype global baseline column computed.")

# ==============================================================================
# Negative-value detection and floor (kglut)  [Thread 2b, pre-C12]
# ==============================================================================
# User decision: negative percent_wildtype values come from background
# subtraction in low-signal lanes and are QUANTIFICATION ARTIFACTS. Handling is
# FLOOR-AFTER-DETECTION (not silent floor, not exclusion): print every offending
# row so the user can confirm they are small artifacts, assert they are FEW and
# SMALL (a systematic failure must STOP the script loudly), then floor into NEW
# columns. Raw columns are retained for transparency and the methods doc.
# Flooring can only move a value UP toward zero, so it cannot manufacture a
# restorer effect or exaggerate a WT-mutant gap -> this is the safety argument.

# Membership-checked guard (carry-forward #1): NULL columns would make the
# value scan vacuously pass, so confirm the columns exist first.
stopifnot(
    "Negative-detection: percent_wildtype column missing." =
        "percent_wildtype" %in% names(loading_data),
    "Negative-detection: percent_wildtype_global column missing." =
        "percent_wildtype_global" %in% names(loading_data),
    "Negative-detection: net_intensity column missing." =
        "net_intensity" %in% names(loading_data)
)

negative_value_rows <- loading_data %>%
    filter(net_intensity < 0 | percent_wildtype < 0 | percent_wildtype_global < 0) %>%
    select(label, replicate, kglut, net_intensity, relative_to_input,
           percent_wildtype, percent_wildtype_global)

negative_value_count <- nrow(negative_value_rows)
message("Negative-value detection: found ", negative_value_count,
        " row(s) with net_intensity < 0 or a negative percent value.")
if (negative_value_count > 0) {
    message("Full list of negative (floored) rows:")
    print(as.data.frame(negative_value_rows))
}

# Documented envelope thresholds. A breach means a systematic quantification
# failure, which MUST stop the script rather than being floored away.
NEGATIVE_VALUE_MAX_COUNT <- 3          # stray artifacts only (of 36 rows)
NEGATIVE_VALUE_MAX_ABS_PERCENT <- 50   # small relative to WT = 100

stopifnot(
    "Negative envelope: too many negative rows -> systematic failure, not a stray artifact. STOP." =
        negative_value_count <= NEGATIVE_VALUE_MAX_COUNT,
    "Negative envelope: a negative value is too large in magnitude -> not a low-signal artifact. STOP." =
        (negative_value_count == 0 ||
         max(abs(negative_value_rows$percent_wildtype)) < NEGATIVE_VALUE_MAX_ABS_PERCENT)
)
message("Negative envelope assertions passed (count <= ", NEGATIVE_VALUE_MAX_COUNT,
        ", |percent_wildtype| < ", NEGATIVE_VALUE_MAX_ABS_PERCENT, ").")

# Floor into NEW columns. Raw columns retained. Decision (carry-forward #1d):
# the global-baseline column used by C13 is ALSO derived from floored values so
# a single artifact cannot distort the salt-sensitivity baseline. For the test
# genotypes (WT/ORC4R/+4sofa) NO value is negative, so floored == raw for every
# value entering C12/C13; the only floored row is +5sofa (rotating, filtered).
loading_data <- loading_data %>%
    mutate(
        percent_wildtype_floored = pmax(percent_wildtype, 0),
        percent_wildtype_global_floored = pmax(percent_wildtype_global, 0)
    )
message("Floored columns computed (percent_wildtype_floored, percent_wildtype_global_floored).")

# Persist the floored-row list (empty CSV if none) for the methods record.
floored_rows_csv_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_floored-rows.csv")
if (!file.exists(floored_rows_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(as.data.frame(negative_value_rows),
              floored_rows_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(floored_rows_csv_output_filepath))
}

# Derived calculations: paired within replicate and kglut so each condition
# is compared to the WT and ORC4R values from the same experiment and salt
# concentration. Formulas match Script 1 for cross-script comparability.
loading_data <- loading_data %>%
    group_by(replicate, kglut) %>%
    mutate(
        # Percent difference: |A - B| / ((A + B) / 2) * 100
        percent_difference_from_wildtype = abs(percent_wildtype - percent_wildtype[label == "WT"]) /
            ((percent_wildtype + percent_wildtype[label == "WT"]) / 2) * 100,
        percent_difference_from_orc4r = abs(percent_wildtype - percent_wildtype[label == "ORC4R"]) /
            ((percent_wildtype + percent_wildtype[label == "ORC4R"]) / 2) * 100,
        # Percent change: (A - reference) / reference * 100
        # Directional: negative means condition loads less than reference.
        percent_change_from_wildtype = (percent_wildtype - percent_wildtype[label == "WT"]) /
            percent_wildtype[label == "WT"] * 100,
        percent_change_from_orc4r = (percent_wildtype - percent_wildtype[label == "ORC4R"]) /
            percent_wildtype[label == "ORC4R"] * 100,
        # Fold change: condition / ORC4R. Values > 1 indicate rescue.
        fold_change_from_orc4r = percent_wildtype / percent_wildtype[label == "ORC4R"]
    ) %>%
    ungroup()

message("Percent difference, percent change, and fold change columns computed.")

# Guard against division by zero in derived calculations.
derived_cols <- c("fold_change_from_orc4r", "percent_change_from_wildtype",
                  "percent_change_from_orc4r")
for (col in derived_cols) {
    stopifnot(
        !any(is.infinite(loading_data[[col]])),
        !any(is.nan(loading_data[[col]]))
    )
}
message("Inf/NaN guard passed on derived calculations.")

# Derived calculations on global baseline: same formulas as per-kglut versions
# but operating on percent_wildtype_global. Paired within replicate and kglut.
loading_data <- loading_data %>%
    group_by(replicate, kglut) %>%
    mutate(
        # Percent difference: |A - B| / ((A + B) / 2) * 100
        percent_difference_from_wildtype_global = abs(percent_wildtype_global - percent_wildtype_global[label == "WT"]) /
            ((percent_wildtype_global + percent_wildtype_global[label == "WT"]) / 2) * 100,
        percent_difference_from_orc4r_global = abs(percent_wildtype_global - percent_wildtype_global[label == "ORC4R"]) /
            ((percent_wildtype_global + percent_wildtype_global[label == "ORC4R"]) / 2) * 100,
        # Percent change: (A - reference) / reference * 100
        # Directional: negative means condition loads less than reference.
        percent_change_from_wildtype_global = (percent_wildtype_global - percent_wildtype_global[label == "WT"]) /
            percent_wildtype_global[label == "WT"] * 100,
        percent_change_from_orc4r_global = (percent_wildtype_global - percent_wildtype_global[label == "ORC4R"]) /
            percent_wildtype_global[label == "ORC4R"] * 100,
        # Fold change: condition / ORC4R. Values > 1 indicate rescue.
        fold_change_from_orc4r_global = percent_wildtype_global / percent_wildtype_global[label == "ORC4R"]
    ) %>%
    ungroup()
message("Global baseline derived calculations computed.")

# Guard against division by zero in global baseline derived calculations.
derived_cols_global <- c("fold_change_from_orc4r_global",
                         "percent_change_from_wildtype_global",
                         "percent_change_from_orc4r_global")
for (col in derived_cols_global) {
    stopifnot(
        !any(is.infinite(loading_data[[col]])),
        !any(is.nan(loading_data[[col]]))
    )
}
message("Inf/NaN guard passed on global baseline derived calculations.")


# ==============================================================================
# Summary Statistics
# ==============================================================================

summary_loading_data <- loading_data %>%
  filter(!(label %in% c("+1sofa", "+3sofa", "+5sofa"))) %>%
  droplevels() %>%
  group_by(kglut, label) %>%
  summarise(
    mean_percent_wildtype = mean(percent_wildtype, na.rm = TRUE),
    sd_percent_wildtype = sd(percent_wildtype, na.rm = TRUE),
    mean_percent_difference_from_wildtype = mean(percent_difference_from_wildtype, na.rm = TRUE),
    sd_percent_difference_from_wildtype = sd(percent_difference_from_wildtype, na.rm = TRUE),
    mean_percent_difference_from_orc4r = mean(percent_difference_from_orc4r, na.rm = TRUE),
    sd_percent_difference_from_orc4r = sd(percent_difference_from_orc4r, na.rm = TRUE),
    mean_percent_change_from_wildtype = mean(percent_change_from_wildtype, na.rm = TRUE),
    sd_percent_change_from_wildtype = sd(percent_change_from_wildtype, na.rm = TRUE),
    mean_percent_change_from_orc4r = mean(percent_change_from_orc4r, na.rm = TRUE),
    sd_percent_change_from_orc4r = sd(percent_change_from_orc4r, na.rm = TRUE),
    mean_fold_change_from_orc4r = mean(fold_change_from_orc4r, na.rm = TRUE),
    sd_fold_change_from_orc4r = sd(fold_change_from_orc4r, na.rm = TRUE),
    # Retain raw and global baseline summaries for reference
    mean_relative_to_input = mean(relative_to_input, na.rm = TRUE),
    sd_relative_to_input = sd(relative_to_input, na.rm = TRUE),
    mean_normalized_to_wildtype_250mm = mean(normalized_to_wildtype_250mm, na.rm = TRUE),
    sd_normalized_to_wildtype_250mm = sd(normalized_to_wildtype_250mm, na.rm = TRUE),
    replicate_count = n(),
    .groups = "drop"
  )

message("Summary statistics computed.")

stopifnot(
    "Not all summary groups have 3 replicates." =
        all(summary_loading_data$replicate_count == 3)
)

# Global baseline summary: same pipeline as summary_loading_data but
# summarizing percent_wildtype_global and its derived calculations.
# Normalization reference is WT at 250 mM across all panels.
summary_loading_data_global <- loading_data %>%
    filter(!(label %in% c("+1sofa", "+3sofa", "+5sofa"))) %>%
    droplevels() %>%
    group_by(kglut, label) %>%
    summarise(
        mean_percent_wildtype_global = mean(percent_wildtype_global, na.rm = TRUE),
        sd_percent_wildtype_global = sd(percent_wildtype_global, na.rm = TRUE),
        mean_percent_difference_from_wildtype_global = mean(percent_difference_from_wildtype_global, na.rm = TRUE),
        sd_percent_difference_from_wildtype_global = sd(percent_difference_from_wildtype_global, na.rm = TRUE),
        mean_percent_difference_from_orc4r_global = mean(percent_difference_from_orc4r_global, na.rm = TRUE),
        sd_percent_difference_from_orc4r_global = sd(percent_difference_from_orc4r_global, na.rm = TRUE),
        mean_percent_change_from_wildtype_global = mean(percent_change_from_wildtype_global, na.rm = TRUE),
        sd_percent_change_from_wildtype_global = sd(percent_change_from_wildtype_global, na.rm = TRUE),
        mean_percent_change_from_orc4r_global = mean(percent_change_from_orc4r_global, na.rm = TRUE),
        sd_percent_change_from_orc4r_global = sd(percent_change_from_orc4r_global, na.rm = TRUE),
        mean_fold_change_from_orc4r_global = mean(fold_change_from_orc4r_global, na.rm = TRUE),
        sd_fold_change_from_orc4r_global = sd(fold_change_from_orc4r_global, na.rm = TRUE),
        replicate_count = n(),
        .groups = "drop"
    )
message("Global baseline summary statistics computed.")

stopifnot(
    "Not all global summary groups have 3 replicates." =
        all(summary_loading_data_global$replicate_count == 3)
)

# ==============================================================================
# C11: Assumption diagnostics (kglut, per-salt + pooled residuals)
# ==============================================================================
# Mirrors the 350mm diagnostics (Thread 2a C6): per-cell {n, mean, sd, variance}
# plus pooled-residual Shapiro. ADDS a per-salt split (per-salt residuals AND
# pooled residuals) because kglut spans three salts.
#
# I4 CAVEAT: at n=3 per cell, normality CANNOT be meaningfully tested. The
# Shapiro statistics below are DECORATIVE. Justification for the parametric
# paired t-tests rests on the paired/blocked design and the assay's established
# track record, NOT on a passing normality test.
#
# Bartlett/Levene are OMITTED: under per-kglut normalization WT is pinned to 100
# at every salt (zero within-cell variance), making those tests undefined
# (mirrors 2a). WT residuals below are therefore exactly 0 -- a normalization
# artifact, not biology.
#
# Diagnostics computed on percent_wildtype_floored (the column the C12 per-salt
# tests use).

diagnostic_data <- loading_data %>%
    filter(label %in% c("WT", "ORC4R", "+4sofa")) %>%
    droplevels()

stopifnot(
    "C11 diagnostics: percent_wildtype_floored column missing." =
        "percent_wildtype_floored" %in% names(diagnostic_data)
)

diagnostics_cell_stats <- diagnostic_data %>%
    group_by(kglut, label) %>%
    summarise(
        n = n(),
        mean_percent_wildtype_floored = mean(percent_wildtype_floored),
        sd_percent_wildtype_floored = sd(percent_wildtype_floored),
        variance_percent_wildtype_floored = var(percent_wildtype_floored),
        .groups = "drop"
    )
message("C11: per-cell descriptive stats computed.")

# Residual = value - (label x salt) cell mean.
diagnostic_data <- diagnostic_data %>%
    group_by(kglut, label) %>%
    mutate(cell_residual_percent_wildtype_floored =
               percent_wildtype_floored - mean(percent_wildtype_floored)) %>%
    ungroup()

pooled_residual_vector <- diagnostic_data$cell_residual_percent_wildtype_floored
pooled_residual_shapiro <- shapiro.test(pooled_residual_vector)

diagnostics_shapiro_rows <- list()
diagnostics_shapiro_rows[["pooled_all_salts"]] <- data.frame(
    residual_scope = "pooled_all_salts",
    shapiro_W = unname(pooled_residual_shapiro$statistic),
    shapiro_p = pooled_residual_shapiro$p.value,
    n_residuals = length(pooled_residual_vector),
    caveat = "n=3 per cell; normality untestable; DECORATIVE only",
    stringsAsFactors = FALSE
)
for (one_salt in levels(diagnostic_data$kglut)) {
    salt_residuals <- diagnostic_data$cell_residual_percent_wildtype_floored[
        diagnostic_data$kglut == one_salt]
    salt_residual_shapiro <- shapiro.test(salt_residuals)
    diagnostics_shapiro_rows[[one_salt]] <- data.frame(
        residual_scope = paste0("salt_", one_salt),
        shapiro_W = unname(salt_residual_shapiro$statistic),
        shapiro_p = salt_residual_shapiro$p.value,
        n_residuals = length(salt_residuals),
        caveat = "n=3 per cell; normality untestable; DECORATIVE only",
        stringsAsFactors = FALSE
    )
}
diagnostics_shapiro <- do.call(rbind, diagnostics_shapiro_rows)
rownames(diagnostics_shapiro) <- NULL
message("C11: pooled and per-salt Shapiro computed (DECORATIVE; see I4 caveat).")

diagnostics_cellstats_csv_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_diagnostics-cellstats.csv")
diagnostics_shapiro_csv_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_diagnostics-shapiro.csv")
if (!file.exists(diagnostics_cellstats_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(as.data.frame(diagnostics_cell_stats),
              diagnostics_cellstats_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(diagnostics_cellstats_csv_output_filepath))
}
if (!file.exists(diagnostics_shapiro_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(diagnostics_shapiro,
              diagnostics_shapiro_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(diagnostics_shapiro_csv_output_filepath))
}

# ==============================================================================
# C12: Per-salt paired tests (kglut) -- SIMPLE FRAMING, direct reviewer answer
# ==============================================================================
# For each salt (250/300/350), paired t-tests blocked on replicate:
#   WT vs ORC4R     (DIRECT reviewer answer: is loading reduced in ORC4R?)
#   WT vs +4sofa    (is +4sofa restored toward WT?)
#   ORC4R vs +4sofa (NARRATIVE-FIREWALL ordinary effect size -- see below)
#
# HOLM FAMILY (I2): the correction family is the WT-CONTRASTS within a SINGLE
# salt = { WT-vs-ORC4R, WT-vs-+4sofa }. Holm is applied WITHIN each salt across
# exactly these two contrasts. NEVER pooled across the three salts -- each salt
# is a separate question, hence a separate family. ORC4R-vs-+4sofa is NOT a WT
# contrast and is reported with its raw p, OUTSIDE the Holm family.
#
# NARRATIVE FIREWALL (handoff F6): the (+4sofa - ORC4R) comparison is an
# ORDINARY loading effect size. It is NOT the ATPase floor-proximity argument.
# In a loading assay, restorers/reducers deviate from ORC4R by design.
#
# I4 CAVEAT: parametric paired t at n=3; normality untestable; justification =
# paired/blocked design + assay track record, not a normality test. Wilcoxon
# companions are DECORATIVE -- at n=3 paired the two-sided exact p floors at
# ~0.25 and can never reach 0.05; reported for transparency, never leaned on.
#
# Tests run on percent_wildtype_floored (per-kglut normalization). No test
# genotype carries a negative value, so floored == raw for every value here.

per_salt_contrast_rows <- list()
salts_to_test <- levels(diagnostic_data$kglut)
contrast_definitions <- list(
    c("WT", "ORC4R"),
    c("WT", "+4sofa"),
    c("ORC4R", "+4sofa")
)

for (one_salt in salts_to_test) {
    salt_subset <- diagnostic_data %>% filter(kglut == one_salt) %>% arrange(replicate)
    for (one_contrast in contrast_definitions) {
        group_a_label <- one_contrast[1]
        group_b_label <- one_contrast[2]
        group_a_values <- salt_subset$percent_wildtype_floored[salt_subset$label == group_a_label]
        group_b_values <- salt_subset$percent_wildtype_floored[salt_subset$label == group_b_label]
        group_a_reps   <- salt_subset$replicate[salt_subset$label == group_a_label]
        group_b_reps   <- salt_subset$replicate[salt_subset$label == group_b_label]

        stopifnot(
            "C12 pairing: unequal pair lengths." =
                length(group_a_values) == length(group_b_values),
            "C12 pairing: expected n=3 per group." =
                length(group_a_values) == 3,
            "C12 pairing: replicate alignment broken (pairs not in same block order)." =
                identical(group_a_reps, group_b_reps)
        )

        paired_t <- t.test(group_a_values, group_b_values, paired = TRUE)
        paired_wilcox <- suppressWarnings(
            wilcox.test(group_a_values, group_b_values, paired = TRUE))
        is_wt_contrast <- group_a_label == "WT"

        per_salt_contrast_rows[[paste(one_salt, group_a_label, group_b_label, sep = "_")]] <-
            data.frame(
                salt_kglut = one_salt,
                contrast = paste0(group_a_label, "_vs_", group_b_label),
                group_a = group_a_label,
                group_b = group_b_label,
                n_pairs = length(group_a_values),
                mean_difference_a_minus_b = unname(paired_t$estimate),
                ci_low = paired_t$conf.int[1],
                ci_high = paired_t$conf.int[2],
                t_statistic = unname(paired_t$statistic),
                df = unname(paired_t$parameter),
                raw_p_paired_t = paired_t$p.value,
                wilcoxon_p_decorative = paired_wilcox$p.value,
                is_wt_contrast = is_wt_contrast,
                holm_family = if (is_wt_contrast)
                    paste0("WT-contrasts_salt_", one_salt) else "not_in_WT_family",
                test = "paired t-test (block = replicate)",
                normalization_column = "percent_wildtype_floored (per-kglut, WT=100/salt)",
                stringsAsFactors = FALSE
            )
    }
}
per_salt_contrasts <- do.call(rbind, per_salt_contrast_rows)
rownames(per_salt_contrasts) <- NULL

# Holm WITHIN each salt across the WT-contrasts ONLY (I2). Non-WT contrasts: NA.
per_salt_contrasts$holm_p_within_salt_wt_family <- NA_real_
for (one_salt in salts_to_test) {
    wt_family_index <- which(per_salt_contrasts$salt_kglut == one_salt &
                             per_salt_contrasts$is_wt_contrast)
    per_salt_contrasts$holm_p_within_salt_wt_family[wt_family_index] <-
        p.adjust(per_salt_contrasts$raw_p_paired_t[wt_family_index], method = "holm")
}
message("C12: per-salt paired tests + within-salt Holm (WT family) computed.")

# ==============================================================================
# C13: Salt-sensitivity interaction (kglut) -- STRUCTURED FRAMING, REPORTED
# ==============================================================================
# REPORTED result, NOT a diagnostic (I5 note): it answers the reviewer's salt-
# sensitivity question by quantifying how the WT-vs-genotype gap CHANGES with
# salt. C12 is the simple framing; this is the structured one. BOTH are reported.
#
# Model: percent_wildtype_global ~ label * kglut, random = ~1 | replicate (lme).
# salt (kglut) as FACTOR (3 levels, no monotonic/trend assumption) -- DEFAULT.
#   Alternative considered: salt as numeric/ordered (linear trend across
#   250/300/350). NOT adopted: only 3 levels, the decline need not be linear,
#   and the factor model is the assumption-light conservative choice.
#
# I3 NORMALIZATION FIREWALL: the interaction MUST use the GLOBAL-baseline column
# (WT@250mM). Under per-kglut normalization WT is pinned to 100 at every salt
# (within-WT sd across salts = 0) so WT cannot show a salt effect -> the
# genotype:salt interaction is MECHANICALLY DEGENERATE. The firewall below
# confirms (membership + VALUE + name) the response is the global column and
# that WT is NOT pinned flat across salts.

# *** SANCTIONED DELIBERATE RED (C13 -- the ONLY one in the project) ***
# GREEN (committed): the global-baseline floored column.
# To WITNESS the firewall fire, comment the GREEN line and uncomment the RED
# line (per-kglut column, WT pinned to 100). The VALUE check then fires
# (WT sd = 0), proving the guard is real -- not merely a name check. Restore
# GREEN to proceed.
interaction_response_column <- "percent_wildtype_global_floored"   # GREEN (committed)
# interaction_response_column <- "percent_wildtype_floored"        # RED (witness guard)

interaction_data <- loading_data %>%
    filter(label %in% c("WT", "ORC4R", "+4sofa")) %>%
    droplevels() %>%
    mutate(replicate_factor = factor(replicate))

# Membership FIRST (carry-forward #1): a missing column makes [[ ]] return NULL.
stopifnot(
    "C13 firewall: chosen response column absent from interaction data." =
        interaction_response_column %in% names(interaction_data)
)
wt_values_in_chosen_column <- interaction_data[[interaction_response_column]][
    interaction_data$label == "WT"]
# VALUE check is the load-bearing one; NAME check is belt-and-suspenders.
stopifnot(
    "C13 firewall (VALUE): WT is pinned flat across salts -> per-kglut column fed to the interaction; use the global-baseline column (I3)." =
        sd(wt_values_in_chosen_column) > 1,
    "C13 firewall (NAME): response column is not a global-baseline column (I3)." =
        interaction_response_column %in%
            c("percent_wildtype_global", "percent_wildtype_global_floored")
)
message("C13 firewall passed: interaction runs on '", interaction_response_column,
        "' (WT sd across salts = ", round(sd(wt_values_in_chosen_column), 2), ").")

interaction_model_formula <- as.formula(paste0(interaction_response_column, " ~ label * kglut"))

interaction_model_fit <- tryCatch(
    nlme::lme(interaction_model_formula,
              random = ~ 1 | replicate_factor,
              data = interaction_data, method = "REML"),
    error = function(e) {
        message("C13: nlme::lme did NOT converge -> reporting as a finding, not a code bug:\n  ",
                conditionMessage(e))
        NULL
    }
)

if (is.null(interaction_model_fit)) {
    interaction_results <- data.frame(
        term = NA_character_, estimate = NA_real_, std_error = NA_real_,
        df = NA_real_, t_value = NA_real_, p_value = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_,
        normalization_column = "percent_wildtype_global_floored (global, WT@250mM)",
        note = "nlme::lme did not converge; see console.", stringsAsFactors = FALSE
    )
} else {
    interaction_fixed_table <- summary(interaction_model_fit)$tTable
    interaction_intervals_fixed <- tryCatch(
        nlme::intervals(interaction_model_fit, which = "fixed")$fixed,
        error = function(e) {
            message("C13: intervals() failed (", conditionMessage(e), "); CI columns set NA.")
            NULL
        }
    )
    interaction_results <- data.frame(
        term = rownames(interaction_fixed_table),
        estimate = interaction_fixed_table[, "Value"],
        std_error = interaction_fixed_table[, "Std.Error"],
        df = interaction_fixed_table[, "DF"],
        t_value = interaction_fixed_table[, "t-value"],
        p_value = interaction_fixed_table[, "p-value"],
        stringsAsFactors = FALSE
    )
    rownames(interaction_results) <- NULL
    if (is.null(interaction_intervals_fixed)) {
        interaction_results$ci_low <- NA_real_
        interaction_results$ci_high <- NA_real_
    } else {
        interaction_results$ci_low <-
            interaction_intervals_fixed[interaction_results$term, "lower"]
        interaction_results$ci_high <-
            interaction_intervals_fixed[interaction_results$term, "upper"]
    }
    interaction_results$normalization_column <-
        "percent_wildtype_global_floored (global, WT@250mM)"
    interaction_results$note <-
        "salt as factor; genotype:salt term = salt-sensitivity of the WT gap"
}
message("C13: salt-sensitivity interaction model fit; fixed-effects table extracted.")

# ==============================================================================
# C14a: Results CSVs -- BOTH framings, with normalization-column lineage
# ==============================================================================
per_salt_contrasts_csv_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_per-salt-contrasts.csv")
salt_interaction_csv_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_salt-interaction.csv")

if (!file.exists(per_salt_contrasts_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(per_salt_contrasts, per_salt_contrasts_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(per_salt_contrasts_csv_output_filepath))
}
if (!file.exists(salt_interaction_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(interaction_results, salt_interaction_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(salt_interaction_csv_output_filepath))
}

# ==============================================================================
# C14b prep: annotation data + interaction caption (built before plot)
# ==============================================================================
wt_vs_orc4r_per_salt_annotation <- per_salt_contrasts %>%
    filter(contrast == "WT_vs_ORC4R") %>%
    transmute(
        kglut = factor(salt_kglut, levels = levels(loading_data$kglut)),
        annotation_text = paste0(
            "WT vs ORC4R\np = ", formatC(raw_p_paired_t, format = "f", digits = 4),
            "\nHolm p = ", formatC(holm_p_within_salt_wt_family, format = "f", digits = 4)
        )
    )

salt_sensitivity_term_for_caption <- "labelORC4R:kglut350"
if (!is.null(interaction_model_fit) &&
    salt_sensitivity_term_for_caption %in% interaction_results$term) {
    caption_estimate <- interaction_results$estimate[interaction_results$term == salt_sensitivity_term_for_caption]
    caption_ci_low   <- interaction_results$ci_low[interaction_results$term == salt_sensitivity_term_for_caption]
    caption_ci_high  <- interaction_results$ci_high[interaction_results$term == salt_sensitivity_term_for_caption]
    caption_p        <- interaction_results$p_value[interaction_results$term == salt_sensitivity_term_for_caption]
    salt_sensitivity_caption_text <- paste0(
        "Per-panel: WT vs ORC4R exact paired-t p (Holm-adjusted within salt; family = WT contrasts).\n",
        "Salt-sensitivity interaction (nlme genotype x salt, global baseline WT@250mM): ",
        "WT-ORC4R gap change 250->350 mM = ",
        formatC(caption_estimate, format = "f", digits = 2), " %WT-units, 95% CI [",
        formatC(caption_ci_low, format = "f", digits = 2), ", ",
        formatC(caption_ci_high, format = "f", digits = 2), "], p = ",
        formatC(caption_p, format = "f", digits = 4), "."
    )
} else {
    salt_sensitivity_caption_text <-
        "Salt-sensitivity interaction: model did not converge (see console)."
}
message("C14b: per-salt annotation and interaction caption prepared.")

# ==============================================================================
# Plots
# ==============================================================================

# Primary plot. Exploratory plots are in
# quantification_kgluttitr_wt-4r-ps_exploratory-plots.R
# C14b: per-salt WT-vs-ORC4R EXACT p annotated above each panel; salt-sensitivity
# interaction estimate/CI/p bound into the caption. EXACT p-values, never stars.
faceted_by_kglut_plot <- ggplot(summary_loading_data, aes(x = label, y = mean_percent_wildtype, fill = label)) +
    geom_col(
        width = PLOT_CONFIG$bar$width,
        color = PLOT_CONFIG$bar$color,
        linewidth = PLOT_CONFIG$bar$linewidth
    ) +
    geom_errorbar(
        aes(
            # pmax(0, ...): negative MCM loading is biologically meaningless.
            # Clamps lower error bar at 0 to prevent SD overshoot visual artifacts.
            ymin = pmax(0, mean_percent_wildtype - sd_percent_wildtype),
            ymax = mean_percent_wildtype + sd_percent_wildtype
        ),
        width = PLOT_CONFIG$errorbar$width,
        linewidth = PLOT_CONFIG$errorbar$linewidth
    ) +
    geom_point(
        data = loading_data %>%
            filter(label %in% c("WT", "ORC4R", "+4sofa")) %>%
            droplevels(),
        aes(x = label, y = percent_wildtype, shape = factor(.data$replicate)),
        position = position_jitter(
            width = PLOT_CONFIG$point$jitter_width,
            # seed = 42: deterministic jitter for reproducible figures across runs.
            seed = PLOT_CONFIG$point$jitter_seed
        ),
        size = PLOT_CONFIG$point$size,
        fill = PLOT_CONFIG$point$fill,
        color = PLOT_CONFIG$point$color,
        stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE
    ) +
    geom_text(
        data = wt_vs_orc4r_per_salt_annotation,
        aes(x = 1.5, y = Inf, label = annotation_text),
        inherit.aes = FALSE, vjust = 1.15, hjust = 0.5,
        size = 3, family = PLOT_CONFIG$theme$base_family, lineheight = 0.9
    ) +
    facet_wrap(~kglut, nrow = 1, labeller = label_both) +
    scale_fill_manual(values = PLOT_CONFIG$fill_colors) +
    scale_shape_manual(
        values = PLOT_CONFIG$replicate_shapes,
        labels = c("1" = "Replicate 1", "2" = "Replicate 2", "3" = "Replicate 3"),
        name = "Replicate"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        x = "Sample",
        y = "MCM Loading (% WT)",
        title = "MCM Loading Across KGlut Concentrations",
        fill = "Sample",
        caption = salt_sensitivity_caption_text
    ) +
    theme_classic(
        base_size = PLOT_CONFIG$theme$base_size,
        base_family = PLOT_CONFIG$theme$base_family
    ) +
    theme(
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = PLOT_CONFIG$theme$legend_position,
        plot.caption = element_text(hjust = 0, size = 8),
        panel.spacing = unit(1, "lines")
    )
message("Faceted-by-kglut plot constructed.")

# Global baseline plot: same layout as faceted_by_kglut_plot but mapping
# percent_wildtype_global. WT bars decrease across panels (not flat at 100%).
# Figure legend note: "All values normalized to WT MCM loading at 250 mM KGlut.
# Error bars represent +/-1 SD of three biological replicates. Individual
# replicates shown as distinct shapes."
global_baseline_plot <- ggplot(summary_loading_data_global,
    aes(x = label, y = mean_percent_wildtype_global, fill = label)) +
    geom_col(
        width = PLOT_CONFIG$bar$width,
        color = PLOT_CONFIG$bar$color,
        linewidth = PLOT_CONFIG$bar$linewidth
    ) +
    geom_errorbar(
        aes(
            # pmax(0, ...): negative MCM loading is biologically meaningless.
            # Clamps lower error bar at 0 to prevent SD overshoot visual artifacts.
            ymin = pmax(0, mean_percent_wildtype_global - sd_percent_wildtype_global),
            ymax = mean_percent_wildtype_global + sd_percent_wildtype_global
        ),
        width = PLOT_CONFIG$errorbar$width,
        linewidth = PLOT_CONFIG$errorbar$linewidth
    ) +
    geom_point(
        data = loading_data %>%
            filter(label %in% c("WT", "ORC4R", "+4sofa")) %>%
            droplevels(),
        aes(x = label, y = percent_wildtype_global, shape = factor(.data$replicate)),
        position = position_jitter(
            width = PLOT_CONFIG$point$jitter_width,
            # seed = 42: deterministic jitter for reproducible figures across runs.
            seed = PLOT_CONFIG$point$jitter_seed
        ),
        size = PLOT_CONFIG$point$size,
        fill = PLOT_CONFIG$point$fill,
        color = PLOT_CONFIG$point$color,
        stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE
    ) +
    facet_wrap(~kglut, nrow = 1, labeller = label_both) +
    scale_fill_manual(values = PLOT_CONFIG$fill_colors) +
    scale_shape_manual(
        values = PLOT_CONFIG$replicate_shapes,
        labels = c("1" = "Replicate 1", "2" = "Replicate 2", "3" = "Replicate 3"),
        name = "Replicate"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Sample",
        y = "MCM Loading (% WT at 250 mM)",
        title = "MCM Loading Across KGlut Concentrations (Global Baseline)",
        fill = "Sample"
    ) +
    theme_classic(
        base_size = PLOT_CONFIG$theme$base_size,
        base_family = PLOT_CONFIG$theme$base_family
    ) +
    theme(
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = PLOT_CONFIG$theme$legend_position,
        panel.spacing = unit(1, "lines")
    )
message("Global baseline plot constructed.")

# ==============================================================================
# Output
# ==============================================================================

plot_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_kglut-titration_per-kglut-plot.pdf")
if (!file.exists(plot_output_filepath) || OVERWRITE_PLOTS) {
  ggsave(
       plot_output_filepath,
       faceted_by_kglut_plot,
       device = PLOT_CONFIG$output$device,
       width = PLOT_CONFIG$output$width,
       height = PLOT_CONFIG$output$height
   )
  message("Saved plot: ", basename(plot_output_filepath))
} else {
  message("Skipped plot (already exists): ", basename(plot_output_filepath))
}

# Save Summary CSV
if (!file.exists(summary_csv_output_filepath) || OVERWRITE_CSVS) {
  write.csv(summary_loading_data, summary_csv_output_filepath, row.names = FALSE)
  message("Saved CSV: ", basename(summary_csv_output_filepath))
} else {
  message("Skipped CSV (already exists): ", basename(summary_csv_output_filepath))
}

# Save Full Data CSV
if (!file.exists(full_data_csv_output_filepath) || OVERWRITE_CSVS) {
  write.csv(loading_data, full_data_csv_output_filepath, row.names = FALSE)
  message("Saved CSV: ", basename(full_data_csv_output_filepath))
} else {
  message("Skipped CSV (already exists): ", basename(full_data_csv_output_filepath))
}

# Global baseline plot output
global_baseline_plot_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_global-baseline-plot.pdf"
)
if (!file.exists(global_baseline_plot_output_filepath) || OVERWRITE_PLOTS) {
    ggsave(
        global_baseline_plot_output_filepath,
        global_baseline_plot,
        device = PLOT_CONFIG$output$device,
        width = PLOT_CONFIG$output$width,
        height = PLOT_CONFIG$output$height
    )
    message("Saved plot: ", basename(global_baseline_plot_output_filepath))
} else {
    message("Skipped plot (already exists): ", basename(global_baseline_plot_output_filepath))
}

# Global baseline summary CSV
global_baseline_summary_csv_output_filepath <- file.path(
    OUTPUT_DIRECTORY, "loading_kglut-titration_global-baseline-summary.csv"
)
if (!file.exists(global_baseline_summary_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(summary_loading_data_global, global_baseline_summary_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(global_baseline_summary_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(global_baseline_summary_csv_output_filepath))
}

message("Script complete.")
