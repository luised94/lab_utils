# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., readxl::read_excel).
# Date updated: 2026-04-22
# Data of the xlsx was produced by analyzing tiff files using imagej.
# Usage: source("orc4r-screen_loading-kglut-titration.R")
# The output is the plot called faceted_by_kglut_plot.
# All other plots kept for reference.
# Prerequisites:
#   1. Environment variable MC_DROPBOX_PATH must be set in shell profile.
#      Points to the root of the shared Dropbox folder.
#   2. Font setup (one-time per machine):
#        renv::install(c("extrafont", "ragg"))
#        library(extrafont); font_import(prompt = FALSE)
#      Per-session: extrafont::loadfonts(device = "pdf")
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
# Configuration
# ==============================================================================

MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")

if (nchar(MC_DROPBOX_PATH) == 0) {
  stop(
    "MC_DROPBOX_PATH not defined in bash environment.\n",
    "Set manually via command line if necessary."
  )

}

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
        "Run once: library(extrafont); font_import(prompt = FALSE)\n",
        "Then restart R and re-source this script."
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
INPUT_FILEPATH <- file.path(MC_DROPBOX_PATH, EXPERIMENT_DIRECTORY, INPUT_FILENAME)

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
# Plots
# ==============================================================================

# Primary plot. Exploratory plots are in
# quantification_kgluttitr_wt-4r-ps_exploratory-plots.R
faceted_by_kglut_plot <- ggplot(summary_loading_data, aes(x = label, y = mean_percent_wildtype, fill = label)) +
    geom_col(
        width = PLOT_CONFIG$bar$width,
        color = PLOT_CONFIG$bar$color,
        linewidth = PLOT_CONFIG$bar$linewidth
    ) +
    geom_errorbar(
        aes(
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
        y = "MCM Loading (% WT)",
        title = "MCM Loading Across KGlut Concentrations",
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
message("Plotting completed...")

# Global baseline plot: same layout as faceted_by_kglut_plot but mapping
# percent_wildtype_global. WT bars decrease across panels (not flat at 100%).
# Figure legend note: "All values normalized to WT MCM loading at 250 mM KGlut.
# Error bars represent ñ1 SD of three biological replicates. Individual
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
