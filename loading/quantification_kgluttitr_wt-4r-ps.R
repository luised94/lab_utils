# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., readxl::read_excel).
# Date updated: 2026-04-22
# Data of the xlsx was produced by analyzing tiff files using imagej.
# Usage: source("loading-quantification_kgluttitr_wt-4r-ps.R")
# The output is the plot called faceted_by_kglut_plot.
# All other plots kept for reference.

library(readxl)
library(tidyverse)

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

FILENAME <- "Analysis.xlsx"
OUTPUT_DIRECTORY <- "~/data/loading_analysis"
EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication"
# Define CSV file paths
summary_csv_path <- file.path(OUTPUT_DIRECTORY, "loading_summary_statistics.csv")
full_data_csv_path <- file.path(OUTPUT_DIRECTORY, "loading_processed_full_data.csv")

sheet_indices <- c(2, 3, 4)
EXPECTED_NUMBER_OF_ROWS <- 14
EXPECTED_NUMBER_OF_COLS <- 7
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

MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")
FILE_PATH <- file.path(MC_DROPBOX_PATH, EXPERIMENT_DIRECTORY, FILENAME)

if (nchar(MC_DROPBOX_PATH) == 0) {
  stop(
    "MC_DROPBOX_PATH not defined in bash environment.\n",
    "Set manually via command line if necessary."
  )

}
if (!dir.exists(OUTPUT_DIRECTORY)) {
  dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)

}
if (!file.exists(FILE_PATH)) {
  stop("FILE_PATH does not exist: ", FILE_PATH)

}

message("Paths for input and output set...")
df_lst <- vector(mode = "list", length = length(sheet_indices))
temp_df <- data.frame(matrix(
  NA,
  nrow = EXPECTED_NUMBER_OF_ROWS,
  ncol = EXPECTED_NUMBER_OF_COLS
))
loading_df <- data.frame(matrix(
  NA,
  nrow = EXPECTED_NUMBER_OF_ROWS * length(sheet_indices),
  ncol = EXPECTED_NUMBER_OF_COLS
))

message("Reading in loading assay sheets...")

df_count <- 0
for (sheet_idx in sheet_indices){
  df_count <- df_count + 1

  temp_df <- readxl::read_excel(FILE_PATH, sheet = sheet_idx)

  stopifnot(
    "Number of rows in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      nrow(temp_df) == EXPECTED_NUMBER_OF_ROWS,
    "Number of columns in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      ncol(temp_df) == EXPECTED_NUMBER_OF_COLS,
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
      net_intensity = intensity - intensity[2],
      replicate = df_count
    ) %>%
    # Correct for input volume (0.5 factor) using the net_intensity
    mutate(relative_to_input = (net_intensity / net_intensity[input == "yes"]) * 0.5) %>%
    filter(input == "no")
  df_lst[[df_count]] <- temp_df

} # end read-data for loop

loading_df <- do.call(rbind, df_lst)
loading_df <- loading_df[, names(loading_df) != COLUMN_TO_REMOVE]

message("loading_df preparation complete...")

# Remove negative control rows (orc == "None") before label assignment.
# These rows have no meaningful label and would cause indexing errors in
# downstream paired calculations that subset by label.
loading_df <- loading_df %>%
    filter(orc != "None")
message("Negative control rows removed (orc == 'None').")


loading_df <- loading_df %>%
  mutate(label = case_when(
    orc == "WT" & suppressor == "None" ~ "WT",
    orc == "RA" & suppressor == "None" ~ "ORC4R",
    orc == "RA" & suppressor == "4PS" ~ "+4sofa",
    orc == "RA" & suppressor == "1EK" ~ "+1sofa",
    orc == "RA" & suppressor == "3PL" ~ "+3sofa",
    orc == "RA" & suppressor == "5EK" ~ "+5sofa"
  ))

message("Adjust loading_df to factor order...")
for (col_name in names(factor_order)) {

  loading_df[[col_name]] <- factor(
    loading_df[[col_name]],
    levels = factor_order[[col_name]],
    ordered = TRUE
  )
}

loading_df <- loading_df %>%
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
loading_df <- loading_df %>%
    mutate(percent_wildtype = normalized_to_wildtype_per_kglut * 100)

message("Percent wildtype column computed.")

# Derived calculations: paired within replicate and kglut so each condition
# is compared to the WT and ORC4R values from the same experiment and salt
# concentration. Formulas match Script 1 for cross-script comparability.
loading_df <- loading_df %>%
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


df_summary <- loading_df %>%
  filter(!(label %in% c("+1sofa", "+3sofa", "+5sofa"))) %>%
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

# Primary plot. Exploratory plots are in
# quantification_kgluttitr_wt-4r-ps_exploratory-plots.R
faceted_by_kglut_plot <- ggplot(df_summary, aes(x = label, y = mean_percent_wildtype, fill = label)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_percent_wildtype - sd_percent_wildtype, ymax = mean_percent_wildtype + sd_percent_wildtype), 
                width = 0.25, linewidth = 0.6) +
  facet_wrap(~kglut, nrow = 1, labeller = label_both) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sample",
       y = "MCM Loading (% WT)",
       title = "MCM Loading Across KGlut Concentrations",
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

message("Plotting completed...")

plot_filepath <- file.path(OUTPUT_DIRECTORY, "faceted_by_kglut_plot.pdf")
if (!file.exists(plot_filepath) || OVERWRITE_PLOTS) {
  ggsave(plot_filepath, faceted_by_kglut_plot, width = 8, height = 5)
  message("Saved plot: ", basename(plot_filepath))
} else {
  message("Skipped plot (already exists): ", basename(plot_filepath))
}

# Save Summary CSV
if (!file.exists(summary_csv_path) || OVERWRITE_CSVS) {
  write.csv(df_summary, summary_csv_path, row.names = FALSE)
  message("Saved CSV: ", basename(summary_csv_path))
} else {
  message("Skipped CSV (already exists): ", basename(summary_csv_path))
}

# Save Full Data CSV
if (!file.exists(full_data_csv_path) || OVERWRITE_CSVS) {
  write.csv(loading_df, full_data_csv_path, row.names = FALSE)
  message("Saved CSV: ", basename(full_data_csv_path))
} else {
  message("Skipped CSV (already exists): ", basename(full_data_csv_path))
}

message("Script complete.")
