# Date updated: 2026-01-12
# Data of the xlsx was produced by analyzing tiff files using imagej.
# Usage: source("loading-quantification_kgluttitr_wt-4r-ps.R")
# The output is the plot called faceted_by_kglut_plot.
# All other plots kept for reference.

library(xlsx)
library(tidyverse)

OVERWRITE_PLOTS <- FALSE
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
COLUMN_TO_REMOVE <- "NA." # Column has row numbers from imagej export/copy-paste
# @NOTE: An additional space in the excel sheet cell with cause column name to have dot at the end.
#   > "ORC " would be read in as "ORC."
REQUIRED_COLUMNS <- c(
  "NA.", "Intensity", "Lane",
  "ORC", "Suppressor",
  "kGlut", "Input"
)
# Define custom orderings
factor_order <- list(
  "Suppressor" = c("None", "1EK", "3PL", "4PS"),
  "kGlut" = c("250", "300", "350"),
  "Label" = c("WT", "ORC4R",  "+1sofr",  "+3sofr", "+4sofr")
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
  temp_df <- read.xlsx(FILE_PATH, header = TRUE, sheetIndex = sheet_idx)

  stopifnot(
    "Number of rows in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      nrow(temp_df) == EXPECTED_NUMBER_OF_ROWS,
    "Number of columns in temp_df does not match EXPECTED_NUMBER_OF_ROWS." =
      ncol(temp_df) == EXPECTED_NUMBER_OF_COLS,
    "df_lst[[df_count]] colnames not identical to REQUIRED_COLUMNS." =
      identical(REQUIRED_COLUMNS, colnames(temp_df))
  )


temp_df <- temp_df %>%
  mutate(
    # Subtract background (Row 2) from all rows
    Net_Intensity = Intensity - Intensity[2],
    Experiment = paste0("Exp_", df_count)
  ) %>%
  # Correct for input volume (0.5 factor) using the Net_Intensity
  mutate(Rel_to_Input = (Net_Intensity / Net_Intensity[Input == "yes"]) * 0.5) %>%
  filter(Input == "no")
  df_lst[[df_count]] <- temp_df


} # end read-data for loop

loading_df <- do.call(rbind, df_lst)
loading_df <- loading_df[, names(loading_df) != COLUMN_TO_REMOVE]

message("loading_df preparation complete...")

loading_df <- loading_df %>%
  mutate(Label = case_when(
    ORC == "None" & Suppressor == "None" ~ "None",
    ORC == "WT" & Suppressor == "None" ~ "WT",
    ORC == "RA" & Suppressor == "None" ~ "ORC4R",
    ORC == "RA" & Suppressor == "4PS" ~ "+4sofr",
    ORC == "RA" & Suppressor == "1EK" ~ "+1sofr",
    ORC == "RA" & Suppressor == "3PL" ~ "+3sofr"
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
  group_by(Experiment, kGlut) %>%
  mutate(
    Norm_to_WT_kGlut = Rel_to_Input / Rel_to_Input[ORC == "WT"][1]
  ) %>%
  
  # Step 2: Normalize to the WT of the 250 mM concentration (Global Baseline)
  group_by(Experiment) %>%
  mutate(
    Norm_to_WT_250 = Rel_to_Input / Rel_to_Input[ORC == "WT" & kGlut == "250"][1]
  ) %>%
  ungroup()

df_summary <- loading_df %>%
  filter(
    ORC != "None",                       # Exclude negative control from plots
    !(Label %in% c("+1sofr", "+3sofr")), # Exclude other specific mutants
    !is.na(Label)                        # Remove any remaining NA labels
  ) %>%
  group_by(kGlut, Label) %>%
  summarise(
    across(
      c(Rel_to_Input, Norm_to_WT_kGlut, Norm_to_WT_250),
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Calculate fold change relative to WT
df_normalized <- df_summary %>%
  group_by(kGlut) %>%
  mutate(
    # Ensure a single value is pulled for WT reference
    wt_val = Rel_to_Input_mean[Label == "WT"][1],
    wt_sd = Rel_to_Input_sd[Label == "WT"][1],
    fold_change = Rel_to_Input_mean / wt_val,
    # Proper error propagation
    fold_change_sd = fold_change * sqrt((Rel_to_Input_sd/Rel_to_Input_mean)^2 + (wt_sd/wt_val)^2)
  ) %>%
  ungroup()

intensity_vs_kglut_plot <- ggplot(df_summary, aes(x = kGlut, y = Rel_to_Input_mean, color = Label, fill = Label)) +
  geom_line(linewidth = 1.2, aes(group = Label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.15, linewidth = 0.6) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "kGlut Concentration (mM)", 
       y = "MCM (pmol)", 
       title = "Intensity vs kGlut",
       color = "Protein",
       fill = "Protein") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

all_samples_dodged_plot <- ggplot(df_summary, aes(x = Label, y = Rel_to_Input_mean, fill = Label, group = kGlut)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, 
           color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd),
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sample Type", 
       y = "MCM (pmol)", 
       title = "MCM Levels by Sample Type and kGlut Concentration", 
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold")
  )

# This is the correct plot. All other plots kept for reference.
faceted_by_kglut_plot <- ggplot(df_summary, aes(x = Label, y = Rel_to_Input_mean, fill = Label)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.25, linewidth = 0.6) +
  facet_wrap(~kGlut, nrow = 1, labeller = label_both) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sample Type", 
       y = "MCM (pmol)", 
       title = "Label Effects Across kGlut Concentrations",
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

df_summary_300 <- df_summary %>% filter(kGlut == "300")

kglut_300_only_plot <- ggplot(df_summary_300, aes(x = Label, y = Rel_to_Input_mean, fill = Label)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.25, linewidth = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Sample Type", 
       y = "MCM (pmol)", 
       title = "300 mM kGlut Only", 
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold")
  )

normalized_to_wt_plot <- ggplot(df_normalized, aes(x = kGlut, y = fold_change, color = Label, fill = Label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1.2, aes(group = Label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_errorbar(aes(ymin = fold_change - fold_change_sd, 
                    ymax = fold_change + fold_change_sd), 
                width = 0.15, linewidth = 0.6) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "kGlut Concentration (mM)", 
       y = "Fold Change (relative to WT)", 
       title = "Normalized to Wildtype",
       color = "Sample",
       fill = "Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )

faceted_by_label_plot <- ggplot(df_summary, aes(x = kGlut, y = Rel_to_Input_mean, fill = kGlut)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Rel_to_Input_mean - Rel_to_Input_sd, ymax = Rel_to_Input_mean + Rel_to_Input_sd), 
                width = 0.25, linewidth = 0.6) +
  facet_wrap(~Label, nrow = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "kGlut Concentration (mM)", 
       y = "MCM (pmol)", 
       title = "kGlut Dose-Response by Sample Type",
       fill = "kGlut") +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  )

message("Plotting completed...")

plot_object_names <- ls(pattern = "_plot$", envir = .GlobalEnv)
file_paths <- normalizePath(
  file.path(OUTPUT_DIRECTORY, paste0(plot_object_names, ".pdf")),
  mustWork = FALSE
)
names(file_paths) <- plot_object_names

# Save using names
# Update the plot saving loop
for (plot_name in names(file_paths)) {
  if (!file.exists(file_paths[plot_name]) || OVERWRITE_PLOTS) {
    current_plot <- get(plot_name, envir = .GlobalEnv)
    ggsave(
      file_paths[plot_name],
      current_plot,
      width = 8,
      height = 5
    )
    cat("Saved plot:", basename(file_paths[plot_name]), "\n")
  } else {
    cat("Skipped plot (already exists):", basename(file_paths[plot_name]), "\n")
  }
}

# Save Summary CSV
if (!file.exists(summary_csv_path) || OVERWRITE_CSVS) {
  write.csv(df_summary, summary_csv_path, row.names = FALSE)
  cat("Saved CSV:", basename(summary_csv_path), "\n")
} else {
  cat("Skipped CSV (already exists):", basename(summary_csv_path), "\n")
}

# Save Full Data CSV
if (!file.exists(full_data_csv_path) || OVERWRITE_CSVS) {
  write.csv(loading_df, full_data_csv_path, row.names = FALSE)
  cat("Saved CSV:", basename(full_data_csv_path), "\n")
} else {
  cat("Skipped CSV (already exists):", basename(full_data_csv_path), "\n")
}

message("Script complete.")
