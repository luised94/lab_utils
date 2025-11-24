library(xlsx)
library(tidyverse)

FILENAME <- "Analysis.xlsx"
OUTPUT_DIRECTORY <- "~/data/loading_analysis"
EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication"
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
  "kGlut" = c("250", "300", "350"), # In inverse order (as factor, alphabetical would be "250", ...)
  "Label" = c("WT", "ORC4R",  "+1sofr",  "+3sofr", "+4sofr") # Adjust this as needed
)

DROPBOX_PATH <- Sys.getenv("DROPBOX_PATH")
FILE_PATH <- file.path(DROPBOX_PATH, EXPERIMENT_DIRECTORY, FILENAME)

if (nchar(DROPBOX_PATH) == 0) {
  stop("DROPBOX_PATH not defined: ")

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
      Intensity = Intensity - Intensity[2],
      Experiment = paste0("Exp_", df_count)
    ) %>%
    slice(-2) %>%
    mutate(Intensity = (Intensity / Intensity[Input == "yes"]) * 0.5) %>%
    filter(Input == "no")

  df_lst[[df_count]] <- temp_df


} # end read-data for loop

loading_df <- do.call(rbind, df_lst)
loading_df <- loading_df[, names(loading_df) != COLUMN_TO_REMOVE]

message("loading_df preparation complete...")

loading_df <- loading_df %>%
  mutate(Label = case_when(
    ORC == "WT" & Suppressor == "None" ~ "WT",
    ORC == "RA" & Suppressor == "None" ~ "ORC4R",
    ORC == "RA" & Suppressor == "4PS" ~ "+4sofr",
    ORC == "RA" & Suppressor == "1EK" ~ "+1sofr",
    ORC == "RA" & Suppressor == "3PL" ~ "+3sofr"
  ))

for (col_name in names(factor_order)) {

  loading_df[[col_name]] <- factor(
    loading_df[[col_name]],
    levels = factor_order[[col_name]],
    ordered = TRUE
  )
}

message("Adjust loading_df to factor order...")
#sapply(loading_df, function(x) if(!is.numeric(x)) unique(x))

# Calculate mean and sd for each kGlut
df_summary <- loading_df %>% 
  filter(!(Label %in% c("+1sofr","+3sofr"))) %>%
  group_by(kGlut, Label) %>%
  summarise(pmol_MCM = mean(Intensity), 
            sd = sd(Intensity))

intensity_vs_kglut_plot <- ggplot(df_summary, aes(x = kGlut, y = pmol_MCM, color = Label, fill = Label)) +
  geom_line(size = 1.2, aes(group = Label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd), 
                width = 0.15, size = 0.6) +
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

all_samples_dodged_plot <- ggplot(df_summary, aes(x = Label, y = pmol_MCM, fill = Label, group = kGlut)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, 
           color = "black", size = 0.4) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd),
                position = position_dodge(width = 0.8), 
                width = 0.25, size = 0.6) +
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

faceted_by_kglut_plot <- ggplot(df_summary, aes(x = Label, y = pmol_MCM, fill = Label)) +
  geom_col(width = 0.7, color = "black", size = 0.4) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd), 
                width = 0.25, size = 0.6) +
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

kglut_300_only_plot <- ggplot(df_summary_300, aes(x = Label, y = pmol_MCM, fill = Label)) +
  geom_col(width = 0.7, color = "black", size = 0.4) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd), 
                width = 0.25, size = 0.6) +
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

# Calculate fold change relative to WT
df_normalized <- df_summary %>%
  group_by(kGlut) %>%
  mutate(
    wt_value = pmol_MCM[Label == "WT"],
    fold_change = pmol_MCM / wt_value,
    # Propagate error for fold change (approximate)
    fold_change_sd = fold_change * sqrt((sd/pmol_MCM)^2 + (sd[Label == "WT"]/wt_value)^2)
  ) %>%
  ungroup()

normalized_to_wt_plot <- ggplot(df_normalized, aes(x = kGlut, y = fold_change, color = Label, fill = Label)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", size = 0.8) +
  geom_line(size = 1.2, aes(group = Label)) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_errorbar(aes(ymin = fold_change - fold_change_sd, 
                    ymax = fold_change + fold_change_sd), 
                width = 0.15, size = 0.6) +
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

faceted_by_label_plot <- ggplot(df_summary, aes(x = kGlut, y = pmol_MCM, fill = kGlut)) +
  geom_col(width = 0.7, color = "black", size = 0.4) +
  geom_errorbar(aes(ymin = pmol_MCM - sd, ymax = pmol_MCM + sd), 
                width = 0.25, size = 0.6) +
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
# Your pattern-based saving approach
plot_object_names <- ls(pattern = "_plot$", envir = .GlobalEnv)
file_paths <- normalizePath(
  file.path(OUTPUT_DIRECTORY, paste0(plot_object_names, ".pdf")),
  mustWork = FALSE
)
names(file_paths) <- plot_object_names

# Save using names
for (plot_name in names(file_paths)) {
  current_plot <- get(plot_name, envir = .GlobalEnv)
  ggsave(file_paths[plot_name], current_plot, width = 8, height = 5)
  cat("Saved:", basename(file_paths[plot_name]), "\n")
}
message("Plots saved to OUTPUT_DIRECTORY...")
message("Script complete.")
