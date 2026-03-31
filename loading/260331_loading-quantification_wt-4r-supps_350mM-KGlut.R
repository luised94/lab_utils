# Date created: 2026-03-31
# Data produced by analyzing tiff files using ImageJ. Manual gel processing
# due to noisy gels. Results are consistent across replicates.
# Usage: source("260331_loading-quantification_wt-4r-supps_350mM-KGlut.R")
# Output: Bar chart of MCM loading (% WT) for WT, ORC4R, and sofa suppressors
# at 350 mM KGlut, saved as PDF.

# ==============================================================================
# Configuration
# ==============================================================================
library(readxl)
library(tidyverse)

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

EXPERIMENT_DIRECTORY <- "/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication"
INPUT_FILENAME <- "260331_aggregate-analysis_load_wt-4r-supps_350mM-KGlut.xlsx"
INPUT_FILEPATH <- file.path(EXPERIMENT_DIRECTORY, "analysis", INPUT_FILENAME)
OUTPUT_DIRECTORY <- file.path(EXPERIMENT_DIRECTORY, "analysis")

SHEET_NAME <- "Sheet1"
EXPECTED_NUMBER_OF_ROWS <- 21
EXPECTED_NUMBER_OF_COLUMNS <- 4
REQUIRED_COLUMNS <- c("orc4", "sofa", "repeat", "Percent Wildtype")

# ==============================================================================
# File validation
# ==============================================================================
if (!dir.exists(OUTPUT_DIRECTORY)) {
    dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)
}
if (!file.exists(INPUT_FILEPATH)) {
    stop("INPUT_FILEPATH does not exist: ", INPUT_FILEPATH)
}
message("Input and output paths validated.")

# ==============================================================================
# Data loading
# ==============================================================================
raw_loading_data <- read_excel(INPUT_FILEPATH, sheet = SHEET_NAME)

stopifnot(
    "Number of rows does not match EXPECTED_NUMBER_OF_ROWS." =
        nrow(raw_loading_data) == EXPECTED_NUMBER_OF_ROWS,
    "Number of columns does not match EXPECTED_NUMBER_OF_COLUMNS." =
        ncol(raw_loading_data) == EXPECTED_NUMBER_OF_COLUMNS,
    "Column names do not match REQUIRED_COLUMNS." =
        identical(colnames(raw_loading_data), REQUIRED_COLUMNS)
)
message("Data loaded and validated: ", nrow(raw_loading_data), " rows x ", ncol(raw_loading_data), " columns.")

# ==============================================================================
# Preprocessing
# ==============================================================================
LABEL_FACTOR_ORDER <- c("WT", "ORC4R", "+1sofa", "+3sofa", "+4sofa", "+5sofa", "+6sofa")

loading_data <- raw_loading_data %>%
    mutate(label = case_when(
        orc4 == "WT" & sofa == "none" ~ "WT",
        orc4 == "RA" & sofa == "none" ~ "ORC4R",
        orc4 == "RA" & sofa == "orc1" ~ "+1sofa",
        orc4 == "RA" & sofa == "orc3" ~ "+3sofa",
        orc4 == "RA" & sofa == "orc4" ~ "+4sofa",
        orc4 == "RA" & sofa == "orc5" ~ "+5sofa",
        orc4 == "RA" & sofa == "orc6" ~ "+6sofa"
    )) %>%
    mutate(label = factor(label, levels = LABEL_FACTOR_ORDER, ordered = TRUE))

stopifnot(
    "NA labels found after case_when mapping." =
        sum(is.na(loading_data$label)) == 0
)
message("Labels mapped and factor order applied.")

# ==============================================================================
# Summary statistics
# ==============================================================================
summary_loading_data <- loading_data %>%
    group_by(label) %>%
    summarise(
        mean_percent_wildtype = mean(`Percent Wildtype`, na.rm = TRUE),
        sd_percent_wildtype = sd(`Percent Wildtype`, na.rm = TRUE),
        replicate_count = n(),
        .groups = "drop"
    )
message("Summary statistics computed.")

# ==============================================================================
# Plot
# ==============================================================================
loading_bar_chart <- ggplot(summary_loading_data, aes(x = label, y = mean_percent_wildtype, fill = label)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.4) +
    geom_errorbar(
        aes(
            ymin = mean_percent_wildtype - sd_percent_wildtype,
            ymax = mean_percent_wildtype + sd_percent_wildtype
        ),
        width = 0.25, linewidth = 0.6
    ) +
    geom_jitter(
        data = loading_data,
        aes(x = label, y = `Percent Wildtype`),
        width = 0.15, size = 2, shape = 21, fill = "grey30", color = "black", stroke = 0.5,
        inherit.aes = FALSE
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Sample",
        y = "MCM Loading (% WT)",
        title = "MCM Loading at 350 mM KGlut"
    ) +
    theme_classic(base_size = 13) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(face = "bold")
    )
message("Plot constructed.")

# ==============================================================================
# Output
# ==============================================================================
plot_output_filepath <- file.path(OUTPUT_DIRECTORY, "260331_loading-bar-chart_wt-4r-supps_350mM-KGlut.pdf")
summary_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "260331_loading-summary_wt-4r-supps_350mM-KGlut.csv")

if (!file.exists(plot_output_filepath) || OVERWRITE_PLOTS) {
    ggsave(plot_output_filepath, loading_bar_chart, width = 8, height = 5)
    message("Saved plot: ", basename(plot_output_filepath))
} else {
    message("Skipped plot (already exists): ", basename(plot_output_filepath))
}

if (!file.exists(summary_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(summary_loading_data, summary_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(summary_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(summary_csv_output_filepath))
}

message("Script complete.")
