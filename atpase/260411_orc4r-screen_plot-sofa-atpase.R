# ==============================================================================
# ATPase Percent Hydrolysis Plot - Standalone
# ==============================================================================
# Date created: 2026-04-11
# Reads processed_data.csv produced by 260411_orc4r-screen_analyze-sofa-atpase.R
# and generates the percent hydrolysis timecourse plot.
# Usage: Source in R REPL or run via Rscript.
# ==============================================================================

library(tidyverse)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
message("=== CONFIGURATION ===")

OVERWRITE_PLOTS <- TRUE

TIMEPOINTS <- c(0, 15, 45, 90)

# ------------------------------------------------------------------------------
# Sample ordering and colors
# ------------------------------------------------------------------------------
# Factor order determines legend and axis ordering.
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

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
# Expects processed_data.csv in the same directory as this script.
SCRIPT_DIRECTORY <- getwd()

INPUT_DIRECTORY <- SCRIPT_DIRECTORY
OUTPUT_DIRECTORY <- SCRIPT_DIRECTORY

PROCESSED_DATA_FILENAME <- "processed_data.csv"
PROCESSED_DATA_FILEPATH <- file.path(INPUT_DIRECTORY, PROCESSED_DATA_FILENAME)

if (!file.exists(PROCESSED_DATA_FILEPATH)) {
    stop(
        "Processed data file not found: ", PROCESSED_DATA_FILEPATH, "\n",
        "Run 260411_orc4r-screen_analyze-sofa-atpase.R first."
    )
}

message("Input file: ", PROCESSED_DATA_FILEPATH)
message("Output directory: ", OUTPUT_DIRECTORY)
message("=== CONFIGURATION COMPLETE ===")

# ==============================================================================
# DATA LOADING AND VALIDATION
# ==============================================================================
message("=== DATA LOADING ===")

processed_data <- read.csv(PROCESSED_DATA_FILEPATH, stringsAsFactors = FALSE)
message("Loaded ", nrow(processed_data), " rows x ", ncol(processed_data), " columns.")

REQUIRED_COLUMNS <- c(
    "adp_intensity", "atp_intensity", "timepoint", "sample",
    "experiment_label", "percent_adp", "percent_adp_corrected"
)
missing_columns <- setdiff(REQUIRED_COLUMNS, colnames(processed_data))
if (length(missing_columns) > 0) {
    stop(
        "Missing columns in processed data: ",
        paste(missing_columns, collapse = ", ")
    )
}

if (anyNA(processed_data$percent_adp_corrected)) {
    stop("NA values found in percent_adp_corrected column.")
}

# Verify all expected samples are present in the data.
plotted_samples <- setdiff(unique(processed_data$sample), c("No_ORC", "WT_DNA"))
unknown_samples <- setdiff(plotted_samples, SAMPLE_DISPLAY_ORDER)
if (length(unknown_samples) > 0) {
    stop(
        "Unknown sample names in data: ",
        paste(unknown_samples, collapse = ", "),
        "\nUpdate SAMPLE_DISPLAY_ORDER and SAMPLE_COLORS."
    )
}

message("All required columns present. All sample names recognized.")
message("=== DATA LOADING COMPLETE ===")

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
print(
    summary_data[summary_data$timepoint == 0, c("sample", "replicate_count")]
)

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
for (output_file in list.files(OUTPUT_DIRECTORY, pattern = "\\.pdf$")) {
    message("  ", output_file)
}
