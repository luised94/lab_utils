# ==============================================================================
# ATPase Percent Hydrolysis Plot - Standalone
# ==============================================================================
# Date created: 2026-04-11
# Reads processed_data.csv produced by 260411_orc4r-screen_analyze-sofa-atpase.R
# and generates the percent hydrolysis timecourse plot.
# Usage: source("260411_orc4r-screen_plot-sofa-atpase.R")
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
        "This script must be run via source(\"260411_orc4r-screen_plot-sofa-atpase.R\").\n",
        "Rscript and interactive invocation are not supported ",
        "(no script path is available to resolve data locations)."
    )
}
SCRIPT_DIRECTORY <- dirname(normalizePath(script_path_under_source))

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
# MC_DROPBOX_PATH is the original (Dropbox) data home; may be unset in a
# Zenodo deposit where processed_data.csv sits alongside this script.
MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")

PROCESSED_DATA_FILENAME <- "processed_data.csv"

# Path resolution (C1): script-relative (Zenodo co-located) -> MC_DROPBOX_PATH
# (consolidated_analysis, where the analyze script writes it) -> stop().
script_relative_processed_data_filepath <- file.path(
    SCRIPT_DIRECTORY, PROCESSED_DATA_FILENAME
)
if (nchar(MC_DROPBOX_PATH) > 0) {
    dropbox_processed_data_filepath <- file.path(
        MC_DROPBOX_PATH, "Lab", "Experiments", "ATPase",
        "2020_09_03 ATPase Analysis of 4R supps", "consolidated_analysis",
        PROCESSED_DATA_FILENAME
    )
} else {
    dropbox_processed_data_filepath <- NA_character_
}

if (file.exists(script_relative_processed_data_filepath)) {
    PROCESSED_DATA_FILEPATH <- script_relative_processed_data_filepath
} else if (!is.na(dropbox_processed_data_filepath) &&
    file.exists(dropbox_processed_data_filepath)) {
    PROCESSED_DATA_FILEPATH <- dropbox_processed_data_filepath
} else {
    stop(
        "processed_data.csv not found in either supported location:\n",
        "  script-relative: ", script_relative_processed_data_filepath, "\n",
        "  MC_DROPBOX_PATH:  ",
        if (is.na(dropbox_processed_data_filepath)) {
            "<MC_DROPBOX_PATH not set>"
        } else {
            dropbox_processed_data_filepath
        }, "\n",
        "Run 260411_orc4r-screen_analyze-sofa-atpase.R first."
    )
}

INPUT_DIRECTORY <- dirname(PROCESSED_DATA_FILEPATH)
OUTPUT_DIRECTORY <- INPUT_DIRECTORY

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
