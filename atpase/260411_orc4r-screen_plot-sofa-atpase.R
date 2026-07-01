# ==============================================================================
# ATPase Percent Hydrolysis Plot - Standalone
# ==============================================================================
# Date created: 2026-04-11
# Reads processed_data.csv AND atpase_per_timepoint_wt_contrasts.csv produced by
# 260411_orc4r-screen_analyze-sofa-atpase.R, and generates the percent
# hydrolysis timecourse plot WITH exact per-timepoint WT-contrast p-value
# annotations.
# ANTI-DRIFT (C18b): all inferential statistics (p-values) are READ from the
# analyze script's results CSV; this script never recomputes a test. The line
# geometry (means/SD) is recomputed deterministically from processed_data for
# display only.
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
# Script location (C1)
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
MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")

PROCESSED_DATA_FILENAME <- "processed_data.csv"

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

# C18b: per-timepoint WT-contrast results CSV (written by the analyze script).
# This is the SOLE source of plotted p-values (anti-drift).
WT_CONTRASTS_FILENAME <- "atpase_per_timepoint_wt_contrasts.csv"
WT_CONTRASTS_FILEPATH <- file.path(INPUT_DIRECTORY, WT_CONTRASTS_FILENAME)

message("Input file: ", PROCESSED_DATA_FILEPATH)
message("Stats CSV:  ", WT_CONTRASTS_FILEPATH)
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

plotted_samples <- setdiff(unique(processed_data$sample), c("No_ORC", "WT_DNA"))
unknown_samples <- setdiff(plotted_samples, SAMPLE_DISPLAY_ORDER)
if (length(unknown_samples) > 0) {
    stop(
        "Unknown sample names in data: ",
        paste(unknown_samples, collapse = ", "),
        "\nUpdate SAMPLE_DISPLAY_ORDER and SAMPLE_COLORS."
    )
}

# C18b: stats CSV is REQUIRED (anti-drift: we read, never recompute).
if (!file.exists(WT_CONTRASTS_FILEPATH)) {
    stop(
        "Per-timepoint WT-contrasts CSV not found:\n  ", WT_CONTRASTS_FILEPATH,
        "\nRun 260411_orc4r-screen_analyze-sofa-atpase.R first (it computes and ",
        "writes the statistics this plot annotates)."
    )
}
wt_contrasts <- read.csv(WT_CONTRASTS_FILEPATH, stringsAsFactors = FALSE)
REQUIRED_CONTRAST_COLUMNS <- c(
    "timepoint", "comparison", "group_label", "n_pairs",
    "holm_adjusted_p_value", "computable"
)
missing_contrast_columns <- setdiff(REQUIRED_CONTRAST_COLUMNS, colnames(wt_contrasts))
if (length(missing_contrast_columns) > 0) {
    stop(
        "Missing columns in WT-contrasts CSV: ",
        paste(missing_contrast_columns, collapse = ", "),
        "\nRe-run the analyze script (C18a writes these)."
    )
}
message("Loaded ", nrow(wt_contrasts), " WT-contrast rows from stats CSV.")

message("All required columns present. All sample names recognized.")
message("=== DATA LOADING COMPLETE ===")

# ==============================================================================
# SUMMARY STATISTICS (display geometry only -- means/SD for lines + error bars)
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
# ANNOTATION TABLE (C18b: exact p-values READ from the stats CSV, never stars)
# ==============================================================================
message("=== ANNOTATION (exact WT-contrast p-values, from stats CSV) ===")

# WT-vs-{4R + each suppressor}, Holm-adjusted WITHIN each timepoint (I2).
annotation_data <- wt_contrasts
annotation_data$short_label <- sub("_.*", "", annotation_data$group_label)
annotation_data$p_text <- ifelse(
    annotation_data$computable %in% c(TRUE, "TRUE", "True"),
    paste0(
        annotation_data$short_label, " p=",
        sprintf("%.4f", annotation_data$holm_adjusted_p_value),
        " (n=", annotation_data$n_pairs, ")"
    ),
    paste0(annotation_data$short_label, " p=n.c. (n=", annotation_data$n_pairs, ")")
)

# Stack the (up to) six labels above each timepoint cluster.
annotation_data <- annotation_data[
    order(annotation_data$timepoint, annotation_data$group_label),
]
ANNOTATION_Y_TOP <- 0.62
ANNOTATION_Y_STEP <- 0.045
annotation_data$y_pos <- NA_real_
for (current_timepoint in unique(annotation_data$timepoint)) {
    idx <- which(annotation_data$timepoint == current_timepoint)
    annotation_data$y_pos[idx] <- ANNOTATION_Y_TOP - (seq_along(idx) - 1) * ANNOTATION_Y_STEP
}

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
    geom_text(
        data = annotation_data,
        aes(x = timepoint, y = y_pos, label = p_text, color = group_label),
        hjust = 0.5,
        size = 2.3,
        show.legend = FALSE
    ) +
    scale_color_manual(values = SAMPLE_COLORS) +
    coord_cartesian(xlim = c(-8, 104), ylim = c(-0.08, 0.66)) +
    labs(
        title = "ORC ATPase Timecourse",
        subtitle = "Annotations: WT vs sample, paired t-test (block=experiment), Holm-adjusted WITHIN timepoint; exact p, n.c.=not computable",
        x = "Time (min)",
        y = "Percent hydrolysis (background-corrected)",
        color = "Sample",
        caption = paste0(
            "Per-timepoint WT-contrasts read from ", WT_CONTRASTS_FILENAME,
            " (not recomputed). Holm family = WT-contrasts within a single ",
            "timepoint (NEVER across timepoints). WT subset to each sample's ",
            "shared experiments; n shown per label. Pooled suppressor-vs-4R ",
            "contrasts + CIs and the sample:timepoint interaction are in the ",
            "mixed-model CSVs."
        )
    ) +
    theme_classic(base_size = 13) +
    theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        plot.caption = element_text(hjust = 0, size = 6),
        legend.position = "right"
    )

plot_pdf_path <- file.path(OUTPUT_DIRECTORY, "atpase_timecourse_plot.pdf")
if (!file.exists(plot_pdf_path) || OVERWRITE_PLOTS) {
    ggsave(plot_pdf_path, atpase_timecourse_plot, width = 9, height = 6)
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
