# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., readxl::read_excel).
# Date created: 2026-03-31
# Data produced by analyzing tiff files using ImageJ. Manual gel processing
# due to noisy gels. Results are consistent across replicates.
# Usage: source("260331_loading-quantification_wt-4r-supps_350mM-KGlut.R")
# Output: Bar chart of MCM loading (% WT) for WT, ORC4R, and sofa suppressors
# at 350 mM KGlut, saved as PDF.
#
# Experiment context:
#   Loading assay measures how much MCM helicase is loaded onto chromatin.
#   ORC4R = ORC4 subunit with an R-to-A mutation that impairs loading.
#   sofa = "suppressor of four R/A" - second-site mutations in other ORC
#   subunits (orc1, orc3, orc4, orc5, orc6) that partially rescue ORC4R.
#   This dataset is the 350 mM KGlut salt condition only. A separate script
#   (quantification_kgluttitr_wt-4r-ps.R) handles the multi-salt titration.
#
# Data layout (Sheet1, 21 rows = 7 conditions x 3 biological replicates):
#   orc4: WT or RA (the ORC4R mutant)
#   sofa: none, orc1, orc3, orc4, orc5, orc6
#   repeat: biological replicate number (1, 2, 3)
#   Percent Wildtype: already normalized to WT = 100 per replicate in Excel
#   Sheet2 contains source image filenames only - not used here.

MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")
if (nchar(MC_DROPBOX_PATH) == 0) {
  stop(
    "MC_DROPBOX_PATH not defined in bash environment.\n",
    "Set manually via command line if necessary."
  )

}

# ==============================================================================
# Configuration
# ==============================================================================
library(readxl)
library(tidyverse)

OVERWRITE_PLOTS <- TRUE
OVERWRITE_CSVS <- TRUE

EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication/analysis"
INPUT_FILENAME <- "260331_aggregate-analysis_load_wt-4r-supps_350mM-KGlut.xlsx"
INPUT_FILEPATH <- file.path(MC_DROPBOX_PATH, EXPERIMENT_DIRECTORY, INPUT_FILENAME)
OUTPUT_DIRECTORY <- "~/data/loading_analysis"
SHEET_NAME <- "Sheet1"
EXPECTED_NUMBER_OF_ROWS <- 21
EXPECTED_NUMBER_OF_COLUMNS <- 4
REQUIRED_COLUMNS <- c("orc4", "sofa", "repeat", "Percent Wildtype")

# Centralized plot configuration. Consumed by ggplot calls in the plot section.
# Any parameter set to NULL means "use ggplot default."
# fill_colors: +1sofa and +4sofa colors are swapped relative to default Set1
# so +4sofa keeps its color in the companion faceted figure.
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
        width = 5.2,
        height = 4.5
    )
)

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

# Structural validation against original Excel column names.
stopifnot(
    "Number of rows does not match EXPECTED_NUMBER_OF_ROWS." =
        nrow(raw_loading_data) == EXPECTED_NUMBER_OF_ROWS,
    "Number of columns does not match EXPECTED_NUMBER_OF_COLUMNS." =
        ncol(raw_loading_data) == EXPECTED_NUMBER_OF_COLUMNS,
    "Column names do not match REQUIRED_COLUMNS." =
        identical(colnames(raw_loading_data), REQUIRED_COLUMNS)
)

# Rename to snake_case immediately after structural validation.
# `repeat` is an R reserved word; `Percent Wildtype` contains a space.
raw_loading_data <- dplyr::rename(
    raw_loading_data,
    replicate = `repeat`,
    percent_wildtype = `Percent Wildtype`
)

# Value validation using renamed columns.
stopifnot(
    "NA values found in percent_wildtype column." =
        !anyNA(raw_loading_data[["percent_wildtype"]]),
    "Negative values found in percent_wildtype column." =
        all(raw_loading_data[["percent_wildtype"]] >= 0),
    "WT rows are not exactly 100 in percent_wildtype column." =
        all(raw_loading_data[["percent_wildtype"]][raw_loading_data[["orc4"]] == "WT"] == 100)
)

message("Data loaded and validated: ", nrow(raw_loading_data), " rows x ", ncol(raw_loading_data), " columns.")

# ==============================================================================
# Preprocessing
# ==============================================================================
# Define custom orderings. Named list pattern shared with Script 2.
# suppressor and kglut keys are no-ops in Script 1 (columns don't exist)
# but are included for cross-script consistency.
factor_order <- list(
    "suppressor" = c("None", "1EK", "3PL", "4PS", "5EK"),
    "kglut" = c("250", "300", "350"),
    "label" = c("WT", "ORC4R", "+1sofa", "+3sofa", "+4sofa", "+5sofa", "+6sofa")
)

# Map orc4/sofa column pairs to display labels.
# WT + none = wildtype control. RA + none = ORC4R mutant alone.
# RA + orc<N> = ORC4R rescued by suppressor in subunit N.
loading_data <- raw_loading_data %>%
    mutate(label = case_when(
        orc4 == "WT" & sofa == "none" ~ "WT",
        orc4 == "RA" & sofa == "none" ~ "ORC4R",
        orc4 == "RA" & sofa == "orc1" ~ "+1sofa",
        orc4 == "RA" & sofa == "orc3" ~ "+3sofa",
        orc4 == "RA" & sofa == "orc4" ~ "+4sofa",
        orc4 == "RA" & sofa == "orc5" ~ "+5sofa",
        orc4 == "RA" & sofa == "orc6" ~ "+6sofa"
    ))
# If any row fails to match a case_when branch, it gets NA - catch that here.
stopifnot(
    "NA labels found after case_when mapping." =
        sum(is.na(loading_data$label)) == 0
)

# Apply factor ordering. Column-existence guard skips keys not in data
# (suppressor, kglut are no-ops in Script 1).
for (col_name in names(factor_order)) {
    if (col_name %in% names(loading_data)) {
        loading_data[[col_name]] <- factor(
            loading_data[[col_name]],
            levels = factor_order[[col_name]],
            ordered = FALSE
        )
    }
}

message("Factor ordering applied.")
message("Labels mapped and factor order applied.")


# Percent difference: |A - B| / ((A + B) / 2) * 100
# Computed per-replicate so we can derive mean +/- SD downstream.
# Paired within replicate: each condition compared to the WT and ORC4R
# values from the same experiment.
message("Percent difference, percent change, and fold change columns computed.")
loading_data <- loading_data %>%
    group_by(replicate) %>%
    mutate(
        percent_difference_from_wildtype = abs(percent_wildtype - percent_wildtype[label == "WT"]) /
            ((percent_wildtype + percent_wildtype[label == "WT"]) / 2) * 100,
        percent_difference_from_orc4r = abs(percent_wildtype - percent_wildtype[label == "ORC4R"]) /
            ((percent_wildtype + percent_wildtype[label == "ORC4R"]) / 2) * 100,
        percent_change_from_wildtype = (percent_wildtype - percent_wildtype[label == "WT"]) /
            percent_wildtype[label == "WT"] * 100,
        percent_change_from_orc4r = (percent_wildtype - percent_wildtype[label == "ORC4R"]) /
            percent_wildtype[label == "ORC4R"] * 100,
        fold_change_from_orc4r = percent_wildtype / percent_wildtype[label == "ORC4R"]
    ) %>%
    ungroup()

# ==============================================================================
# Summary statistics
# ==============================================================================
summary_loading_data <- loading_data %>%
    group_by(label) %>%
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
        replicate_count = n(),
        .groups = "drop"
    )
message("Summary statistics computed.")

# ==============================================================================
# Plot
# ==============================================================================

loading_bar_chart <- ggplot(summary_loading_data, aes(x = label, y = mean_percent_wildtype, fill = label)) +
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
        data = loading_data,
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
        title = "MCM Loading at 350 mM KGlut"
    ) +
    theme_classic(
        base_size = PLOT_CONFIG$theme$base_size,
        base_family = PLOT_CONFIG$theme$base_family
    ) +
    theme(
        legend.position = PLOT_CONFIG$theme$legend_position,
        axis.text.x = element_text(face = "bold")
    )
message("Plot constructed.")

# ==============================================================================
# Output
# ==============================================================================
plot_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_bar-chart.pdf")
summary_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_summary.csv")
full_data_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_full-data.csv")

if (!file.exists(plot_output_filepath) || OVERWRITE_PLOTS) {
    ggsave(
        plot_output_filepath,
        loading_bar_chart,
        device = PLOT_CONFIG$output$device,
        width = PLOT_CONFIG$output$width,
        height = PLOT_CONFIG$output$height
    )
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

if (!file.exists(full_data_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(loading_data, full_data_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(full_data_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(full_data_csv_output_filepath))
}

message("Summary statistics:")
print(as.data.frame(summary_loading_data))

message("Script complete.")
