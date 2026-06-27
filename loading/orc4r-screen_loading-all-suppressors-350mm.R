# Namespace convention: tidyverse verbs used unqualified via library(tidyverse).
# All other packages namespaced (e.g., readxl::read_excel).
# Date created: 2026-03-31
# Data produced by analyzing tiff files using ImageJ. Manual gel processing
# due to noisy gels. Results are consistent across replicates.
# Usage: source("orc4r-screen_loading-all-suppressors-350mm.R")
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
#   5. Input: single Excel file, Sheet1, 21 rows = 7 conditions x 3 replicates.
#      Columns: orc4, sofa, repeat, Percent Wildtype.
#      Data is pre-normalized (WT = 100 per replicate in Excel).
# Output: Bar chart of MCM loading (% WT) for WT, ORC4R, and sofa suppressors
# at 350 mM KGlut, saved as PDF. Plus statistical-analysis CSVs (Thread 2a):
# WT-contrast paired t-tests, (suppressor - ORC4R) effect sizes, assumption
# diagnostics, and a mixed-model consistency diagnostic.
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

# ==============================================================================
# GIT STATE REFERENCE (manual-fill at deposit time; no runtime git calls)
# ==============================================================================
# Commit hash:    ____________________________________________
# Branch:         ____________________________________________
# Tag / release:  ____________________________________________
# Snapshot date:  ____________________________________________
# Repository URL: ____________________________________________
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
        "This script must be run via source(\"orc4r-screen_loading-all-suppressors-350mm.R\").\n",
        "Rscript and interactive invocation are not supported ",
        "(no script path is available to resolve data locations)."
    )
}
SCRIPT_DIRECTORY <- dirname(normalizePath(script_path_under_source))

# MC_DROPBOX_PATH is the original (Dropbox) data home; may be unset in a
# Zenodo deposit where the input Excel sits alongside this script.
MC_DROPBOX_PATH <- Sys.getenv("MC_DROPBOX_PATH")

# ==============================================================================
# Configuration
# ==============================================================================
library(readxl)
library(tidyverse)
library(nlme)  # Thread 2a / C9: mixed-model consistency diagnostic only.

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

EXPERIMENT_DIRECTORY <- "Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication/analysis"
INPUT_FILENAME <- "260331_aggregate-analysis_load_wt-4r-supps_350mM-KGlut.xlsx"

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

# Row-count assertion: each label should appear exactly 3 times (one per replicate).
label_counts <- loading_data %>% count(label)
bad_labels <- label_counts %>% filter(n != 3)
if (nrow(bad_labels) > 0) {
    stop("Unexpected row counts per label:\n",
         paste(capture.output(print(as.data.frame(bad_labels))), collapse = "\n"))
}

# WT uniqueness assertion: exactly one WT per replicate.
wt_counts <- loading_data %>% filter(label == "WT") %>% count(replicate)
bad_wt <- wt_counts %>% filter(n != 1)
if (nrow(bad_wt) > 0) {
    stop("Expected exactly 1 WT per replicate. Offending replicates:\n",
         paste(capture.output(print(as.data.frame(bad_wt))), collapse = "\n"))
}

# ORC4R uniqueness assertion: exactly one ORC4R per replicate.
orc4r_counts <- loading_data %>% filter(label == "ORC4R") %>% count(replicate)
bad_orc4r <- orc4r_counts %>% filter(n != 1)
if (nrow(bad_orc4r) > 0) {
    stop("Expected exactly 1 ORC4R per replicate. Offending replicates:\n",
         paste(capture.output(print(as.data.frame(bad_orc4r))), collapse = "\n"))
}

message("Row-count and uniqueness assertions passed.")

# Percent difference / change / fold change, paired within replicate: each
# condition compared to the WT and ORC4R values from the same experiment.
# C4: this derivation MUST run before the Inf/NaN guard below, which inspects
# the columns it creates.
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
message("Percent difference, percent change, and fold change columns computed.")

# Guard against division by zero in derived calculations.
# C4: runs AFTER the mutate above so it inspects existing columns (previously
# it ran first and passed vacuously against NULL columns).
derived_cols <- c("fold_change_from_orc4r", "percent_change_from_wildtype",
                  "percent_change_from_orc4r")
for (col in derived_cols) {
    stopifnot(
        !any(is.infinite(loading_data[[col]])),
        !any(is.nan(loading_data[[col]]))
    )
}
message("Inf/NaN guard passed on derived calculations.")

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

stopifnot(
    "Not all summary groups have 3 replicates." =
        all(summary_loading_data$replicate_count == 3)
)

# ==============================================================================
# Statistical analysis (Thread 2a: C6, C7, C8, C9)
# ------------------------------------------------------------------------------
# ALL statistical COMPUTATION + augmented-frame construction lives in THIS
# section, which is UPSTREAM of the font/Cairo-dependent Plot section. If
# plotting later fails, the stats CSVs (written in Output) are still produced.
# (Thread 1 carry-forward #2.)
#
# DESIGN FACT confirmed against the real 350 mM data (Thread 2a verification):
#   - block variable = `replicate`, values 1/2/3.
#   - 7 labels (WT, ORC4R, +1/+3/+4/+5/+6sofa), n = 3 each, 21 rows total.
#   - pairing is COMPLETE: WT, ORC4R and every suppressor co-occur in every
#     replicate, so all 6 WT-contrasts use the full n = 3 (no incomplete block).
#   - percent_wildtype pins WT = 100 in every replicate -> WT within-group
#     variance is exactly 0 (relevant to C6 and C9 below).
# ==============================================================================

# ------------------------------------------------------------------------------
# C6: Assumption diagnostics for the 350 mM loading contrasts.
# CLAIM SUPPORTED: none directly -- these are assumption checks feeding C7/C8.
# I4 n=3 CAVEAT: with 3 biological replicates per condition, normality CANNOT
# be meaningfully tested (Shapiro-Wilk at n=3 has almost no power and is NOT a
# license to assume normality). Justification for the parametric paired test
# rests on the PAIRED/BLOCKED design (WT, ORC4R and every suppressor co-run in
# each replicate) plus the assay's established track record -- NOT on the
# Shapiro number below passing. Reported for transparency only.
# DATA FACT: percent_wildtype pins WT = 100 in every replicate, so WT has
# exactly zero within-group variance. Standard variance-homogeneity tests
# (Bartlett/Levene) are undefined with a zero-variance group; we therefore
# report per-group variances directly and lean on the within-replicate paired
# differences (which is what the t-test actually consumes). QQ is reported as
# the numeric Shapiro statistic rather than a font-gated plot, so diagnostics
# survive a font-less bare run.
# ------------------------------------------------------------------------------

# Membership-checked guard (carry-forward #1): each referenced column must
# exist BEFORE its values are tested, else the guard passes vacuously.
stopifnot(
    "label column missing for diagnostics." = "label" %in% names(loading_data),
    "replicate column missing for diagnostics." = "replicate" %in% names(loading_data),
    "percent_wildtype column missing for diagnostics." = "percent_wildtype" %in% names(loading_data)
)

loading_diagnostics_per_label <- loading_data %>%
    group_by(label) %>%
    summarise(
        n_replicates = n(),
        mean_percent_wildtype = mean(percent_wildtype),
        sd_percent_wildtype = sd(percent_wildtype),
        variance_percent_wildtype = var(percent_wildtype),
        .groups = "drop"
    )

# Pooled within-group-centered residuals (each value minus its label mean).
loading_residuals <- loading_data %>%
    group_by(label) %>%
    mutate(within_label_residual = percent_wildtype - mean(percent_wildtype)) %>%
    ungroup()
pooled_residual_shapiro <- shapiro.test(loading_residuals$within_label_residual)

loading_diagnostics_normality <- data.frame(
    diagnostic = "shapiro_wilk_pooled_within_label_residuals",
    shapiro_w = unname(pooled_residual_shapiro$statistic),
    shapiro_p_value = pooled_residual_shapiro$p.value,
    n_residuals = nrow(loading_residuals),
    caveat = "n=3 per group: normality untestable; reported for transparency only.",
    stringsAsFactors = FALSE
)
message("C6 diagnostics computed (per-label variances + pooled-residual Shapiro).")

# ------------------------------------------------------------------------------
# C7: WT-contrast tests -- the "non-restoration-to-WT" claim.
# For each of ORC4R and the 5 suppressors: a PAIRED t-test (block = replicate)
# against WT is the REPORTED result. A paired Wilcoxon signed-rank test is run
# as a DECORATIVE companion only (I4): at n=3 paired, the two-sided signed-rank
# p-value FLOORS at ~0.25, so it can never reach significance; shown for
# transparency, never leaned on.
# HOLM FAMILY (I2): the SIX WT-contrasts within the SINGLE 350 mM KGlut
# condition form ONE family; Holm is applied across exactly these six and
# NOTHING else (never across salts/timepoints -- there are none in this script).
# I4 n=3 CAVEAT: normality of the paired differences is untestable at n=3;
# justification is the paired/blocked design + assay track record (see C6).
# ------------------------------------------------------------------------------

stopifnot(
    "label column missing for WT contrasts." = "label" %in% names(loading_data),
    "percent_wildtype column missing for WT contrasts." = "percent_wildtype" %in% names(loading_data),
    "replicate column missing for WT contrasts." = "replicate" %in% names(loading_data)
)

# WT reference values, ordered by replicate so pairing aligns by block.
wildtype_values_by_replicate <- loading_data %>%
    filter(label == "WT") %>%
    arrange(replicate) %>%
    pull(percent_wildtype)
stopifnot(
    "WT does not have exactly 3 paired values." = length(wildtype_values_by_replicate) == 3
)

wt_contrast_comparison_labels <- c("ORC4R", "+1sofa", "+3sofa", "+4sofa", "+5sofa", "+6sofa")

wt_contrast_label <- character(0)
wt_contrast_n_pairs <- integer(0)
wt_contrast_estimate_wt_minus_group <- numeric(0)
wt_contrast_ci_lower <- numeric(0)
wt_contrast_ci_upper <- numeric(0)
wt_contrast_t_statistic <- numeric(0)
wt_contrast_df <- numeric(0)
wt_contrast_raw_p_value <- numeric(0)
wt_contrast_wilcoxon_p_value <- numeric(0)

for (comparison_label in wt_contrast_comparison_labels) {
    group_values_by_replicate <- loading_data %>%
        filter(label == comparison_label) %>%
        arrange(replicate) %>%
        pull(percent_wildtype)

    # Pairing-completeness guard: each comparison group must have exactly one
    # value per replicate (complete blocks), or the paired test is invalid.
    stopifnot(
        "Comparison group lacks 3 paired values (incomplete block)." =
            length(group_values_by_replicate) == 3
    )

    paired_t <- t.test(wildtype_values_by_replicate, group_values_by_replicate, paired = TRUE)
    paired_wilcox <- suppressWarnings(
        wilcox.test(wildtype_values_by_replicate, group_values_by_replicate, paired = TRUE)
    )

    wt_contrast_label <- c(wt_contrast_label, comparison_label)
    wt_contrast_n_pairs <- c(wt_contrast_n_pairs, length(group_values_by_replicate))
    wt_contrast_estimate_wt_minus_group <- c(wt_contrast_estimate_wt_minus_group, unname(paired_t$estimate))
    wt_contrast_ci_lower <- c(wt_contrast_ci_lower, paired_t$conf.int[1])
    wt_contrast_ci_upper <- c(wt_contrast_ci_upper, paired_t$conf.int[2])
    wt_contrast_t_statistic <- c(wt_contrast_t_statistic, unname(paired_t$statistic))
    wt_contrast_df <- c(wt_contrast_df, unname(paired_t$parameter))
    wt_contrast_raw_p_value <- c(wt_contrast_raw_p_value, paired_t$p.value)
    wt_contrast_wilcoxon_p_value <- c(wt_contrast_wilcoxon_p_value, paired_wilcox$p.value)
}

# Holm correction WITHIN the six-member WT-contrast family (I2).
wt_contrast_holm_p_value <- p.adjust(wt_contrast_raw_p_value, method = "holm")

wt_contrast_results <- data.frame(
    comparison = paste0("WT_vs_", wt_contrast_label),
    group_label = wt_contrast_label,
    claim = "non-restoration-to-WT (rejectable WT-vs-group difference)",
    test = "paired t-test, block = replicate",
    holm_family = "6 WT-contrasts within 350mM KGlut condition",
    n_pairs = wt_contrast_n_pairs,
    estimate_wt_minus_group = wt_contrast_estimate_wt_minus_group,
    ci_lower_wt_minus_group = wt_contrast_ci_lower,
    ci_upper_wt_minus_group = wt_contrast_ci_upper,
    t_statistic = wt_contrast_t_statistic,
    df = wt_contrast_df,
    raw_p_value = wt_contrast_raw_p_value,
    holm_adjusted_p_value = wt_contrast_holm_p_value,
    wilcoxon_p_value_decorative = wt_contrast_wilcoxon_p_value,
    stringsAsFactors = FALSE
)
message("C7 WT-contrast paired t-tests + decorative Wilcoxon + Holm computed.")
print(wt_contrast_results)

# ------------------------------------------------------------------------------
# C8: (suppressor - ORC4R) ORDINARY EFFECT SIZE with 95% CI.
# NARRATIVE FIREWALL (handoff F6): this is a PLAIN paired effect size measuring
# how far each suppressor sits from ORC4R on the loading axis. It is NOT, and
# must not be read as, the ATPase "4R-floor-proximity" argument. In LOADING the
# reducer suppressors (+3sofa/Orc3, +5sofa/Orc5, +6sofa/Orc6) fall BELOW ORC4R
# BY DESIGN, so a negative (suppressor - ORC4R) estimate is an EXPECTED
# biological result, not evidence of a floor. The ATPase floor claim rests on a
# shared catalytic mutation and lives ONLY in the ATPase thread (Thread 3).
# Do not conflate the two narratives.
# I4 n=3 CAVEAT: df = 2, so these CIs are wide; reported as effect sizes for
# magnitude/direction, not as significance tests.
# ------------------------------------------------------------------------------

stopifnot(
    "label column missing for effect-size-vs-ORC4R." = "label" %in% names(loading_data),
    "percent_wildtype column missing for effect-size-vs-ORC4R." = "percent_wildtype" %in% names(loading_data),
    "replicate column missing for effect-size-vs-ORC4R." = "replicate" %in% names(loading_data)
)

orc4r_values_by_replicate <- loading_data %>%
    filter(label == "ORC4R") %>%
    arrange(replicate) %>%
    pull(percent_wildtype)
stopifnot(
    "ORC4R does not have exactly 3 paired values." = length(orc4r_values_by_replicate) == 3
)

effect_size_suppressor_labels <- c("+1sofa", "+3sofa", "+4sofa", "+5sofa", "+6sofa")

effect_size_label <- character(0)
effect_size_n_pairs <- integer(0)
effect_size_estimate_group_minus_orc4r <- numeric(0)
effect_size_ci_lower <- numeric(0)
effect_size_ci_upper <- numeric(0)

for (suppressor_label in effect_size_suppressor_labels) {
    suppressor_values_by_replicate <- loading_data %>%
        filter(label == suppressor_label) %>%
        arrange(replicate) %>%
        pull(percent_wildtype)
    stopifnot(
        "Suppressor lacks 3 paired values (incomplete block)." =
            length(suppressor_values_by_replicate) == 3
    )

    effect_t <- t.test(suppressor_values_by_replicate, orc4r_values_by_replicate, paired = TRUE)

    effect_size_label <- c(effect_size_label, suppressor_label)
    effect_size_n_pairs <- c(effect_size_n_pairs, length(suppressor_values_by_replicate))
    effect_size_estimate_group_minus_orc4r <- c(effect_size_estimate_group_minus_orc4r, unname(effect_t$estimate))
    effect_size_ci_lower <- c(effect_size_ci_lower, effect_t$conf.int[1])
    effect_size_ci_upper <- c(effect_size_ci_upper, effect_t$conf.int[2])
}

effect_size_from_orc4r_results <- data.frame(
    suppressor = effect_size_label,
    quantity = "ORDINARY effect size (suppressor - ORC4R); NOT an ATPase floor-proximity claim",
    n_pairs = effect_size_n_pairs,
    estimate_group_minus_orc4r = effect_size_estimate_group_minus_orc4r,
    ci_lower = effect_size_ci_lower,
    ci_upper = effect_size_ci_upper,
    stringsAsFactors = FALSE
)
message("C8 (suppressor - ORC4R) ordinary effect sizes + 95% CI computed.")
print(effect_size_from_orc4r_results)

# ------------------------------------------------------------------------------
# C9: PRE-COMMITTED mixed-model consistency DIAGNOSTIC (I5).
# The REPORTED result for every WT-contrast is the paired t-test in C7, fixed
# in advance. This nlme random-intercept model (percent_wildtype ~ label,
# random = ~1|replicate) is a consistency CHECK ONLY. If it disagrees
# MATERIALLY with C7, that triggers INVESTIGATION -- it NEVER licenses swapping
# in the friendlier p-value.
# DATA FACT / why this is only a diagnostic: WT is pinned to 100 with zero
# within-group variance, so residuals are heteroscedastic and lme's single
# residual-variance assumption is violated; the within-replicate paired t-test
# (which differences WT out) is the honest frame. The fit is wrapped so a
# convergence failure FLAGS the diagnostic rather than crashing the script
# (a non-converging diagnostic is itself an investigation flag, per I5).
# WT is the factor reference level, so each label coefficient = (group - WT).
# ------------------------------------------------------------------------------

stopifnot(
    "label column missing for mixed-model diagnostic." = "label" %in% names(loading_data),
    "replicate column missing for mixed-model diagnostic." = "replicate" %in% names(loading_data),
    "percent_wildtype column missing for mixed-model diagnostic." = "percent_wildtype" %in% names(loading_data)
)

mixed_model_diagnostic_fit <- tryCatch(
    nlme::lme(percent_wildtype ~ label, random = ~ 1 | replicate, data = loading_data),
    error = function(e) e
)

if (inherits(mixed_model_diagnostic_fit, "error")) {
    mixed_model_diagnostic_results <- data.frame(
        term = NA_character_,
        estimate_group_minus_wt = NA_real_,
        std_error = NA_real_,
        t_value = NA_real_,
        model_p_value = NA_real_,
        converged = FALSE,
        note = paste0("lme did not converge: ", conditionMessage(mixed_model_diagnostic_fit),
                      " -- C7 paired t-test stands as reported (I5)."),
        stringsAsFactors = FALSE
    )
    message("C9 mixed-model diagnostic DID NOT CONVERGE; flagged. C7 stands (I5).")
} else {
    mixed_model_fixed_effects <- summary(mixed_model_diagnostic_fit)$tTable
    # Drop the intercept row (= WT itself); keep label coefficients.
    label_coefficient_rows <- rownames(mixed_model_fixed_effects) != "(Intercept)"
    mixed_model_diagnostic_results <- data.frame(
        term = rownames(mixed_model_fixed_effects)[label_coefficient_rows],
        estimate_group_minus_wt = mixed_model_fixed_effects[label_coefficient_rows, "Value"],
        std_error = mixed_model_fixed_effects[label_coefficient_rows, "Std.Error"],
        t_value = mixed_model_fixed_effects[label_coefficient_rows, "t-value"],
        model_p_value = mixed_model_fixed_effects[label_coefficient_rows, "p-value"],
        converged = TRUE,
        note = "DIAGNOSTIC ONLY (I5); reported result is the C7 paired t-test.",
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    message("C9 mixed-model diagnostic converged; flagged as consistency check (I5).")
}
print(mixed_model_diagnostic_results)

# ------------------------------------------------------------------------------
# C10a (compute half): build augmented summary + full-data frames.
# Column lineage: every WT-contrast/effect column derives from percent_wildtype
# (the per-replicate-normalized, WT=100 column). CSV writes happen in Output.
# ------------------------------------------------------------------------------
summary_loading_data_augmented <- summary_loading_data %>%
    mutate(label = as.character(label)) %>%
    left_join(
        data.frame(
            label = wt_contrast_results$group_label,
            wt_contrast_test = wt_contrast_results$test,
            wt_contrast_n_pairs = wt_contrast_results$n_pairs,
            wt_contrast_estimate_wt_minus_group = wt_contrast_results$estimate_wt_minus_group,
            wt_contrast_ci_lower = wt_contrast_results$ci_lower_wt_minus_group,
            wt_contrast_ci_upper = wt_contrast_results$ci_upper_wt_minus_group,
            wt_contrast_raw_p_value = wt_contrast_results$raw_p_value,
            wt_contrast_holm_adjusted_p_value = wt_contrast_results$holm_adjusted_p_value,
            wt_contrast_wilcoxon_p_value_decorative = wt_contrast_results$wilcoxon_p_value_decorative,
            stringsAsFactors = FALSE
        ),
        by = "label"
    ) %>%
    left_join(
        data.frame(
            label = effect_size_from_orc4r_results$suppressor,
            effect_estimate_group_minus_orc4r = effect_size_from_orc4r_results$estimate_group_minus_orc4r,
            effect_ci_lower = effect_size_from_orc4r_results$ci_lower,
            effect_ci_upper = effect_size_from_orc4r_results$ci_upper,
            stringsAsFactors = FALSE
        ),
        by = "label"
    )

loading_data_augmented <- loading_data %>%
    mutate(label = as.character(label)) %>%
    left_join(
        data.frame(
            label = wt_contrast_results$group_label,
            wt_contrast_raw_p_value = wt_contrast_results$raw_p_value,
            wt_contrast_holm_adjusted_p_value = wt_contrast_results$holm_adjusted_p_value,
            stringsAsFactors = FALSE
        ),
        by = "label"
    )

message("C10a augmented summary + full-data frames constructed.")

# ==============================================================================
# Plot
# ==============================================================================

# C10b (display compute): exact Holm-adjusted p annotations for the WT-contrast
# family. EXACT p-values, never stars (style contract). Display-only; the
# numbers come from C7 upstream (font-gate-independent).
wt_contrast_plot_annotations <- summary_loading_data %>%
    mutate(label_character = as.character(label)) %>%
    inner_join(
        data.frame(
            label_character = wt_contrast_results$group_label,
            holm_adjusted_p_value = wt_contrast_results$holm_adjusted_p_value,
            stringsAsFactors = FALSE
        ),
        by = "label_character"
    ) %>%
    mutate(
        annotation_y = mean_percent_wildtype + sd_percent_wildtype + 6,
        annotation_text = paste0(
            "p = ", formatC(holm_adjusted_p_value, format = "g", digits = 3), "\n(Holm vs WT)"
        )
    )

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
    geom_text(
        data = wt_contrast_plot_annotations,
        aes(x = label, y = annotation_y, label = annotation_text),
        inherit.aes = FALSE,
        size = 3,
        family = PLOT_CONFIG$theme$base_family,
        lineheight = 0.9
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
        title = "MCM Loading at 350 mM KGlut",
        caption = paste0(
            "Exact Holm-adjusted p-values vs WT (paired t-test, block = replicate; ",
            "Holm family = 6 WT-contrasts within 350 mM). n = 3 per condition."
        )
    ) +
    theme_classic(
        base_size = PLOT_CONFIG$theme$base_size,
        base_family = PLOT_CONFIG$theme$base_family
    ) +
    theme(
        legend.position = PLOT_CONFIG$theme$legend_position,
        axis.text.x = element_text(face = "bold"),
        plot.caption = element_text(size = 7, hjust = 0)
    )
message("Plot constructed (with exact Holm-adjusted p-value annotations).")

# ==============================================================================
# Output
# ==============================================================================
plot_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_bar-chart.pdf")
summary_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_summary.csv")
full_data_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_full-data.csv")
wt_contrast_results_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_wt-contrasts.csv")
effect_size_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_effect-size-vs-orc4r.csv")
diagnostics_per_label_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_diagnostics-per-label.csv")
diagnostics_normality_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_diagnostics-normality.csv")
mixed_model_diagnostic_csv_output_filepath <- file.path(OUTPUT_DIRECTORY, "loading_350mm-all-suppressors_mixed-model-diagnostic.csv")

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
    write.csv(summary_loading_data_augmented, summary_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(summary_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(summary_csv_output_filepath))
}

if (!file.exists(full_data_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(loading_data_augmented, full_data_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(full_data_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(full_data_csv_output_filepath))
}

if (!file.exists(wt_contrast_results_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(wt_contrast_results, wt_contrast_results_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(wt_contrast_results_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(wt_contrast_results_csv_output_filepath))
}

if (!file.exists(effect_size_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(effect_size_from_orc4r_results, effect_size_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(effect_size_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(effect_size_csv_output_filepath))
}

if (!file.exists(diagnostics_per_label_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(loading_diagnostics_per_label, diagnostics_per_label_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(diagnostics_per_label_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(diagnostics_per_label_csv_output_filepath))
}

if (!file.exists(diagnostics_normality_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(loading_diagnostics_normality, diagnostics_normality_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(diagnostics_normality_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(diagnostics_normality_csv_output_filepath))
}

if (!file.exists(mixed_model_diagnostic_csv_output_filepath) || OVERWRITE_CSVS) {
    write.csv(mixed_model_diagnostic_results, mixed_model_diagnostic_csv_output_filepath, row.names = FALSE)
    message("Saved CSV: ", basename(mixed_model_diagnostic_csv_output_filepath))
} else {
    message("Skipped CSV (already exists): ", basename(mixed_model_diagnostic_csv_output_filepath))
}

message("Summary statistics (augmented):")
print(as.data.frame(summary_loading_data_augmented))

message("Script complete.")
