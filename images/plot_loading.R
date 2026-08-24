# ==============================================================================
# plot_loading.R
#
# Aggregate loading per-lane CSVs across replicate gels, normalize each gel to its
# own reference lane, write the summary and all-against-all matrices the paper
# references, and render a genotype bar figure (finished in Illustrator).
#
# Input is a MANIFEST CSV (not a directory, not a TIFF): one row per replicate gel,
# columns
#     replicate,loading_csv
# where loading_csv is a path to a loading_value_*.csv produced per gel by
# calculate_loading_value.py. Listing the CSVs directly avoids re-deriving the
# per-gel filename, whose millimetre window differs between gels.
#
# Usage (Rscript, headless):
#     Rscript plot_loading.R <manifest.csv> <output_directory>
#
# The reported quantity is a loading ratio, but -- unlike gel shift's self-
# normalized bound_fraction -- it is lane-to-REFERENCE-lane: within each gel, every
# lane's value is divided by the analysis_role == "reference" lane's value, so the
# reference reads 1.0 and every other lane reads relative to it. That per-gel
# normalization is done HERE (matching the loading reference scripts, which compute
# "percent of reference"), because it needs each gel's reference lane, which the
# per-gel calc deliberately left in the CSV rather than dividing away.
#
# Dependencies: tidyverse only. The loading reference scripts also load readxl
# (Excel input; not needed here, we read CSVs) and nlme (a mixed-model diagnostic
# not requested here); both are omitted so this stays a strict subset. extrafont is
# intentionally not loaded: base_family = "Arial" is requested in the theme and R
# falls back cleanly if the font is not registered, so the script runs headless.
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# Per-run analysis choices. These are the only knobs; edit and re-run.
#
# NORMALIZATION_VALUE_COLUMN is the per-lane quantity divided by the WT lane.
# "value_corrected" is the blank-corrected value from the calc; "value" is the raw
# window sum. The reference model wants the corrected value by default.
#
# The normalization reference is the WT lane (it reads 1.0 = 100%), found by derived
# genotype == "WT", not by an analysis_role marker; see the per-gel loop below.
#
# GROUPING: two factor columns, orc4-R267 and suppressor, combined into a genotype
# key exactly as gel shift did, but with NO ATP axis (loading has no ATP contrast),
# so there is one bar per genotype. If a future loading experiment varies different
# columns, change derive_genotype_label and GENOTYPE_LEVELS_IN_ORDER together.
# ------------------------------------------------------------------------------
NORMALIZATION_VALUE_COLUMN <- "value_corrected"

# ------------------------------------------------------------------------------
# Centralized plot configuration. fill_colors is the canonical palette copied
# verbatim from the loading reference scripts (orc4r-screen_loading-*.R) so a
# genotype keeps one colour across every figure. Label form is "+Nsofa" (the
# reference convention), NOT "4R +orcNsofa". "No ORC" has NO colour key because it
# is dropped from the plot entirely (used only for background subtraction upstream,
# never a bar); all seven canonical labels are listed so this one script serves
# both the orc1/3/4 and orc5/6 loading manifests without edits.
# ------------------------------------------------------------------------------
PLOT_CONFIG <- list(
    fill_colors = c(
        "WT" = "#E41A1C",
        "ORC4R" = "#377EB8",
        "+1sofa" = "#FF7F00",
        "+3sofa" = "#984EA3",
        "+4sofa" = "#4DAF4A",
        "+5sofa" = "#FFFF33",
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

# Genotype axis order, fixed so bars read left-to-right in experiment order. WT is
# the leftmost bar (the normalization reference, reads 1.0), then the ORC4R mutant,
# then the suppressors by orcN. The list is the union across BOTH loading datasets
# (orc1/3/4 and orc5/6) so one script plots either manifest; an absent genotype
# yields no bar. "No ORC" is deliberately NOT a level: it is dropped before plotting
# (background-subtraction control only), and including it would create an empty
# phantom bar. A derived label not in this list is a hard failure (see the NA
# assertion after derive_genotype_label), never a silent drop -- that silent drop is
# what collapsed the orc5/6 bars in the gel-shift script.
GENOTYPE_LEVELS_IN_ORDER <- c(
    "WT", "ORC4R", "+1sofa", "+3sofa", "+4sofa", "+5sofa", "+6sofa"
)

# ------------------------------------------------------------------------------
# Arguments.
# ------------------------------------------------------------------------------
command_arguments <- commandArgs(trailingOnly = TRUE)
if (length(command_arguments) != 2) {
    stop(
        "usage: Rscript plot_loading.R <manifest.csv> <output_directory>",
        call. = FALSE
    )
}
manifest_filepath <- command_arguments[[1]]
output_directory <- command_arguments[[2]]
if (!file.exists(manifest_filepath)) {
    stop("manifest not found: ", manifest_filepath, call. = FALSE)
}
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Read the manifest and load every referenced loading CSV, tagging each row with
# its replicate so within-gel pairing survives into the plots.
# ------------------------------------------------------------------------------
manifest_data <- read_csv(manifest_filepath, show_col_types = FALSE)
missing_manifest_columns <- setdiff(c("replicate", "loading_csv"), names(manifest_data))
if (length(missing_manifest_columns) > 0) {
    stop(
        "manifest is missing column(s): ",
        paste(missing_manifest_columns, collapse = ", "),
        call. = FALSE
    )
}

# DATA-MODEL SHIM (temporary; scoped to the orc1/3/4 and orc5/6 loading datasets).
# This pipeline spells the factors differently from the canonical loading scripts
# (orc4r-screen_loading-*.R): here orc4-R267 is {WT, 4R, none} and suppressor is
# {none, orcN-sofa}; canon uses orc4 {WT, RA} and sofa {none, orcN}. The spellings
# are mapped to canonical labels HERE rather than by rewriting every ratio CSV. The
# long-term fix is a canonical data model; this shim exists so the current datasets
# plot now. Loading has NO ATP factor, so there is no atp label here at all.
#
# Branch ORDER matters: No ORC (orc4-R267 == "none") and ORC4R (orc4-R267 == "4R",
# suppressor == "none") both have suppressor == "none" and are separable ONLY by
# orc4-R267, so orc4-R267 is tested first. A suppressor-first rule merges the two.
# No ORC is labelled here so the normalization step can find and exclude it, then it
# is dropped before plotting.
derive_genotype_label <- function(orc4_r267_value, suppressor_value) {
    if (orc4_r267_value == "none") {
        return("No ORC")
    }
    if (orc4_r267_value == "WT") {
        return("WT")
    }
    # orc4-R267 == "4R" (canon ORC4R). No suppressor -> ORC4R; else the suppressor.
    if (suppressor_value == "none") {
        return("ORC4R")
    }
    # Only KNOWN suppressors are accepted; a well-formed but unknown value (e.g.
    # "orc9-sofa") must NOT become a plausible "+9sofa" that then vanishes at
    # factor(). Unknown -> NA, caught by the assertion below. Add new suppressors
    # here, to GENOTYPE_LEVELS_IN_ORDER, and to fill_colors together.
    KNOWN_SUPPRESSOR_LABELS <- c(
        "orc1-sofa" = "+1sofa",
        "orc3-sofa" = "+3sofa",
        "orc4-sofa" = "+4sofa",
        "orc5-sofa" = "+5sofa",
        "orc6-sofa" = "+6sofa"
    )
    if (!(suppressor_value %in% names(KNOWN_SUPPRESSOR_LABELS))) {
        return(NA_character_)
    }
    return(unname(KNOWN_SUPPRESSOR_LABELS[[suppressor_value]]))
}

# ------------------------------------------------------------------------------
# Per-gel normalization to the WT lane. Each gel is normalized to its OWN WT lane
# BEFORE gels are combined, because loading exposure differs per gel; pooling raw
# values across gels then dividing would mix exposures. The reference is WT (it
# reads 1.0 = 100%), NOT the analysis_role == "reference" marker: "reference is WT"
# is expressed literally here so the result does not depend on where the reference
# role happens to be marked in the sheet. Genotype must be derived first (to find
# the WT lane), so derive_genotype_label is applied inside this loop. Exactly one WT
# lane with a positive value is required per gel, or the gel is a hard failure --
# normalizing to zero or an ambiguous WT is meaningless.
#
# The No ORC lane is NOT the reference here: it is the low background control
# (subtracted upstream by the calc), reads low, and is dropped before plotting.
# ------------------------------------------------------------------------------
per_gel_frames <- vector(mode = "list", length = nrow(manifest_data))
for (manifest_row_index in seq_len(nrow(manifest_data))) {
    replicate_value <- manifest_data$replicate[[manifest_row_index]]
    loading_csv_path <- manifest_data$loading_csv[[manifest_row_index]]
    if (!file.exists(loading_csv_path)) {
        stop(
            "loading CSV listed in manifest not found: ", loading_csv_path,
            call. = FALSE
        )
    }
    one_gel <- read_csv(loading_csv_path, show_col_types = FALSE)

    required_columns <- c("orc4-R267", "suppressor", NORMALIZATION_VALUE_COLUMN)
    missing_columns <- setdiff(required_columns, names(one_gel))
    if (length(missing_columns) > 0) {
        stop(
            "loading CSV ", loading_csv_path, " is missing column(s): ",
            paste(missing_columns, collapse = ", "),
            " (produced by calculate_loading_value.py? sample sheet must carry ",
            "orc4-R267 and suppressor)",
            call. = FALSE
        )
    }

    # Derive genotype for this gel so the WT lane can be found by label.
    one_gel$genotype <- map2_chr(
        one_gel$`orc4-R267`, one_gel$suppressor, derive_genotype_label
    )
    wt_rows <- one_gel %>% filter(.data$genotype == "WT")
    if (nrow(wt_rows) != 1) {
        stop(
            "gel ", loading_csv_path, " has ", nrow(wt_rows),
            " WT lane(s) (orc4-R267 == 'WT', suppressor == 'none'); exactly one is ",
            "required to normalize against (reference is WT).",
            call. = FALSE
        )
    }
    wt_value <- as.numeric(wt_rows[[NORMALIZATION_VALUE_COLUMN]][[1]])
    if (is.na(wt_value) || wt_value <= 0) {
        stop(
            "gel ", loading_csv_path, " WT lane has a non-positive ",
            NORMALIZATION_VALUE_COLUMN, " (", wt_value,
            "); cannot normalize against it.",
            call. = FALSE
        )
    }

    one_gel$replicate <- as.character(replicate_value)
    one_gel$loading_ratio <-
        as.numeric(one_gel[[NORMALIZATION_VALUE_COLUMN]]) / wt_value
    per_gel_frames[[manifest_row_index]] <- one_gel
}
combined_data <- bind_rows(per_gel_frames)

# ------------------------------------------------------------------------------
# Derive the genotype factor. Genotype is the combined (orc4-R267, suppressor) key
# the experiment varies -- orc4-R267 alone would collapse the three suppressors into
# one label. This is the same derivation gel shift used; loading differs only in
# having no ATP axis. A lane with an undefined loading_ratio (no signal, so the calc
# clipped it to zero and the ratio is 0, or a non-analyte lane) is handled below.
# Non-analyte lanes (ladder/empty, analysis_role == "excluded") are dropped: they
# carry no genotype and were emitted by the calc only so the reference lane and the
# full gel were available for normalization.
# ------------------------------------------------------------------------------
# Drop non-analyte lanes FIRST (ladder/empty/control marked analysis_role ==
# "excluded"): they carry non-genotype factor values (e.g. not_applicable) that do
# not map to a genotype, and were emitted by the calc only so the full gel was
# available for normalization. Filtering them before the NA assertion below means
# the assertion checks only true analyte lanes, so a ladder lane cannot trip it.
analytes_only <- combined_data %>%
    filter(.data$analysis_role != "excluded")

# Genotype was already derived per gel (needed to find the WT lane); coerce the
# ratio to numeric here.
analysis_data <- analytes_only %>%
    mutate(loading_ratio = as.numeric(.data$loading_ratio))

# NA-label assertion, ported from the canonical loading scripts and the gel-shift
# fix. This is the guard whose ABSENCE let the orc5/6 bars vanish: an unmapped
# genotype label silently became NA at factor() and pooled. Catch it loudly here,
# naming the exact (orc4-R267, suppressor) combinations that failed to map, so an
# unrecognized factor value stops the run instead of corrupting the figure.
unmapped_genotype_rows <- analytes_only[is.na(analysis_data$genotype), ]
if (nrow(unmapped_genotype_rows) > 0) {
    offending_combinations <- unique(
        paste0(
            "orc4-R267=", unmapped_genotype_rows$`orc4-R267`,
            ", suppressor=", unmapped_genotype_rows$suppressor
        )
    )
    stop(
        "derive_genotype_label produced NA for ", nrow(unmapped_genotype_rows),
        " analyte lane-row(s); unmapped factor combination(s):\n  ",
        paste(offending_combinations, collapse = "\n  "),
        "\nAdd the mapping to derive_genotype_label (and the label to ",
        "GENOTYPE_LEVELS_IN_ORDER and fill_colors) before plotting.",
        call. = FALSE
    )
}

undefined_row_count <- sum(is.na(analysis_data$loading_ratio))
if (undefined_row_count > 0) {
    message(
        "dropping ", undefined_row_count,
        " lane-row(s) with an undefined loading_ratio (no value to normalize)."
    )
    analysis_data <- analysis_data %>% filter(!is.na(.data$loading_ratio))
}

# Drop No ORC from the plot. It is the low background control (subtracted upstream),
# not an analyte bar; you asked it never appear. It is dropped AFTER normalization
# (where it played no part -- WT was the divisor) and AFTER the NA guard (it maps to
# a real label, "No ORC", so it is not an error). It is also absent from
# GENOTYPE_LEVELS_IN_ORDER, so it would drop to NA at factor() anyway; removing it
# explicitly here keeps the intent visible and avoids a spurious "dropped to NA"
# from the post-factor guard below.
no_orc_row_count <- sum(analysis_data$genotype == "No ORC", na.rm = TRUE)
if (no_orc_row_count > 0) {
    message(
        "dropping ", no_orc_row_count,
        " No ORC lane-row(s) from the plot (background-subtraction control only)."
    )
    analysis_data <- analysis_data %>% filter(.data$genotype != "No ORC")
}

# Capture labels before factoring so a valid label absent from the levels list (a
# mapping/levels drift) is caught as a loud stop rather than a silent NA drop -- the
# original bug, one layer down. No ORC has already been removed above, so it will
# not appear here.
genotype_labels_before_factoring <- analysis_data$genotype
labels_dropped_by_factoring <- setdiff(
    unique(genotype_labels_before_factoring[!is.na(genotype_labels_before_factoring)]),
    GENOTYPE_LEVELS_IN_ORDER
)
if (length(labels_dropped_by_factoring) > 0) {
    stop(
        "genotype label(s) not in GENOTYPE_LEVELS_IN_ORDER, would drop to NA at ",
        "factor(): ", paste(labels_dropped_by_factoring, collapse = ", "),
        "\nAdd them to GENOTYPE_LEVELS_IN_ORDER and fill_colors.",
        call. = FALSE
    )
}

analysis_data <- analysis_data %>%
    mutate(genotype = factor(.data$genotype, levels = GENOTYPE_LEVELS_IN_ORDER))

# A single condition key (genotype) used as the grouping unit for the summary and
# the all-against-all matrices. Loading has no second factor, so condition ==
# genotype; the column is kept named "condition" for symmetry with the gel-shift
# outputs the paper cross-references.
analysis_data <- analysis_data %>%
    mutate(condition = as.character(.data$genotype))

# ------------------------------------------------------------------------------
# Per-condition summary: mean, sd, cv, and standard error. n is the number of
# replicate-gel lanes contributing; se is reported but with n this small it is a
# rough quantity, so it is labelled as such and the plot uses sd, not se.
# ------------------------------------------------------------------------------
summary_data <- analysis_data %>%
    group_by(.data$genotype, .data$condition) %>%
    summarise(
        replicate_count = n(),
        mean_loading_ratio = mean(.data$loading_ratio),
        sd_loading_ratio = sd(.data$loading_ratio),
        cv_loading_ratio = if_else(
            mean(.data$loading_ratio) != 0,
            sd(.data$loading_ratio) / abs(mean(.data$loading_ratio)),
            NA_real_
        ),
        # Standard error of the mean; rough at n=3, provided for completeness.
        se_loading_ratio = sd(.data$loading_ratio) / sqrt(n()),
        .groups = "drop"
    ) %>%
    arrange(.data$genotype)

summary_csv_path <- file.path(output_directory, "loading_summary.csv")
write.csv(summary_data, summary_csv_path, row.names = FALSE)
message("wrote ", summary_csv_path)

# ------------------------------------------------------------------------------
# All-against-all matrices on the per-condition MEANS (one value per condition,
# pooled over replicate gels). Fold change is mean_row / mean_col; percent
# difference is 100 * (mean_row - mean_col) / mean_col. Both are directional and
# read "row relative to column". Referenced in the paper, not plotted. This is the
# orphaned capability the R-script triage flagged (all-against-all percent/fold),
# reused here.
# ------------------------------------------------------------------------------
condition_order <- summary_data$condition
mean_by_condition <- summary_data$mean_loading_ratio
names(mean_by_condition) <- condition_order

fold_change_matrix <- outer(
    mean_by_condition, mean_by_condition,
    FUN = function(row_mean, col_mean) row_mean / col_mean
)
percent_difference_matrix <- outer(
    mean_by_condition, mean_by_condition,
    FUN = function(row_mean, col_mean) 100 * (row_mean - col_mean) / col_mean
)
# Write with an explicit leading "condition" column so the row identity survives in
# the CSV (row.names alone are easy to lose on re-import).
fold_change_frame <- data.frame(
    condition = condition_order, fold_change_matrix, check.names = FALSE
)
percent_difference_frame <- data.frame(
    condition = condition_order, percent_difference_matrix, check.names = FALSE
)
fold_change_csv_path <- file.path(
    output_directory, "loading_fold_change_all_vs_all.csv"
)
percent_difference_csv_path <- file.path(
    output_directory, "loading_percent_difference_all_vs_all.csv"
)
write.csv(fold_change_frame, fold_change_csv_path, row.names = FALSE)
write.csv(percent_difference_frame, percent_difference_csv_path, row.names = FALSE)
message("wrote ", fold_change_csv_path)
message("wrote ", percent_difference_csv_path)

# ------------------------------------------------------------------------------
# A reusable theme matching the reference scripts.
# ------------------------------------------------------------------------------
house_theme <- theme_classic(
    base_size = PLOT_CONFIG$theme$base_size,
    base_family = PLOT_CONFIG$theme$base_family
) +
    theme(
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = PLOT_CONFIG$theme$legend_position,
        panel.spacing = unit(1, "lines")
    )

# ==============================================================================
# The loading figure: x = genotype, one bar per genotype, fill = genotype, SD error
# bars, per-replicate dots. No ATP dodge (loading has no ATP contrast), so this is
# the gel-shift P1 stripped to a single bar per genotype.
# ==============================================================================
loading_plot <- ggplot(
    summary_data,
    aes(
        x = .data$genotype, y = .data$mean_loading_ratio,
        fill = .data$genotype
    )
) +
    geom_col(
        width = PLOT_CONFIG$bar$width,
        color = PLOT_CONFIG$bar$color,
        linewidth = PLOT_CONFIG$bar$linewidth
    ) +
    geom_errorbar(
        aes(
            ymin = pmax(0, .data$mean_loading_ratio - .data$sd_loading_ratio),
            ymax = .data$mean_loading_ratio + .data$sd_loading_ratio
        ),
        width = PLOT_CONFIG$errorbar$width,
        linewidth = PLOT_CONFIG$errorbar$linewidth
    ) +
    geom_point(
        data = analysis_data,
        aes(
            x = .data$genotype, y = .data$loading_ratio,
            shape = factor(.data$replicate)
        ),
        position = position_jitter(
            width = PLOT_CONFIG$point$jitter_width,
            height = 0,
            seed = PLOT_CONFIG$point$jitter_seed
        ),
        size = PLOT_CONFIG$point$size,
        fill = PLOT_CONFIG$point$fill,
        color = PLOT_CONFIG$point$color,
        stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE
    ) +
    scale_fill_manual(
        values = PLOT_CONFIG$fill_colors,
        name = "Genotype",
        breaks = GENOTYPE_LEVELS_IN_ORDER
    ) +
    scale_shape_manual(
        values = PLOT_CONFIG$replicate_shapes,
        labels = c("1" = "Replicate 1", "2" = "Replicate 2", "3" = "Replicate 3"),
        name = "Replicate"
    ) +
    # Y-axis is intentionally AUTO-SCALED (limits floated, 18% headroom above the
    # tallest bar). The loading ratio is lane / WT, which is UNBOUNDED above -- a
    # genotype that loads more than WT reads > 1.0 -- so a fixed upper limit like the
    # gel-shift plot's 1.2 would silently clip real data to NA. Do not add fixed
    # limits here without confirming no genotype exceeds them.
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        x = "Genotype", y = "Loading  [ lane / WT ]",
        title = "Loading by genotype (normalized to WT)"
    ) +
    guides(
        fill = guide_legend(nrow = 1, order = 1),
        shape = guide_legend(nrow = 1, order = 2)
    ) +
    house_theme
loading_plot_path <- file.path(output_directory, "loading_by_genotype.pdf")
ggsave(
    loading_plot_path, plot = loading_plot,
    device = PLOT_CONFIG$output$device,
    width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height
)
message("wrote ", loading_plot_path)

message("done.")
