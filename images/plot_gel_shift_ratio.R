# ==============================================================================
# plot_gel_shift_ratio.R
#
# Aggregate gel shift bound-fraction CSVs across replicate gels, write the
# summary and all-against-all matrices the paper references, and render three
# candidate figures (one is chosen and finished in Illustrator).
#
# Input is a MANIFEST CSV (not a directory, not a TIFF): one row per replicate
# gel, columns
#     replicate,ratio_csv
# where ratio_csv is a path to a gel_shift_ratio_*.csv produced per gel by
# gel_shift_ratio.py. Listing the CSVs directly avoids re-deriving the per-gel
# filename, whose millimetre window differs between gels.
#
# Usage (Rscript, headless):
#     Rscript plot_gel_shift_ratio.R <manifest.csv> <output_directory>
#
# The reported quantity is bound_fraction = bound / (bound + free), already a
# self-normalized proportion in [0, 1]; there is deliberately no normalization to
# WT (unlike the loading analysis), because a proportion needs no reference.
#
# Dependencies: tidyverse only. The loading reference scripts also load readxl
# (Excel input; not needed here, we read CSVs) and nlme (a mixed-model diagnostic
# not requested here); both are omitted so this stays a strict subset. extrafont
# is intentionally not loaded: base_family = "Arial" is requested in the theme and
# R falls back cleanly if the font is not registered, so the script runs headless.
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# Centralized plot configuration, carried over verbatim from the loading
# reference scripts so the gel shift figure matches the companion figures. Any
# parameter set to NULL means "use ggplot default." fill_colors is keyed by the
# gel shift genotype display names (below), reusing the reference Set1 hues so a
# genotype keeps one colour across figures.
# ------------------------------------------------------------------------------
PLOT_CONFIG <- list(
    fill_colors = c(
        "WT" = "#E41A1C",
        "4R +orc1sofa" = "#FF7F00",
        "4R +orc3sofa" = "#984EA3",
        "4R +orc4sofa" = "#4DAF4A",
        "No ORC" = "#999999"
    ),
    bar = list(width = 0.7, color = "black", linewidth = 0.4),
    errorbar = list(width = 0.25, linewidth = 0.6),
    # noATP bars are the genotype colour lightened by this fraction toward white,
    # so genotype reads by hue and ATP state reads by shade (plusATP = full colour,
    # the stronger-binding condition, so the fuller colour). No extra dependency:
    # the lighten helper below is base grDevices.
    noatp_lighten_amount = 0.55,
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

# The ATP axis order and the genotype axis order are fixed here so bars and facets
# read left-to-right in experiment order rather than alphabetically.
ATP_LEVELS_IN_ORDER <- c("noATP", "plusATP")
GENOTYPE_LEVELS_IN_ORDER <- c(
    "No ORC", "WT", "4R +orc1sofa", "4R +orc3sofa", "4R +orc4sofa"
)

# ------------------------------------------------------------------------------
# Arguments.
# ------------------------------------------------------------------------------
command_arguments <- commandArgs(trailingOnly = TRUE)
if (length(command_arguments) != 2) {
    stop(
        "usage: Rscript plot_gel_shift_ratio.R <manifest.csv> <output_directory>",
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
# Read the manifest and load every referenced ratio CSV, tagging each row with
# its replicate so within-gel pairing survives into the plots.
# ------------------------------------------------------------------------------
manifest_data <- read_csv(manifest_filepath, show_col_types = FALSE)
missing_manifest_columns <- setdiff(c("replicate", "ratio_csv"), names(manifest_data))
if (length(missing_manifest_columns) > 0) {
    stop(
        "manifest is missing column(s): ",
        paste(missing_manifest_columns, collapse = ", "),
        call. = FALSE
    )
}

per_gel_frames <- vector(mode = "list", length = nrow(manifest_data))
for (manifest_row_index in seq_len(nrow(manifest_data))) {
    replicate_value <- manifest_data$replicate[[manifest_row_index]]
    ratio_csv_path <- manifest_data$ratio_csv[[manifest_row_index]]
    if (!file.exists(ratio_csv_path)) {
        stop(
            "ratio CSV listed in manifest not found: ", ratio_csv_path,
            call. = FALSE
        )
    }
    one_gel <- read_csv(ratio_csv_path, show_col_types = FALSE)
    one_gel$replicate <- as.character(replicate_value)
    per_gel_frames[[manifest_row_index]] <- one_gel
}
combined_data <- bind_rows(per_gel_frames)

# ------------------------------------------------------------------------------
# Derive the two design factors. Genotype is the combined (orc4-R267, suppressor)
# key the experiment actually varies: orc4-R267 alone only separates WT from 4R
# and would collapse the three suppressors into one colour. ATP comes from the
# atp_presence column (yes/no), which the sample sheet is the source of truth for.
# A lane with an undefined bound_fraction (both bands absent) is dropped from the
# analysis with a message; it carries no proportion to average.
# ------------------------------------------------------------------------------
derive_genotype_label <- function(orc4_r267_value, suppressor_value) {
    # not_applicable/not_applicable is the No ORC control lane.
    if (orc4_r267_value == "not_applicable") {
        return("No ORC")
    }
    if (orc4_r267_value == "WT") {
        return("WT")
    }
    # 4R carries a suppressor; name it by the suppressor's orcN identity.
    if (suppressor_value == "orc1-sofa") return("4R +orc1sofa")
    if (suppressor_value == "orc3-sofa") return("4R +orc3sofa")
    if (suppressor_value == "orc4-sofa") return("4R +orc4sofa")
    return(paste0("4R +", suppressor_value))
}
derive_atp_label <- function(atp_presence_value) {
    if (atp_presence_value == "yes") return("plusATP")
    if (atp_presence_value == "no") return("noATP")
    return(NA_character_)  # not_applicable (the No ORC control)
}

analysis_data <- combined_data %>%
    mutate(
        genotype = map2_chr(.data$`orc4-R267`, .data$suppressor, derive_genotype_label),
        atp = map_chr(.data$atp_presence, derive_atp_label),
        bound_fraction = as.numeric(.data$bound_fraction)
    )

undefined_row_count <- sum(is.na(analysis_data$bound_fraction))
if (undefined_row_count > 0) {
    message(
        "dropping ", undefined_row_count,
        " lane-row(s) with an undefined bound_fraction (no signal in either band)."
    )
    analysis_data <- analysis_data %>% filter(!is.na(.data$bound_fraction))
}

# The No ORC control has no ATP condition; keep it as its own single-bar genotype
# for the baseline, and give it an atp level so it can sit on the axis. It is
# assigned to noATP purely for placement and is excluded from the ATP contrast.
analysis_data <- analysis_data %>%
    mutate(
        atp = if_else(.data$genotype == "No ORC" & is.na(.data$atp), "noATP", .data$atp),
        genotype = factor(.data$genotype, levels = GENOTYPE_LEVELS_IN_ORDER),
        atp = factor(.data$atp, levels = ATP_LEVELS_IN_ORDER)
    )

# A single condition key (genotype x ATP) used as the grouping unit for the
# summary and for the all-against-all matrices.
analysis_data <- analysis_data %>%
    mutate(
        condition = if_else(
            .data$genotype == "No ORC",
            "No ORC",
            paste0(as.character(.data$genotype), " ", as.character(.data$atp))
        )
    )

# ------------------------------------------------------------------------------
# Per-condition summary: mean, sd, cv, and standard error. n is the number of
# replicate gels contributing; se is reported but with n this small it is a rough
# quantity, so it is labelled as such and the plots use sd, not se.
# ------------------------------------------------------------------------------
summary_data <- analysis_data %>%
    group_by(.data$genotype, .data$atp, .data$condition) %>%
    summarise(
        replicate_count = n(),
        mean_bound_fraction = mean(.data$bound_fraction),
        sd_bound_fraction = sd(.data$bound_fraction),
        cv_bound_fraction = if_else(
            mean(.data$bound_fraction) != 0,
            sd(.data$bound_fraction) / abs(mean(.data$bound_fraction)),
            NA_real_
        ),
        # Standard error of the mean; rough at n=3, provided for completeness.
        se_bound_fraction = sd(.data$bound_fraction) / sqrt(n()),
        .groups = "drop"
    ) %>%
    arrange(.data$genotype, .data$atp)

summary_csv_path <- file.path(output_directory, "gel_shift_summary.csv")
write.csv(summary_data, summary_csv_path, row.names = FALSE)
message("wrote ", summary_csv_path)

# ------------------------------------------------------------------------------
# All-against-all matrices on the per-condition MEANS (one value per condition,
# pooled over replicate gels). Fold change is mean_row / mean_col; percent
# difference is 100 * (mean_row - mean_col) / mean_col. Both are directional and
# read "row relative to column". Referenced in the paper, not plotted.
# ------------------------------------------------------------------------------
condition_order <- summary_data$condition
mean_by_condition <- summary_data$mean_bound_fraction
names(mean_by_condition) <- condition_order

fold_change_matrix <- outer(
    mean_by_condition, mean_by_condition,
    FUN = function(row_mean, col_mean) row_mean / col_mean
)
percent_difference_matrix <- outer(
    mean_by_condition, mean_by_condition,
    FUN = function(row_mean, col_mean) 100 * (row_mean - col_mean) / col_mean
)
# Write with an explicit leading "condition" column so the row identity survives
# in the CSV (row.names alone are easy to lose on re-import).
fold_change_frame <- data.frame(
    condition = condition_order, fold_change_matrix, check.names = FALSE
)
percent_difference_frame <- data.frame(
    condition = condition_order, percent_difference_matrix, check.names = FALSE
)
fold_change_csv_path <- file.path(
    output_directory, "gel_shift_fold_change_all_vs_all.csv"
)
percent_difference_csv_path <- file.path(
    output_directory, "gel_shift_percent_difference_all_vs_all.csv"
)
write.csv(fold_change_frame, fold_change_csv_path, row.names = FALSE)
write.csv(percent_difference_frame, percent_difference_csv_path, row.names = FALSE)
message("wrote ", fold_change_csv_path)
message("wrote ", percent_difference_csv_path)

# ------------------------------------------------------------------------------
# Data subsets for plotting. The ATP contrast plots exclude the No ORC control
# (it has no ATP pair); it appears only where a baseline is drawn.
# ------------------------------------------------------------------------------
summary_with_atp <- summary_data %>%
    filter(.data$genotype != "No ORC") %>%
    mutate(genotype = droplevels(.data$genotype))
points_with_atp <- analysis_data %>%
    filter(.data$genotype != "No ORC") %>%
    mutate(genotype = droplevels(.data$genotype))
no_orc_mean <- summary_data %>%
    filter(.data$genotype == "No ORC") %>%
    pull(.data$mean_bound_fraction)
no_orc_baseline <- if (length(no_orc_mean) == 1) no_orc_mean else NA_real_

# A reusable theme matching the reference scripts.
house_theme <- theme_classic(
    base_size = PLOT_CONFIG$theme$base_size,
    base_family = PLOT_CONFIG$theme$base_family
) +
    theme(
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = PLOT_CONFIG$theme$legend_position,
        # Stack the fill and shape legends as separate boxes. With the eight-entry
        # genotype/ATP fill legend laid out in two rows, a side-by-side shape legend
        # gets consumed by ggplot's guide layout and collapses to one key; stacking
        # them vertically keeps the three replicate keys intact.
        legend.box = "vertical",
        panel.spacing = unit(1, "lines")
    )

# Fill is genotype x ATP: each genotype keeps its palette hue, and noATP is that
# hue lightened toward white so the ATP state reads by shade. lighten_colour mixes
# the colour toward white by a fraction; base grDevices only, no new dependency.
lighten_colour <- function(hex_colour, amount) {
    rgb_fraction <- col2rgb(hex_colour) / 255
    lightened_fraction <- rgb_fraction + (1 - rgb_fraction) * amount
    rgb(lightened_fraction[1], lightened_fraction[2], lightened_fraction[3])
}
genotype_atp_fill_values <- c()
genotype_atp_fill_labels <- c()
for (genotype_name in GENOTYPE_LEVELS_IN_ORDER) {
    if (genotype_name == "No ORC") {
        next  # No ORC has no ATP pair; it is not a filled bar in the ATP plots
    }
    full_colour <- PLOT_CONFIG$fill_colors[[genotype_name]]
    plus_key <- paste0(genotype_name, ".plusATP")
    no_key <- paste0(genotype_name, ".noATP")
    genotype_atp_fill_values[[plus_key]] <- full_colour
    genotype_atp_fill_values[[no_key]] <- lighten_colour(
        full_colour, PLOT_CONFIG$noatp_lighten_amount
    )
    genotype_atp_fill_labels[[plus_key]] <- paste0(genotype_name, " +ATP")
    genotype_atp_fill_labels[[no_key]] <- paste0(genotype_name, " -ATP")
}
# The interaction key ties each bar to its genotype-and-ATP fill. The factor level
# order (genotype outer, ATP inner with noATP first) matches the dodge order so the
# lighter noATP bar sits left of the fuller plusATP bar within each genotype.
genotype_atp_level_order <- unlist(lapply(
    setdiff(GENOTYPE_LEVELS_IN_ORDER, "No ORC"),
    function(genotype_name) c(
        paste0(genotype_name, ".noATP"), paste0(genotype_name, ".plusATP")
    )
))
summary_with_atp <- summary_with_atp %>%
    mutate(
        genotype_atp = factor(
            paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
            levels = genotype_atp_level_order
        )
    )
points_with_atp <- points_with_atp %>%
    mutate(
        genotype_atp = factor(
            paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
            levels = genotype_atp_level_order
        )
    )

# ==============================================================================
# P1: x = genotype, dodged noATP/plusATP bars, SD error bars, replicate dots.
# One panel; the ATP pair sits adjacent within each genotype for direct contrast.
# Fill encodes genotype (hue) and ATP (shade): noATP lighter, plusATP full.
# ==============================================================================
dodge_width <- 0.8
plot_p1 <- ggplot(
    summary_with_atp,
    aes(
        x = .data$genotype, y = .data$mean_bound_fraction,
        fill = .data$genotype_atp, group = .data$atp
    )
) +
    geom_col(
        position = position_dodge(width = dodge_width),
        width = PLOT_CONFIG$bar$width,
        color = PLOT_CONFIG$bar$color,
        linewidth = PLOT_CONFIG$bar$linewidth
    ) +
    geom_errorbar(
        aes(
            ymin = pmax(0, .data$mean_bound_fraction - .data$sd_bound_fraction),
            ymax = .data$mean_bound_fraction + .data$sd_bound_fraction
        ),
        position = position_dodge(width = dodge_width),
        width = PLOT_CONFIG$errorbar$width,
        linewidth = PLOT_CONFIG$errorbar$linewidth
    ) +
    geom_point(
        data = points_with_atp,
        aes(
            x = .data$genotype, y = .data$bound_fraction,
            group = .data$atp, shape = factor(.data$replicate)
        ),
        position = position_jitterdodge(
            jitter.width = PLOT_CONFIG$point$jitter_width,
            dodge.width = dodge_width,
            seed = PLOT_CONFIG$point$jitter_seed
        ),
        size = PLOT_CONFIG$point$size,
        fill = PLOT_CONFIG$point$fill,
        color = PLOT_CONFIG$point$color,
        stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE
    ) +
    scale_fill_manual(
        values = genotype_atp_fill_values,
        labels = genotype_atp_fill_labels,
        name = "Genotype / ATP",
        breaks = genotype_atp_level_order
    ) +
    scale_shape_manual(
        values = PLOT_CONFIG$replicate_shapes,
        labels = c("1" = "Replicate 1", "2" = "Replicate 2", "3" = "Replicate 3"),
        name = "Replicate"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        x = "Genotype", y = "Bound fraction  [ bound / (bound + free) ]",
        title = "Gel shift bound fraction by genotype and ATP"
    ) +
    guides(
        # The fill scale has eight entries (four genotypes x two ATP states); lay it
        # out in two rows so the legend stays compact. The shape guide lists the
        # three replicates; keep it a single row beside the fill legend.
        fill = guide_legend(nrow = 2, order = 1),
        shape = guide_legend(nrow = 1, order = 2)
    ) +
    house_theme
p1_path <- file.path(output_directory, "gel_shift_P1_genotype_dodged_ATP.pdf")
ggsave(
    p1_path, plot = plot_p1,
    device = PLOT_CONFIG$output$device,
    width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height
)
message("wrote ", p1_path)

# ==============================================================================
# P6: P1 plus the No ORC control drawn as a horizontal baseline, so "how far
# above background" is legible for every bar.
# ==============================================================================
plot_p6 <- plot_p1
if (!is.na(no_orc_baseline)) {
    plot_p6 <- plot_p6 +
        geom_hline(
            yintercept = no_orc_baseline,
            linetype = "dashed", color = "grey40", linewidth = 0.5
        ) +
        annotate(
            "text",
            x = 0.6, y = no_orc_baseline, vjust = -0.4, hjust = 0,
            label = "No ORC", size = 3, color = "grey40",
            family = PLOT_CONFIG$theme$base_family
        ) +
        labs(title = "Gel shift bound fraction by genotype and ATP (No-ORC baseline)")
}
p6_path <- file.path(output_directory, "gel_shift_P6_genotype_dodged_ATP_baseline.pdf")
ggsave(
    p6_path, plot = plot_p6,
    device = PLOT_CONFIG$output$device,
    width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height
)
message("wrote ", p6_path)

# ==============================================================================
# P5: per-gel slopegraph. Each replicate gel is a line from noATP to plusATP,
# faceted by genotype, so within-gel reproducibility of the ATP effect is visible
# rather than averaged away. Uses the raw per-lane points, not the summary.
# ==============================================================================
plot_p5 <- ggplot(
    points_with_atp,
    aes(
        x = .data$atp, y = .data$bound_fraction,
        group = .data$replicate, shape = factor(.data$replicate)
    )
) +
    geom_line(color = "grey50", linewidth = 0.5) +
    geom_point(
        size = PLOT_CONFIG$point$size,
        fill = PLOT_CONFIG$point$fill,
        color = PLOT_CONFIG$point$color,
        stroke = PLOT_CONFIG$point$stroke
    ) +
    facet_wrap(~ .data$genotype, nrow = 1) +
    scale_shape_manual(
        values = PLOT_CONFIG$replicate_shapes,
        labels = c("1" = "Replicate 1", "2" = "Replicate 2", "3" = "Replicate 3"),
        name = "Replicate"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        x = "ATP", y = "Bound fraction  [ bound / (bound + free) ]",
        title = "Gel shift ATP effect per replicate gel"
    ) +
    house_theme
p5_path <- file.path(output_directory, "gel_shift_P5_slopegraph_per_gel.pdf")
ggsave(
    p5_path, plot = plot_p5,
    device = PLOT_CONFIG$output$device,
    width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height
)
message("wrote ", p5_path)

message("done.")
