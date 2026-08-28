# ==============================================================================
# plot_gel_shift_combined.R  (SUPPLEMENTARY -- does not replace plot_gel_shift_ratio.R)
#
# Combine two gel-shift screens (orc1/3/4 and orc5/6) into one figure set in which
# the shared controls are FOLDED into a single bar each, so a reader does not see
# two different-looking WT bars and wonder which is right. The screens share WT and
# the No ORC control (each run in both screens, so six WT gels and six No ORC gels
# total); ORC4R and each suppressor are present in only one screen. This script does
# not hard-code any of those counts: it pools bound_fraction by (genotype, ATP)
# across whatever screens contain a genotype, so WT folds to n=6, No ORC to n=6, and
# a genotype that appears in one screen stays at its own n. Provenance is preserved
# by colouring each replicate point by its screen, so the folded WT bar visibly
# draws from both screens.
#
# INPUT: two (or more) screen-tagged manifests of RAW per-gel ratio CSVs, NOT the
# summary CSVs. Raw is required because folding a control across screens means
# recomputing the mean and SD from the pooled replicate values; a summary CSV has
# already collapsed to mean/SD/n and cannot be re-pooled to a correct SD. Each
# manifest is the same "replicate,ratio_csv" file the per-screen script consumes.
#
# Usage (Rscript, headless):
#     Rscript plot_gel_shift_combined.R <output_directory> \
#         <screen_label_1>=<manifest_1.csv> <screen_label_2>=<manifest_2.csv> [...]
# e.g.
#     Rscript plot_gel_shift_combined.R out \
#         orc1-3-4=manifest_134.csv orc5-6=manifest_56.csv
#
# OUTPUTS (all in <output_directory>):
#   combined_summary.csv                  pooled mean/sd/cv/se/n per (genotype, ATP),
#                                         with the screen composition per genotype
#   combined_screen_composition.csv       which screens and how many gels feed each
#                                         (genotype, ATP) cell -- the check that WT
#                                         drew from both screens and ORC4R from one
#   combined_noATP.pdf                    one ATP state, WT once, points by screen
#   combined_plusATP.pdf                  the other ATP state, WT once
#   combined_faceted_by_ATP.pdf           both ATP states in one figure, faceted
#
# The per-screen plot_gel_shift_ratio.R is intentionally NOT modified or imported.
# Per house rule the derivation, palette, and guards below are DUPLICATED from it
# verbatim rather than shared, so this script runs standalone and drift fails loudly.
#
# Dependencies: tidyverse only.
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# Palette and level order, DUPLICATED verbatim from plot_gel_shift_ratio.R (house
# rule: no shared module). Canonical hues from orc4r-screen_loading-*.R.
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
    # Screen point colours. Distinct from the genotype fills so a point's screen
    # reads against any bar. Extend when a third screen is added.
    screen_colors = c(
        "orc1-3-4" = "#1B9E77",
        "orc5-6" = "#7570B3"
    ),
    bar = list(width = 0.7, color = "black", linewidth = 0.4),
    errorbar = list(width = 0.25, linewidth = 0.6),
    point = list(size = 2, color = "black", stroke = 0.5,
                 jitter_width = 0.15, jitter_seed = 42),
    # noATP fill is the genotype hue lightened toward white so ATP reads by shade,
    # matching the per-screen P1 (fill = genotype x ATP). Used by the dodged panel.
    noatp_lighten_amount = 0.55,
    # Replicate point shapes, matching the per-screen P1. Replicate numbering is
    # per-screen (each screen has its own 1/2/3), so in the color=screen variant a
    # shape identifies the within-screen replicate and the colour identifies which
    # screen; the two together are unique. Three shapes cover the max replicates per
    # screen; extend if a screen has more than three gels.
    replicate_shapes = c("1" = 21, "2" = 24, "3" = 22),
    # Fixed y-axis carried over from the per-screen P1. bound_fraction is bounded
    # [0, 1] by construction, so 1.2 leaves headroom without ever clipping real data,
    # and the explicit breaks guarantee the 1.0 tick appears (auto-breaks can skip
    # it). Unlike the loading ratio, this quantity is bounded, so a fixed max is safe.
    y_axis_limits = c(0, 1.2),
    y_axis_breaks = seq(0, 1.2, 0.2),
    y_axis_expand = expansion(mult = c(0, 0), add = c(0, 0.02)),
    theme = list(base_family = "Arial", base_size = 12, legend_position = "bottom"),
    output = list(device = cairo_pdf, width = 8.5, height = 4.5)
)

# Genotype axis order, CLUSTERED BY SCREEN/SET so the figure mirrors how the sets
# sit on the gel. WT is placed first as the shared reference control (it is folded
# across BOTH screens, so it belongs to neither set exclusively and is pulled out
# rather than forced into one). Then each set clusters, and within a set the ORC4R
# bar (orc4 = 4R, no suppressor) precedes its suppressors, which sort by number:
#   WT                          shared control (both screens, folded n=6)
#   +1sofa, +3sofa, +4sofa      SET 1  (orc1-3-4 screen)
#   ORC4R, +5sofa, +6sofa       SET 2  (orc5-6 screen; ORC4R is orc5-6-only)
# This is DISPLAY ORDER only -- it does not change the fold, the counts, or any
# value; it only moves where each bar sits. When a third screen is added, insert its
# block here (ORC4R/orc4 bar first if present, then its suppressors by number) and
# add the new suppressor labels to fill_colors and derive_genotype_label together.
# No ORC is deliberately absent (drawn as a baseline line, not a bar). The list must
# contain exactly the canonical analyte labels -- a derived label not here is a hard
# failure at the post-factor guard, and a label here with no data simply draws no bar.
GENOTYPE_LEVELS_IN_ORDER <- c(
    "WT",
    "ORC4R", "+5sofa", "+6sofa",
    "+1sofa", "+3sofa", "+4sofa"
)
ATP_LEVELS_IN_ORDER <- c("noATP", "plusATP")

# ------------------------------------------------------------------------------
# Arguments: first is the output directory, the rest are label=manifest pairs.
# ------------------------------------------------------------------------------
command_arguments <- commandArgs(trailingOnly = TRUE)
if (length(command_arguments) < 2) {
    stop(
        "usage: Rscript plot_gel_shift_combined.R <output_directory> ",
        "<screen_label>=<manifest.csv> [<screen_label>=<manifest.csv> ...]",
        call. = FALSE
    )
}
output_directory <- command_arguments[[1]]
screen_manifest_arguments <- command_arguments[-1]
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# Parse label=manifest pairs. A screen label with an "=" in it would split wrong,
# so split on the FIRST "=" only.
screen_labels <- character(0)
screen_manifest_paths <- character(0)
for (pair_text in screen_manifest_arguments) {
    equals_position <- regexpr("=", pair_text, fixed = TRUE)
    if (equals_position < 1) {
        stop(
            "screen argument is not in <label>=<manifest.csv> form: ", pair_text,
            call. = FALSE
        )
    }
    screen_label <- substr(pair_text, 1, equals_position - 1)
    manifest_path <- substr(pair_text, equals_position + 1, nchar(pair_text))
    if (!file.exists(manifest_path)) {
        stop("manifest not found for screen '", screen_label, "': ", manifest_path,
             call. = FALSE)
    }
    screen_labels <- c(screen_labels, screen_label)
    screen_manifest_paths <- c(screen_manifest_paths, manifest_path)
}
if (any(duplicated(screen_labels))) {
    stop("duplicate screen label(s): ",
         paste(screen_labels[duplicated(screen_labels)], collapse = ", "),
         call. = FALSE)
}
missing_screen_colours <- setdiff(screen_labels, names(PLOT_CONFIG$screen_colors))
if (length(missing_screen_colours) > 0) {
    stop("no point colour configured for screen(s): ",
         paste(missing_screen_colours, collapse = ", "),
         " -- add them to PLOT_CONFIG$screen_colors.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Genotype derivation, DUPLICATED verbatim from plot_gel_shift_ratio.R. Branch on
# orc4-R267 first so No ORC (none) and ORC4R (4R, none) do not collide; only known
# suppressors map, an unknown one returns NA for the loud assertion below.
# ------------------------------------------------------------------------------
derive_genotype_label <- function(orc4_r267_value, suppressor_value) {
    if (orc4_r267_value == "none") {
        return("No ORC")
    }
    if (orc4_r267_value == "WT") {
        return("WT")
    }
    if (suppressor_value == "none") {
        return("ORC4R")
    }
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

derive_atp_label <- function(atp_presence_value) {
    if (atp_presence_value == "yes") return("plusATP")
    if (atp_presence_value == "no") return("noATP")
    return(NA_character_)
}

# ------------------------------------------------------------------------------
# Read every gel from every screen, tagging each row with its screen label and its
# within-screen replicate. Replicate is kept screen-scoped (screen + replicate) so
# a point's identity survives when two screens both have a "replicate 1".
# ------------------------------------------------------------------------------
per_gel_frames <- list()
for (screen_index in seq_along(screen_labels)) {
    screen_label <- screen_labels[[screen_index]]
    manifest_path <- screen_manifest_paths[[screen_index]]
    manifest_data <- read_csv(manifest_path, show_col_types = FALSE)
    missing_manifest_columns <- setdiff(c("replicate", "ratio_csv"), names(manifest_data))
    if (length(missing_manifest_columns) > 0) {
        stop("manifest ", manifest_path, " missing column(s): ",
             paste(missing_manifest_columns, collapse = ", "), call. = FALSE)
    }
    for (manifest_row_index in seq_len(nrow(manifest_data))) {
        replicate_value <- manifest_data$replicate[[manifest_row_index]]
        ratio_csv_path <- manifest_data$ratio_csv[[manifest_row_index]]
        if (!file.exists(ratio_csv_path)) {
            stop("ratio CSV listed in manifest '", manifest_path, "' not found: ",
                 ratio_csv_path, call. = FALSE)
        }
        one_gel <- read_csv(ratio_csv_path, show_col_types = FALSE,
                            col_types = cols(.default = col_character()))
        # Read all columns as character. Two screens' CSVs can type a shared column
        # differently (e.g. prep_source all-numeric in one screen, alphanumeric in
        # another), and bind_rows refuses to combine a <double> column with a
        # <character> one. Only bound_fraction is used numerically, and it is coerced
        # after binding; every other column is carried as text.
        one_gel$screen <- screen_label
        one_gel$replicate <- as.character(replicate_value)
        one_gel$screen_replicate <- paste0(screen_label, ":", replicate_value)
        per_gel_frames[[length(per_gel_frames) + 1]] <- one_gel
    }
}
combined_data <- bind_rows(per_gel_frames)

# ------------------------------------------------------------------------------
# Derive genotype and ATP, then the loud guards (DUPLICATED from the per-screen
# script). The guard whose absence ate the orc5/6 bars: an unmapped label silently
# becomes NA at factor() and pools; catch it here naming the offending combination.
# ------------------------------------------------------------------------------
combined_data <- combined_data %>%
    mutate(
        genotype = map2_chr(.data$`orc4-R267`, .data$suppressor, derive_genotype_label),
        atp = map_chr(.data$atp_presence, derive_atp_label),
        bound_fraction = as.numeric(.data$bound_fraction)
    )

unmapped_genotype_rows <- combined_data[is.na(combined_data$genotype), ]
if (nrow(unmapped_genotype_rows) > 0) {
    offending_combinations <- unique(paste0(
        "orc4-R267=", unmapped_genotype_rows$`orc4-R267`,
        ", suppressor=", unmapped_genotype_rows$suppressor))
    stop("derive_genotype_label produced NA for ", nrow(unmapped_genotype_rows),
         " lane-row(s); unmapped factor combination(s):\n  ",
         paste(offending_combinations, collapse = "\n  "),
         "\nAdd the mapping to derive_genotype_label, GENOTYPE_LEVELS_IN_ORDER, ",
         "and fill_colors.", call. = FALSE)
}
unmapped_atp_rows <- combined_data[is.na(combined_data$atp), ]
if (nrow(unmapped_atp_rows) > 0) {
    stop("derive_atp_label produced NA for ", nrow(unmapped_atp_rows),
         " lane-row(s); unmapped atp_presence value(s): ",
         paste(unique(unmapped_atp_rows$atp_presence), collapse = ", "), call. = FALSE)
}

# Drop lanes with an undefined bound_fraction (both bands absent): no proportion to
# pool. Message rather than fail; a blank analyte lane is legal.
undefined_row_count <- sum(is.na(combined_data$bound_fraction))
if (undefined_row_count > 0) {
    message("dropping ", undefined_row_count,
            " lane-row(s) with an undefined bound_fraction.")
    combined_data <- combined_data %>% filter(!is.na(.data$bound_fraction))
}

# ------------------------------------------------------------------------------
# Screen composition per (genotype, ATP): the CHECK that folding did what we think.
# It reports, for each cell, which screens contributed and how many gels each gave,
# so WT should show both screens (3 + 3 = 6) and ORC4R only one (3). This is written
# out AND a soft message is printed; it does not stop the run (an asymmetric design
# is legitimate), but it makes the fold auditable instead of trusted.
# ------------------------------------------------------------------------------
screen_composition <- combined_data %>%
    group_by(.data$genotype, .data$atp, .data$screen) %>%
    summarise(gel_count = n(), .groups = "drop") %>%
    arrange(.data$genotype, .data$atp, .data$screen)
write.csv(screen_composition,
          file.path(output_directory, "combined_screen_composition.csv"),
          row.names = FALSE)
message("wrote ", file.path(output_directory, "combined_screen_composition.csv"))
folded_genotypes <- screen_composition %>%
    group_by(.data$genotype, .data$atp) %>%
    summarise(screen_count = n_distinct(.data$screen),
              total_gels = sum(.data$gel_count), .groups = "drop") %>%
    filter(.data$screen_count > 1)
if (nrow(folded_genotypes) > 0) {
    for (folded_row_index in seq_len(nrow(folded_genotypes))) {
        message("folded across screens: ",
                folded_genotypes$genotype[[folded_row_index]], " / ",
                folded_genotypes$atp[[folded_row_index]], " -> n=",
                folded_genotypes$total_gels[[folded_row_index]],
                " from ", folded_genotypes$screen_count[[folded_row_index]], " screens")
    }
}

# ------------------------------------------------------------------------------
# Split No ORC out as the baseline (pooled across screens), same as the per-screen
# plots. It is not an ATP-contrast bar.
# ------------------------------------------------------------------------------
no_orc_data <- combined_data %>% filter(.data$genotype == "No ORC")
no_orc_baseline <- if (nrow(no_orc_data) > 0) mean(no_orc_data$bound_fraction) else NA_real_
if (!is.na(no_orc_baseline)) {
    message("No ORC baseline (pooled n=", nrow(no_orc_data), "): ",
            round(no_orc_baseline, 4))
}

analyte_data <- combined_data %>%
    filter(.data$genotype != "No ORC")

# Loud post-factor guard: a mapped label absent from the level list would drop to
# NA at factor(). No ORC already removed.
labels_dropped_by_factoring <- setdiff(
    unique(analyte_data$genotype[!is.na(analyte_data$genotype)]),
    GENOTYPE_LEVELS_IN_ORDER)
if (length(labels_dropped_by_factoring) > 0) {
    stop("genotype label(s) not in GENOTYPE_LEVELS_IN_ORDER, would drop to NA at ",
         "factor(): ", paste(labels_dropped_by_factoring, collapse = ", "),
         "\nAdd them to GENOTYPE_LEVELS_IN_ORDER and fill_colors.", call. = FALSE)
}
analyte_data <- analyte_data %>%
    mutate(
        genotype = factor(.data$genotype, levels = GENOTYPE_LEVELS_IN_ORDER),
        atp = factor(.data$atp, levels = ATP_LEVELS_IN_ORDER),
        screen = factor(.data$screen, levels = screen_labels)
    )

# ------------------------------------------------------------------------------
# Pooled summary per (genotype, ATP): the fold. Mean and SD are over the POOLED
# replicate bound_fraction values, so WT's n=6 SD reflects all six gels (including
# between-screen spread), which is exactly the point of folding.
# ------------------------------------------------------------------------------
summary_data <- analyte_data %>%
    group_by(.data$genotype, .data$atp) %>%
    summarise(
        replicate_count = n(),
        mean_bound_fraction = mean(.data$bound_fraction),
        sd_bound_fraction = sd(.data$bound_fraction),
        cv_bound_fraction = if_else(mean(.data$bound_fraction) != 0,
                                    sd(.data$bound_fraction) / abs(mean(.data$bound_fraction)),
                                    NA_real_),
        se_bound_fraction = sd(.data$bound_fraction) / sqrt(n()),
        .groups = "drop") %>%
    arrange(.data$genotype, .data$atp)
write.csv(summary_data, file.path(output_directory, "combined_summary.csv"),
          row.names = FALSE)
message("wrote ", file.path(output_directory, "combined_summary.csv"))

# Assert the fold produced exactly one bar per (genotype, ATP): screen must be a
# point aesthetic only and must NOT have split a genotype into two bars. If this
# fires, screen leaked into the bar grouping somewhere above.
bars_per_cell <- summary_data %>%
    group_by(.data$genotype, .data$atp) %>%
    summarise(row_count = n(), .groups = "drop")
if (any(bars_per_cell$row_count != 1)) {
    stop("a (genotype, ATP) cell produced more than one summary row; screen leaked ",
         "into the bar grouping.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# A house theme matching the per-screen plots.
# ------------------------------------------------------------------------------
house_theme <- theme_classic(base_size = PLOT_CONFIG$theme$base_size,
                             base_family = PLOT_CONFIG$theme$base_family) +
    theme(strip.background = element_rect(fill = "gray90", color = "black"),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = PLOT_CONFIG$theme$legend_position,
          legend.box = "vertical",
          panel.spacing = unit(1, "lines"))

# ------------------------------------------------------------------------------
# Fill = genotype x ATP: each genotype keeps its palette hue, and noATP is that hue
# lightened toward white so the ATP state reads by shade. DUPLICATED from the
# per-screen P1 (house rule: no shared module). lighten_colour mixes toward white by
# a fraction using base grDevices only.
# ------------------------------------------------------------------------------
lighten_colour <- function(hex_colour, amount) {
    rgb_fraction <- col2rgb(hex_colour) / 255
    lightened_fraction <- rgb_fraction + (1 - rgb_fraction) * amount
    rgb(lightened_fraction[1], lightened_fraction[2], lightened_fraction[3])
}
genotype_atp_fill_values <- c()
genotype_atp_fill_labels <- c()
for (genotype_name in GENOTYPE_LEVELS_IN_ORDER) {
    full_colour <- PLOT_CONFIG$fill_colors[[genotype_name]]
    plus_key <- paste0(genotype_name, ".plusATP")
    no_key <- paste0(genotype_name, ".noATP")
    genotype_atp_fill_values[[plus_key]] <- full_colour
    genotype_atp_fill_values[[no_key]] <- lighten_colour(
        full_colour, PLOT_CONFIG$noatp_lighten_amount)
    genotype_atp_fill_labels[[plus_key]] <- paste0(genotype_name, " +ATP")
    genotype_atp_fill_labels[[no_key]] <- paste0(genotype_name, " -ATP")
}
# Level order: genotype outer, ATP inner with noATP first, so within each genotype
# the lighter noATP bar dodges to the left of the fuller plusATP bar -- the adjacent
# ATP pair the professor asked for.
genotype_atp_level_order <- unlist(lapply(
    GENOTYPE_LEVELS_IN_ORDER,
    function(genotype_name) c(paste0(genotype_name, ".noATP"),
                             paste0(genotype_name, ".plusATP"))))
summary_data <- summary_data %>%
    mutate(genotype_atp = factor(
        paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
        levels = genotype_atp_level_order))
analyte_data <- analyte_data %>%
    mutate(genotype_atp = factor(
        paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
        levels = genotype_atp_level_order))

# ==============================================================================
# PRIMARY OUTPUT: the dodged P1 panel. x = genotype; within each genotype the noATP
# and plusATP bars sit adjacent for direct contrast; fill = genotype hue with ATP as
# shade (noATP lighter, plusATP full); SD error bars; replicate dots. WT and No ORC
# are folded across screens (one bar each). Emitted in TWO point-encoding variants
# so you can choose:
#   _shape_rep_color_screen : shape = replicate, colour = screen (keeps provenance,
#                             so the folded WT dots visibly come from both screens)
#   _shape_rep_plain        : plain grey dots, shape = replicate, no screen colour
# Everything else about the two figures is identical.
# ==============================================================================
dodge_width <- 0.8

# The bars and error bars are identical across both variants; only the geom_point
# layer differs. Build the shared base once, then add the two point layers.
combined_dodged_base <- ggplot(
    summary_data,
    aes(x = .data$genotype, y = .data$mean_bound_fraction,
        fill = .data$genotype_atp, group = .data$atp)) +
    geom_col(position = position_dodge(width = dodge_width),
             width = PLOT_CONFIG$bar$width, color = PLOT_CONFIG$bar$color,
             linewidth = PLOT_CONFIG$bar$linewidth) +
    geom_errorbar(
        aes(ymin = pmax(0, .data$mean_bound_fraction - .data$sd_bound_fraction),
            ymax = .data$mean_bound_fraction + .data$sd_bound_fraction),
        position = position_dodge(width = dodge_width),
        width = PLOT_CONFIG$errorbar$width,
        linewidth = PLOT_CONFIG$errorbar$linewidth) +
    scale_fill_manual(values = genotype_atp_fill_values,
                      labels = genotype_atp_fill_labels, name = "Genotype / ATP",
                      breaks = genotype_atp_level_order) +
    scale_y_continuous(limits = PLOT_CONFIG$y_axis_limits,
                       breaks = PLOT_CONFIG$y_axis_breaks,
                       expand = PLOT_CONFIG$y_axis_expand) +
    labs(x = "Genotype", y = "Bound fraction  [ bound / (bound + free) ]",
         title = "Gel shift bound fraction by genotype and ATP (WT folded across screens)")

if (!is.na(no_orc_baseline)) {
    combined_dodged_base <- combined_dodged_base +
        geom_hline(yintercept = no_orc_baseline, linetype = "dashed",
                   color = "grey40", linewidth = 0.5) +
        annotate("text", x = 0.6, y = no_orc_baseline, vjust = -0.4, hjust = 0,
                 label = "No ORC", size = 3, color = "grey40")
}

# Variant A: shape = replicate, colour = screen.
plot_dodged_screen <- combined_dodged_base +
    geom_point(
        data = analyte_data,
        aes(x = .data$genotype, y = .data$bound_fraction, group = .data$atp,
            shape = factor(.data$replicate), color = .data$screen),
        position = position_jitterdodge(
            jitter.width = PLOT_CONFIG$point$jitter_width, dodge.width = dodge_width,
            seed = PLOT_CONFIG$point$jitter_seed),
        size = PLOT_CONFIG$point$size, stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE) +
    scale_shape_manual(values = PLOT_CONFIG$replicate_shapes,
                       labels = c("1" = "Replicate 1", "2" = "Replicate 2",
                                  "3" = "Replicate 3"), name = "Replicate") +
    scale_color_manual(values = PLOT_CONFIG$screen_colors, name = "Screen",
                       drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, order = 1),
           shape = guide_legend(nrow = 1, order = 2),
           color = guide_legend(nrow = 1, order = 3)) +
    house_theme
plot_dodged_screen_path <- file.path(
    output_directory, "combined_dodged_ATP_shape_rep_color_screen.pdf")
ggsave(plot_dodged_screen_path, plot = plot_dodged_screen,
       device = PLOT_CONFIG$output$device,
       width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
message("wrote ", plot_dodged_screen_path, " (PRIMARY)")

# Variant B: plain grey dots, shape = replicate, no screen colour.
plot_dodged_plain <- combined_dodged_base +
    geom_point(
        data = analyte_data,
        aes(x = .data$genotype, y = .data$bound_fraction, group = .data$atp,
            shape = factor(.data$replicate)),
        position = position_jitterdodge(
            jitter.width = PLOT_CONFIG$point$jitter_width, dodge.width = dodge_width,
            seed = PLOT_CONFIG$point$jitter_seed),
        size = PLOT_CONFIG$point$size, fill = "grey30",
        color = PLOT_CONFIG$point$color, stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE) +
    scale_shape_manual(values = PLOT_CONFIG$replicate_shapes,
                       labels = c("1" = "Replicate 1", "2" = "Replicate 2",
                                  "3" = "Replicate 3"), name = "Replicate") +
    guides(fill = guide_legend(nrow = 2, order = 1),
           shape = guide_legend(nrow = 1, order = 2)) +
    house_theme
plot_dodged_plain_path <- file.path(
    output_directory, "combined_dodged_ATP_shape_rep_plain.pdf")
ggsave(plot_dodged_plain_path, plot = plot_dodged_plain,
       device = PLOT_CONFIG$output$device,
       width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
message("wrote ", plot_dodged_plain_path, " (PRIMARY, plain-dot variant)")

# ==============================================================================
# NON-PRIMARY OUTPUTS below. The dodged panel above is the primary figure; the two
# single-ATP files and the ATP facet are kept only as alternate views and can be
# ignored for the main figure. They are not maintained to the same standard.
# ==============================================================================

# ------------------------------------------------------------------------------
# NON-PRIMARY: the two single-ATP figures (one bar per genotype, flat hue).
# ------------------------------------------------------------------------------
for (atp_state in ATP_LEVELS_IN_ORDER) {
    summary_this_atp <- summary_data %>% filter(.data$atp == atp_state)
    points_this_atp <- analyte_data %>% filter(.data$atp == atp_state)
    if (nrow(summary_this_atp) == 0) {
        message("no data for ATP state ", atp_state, "; skipping its figure.")
        next
    }
    single_atp_plot <- ggplot(
        summary_this_atp,
        aes(x = .data$genotype, y = .data$mean_bound_fraction, fill = .data$genotype)) +
        geom_col(width = PLOT_CONFIG$bar$width, color = PLOT_CONFIG$bar$color,
                 linewidth = PLOT_CONFIG$bar$linewidth) +
        geom_errorbar(
            aes(ymin = pmax(0, .data$mean_bound_fraction - .data$sd_bound_fraction),
                ymax = .data$mean_bound_fraction + .data$sd_bound_fraction),
            width = PLOT_CONFIG$errorbar$width, linewidth = PLOT_CONFIG$errorbar$linewidth) +
        geom_point(
            data = points_this_atp,
            aes(x = .data$genotype, y = .data$bound_fraction, color = .data$screen),
            position = position_jitter(width = PLOT_CONFIG$point$jitter_width, height = 0,
                                       seed = PLOT_CONFIG$point$jitter_seed),
            size = PLOT_CONFIG$point$size, stroke = PLOT_CONFIG$point$stroke,
            inherit.aes = FALSE) +
        scale_fill_manual(values = PLOT_CONFIG$fill_colors,
                          breaks = GENOTYPE_LEVELS_IN_ORDER, name = "Genotype",
                          guide = "none") +
        scale_color_manual(values = PLOT_CONFIG$screen_colors, name = "Screen",
                           drop = FALSE) +
        scale_y_continuous(limits = PLOT_CONFIG$y_axis_limits,
                           breaks = PLOT_CONFIG$y_axis_breaks,
                           expand = PLOT_CONFIG$y_axis_expand) +
        labs(x = "Genotype", y = "Bound fraction",
             title = paste0("Gel shift bound fraction (", atp_state,
                            ", WT folded across screens) [non-primary]")) +
        house_theme
    if (!is.na(no_orc_baseline)) {
        single_atp_plot <- single_atp_plot +
            geom_hline(yintercept = no_orc_baseline, linetype = "dashed",
                       color = "grey40", linewidth = 0.5) +
            annotate("text", x = 0.6, y = no_orc_baseline, vjust = -0.4, hjust = 0,
                     label = "No ORC", size = 3, color = "grey40")
    }
    single_atp_path <- file.path(output_directory,
                                 paste0("combined_", atp_state, ".pdf"))
    ggsave(single_atp_path, plot = single_atp_plot, device = PLOT_CONFIG$output$device,
           width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
    message("wrote ", single_atp_path, " (non-primary)")
}

# ------------------------------------------------------------------------------
# NON-PRIMARY: the ATP-faceted figure. Superseded by the dodged panel above, which
# puts the ATP pair adjacent; kept only as an alternate view.
# ------------------------------------------------------------------------------
faceted_plot <- ggplot(
    summary_data,
    aes(x = .data$genotype, y = .data$mean_bound_fraction, fill = .data$genotype)) +
    geom_col(width = PLOT_CONFIG$bar$width, color = PLOT_CONFIG$bar$color,
             linewidth = PLOT_CONFIG$bar$linewidth) +
    geom_errorbar(
        aes(ymin = pmax(0, .data$mean_bound_fraction - .data$sd_bound_fraction),
            ymax = .data$mean_bound_fraction + .data$sd_bound_fraction),
        width = PLOT_CONFIG$errorbar$width, linewidth = PLOT_CONFIG$errorbar$linewidth) +
    geom_point(
        data = analyte_data,
        aes(x = .data$genotype, y = .data$bound_fraction, color = .data$screen),
        position = position_jitter(width = PLOT_CONFIG$point$jitter_width, height = 0,
                                   seed = PLOT_CONFIG$point$jitter_seed),
        size = PLOT_CONFIG$point$size, stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE) +
    scale_fill_manual(values = PLOT_CONFIG$fill_colors,
                      breaks = GENOTYPE_LEVELS_IN_ORDER, name = "Genotype",
                      guide = "none") +
    scale_color_manual(values = PLOT_CONFIG$screen_colors, name = "Screen",
                       drop = FALSE) +
    scale_y_continuous(limits = PLOT_CONFIG$y_axis_limits,
                       breaks = PLOT_CONFIG$y_axis_breaks,
                       expand = PLOT_CONFIG$y_axis_expand) +
    facet_wrap(~ .data$atp, nrow = 1) +
    labs(x = "Genotype", y = "Bound fraction",
         title = "Gel shift bound fraction by genotype and ATP (WT folded) [non-primary]") +
    house_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (!is.na(no_orc_baseline)) {
    faceted_plot <- faceted_plot +
        geom_hline(yintercept = no_orc_baseline, linetype = "dashed",
                   color = "grey40", linewidth = 0.5)
}
faceted_path <- file.path(output_directory, "combined_faceted_by_ATP.pdf")
ggsave(faceted_path, plot = faceted_plot, device = PLOT_CONFIG$output$device,
       width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
message("wrote ", faceted_path, " (non-primary)")

message("done.")
