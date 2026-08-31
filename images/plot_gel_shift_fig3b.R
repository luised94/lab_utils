# ==============================================================================
# plot_gel_shift_fig3c.R
#
# SELF-CONTAINED REPRODUCTION script for Figure 3C and its raw supplement.
# Run with no arguments, from anywhere:
#
#     Rscript plot_gel_shift_fig3c.R
#
# It resolves its inputs relative to ITS OWN location (not the current working
# directory), so it runs correctly however it is invoked, as long as the six
# ratio CSVs sit in the SAME directory as this script (flat layout, by design).
# The output directory is created beside the script on run.
#
# This is the trimmed, hardcoded companion to the exploratory
# plot_gel_shift_combined.R. That script takes CLI arguments and emits every
# exploratory view; THIS script hardcodes the six inputs and emits ONLY the
# figures that back the Figure 3C claims. Trimmed deliberately: no argument
# parsing, no manifests, no non-claim figures.
#
# ------------------------------------------------------------------------------
# ILLUSTRATOR / FIGURE-EDITING WARNING -- READ BEFORE EDITING ANY OUTPUT PDF
# ------------------------------------------------------------------------------
# The output PDFs are hand-edited in Illustrator for LABELS and ANNOTATION ONLY
# (text, arrows, panel letters, legend placement). Bar heights, error-bar
# lengths, point positions, axis limits, and the aspect ratio ARE DATA. They must
# NOT be rescaled, stretched, or reproportioned. Stretching one axis
# independently, or resizing a bar, silently misrepresents the underlying values.
# When resizing the whole panel, scale it UNIFORMLY (lock aspect ratio) so every
# encoded quantity keeps its meaning. If a value needs to change, change the data
# and re-run this script -- do not edit the geometry in Illustrator.
# ------------------------------------------------------------------------------
#
# What it produces (in ./fig3c_output/ beside this script):
#   fig3c_separate_gelsets_normalized_plusATP.pdf   the Figure 3C source: gel sets
#                                                   separate, each lane normalized
#                                                   to its own gel's WT plusATP
#                                                   (WT plusATP = 1.0 anchor).
#   fig3c_supplement_raw_combined.pdf               raw (unnormalized) combined
#                                                   view, WT folded, INCLUDING the
#                                                   low wt-1-3-4 rep-2 replicate --
#                                                   the honest "nothing hidden"
#                                                   supplement.
#   fig3c_summary_normalized_plusATP.csv            the numbers behind the bars.
#
# Determinism: fixed jitter seed and fixed factor orders; no RNG elsewhere. The
# same inputs always produce the same figure.
#
# Dependencies: tidyverse only.
# ==============================================================================

library(tidyverse)

# Print the editing warning to the console on EVERY run, so it is seen even by
# someone who never opens the source.
message("--------------------------------------------------------------------------------")
message("FIGURE-EDITING WARNING: output PDFs may be edited in Illustrator for labels and")
message("annotation ONLY. Bar heights, error bars, point positions, axis limits and aspect")
message("ratio are DATA -- do not rescale, stretch, or reproportion. Scale the panel")
message("uniformly (locked aspect) only. To change a value, change the data and re-run.")
message("--------------------------------------------------------------------------------")

# ------------------------------------------------------------------------------
# Resolve this script's own directory, so relative input paths do not depend on
# the caller's working directory. Rscript exposes the path via --file= in the
# command line; fall back to the current directory if run interactively.
# ------------------------------------------------------------------------------
command_arguments_full <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_arguments_full, value = TRUE)
if (length(file_argument) == 1) {
    script_path <- sub("^--file=", "", file_argument)
    script_directory <- dirname(normalizePath(script_path))
} else {
    # Interactive / sourced: assume the working directory holds the inputs.
    script_directory <- normalizePath(getwd())
}

# ------------------------------------------------------------------------------
# HARDCODED INPUTS. The six ratio CSVs, expected in the SAME directory as this
# script, tagged with screen label and within-screen replicate. This replaces the
# two manifests and the CLI arguments of the exploratory script. Screen labels
# must match PLOT_CONFIG$screen_colors below.
#
# The md5 checksums are the baseline this figure was produced from. On run, each
# input is checksummed and compared; a mismatch is WARNED (not fatal), because a
# legitimate re-export can change bytes without changing values -- but a reader
# then knows their data differs from the published baseline.
# ------------------------------------------------------------------------------
INPUT_FILES <- tibble::tribble(
    ~screen,      ~replicate, ~file_name,                                                        ~md5_baseline,
    "orc1-3-4",   1L,         "gel_shift_ratio_16.8-23.6mm_bound_over_46.6-56mm_free.csv",        "e6933abe17aed4c15dc9866f60620baa",
    "orc1-3-4",   2L,         "gel_shift_ratio_23.6-27.4mm_bound_over_38.9-45.6mm_free.csv",      "5d36d89d325f949062f24425ce863be7",
    "orc1-3-4",   3L,         "gel_shift_ratio_23.8-29.4mm_bound_over_39.8-45.6mm_free.csv",      "0a2a8026415cc007425d1d131f0604de",
    "orc5-6",     1L,         "gel_shift_ratio_16.8-21.4mm_bound_over_42.8-52.4mm_free.csv",      "de26db7e3c3617af43c717cb59ade943",
    "orc5-6",     2L,         "gel_shift_ratio_15.6-22.2mm_bound_over_43.2-54.6mm_free.csv",      "e3a44c91fa51466416ee4ef057449fcf",
    "orc5-6",     3L,         "gel_shift_ratio_19-23.8mm_bound_over_44.2-53.4mm_free.csv",        "e0e7ed49a210417b8dd936d2624d180a"
)

output_directory <- file.path(script_directory, "fig3c_output")
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
message("output directory: ", output_directory)

# Loud missing-input guard: name every input that is not found beside the script,
# so a stranger with the wrong layout gets a precise error, not a cryptic read
# failure deeper in.
INPUT_FILES$full_path <- file.path(script_directory, INPUT_FILES$file_name)
missing_inputs <- INPUT_FILES$file_name[!file.exists(INPUT_FILES$full_path)]
if (length(missing_inputs) > 0) {
    stop("expected these input CSV(s) in the same directory as the script (",
         script_directory, ") but they are missing:\n  ",
         paste(missing_inputs, collapse = "\n  "),
         "\nThis script uses a flat layout: the six ratio CSVs sit beside it.",
         call. = FALSE)
}

# Checksum verification (warn, do not stop).
for (input_index in seq_len(nrow(INPUT_FILES))) {
    actual_md5 <- unname(tools::md5sum(INPUT_FILES$full_path[input_index]))
    if (actual_md5 != INPUT_FILES$md5_baseline[input_index]) {
        warning("checksum mismatch for ", INPUT_FILES$file_name[input_index],
                ": baseline ", INPUT_FILES$md5_baseline[input_index],
                ", got ", actual_md5,
                ". Your data differs from the published Figure 3C baseline.",
                call. = FALSE)
    }
}

# ------------------------------------------------------------------------------
# Palette and level order, DUPLICATED verbatim from plot_gel_shift_combined.R
# (house rule: no shared module, so this repro script stands alone and drift
# fails loudly). Only the entries actually used by the two kept figures matter,
# but the full set is carried so the guards behave identically to the source.
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
    screen_colors = c(
        "orc1-3-4" = "#1B9E77",
        "orc5-6" = "#7570B3"
    ),
    bar = list(width = 0.7, color = "black", linewidth = 0.4),
    errorbar = list(width = 0.25, linewidth = 0.6),
    point = list(size = 2, color = "black", stroke = 0.5,
                 jitter_width = 0.15, jitter_seed = 42),
    noatp_lighten_amount = 0.55,
    replicate_shapes = c("1" = 21, "2" = 24, "3" = 22),
    y_axis_limits = c(0, 1.2),
    y_axis_breaks = seq(0, 1.2, 0.2),
    y_axis_expand = expansion(mult = c(0, 0), add = c(0, 0.02)),
    theme = list(base_family = "Arial", base_size = 12, legend_position = "bottom"),
    output = list(device = cairo_pdf, width = 8.5, height = 4.5)
)

GENOTYPE_LEVELS_IN_ORDER <- c(
    "WT",
    "ORC4R", "+5sofa", "+6sofa",
    "+1sofa", "+3sofa", "+4sofa"
)
ATP_LEVELS_IN_ORDER <- c("noATP", "plusATP")

# Screen colour guard, DUPLICATED from source: a screen without a configured
# point colour is a hard failure.
missing_screen_colours <- setdiff(unique(INPUT_FILES$screen),
                                  names(PLOT_CONFIG$screen_colors))
if (length(missing_screen_colours) > 0) {
    stop("no point colour configured for screen(s): ",
         paste(missing_screen_colours, collapse = ", "),
         " -- add them to PLOT_CONFIG$screen_colors.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Derivations, DUPLICATED verbatim from plot_gel_shift_combined.R.
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
# Read every gel, tag with screen and within-screen replicate. All columns read
# as character (two screens can type a shared column differently and bind_rows
# refuses to combine <double> with <character>); only bound_fraction is coerced
# numeric after binding. DUPLICATED logic from source.
# ------------------------------------------------------------------------------
per_gel_frames <- list()
for (input_index in seq_len(nrow(INPUT_FILES))) {
    one_gel <- read_csv(INPUT_FILES$full_path[input_index], show_col_types = FALSE,
                        col_types = cols(.default = col_character()))
    one_gel$screen <- INPUT_FILES$screen[input_index]
    one_gel$replicate <- as.character(INPUT_FILES$replicate[input_index])
    one_gel$screen_replicate <- paste0(INPUT_FILES$screen[input_index], ":",
                                       INPUT_FILES$replicate[input_index])
    per_gel_frames[[length(per_gel_frames) + 1]] <- one_gel
}
combined_data <- bind_rows(per_gel_frames)

combined_data <- combined_data %>%
    mutate(
        genotype = map2_chr(.data$`orc4-R267`, .data$suppressor, derive_genotype_label),
        atp = map_chr(.data$atp_presence, derive_atp_label),
        bound_fraction = as.numeric(.data$bound_fraction)
    )

# Loud guards, DUPLICATED from source: an unmapped genotype or ATP would silently
# become NA and pool.
unmapped_genotype_rows <- combined_data[is.na(combined_data$genotype), ]
if (nrow(unmapped_genotype_rows) > 0) {
    offending_combinations <- unique(paste0(
        "orc4-R267=", unmapped_genotype_rows$`orc4-R267`,
        ", suppressor=", unmapped_genotype_rows$suppressor))
    stop("derive_genotype_label produced NA for ", nrow(unmapped_genotype_rows),
         " lane-row(s); unmapped factor combination(s):\n  ",
         paste(offending_combinations, collapse = "\n  "), call. = FALSE)
}
unmapped_atp_rows <- combined_data[is.na(combined_data$atp), ]
if (nrow(unmapped_atp_rows) > 0) {
    stop("derive_atp_label produced NA for ", nrow(unmapped_atp_rows),
         " lane-row(s); unmapped atp_presence value(s): ",
         paste(unique(unmapped_atp_rows$atp_presence), collapse = ", "), call. = FALSE)
}
undefined_row_count <- sum(is.na(combined_data$bound_fraction))
if (undefined_row_count > 0) {
    message("dropping ", undefined_row_count,
            " lane-row(s) with an undefined bound_fraction.")
    combined_data <- combined_data %>% filter(!is.na(.data$bound_fraction))
}

# ------------------------------------------------------------------------------
# No ORC baseline (pooled) -- used only by the raw supplement figure, as in the
# source. Not drawn on the normalized axis.
# ------------------------------------------------------------------------------
no_orc_data <- combined_data %>% filter(.data$genotype == "No ORC")
no_orc_baseline <- if (nrow(no_orc_data) > 0) mean(no_orc_data$bound_fraction) else NA_real_

analyte_data <- combined_data %>% filter(.data$genotype != "No ORC")

labels_dropped_by_factoring <- setdiff(
    unique(analyte_data$genotype[!is.na(analyte_data$genotype)]),
    GENOTYPE_LEVELS_IN_ORDER)
if (length(labels_dropped_by_factoring) > 0) {
    stop("genotype label(s) not in GENOTYPE_LEVELS_IN_ORDER, would drop to NA at ",
         "factor(): ", paste(labels_dropped_by_factoring, collapse = ", "),
         call. = FALSE)
}
analyte_data <- analyte_data %>%
    mutate(
        genotype = factor(.data$genotype, levels = GENOTYPE_LEVELS_IN_ORDER),
        atp = factor(.data$atp, levels = ATP_LEVELS_IN_ORDER),
        screen = factor(.data$screen, levels = unique(INPUT_FILES$screen))
    )

# ------------------------------------------------------------------------------
# Shared theme, fill map, and level order, DUPLICATED from source.
# ------------------------------------------------------------------------------
house_theme <- theme_classic(base_size = PLOT_CONFIG$theme$base_size,
                             base_family = PLOT_CONFIG$theme$base_family) +
    theme(strip.background = element_rect(fill = "gray90", color = "black"),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = PLOT_CONFIG$theme$legend_position,
          legend.box = "vertical",
          panel.spacing = unit(1, "lines"))

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
genotype_atp_level_order <- unlist(lapply(
    GENOTYPE_LEVELS_IN_ORDER,
    function(genotype_name) c(paste0(genotype_name, ".noATP"),
                             paste0(genotype_name, ".plusATP"))))

dodge_width <- 0.8

# ==============================================================================
# RAW SUPPLEMENT: combined / folded, unnormalized, WT n=6 including the low
# wt-1-3-4 rep-2 replicate. This is the honest "nothing hidden" panel. One
# variant (shape = replicate, colour = screen) is kept; the exploratory script's
# other raw variants are intentionally omitted here.
# ==============================================================================
summary_data <- analyte_data %>%
    group_by(.data$genotype, .data$atp) %>%
    summarise(
        replicate_count = n(),
        mean_bound_fraction = mean(.data$bound_fraction),
        sd_bound_fraction = sd(.data$bound_fraction),
        .groups = "drop") %>%
    arrange(.data$genotype, .data$atp) %>%
    mutate(genotype_atp = factor(
        paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
        levels = genotype_atp_level_order))
analyte_data <- analyte_data %>%
    mutate(genotype_atp = factor(
        paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
        levels = genotype_atp_level_order))

raw_supplement_plot <- ggplot(
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
    geom_point(
        data = analyte_data,
        aes(x = .data$genotype, y = .data$bound_fraction, group = .data$atp,
            shape = factor(.data$replicate), color = .data$screen),
        position = position_jitterdodge(
            jitter.width = PLOT_CONFIG$point$jitter_width, dodge.width = dodge_width,
            seed = PLOT_CONFIG$point$jitter_seed),
        size = PLOT_CONFIG$point$size, stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE) +
    scale_fill_manual(values = genotype_atp_fill_values,
                      labels = genotype_atp_fill_labels, name = "Genotype / ATP",
                      breaks = genotype_atp_level_order) +
    scale_shape_manual(values = PLOT_CONFIG$replicate_shapes,
                       labels = c("1" = "Replicate 1", "2" = "Replicate 2",
                                  "3" = "Replicate 3"), name = "Replicate") +
    scale_color_manual(values = PLOT_CONFIG$screen_colors, name = "Screen",
                       drop = FALSE) +
    scale_y_continuous(limits = PLOT_CONFIG$y_axis_limits,
                       breaks = PLOT_CONFIG$y_axis_breaks,
                       expand = PLOT_CONFIG$y_axis_expand) +
    labs(x = "Genotype", y = "Bound fraction  [ bound / (bound + free) ]",
         title = "Raw bound fraction by genotype and ATP (WT folded; includes low rep) [supplement]") +
    guides(fill = guide_legend(nrow = 2, order = 1),
           shape = guide_legend(nrow = 1, order = 2),
           color = guide_legend(nrow = 1, order = 3)) +
    house_theme
if (!is.na(no_orc_baseline)) {
    raw_supplement_plot <- raw_supplement_plot +
        geom_hline(yintercept = no_orc_baseline, linetype = "dashed",
                   color = "grey40", linewidth = 0.5) +
        annotate("text", x = 0.6, y = no_orc_baseline, vjust = -0.4, hjust = 0,
                 label = "No ORC", size = 3, color = "grey40")
}
raw_supplement_path <- file.path(output_directory, "fig3c_supplement_raw_combined.pdf")
ggsave(raw_supplement_path, plot = raw_supplement_plot,
       device = PLOT_CONFIG$output$device,
       width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
message("wrote ", raw_supplement_path, " (raw supplement)")

# ==============================================================================
# FIGURE 3C: gel sets separate, each lane normalized to its own gel's WT PLUSATP
# (WT plusATP = 1.0 anchor per gel; WT noATP floats below, preserving the ATP
# contrast). Divisor lookup is (screen_replicate, WT, plusATP) applied to every
# lane regardless of its own ATP state.
# ==============================================================================
wt_plusatp_divisor_table <- analyte_data %>%
    filter(.data$genotype == "WT", .data$atp == "plusATP") %>%
    group_by(.data$screen_replicate) %>%
    summarise(wt_plusatp_bound_fraction = mean(.data$bound_fraction),
              wt_plusatp_row_count = n(), .groups = "drop")
if (any(wt_plusatp_divisor_table$wt_plusatp_row_count != 1)) {
    offending <- wt_plusatp_divisor_table %>% filter(.data$wt_plusatp_row_count != 1)
    stop("per-gel WT plusATP divisor is not exactly one lane for ",
         nrow(offending), " gel(s); e.g. ",
         paste(paste0(offending$screen_replicate,
                      " (n=", offending$wt_plusatp_row_count, ")"), collapse = ", "),
         ". Need exactly one WT plusATP lane per gel.", call. = FALSE)
}
WT_DIVISOR_FLOOR <- 0.02
tiny_wt_plusatp <- wt_plusatp_divisor_table %>%
    filter(.data$wt_plusatp_bound_fraction < WT_DIVISOR_FLOOR)
if (nrow(tiny_wt_plusatp) > 0) {
    stop("per-gel WT plusATP divisor below floor ", WT_DIVISOR_FLOOR, " for ",
         nrow(tiny_wt_plusatp), " gel(s): ",
         paste(paste0(tiny_wt_plusatp$screen_replicate, " (WT+ATP=",
                      round(tiny_wt_plusatp$wt_plusatp_bound_fraction, 4), ")"),
               collapse = ", "),
         ". A near-zero WT plusATP would produce explosive ratios.", call. = FALSE)
}

analyte_data_normalized_plusatp <- analyte_data %>%
    left_join(wt_plusatp_divisor_table %>% select(-.data$wt_plusatp_row_count),
              by = "screen_replicate") %>%
    mutate(normalized_bound_fraction =
               .data$bound_fraction / .data$wt_plusatp_bound_fraction)
unjoined_plusatp_count <- sum(is.na(analyte_data_normalized_plusatp$wt_plusatp_bound_fraction))
if (unjoined_plusatp_count > 0) {
    stop(unjoined_plusatp_count, " analyte lane-row(s) found no per-gel WT plusATP ",
         "divisor; a gel is missing its WT plusATP lane.", call. = FALSE)
}

# Data-driven y-limit with a 1.2 floor, rounded up to the next 0.2, so a value
# above 1.0 cannot clip silently. Reported on run.
normalized_plusatp_max_value <- max(analyte_data_normalized_plusatp$normalized_bound_fraction)
normalized_plusatp_y_upper <- max(1.2, ceiling(normalized_plusatp_max_value / 0.2) * 0.2)
message("Figure 3C normalized max value = ", round(normalized_plusatp_max_value, 4),
        "; y-limit set to ", normalized_plusatp_y_upper, ".")
normalized_plusatp_y_breaks <- seq(0, normalized_plusatp_y_upper, 0.2)

summary_data_normalized_plusatp_by_screen <- analyte_data_normalized_plusatp %>%
    group_by(.data$screen, .data$genotype, .data$atp) %>%
    summarise(
        replicate_count = n(),
        mean_normalized = mean(.data$normalized_bound_fraction),
        sd_normalized = sd(.data$normalized_bound_fraction),
        .groups = "drop") %>%
    arrange(.data$screen, .data$genotype, .data$atp) %>%
    mutate(genotype_atp = factor(
        paste0(as.character(.data$genotype), ".", as.character(.data$atp)),
        levels = genotype_atp_level_order))
summary_csv_path <- file.path(output_directory, "fig3c_summary_normalized_plusATP.csv")
write.csv(summary_data_normalized_plusatp_by_screen, summary_csv_path, row.names = FALSE)
message("wrote ", summary_csv_path)

fig3c_plot <- ggplot(
    summary_data_normalized_plusatp_by_screen,
    aes(x = .data$genotype, y = .data$mean_normalized,
        fill = .data$genotype_atp, group = .data$atp)) +
    geom_col(position = position_dodge(width = dodge_width),
             width = PLOT_CONFIG$bar$width, color = PLOT_CONFIG$bar$color,
             linewidth = PLOT_CONFIG$bar$linewidth) +
    geom_errorbar(
        aes(ymin = pmax(0, .data$mean_normalized - .data$sd_normalized),
            ymax = .data$mean_normalized + .data$sd_normalized),
        position = position_dodge(width = dodge_width),
        width = PLOT_CONFIG$errorbar$width,
        linewidth = PLOT_CONFIG$errorbar$linewidth) +
    geom_point(
        data = analyte_data_normalized_plusatp,
        aes(x = .data$genotype, y = .data$normalized_bound_fraction, group = .data$atp,
            shape = factor(.data$replicate)),
        position = position_jitterdodge(
            jitter.width = PLOT_CONFIG$point$jitter_width, dodge.width = dodge_width,
            seed = PLOT_CONFIG$point$jitter_seed),
        size = PLOT_CONFIG$point$size, fill = "grey30",
        color = PLOT_CONFIG$point$color, stroke = PLOT_CONFIG$point$stroke,
        inherit.aes = FALSE) +
    scale_fill_manual(values = genotype_atp_fill_values,
                      labels = genotype_atp_fill_labels, name = "Genotype / ATP",
                      breaks = genotype_atp_level_order) +
    scale_shape_manual(values = PLOT_CONFIG$replicate_shapes,
                       labels = c("1" = "Replicate 1", "2" = "Replicate 2",
                                  "3" = "Replicate 3"), name = "Replicate") +
    scale_y_continuous(limits = c(0, normalized_plusatp_y_upper),
                       breaks = normalized_plusatp_y_breaks,
                       expand = PLOT_CONFIG$y_axis_expand) +
    facet_wrap(~ .data$screen, nrow = 1) +
    labs(x = "Genotype", y = "Bound fraction, normalized  (per-gel WT +ATP = 1.0)",
         title = "Figure 3C: gel shift normalized to per-gel WT plusATP, gel sets separate") +
    guides(fill = guide_legend(nrow = 2, order = 1),
           shape = guide_legend(nrow = 1, order = 2)) +
    house_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig3c_path <- file.path(output_directory, "fig3c_separate_gelsets_normalized_plusATP.pdf")
ggsave(fig3c_path, plot = fig3c_plot, device = PLOT_CONFIG$output$device,
       width = PLOT_CONFIG$output$width, height = PLOT_CONFIG$output$height)
message("wrote ", fig3c_path, " (FIGURE 3C)")

message("done.")
