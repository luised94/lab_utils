#!/usr/bin/env Rscript
# generate_color_legend.R
# Extract Chemistry_AA color scheme from ggmsa and generate a legend figure
#
# Approach: Use ggplot_build() to extract actual hex codes from a rendered
# ggmsa plot, ensuring we capture the exact colors the package uses.

library(ggmsa)
library(ggplot2)
library(Biostrings)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

OUTPUT_DIR <- "~/data/protein_files/alignments/plots"
OUTPUT_BASE <- "chemistry_aa_legend"

PLOT_WIDTH <- 8
PLOT_HEIGHT <- 6

SAVE_PNG <- TRUE
SAVE_SVG <- TRUE

# Chemical property groupings (standard biochemistry classification)
AA_GROUPS <- list(
    "Hydrophobic" = c("G", "A", "V", "L", "I", "P"),
    "Aromatic" = c("F", "W", "Y"),
    "Polar" = c("S", "T", "N", "Q"),
    "Sulfur" = c("C", "M"),
    "Acidic" = c("D", "E"),
    "Basic" = c("K", "R", "H")
)

# ==============================================================================
# SETUP
# ==============================================================================

OUTPUT_DIR <- path.expand(OUTPUT_DIR)

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(path = OUTPUT_DIR, recursive = TRUE)
}

cat("=== CHEMISTRY_AA LEGEND GENERATOR ===\n\n")


# ==============================================================================
# EXTRACT COLORS FROM GGMSA
# ==============================================================================

cat("Extracting colors from ggmsa...\n")

# Create a dummy sequence containing all 20 standard amino acids
aa_sequence <- "GAVLIPFWYSTNQCMDEKRH"
dummy_seq <- Biostrings::AAStringSet(aa_sequence)
names(dummy_seq) <- "dummy"

# Generate plot with Chemistry_AA scheme
p_dummy <- suppressWarnings(
    ggmsa::ggmsa(msa = dummy_seq, color = "Chemistry_AA", seq_name = FALSE)
)

# Extract built plot data
build <- ggplot2::ggplot_build(p_dummy)

# Layer 1 contains tile data with fill colors
tile_data <- build$data[[1]]

# Get unique x positions and their fill colors
color_by_position <- unique(tile_data[, c("x", "fill")])
color_by_position <- color_by_position[order(color_by_position$x), ]

# Map position to amino acid from our known sequence
aa_chars <- strsplit(aa_sequence, "")[[1]]

stopifnot(
    "Position count mismatch" = nrow(color_by_position) == length(aa_chars)
)

color_map <- data.frame(
    aa = aa_chars,
    color = color_by_position$fill,
    stringsAsFactors = FALSE
)

cat("  Extracted", nrow(color_map), "amino acid colors:\n")
for (i in seq_len(nrow(color_map))) {
    cat("    ", color_map$aa[i], " -> ", color_map$color[i], "\n", sep = "")
}
cat("\n")

# ==============================================================================
# BUILD LEGEND DATA
# ==============================================================================

cat("Building legend data...\n")

# Create data frame with group assignments
legend_data <- data.frame(
    aa = character(),
    group = character(),
    color = character(),
    stringsAsFactors = FALSE
)

for (group_name in names(AA_GROUPS)) {
    aas <- AA_GROUPS[[group_name]]
    for (aa in aas) {
        color_val <- color_map$color[color_map$aa == aa]
        if (length(color_val) == 1) {
            legend_data <- rbind(legend_data, data.frame(
                aa = aa,
                group = group_name,
                color = color_val,
                stringsAsFactors = FALSE
            ))
        } else {
            cat("  WARNING: No color found for", aa, "\n")
        }
    }
}

# Set factor levels for ordering
legend_data$group <- factor(
    x = legend_data$group,
    levels = names(AA_GROUPS)
)

# Order within groups
legend_data <- legend_data[order(legend_data$group, legend_data$aa), ]

cat("  Legend entries:", nrow(legend_data), "\n\n")

# ==============================================================================
# CREATE LEGEND FIGURE
# ==============================================================================

cat("Creating legend figure...\n")

# Calculate grid positions
legend_data$group_idx <- as.integer(legend_data$group)
legend_data$aa_idx <- ave(
    seq_len(nrow(legend_data)),
    legend_data$group,
    FUN = seq_along
)

p_legend <- ggplot2::ggplot(
    data = legend_data,
    mapping = ggplot2::aes(x = group, y = aa_idx)
) +
    ggplot2::geom_tile(
        mapping = ggplot2::aes(fill = color),
        color = "white",
        linewidth = 0.5,
        width = 0.9,
        height = 0.9
    ) +
    ggplot2::geom_text(
        mapping = ggplot2::aes(label = aa),
        size = 5,
        fontface = "bold"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
        title = "Chemistry_AA Color Scheme",
        subtitle = "Amino acids grouped by chemical properties",
        x = NULL,
        y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
    )


# ==============================================================================
# CREATE CLASSIFICATION FIGURE
# ==============================================================================

cat("Creating classification figure...\n")

# Map classification to representative color (from ggmsa extraction)
class_colors <- data.frame(
    group = c("Hydrophobic", "Aromatic", "Polar", "Acidic", "Basic"),
    color = c(
        color_map$color[color_map$aa == "A"],  # Hydrophobic
        color_map$color[color_map$aa == "F"],  # Aromatic
        color_map$color[color_map$aa == "S"],  # Polar
        color_map$color[color_map$aa == "D"],  # Acidic
        color_map$color[color_map$aa == "K"]   # Basic
    ),
    stringsAsFactors = FALSE
)

# Set factor for ordering
class_colors$group <- factor(
    x = class_colors$group,
    levels = class_colors$group
)

p_classification <- ggplot2::ggplot(
    data = class_colors,
    mapping = ggplot2::aes(x = 1, y = group)
) +
    ggplot2::geom_tile(
        mapping = ggplot2::aes(fill = color),
        color = "white",
        linewidth = 0.5,
        width = 0.6,
        height = 0.9
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_discrete(limits = rev(levels(class_colors$group))) +
    ggplot2::labs(
        title = "Chemistry_AA Color Scheme",
        x = NULL,
        y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
    )
# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

if (SAVE_PNG) {
    png_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_BASE, ".png"))
    ggplot2::ggsave(
        filename = png_path,
        plot = p_legend,
        width = PLOT_WIDTH,
        height = PLOT_HEIGHT,
        dpi = 300
    )
    cat("Saved PNG:", png_path, "\n")
}

if (SAVE_SVG) {
    svg_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_BASE, ".svg"))
    ggplot2::ggsave(
        filename = svg_path,
        plot = p_legend,
        width = PLOT_WIDTH,
        height = PLOT_HEIGHT
    )
    cat("Saved SVG:", svg_path, "\n")
}


# Save classification figure
if (SAVE_PNG) {
    png_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_BASE, "_classification.png"))
    ggplot2::ggsave(
        filename = png_path,
        plot = p_classification,
        width = 4,
        height = 3,
        dpi = 300
    )
    cat("Saved PNG:", png_path, "\n")
}

if (SAVE_SVG) {
    svg_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_BASE, "_classification.svg"))
    ggplot2::ggsave(
        filename = svg_path,
        plot = p_classification,
        width = 4,
        height = 3
    )
    cat("Saved SVG:", svg_path, "\n")
}

# ==============================================================================
# SAVE COLOR MAPPING TABLE
# ==============================================================================

color_table_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_BASE, "_colors.tsv"))
utils::write.table(
    x = legend_data[, c("group", "aa", "color")],
    file = color_table_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
cat("Saved color table:", color_table_path, "\n")

# ==============================================================================
# SESSION INFO
# ==============================================================================

session_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_BASE, "_session_info.txt"))
session_info <- utils::capture.output(utils::sessionInfo())
writeLines(text = session_info, con = session_path)
cat("Saved session info:", session_path, "\n")

cat("\n=== COMPLETE ===\n")
