#!/usr/bin/env Rscript
# Visualize Consensus Artifacts with Gviz
# Shows IP/Input tracks, GFF features, and artifact regions

library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BW_DIR <- "/home/luised94/data/250930Bel/coverage"
ARTIFACT_FILE <- "~/consensus_artifacts/consensus_artifacts_annotated.txt"
GFF_FILE <- "/home/luised94/data/feature_files/240830_saccharomyces_cerevisiae.gff"
GENOME_FILE <- "/home/luised94/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
OUT_DIR <- "~/artifact_visualizations"
dir.create(path.expand(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
REPEAT_ZOOM_IN <- FALSE

# Samples to visualize (pick representative ORC, MCM, Input)
ORC_SAMPLE <- "D25-12504"
MCM_SAMPLE <- "D25-12510"
INPUT_SAMPLE <- "D25-12496"

# Genomic context window around artifacts (bp)
CONTEXT_WINDOW <- 20000

# Chromosome lengths (S288C R64)
CHR_LENGTHS <- c(
  chrI = 230218,
  chrII = 813184,
  chrIII = 316620,
  chrIV = 1531933,
  chrV = 576874,
  chrVI = 270161,
  chrVII = 1090940,
  chrVIII = 562643,
  chrIX = 439888,
  chrX = 745751,
  chrXI = 666816,
  chrXII = 1078177,
  chrXIII = 924431,
  chrXIV = 784333,
  chrXV = 1091291,
  chrXVI = 948066
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading artifacts...")
artifacts <- read.table(path.expand(ARTIFACT_FILE), sep = "\t", header = TRUE, 
                       stringsAsFactors = FALSE)

message(paste("Found", nrow(artifacts), "artifacts to visualize"))

# Load GFF once
message("Loading GFF annotations...")
GENOME_FEATURES <- import(GFF_FILE, format = "GFF")
message(paste("Loaded", length(GENOME_FEATURES), "genomic features"))

# ==============================================================================
# CREATE VISUALIZATIONS
# ==============================================================================

if (REPEAT_ZOOM_IN) {

for (i in 1:nrow(artifacts)) {
  artifact <- artifacts[i, ]

  current_chromosome <- artifact$chr
  artifact_start <- artifact$start
  artifact_end <- artifact$end

  # Define region to visualize
  region_start <- max(1, artifact_start - CONTEXT_WINDOW)
  region_end <- artifact_end + CONTEXT_WINDOW

  message(sprintf("\nVisualizing artifact %d/%d: %s:%d-%d", 
                 i, nrow(artifacts), current_chromosome, artifact_start, artifact_end))

  # Create GRanges for this region
  current_genome_range <- GRanges(
    seqnames = current_chromosome,
    ranges = IRanges(start = region_start, end = region_end)
  )

  # Initialize track container
  track_container <- list()

  # 1. Genome axis track
  track_container[[1]] <- GenomeAxisTrack(
    name = sprintf("%s Axis", current_chromosome),
    fontcolor.title = "black",
    cex.title = 0.7,
    background.title = "white"
  )

  # 2. Artifact region highlight
  artifact_gr <- GRanges(
    seqnames = current_chromosome,
    ranges = IRanges(start = artifact_start, end = artifact_end)
  )

  track_container[[2]] <- AnnotationTrack(
    artifact_gr,
    name = "Artifact",
    fill = "pink",
    col = "red",
    alpha = 0.5,
    background.title = "lightcoral",
    fontcolor.title = "black",
    cex.title = 0.7
  )

  # 3. Load and add bigwig tracks
  # ORC track
  orc_file <- file.path(BW_DIR, paste0(ORC_SAMPLE, "_CPM.bw"))
  if (file.exists(orc_file)) {
    message("  Loading ORC bigwig...")
    orc_data <- import(
      orc_file,
      format = "BigWig",
      which = current_genome_range
    )

    if (length(orc_data) > 0) {
      track_container[[length(track_container) + 1]] <- DataTrack(
        range = orc_data,
        name = "ORC",
        showAxis = TRUE,
        showTitle = TRUE,
        type = "histogram",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        cex.title = 0.7,
        col = "darkgreen",
        fill = "darkgreen"
      )
    } else {
      message("  WARNING: No ORC data in this region")
    }
  }

  # MCM track
  mcm_file <- file.path(BW_DIR, paste0(MCM_SAMPLE, "_CPM.bw"))
  if (file.exists(mcm_file)) {
    message("  Loading MCM bigwig...")
    mcm_data <- import(
      mcm_file,
      format = "BigWig",
      which = current_genome_range
    )

    if (length(mcm_data) > 0) {
      track_container[[length(track_container) + 1]] <- DataTrack(
        range = mcm_data,
        name = "MCM",
        showAxis = TRUE,
        showTitle = TRUE,
        type = "histogram",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        cex.title = 0.7,
        col = "purple",
        fill = "purple"
      )
    } else {
      message("  WARNING: No MCM data in this region")
    }
  }

  # Input track
  input_file <- file.path(BW_DIR, paste0(INPUT_SAMPLE, "_CPM.bw"))
  if (file.exists(input_file)) {
    message("  Loading Input bigwig...")
    input_data <- import(
      input_file,
      format = "BigWig",
      which = current_genome_range
    )

    if (length(input_data) > 0) {
      track_container[[length(track_container) + 1]] <- DataTrack(
        range = input_data,
        name = "Input",
        showAxis = TRUE,
        showTitle = TRUE,
        type = "histogram",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        cex.title = 0.7,
        col = "gray60",
        fill = "gray60"
      )
    } else {
      message("  WARNING: No Input data in this region")
    }
  }

  # 4. Add genomic features
  message("  Adding genomic features...")
  current_genome_feature <- keepSeqlevels(
    GENOME_FEATURES,
    current_chromosome,
    pruning.mode = "coarse"
  )

  # Subset to region
  feature_subset <- subsetByOverlaps(current_genome_feature, current_genome_range)

  if (length(feature_subset) > 0) {
    track_container[[length(track_container) + 1]] <- AnnotationTrack(
      feature_subset,
      name = "Features",
      size = 1.0,
      background.title = "lightgray",
      fontcolor.title = "black",
      showAxis = FALSE,
      cex.title = 0.6,
      fill = "#8b4513",
      col = "#8b4513",
      stacking = "squish"
    )
  }

  # Create filename
  features_short <- gsub("[^A-Za-z0-9_]", "_", substr(artifact$overlapping_features, 1, 40))
  filename <- sprintf("%s_%d-%d_%s.pdf", 
                     gsub("chr", "", current_chromosome), 
                     artifact_start, 
                     artifact_end,
                     features_short)

  plot_output_file <- file.path(path.expand(OUT_DIR), filename)

  # Plot
  message("  Creating PDF...")

  tryCatch({
    pdf(
      file = plot_output_file,
      width = 10,
      height = 8,
      bg = "white",
      compress = TRUE,
      colormodel = "srgb",
      useDingbats = FALSE
    )

    plotTracks(
      trackList = track_container,
      chromosome = current_chromosome,
      from = region_start,
      to = region_end,
      margin = 15,
      innerMargin = 5,
      main = sprintf("Artifact: %s:%d-%d\n%s", 
                    current_chromosome, artifact_start, artifact_end,
                    artifact$overlapping_features),
      col.axis = "black",
      cex.axis = 0.8,
      cex.main = 0.7,
      fontface.main = 1,
      background.panel = "white"
    )

    dev.off()
    message(paste("  SUCCESS! Saved:", filename))

  }, error = function(e) {
    dev.off()
    message(paste("  ERROR:", e$message))
  })
}

} # end skip for loop

# ==============================================================================
# CREATE CHROMOSOME-WIDE OVERVIEW PLOTS
# ==============================================================================

message("\n=== Creating Chromosome-Wide Overview Plots ===")
# 1. Define your groupings as a named list.
#    The name is the final category. The value is the regex pattern.
grouped_feature_map <- list(
  "ARS"                = "ARS",
  "Centromere"         = "centromere",
  "Telomere_Associated"= "telomere|X_element|Y_prime",
  "tRNA"               = "tRNA",
  "LTR_Elements"       = "LTR|Ty|retrotransposon",
  "rRNA"               = "rRNA"
)

# 2. Derive the final group names and create the color map.
#    This is now driven directly by your mapping list.
final_group_names <- names(grouped_feature_map)
library(RColorBrewer)
# Add a neutral color for the 'Other' category if we want to show it
number_of_colors <- length(final_group_names)
final_colors <- brewer.pal(min(number_of_colors, 8), "Dark2") # Using 'Dark2' as requested
final_color_map <- setNames(final_colors, final_group_names)

# Get unique chromosomes with artifacts
artifact_chromosomes <- unique(artifacts$chr)

for (current_chromosome in artifact_chromosomes) {
  message(sprintf("\nCreating overview for %s...", current_chromosome))

  # Get all artifacts on this chromosome
  chr_artifacts <- artifacts[artifacts$chr == current_chromosome, ]

  # Create GRanges for all artifacts on this chromosome
  artifact_regions <- GRanges(
    seqnames = current_chromosome,
    ranges = IRanges(
      start = chr_artifacts$start,
      end = chr_artifacts$end
    ),
    name = paste0("Artifact_", seq_len(nrow(chr_artifacts)))
  )

  # Initialize track container
  track_container <- list()

  # 1. Genome axis
   genome_axis_track <- GenomeAxisTrack(
    name = sprintf("%s Position", current_chromosome),
    fontcolor.title = "black",
    cex.title = 0.7,
    background.title = "white"
  )

  # 2. Artifact regions (highlighted)
   artifact_annotation_track <- AnnotationTrack(
    artifact_regions,
    name = "Artifacts",
    fill = "red",
    col = "darkred",
    alpha = 0.4,
    background.title = "lightcoral",
    fontcolor.title = "black",
    cex.title = 0.7,
    stacking = "squish"
  )

  # --- Load datatracks ---
  # 3. Load full chromosome bigwig data
  # ORC track
  orc_file <- file.path(BW_DIR, paste0(ORC_SAMPLE, "_CPM.bw"))
  if (file.exists(orc_file)) {
    message("  Loading ORC bigwig for full chromosome...")

    # Create GRanges for full chromosome (no subsetting)
    orc_data_full <- import(
      orc_file,
      format = "BigWig",
      as = "GRanges"
    )

    # Filter to current chromosome
    orc_data_chr <- orc_data_full[seqnames(orc_data_full) == current_chromosome]

    if (length(orc_data_chr) > 0) {
      track_container[[length(track_container) + 1]] <- DataTrack(
        range = orc_data_chr,
        name = "ORC",
        showAxis = TRUE,
        showTitle = TRUE,
        type = "histogram",
        size = 1.0,
        background.title = "white",
        fontcolor.title = "black",
        cex.title = 0.7,
        col = "darkgreen",
        fill = "darkgreen"
      )
    }
  }

  # MCM track
  mcm_file <- file.path(BW_DIR, paste0(MCM_SAMPLE, "_CPM.bw"))
  if (file.exists(mcm_file)) {
    message("  Loading MCM bigwig for full chromosome...")

    mcm_data_full <- import(
      mcm_file,
      format = "BigWig",
      as = "GRanges"
    )

    mcm_data_chr <- mcm_data_full[seqnames(mcm_data_full) == current_chromosome]

    if (length(mcm_data_chr) > 0) {
      track_container[[length(track_container) + 1]] <- DataTrack(
        range = mcm_data_chr,
        name = "MCM",
        showAxis = TRUE,
        showTitle = TRUE,
        type = "histogram",
        size = 1.0,
        background.title = "white",
        fontcolor.title = "black",
        cex.title = 0.7,
        col = "purple",
        fill = "purple"
      )
    }
  }

  # Input track
  input_file <- file.path(BW_DIR, paste0(INPUT_SAMPLE, "_CPM.bw"))
  if (file.exists(input_file)) {
    message("  Loading Input bigwig for full chromosome...")

    input_data_full <- import(
      input_file,
      format = "BigWig",
      as = "GRanges"
    )

    input_data_chr <- input_data_full[seqnames(input_data_full) == current_chromosome]

    if (length(input_data_chr) > 0) {
      track_container[[length(track_container) + 1]] <- DataTrack(
        range = input_data_chr,
        name = "Input",
        showAxis = TRUE,
        showTitle = TRUE,
        type = "histogram",
        size = 1.0,
        background.title = "white",
        fontcolor.title = "black",
        cex.title = 0.7,
        col = "gray60",
        fill = "gray60"
      )
    }
  }

  # 4. Add genomic features (optional - can be dense)
  message("  Adding genomic features...")
  current_genome_feature <- keepSeqlevels(
    GENOME_FEATURES,
    current_chromosome,
    pruning.mode = "coarse"
  )

  # Start with the data tracks as the base layer to be wrapped
  #highlight_wrapper <- track_container
  major_features <- current_genome_feature[GenomicRanges::width(current_genome_feature) > 0L]
  major_features <- major_features[grepl(feature_patterns_collapsed, major_features$type, ignore.case = TRUE)]
  GenomicRanges::strand(major_features) <- "*"

  # 1. Create a new metadata column for the standardized group name, initializing it.
  major_features$group <- "Other"

  # 2. Loop through your patterns and assign the correct group name.
  message("  Standardizing feature types into groups...")
  for (group_name in final_group_names) {
    # Get the regex pattern for the current group
    pattern <- grouped_feature_map[[group_name]]
    
    # Find features matching the pattern
    matching_indices <- grepl(pattern, major_features$type, ignore.case = TRUE)
    
    # Assign the clean group name to those features
    major_features$group[matching_indices] <- group_name
  }

  # Loop through each feature pattern to build nested highlight layers
  if (length(major_features) > 0 && length(track_container) > 0) {

        highlight_wrapper <- Gviz::HighlightTrack(
          trackList = track_container,
          range = major_features,
          chromosome = current_chromosome,
          inBackground = FALSE,
          alpha = 0.8,
          col = "transparent"
        )

  }

  region_colors <- final_color_map[major_features$group]
  region_colors[is.na(region_colors)] <- "#D3D3D3" # Assign a default gray for any NAs
  displayPars(highlight_wrapper) <- list(fill = region_colors)

  # Assemble the final list of tracks for plotting
  final_plot_list <- list(
    genome_axis_track,
    artifact_annotation_track,
    highlight_wrapper
  )

  # Create filename
  filename <- sprintf("%s_full_chromosome_overview.pdf", 
                     gsub("chr", "", current_chromosome))

  plot_output_file <- file.path(path.expand(OUT_DIR), filename)

  # Plot full chromosome
  message("  Creating chromosome-wide PDF...")

  tryCatch({
    pdf(
      file = plot_output_file,
      width = 16,  # Wider for full chromosome
      height = 8,
      bg = "white",
      compress = TRUE,
      colormodel = "srgb",
      useDingbats = FALSE
    )

    plotTracks(
      trackList = final_plot_list,
      chromosome = current_chromosome,
      margin = 15,
      innerMargin = 5,
      main = sprintf("%s - Full Chromosome with %d Artifacts Highlighted",
                    current_chromosome, nrow(chr_artifacts)),
      col.axis = "black",
      cex.axis = 0.7,
      cex.main = 0.8,
      fontface.main = 2,
      background.panel = "white"
    )

    dev.off()
    message(paste("  SUCCESS! Saved:", filename))

  }, error = function(e) {
    dev.off()
    message(paste("  ERROR:", e$message))
  })
}

# --- Code to generate a graphical legend PDF ---

# 1. Define the names and colors from your map
legend_labels <- names(final_color_map)
legend_colors <- unname(final_color_map)
legend_filename <- "legend_for_plots.pdf"

# 2. Open a PDF device
pdf(file = file.path(path.expand(OUT_DIR), legend_filename), width = 8, height = 6)

# 3. Create an empty plot. The 'plot.new' command creates a blank canvas.
plot.new()

# 4. Add the legend to the blank canvas
legend(
  "center",                     # Position the legend in the center
  legend = legend_labels,       # The text labels
  fill = legend_colors,         # The colors for the boxes
  title = "Feature Types",      # A title for the legend
  bty = "n",                    # No box around the legend (bty = "o" for a box)
  cex = 1.2                     # Adjust text size if needed
)

# 5. Close the PDF device, which saves the file
dev.off()

message("\nLegend has been saved to legend_for_plots.pdf")

message("\n=================================================================")
message("Visualization Complete")
message("=================================================================")
message(paste("Output directory:", path.expand(OUT_DIR)))
message(paste("Generated", nrow(artifacts), "zoomed PDF files"))
message(paste("Generated", length(artifact_chromosomes), "chromosome-wide overview PDF files"))
