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
OUT_DIR <- "~/artifact_visualizations"
dir.create(path.expand(OUT_DIR), showWarnings = FALSE, recursive = TRUE)

# Samples to visualize (pick representative ORC, MCM, Input)
ORC_SAMPLE <- "D25-12504"
MCM_SAMPLE <- "D25-12510"
INPUT_SAMPLE <- "D25-12496"

# Genomic context window around artifacts (bp)
CONTEXT_WINDOW <- 20000

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

message("\n=================================================================")
message("Visualization Complete")
message("=================================================================")
message(paste("Output directory:", path.expand(OUT_DIR)))
message(paste("Generated", nrow(artifacts), "PDF files"))
