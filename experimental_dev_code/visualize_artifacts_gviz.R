#!/usr/bin/env Rscript
# Visualize Consensus Artifacts with Gviz
# Shows IP/Input tracks, GFF features, and artifact regions

library(Gviz)
library(rtracklayer)
library(GenomicRanges)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BW_DIR <- "/home/luised94/data/250930Bel/coverage"
ARTIFACT_FILE <- "~/consensus_artifacts/consensus_artifacts_annotated.txt"
GFF_FILE <- "/home/luised94/data/feature_files/240830_saccharomyces_cerevisiae.gff"
OUT_DIR <- "artifact_visualizations"
dir.create(OUT_DIR, showWarnings = FALSE)

# Samples to visualize (pick representative ORC, MCM, Input)
ORC_SAMPLE <- "D25-12504"  # Adjust to your actual sample IDs
MCM_SAMPLE <- "D25-12510"
INPUT_SAMPLE <- "D25-12496"

# Genomic context window around artifacts (bp)
CONTEXT_WINDOW <- 20000

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading artifacts...")
artifacts <- read.table(ARTIFACT_FILE, sep = "\t", header = TRUE, 
                       stringsAsFactors = FALSE)

message(paste("Found", nrow(artifacts), "artifacts to visualize"))

# Load GFF
message("Loading GFF annotations...")
gff <- import(GFF_FILE, format = "GFF")

# ==============================================================================
# CREATE VISUALIZATIONS
# ==============================================================================

for (i in 1:nrow(artifacts)) {
  artifact <- artifacts[i, ]
  
  chr <- artifact$chr
  start <- max(1, artifact$start - CONTEXT_WINDOW)
  end <- artifact$end + CONTEXT_WINDOW
  
  message(sprintf("\nVisualizing artifact %d/%d: %s:%d-%d", 
                 i, nrow(artifacts), chr, artifact$start, artifact$end))
  
  # Extract feature info
  features <- artifact$overlapping_features
  
  # Create filename
  filename <- sprintf("%s_%d-%d_%s.pdf", 
                     gsub("chr", "", chr), 
                     artifact$start, 
                     artifact$end,
                     gsub("[; ]", "_", substr(features, 1, 30)))
  
  pdf(file.path(OUT_DIR, filename), width = 12, height = 8)
  
  tryCatch({
    # Genome axis
    gtrack <- GenomeAxisTrack()
    
    # Ideogram
    itrack <- IdeogramTrack(genome = "sacCer3", chromosome = chr)
    
    # GFF features in region
    gff_region <- subsetByOverlaps(gff, GRanges(chr, IRanges(start, end)))
    
    if (length(gff_region) > 0) {
      # Create annotation track
      atrack <- AnnotationTrack(gff_region, 
                               name = "Features",
                               fill = "lightblue",
                               col = "darkblue",
                               stacking = "squish",
                               fontsize = 8)
    } else {
      atrack <- NULL
    }
    
    # Artifact region highlight
    artifact_gr <- GRanges(chr, IRanges(artifact$start, artifact$end))
    artifact_track <- AnnotationTrack(artifact_gr,
                                     name = "Artifact",
                                     fill = "red",
                                     col = "darkred",
                                     alpha = 0.3)
    
    # BigWig tracks
    orc_file <- file.path(BW_DIR, paste0(ORC_SAMPLE, "_CPM.bw"))
    mcm_file <- file.path(BW_DIR, paste0(MCM_SAMPLE, "_CPM.bw"))
    input_file <- file.path(BW_DIR, paste0(INPUT_SAMPLE, "_CPM.bw"))
    
    tracks <- list(itrack, gtrack, artifact_track)
    
    if (file.exists(orc_file)) {
      orc_track <- DataTrack(range = orc_file,
                            genome = "sacCer3",
                            chromosome = chr,
                            name = "ORC",
                            type = "histogram",
                            fill = "darkgreen",
                            col = "darkgreen")
      tracks <- c(tracks, list(orc_track))
    }
    
    if (file.exists(mcm_file)) {
      mcm_track <- DataTrack(range = mcm_file,
                            genome = "sacCer3",
                            chromosome = chr,
                            name = "MCM",
                            type = "histogram",
                            fill = "purple",
                            col = "purple")
      tracks <- c(tracks, list(mcm_track))
    }
    
    if (file.exists(input_file)) {
      input_track <- DataTrack(range = input_file,
                              genome = "sacCer3",
                              chromosome = chr,
                              name = "Input",
                              type = "histogram",
                              fill = "gray60",
                              col = "gray60")
      tracks <- c(tracks, list(input_track))
    }
    
    if (!is.null(atrack)) {
      tracks <- c(tracks, list(atrack))
    }
    
    # Plot
    plotTracks(tracks,
              from = start,
              to = end,
              chromosome = chr,
              main = sprintf("Artifact: %s\nFeatures: %s", 
                           artifact$window_id,
                           features),
              cex.main = 0.8)
    
  }, error = function(e) {
    message(paste("  Error plotting:", e$message))
  })
  
  dev.off()
  message(paste("  Saved:", filename))
}

message("\n=== Visualization Complete ===")
message(paste("Output directory:", OUT_DIR))
message(paste("Generated", nrow(artifacts), "PDF files"))
