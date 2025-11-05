#!/usr/bin/env Rscript
# Visualize Consensus Artifacts with Gviz
# Shows IP/Input tracks, GFF features, and artifact regions

library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Directories
BW_DIR <- "/home/luised94/data/250930Bel/coverage"
ARTIFACT_FILE <- "~/consensus_artifacts/consensus_artifacts_annotated.txt"
GFF_FILE <- "/home/luised94/data/feature_files/240830_saccharomyces_cerevisiae.gff"
GENOME_FILE <- "/home/luised94/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"
OUT_DIR <- "~/artifact_visualizations"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Samples to visualize (pick representative ORC, MCM, Input)
ORC_SAMPLE <- "D25-12504"  # Adjust to your actual sample IDs
MCM_SAMPLE <- "D25-12510"
INPUT_SAMPLE <- "D25-12496"

# Genomic context window around artifacts (bp)
CONTEXT_WINDOW <- 20000

# ==============================================================================
# LOAD GENOME AND GET CHROMOSOME LENGTHS
# ==============================================================================

message("Loading genome to get chromosome lengths...")

# Read chromosome lengths from FASTA index
fai_file <- paste0(GENOME_FILE, ".fai")

if (file.exists(fai_file)) {
  fai <- read.table(fai_file, header = FALSE, stringsAsFactors = FALSE)
  chr_lengths <- setNames(fai$V2, fai$V1)
  message(paste("  Loaded", length(chr_lengths), "chromosomes from FAI"))
} else {
  # Create FAI if it doesn't exist
  message("  Creating FASTA index...")
  genome <- readDNAStringSet(GENOME_FILE)
  chr_lengths <- setNames(width(genome), names(genome))
  message(paste("  Indexed", length(chr_lengths), "chromosomes"))
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("\nLoading artifacts...")
artifacts <- read.table(path.expand(ARTIFACT_FILE), sep = "\t", header = TRUE, 
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
  
  # Get chromosome length
  chr_len <- chr_lengths[chr]
  if (is.na(chr_len)) {
    message(sprintf("  WARNING: Chromosome %s not found in genome, skipping", chr))
    next
  }
  
  start <- max(1, artifact$start - CONTEXT_WINDOW)
  end <- min(chr_len, artifact$end + CONTEXT_WINDOW)
  
  message(sprintf("\nVisualizing artifact %d/%d: %s:%d-%d", 
                 i, nrow(artifacts), chr, artifact$start, artifact$end))
  
  # Extract feature info (clean up for filename)
  features <- artifact$overlapping_features
  feature_short <- gsub("[^A-Za-z0-9_]", "_", substr(features, 1, 40))
  
  # Create filename
  filename <- sprintf("%s_%d-%d_%s.pdf", 
                     gsub("chr", "", chr), 
                     artifact$start, 
                     artifact$end,
                     feature_short)
  
  pdf(file.path(path.expand(OUT_DIR), filename), width = 12, height = 8)
  
  tryCatch({
    # Genome axis
    gtrack <- GenomeAxisTrack()
    
    # GFF features in region
    gff_region <- subsetByOverlaps(gff, GRanges(chr, IRanges(start, end)))
    
    track_list <- list(gtrack)
    
    # Artifact region highlight (add first so it's in background)
    artifact_gr <- GRanges(chr, IRanges(artifact$start, artifact$end))
    artifact_track <- AnnotationTrack(artifact_gr,
                                     name = "Artifact",
                                     fill = "pink",
                                     col = "red",
                                     alpha = 0.5)
    track_list <- c(track_list, list(artifact_track))
    
    # BigWig tracks
    orc_file <- file.path(BW_DIR, paste0(ORC_SAMPLE, "_CPM.bw"))
    mcm_file <- file.path(BW_DIR, paste0(MCM_SAMPLE, "_CPM.bw"))
    input_file <- file.path(BW_DIR, paste0(INPUT_SAMPLE, "_CPM.bw"))
    
    if (file.exists(orc_file)) {
      orc_track <- DataTrack(range = orc_file,
                            chromosome = chr,
                            name = "ORC",
                            type = "histogram",
                            fill = "darkgreen",
                            col = "darkgreen",
                            ylim = c(0, NA))
      track_list <- c(track_list, list(orc_track))
    }
    
    if (file.exists(mcm_file)) {
      mcm_track <- DataTrack(range = mcm_file,
                            chromosome = chr,
                            name = "MCM",
                            type = "histogram",
                            fill = "purple",
                            col = "purple",
                            ylim = c(0, NA))
      track_list <- c(track_list, list(mcm_track))
    }
    
    if (file.exists(input_file)) {
      input_track <- DataTrack(range = input_file,
                              chromosome = chr,
                              name = "Input",
                              type = "histogram",
                              fill = "gray60",
                              col = "gray60",
                              ylim = c(0, NA))
      track_list <- c(track_list, list(input_track))
    }
    
    # Add GFF features if present
    if (length(gff_region) > 0) {
      # Color features by type
      feature_colors <- rep("lightblue", length(gff_region))
      feature_colors[grepl("Ty|transpos|LTR", gff_region$type, ignore.case = TRUE)] <- "orange"
      feature_colors[grepl("ARS|origin", gff_region$type, ignore.case = TRUE)] <- "red"
      feature_colors[grepl("tRNA", gff_region$type, ignore.case = TRUE)] <- "green"
      feature_colors[grepl("repeat", gff_region$type, ignore.case = TRUE)] <- "yellow"
      
      atrack <- AnnotationTrack(gff_region, 
                               name = "Features",
                               fill = feature_colors,
                               col = "black",
                               stacking = "squish",
                               fontsize = 6,
                               rotation.title = 0)
      track_list <- c(track_list, list(atrack))
    }
    
    # Plot all tracks
    plotTracks(track_list,
              from = start,
              to = end,
              chromosome = chr,
              main = sprintf("Artifact %d: %s:%d-%d\nFeatures: %s", 
                           i, chr, artifact$start, artifact$end,
                           features),
              cex.main = 0.7,
              cex.axis = 0.7,
              background.title = "darkgray")
    
    message("  Success!")
    
  }, error = function(e) {
    message(paste("  Error plotting:", e$message))
  })
  
  dev.off()
  message(paste("  Saved:", filename))
}

# ==============================================================================
# CREATE SUMMARY FIGURE
# ==============================================================================

message("\n=== Creating Summary Figure ===")

# Pick top 3 artifacts for summary
top_artifacts <- head(artifacts, 3)

pdf(file.path(path.expand(OUT_DIR), "artifact_summary.pdf"), 
    width = 12, height = 10)

par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))

for (i in 1:nrow(top_artifacts)) {
  artifact <- top_artifacts[i, ]
  chr <- artifact$chr
  
  message(sprintf("Adding to summary: %s:%d-%d", chr, artifact$start, artifact$end))
  
  tryCatch({
    # Simple plotting without Gviz for summary
    plot(1, type = "n", 
         xlim = c(artifact$start, artifact$end),
         ylim = c(0, 20),
         xlab = sprintf("%s position (bp)", chr),
         ylab = "Signal (CPM)",
         main = sprintf("Artifact %d: %s\n%s", 
                       i, artifact$window_id, artifact$overlapping_features),
         cex.main = 0.8)
    
    # Add shaded artifact region
    rect(artifact$start, 0, artifact$end, 20, 
         col = rgb(1, 0, 0, 0.2), border = NA)
    
    legend("topright", 
           legend = c("ORC", "MCM", "Input", "Artifact region"),
           col = c("darkgreen", "purple", "gray60", "pink"),
           lty = 1, lwd = 2, cex = 0.7)
    
  }, error = function(e) {
    message(paste("  Error in summary:", e$message))
  })
}

dev.off()

message("\n=== Visualization Complete ===")
message(paste("Output directory:", path.expand(OUT_DIR)))
message(paste("Generated", nrow(artifacts), "individual PDFs"))
message("Generated artifact_summary.pdf")
