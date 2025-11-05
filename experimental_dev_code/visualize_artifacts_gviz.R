#!/usr/bin/env Rscript
# Visualize Consensus Artifacts with Gviz (adjusted)

suppressPackageStartupMessages({
  library(Gviz)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

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
# HELPERS
# ==============================================================================

safe_dev_off <- function() {
  if (length(grDevices::dev.list()) > 0) try(grDevices::dev.off(), silent = TRUE)
}

# Harmonize seqname style to UCSC ("chrI".."chrXVI") and fix artifact chr labels if needed
harmonize_ucsc <- function(x) {
  if (inherits(x, "GRanges")) {
    GenomeInfoDb::seqlevelsStyle(x) <- "UCSC"
  }
  x
}

ensure_chr_prefix <- function(v) {
  v <- as.character(v)
  ifelse(grepl("^chr", v), v, paste0("chr", v))
}

# Get chromosome lengths from a BigWig header; fallback to GenomeInfoDb if missing
get_chr_len <- function(chr, bw_path, genome = "sacCer3") {
  if (file.exists(bw_path)) {
    bwf <- rtracklayer::BigWigFile(bw_path)
    si <- seqinfo(bwf)
    if (!is.null(si) && chr %in% names(seqlengths(si))) {
      L <- seqlengths(si)[chr]
      if (is.finite(L) && !is.na(L) && L > 0) return(as.integer(L))
    }
  }
  # Fallback: try UCSC chrom info (internet required)
  ci <- try(GenomeInfoDb::getChromInfoFromUCSC(genome), silent = TRUE)
  if (!inherits(ci, "try-error") && chr %in% ci$chrom) {
    return(as.integer(ci$size[match(chr, ci$chrom)]))
  }
  NA_integer_
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading artifacts...")
artifacts <- read.table(path.expand(ARTIFACT_FILE), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Normalize artifact chromosome labels to UCSC style
artifacts$chr <- ensure_chr_prefix(artifacts$chr)

message(paste("Found", nrow(artifacts), "artifacts to visualize"))

# Load GFF once and harmonize seqlevels style
message("Loading GFF annotations...")
GENOME_FEATURES <- import(GFF_FILE, format = "GFF")
GENOME_FEATURES <- harmonize_ucsc(GENOME_FEATURES)
message(paste("Loaded", length(GENOME_FEATURES), "genomic features"))

# ==============================================================================
# CREATE VISUALIZATIONS (ZOOMED AROUND EACH ARTIFACT)
# ==============================================================================

for (i in seq_len(nrow(artifacts))) {
  artifact <- artifacts[i, ]

  current_chromosome <- artifact$chr
  artifact_start <- artifact$start
  artifact_end <- artifact$end

  # Define region to visualize
  region_start <- max(1L, as.integer(artifact_start - CONTEXT_WINDOW))
  region_end   <- as.integer(artifact_end + CONTEXT_WINDOW)

  message(sprintf("\nVisualizing artifact %d/%d: %s:%d-%d",
                  i, nrow(artifacts), current_chromosome, artifact_start, artifact_end))

  # Build GRanges for 'which' queries and subsetting; harmonize style
  current_genome_range <- GRanges(seqnames = current_chromosome,
                                  ranges = IRanges(start = region_start, end = region_end))
  current_genome_range <- harmonize_ucsc(current_genome_range)

  # Initialize track container
  track_container <- list()

  # 1) Genome axis
  track_container[[length(track_container) + 1]] <- GenomeAxisTrack(
    name = sprintf("%s Axis", current_chromosome),
    fontcolor.title = "black",
    cex.title = 0.7,
    background.title = "white"
  )

  # 2) Artifact interval highlight
  artifact_gr <- GRanges(seqnames = current_chromosome,
                         ranges = IRanges(start = artifact_start, end = artifact_end))
  artifact_gr <- harmonize_ucsc(artifact_gr)

  track_container[[length(track_container) + 1]] <- AnnotationTrack(
    artifact_gr,
    name = "Artifact",
    fill = "pink",
    col = "red",
    alpha = 0.5,
    background.title = "lightcoral",
    fontcolor.title = "black",
    cex.title = 0.7
  )

  # 3) Add BigWig tracks (slice to window with rtracklayer::import(which=...))
  add_bw_slice <- function(bw_path, nm, col) {
    if (!file.exists(bw_path)) return(NULL)
    msg <- sprintf("  Loading %s bigwig...", nm)
    message(msg)
    dat <- import(bw_path, format = "BigWig", which = current_genome_range)
    if (length(dat) == 0L) {
      message(sprintf("  WARNING: No %s data in this region", nm))
      return(NULL)
    }
    DataTrack(
      range = dat, name = nm,
      showAxis = TRUE, showTitle = TRUE,
      type = "histogram",
      size = 1.2,
      background.title = "white",
      fontcolor.title = "black",
      cex.title = 0.7,
      col = col, fill = col
    )
  }

  orc_file   <- file.path(BW_DIR, paste0(ORC_SAMPLE,   "_CPM.bw"))
  mcm_file   <- file.path(BW_DIR, paste0(MCM_SAMPLE,   "_CPM.bw"))
  input_file <- file.path(BW_DIR, paste0(INPUT_SAMPLE, "_CPM.bw"))

  if (!is.null(trk <- add_bw_slice(orc_file,   "ORC",   "darkgreen"))) track_container[[length(track_container) + 1]] <- trk
  if (!is.null(trk <- add_bw_slice(mcm_file,   "MCM",   "purple")))    track_container[[length(track_container) + 1]] <- trk
  if (!is.null(trk <- add_bw_slice(input_file, "Input", "gray60")))    track_container[[length(track_container) + 1]] <- trk

  # 4) Genomic features for the window
  message("  Adding genomic features...")
  current_genome_feature <- keepSeqlevels(GENOME_FEATURES, current_chromosome, pruning.mode = "coarse")
  feature_subset <- GenomicRanges::subsetByOverlaps(current_genome_feature, current_genome_range)

  if (length(feature_subset) > 0L) {
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

  # Plot zoomed view
  message("  Creating PDF...")
  tryCatch({
    pdf(file = plot_output_file, width = 10, height = 8, bg = "white",
        compress = TRUE, colormodel = "srgb", useDingbats = FALSE)

    plotTracks(
      trackList = track_container,
      chromosome = current_chromosome,
      from = region_start, to = region_end,
      margin = 15, innerMargin = 5,
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
    safe_dev_off()
    message(paste("  ERROR:", e$message))
  })
}

# ==============================================================================
# CREATE CHROMOSOME-WIDE OVERVIEW PLOTS (EXPLICIT from/to + window='auto')
# ==============================================================================

message("\n=== Creating Chromosome-Wide Overview Plots ===")

artifact_chromosomes <- unique(artifacts$chr)

for (current_chromosome in artifact_chromosomes) {
  message(sprintf("\nCreating overview for %s...", current_chromosome))

  # All artifacts on this chromosome
  chr_artifacts <- artifacts[artifacts$chr == current_chromosome, , drop = FALSE]

  artifact_regions <- GRanges(
    seqnames = current_chromosome,
    ranges = IRanges(start = chr_artifacts$start, end = chr_artifacts$end),
    name = paste0("Artifact_", seq_len(nrow(chr_artifacts)))
  )
  artifact_regions <- harmonize_ucsc(artifact_regions)

  # Determine chromosome length from BigWig header (fallback to UCSC if needed)
  orc_file   <- file.path(BW_DIR, paste0(ORC_SAMPLE,   "_CPM.bw"))
  mcm_file   <- file.path(BW_DIR, paste0(MCM_SAMPLE,   "_CPM.bw"))
  input_file <- file.path(BW_DIR, paste0(INPUT_SAMPLE, "_CPM.bw"))

  chr_len <- get_chr_len(current_chromosome, if (file.exists(orc_file)) orc_file else input_file, genome = "sacCer3")
  if (!is.finite(chr_len) || is.na(chr_len) || chr_len <= 0) {
    message(sprintf("  ERROR: Could not resolve length for %s; skipping overview", current_chromosome))
    next
  }

  # Track container
  track_container <- list()

  # 1) Genome axis
  track_container[[length(track_container) + 1]] <- GenomeAxisTrack(
    name = sprintf("%s Position", current_chromosome),
    fontcolor.title = "black",
    cex.title = 0.7,
    background.title = "white"
  )

  # 2) Artifact regions (highlighted)
  track_container[[length(track_container) + 1]] <- AnnotationTrack(
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

  # 3) BigWig tracks using file paths directly; use window='auto' and no NA in ylim
  add_bw_file_track <- function(bw_path, nm, col) {
    if (!file.exists(bw_path)) return(NULL)
    DataTrack(
      range = bw_path,
      genome = "sacCer3",
      chromosome = current_chromosome,
      name = nm,
      showAxis = TRUE,
      showTitle = TRUE,
      type = "histogram",
      size = 1.0,
      background.title = "white",
      fontcolor.title = "black",
      cex.title = 0.7,
      col = col,
      fill = col,
      window = "auto"  # adaptive binning for full-chromosome
    )
  }

  if (!is.null(trk <- add_bw_file_track(orc_file,   "ORC",   "darkgreen"))) track_container[[length(track_container) + 1]] <- trk
  if (!is.null(trk <- add_bw_file_track(mcm_file,   "MCM",   "purple")))    track_container[[length(track_container) + 1]] <- trk
  if (!is.null(trk <- add_bw_file_track(input_file, "Input", "gray60")))    track_container[[length(track_container) + 1]] <- trk

  # 4) Genomic features (use dense stacking)
  message("  Adding genomic features...")
  current_genome_feature <- keepSeqlevels(GENOME_FEATURES, current_chromosome, pruning.mode = "coarse")
  if (length(current_genome_feature) > 0L) {
    major_features <- current_genome_feature[
      grepl("gene|CDS|ARS|Ty|tRNA|rRNA", mcols(current_genome_feature)$type, ignore.case = TRUE)
    ]
    if (length(major_features) > 0L) {
      track_container[[length(track_container) + 1]] <- AnnotationTrack(
        major_features,
        name = "Features",
        size = 0.8,
        background.title = "lightgray",
        fontcolor.title = "black",
        showAxis = FALSE,
        cex.title = 0.6,
        fill = "#8b4513",
        col = "#8b4513",
        stacking = "dense"
      )
    }
  }

  # Filename and plot
  filename <- sprintf("%s_full_chromosome_overview.pdf", gsub("chr", "", current_chromosome))
  plot_output_file <- file.path(path.expand(OUT_DIR), filename)

  message("  Creating chromosome-wide PDF...")
  tryCatch({
    pdf(file = plot_output_file, width = 16, height = 8, bg = "white",
        compress = TRUE, colormodel = "srgb", useDingbats = FALSE)

    plotTracks(
      trackList = track_container,
      chromosome = current_chromosome,
      from = 1L, to = as.integer(chr_len),   # explicit full-length bounds
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
    safe_dev_off()
    message(paste("  ERROR:", e$message))
  })
}

message("\n=================================================================")
message("Visualization Complete")
message("=================================================================")
message(paste("Output directory:", path.expand(OUT_DIR)))
message(paste("Generated", nrow(artifacts), "zoomed PDF files"))
message(paste("Generated", length(unique(artifacts$chr)), "chromosome-wide overview PDF files"))
