#!/usr/bin/env Rscript
# ChIP-seq Artifact Detection and Diagnosis
# Uses SGD annotations and sample metadata for accurate classification
# Scans chr VII, X, XIV for regions with elevated signal in IP vs Input

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# File paths
BW_DIR <- "/home/luised94/data/250930Bel/coverage"
FEATURE_DIR <- "/home/luised94/data/feature_files"
BLACKLIST_BED <- file.path(FEATURE_DIR, "20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed")
SGD_FEATURES <- file.path(FEATURE_DIR, "240830_SGD_features.tab")
SGD_GFF <- file.path(FEATURE_DIR, "240830_saccharomyces_cerevisiae.gff")
METADATA_CSV <- Sys.glob("/home/luised94/data/250930Bel/documentation/*.csv")[1]

# Analysis parameters
SUSPECT_CHRS <- c("chrVII", "chrX", "chrXIV", "chrXII")  # chrXII for rDNA comparison
WINDOW_SIZE <- 10000  # 10kb sliding windows
ARTIFACT_THRESHOLD <- 100  # Minimum signal to consider (arbitrary units)
IP_VS_INPUT_RATIO <- 3  # IP signal must be 3x higher than Input

# Chromosome lengths (S288C)
CHR_LENGTHS <- c(
  chrI = 230218, chrII = 813184, chrIII = 316620, chrIV = 1531933,
  chrV = 576874, chrVI = 270161, chrVII = 1090940, chrVIII = 562643,
  chrIX = 439888, chrX = 745751, chrXI = 666816, chrXII = 1078177,
  chrXIII = 924431, chrXIV = 784333, chrXV = 1091291, chrXVI = 948066
)

# Output directory
OUT_DIR <- "artifact_diagnosis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ==============================================================================
# 1. DISCOVER BIGWIG FILES AND EXTRACT SAMPLE IDs
# ==============================================================================

message("\n=== Discovering BigWig Files ===")

# Get all RAW bigwig files and sort them
bw_files <- list.files(BW_DIR, pattern = "_RAW\\.bw$", full.names = FALSE)
bw_files <- sort(bw_files)  # Sort to ensure consistent order

# Extract sample IDs from filenames (e.g., D25-12496_RAW.bw -> D25-12496)
sample_ids <- gsub("_RAW\\.bw", "", bw_files)

message(paste("Found", length(bw_files), "RAW bigwig files (sorted):"))
for (i in seq_along(sample_ids)) {
  message(paste("  ", i, ":", sample_ids[i]))
}

# ==============================================================================
# 2. LOAD SAMPLE METADATA AND MATCH BY POSITION
# ==============================================================================

message("\n=== Loading Sample Metadata ===")

if (!file.exists(METADATA_CSV)) {
  stop("Sample metadata CSV not found at: ", METADATA_CSV)
}

metadata <- read.csv(METADATA_CSV, stringsAsFactors = FALSE)
message(paste("Loaded metadata for", nrow(metadata), "rows from:"))
message(paste("  ", basename(METADATA_CSV)))

# Display column names
message("\nMetadata columns:")
print(colnames(metadata))

# Check if number of samples matches
if (nrow(metadata) != length(sample_ids)) {
  message("\nWARNING: Mismatch between metadata rows and bigwig files!")
  message(paste("  Metadata rows:", nrow(metadata)))
  message(paste("  Bigwig files:", length(sample_ids)))
  message("\nProceeding with assumption that they match by position...")
  message("Please verify the order is correct!\n")
}

# Identify antibody column (required)
antibody_col <- colnames(metadata)[grep("antibody|target|protein|factor", 
                                         colnames(metadata), ignore.case = TRUE)][1]

if (is.na(antibody_col)) {
  message("\nERROR: Could not find antibody column.")
  message("Available columns: ", paste(colnames(metadata), collapse = ", "))
  stop("Cannot proceed without antibody information")
}

message(paste("\nUsing antibody column:", antibody_col))

# Create sample mapping by matching positions
# Row 1 of metadata -> First sample ID, etc.
n_samples <- min(nrow(metadata), length(sample_ids))

sample_map <- data.frame(
  sample_id = sample_ids[1:n_samples],
  antibody = trimws(metadata[[antibody_col]][1:n_samples]),
  stringsAsFactors = FALSE
)

# Classify as Input or IP based on antibody column
sample_map$class <- ifelse(
  grepl("input|Input|INPUT|control|Control|mock|Mock|IgG", 
        sample_map$antibody, ignore.case = TRUE),
  "Input",
  "IP"
)

# Add any other useful metadata columns
if ("cell_cycle_treatment" %in% colnames(metadata)) {
  sample_map$treatment <- metadata$cell_cycle_treatment[1:n_samples]
}

# Print sample classification with position matching
message("\n=== Sample Classification (Matched by Position) ===")
message("Row | Sample ID    | Antibody | Class")
message("----+-------------+----------+-------")
for (i in 1:nrow(sample_map)) {
  message(sprintf("%3d | %-12s | %-8s | %s", 
                  i, sample_map$sample_id[i], 
                  sample_map$antibody[i], 
                  sample_map$class[i]))
}

message("\nPlease verify this mapping is correct!")
message("If incorrect, stop the script (Ctrl+C) and check file order.\n")

# Pause to let user verify (optional - comment out if running non-interactively)
# Sys.sleep(5)

# Save sample mapping
write.table(sample_map, file.path(OUT_DIR, "sample_mapping.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 3. LOAD BIGWIG FILES
# ==============================================================================

message("\n=== Loading BigWig Files ===")

# Load all samples from sample_map
signal_data <- list()
for (sample_id in sample_map$sample_id) {
  bw_file <- file.path(BW_DIR, paste0(sample_id, "_RAW.bw"))
  if (!file.exists(bw_file)) {
    message(paste("WARNING: File not found:", basename(bw_file)))
    next
  }
  message(paste("Loading:", basename(bw_file)))
  signal_data[[sample_id]] <- import(bw_file, format = "BigWig")
}

message(paste("\nSuccessfully loaded", length(signal_data), "bigwig files"))

# ==============================================================================
# 4. LOAD SGD FEATURES
# ==============================================================================

message("\n=== Loading SGD Features ===")

# Try to load SGD features file
if (file.exists(SGD_FEATURES)) {
  message(paste("Loading:", basename(SGD_FEATURES)))
  
  # SGD features format: tab-delimited with header
  # Columns typically: SGDID, Feature type, Feature qualifier, Feature name, 
  #                    Standard gene name, Alias, Parent feature name, 
  #                    Secondary SGDID, Chromosome, Start coordinate, 
  #                    Stop coordinate, Strand, Genetic position, etc.
  
  sgd_features <- read.table(SGD_FEATURES, sep = "\t", header = TRUE, 
                             quote = "", comment.char = "", 
                             stringsAsFactors = FALSE, fill = TRUE)
  
  message(paste("Loaded", nrow(sgd_features), "features"))
  message("SGD columns:")
  print(colnames(sgd_features)[1:min(10, ncol(sgd_features))])
  
} else if (file.exists(SGD_GFF)) {
  message(paste("SGD features file not found, loading GFF:", basename(SGD_GFF)))
  
  # Load GFF
  sgd_gff <- import(SGD_GFF, format = "GFF")
  
  # Convert to data frame
  sgd_features <- as.data.frame(sgd_gff)
  message(paste("Loaded", nrow(sgd_features), "features from GFF"))
  
} else {
  message("WARNING: Neither SGD features nor GFF file found")
  message("Falling back to blacklist file only")
  sgd_features <- NULL
}

# Load blacklist
message(paste("\nLoading blacklist:", basename(BLACKLIST_BED)))
blacklist <- read.table(BLACKLIST_BED, sep = "\t", header = FALSE,
                        col.names = c("chr", "start", "end", "name", 
                                      "score", "strand", "type", "gene"),
                        fill = TRUE, comment.char = "")

message(paste("Loaded", nrow(blacklist), "blacklist features"))

# Identify key feature types from blacklist
ty_elements <- blacklist[grep("Ty|delta|LTR|transpos", 
                               blacklist$type, ignore.case = TRUE), ]
y_elements <- blacklist[grep("Y_prime|Y'|YRF", 
                              blacklist$name, ignore.case = TRUE), ]
repeats <- blacklist[grep("repeat", blacklist$type, ignore.case = TRUE), ]
ars_elements <- blacklist[grep("ARS", blacklist$type, ignore.case = TRUE), ]
rdna <- blacklist[grep("rDNA|RDN", blacklist$name, ignore.case = TRUE), ]

message("\nFeature counts from blacklist:")
message(paste("  - Ty elements:", nrow(ty_elements)))
message(paste("  - Y' elements:", nrow(y_elements)))
message(paste("  - Repeats:", nrow(repeats)))
message(paste("  - ARS elements:", nrow(ars_elements)))
message(paste("  - rDNA:", nrow(rdna)))

# Convert to GRanges
blacklist_gr <- GRanges(
  seqnames = blacklist$chr,
  ranges = IRanges(start = blacklist$start + 1, end = blacklist$end),  # BED is 0-based
  name = blacklist$name,
  type = blacklist$type
)

ty_gr <- GRanges(
  seqnames = ty_elements$chr,
  ranges = IRanges(start = ty_elements$start + 1, end = ty_elements$end),
  name = ty_elements$name,
  type = ty_elements$type
)

y_gr <- GRanges(
  seqnames = y_elements$chr,
  ranges = IRanges(start = y_elements$start + 1, end = y_elements$end),
  name = y_elements$name
)

ars_gr <- GRanges(
  seqnames = ars_elements$chr,
  ranges = IRanges(start = ars_elements$start + 1, end = ars_elements$end),
  name = ars_elements$name
)

# ==============================================================================
# 5. SCAN CHROMOSOMES IN SLIDING WINDOWS
# ==============================================================================

message("\n=== Scanning Chromosomes for Artifacts ===")
message(paste("Window size:", WINDOW_SIZE, "bp"))
message(paste("Artifact threshold:", ARTIFACT_THRESHOLD, "AU"))
message(paste("IP/Input ratio threshold:", IP_VS_INPUT_RATIO))

# Function to calculate mean signal in a window
calc_window_signal <- function(bw_gr, chr, start, end) {
  window_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  overlaps <- findOverlaps(window_gr, bw_gr)
  
  if (length(overlaps) == 0) return(0)
  
  hits <- bw_gr[subjectHits(overlaps)]
  mean(hits$score, na.rm = TRUE)
}

# Pre-allocate lists for speed (avoid repeated rbind)
all_windows_list <- list()
counter <- 1

# Pre-subset bigwigs by chromosome (do once, not in loop)
bw_by_chr <- list()
for (chr in SUSPECT_CHRS) {
  bw_by_chr[[chr]] <- list()
  for (sample_id in names(signal_data)) {
    bw_by_chr[[chr]][[sample_id]] <- signal_data[[sample_id]][seqnames(signal_data[[sample_id]]) == chr]
  }
}

# Scan each chromosome
for (chr in SUSPECT_CHRS) {
  message(paste("\nScanning", chr, "..."))
  chr_len <- CHR_LENGTHS[chr]
  
  # Create sliding windows
  window_starts <- seq(1, chr_len - WINDOW_SIZE, by = WINDOW_SIZE)
  window_ends <- window_starts + WINDOW_SIZE - 1
  window_ends[length(window_ends)] <- min(window_ends[length(window_ends)], chr_len)
  n_windows <- length(window_starts)
  
  message(paste("  Generated", n_windows, "windows"))
  
  # Calculate signal for each window across all samples
  for (i in seq_along(window_starts)) {
    window_start <- window_starts[i]
    window_end <- window_ends[i]
    window_id <- paste0(chr, ":", window_start, "-", window_end)
    
    # Calculate mean signal for each sample in this window
    for (sample_id in names(signal_data)) {
      bw_chr <- bw_by_chr[[chr]][[sample_id]]  # Pre-subsetted
      
      signal <- calc_window_signal(bw_chr, chr, window_start, window_end)
      
      # Get sample info from metadata (lookup once)
      sample_info <- sample_map[sample_map$sample_id == sample_id, ]
      
      # Store in list (much faster than rbind)
      all_windows_list[[counter]] <- data.frame(
        chr = chr,
        start = window_start,
        end = window_end,
        window_id = window_id,
        sample_id = sample_id,
        sample_class = sample_info$class,
        antibody = sample_info$antibody,
        signal = signal,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
    
    if (i %% 20 == 0 || i == n_windows) {
      message(sprintf("    Progress: %d/%d (%.1f%%)", i, n_windows, i/n_windows*100))
    }
  }
}

# Combine all at once (much faster than repeated rbind)
message("\nCombining results...")
all_windows <- do.call(rbind, all_windows_list)

message("\nWindow scanning complete")

# Save raw window data
write.table(all_windows, file.path(OUT_DIR, "all_window_signals.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 6. IDENTIFY ARTIFACT REGIONS
# ==============================================================================

message("\n=== Identifying Artifact Regions ===")

# Calculate mean signal per window across sample classes
window_summary <- all_windows %>%
  group_by(chr, start, end, window_id, sample_class) %>%
  summarise(
    mean_signal = mean(signal, na.rm = TRUE),
    median_signal = median(signal, na.rm = TRUE),
    max_signal = max(signal, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample_class,
    values_from = c(mean_signal, median_signal, max_signal, n_samples),
    names_sep = "_"
  )

# Handle case where there might be no Input samples
if (!"mean_signal_Input" %in% colnames(window_summary)) {
  window_summary$mean_signal_Input <- 0
  message("WARNING: No Input samples found, using 0 for Input signal")
}

# Calculate IP/Input ratio
window_summary$ip_input_ratio <- window_summary$mean_signal_IP / 
                                  (window_summary$mean_signal_Input + 0.1)

# Print diagnostic info
message("\n=== Signal Distribution Diagnostics ===")
message(sprintf("IP signal range: %.1f - %.1f (median: %.1f)", 
                min(window_summary$mean_signal_IP, na.rm = TRUE),
                max(window_summary$mean_signal_IP, na.rm = TRUE),
                median(window_summary$mean_signal_IP, na.rm = TRUE)))
message(sprintf("Input signal range: %.1f - %.1f (median: %.1f)", 
                min(window_summary$mean_signal_Input, na.rm = TRUE),
                max(window_summary$mean_signal_Input, na.rm = TRUE),
                median(window_summary$mean_signal_Input, na.rm = TRUE)))
message(sprintf("IP/Input ratio range: %.1f - %.1f (median: %.1f)", 
                min(window_summary$ip_input_ratio, na.rm = TRUE),
                max(window_summary$ip_input_ratio, na.rm = TRUE),
                median(window_summary$ip_input_ratio, na.rm = TRUE)))

# Count windows passing each threshold
n_high_ip <- sum(window_summary$mean_signal_IP > ARTIFACT_THRESHOLD, na.rm = TRUE)
n_high_ratio <- sum(window_summary$ip_input_ratio > IP_VS_INPUT_RATIO, na.rm = TRUE)
n_low_input <- sum(window_summary$mean_signal_Input < 50, na.rm = TRUE)

message("\nWindows passing thresholds:")
message(sprintf("  IP signal > %d: %d windows", ARTIFACT_THRESHOLD, n_high_ip))
message(sprintf("  IP/Input ratio > %d: %d windows", IP_VS_INPUT_RATIO, n_high_ratio))
message(sprintf("  Input signal < 50: %d windows", n_low_input))

# Identify artifact windows with relaxed thresholds if none found
artifact_windows <- window_summary %>%
  filter(
    mean_signal_IP > ARTIFACT_THRESHOLD,
    ip_input_ratio > IP_VS_INPUT_RATIO,
    mean_signal_Input < 50
  ) %>%
  arrange(chr, start)

message(sprintf("\nArtifacts (all 3 thresholds): %d windows", nrow(artifact_windows)))

# If no artifacts found, try relaxed criteria
if (nrow(artifact_windows) == 0) {
  message("\n=== No artifacts with strict thresholds. Trying relaxed criteria ===")
  
  # Try: High IP signal OR high ratio
  artifact_windows_relaxed <- window_summary %>%
    filter(
      mean_signal_IP > ARTIFACT_THRESHOLD | 
      (ip_input_ratio > IP_VS_INPUT_RATIO & mean_signal_IP > 50)
    ) %>%
    arrange(desc(mean_signal_IP)) %>%
    head(50)  # Top 50
  
  message(sprintf("Found %d candidate regions with relaxed criteria", 
                  nrow(artifact_windows_relaxed)))
  
  if (nrow(artifact_windows_relaxed) > 0) {
    # Use relaxed criteria
    artifact_windows <- artifact_windows_relaxed
    message("Using relaxed artifact windows for further analysis")
  }
}

message(paste("Identified", nrow(artifact_windows), "artifact windows"))

# Always save top signal regions for inspection
top_signal_windows <- window_summary %>%
  arrange(desc(mean_signal_IP)) %>%
  head(100)

write.table(top_signal_windows, file.path(OUT_DIR, "top_100_signal_windows.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(paste("Saved top 100 signal windows to:", 
              file.path(OUT_DIR, "top_100_signal_windows.txt")))

if (nrow(artifact_windows) > 0) {
  write.table(artifact_windows, file.path(OUT_DIR, "artifact_windows.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Create BED file
  artifact_bed <- artifact_windows %>%
    select(chr, start, end) %>%
    mutate(name = paste0("artifact_", row_number()),
           score = 0,
           strand = ".")
  
  write.table(artifact_bed, file.path(OUT_DIR, "artifact_regions.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message(paste("Saved artifact BED file:", 
                file.path(OUT_DIR, "artifact_regions.bed")))
}

# ==============================================================================
# 7. ANNOTATE ARTIFACTS WITH FEATURES
# ==============================================================================

message("\n=== Annotating Artifacts ===")

if (nrow(artifact_windows) > 0) {
  # Convert to GRanges
  artifact_gr <- GRanges(
    seqnames = artifact_windows$chr,
    ranges = IRanges(start = artifact_windows$start, end = artifact_windows$end),
    window_id = artifact_windows$window_id,
    ip_signal = artifact_windows$mean_signal_IP,
    input_signal = artifact_windows$mean_signal_Input,
    ratio = artifact_windows$ip_input_ratio
  )
  
  # Find overlaps with features
  ty_overlaps <- findOverlaps(artifact_gr, ty_gr)
  y_overlaps <- findOverlaps(artifact_gr, y_gr)
  ars_overlaps <- findOverlaps(artifact_gr, ars_gr)
  
  # Create annotation dataframe
  artifact_annotations <- data.frame(
    window_id = artifact_windows$window_id,
    chr = artifact_windows$chr,
    start = artifact_windows$start,
    end = artifact_windows$end,
    ip_signal = round(artifact_windows$mean_signal_IP, 2),
    input_signal = round(artifact_windows$mean_signal_Input, 2),
    ratio = round(artifact_windows$ip_input_ratio, 2),
    overlaps_ty = seq_along(artifact_gr) %in% queryHits(ty_overlaps),
    overlaps_y_prime = seq_along(artifact_gr) %in% queryHits(y_overlaps),
    overlaps_ars = seq_along(artifact_gr) %in% queryHits(ars_overlaps),
    ty_names = "",
    y_names = "",
    ars_names = ""
  )
  
  # Add feature names
  for (i in seq_along(artifact_gr)) {
    ty_hits <- subjectHits(ty_overlaps)[queryHits(ty_overlaps) == i]
    if (length(ty_hits) > 0) {
      artifact_annotations$ty_names[i] <- paste(ty_gr$name[ty_hits], collapse = "; ")
    }
    
    y_hits <- subjectHits(y_overlaps)[queryHits(y_overlaps) == i]
    if (length(y_hits) > 0) {
      artifact_annotations$y_names[i] <- paste(y_gr$name[y_hits], collapse = "; ")
    }
    
    ars_hits <- subjectHits(ars_overlaps)[queryHits(ars_overlaps) == i]
    if (length(ars_hits) > 0) {
      artifact_annotations$ars_names[i] <- paste(ars_gr$name[ars_hits], collapse = "; ")
    }
  }
  
  write.table(artifact_annotations, 
              file.path(OUT_DIR, "artifact_annotations.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Annotation complete")
  
  # Print feature overlap summary
  message("\n=== Feature Overlap Summary ===")
  message(paste("Artifacts overlapping Ty elements:", 
                sum(artifact_annotations$overlaps_ty)))
  message(paste("Artifacts overlapping Y' elements:", 
                sum(artifact_annotations$overlaps_y_prime)))
  message(paste("Artifacts overlapping ARS:", 
                sum(artifact_annotations$overlaps_ars)))
  message(paste("Artifacts with no annotated features:", 
                sum(!artifact_annotations$overlaps_ty & 
                    !artifact_annotations$overlaps_y_prime & 
                    !artifact_annotations$overlaps_ars)))
}

# ==============================================================================
# 8. VISUALIZATIONS
# ==============================================================================

message("\n=== Generating Plots ===")

# Plot 1: Signal tracks across chromosomes
p1 <- ggplot(all_windows, aes(x = start / 1e6, y = signal, 
                               color = sample_class, group = sample_id)) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  facet_wrap(~ chr, scales = "free_x", ncol = 1) +
  geom_hline(yintercept = ARTIFACT_THRESHOLD, linetype = "dashed", 
             color = "red", alpha = 0.7) +
  labs(title = "ChIP-seq Signal Across Suspect Chromosomes",
       subtitle = paste0("Red line = artifact threshold (", ARTIFACT_THRESHOLD, " AU); ",
                        "chrXII included for rDNA comparison"),
       x = "Position (Mb)", y = "Signal (AU)",
       color = "Sample Class") +
  theme_bw() +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "gray90"))

ggsave(file.path(OUT_DIR, "chromosome_signal_tracks.pdf"), p1, 
       width = 14, height = 12)

# Plot 2: IP vs Input scatter
if (nrow(artifact_windows) > 0 && "mean_signal_Input" %in% colnames(window_summary)) {
  p2 <- ggplot(window_summary, aes(x = mean_signal_Input + 0.1, 
                                    y = mean_signal_IP + 0.1)) +
    geom_point(alpha = 0.2, color = "gray60", size = 1) +
    geom_point(data = artifact_windows, 
               aes(x = mean_signal_Input + 0.1, y = mean_signal_IP + 0.1),
               color = "red", size = 2, alpha = 0.7) +
    geom_abline(slope = IP_VS_INPUT_RATIO, intercept = 0, 
                linetype = "dashed", color = "blue", linewidth = 1) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    facet_wrap(~ chr) +
    labs(title = "IP vs Input Signal Distribution",
         subtitle = paste0("Red points = artifacts; Blue line = ", 
                          IP_VS_INPUT_RATIO, "x ratio threshold"),
         x = "Input Signal (AU, log scale)", 
         y = "IP Signal (AU, log scale)") +
    theme_bw()
  
  ggsave(file.path(OUT_DIR, "ip_vs_input_scatter.pdf"), p2, 
         width = 12, height = 9)
}

# Plot 3: Artifact characterization
if (nrow(artifact_windows) > 0) {
  p3 <- ggplot(artifact_annotations, 
               aes(x = chr, y = ip_signal)) +
    geom_boxplot(outlier.shape = NA, fill = "gray80") +
    geom_jitter(aes(color = overlaps_ty), width = 0.2, alpha = 0.7, size = 2) +
    scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red"),
                       labels = c("No", "Yes")) +
    labs(title = "Artifact Signal by Chromosome",
         subtitle = "Color indicates overlap with Ty transposable elements",
         x = "Chromosome", y = "Mean IP Signal (AU)",
         color = "Overlaps\nTy Element") +
    theme_bw() +
    theme(legend.position = "right")
  
  ggsave(file.path(OUT_DIR, "artifact_characterization.pdf"), p3, 
         width = 10, height = 6)
}

# ==============================================================================
# 9. SUMMARY REPORT
# ==============================================================================

message("\n=== Generating Summary Report ===")

sink(file.path(OUT_DIR, "ANALYSIS_SUMMARY.txt"))

cat("=================================================================\n")
cat("ChIP-seq Artifact Detection Analysis\n")
cat("=================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Input Files:\n")
cat(paste("  BigWig directory:", BW_DIR, "\n"))
cat(paste("  Sample metadata:", basename(METADATA_CSV), "\n"))
cat(paste("  Blacklist:", basename(BLACKLIST_BED), "\n\n"))

cat("Analysis Parameters:\n")
cat(paste("  Chromosomes:", paste(SUSPECT_CHRS, collapse = ", "), "\n"))
cat(paste("  Window size:", WINDOW_SIZE, "bp\n"))
cat(paste("  Artifact threshold:", ARTIFACT_THRESHOLD, "AU\n"))
cat(paste("  IP/Input ratio:", IP_VS_INPUT_RATIO, "\n\n"))

cat("Samples Analyzed:\n")
cat(paste("  Total:", nrow(sample_map), "\n"))
sample_counts <- table(sample_map$class)
for (class in names(sample_counts)) {
  cat(paste("  ", class, ":", sample_counts[class], "\n"))
}
cat("\nSample Details:\n")
print(sample_map, row.names = FALSE)

cat("\n\nResults:\n")
cat(paste("  Total windows scanned:", nrow(window_summary), "\n"))
cat(paste("  Artifact windows detected:", nrow(artifact_windows), "\n"))

if (nrow(artifact_windows) > 0) {
  cat("\n\nArtifacts by Chromosome:\n")
  chr_summary <- artifact_annotations %>%
    group_by(chr) %>%
    summarise(
      n_artifacts = n(),
      mean_ip = round(mean(ip_signal), 1),
      mean_input = round(mean(input_signal), 1),
      mean_ratio = round(mean(ratio), 1),
      with_ty = sum(overlaps_ty),
      with_y_prime = sum(overlaps_y_prime),
      with_ars = sum(overlaps_ars)
    )
  print(chr_summary, row.names = FALSE)
  
  cat("\n\nTop 15 Artifact Regions (by IP signal):\n")
  top_artifacts <- artifact_annotations %>%
    arrange(desc(ip_signal)) %>%
    head(15) %>%
    select(chr, start, end, ip_signal, ratio, overlaps_ty, 
           overlaps_y_prime, overlaps_ars, ty_names)
  print(top_artifacts, row.names = FALSE)
}

cat("\n=================================================================\n")
cat("Output Files:\n")
cat("  - all_window_signals.txt: Raw signal for all windows\n")
cat("  - artifact_windows.txt: Detected artifact windows\n")
cat("  - artifact_regions.bed: BED format (add to blacklist)\n")
cat("  - artifact_annotations.txt: Features overlapping artifacts\n")
cat("  - sample_mapping.txt: Sample classification used\n")
cat("  - *.pdf: Visualization plots\n")
cat("=================================================================\n")

sink()

message("\n=== ANALYSIS COMPLETE ===")
message(paste("Output directory:", OUT_DIR))
message("\nKey files to review:")
message(paste("  1.", file.path(OUT_DIR, "ANALYSIS_SUMMARY.txt")))
message(paste("  2.", file.path(OUT_DIR, "artifact_annotations.txt")))
message(paste("  3.", file.path(OUT_DIR, "artifact_regions.bed")))
