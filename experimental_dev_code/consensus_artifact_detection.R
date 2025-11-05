#!/usr/bin/env Rscript
# Consensus Artifact Detection - Per-Sample Approach
# Finds regions with elevated signal in MOST IP samples but NOT in Input

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Use existing window data from previous run
ARTIFACT_DIR <- "~/artifact_diagnosis"  # Moved outside git repo
WINDOW_DATA <- file.path(ARTIFACT_DIR, "all_window_signals.txt")

# Parameters
SUSPECT_CHRS <- c("chrVII", "chrX", "chrXIV")  # Exclude chrXII (rDNA is expected)
PERCENTILE_THRESHOLD <- 0.95  # Top 5% within each sample
MIN_IP_FRACTION <- 0.5  # Present in at least 50% of IP samples
MAX_INPUT_FRACTION <- 0.2  # Present in at most 20% of Input samples

# Feature annotations
FEATURE_DIR <- "/home/luised94/data/feature_files"
SGD_GFF <- file.path(FEATURE_DIR, "240830_saccharomyces_cerevisiae.gff")

# Output
OUT_DIR <- "consensus_artifacts"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rtracklayer)
  library(GenomicRanges)
})

message("=================================================================")
message("Consensus Artifact Detection - Per-Sample Analysis")
message("=================================================================\n")

# ==============================================================================
# 1. LOAD WINDOW DATA
# ==============================================================================

message("=== Loading Window Data ===")

if (!file.exists(WINDOW_DATA)) {
  stop("Window data not found. Run chipseq_artifact_diagnosis.R first!")
}

all_windows <- read.table(WINDOW_DATA, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)

message(paste("Loaded", nrow(all_windows), "window measurements"))

# Filter to suspect chromosomes only
all_windows <- all_windows %>%
  filter(chr %in% SUSPECT_CHRS)

message(paste("Analyzing", nrow(all_windows),
              "measurements on chromosomes:", paste(SUSPECT_CHRS, collapse = ", ")))

# Get sample info
ip_samples <- unique(all_windows$sample_id[all_windows$sample_class == "IP"])
input_samples <- unique(all_windows$sample_id[all_windows$sample_class == "Input"])

message(paste("\nSamples:"))
message(paste("  IP samples:", length(ip_samples)))
message(paste("  Input samples:", length(input_samples)))

# ==============================================================================
# 2. CALCULATE PER-SAMPLE PERCENTILES
# ==============================================================================

message("\n=== Calculating Per-Sample Signal Thresholds ===")

# For each sample, calculate the threshold for "elevated" signal
sample_thresholds <- all_windows %>%
  group_by(sample_id, sample_class) %>%
  summarise(
    threshold_95 = quantile(signal, PERCENTILE_THRESHOLD, na.rm = TRUE),
    median_signal = median(signal, na.rm = TRUE),
    mean_signal = mean(signal, na.rm = TRUE),
    .groups = "drop"
  )

message("\nPer-sample thresholds (95th percentile):")
print(sample_thresholds, n = Inf)

# Save thresholds
write.table(sample_thresholds, file.path(OUT_DIR, "sample_thresholds.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 3. IDENTIFY ELEVATED WINDOWS PER SAMPLE
# ==============================================================================

message("\n=== Identifying Elevated Windows Per Sample ===")

# Join thresholds to window data
all_windows <- all_windows %>%
  left_join(sample_thresholds %>% select(sample_id, threshold_95),
            by = "sample_id")

# Mark windows as elevated if signal > threshold
all_windows <- all_windows %>%
  mutate(is_elevated = signal > threshold_95)

# Count elevated windows per sample
elevated_counts <- all_windows %>%
  group_by(sample_id, sample_class) %>%
  summarise(
    n_windows = n(),
    n_elevated = sum(is_elevated),
    pct_elevated = n_elevated / n_windows * 100,
    .groups = "drop"
  )

message("\nElevated windows per sample:")
print(elevated_counts, n = Inf)

# ==============================================================================
# 4. FIND CONSENSUS ARTIFACTS
# ==============================================================================

message("\n=== Finding Consensus Artifact Regions ===")

# For each window, count how many samples show elevation
consensus <- all_windows %>%
  group_by(chr, start, end, window_id, sample_class) %>%
  summarise(
    n_samples = n(),
    n_elevated = sum(is_elevated),
    fraction_elevated = n_elevated / n_samples,
    mean_signal = mean(signal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample_class,
    values_from = c(n_elevated, fraction_elevated, mean_signal, n_samples),
    names_sep = "_"
  )

# Identify consensus artifacts:
# - Elevated in >=% of IP samples
# - Elevated in >=% of Input samples
consensus_artifacts <- consensus %>%
  filter(
    fraction_elevated_IP >= MIN_IP_FRACTION,
    fraction_elevated_Input <= MAX_INPUT_FRACTION
  ) %>%
  arrange(desc(fraction_elevated_IP))

message(sprintf("\nFound %d consensus artifact windows:", nrow(consensus_artifacts)))
message(sprintf("  Elevated in >=%.0f%% of IP samples", MIN_IP_FRACTION * 100))
message(sprintf("  Elevated in >=%.0f%% of Input samples", MAX_INPUT_FRACTION * 100))

if (nrow(consensus_artifacts) == 0) {
  message("\nNo consensus artifacts found with current thresholds.")
  message("Try adjusting MIN_IP_FRACTION or MAX_INPUT_FRACTION")
  quit(save = "no")
}

# Save consensus artifacts
write.table(consensus_artifacts, file.path(OUT_DIR, "consensus_artifacts.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Create BED file
artifact_bed <- consensus_artifacts %>%
  select(chr, start, end) %>%
  mutate(
    name = paste0("consensus_artifact_", row_number()),
    score = 0,
    strand = "."
  )

write.table(artifact_bed, file.path(OUT_DIR, "consensus_artifacts.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(paste("Saved:", file.path(OUT_DIR, "consensus_artifacts.bed")))

# ==============================================================================
# 5. ANNOTATE WITH GENOMIC FEATURES FROM GFF
# ==============================================================================

message("\n=== Annotating Artifacts with SGD GFF Features ===")

# Load GFF file
if (!file.exists(SGD_GFF)) {
  message(paste("WARNING: GFF file not found at:", SGD_GFF))
  message("Skipping feature annotation")

  # Add empty annotation columns
  consensus_artifacts$overlapping_features <- ""
  consensus_artifacts$feature_types <- ""
  consensus_artifacts$feature_names <- ""
  consensus_artifacts$has_ty <- FALSE
  consensus_artifacts$has_ars <- FALSE
  consensus_artifacts$has_repeat <- FALSE
  consensus_artifacts$has_y_prime <- FALSE
  consensus_artifacts$has_telomere <- FALSE

} else {
  message(paste("Loading GFF:", basename(SGD_GFF)))

  # Import GFF
  gff <- import(SGD_GFF, format = "GFF")

  message(paste("Loaded", length(gff), "features from GFF"))

  # Convert to data frame for easier manipulation
  gff_df <- as.data.frame(gff)

  # Create GRanges for artifacts
  artifact_gr <- GRanges(
    seqnames = consensus_artifacts$chr,
    ranges = IRanges(start = consensus_artifacts$start, end = consensus_artifacts$end),
    window_id = consensus_artifacts$window_id
  )

  # Find overlaps
  overlaps <- findOverlaps(artifact_gr, gff)

  # Initialize annotation columns
  consensus_artifacts$overlapping_features <- ""
  consensus_artifacts$feature_types <- ""
  consensus_artifacts$feature_names <- ""

  # Annotate each artifact with overlapping features
  for (i in seq_along(artifact_gr)) {
    hits <- subjectHits(overlaps)[queryHits(overlaps) == i]

    if (length(hits) > 0) {
      hit_features <- gff_df[hits, ]

      # Extract feature information
      feature_types <- unique(hit_features$type)

      # Try to get feature names from different possible columns
      feature_names <- character()
      if ("Name" %in% colnames(hit_features)) {
        feature_names <- c(feature_names, na.omit(hit_features$Name))
      }
      if ("ID" %in% colnames(hit_features)) {
        feature_names <- c(feature_names, na.omit(hit_features$ID))
      }
      if ("gene" %in% colnames(hit_features)) {
        feature_names <- c(feature_names, na.omit(hit_features$gene))
      }

      feature_names <- unique(feature_names)

      # Combine all information
      all_info <- unique(c(feature_types, feature_names))

      consensus_artifacts$overlapping_features[i] <- paste(all_info, collapse = "; ")
      consensus_artifacts$feature_types[i] <- paste(feature_types, collapse = "; ")
      consensus_artifacts$feature_names[i] <- paste(feature_names, collapse = "; ")
    }
  }

  # Identify specific feature categories
  consensus_artifacts <- consensus_artifacts %>%
    mutate(
      has_ty = grepl("Ty|delta|transpos|LTR|long_terminal_repeat",
                     overlapping_features, ignore.case = TRUE),
      has_ars = grepl("ARS|autonomously_replicating_sequence|origin",
                      overlapping_features, ignore.case = TRUE),
      has_repeat = grepl("repeat|repetitive",
                         overlapping_features, ignore.case = TRUE),
      has_y_prime = grepl("Y_prime|Y'|YRF",
                          overlapping_features, ignore.case = TRUE),
      has_telomere = grepl("telomere|TEL",
                           overlapping_features, ignore.case = TRUE),
      has_trna = grepl("tRNA|transfer_RNA",
                       overlapping_features, ignore.case = TRUE),
      has_rrna = grepl("rRNA|ribosomal_RNA",
                       overlapping_features, ignore.case = TRUE)
    )

  message("Feature annotation complete")
}

# Save annotated artifacts
write.table(consensus_artifacts, file.path(OUT_DIR, "consensus_artifacts_annotated.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 6. SUMMARIZE FEATURE OVERLAPS
# ==============================================================================

message("\n=== Feature Overlap Summary ===")

feature_summary <- data.frame(
  Feature = c("Ty elements", "ARS elements", "Repeats", "Y' elements",
              "Telomeres", "tRNA", "rRNA", "No annotation"),
  Count = c(
    sum(consensus_artifacts$has_ty),
    sum(consensus_artifacts$has_ars),
    sum(consensus_artifacts$has_repeat),
    sum(consensus_artifacts$has_y_prime),
    sum(consensus_artifacts$has_telomere),
    sum(consensus_artifacts$has_trna),
    sum(consensus_artifacts$has_rrna),
    sum(consensus_artifacts$overlapping_features == "")
  )
)

print(feature_summary)

# ==============================================================================
# 7. DETAILED PER-SAMPLE VIEW
# ==============================================================================

message("\n=== Generating Per-Sample Heatmap Data ===")

# For each artifact window, get signal from each sample
artifact_windows <- consensus_artifacts$window_id

artifact_details <- all_windows %>%
  filter(window_id %in% artifact_windows) %>%
  select(window_id, chr, start, end, sample_id, sample_class,
         antibody, signal, is_elevated) %>%
  arrange(window_id, sample_class, sample_id)

write.table(artifact_details, file.path(OUT_DIR, "artifact_per_sample_signals.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 8. VISUALIZATIONS
# ==============================================================================

message("\n=== Generating Visualizations ===")

# Plot 1: Artifact distribution by chromosome
if (nrow(consensus_artifacts) > 0) {
  p1 <- ggplot(consensus_artifacts, aes(x = chr, fill = has_ty)) +
    geom_bar() +
    labs(title = "Consensus Artifacts by Chromosome",
         subtitle = sprintf("Total: %d regions elevated in >=%.0f%% of IPs, >=%.0f%% of Inputs",
                           nrow(consensus_artifacts),
                           MIN_IP_FRACTION * 100,
                           MAX_INPUT_FRACTION * 100),
         x = "Chromosome", y = "Number of Artifact Windows",
         fill = "Contains\nTy Element") +
    theme_bw() +
    theme(legend.position = "right")

  ggsave(file.path(OUT_DIR, "artifacts_by_chromosome.pdf"), p1,
         width = 8, height = 6)
}

# Plot 2: IP fraction vs Input fraction
p2 <- ggplot(consensus, aes(x = fraction_elevated_Input,
                             y = fraction_elevated_IP)) +
  geom_point(alpha = 0.3, color = "gray60") +
  geom_point(data = consensus_artifacts,
             aes(x = fraction_elevated_Input, y = fraction_elevated_IP),
             color = "red", size = 2) +
  geom_hline(yintercept = MIN_IP_FRACTION, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = MAX_INPUT_FRACTION, linetype = "dashed", color = "blue") +
  labs(title = "Consensus Artifact Detection",
       subtitle = "Red points = artifacts (high in IPs, low in Inputs)",
       x = "Fraction of Input Samples with Elevated Signal",
       y = "Fraction of IP Samples with Elevated Signal") +
  theme_bw()

ggsave(file.path(OUT_DIR, "consensus_scatter.pdf"), p2,
       width = 8, height = 6)

# Plot 3: Signal heatmap for top artifacts
if (nrow(consensus_artifacts) > 0) {
  top_artifacts <- head(consensus_artifacts$window_id, 20)

  heatmap_data <- artifact_details %>%
    filter(window_id %in% top_artifacts) %>%
    mutate(
      window_label = paste0(chr, ":", start/1000, "kb"),
      sample_label = paste0(sample_id, " (", antibody, ")")
    )

  p3 <- ggplot(heatmap_data, aes(x = sample_label, y = window_label, fill = signal)) +
    geom_tile() +
    geom_point(data = heatmap_data %>% filter(is_elevated),
               aes(x = sample_label, y = window_label),
               shape = 4, size = 1, color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    facet_grid(~ sample_class, scales = "free_x", space = "free_x") +
    labs(title = "Top 20 Consensus Artifacts - Per-Sample Signal",
         subtitle = "X marks = elevated within that sample",
         x = "", y = "Genomic Region",
         fill = "Signal (AU)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8))

  ggsave(file.path(OUT_DIR, "artifact_heatmap.pdf"), p3,
         width = 14, height = 10)
}

# ==============================================================================
# 9. GENERATE REPORT
# ==============================================================================

message("\n=== Generating Report ===")

sink(file.path(OUT_DIR, "CONSENSUS_ANALYSIS_REPORT.txt"))

cat("=================================================================\n")
cat("Consensus Artifact Detection Report\n")
cat("=================================================================\n\n")

cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Analysis Parameters:\n")
cat(paste("  Chromosomes analyzed:", paste(SUSPECT_CHRS, collapse = ", "), "\n"))
cat(paste("  Per-sample threshold: Top", (1 - PERCENTILE_THRESHOLD) * 100, "% of signal\n"))
cat(paste("  Minimum IP fraction:", MIN_IP_FRACTION,
          sprintf("(>=%.0f%% of IP samples)\n", MIN_IP_FRACTION * 100)))
cat(paste("  Maximum Input fraction:", MAX_INPUT_FRACTION,
          sprintf("(>=%.0f%% of Input samples)\n", MAX_INPUT_FRACTION * 100)))

cat("\nSamples:\n")
cat(paste("  IP samples:", length(ip_samples), "\n"))
cat(paste("  Input samples:", length(input_samples), "\n"))

cat("\nResults:\n")
cat(paste("  Total windows analyzed:", length(unique(all_windows$window_id)), "\n"))
cat(paste("  Consensus artifacts found:", nrow(consensus_artifacts), "\n\n"))

cat("Feature Overlap Summary:\n")
print(feature_summary, row.names = FALSE)

cat("\n\nTop 20 Consensus Artifacts:\n")
top20 <- consensus_artifacts %>%
  select(chr, start, end, fraction_elevated_IP, fraction_elevated_Input,
         mean_signal_IP, mean_signal_Input, overlapping_features) %>%
  head(20)
print(top20, row.names = FALSE)

cat("\n\nArtifacts by Chromosome:\n")
chr_summary <- consensus_artifacts %>%
  group_by(chr) %>%
  summarise(
    n_artifacts = n(),
    mean_ip_fraction = mean(fraction_elevated_IP),
    with_ty = sum(has_ty),
    with_ars = sum(has_ars),
    with_repeats = sum(has_repeat)
  )
print(chr_summary, row.names = FALSE)

cat("\n=================================================================\n")
cat("Output Files:\n")
cat("  - consensus_artifacts.txt: Summary of all artifacts\n")
cat("  - consensus_artifacts_annotated.txt: With genomic features\n")
cat("  - consensus_artifacts.bed: BED file for IGV/blacklist\n")
cat("  - artifact_per_sample_signals.txt: Detailed per-sample view\n")
cat("  - sample_thresholds.txt: Per-sample signal thresholds\n")
cat("  - *.pdf: Visualization plots\n")
cat("=================================================================\n")

sink()

message("\n=== ANALYSIS COMPLETE ===")
message(paste("Output directory:", OUT_DIR))
message("\nKey findings:")
message(paste("  Consensus artifacts:", nrow(consensus_artifacts)))
message(paste("  With Ty elements:", sum(consensus_artifacts$has_ty)))
message(paste("  With ARS elements:", sum(consensus_artifacts$has_ars)))
message(paste("  With repeats:", sum(consensus_artifacts$has_repeat)))
message(paste("  With tRNA:", sum(consensus_artifacts$has_trna)))
message(paste("  With rRNA:", sum(consensus_artifacts$has_rrna)))
message(paste("  No annotation:", sum(consensus_artifacts$overlapping_features == "")))
message("\nReview files:")
message(paste("  1.", file.path(OUT_DIR, "CONSENSUS_ANALYSIS_REPORT.txt")))
message(paste("  2.", file.path(OUT_DIR, "consensus_artifacts_annotated.txt")))
message(paste("  3.", file.path(OUT_DIR, "artifact_heatmap.pdf")))
