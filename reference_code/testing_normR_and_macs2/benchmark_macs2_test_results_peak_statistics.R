# Verify required libraries are available.
library(ggplot2)
required_packages <- c("GenomeInfoDb", "IRanges", "GenomicRanges", "rtracklayer", "ggplot2", "Gviz", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "UpSetR")
# add cosmo after manually installing
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

# File paths and configuration
INPUT_FILES <- list(
    narrow_peaks = "~/macs2_test_results/macs2_test_241010Bel_peaks.narrowPeak",
    broad_peaks = "~/macs2_test_results/macs2_test_100303Bel_peaks.broadPeak",
    narrow_signal = "~/macs2_test_results/macs2_test_241010Bel_treat_pileup.bdg",
    broad_signal = "~/macs2_test_results/macs2_test_100303Bel_treat_pileup.bdg",
    reference_peaks = list.files(
        file.path(Sys.getenv("HOME"), "data", "feature_files"),
        pattern = "eaton_peaks",
        full.names = TRUE
    )[1],
    ref_genome_file <- list.files(
        file.path(Sys.getenv("HOME"), "data", "REFGENS"),
        pattern = "S288C_refgenome.fna",
        full.names = TRUE,
        recursive = TRUE
    )[1]
)
output_dir <- "~/macs2_test_results"
PEAK_WIDTH <- 100  # Width for sequence extraction
MOTIF_MIN_LENGTH <- 6
MOTIF_MAX_LENGTH <- 12

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# Validate file existence
for (file in INPUT_FILES) {
    if (!file.exists(file)) {
        stop(sprintf("File not found: %s", file))
    }
}

# 1. Import data
message("Importing peak files...")
narrow_peaks <- rtracklayer::import(INPUT_FILES$narrow_peaks)
broad_peaks <- rtracklayer::import(INPUT_FILES$broad_peaks)
reference_peaks <- rtracklayer::import(INPUT_FILES$reference_peaks)
# Import signal tracks for visualization
message("Importing signal tracks...")
treat_signal1 <- rtracklayer::import(INPUT_FILES$narrow_signal, format = "BedGraph")
treat_signal2 <- rtracklayer::import(INPUT_FILES$broad_signal, format = "BedGraph")

if (!is.null(reference_peaks)) {
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(reference_peaks) <- paste0(
        "chr",
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(reference_peaks)))
    )
}

# Find a reference origin on chrX for focused view
chromosome_to_plot <- "chrX"
chrX_refs <- reference_peaks[GenomicRanges::seqnames(reference_peaks) == chromosome_to_plot]
example_origin <- chrX_refs[1]
origin_start <- example_origin@ranges@start
origin_end <- origin_start + example_origin@ranges@width

# Expand region around origin for visualization
view_window <- 5000  # 5kb window
region_start <- max(0, origin_start - view_window)
region_end <- origin_end + view_window

################################################################################
# Load reference genome
################################################################################
stopifnot(
    "No reference genome found. One expected." = length(ref_genome_file) == 1
)

genome_data <- Biostrings::readDNAStringSet(INPUT_FILES$ref_genome_file)

# Convert genome data to required format
genome_info <- data.frame(
    chrom = names(genome_data),
    size = genome_data@ranges@width
)

# Validation
stopifnot(
    "genome_info must be a data frame" = is.data.frame(genome_info),
    "genome_info must have exactly 16 rows and 2 columns" = all(dim(genome_info) == c(16, 2)),
    "Column names must be 'chrom' and 'size'" = all(colnames(genome_info) == c("chrom", "size")),
    "Chromosome sizes must be positive" = all(genome_info$size > 0),
    "Chromosome names must start with 'chr'" = all(grepl("^chr", genome_info$chrom)),
    "No missing values allowed" = !any(is.na(genome_info))
)

################################################################################
# Basic statistics
################################################################################
message("Calculating basic statistics...")
peak_stats <- data.frame(
    Dataset = c("With_Control", "No_Control"),
    Total_Peaks = c(length(narrow_peaks), length(broad_peaks)),
    Mean_Width = c(mean(narrow_peaks@ranges@width), mean(broad_peaks@ranges@width)),
    Median_Score = c(median(narrow_peaks$score), median(broad_peaks$score))
)

# Chromosome distribution
chr_stats_narrow <- table(GenomicRanges::seqnames(narrow_peaks))
chr_stats_broad <- table(GenomicRanges::seqnames(broad_peaks))

# 6. Overlap Analysis (Improved with UpSet plot)
message("Calculating overlaps...")
narrow_overlaps <- GenomicRanges::findOverlaps(narrow_peaks, reference_peaks)
broad_overlaps <- GenomicRanges::findOverlaps(broad_peaks, reference_peaks)

# Plot score distribution for non-overlapping and overlapping peaks.
narrow_peaks$category <- "Non-overlapping"
narrow_peaks$category[S4Vectors::queryHits(narrow_overlaps)] <- "Overlapping with reference"
score_df <- data.frame(
    score = narrow_peaks$score,
    category = narrow_peaks$category,
    dataset = "241010Bel"
)

pdf(file.path(output_dir, "narrow_peaks_peak_score_dsitribution_by_reference_overlap.pdf"))
ggplot(score_df, aes(x = score, fill = category)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~dataset) +
    theme_minimal() +
    labs(title = "Peak Score Distribution by Reference Overlap")
dev.off()

broad_peaks$category <- "Non-overlapping"
broad_peaks$category[S4Vectors::queryHits(broad_overlaps)] <- "Overlapping with reference"
score_df <- data.frame(
    score = broad_peaks$score,
    category = broad_peaks$category,
    dataset = "241010Bel"
)

pdf(file.path(output_dir, "broad_peaks_peak_score_dsitribution_by_reference_overlap.pdf"))
ggplot(score_df, aes(x = score, fill = category)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~dataset) +
    theme_minimal() +
    labs(title = "Peak Score Distribution by Reference Overlap")
dev.off()

overlap_stats <- data.frame(
    Dataset = c("With_Control", "No_Control"),
    Total_Peaks = c(length(narrow_peaks), length(broad_peaks)),
    Reference_Peaks = c(length(reference_peaks), length(reference_peaks)),
    Overlapping_Peaks = c(
        length(unique(S4Vectors::queryHits(narrow_overlaps))),
        length(unique(S4Vectors::queryHits(broad_overlaps)))
    ),
    Percent_Overlap = c(
        100 * length(unique(S4Vectors::queryHits(narrow_overlaps))) / length(narrow_peaks),
        100 * length(unique(S4Vectors::queryHits(broad_overlaps))) / length(broad_peaks)
    )
)

# UpSet plot for overlap visualization (Narrow Peaks)
narrow_peaks_in_control <- rep(FALSE, length(narrow_peaks))
narrow_peaks_in_control[S4Vectors::queryHits(narrow_overlaps)] <- TRUE
narrow_peak_df <- data.frame(Peak = seq_along(narrow_peaks), InControl = narrow_peaks_in_control)

narrow_peak_df$InControl <- as.character(narrow_peak_df$InControl)

# Create a list of sets for UpSetR
peak_sets <- list(
    "In Control" = narrow_peak_df$Peak[narrow_peak_df$InControl == "TRUE"],
    "Not In Control" = narrow_peak_df$Peak[narrow_peak_df$InControl == "FALSE"]
)

# Create the UpSetR plot
pdf(file.path(output_dir, "narrow_peaks_upset.pdf"))
UpSetR::upset(
    UpSetR::fromList(peak_sets),
    main.bar.color = "#CC0000",
    sets.bar.color = "#CC0000",
    matrix.color = "#CC0000",
    mainbar.y.label = "Number of Peaks",
    sets.x.label = "Total Number of Peaks",
    text.scale = c(1.5, 1.5, 1.3, 1, 1, 0.75)
)
dev.off()

# UpSet plot for overlap visualization (Broad Peaks)
broad_peaks_in_control <- rep(FALSE, length(broad_peaks))
broad_peaks_in_control[S4Vectors::queryHits(broad_overlaps)] <- TRUE
broad_peak_df <- data.frame(Peak = seq_along(broad_peaks), InControl = broad_peaks_in_control)
# Create a proper format for UpSetR
# Convert the logical column into a factor or character
broad_peak_df$InControl <- as.character(broad_peak_df$InControl)

# Create a list of sets for UpSetR
peak_sets <- list(
    "In Control" = broad_peak_df$Peak[broad_peak_df$InControl == "TRUE"],
    "Not In Control" = broad_peak_df$Peak[broad_peak_df$InControl == "FALSE"]
)

# Create the UpSetR plot
pdf(file.path(output_dir, "broad_peaks_upset.pdf"))
UpSetR::upset(
    UpSetR::fromList(peak_sets),
    main.bar.color = "#CC0000",
    sets.bar.color = "#CC0000",
    matrix.color = "#CC0000",
    mainbar.y.label = "Number of Peaks",
    sets.x.label = "Total Number of Peaks",
    text.scale = c(1.5, 1.5, 1.3, 1, 1, 0.75)
)
dev.off()

# Print results
message("Results summary:")
print(peak_stats)
print(chr_stats_narrow)
print(chr_stats_broad)
print(overlap_stats)
