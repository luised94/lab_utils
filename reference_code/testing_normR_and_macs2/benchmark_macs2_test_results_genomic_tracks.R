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
# Genomic Tracks
################################################################################
# 4. Create visualization tracks
# Define Ideogram using local genome information.
# Create custom ideogram track
# First, create a proper bands data frame
bands_df <- data.frame(
    chrom = genome_info$chrom,
    chromStart = 0,
    chromEnd = genome_info$size,
    name = paste0(genome_info$chrom, "_band"),
    gieStain = "gpos50"  # Default staining
)

# Create chromosome-specific subset
chromosome_bands <- bands_df[bands_df$chrom == chromosome_to_plot, ]

# Create ideogram track with explicit bands data
itrack <- Gviz::IdeogramTrack(
    chromosome = chromosome_to_plot[1],  # Ensure single value
    genome = "sacCer3",
    bands = bands_df,
    name = "Ideogram"
)

# Plot and save focused views (Narrow Peaks)
pdf(file.path(output_dir, "narrow_peaks_origin_region.pdf"), width=10, height=6)
Gviz::plotTracks(list(
    itrack,
    Gviz::GenomeAxisTrack(),
    Gviz::DataTrack(treat_signal1, name="With Control Signal", type="h", col.histogram="darkblue", fill.histogram="lightblue"),
    Gviz::AnnotationTrack(narrow_peaks, name="Narrow Peaks", col="darkblue"),
    Gviz::AnnotationTrack(reference_peaks, name="Reference Origins", col="darkred")
), chromosome=chromosome_to_plot, from=region_start, to=region_end, main="Narrow Peaks - Origin Region")
dev.off()

# Plot and save whole chromosome views (Narrow Peaks)
pdf(file.path(output_dir, "narrow_peaks_chrX.pdf"), width=15, height=6)
Gviz::plotTracks(list(
    itrack,
    Gviz::GenomeAxisTrack(),
    Gviz::DataTrack(treat_signal1, name="With Control Signal", type="h", col.histogram="darkblue", fill.histogram="lightblue"),
    Gviz::AnnotationTrack(narrow_peaks, name="Narrow Peaks", col="darkblue"),
    Gviz::AnnotationTrack(reference_peaks, name="Reference Origins", col="darkred")
), chromosome=chromosome_to_plot, main="Narrow Peaks - Chromosome X")
dev.off()

# Plot and save focused views (Broad Peaks)
pdf(file.path(output_dir, "broad_peaks_origin_region.pdf"), width=10, height=6)
Gviz::plotTracks(list(
    itrack,
    Gviz::GenomeAxisTrack(),
    Gviz::DataTrack(treat_signal2, name="No Control Signal", type="h", col.histogram="darkgreen", fill.histogram="lightgreen"),
    Gviz::AnnotationTrack(broad_peaks, name="Broad Peaks", col="darkgreen"),
    Gviz::AnnotationTrack(reference_peaks, name="Reference Origins", col="darkred")
), chromosome=chromosome_to_plot, from=region_start, to=region_end, main="Broad Peaks - Origin Region")
dev.off()

# Plot and save whole chromosome views (Broad Peaks)
pdf(file.path(output_dir, "broad_peaks_chrX.pdf"), width=15, height=6)
Gviz::plotTracks(list(
    itrack,
    Gviz::GenomeAxisTrack(),
    Gviz::DataTrack(treat_signal2, name="No Control Signal", type="h", col.histogram="darkgreen", fill.histogram="lightgreen"),
    Gviz::AnnotationTrack(broad_peaks, name="Broad Peaks", col="darkgreen"),
    Gviz::AnnotationTrack(reference_peaks, name="Reference Origins", col="darkred")
), chromosome=chromosome_to_plot, main="Broad Peaks - Chromosome X")
dev.off()
