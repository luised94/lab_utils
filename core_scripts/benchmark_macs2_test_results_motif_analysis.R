path_to_meme_binary <- file.path(Sys.getenv("HOME"), "meme/bin/meme")
options(meme.bin = path_to_meme_binary)
system("export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.3.3:$PATH")
# Test MEME installation and setup
meme_test <- system("meme -version", intern = TRUE)
message("MEME version: ", meme_test)

# Verify MEME setup
if (is.null(getOption("meme.bin"))) {
    stop("MEME binary location not set. Use options(meme.bin = path_to_meme)")
}
# Verify required libraries are available.
# Environment validation
if (system("hostname", intern = TRUE) == "luria") {
    message("universalmotif package is used in this script. To install universalmotif, some packages\n")
    message("are not available on luria.")
    stop("This script should not be run on luria cluster")
}

# MEME parameters specific for ORC analysis
MEME_PARAMS <- list(
    minw = 8,          # ORC binding sites typically 8-15bp
    maxw = 15,
    nmotifs = 3,       # Start with top 3 motifs
    mod = "zoops",     # Zero or one occurrence per sequence
    markov_order = 2   # Background model order
)

required_packages <- c("GenomeInfoDb", "IRanges", "GenomicRanges", "rtracklayer", "ggplot2", "Gviz", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "UpSetR", "universalmotif")
# add cosmo after manually installing
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}
library(universalmotif)
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
    ref_genome_file = list.files(
        file.path(Sys.getenv("HOME"), "data", "REFGENS"),
        pattern = "S288C_refgenome.fna",
        full.names = TRUE,
        recursive = TRUE
    )[1]
)
output_dir <- file.path(Sys.getenv("HOME"), "macs2_test_results")
# Create MEME-specific output directories
meme_overlap_dir <- file.path(output_dir, "meme_overlap_output")
meme_ref_dir <- file.path(output_dir, "meme_reference_output")
meme_all_dir <- file.path(output_dir, "meme_all_output")

for (dir in c(meme_overlap_dir, meme_ref_dir)) {
    if (!dir.exists(dir)) {
        dir.create(dir)
    }
}

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

################################################################################
# Load reference genome
################################################################################
stopifnot(
    "No reference genome found. One expected." = length(INPUT_FILES$ref_genome_file) == 1
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
# Sequence Preparation for Motif Analysis
################################################################################
PEAK_WIDTH <- 100  # Width for sequence extraction

# Get overlapping peaks with reference
narrow_overlaps <- GenomicRanges::findOverlaps(narrow_peaks, reference_peaks)
overlapping_peaks <- narrow_peaks[S4Vectors::queryHits(narrow_overlaps)]

message("Preparing sequences for motif analysis...")
# Function to extract centered sequences
get_centered_sequences <- function(peaks, width = PEAK_WIDTH) {
    peak_centers <- peaks@ranges@start + peaks@ranges@width %/% 2
    ranges <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(peaks),
        IRanges::IRanges(peak_centers - (width %/% 2), peak_centers + (width %/% 2)),
        strand = "*"
    )
    return(Biostrings::getSeq(genome_data, ranges))
}

# Extract sequences for different peak sets
overlapping_seqs <- get_centered_sequences(overlapping_peaks)
reference_seqs <- get_centered_sequences(reference_peaks)
all_peaks_seqs <- get_centered_sequences(narrow_peaks)

################################################################################
# MEME Motif Discovery
################################################################################
message("Running MEME motif discovery...")

# Run MEME on overlapping peaks (highest confidence)
meme_results_overlap <- run_meme(
    overlapping_seqs,
    output = meme_overlap_dir,
    overwrite.dir = TRUE,
    minw = MEME_PARAMS$minw,
    maxw = MEME_PARAMS$maxw,
    nmotifs = MEME_PARAMS$nmotifs,
    mod = MEME_PARAMS$mod,
    markov_order = MEME_PARAMS$markov_order,
    revcomp = TRUE,
    readsites = TRUE,
    verbose = 2
)

# Run MEME on reference peaks for comparison
meme_results_ref <- run_meme(
    reference_seqs,
    output = meme_ref_dir,
    overwrite.dir = TRUE,
    minw = MEME_PARAMS$minw,
    maxw = MEME_PARAMS$maxw,
    nmotifs = MEME_PARAMS$nmotifs,
    mod = MEME_PARAMS$mod,
    markov_order = MEME_PARAMS$markov_order,
    revcomp = TRUE,
    readsites = TRUE,
    verbose = 2
)

meme_results_overlap <- run_meme(
    all_peaks_seqs,
    output = meme_all_dir,
    overwrite.dir = TRUE,
    minw = MEME_PARAMS$minw,
    maxw = MEME_PARAMS$maxw,
    nmotifs = MEME_PARAMS$nmotifs,
    mod = MEME_PARAMS$mod,
    markov_order = MEME_PARAMS$markov_order,
    revcomp = TRUE,
    readsites = TRUE,
    verbose = 2
)

# Check MEME output files
meme_output_files <- list.files(meme_overlap_dir, full.names = TRUE)
message("MEME output generated: ", length(meme_output_files), " files")

# Read MEME results using universalmotif
meme_motifs <- universalmotif::read_meme(
    file.path(meme_overlap_dir, "meme.txt")
)

# Basic validation of motif discovery
if (length(meme_motifs) == 0) {
    warning("No motifs were discovered")
} else {
    message("Successfully discovered ", length(meme_motifs), " motifs")
}

# Save motif plots
pdf(file.path(output_dir, "overlap_discovered_motifs.pdf"))
universalmotif::view_motifs(
    meme_motifs,
    show.positions = TRUE
)
dev.off()


# Check MEME output files
meme_output_files <- list.files(meme_ref_dir, full.names = TRUE)
message("MEME output generated: ", length(meme_ref_dir), " files")

# Read MEME results using universalmotif
meme_motifs <- universalmotif::read_meme(
    file.path(meme_ref_dir, "meme.txt")
)

# Basic validation of motif discovery
if (length(meme_motifs) == 0) {
    warning("No motifs were discovered")
} else {
    message("Successfully discovered ", length(meme_motifs), " motifs")
}

# Save motif plots
pdf(file.path(output_dir, "reference_discovered_motifs.pdf"))
universalmotif::view_motifs(
    meme_motifs,
    show.positions = TRUE
)
dev.off()

# Check MEME output files
meme_output_files <- list.files(meme_all_dir, full.names = TRUE)
message("MEME output generated: ", length(meme_all_dir), " files")

# Read MEME results using universalmotif
meme_motifs <- universalmotif::read_meme(
    file.path(meme_allf_dir, "meme.txt")
)

# Basic validation of motif discovery
if (length(meme_motifs) == 0) {
    warning("No motifs were discovered")
} else {
    message("Successfully discovered ", length(meme_motifs), " motifs")
}

# Save motif plots
pdf(file.path(output_dir, "allpeaks_discovered_motifs.pdf"))
universalmotif::view_motifs(
    meme_motifs,
    show.positions = TRUE
)
dev.off()
