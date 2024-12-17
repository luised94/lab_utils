################################################################################
# Plot eaton genomic tracks
################################################################################
# PURPOSE: View control data using same few to see noise.
# Conclusion: See HM1108 of Eaton 2010 paper to see noise and alignment to reference bed file.
# Doesnt look particularly different but less noise definitely although mine seem to have more signal.
# USAGE: 
# !! ----> REQUIRED UPDATES: experiment_id, input files paths
# STRUCTURE: Configuration -> Validation -> Processing -> Output
# VALIDATION: Input file existence, format compatibility, memory requirements
# DEPENDENCIES: normR, GenomicRanges, rtracklayer
# COMMON ISSUES: Memory limitations for large BAM files, sparse control samples
# OUTPUT: BED files with called peaks, quality metrics
# AUTHOR: [Your Name]
# DATE: [Current Date]
# VERSION: 1.0.0
################################################################################
################################################################################
# Configuration and Debug Settings
################################################################################
DEBUG_CONFIG <- list(
    single_file_mode = TRUE,
    verbose = TRUE,
    dry_run = TRUE,
    files_to_process_idx = 7
)

NORMR_CONFIG <- list(
    # Core analysis parameters
    bin_size = 100L,    # Base pairs, adjusted for yeast genome
    min_mapq = 30L,      # Minimum mapping quality for high-quality alignments
    iterations = 10L,    # Number of EM iterations
    processors = 1L,     # Number of processors to use
    paired_end = FALSE,  # Specify single-end data
    
    # Peak calling thresholds
    default_fdr = 10^-6,                    # Default threshold for significant peaks
    min_enrichment_score = 1.5,            # Minimum fold change
    min_reads_per_bin = 1,                 # Minimum reads required per bin
    
    # S. cerevisiae ORC-specific parameters
    expected_peak_range = list(
        min = 100,
        max = 500,
        typical = 250:400
    ),
    
    # Genome specifications
    expected_chromosomes = 16,
    chromosome_prefix = "chr",
    chromosome_adjustment = TRUE,  # Round chromosome sizes to bin size multiples
    
    # File patterns and templates
    bam_pattern = "consolidated_([0-9]{5,6})_sequence_to_S288C_sorted\\.bam$",
    genome_pattern = "S288C_refgenome.fna",
    
    # Output templates
    output_name_template = "%s_peaks_%s_vs_%s_%s.bed",        # timestamp, chip, input, package
    region_file_template = "%s_regions_%s_vs_%s_%s.tsv",
    bedgraph_template = "%s_enrichment_%s_vs_%s_%s.bedGraph",
    
    # Output specifications
    region_columns = c(
        "chromosome", "start", "end",
        "treatment_count", "control_count",
        "enrichment", "qvalue", "peak_class"
    )
)

# Time formatting configuration
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d"               # YYYYMMDD
)

# Generate timestamps once at script start
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)

################################################################################
# Load Required Libraries
################################################################################
required_packages <- c("IRanges", "normr", "GenomicRanges", "rtracklayer", "ggplot2", "Gviz", "stats", "boot")
# add cosmo after manually installing
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

################################################################################
# Directory Setup
################################################################################
experiment_id <- "100303Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
bam_dir <- file.path(base_dir, "alignment")
peak_dir <- file.path(base_dir, "peak")
metadata_path <- file.path(base_dir, "documentation",
                          paste0(experiment_id, "_sample_grid.csv"))

config_path <- file.path(base_dir, "documentation", paste0(experiment_id, "_bmc_config.R"))
# Validate directories
stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "BAM directory does not exist" = dir.exists(bam_dir),
    "Metadata path file does not exist" = file.exists(metadata_path)
)

# Create output directory if it doesn't exist
if (!dir.exists(peak_dir)) {
    dir.create(peak_dir, recursive = TRUE, showWarnings = FALSE)
}

################################################################################
# Load and Validate Experiment Configuration and Dependencies
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/core_scripts/functions_for_file_operations.R",
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)

# Define required dependencies
required_modules <- list(
    list(
        path = "~/lab_utils/core_scripts/functions_for_logging.R",
        description = "Logging functions",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_metadata_processing.R",
        description = "Process metadata grid for downstream analysis.",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_peak_calling.R",
        description = "Load functions required for peak calling and benchmarking script",
        required = TRUE
    ),
    list(
        path = config_path,
        description = "Load EXPERIMENT CONFIG object.",
        required = TRUE
    )
)

# Validate module structure
stopifnot(
    "modules must have required fields" = all(sapply(required_modules, function(m) {
        all(c("path", "description", "required") %in% names(m))
    }))
)

# Load dependencies with status tracking
load_status <- lapply(required_modules, function(module) {
    if (DEBUG_CONFIG$verbose) {
        cat(sprintf("\n[LOADING] %s\n", module$description))
    }

    success <- safe_source(module$path, verbose = TRUE)

    if (!success && module$required) {
        stop(sprintf(
            "[FATAL] Failed to load required module: %s\n  Path: %s",
            module$description, module$path
        ))
    } else if (!success) {
        warning(sprintf(
            "[WARNING] Optional module not loaded: %s\n  Path: %s",
            module$description, module$path
        ))
    }

    return(list(
        module = module$description,
        path = module$path,
        loaded = success
    ))
})

# Display loading summary using ASCII
if (DEBUG_CONFIG$verbose) {
    cat("\n=== Module Loading Summary ===\n")
    invisible(lapply(load_status, function(status) {
        cat(sprintf(
            "[%s] %s\n    Path: %s\n",
            if(status$loaded) "+" else "-",
            status$module,
            status$path
        ))
    }))
}

if (DEBUG_CONFIG$verbose) {
    print_config_settings(DEBUG_CONFIG)
}

################################################################################
# File Discovery and Validation
################################################################################
bam_files <- list.files(
    path = bam_dir,
    pattern = NORMR_CONFIG$bam_pattern,
    recursive = TRUE,
    full.names = TRUE
)

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = NORMR_CONFIG$bam_pattern,
    replacement = "\\1",
    x = basename(bam_files)
)

# Find bigwig files
bigwig_pattern <- sprintf("_%s\\.bw$", EXPERIMENT_CONFIG$NORMALIZATION$active)
bigwig_files <- list.files(
    file.path(base_dir, "coverage"),
    pattern = bigwig_pattern,
    full.names = TRUE
)
normalization_method <- sub(".*_([^_]+)\\.bw$", "\\1",
                          basename(bigwig_files[1]))

stopifnot(
    "No BAM files found" = length(bam_files) > 0,
    "No bigwig files found" = length(bigwig_files) > 0,
    "No sample_ids found" = length(sample_ids) > 0,
    "Number of sample ids must match number of files" = length(bam_files) == length(sample_ids)
)

################################################################################
# Load and process sample metadata
################################################################################
# Load and process metadata
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Enforce factor levels from config
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
    if (col_name %in% colnames(metadata)) {
        metadata[[col_name]] <- factor(
            metadata[[col_name]],
            levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
            ordered = TRUE
        )
    }
}

# Sort metadata using config column order
sorted_metadata <- metadata[do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]

# Add sample IDs to metadata
sorted_metadata$sample_id <- sample_ids


# Add after processing metadata but before track creation
short_sample_ids <- create_minimal_identifiers(
    sorted_metadata$sample_id,
    verbose = DEBUG_CONFIG$verbose
)

# Create mapping between full and short IDs
sample_id_mapping <- setNames(short_sample_ids, sorted_metadata$sample_id)

################################################################################
# Load reference genome
################################################################################
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = NORMR_CONFIG$genome_pattern,
    full.names = TRUE,
    recursive = TRUE
)[1]

stopifnot(
    "No reference genome found. One expected." = length(ref_genome_file) == 1
)

genome_data <- Biostrings::readDNAStringSet(ref_genome_file)

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
# Load feature file (annotation)
################################################################################
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]

if (!is.null(feature_file)) {
    features <- rtracklayer::import(feature_file)
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(features) <- paste0(
        "chr",
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(features)))
    )
}

################################################################################
# Process given control file.
################################################################################
# Set the appropriate control for the given experiment.
# Eaton data is misordered compared to metadata I remembered. The correct number is 4.
control_idx <- 4
chromosome_to_plot <- 10
chromosome_width <- genome_data[chromosome_to_plot]@ranges@width
chromosome_roman <- paste0("chr", utils::as.roman(chromosome_to_plot))
#TODO: Add interactive that displays the chosen control to confirm benchmarking.
print("Sample information as positive control")
print(sorted_metadata[control_idx, ])
# Find control sample
#control_sample <- find_control_sample(
#    experimental_sample = sorted_metadata[control_idx, ],
#    metadata = sorted_metadata,
#    control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS
#)
# Get sample IDs
chip_id <- sorted_metadata[control_idx, "sample_id"]
track_name <- sprintf(
    "%s: %s - %s",
    sample_id_mapping[chip_id],
    sorted_metadata$short_name[control_idx],
    sorted_metadata$antibody[control_idx]
)
if (DEBUG_CONFIG$verbose) {
    message(sprintf("Processing sample: %s", chip_id))
    message(sprintf("Track name for visualization: %s", track_name))
}
# Find BAM files using grepl for more robust matching
chip_bam <- bam_files[grepl(paste0("consolidated_", chip_id, "_sequence_to_S288C_sorted\\.bam$"), bam_files)]
bigwig_file <- bigwig_files[grepl(chip_id, bigwig_files)][1]

# Validate BAM files
stopifnot(
    "ChIP BAM file not found" = length(chip_bam) == 1,
    "ChIP BAM file does not exist" = file.exists(chip_bam),
    "ChIP bigwig file does not exist" = file.exists(bigwig_file)
)
# Generate output filename using short IDs
output_filename <- sprintf(
    NORMR_CONFIG$output_name_template,
    TIMESTAMPS$full,
    sample_id_mapping[chip_id],
    "none",
    "normr"
)
output_path <- file.path(peak_dir, output_filename)
if (DEBUG_CONFIG$verbose) {
    message(sprintf("Output will be written to: %s", output_path))
}
# Perform peak calling
tryCatch({
    if (DEBUG_CONFIG$verbose) {
        message("\nStarting peak calling...")
        message(sprintf("Genome size: %d bp across %d chromosomes",
                    sum(genome_info$size),
                    nrow(genome_info)))
    }
    if (DEBUG_CONFIG$verbose) {
        message("\nPre-processing count data...")
    }
    tryCatch({

    track_data <- rtracklayer::import(bigwig_file)

    filtered_tracks <- list(
        Gviz::GenomeAxisTrack(),
        Gviz::DataTrack(
            track_data[GenomicRanges::seqnames(track_data) == chromosome_roman], 
            name = track_name, 
            type = "l"
        ),
        Gviz::AnnotationTrack(
            features, 
            name = "eaton peaks"
        )
    )
    png("my_plot_track.png")
    Gviz::plotTracks(
        filtered_tracks,
        chromosome = chromosome_roman,
    )
    dev.off()
    }, error = function(e) {
      print(paste("Error plotting tracks:", e$message))
    })
}, error = function(e) {
    message(sprintf("Error processing sample %s: %s",
                    chip_id,
                    e$message))
})
