################################################################################
# Peak Calling using normR for ChIP-seq Analysis
################################################################################
# PURPOSE: Perform peak calling on ChIP-seq data using normR
# USAGE: source("peak_calling_normr.R")
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
    files_to_process_idx = 1
)

PEAK_CALLING_CONFIG <- list(
    min_mapping_quality = 30,
    fdr_threshold = 0.01,
    bin_size = 100,
    peak_merge_distance = 250,
    min_read_count = 10
)

################################################################################
# Load Required Libraries
################################################################################
required_packages <- c("normr", "GenomicRanges", "rtracklayer")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
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
        path = "~/lab_utils/core_scripts/bmc_config.R",
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

################################################################################
# Directory Setup
################################################################################
experiment_id <- "241010Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
bam_dir <- file.path(base_dir, "alignment")
peak_dir <- file.path(base_dir, "peak")
metadata_path <- file.path(base_dir, "documentation",
                          paste0(experiment_id, "_sample_grid.csv"))

# Validate directories
stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "BAM directory does not exist" = dir.exists(bam_dir)
)

# Create output directory if it doesn't exist
if (!dir.exists(peak_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}


################################################################################
# File Discovery and Validation
################################################################################
bam_files <- list.files(
    path = bam_dir,
    pattern = "consolidated_([0-9]{5,6})_sequence_to_S288C_sorted\\.bam$",
    recursive = TRUE,
    full.names = TRUE
)

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = "consolidated_([0-9]{5,6})_sequence_to_S288C_sorted\\.bam",
    replacement = "\\1",
    x = basename(bam_files)
)

stopifnot(
    "No BAM files found" = length(bam_files) > 0,
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
    pattern = "S288C_refgenome.fna",
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
# Process Files
################################################################################
files_to_process <- if (DEBUG_CONFIG$single_file_mode) {
    DEBUG_CONFIG$files_to_process_idx
} else {
    seq_along(bam_files)
}

for (file in files_to_process){
    # find control sample.
    control_sample <- find_control_sample(

        experimental_sample = sorted_metadata[file, ]
        metadta = sorted_metadata,
        control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS
    )
    if (DEBUG_CONFIG$verbose) {
        message(sprintf("Adding control track: %s",
                control_sample$sample_id))
    }
    # run normR with default values.
}
