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
    single_file_mode = FALSE,
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
# Load and Validate Experiment Configuration
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
        path = "~/lab_utils/failsafe_scripts/functions_for_logging.R",
        description = "BMC Configuration",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/failsafe_scripts/functions_for_logging.R",
        description = "BMC Configuration",
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
# Directory Setup and Validation
################################################################################
experiment_id <- "241010Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
bam_dir <- file.path(base_dir, "alignment")
peak_dir <- file.path(base_dir, "peak")

# Validate directories
stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "BAM directory does not exist" = dir.exists(bam_dir)
)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}


################################################################################
# File Discovery and Validation
################################################################################
bam_files <- list.files(
    path = bam_dir,
    pattern = "\\.bam$",
    recursive = TRUE,
    full.names = TRUE
)

stopifnot(
    "No BAM files found" = length(bam_files) > 0,
)

# Load reference genome
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]
stopifnot(
    "No reference genome found. One expected." = length(ref_genome_file) < 1,
    "More than one reference genome found. Only one expected." = length(ref_genome_file) > 1,
)
genome_data <- Biostrings::readDNAStringSet(ref_genome_file)

# Convert genome data to required format
genome_info <- data.frame(
    chrom = names(genome_data),
    size = width(genome_data)
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

}
