################################################################################
# SCRIPT: align_subsampled_fastq.R
# PURPOSE: Align subsampled FASTQ and create indexed BAM file
# USAGE: Run after genomic_subsampling.R
# VALIDATION: Tested with single-end FASTQ files
# REQUIREMENTS: Reference genome FASTA file in home directory
# UPDATES: 
#   2024-12-03: Initial version
################################################################################

#===============================================================================
# DEBUG CONFIGURATION
#===============================================================================
DEBUG <- TRUE                    # Enable/disable debug messages
VERBOSE <- TRUE                  # Enable/disable verbose output
DRY_RUN <- FALSE                # If TRUE, shows actions without execution

#===============================================================================
# PROCESSING CONFIGURATION
#===============================================================================
# File patterns
INPUT_PATTERN <- "^subsampled_.*\\.fastq$"
REFERENCE_GENOME <- "sacCer3.fa"  # Update with your reference genome filename

# Output configuration
OUTPUT_PREFIX <- "aligned"
THREADS <- 4                      # Number of threads for alignment

#===============================================================================
# ENVIRONMENT SETUP AND VALIDATION
#===============================================================================
if (DEBUG) message("Setting up environment...")

# Get home directory
HOME_DIR <- path.expand("~")

# Validate directory and files
stopifnot(
    "Home directory does not exist" = dir.exists(HOME_DIR),
    "Reference genome not found" = file.exists(file.path(HOME_DIR, REFERENCE_GENOME))
)

#===============================================================================
# FILE DISCOVERY
#===============================================================================
if (DEBUG) message("Searching for subsampled FASTQ file...")

# Find subsampled FASTQ file
FASTQ_FILE <- list.files(
    path = HOME_DIR,
    pattern = INPUT_PATTERN,
    full.names = TRUE
)[1]

stopifnot("No subsampled FASTQ file found" = !is.na(FASTQ_FILE))

if (VERBOSE) {
    message("Found FASTQ file: ", basename(FASTQ_FILE))
}

#===============================================================================
# REFERENCE INDEXING
#===============================================================================
if (DEBUG) message("\nChecking reference index...")

REF_PATH <- file.path(HOME_DIR, REFERENCE_GENOME)
INDEX_PREFIX <- file.path(HOME_DIR, "reference_index")

if (!file.exists(paste0(INDEX_PREFIX, ".00.b.tab"))) {
    if (VERBOSE) message("Building reference index...")
    if (!DRY_RUN) {
        Rsubread::buildindex(
            basename = INDEX_PREFIX,
            reference = REF_PATH,
            memory = 8000
        )
    }
}

#===============================================================================
# ALIGNMENT
#===============================================================================
if (DEBUG) message("\nPerforming alignment...")

# Set output names
SAM_OUTPUT <- file.path(HOME_DIR, paste0(OUTPUT_PREFIX, ".sam"))
BAM_OUTPUT <- file.path(HOME_DIR, paste0(OUTPUT_PREFIX, ".bam"))
SORTED_BAM <- file.path(HOME_DIR, paste0(OUTPUT_PREFIX, "_sorted.bam"))

if (!DRY_RUN) {
    # Perform alignment
    Rsubread::align(
        index = INDEX_PREFIX,
        readfile1 = FASTQ_FILE,
        output_file = SAM_OUTPUT,
        nthreads = THREADS,
        unique = TRUE,
        indels = 5
    )
    
    # Convert SAM to BAM
    if (VERBOSE) message("Converting SAM to BAM...")
    
    Rsamtools::asBam(
        file = SAM_OUTPUT,
        destination = BAM_OUTPUT,
        indexDestination = FALSE
    )
    
    # Sort BAM
    if (VERBOSE) message("Sorting BAM file...")
    
    Rsamtools::sortBam(
        file = BAM_OUTPUT,
        destination = sub("\\.bam$", "", SORTED_BAM)
    )
    
    # Index sorted BAM
    if (VERBOSE) message("Indexing sorted BAM...")
    
    Rsamtools::indexBam(SORTED_BAM)
    
    # Clean up intermediate files
    if (file.exists(SAM_OUTPUT)) file.remove(SAM_OUTPUT)
    if (file.exists(BAM_OUTPUT)) file.remove(BAM_OUTPUT)
}

#===============================================================================
# COMPLETION
#===============================================================================
if (DEBUG) message("\nAlignment pipeline completed successfully")
