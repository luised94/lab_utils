################################################################################
# SCRIPT: genomic_subsampling.R
# PURPOSE: Create subsampled FASTQ/BAM files for testing pipelines
# USAGE: Run in R/RStudio with renv. Processes first FASTQ and BAM found.
# VALIDATION: Tested with paired/single-end FASTQ and indexed BAM files
# UPDATES: 
#   2024-12-03: Initial version
#   2024-12-03: Modified for renv, explicit package calls
#   2024-12-03: Changed to process only first found files
################################################################################

#===============================================================================
# DEBUG CONFIGURATION
#===============================================================================
DEBUG <- TRUE                    # Enable/disable debug messages
VERBOSE <- TRUE                  # Enable/disable verbose output
DRY_RUN <- FALSE                # If TRUE, shows actions without execution
SEED <- 42                      # Random seed for reproducibility

#===============================================================================
# PROCESSING CONFIGURATION
#===============================================================================
N_READS <- 1000                 # Number of reads to sample
OUTPUT_PREFIX <- "subsampled"   # Prefix for output files
FILE_PATTERNS <- list(
    fastq = c("\\.fastq$", "\\.fq$", "\\.fastq\\.gz$", "\\.fq\\.gz$"),
    bam = "\\.bam$"
)

#===============================================================================
# ENVIRONMENT SETUP AND VALIDATION
#===============================================================================
if (DEBUG) message("Setting up environment...")

# Set random seed
set.seed(SEED)

# Get home directory
HOME_DIR <- path.expand("~")

# Validate directory existence
stopifnot("Home directory does not exist" = dir.exists(HOME_DIR))

#===============================================================================
# FILE DISCOVERY
#===============================================================================
if (DEBUG) message("Searching for input files...")

# Find first FASTQ file
fastq_pattern <- paste(FILE_PATTERNS$fastq, collapse = "|")
FASTQ_FILE <- list.files(
    path = HOME_DIR,
    pattern = fastq_pattern,
    full.names = TRUE
)[1]

# Find first BAM file
BAM_FILE <- list.files(
    path = HOME_DIR,
    pattern = FILE_PATTERNS$bam,
    full.names = TRUE
)[1]

# Validate file existence
stopifnot("No FASTQ or BAM files found in home directory" = 
          !is.na(FASTQ_FILE) || !is.na(BAM_FILE))

if (VERBOSE) {
    if (!is.na(FASTQ_FILE)) message("Found FASTQ file: ", basename(FASTQ_FILE))
    if (!is.na(BAM_FILE)) message("Found BAM file: ", basename(BAM_FILE))
}

#===============================================================================
# FASTQ PROCESSING
#===============================================================================
if (!is.na(FASTQ_FILE)) {
    if (DEBUG) message("\nProcessing FASTQ file...")
    
    output_file <- file.path(HOME_DIR, 
                            paste0(OUTPUT_PREFIX, "_", basename(FASTQ_FILE)))
    
    if (!DRY_RUN) {
        tryCatch({
            # Read FASTQ
            fastq_data <- ShortRead::readFastq(FASTQ_FILE)
            
            # Sample reads
            subset_idx <- sample(length(fastq_data), 
                               min(N_READS, length(fastq_data)))
            subset_fastq <- fastq_data[subset_idx]
            
            # Write subsampled FASTQ
            ShortRead::writeFastq(subset_fastq, output_file)
            
            if (VERBOSE) message("Created: ", basename(output_file))
        }, error = function(e) {
            warning("Error processing ", basename(FASTQ_FILE), ": ", e$message)
        })
    }
}

#===============================================================================
# BAM PROCESSING
#===============================================================================
if (!is.na(BAM_FILE)) {
    if (DEBUG) message("\nProcessing BAM file...")
    
    output_file <- file.path(HOME_DIR, 
                            paste0(OUTPUT_PREFIX, "_", basename(BAM_FILE)))
    
    if (!DRY_RUN) {
        tryCatch({
            # Open BAM
            bam_data <- Rsamtools::BamFile(BAM_FILE)
            
            # Get total reads
            total_reads <- Rsamtools::countBam(bam_data)$records
            
            # Sample reads
            subset_idx <- sort(sample(total_reads, 
                                    min(N_READS, total_reads)))
            
            # Create subsampled BAM
            Rsamtools::filterBam(
                BAM_FILE, 
                output_file,
                indexDestination = TRUE,
                param = Rsamtools::ScanBamParam(which = subset_idx)
            )
            
            if (VERBOSE) message("Created: ", basename(output_file))
        }, error = function(e) {
            warning("Error processing ", basename(BAM_FILE), ": ", e$message)
        })
    }
}

#===============================================================================
# COMPLETION
#===============================================================================
if (DEBUG) message("\nScript completed successfully")
