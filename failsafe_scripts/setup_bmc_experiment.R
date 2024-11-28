################################################################################
# BMC Experiment Setup Script
################################################################################
#
# PURPOSE:
#   Creates standardized directory structure and metadata files for BMC ChIP-seq
#   experiments, including sample tracking and submission documents.
#
# USAGE:
#   1. Use '/!!' in vim/neovim to jump to required updates
#   2. Update experiment_id (format: YYMMDD'Bel', e.g., "241122Bel")
#   3. Set DEBUG_CONFIG options as needed
#   4. Source or run script
#
# !! ----> REQUIRED UPDATES:
# !! experiment_id <- "241122Bel"
# !! DEBUG_CONFIG <- list(
# !!      enabled = FALSE,
# !!      verbose = TRUE,
# !!      interactive = TRUE,
# !!      dry_run = FALSE
# !!  )

# INPUTS:
#   - experiment_id: 9-character experiment identifier
#   - bmc_config.R: External configuration defining experimental design
#
# OUTPUTS:
#   1. Directory Structure:
#      ~/data/[experiment_id]/
#      +-- peak/
#      +-- fastq/
#      |   +-- raw/
#      |   +-- processed/
#      +-- alignment/
#      +-- bigwig/
#      +-- plots/
#      +-- documentation/
#
#   2. Files:
#      - [experiment_id]_sample_grid.csv: Complete experimental design
#      - [experiment_id]_bmc_table.tsv: BMC submission metadata
#      - [experiment_id]_bmc_config.R: Configuration snapshot
#
# CONTROLS:
#   DEBUG_CONFIG$dry_run    = TRUE   # Preview without creating files
#   DEBUG_CONFIG$verbose    = TRUE   # Show detailed progress
#   DEBUG_CONFIG$interactive = TRUE  # Confirm before proceeding
#
# DEPENDENCIES:
#   - R base packages only
#   - ~/lab_utils/failsafe_scripts/bmc_config.R
#
# COMMON ISSUES:
#   1. Wrong experiment ID format -> Check YYMMDD pattern
#   2. Unexpected sample count -> Review antibody distribution
#   3. File access denied -> Check ~/data permissions
#
# AUTHOR: Luis
# DATE: 2024-11-27
# VERSION: 2.0.0
#
################################################################################

################################################################################
# Configuration and Debug Settings
################################################################################
# !! Review debug configuration
DEBUG_CONFIG <- list(
    enabled = FALSE,
    verbose = TRUE,
    interactive = TRUE,
    dry_run = FALSE
)

################################################################################
# Experiment ID Validation
################################################################################
# !! Update experiment ID
experiment_id <- "241122Bel"
stopifnot(
    "Experiment ID must be a character string" = is.character(experiment_id),
    "Invalid experiment ID format. Expected: YYMMDD'Bel'" = grepl("^\\d{6}Bel$", experiment_id)
)

################################################################################
# Load and Validate Experiment Configuration
################################################################################
bmc_configuration_definition_path <- "~/lab_utils/failsafe_scripts/bmc_config.R"
source(bmc_configuration_definition_path)
stopifnot("Script experiment_id is not the same as CONFIG EXPERIMENT_ID" = experiment_id == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID)

################################################################################
# Directory Setup and User Confirmation
################################################################################
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
if (DEBUG_CONFIG$interactive) {
    cat(sprintf("\nExperiment ID: %s\n", experiment_id))
    cat("Base directory will be:", base_dir, "\n")

    user_response <- readline("Continue with this experiment? (y/n): ")
    if (tolower(user_response) != "y") {
        stop("Script terminated by user")
    }
    cat("Proceeding with directory creation...\n\n")
}

################################################################################
# Directory Structure Definition and Creation
################################################################################
# Define directory structure
data_directories <- c(
    "peak",
    "fastq/raw",
    "fastq/processed",
    "alignment",
    "bigwig",
    "plots/genome_tracks/overview",
    "plots/genome_tracks/experimental_comparisons",
    "documentation/dna_qc_traces"
)

# Create directory structure
full_paths <- file.path(base_dir, data_directories)
invisible(lapply(full_paths, function(path) {
    if (DEBUG_CONFIG$dry_run) {
        cat(sprintf("[DRY RUN] Would create directory: %s\n", path))
    } else {
        dir_created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (DEBUG_CONFIG$verbose) {
            status <- if (dir_created) "Created" else "Already exists"
            cat(sprintf("[%s] %s\n", status, path))
        }
    }
}))

# Report directory creation status
if (DEBUG_CONFIG$verbose) {
    mode <- if (DEBUG_CONFIG$dry_run) "DRY RUN" else "LIVE RUN"
    cat(sprintf("\n[%s] Directory structure for experiment: %s\n", mode, experiment_id))
    cat(sprintf("[%s] Base directory: %s\n", mode, base_dir))
}

cat("Directories created successfully!\n")

################################################################################
# Sample Metadata Generation and Validation
################################################################################
# Generate experimental combinations
metadata <- do.call(expand.grid, EXPERIMENT_CONFIG$CATEGORIES)

# Filter invalid combinations
invalid_idx <- Reduce(
    `|`,
    lapply(EXPERIMENT_CONFIG$INVALID_COMBINATIONS, eval, envir = metadata)
)
metadata <- subset(metadata, !invalid_idx)

# Apply experimental conditions
valid_idx <- Reduce(
    `|`,
    lapply(EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS, eval, envir = metadata)
)
metadata <- subset(metadata, valid_idx)

# Verify sample count
n_samples <- nrow(metadata)
expected <- EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
if (n_samples != expected) {
    # Print diagnostic information
    cat("\nDiagnostic Information:\n")
    cat("----------------------\n")
    print(table(metadata$antibody))  # Show antibody distribution
    cat("\nFull sample breakdown:\n")
    print(summary(metadata))         # Show all category distributions
    cat("\n")

    stop(sprintf("Expected %d samples, got %d", expected, n_samples))
}

################################################################################
# Metadata Formatting and Organization
################################################################################
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

# Sort metadata according to column order
metadata <- metadata[do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]

# Generate sample names
metadata$full_name <- apply(metadata, 1, paste, collapse = "_")
metadata$short_name <- apply(metadata[, EXPERIMENT_CONFIG$COLUMN_ORDER], 1,
    function(x) paste0(substr(x, 1, 1), collapse = ""))

################################################################################
# BMC Metadata Generation
################################################################################
bmc_metadata <- data.frame(
    SampleName = metadata$full_name,
    Vol_uL = 10,
    Conc = 0,
    Type = "ChIP",
    Genome = "Saccharomyces cerevisiae",
    Notes = ifelse(
        metadata$antibody == "Input",
        "Run on fragment analyzer.",
        "Run on femto pulse."
    ),
    Pool = "A",
    stringsAsFactors = FALSE
)

################################################################################
# File Output Generation
################################################################################
# Define output file paths
sample_grid_path <- file.path(base_dir, "documentation",
                             paste0(experiment_id,"_", "sample_grid.csv"))

bmc_table_path <- file.path(base_dir, "documentation",
                           paste0(experiment_id,"_", "bmc_table.tsv"))

bmc_experiment_config_path <- file.path(base_dir, "documentation",
                           paste0(experiment_id,"_", "bmc_config.R"))

# Handle file writing with dry run checks
if (DEBUG_CONFIG$dry_run) {
    cat(sprintf("[DRY RUN] Would write sample grid to: %s\n", sample_grid_path))
    cat(sprintf("[DRY RUN] Would write BMC table to: %s\n", bmc_table_path))
    cat(sprintf("[DRY RUN] Would write BMC config script: %s\n", bmc_experiment_config_path))
} else {
    # Write sample grid file
    if (file.exists(sample_grid_path)) {
        if (DEBUG_CONFIG$verbose) {
            cat(sprintf("[SKIP] Sample grid file already exists: %s\n", sample_grid_path))
        }
    } else {
        write.csv(
            metadata,
            file = sample_grid_path,
            row.names = FALSE
        )
        if (DEBUG_CONFIG$verbose) {
            cat(sprintf("[WROTE] Sample grid to: %s\n", sample_grid_path))
        }
    }

    # Write BMC table file
    if (file.exists(bmc_table_path)) {
        if (DEBUG_CONFIG$verbose) {
            cat(sprintf("[SKIP] BMC table file already exists: %s\n", bmc_table_path))
        }
    } else {
        write.table(
            bmc_metadata,
            file = bmc_table_path,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )
        if (DEBUG_CONFIG$verbose) {
            cat(sprintf("[WROTE] BMC table to: %s\n", bmc_table_path))
        }
    }

    # Copy BMC Experiment Config
    if (file.exists(bmc_experiment_config_path)) {
        if (DEBUG_CONFIG$verbose) {
            cat(sprintf("[SKIP] BMC config file already exists: %s\n", bmc_configuration_definition_path))
        }
    } else {
        file.copy(bmc_configuration_definition_path, to = bmc_experiment_config_path)
        if (DEBUG_CONFIG$verbose) {
            cat(sprintf("[WROTE] BMC table to: %s\n", bmc_table_path))
        }
    }
}

##' Add to existing configurations
#CONFIG <- list(
#    EXPERIMENT = list(
#        PATHS = list(
#            DROPBOX = "/mnt/c/Users/%s/Dropbox (MIT)/",
#            CONFIG = "sampleGridConfig.R"
#        ),
#
#        DIRECTORIES = c(
#            "peak",
#            "fastq",
#            "alignment",
#            "qualityControl",
#            "bigwig",
#            "plots",
#            "documentation"
#        ),
#
#        OUTPUT = list(
#            ENABLED = FALSE,
#            CONFIG_TEMPLATE = "%s_%s.R"
#        )
#    ),
#
#    ENVIRONMENT = list(
#        REQUIRED_VARS = c(
#            "WINDOWS_USER"
#        )
#    )
#)
#CONFIG <- list(
#    ID_REGEX = "\\d{5,6}",
#    FILE_PATTERNS = list(
#        FASTQ = "*.fastq",
#        SAMPLE_TABLE = "sample_table"
#    ),
#    PATHS = list(
#        BASE_DATA = file.path(Sys.getenv("HOME"), "data"),
#        DOCUMENTATION = "documentation",
#        FASTQ = "fastq"
#    )
#)
