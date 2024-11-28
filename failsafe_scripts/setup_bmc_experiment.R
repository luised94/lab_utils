DEBUG_CONFIG <- list(
    enabled = FALSE,
    verbose = TRUE,
    interactive = TRUE,
    dry_run = TRUE
)

experiment_id <- "241122Bel"
stopifnot(
    "Experiment ID must be a character string" = is.character(experiment_id),
    "Invalid experiment ID format. Expected: YYMMDD'Bel'" = grepl("^\\d{6}Bel$", experiment_id)
)

source("~/lab_utils/failsafe_scripts/bmc_config.R")
stopifnot("Script experiment_id is not the same as CONFIG EXPERIMENT_ID" = experiment_id == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID)
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

# Define the paths directly in a character vector
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

# Create directories with debug output
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

if (DEBUG_CONFIG$verbose) {
    mode <- if (DEBUG_CONFIG$dry_run) "DRY RUN" else "LIVE RUN"
    cat(sprintf("\n[%s] Directory structure for experiment: %s\n", mode, experiment_id))
    cat(sprintf("[%s] Base directory: %s\n", mode, base_dir))
}

cat("Directories created successfully!\n")

# Generate combinations
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
    stop(sprintf("Expected %d samples, got %d", expected, n_samples))
}
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

metadata <- metadata[do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]
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
