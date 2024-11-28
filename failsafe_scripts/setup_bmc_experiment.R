#' Add to existing configurations
CONFIG <- list(
    EXPERIMENT = list(
        PATHS = list(
            DROPBOX = "/mnt/c/Users/%s/Dropbox (MIT)/",
            CONFIG = "sampleGridConfig.R"
        ),
        
        DIRECTORIES = c(
            "peak",
            "fastq",
            "alignment",
            "qualityControl",
            "bigwig",
            "plots",
            "documentation"
        ),
        
        OUTPUT = list(
            ENABLED = FALSE,
            CONFIG_TEMPLATE = "%s_%s.R"
        )
    ),
    
    ENVIRONMENT = list(
        REQUIRED_VARS = c(
            "WINDOWS_USER"
        )
    )
)
CONFIG <- list(
    ID_REGEX = "\\d{5,6}",
    FILE_PATTERNS = list(
        FASTQ = "*.fastq",
        SAMPLE_TABLE = "sample_table"
    ),
    PATHS = list(
        BASE_DATA = file.path(Sys.getenv("HOME"), "data"),
        DOCUMENTATION = "documentation",
        FASTQ = "fastq"
    )
)
DEBUG_CONFIG <- list(
    enabled = FALSE,           # TRUE for testing single group, FALSE for all
    group = 10,               # Which group to process when in debug mode
    samples_per_group = 4,    # Samples per plot
    save_plots = TRUE,       # Whether to save plots to files
    verbose = TRUE,           # Print debug information
    chromosome = 10,
    interactive = FALSE,
    display_time = 2
)

experiment_id <- "241122Bel"
if (!grepl("^\\d{6}Bel$", experiment_id)) {
    stop("Invalid experiment ID format. Expected: YYMMDD'Bel'")
}

base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)

# Define the paths directly in a character vector
data_directories <- c(
    "peak",
    "fastq/raw",
    "fastq/processed",
    "alignment",
    "bigwig",
    "plots/genome_tracks/overview",
    "plots/genome_tracks/experimental_comparisons",
    "documentation/dna_qc_traces",
)

full_paths <- file.path(base_dir, data_directories)
sapply(full_paths, dir.create, recursive = TRUE, showWarnings = FALSE)
cat("Directories created successfully!\n")
