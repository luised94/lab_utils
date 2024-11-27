
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
create_experiment_dir <- function(directory_path, subdirectories){
    if (!dir.exists(directory_path)) {
        dir.create(directory_path, recursive = TRUE)
    } else {
        cat(sprintf("Directory %s already exists.\n", directory_path))
    }
    sapply(file.path(directory_path, subdirectories), function(path_to_create) {
        if (!dir.exists(path_to_create)) {
            dir.create(path_to_create)
        } else {
            cat(sprintf("Directory %s already exists.\n", path_to_create))
        }
    })
}
