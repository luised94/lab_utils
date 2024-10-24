#' Script Validation Configuration
CONFIG <- list(
    VALIDATION = list(
        TYPES = list(
            STRING = "character",
            INTEGER = "integer",
            NUMERIC = "double",
            LOGICAL = "logical"
        ),
        
        PATHS = list(
            DATA_DIR = file.path(Sys.getenv('HOME'), 'data')
        ),
        
        LIMITS = list(
            CHROMOSOME = list(
                MIN = 1,
                MAX = 16
            )
        )
    ),
    
    SCRIPTS = list(
        "setup_experiment" = list(
            args = list(
                directory = list(
                    type = "character",
                    required = TRUE,
                    description = "Experiment directory name",
                    validation = function(value) {
                        dir.exists(file.path(CONFIG$VALIDATION$PATHS$DATA_DIR, value))
                    },
                    error_message = "Directory does not exist in data directory"
                )
            )
        ),
        
        "update_sample_grid" = list(
            args = list(
                directory = list(
                    type = "character",
                    required = TRUE,
                    description = "Experiment directory name",
                    validation = function(value) {
                        dir.exists(file.path(CONFIG$VALIDATION$PATHS$DATA_DIR, value))
                    },
                    error_message = "Directory does not exist in data directory"
                )
            )
        ),
        
        "plot_sample_tracks" = list(
            args = list(
                directory = list(
                    type = "character",
                    required = TRUE,
                    description = "Experiment directory name",
                    validation = function(value) {
                        dir.exists(file.path(CONFIG$VALIDATION$PATHS$DATA_DIR, value))
                    },
                    error_message = "Directory does not exist in data directory"
                ),
                chromosome = list(
                    type = "integer",
                    required = FALSE,
                    default = 10,
                    description = "Chromosome number to plot",
                    validation = function(value) {
                        value >= CONFIG$VALIDATION$LIMITS$CHROMOSOME$MIN && 
                        value <= CONFIG$VALIDATION$LIMITS$CHROMOSOME$MAX
                    },
                    error_message = "Chromosome must be between 1 and 16"
                )
            )
        )
    )
)
