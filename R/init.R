# R/init.R

# Load all utility functions
source("R/functions/validation_utils.R")
# source other utility files here

# Function to load all scripts in a directory
load_scripts <- function(dir_path) {
    scripts <- list.files(dir_path, pattern = "\\.R$", full.names = TRUE)
    for (script in scripts) {
        source(script)
    }
}

# Load all function scripts
load_scripts("R/functions")

## Load all analysis scripts
#load_scripts("R/analysis")

# Any other initialization code can go here

# Log that initialization is complete
log_info("Initialization complete. All utility functions and scripts loaded.")
#!/usr/bin/env Rscript

source("config/package_config.R")
source("functions/package_manager.R")
source("functions/environment_validator.R")

#' Main initialization function
main <- function() {
    log_info("Starting environment setup")
    
    tryCatch({
        # Initialize environment
        initialize_environment()
        
        # Validate setup
        validate_environment()
        
        log_info("Environment setup completed successfully")
    }, error = function(e) {
        log_error("Environment setup failed:", e$message)
        quit(status = 1)
    })
}

if (!interactive()) {
    main()
}
