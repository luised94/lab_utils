#!/usr/bin/env Rscript

#' Project Initialization
initialize_project <- function() {
    # Setup error handling
    options(error = function() {
        cat("Error occurred in initialization\n")
        if (interactive()) recover()
    })
    
    tryCatch({
        # Load basic configuration
        source("config/init_config.R")
        
        # Check environment
        check_environment()
        
        # Load priority scripts first
        load_priority_scripts()
        
        # Load remaining function scripts
        load_directory_scripts(CONFIG$PATHS$FUNCTIONS)
        
        # Load script directories
        for (dir in c(CONFIG$PATHS$SCRIPTS, CONFIG$PATHS$CONFIG)) {
            load_directory_scripts(dir)
        }
        
        log_info("Project initialization completed successfully")
        
    }, error = function(e) {
        cat("Critical initialization error:", e$message, "\n")
        if (interactive()) {
            recover()
        } else {
            quit(status = 1)
        }
    })
}

initialize_project()
