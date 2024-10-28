#!/usr/bin/env Rscript

source("../functions/environment_validator.R")
source("../functions/experiment_setup.R")

#' Main experiment setup function
main <- function() {
    log_info("Starting experiment setup")
    
    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) != 1) {
        log_error("Invalid number of arguments")
        show_usage()
        stop("Invalid input")
    }
    
    # Validate environment
    user_info <- validate_environment()
    
    # Setup experiment
    setup_experiment(args[1], user_info)
    
    # Output sync command
    output_sync_command()
    
    log_info("Experiment setup completed")
}

show_usage <- function() {
    cat("Usage: Rscript setup_experiment.R <experiment_name>\n")
    cat("Example: Rscript setup_experiment.R 240808Bel\n")
    cat("\nNote: Ensure WINDOWS_USER is defined in bashrc\n")
}

output_sync_command <- function() {
    log_info("To sync files:")
    log_info("Use: scp command (see documentation)")
}

if (!interactive()) {
    main()
}
