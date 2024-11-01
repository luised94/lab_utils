#' Experiment Setup Functions
setup_experiment <- function(experiment_name,
                           user_info) {
    log_info("Setting up experiment:", experiment_name)
    
    # Setup directories
    dirs <- setup_experiment_directories(
        experiment_name,
        user_info
    )
    
    # Setup configuration
    setup_experiment_config(
        experiment_name,
        dirs
    )
    
    dirs
}

setup_experiment_directories <- function(experiment_name,
                                      user_info) {
    log_info("Setting up experiment directories")
    
    # Create base directory
    base_dir <- file.path(
        sprintf(CONFIG$EXPERIMENT$PATHS$DROPBOX, 
                user_info$username),
        experiment_name
    )
    
    # Create subdirectories
    create_directory_structure(
        base_dir,
        CONFIG$EXPERIMENT$DIRECTORIES
    )
    
    list(
        base = base_dir,
        config = file.path(base_dir, "documentation")
    )
}

setup_experiment_config <- function(experiment_name,
                                  directories) {
    log_info("Setting up experiment configuration")
    
    # Load config template
    config_template <- load_config_template(experiment_name)
    
    # Validate configuration
    validate_experiment_config(
        config_template,
        experiment_name
    )
    
    # Output configuration
    if (CONFIG$EXPERIMENT$OUTPUT$ENABLED) {
        output_experiment_config(
            config_template,
            directories$config,
            experiment_name
        )
    }
}

validate_experiment_config <- function(config,
                                     experiment_name) {
    if (config$current_experiment != experiment_name) {
        log_error("Configuration mismatch:")
        log_error("Config experiment:", config$current_experiment)
        log_error("Provided experiment:", experiment_name)
        stop("Invalid configuration")
    }
    
    log_info("Configuration validated")
}
