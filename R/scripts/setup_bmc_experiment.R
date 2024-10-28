#!/usr/bin/env Rscript
# R/scripts/setup_bmc_experiment.R

#' Setup BMC Sequencing Experiment
#' @description Initialize experiment directory structure and configuration
#' @export

#' Initialize Experiment
#' @param experiment_id Character Experiment identifier
#' @param config_path Character Path to project configuration
#' @return List Experiment directories and status
initialize_bmc_experiment <- function(
    experiment_id,
    config_path = "~/lab_utils/R/config/project_config.R"
) {
    tryCatch({
        # Load configurations
        source(config_path)
        source("~/lab_utils/R/config/bmc_sample_grid_config.R")
        
        # Initialize logging
        log_file <- initialize_logging(
            script_name = sprintf("setup_experiment_%s", experiment_id)
        )
        log_info("Initializing experiment setup", log_file)
        
        # Validate experiment ID format
        if (!grepl("^\\d{6}Bel$", experiment_id)) {
            stop("Invalid experiment ID format. Expected: YYMMDD'Bel'")
        }
        
        # Setup directory structure
        experiment_dirs <- create_experiment_directories(
            experiment_id = experiment_id,
            required_dirs = PROJECT_CONFIG$VALIDATION$REQUIRED_DIRS,
            base_path = "~/data"
        )
        
        # Generate and validate sample grid
        sample_grid <- generate_experiment_grid(
            experiment_config = EXPERIMENT_CONFIG
        )
        
        # Save configurations
        save_experiment_files(
            experiment_id = experiment_id,
            dirs = experiment_dirs,
            sample_grid = sample_grid
        )
        
        log_info("Experiment setup completed successfully", log_file)
        invisible(experiment_dirs)
        
    }, error = function(e) {
        log_error(sprintf("Experiment setup failed: %s", e$message))
        stop(e)
    })
}

#' Create Experiment Directories
#' @param experiment_id Character Experiment identifier
#' @param required_dirs Character vector Required directories
#' @param base_path Character Base path for experiment
#' @return List Directory paths
create_experiment_directories <- function(
    experiment_id,
    required_dirs,
    base_path
) {
    experiment_path <- file.path(base_path, experiment_id)
    
    # Create main directory
    if (!dir.exists(experiment_path)) {
        dir.create(experiment_path, recursive = TRUE)
    }
    
    # Create subdirectories
    for (dir_name in required_dirs) {
        dir_path <- file.path(experiment_path, dir_name)
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
    }
    
    list(
        root = experiment_path,
        documentation = file.path(experiment_path, "documentation"),
        logs = file.path(experiment_path, "logs")
    )
}

#' Save Experiment Files
#' @param experiment_id Character Experiment identifier
#' @param dirs List Directory paths
#' @param sample_grid data.frame Experiment sample grid
#' @return None
save_experiment_files <- function(
    experiment_id,
    dirs,
    sample_grid
) {
    # Save sample grid
    write.csv(
        sample_grid,
        file = file.path(dirs$documentation, "sample_grid.csv"),
        row.names = FALSE
    )
    
    # Copy configuration
    file.copy(
        from = "~/lab_utils/R/config/bmc_sample_grid_config.R",
        to = file.path(dirs$documentation, "experiment_config.R")
    )
}

# Execute if run as script
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
        stop("Usage: Rscript setup_bmc_experiment.R <experiment_id>")
    }
    initialize_bmc_experiment(args[1])
}
