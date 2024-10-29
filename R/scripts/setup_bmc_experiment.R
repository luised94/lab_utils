#!/usr/bin/env Rscript
# R/scripts/setup_bmc_experiment.R

#' Setup BMC Experiment
#' @param experiment_id Character Experiment identifier
#' @param config_path Character Path to project configuration
#' @return List Experiment setup results
setup_bmc_experiment <- function(
    experiment_id,
    bmc_sample_grid_config_path = "~/lab_utils/R/config/bmc_sample_grid_config.R",
    initialization_script_path = "~/lab_utils/R/project_init.R"
) {
    tryCatch({
        # Validate experiment ID format
        if (!grepl("^\\d{6}Bel$", experiment_id)) {
            stop("Invalid experiment ID format. Expected: YYMMDD'Bel'")
        }

        # Initialize project
        # Load project_config.R and logging_utils.R
        source(initialization_script_path)
        source(bmc_sample_grid_config_path)
        
        # Setup logging
        log_file <- initialize_logging(
            script_name = sprintf("setup_experiment_%s", experiment_id)
        )
        
        # Generate and process sample grid
        sample_grid <- generate_experiment_grid(
            experiment_config = EXPERIMENT_CONFIG,
            init_logging = TRUE
        )
        processed_grid <- process_sample_grid(
            grid = sample_grid,
            experiment_id = experiment_id,
            experiment_config = EXPERIMENT_CONFIG
        )
        
        # Create experiment directories
        dirs <- create_experiment_directories(
            experiment_id = experiment_id,
            required_dirs = PROJECT_CONFIG$VALIDATION$REQUIRED_DIRS,
            base_path = "~/data"
        )
        
        # Generate BMC submission table
        bmc_table <- generate_bmc_table(processed_grid)
        
        # Save outputs
        save_experiment_files(
            experiment_id = experiment_id,
            dirs = dirs,
            sample_grid = processed_grid,
            bmc_table = bmc_table
        )
        
        log_info("Experiment setup completed successfully", log_file)
        invisible(list(
            dirs = dirs,
            sample_grid = processed_grid,
            bmc_table = bmc_table
        ))
        
    }, error = function(e) {
        log_error(sprintf("Experiment setup failed: %s", e$message))
        stop(e)
    })
}

#' Process Sample Grid
#' @param grid data.frame Raw experiment grid
#' @param experiment_id Character Experiment identifier
#' @return data.frame Processed grid
process_sample_grid <- function(
    grid,
    experiment_id,
    experiment_config
) {
    # Add sample names
    grid$full_name <- apply(grid, 1, paste, collapse = "_")
    grid$short_name <- apply(grid[, experiment_config$COLUMN_ORDER], 1, 
        function(x) paste0(substr(x, 1, 1), collapse = ""))
    
    # Add control columns
    for (factor in names(experiment_config$CONTROL_FACTORS)) {
        col_name <- sprintf("X__cf_%s", factor)
        grid[[col_name]] <- paste(
            experiment_config$CONTROL_FACTORS[[factor]], 
            collapse = ","
        )
    }
    
    return(grid)
}

#' Generate BMC Table
#' @param sample_grid data.frame Processed sample grid
#' @return data.frame BMC submission table
generate_bmc_table <- function(sample_grid) {
    data.frame(
        SampleName = sample_grid$full_name,
        Vol_uL = 10,
        Conc = 0,
        Type = "ChIP",
        Genome = "Saccharomyces cerevisiae",
        Notes = ifelse(
            sample_grid$antibody == "Input", 
            "Run on fragment analyzer.", 
            "Run on femto pulse."
        ),
        Pool = "A",
        stringsAsFactors = FALSE
    )
}

#' Save Experiment Files
#' @param experiment_id Character Experiment identifier
#' @param dirs List Directory paths
#' @param sample_grid data.frame Processed sample grid
#' @param bmc_table data.frame BMC submission table
#' @return None
save_experiment_files <- function(
    experiment_id,
    dirs,
    sample_grid,
    bmc_table
) {
    # Save sample grid
    write.csv(
        sample_grid,
        file = file.path(dirs$documentation, paste0(experiment_id,"_", "sample_grid.csv")),
        row.names = FALSE
    )
    
    # Save BMC table
    write.table(
        bmc_table,
        file = file.path(dirs$documentation, paste0(experiment_id,"_", "bmc_table.tsv")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
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

# Execute if run as script
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
        stop("Usage: Rscript setup_bmc_experiment.R <experiment_id>")
    }
    initialize_bmc_experiment(args[1])
}
