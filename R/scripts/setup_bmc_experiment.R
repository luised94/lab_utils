#!/usr/bin/env Rscript
# R/scripts/setup_bmc_experiment.R

#' Setup BMC Experiment
#' @param experiment_id Character Experiment identifier
#' @param config_path Character Path to project configuration
#' @return List Experiment setup results
setup_bmc_experiment <- function(
    experiment_id,
    config_path = "~/lab_utils/R/config/project_config.R"
) {
    tryCatch({
        # Initialize project
        source(config_path)
        source("~/lab_utils/R/config/bmc_sample_grid_config.R")
        
        # Setup logging
        log_file <- initialize_logging(
            script_name = sprintf("setup_experiment_%s", experiment_id)
        )
        
        # Generate and process sample grid
        sample_grid <- generate_experiment_grid()
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
