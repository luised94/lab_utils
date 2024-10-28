#!/usr/bin/env Rscript
# R/config/experiment_config.R

#' Experiment Configuration
#' @description Define experimental setup and validation
#' @export EXPERIMENT_CONFIG

EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "240808Bel",
        EXPECTED_SAMPLES = 33,
        VERSION = "1.0.0"
    ),
    
    CATEGORIES = list(
        strain_source = c("lemr", "oa"),
        rescue_allele = c("none", "wt"),
        mcm_tag = c("none", "2", "7"),
        auxin_treatment = c("no", "yes"),
        cell_cycle = c("G1", "M"),
        antibody = c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")
    ),
    
    INVALID_COMBINATIONS = list(
        lemr_mcm7 = quote(strain_source == "lemr" & mcm_tag == "7"),
        lemr_mcm2 = quote(strain_source == "lemr" & mcm_tag == "2"),
        oa_mcm2 = quote(strain_source == "oa" & mcm_tag == "2"),
        oa_wt = quote(strain_source == "oa" & rescue_allele == "wt")
    ),
    
    EXPERIMENTAL_CONDITIONS = list(
        is_input = quote(rescue_allele == "none" & cell_cycle == "M" & 
                        antibody == "Input" & 
                        ((strain_source == "oa" & auxin_treatment == "no") | 
                         (strain_source == "lemr" & auxin_treatment == "no"))),
        is_protg = quote(rescue_allele == "wt" & mcm_tag == "none" & 
                        cell_cycle == "M" & antibody == "ProtG" & 
                        strain_source == "lemr" & auxin_treatment == "no")
        # [Additional conditions...]
    ),
    
    COMPARISONS = list(
        comp_cellCycle = quote(strain_source == "lemr" & antibody == "74"),
        comp_alfa = quote(antibody == "ALFA")
        # [Additional comparisons...]
    ),
    
    CONTROL_FACTORS = list(
        genotype = c("strain_source", "rescue_allele", "mcm_tag")
    ),
    
    COLUMN_ORDER = c(
        "antibody", "strain_source", "rescue_allele", 
        "mcm_tag", "auxin_treatment", "cell_cycle"
    )
)

#' Generate Experiment Sample Grid
#' @param project_config List Project configuration
#' @param experiment_config List Experiment configuration
#' @param init_logging Logical Whether to initialize logging
#' @return data.frame Filtered experiment grid
generate_experiment_grid <- function(
    project_config = NULL,
    experiment_config = EXPERIMENT_CONFIG,
    init_logging = TRUE
) {
    # Handle project configuration
    if (is.null(project_config)) {
        if (!exists("PROJECT_CONFIG", envir = .GlobalEnv)) {
            warning("PROJECT_CONFIG not found. Some features may be limited.")
        } else {
            project_config <- PROJECT_CONFIG
        }
    }
    
    # Initialize logging if requested
    if (init_logging) {
        if (file.exists("~/lab_utils/R/functions/logging_utils.R")) {
            source("~/lab_utils/R/functions/logging_utils.R")
            log_file <- initialize_logging(script_name = "experiment_grid")
            log_info("Starting experiment grid generation", log_file)
        } else {
            warning("logging_utils.R not found. Proceeding without logging.")
        }
    }
    
    tryCatch({
        # Generate combinations
        grid <- do.call(expand.grid, experiment_config$CATEGORIES)
        
        # Filter invalid combinations
        invalid_idx <- Reduce(
            `|`, 
            lapply(experiment_config$INVALID_COMBINATIONS, eval, envir = grid)
        )
        grid <- subset(grid, !invalid_idx)
        
        # Apply experimental conditions
        valid_idx <- Reduce(
            `|`, 
            lapply(experiment_config$EXPERIMENTAL_CONDITIONS, eval, envir = grid)
        )
        grid <- subset(grid, valid_idx)
        
        # Verify sample count
        n_samples <- nrow(grid)
        if (n_samples != experiment_config$METADATA$EXPECTED_SAMPLES) {
            warning(sprintf(
                "Expected %d samples, got %d", 
                experiment_config$METADATA$EXPECTED_SAMPLES, 
                n_samples
            ))
        }
        
        # Sort columns
        grid <- grid[, experiment_config$COLUMN_ORDER]
        
        # Add attributes
        attr(grid, "control_factors") <- experiment_config$CONTROL_FACTORS
        attr(grid, "experiment_id") <- experiment_config$METADATA$EXPERIMENT_ID
        
        if (init_logging) {
            log_info(sprintf("Generated grid with %d samples", nrow(grid)))
        }
        
        return(grid)
        
    }, error = function(e) {
        msg <- sprintf("Failed to generate experiment grid: %s", e$message)
        if (init_logging) {
            log_error(msg)
        }
        stop(msg)
    })
}
