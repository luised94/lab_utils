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

#' Generate Experiment Grid
#' @param config List Experiment configuration
#' @return data.frame Filtered experiment grid
generate_experiment_grid <- function(
    project_config = PROJECT_CONFIG
    experimet_config = EXPERIMENT_CONFIG
) {
    # Check project configuration
    if (!exists("PROJECT_CONFIG")) {
        stop("Project configuration not loaded")
    }
    
    # Generate combinations
    grid <- do.call(expand.grid, config$CATEGORIES)
    
    # Filter invalid combinations
    invalid_idx <- Reduce(
        `|`, 
        lapply(config$INVALID_COMBINATIONS, eval, envir = grid)
    )
    grid <- subset(grid, !invalid_idx)
    
    # Apply experimental conditions
    valid_idx <- Reduce(
        `|`, 
        lapply(config$EXPERIMENTAL_CONDITIONS, eval, envir = grid)
    )
    grid <- subset(grid, valid_idx)
    
    # Verify sample count
    n_samples <- nrow(grid)
    if (n_samples != config$METADATA$EXPECTED_SAMPLES) {
        log_warning(sprintf(
            "Expected %d samples, got %d", 
            config$METADATA$EXPECTED_SAMPLES, 
            n_samples
        ))
    }
    
    # Sort columns
    grid <- grid[, config$COLUMN_ORDER]
    
    # Add attributes
    attr(grid, "control_factors") <- config$CONTROL_FACTORS
    attr(grid, "experiment_id") <- config$METADATA$EXPERIMENT_ID
    
    return(grid)
}
