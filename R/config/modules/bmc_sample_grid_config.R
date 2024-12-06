#!/usr/bin/env Rscript
# TODO: Do we need to validate the experiment_id format?
# TODO: Should we add validation for category value uniqueness?
# R/config/experiment_config.R

#' Experiment Configuration
#' @description Define experimental setup and validation
#' @export EXPERIMENT_CONFIG


EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "241010Bel",
        EXPECTED_SAMPLES = 65,
        VERSION = "1.0.0"
    ),
    
    CATEGORIES = list(
        rescue_allele = c("NONE", "WT", "4R", "PS"),
        auxin_treatment = c("NO", "YES"),
        time_after_release = c("0", "1", "2"),
        antibody = c("Input", "ProtG", "HM1108", "V5", "ALFA", "UM174")
    ),
    
    INVALID_COMBINATIONS = list(
        rescue_allele_auxin_treatment = quote(rescue_allele %in% c("4R", "PS") & auxin_treatment == "NO"),
        protg_time_after_release = quote(antibody == "ProtG" & time_after_release %in% c("1", "2")),
        input_time_after_release = quote(antibody == "Input" & time_after_release %in% c("1", "2")),
        input_rescue_allele_auxin_treatment = quote(antibody == "Input" & rescue_allele %in% c("NONE", "WT") & auxin_treatment == "YES")
    ),
    
    EXPERIMENTAL_CONDITIONS = list(
        is_input = quote(time_after_release == "0" & antibody == "Input"),
        is_protg = quote(rescue_allele == "WT" & time_after_release == "0" & antibody == "ProtG" & auxin_treatment == "NO"),
        is_v5 = quote(antibody == "V5"),
        is_alfa = quote(antibody == "ALFA"),
        is_1108 = quote(antibody == "HM1108" & time_after_release == "0"),
        is_174 = quote(antibody == "UM174")
    ),
    
    COMPARISONS = list(
        comp_1108forNoneAndWT = quote(antibody == "HM1108" & rescue_allele %in% c("NONE", "WT")),
        comp_1108forNoneAndWT_auxin = quote(antibody == "HM1108" & auxin_treatment == "YES"),
        comp_timeAfterReleaseV5WT = quote(antibody == "V5" & rescue_allele == "WT" & auxin_treatment == "YES"),
        comp_timeAfterReleaseV5NoTag = quote(antibody == "V5" & rescue_allele == "NONE" & auxin_treatment == "YES"),
        comp_V5atTwoHours = quote(antibody == "V5" & time_after_release == "2" & auxin_treatment == "YES"),
        comp_UM174atTwoHours = quote(antibody == "UM174" & time_after_release == "2" & auxin_treatment == "YES"),
        comp_ALFAforNoRescueNoTreat = quote(antibody == "ALFA" & rescue_allele == "NONE" & auxin_treatment == "NO"),
        comp_ALFAforNoRescueWithTreat = quote(antibody == "ALFA" & rescue_allele == "NONE" & auxin_treatment == "YES"),
        comp_ALFAatTwoHoursForAllAlleles = quote(antibody == "ALFA" & time_after_release == "2" & auxin_treatment == "YES"),
        comp_UM174atZeroHoursForAllAlleles = quote(antibody == "UM174" & time_after_release == "0" & auxin_treatment == "YES"),
        comp_AuxinEffectOnUM174 = quote(antibody == "UM174" & time_after_release == "2" & rescue_allele %in% c("NONE", "WT"))
    ),
    
    CONTROL_FACTORS = list(
        genotype = c("rescue_allele")
    ),
    
    COLUMN_ORDER = c("antibody", "rescue_allele", "auxin_treatment", "time_after_release")
)

#' Validate Category Values
validate_experiment_categories <- function(experiment_config) {
    for (category_name in names(experiment_config$CATEGORIES)) {
        values <- experiment_config$CATEGORIES[[category_name]]
        if (!is.character(values)) {
            stop(sprintf(
                "Category '%s' values must be character vectors, got %s",
                category_name, class(values)
            ))
        }
    }
    invisible(TRUE)
}

#' Validate Column References
validate_experiment_column_references <- function(experiment_config) {
    # Get valid column names
    valid_columns <- names(experiment_config$CATEGORIES)
    
    # Check comparison groups
    for (comp_name in names(experiment_config$COMPARISONS)) {
        comp_expr <- experiment_config$COMPARISONS[[comp_name]]
        comp_vars <- all.vars(comp_expr)
        invalid_cols <- setdiff(comp_vars, valid_columns)
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in comparison '%s': %s",
                comp_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    # Check control factors
    for (factor_name in names(experiment_config$CONTROL_FACTORS)) {
        invalid_cols <- setdiff(
            experiment_config$CONTROL_FACTORS[[factor_name]],
            valid_columns
        )
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in control factor '%s': %s",
                factor_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    # Check experimental conditions
    for (cond_name in names(experiment_config$EXPERIMENTAL_CONDITIONS)) {
        cond_expr <- experiment_config$EXPERIMENTAL_CONDITIONS[[cond_name]]
        cond_vars <- all.vars(cond_expr)
        invalid_cols <- setdiff(cond_vars, valid_columns)
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in condition '%s': %s",
                cond_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    invisible(TRUE)
}

#' Validate Column Order
validate_experiment_column_order <- function(experiment_config) {
    if (!identical(
        sort(names(experiment_config$CATEGORIES)),
        sort(experiment_config$COLUMN_ORDER)
    )) {
        stop("Column order must include all category columns")
    }
    invisible(TRUE)
}

#' Generate Experiment Sample Grid
#' @param project_config List Project configuration
#' @param experiment_config List Experiment configuration
#' @param init_logging Logical Whether to initialize logging
#' @return data.frame Filtered experiment grid
generate_experiment_grid <- function(
    experiment_config = EXPERIMENT_CONFIG,
    init_logging = TRUE
) {
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
        # Run validations
        validate_experiment_categories(experiment_config)
        validate_experiment_column_references(experiment_config)
        validate_experiment_column_order(experiment_config)

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
