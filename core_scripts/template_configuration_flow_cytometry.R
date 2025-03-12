
EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "Exp_20250310_1",
        EXPECTED_SAMPLES = 42,
        VERSION = "1.0.0",
        PROJECT_ID = "project_001"
    ),

    CATEGORIES = list(
        rescue_allele = c("NONE", "WT", "4R"),
        suppressor_allele = c("NONE", "4PS"),
        auxin_treatment = c("NO", "YES"),
        timepoints = c("0", "20", "40", "60", "80", "100", "120")
    ),

    INVALID_COMBINATIONS = list(
        auxin_treatment = quote(
            auxin_treatment == "NO" &
            (suppressor_allele == "4PS" | rescue_allele == "4R")
        )
    ),

    #EXPERIMENTAL_CONDITIONS = list(
    #    is_orc = quote(orc_phenotype == "WT" & cell_cycle == "nocodazole" & antibody == "ORC"),
    #    is_nucleosomes = quote(antibody == "Nucleosomes"),
    #    is_wt_37c = quote(orc_phenotype == "WT" & temperature == "37"),
    #    is_async = quote(orc_phenotype == "WT" & cell_cycle == "async" & temperature == "23")
    #),

    #SAMPLE_CLASSIFICATIONS = list(
    #    is_positive = quote(orc_phenotype == "WT" & antibody == "ORC"),
    #    is_negative = quote(orc_phenotype == "orc1-161" & antibody == "ORC" & temperature == "37")
    #),

    COLUMN_ORDER = c(
        "rescue_allele", "suppressor_allele", "auxin_treatment", "timepoints"
    )
)

#########
# Troubleshooting flags: Control what is displayed during setup_flow_cytometry_experiment.R
# Distinct from debug configuration. 
# This is relevant when setting up the experiment only.
#########
show_all_metadata <- TRUE
show_particular_metadata <- TRUE
category_to_show <- "auxin_treatment"
value_to_show <- "NO"
values_in_category <- unlist(EXPERIMENT_CONFIG$CATEGORIES[category_to_show])

stopifnot(
    "show_all_metadata has to be logical" = is.logical(show_all_metadata),
    "flag must be logical type" = is.logical(show_particular_metadata),
    "category to show must be in grid" = category_to_show %in% names(EXPERIMENT_CONFIG$CATEGORIES),
    "Value must be in category" = value_to_show %in% values_in_category
)

################################################################################
# Time Configurations
################################################################################
TIME_CONFIG <- list(
    # Format specifications
    timestamp_format = "%Y%m%d_%H%M%S",    # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d",                # YYYYMMDD

    # Current values
    current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
    current_date = format(Sys.Date(), "%Y%m%d")
)

################################################################################
# DEBUG CONFIGURATIONS
################################################################################
RUNTIME_CONFIG <- list(
    # Core control flags
    debug_interactive = FALSE,
    debug_verbose = TRUE,
    debug_validate = TRUE,
    # Processing control
    process_single_file = FALSE,
    process_file_index = 1,
    # Output control
    output_save_plots = FALSE,
    output_dry_run = TRUE
)

################################################################################
# Configuration Validation
################################################################################
experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
expected_format_for_experiment_id <- "Exp_\\d{8}_\\d{1,6}"
stopifnot(
    "Experiment ID must be a character string" = is.character(experiment_id),
    "Invalid experiment ID format. Expected: Exp_YYYMMDD_[0-9]{6}" = grepl(expected_format_for_experiment_id, experiment_id)
)
source("~/lab_utils/core_scripts/functions_for_bmc_config_validation.R")

# !! Update if you want thourough messages during validation.
validation_verbose <- FALSE

# Validate configuration structure
if (validation_verbose) cat("\nValidating EXPERIMENT_CONFIG structure...\n")

required_sections <- c("METADATA", "CATEGORIES", "INVALID_COMBINATIONS",
                       "COLUMN_ORDER")

missing_sections <- setdiff(required_sections, names(EXPERIMENT_CONFIG))
if (length(missing_sections) > 0) {
    stop(sprintf("Missing required config sections: %s",
                paste(missing_sections, collapse = ", ")))
}

if (validation_verbose) cat("[PASS] All required sections present\n\n")

# Validate each section
validate_category_values(
    EXPERIMENT_CONFIG$CATEGORIES,
    verbose = validation_verbose
)

validate_column_references(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    comparisons = EXPERIMENT_CONFIG$COMPARISONS,
    control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS,
    #conditions = EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS,
    verbose = validation_verbose
)

validate_column_order(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
    verbose = validation_verbose
)

cat("\n[VALIDATED] Experiment configuration loaded successfully\n")
