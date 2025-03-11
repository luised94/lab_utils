
EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "100303Bel",
        EXPECTED_SAMPLES = 5,
        VERSION = "1.0.0",
        PROJECT_ID = "project_001"
    ),

    CATEGORIES = list(
        orc_phenotype = c("WT", "orc1-161"),
        cell_cycle = c("alpha", "nocodazole", "async"),
        temperature = c("23", "37"),
        antibody = c("ORC", "Nucleosomes")
    ),

    INVALID_COMBINATIONS = list(

    ),

    EXPERIMENTAL_CONDITIONS = list(
        is_orc = quote(orc_phenotype == "WT" & cell_cycle == "nocodazole" & antibody == "ORC"),
        is_nucleosomes = quote(antibody == "Nucleosomes"),
        is_wt_37c = quote(orc_phenotype == "WT" & temperature == "37"),
        is_async = quote(orc_phenotype == "WT" & cell_cycle == "async" & temperature == "23")
    ),

    SAMPLE_CLASSIFICATIONS = list(
        is_positive = quote(orc_phenotype == "WT" & antibody == "ORC"),
        is_negative = quote(orc_phenotype == "orc1-161" & antibody == "ORC" & temperature == "37")
    ),

    COMPARISONS = list(
        wt_vs_mutant = quote(orc_phenotype == "WT" | orc_phenotype == "orc1-161"),
        temp_effect = quote(temperature == "23" | temperature == "37"),
        cell_cycle_effect = quote(cell_cycle %in% c("alpha", "nocodazole", "async"))
    ),

    CONTROL_FACTORS = list(
        genotype = c("orc_phenotype")
    ),

    COLUMN_ORDER = c("antibody", "orc_phenotype", "cell_cycle", "temperature"),

    NORMALIZATION = list(
        methods = c("CPM", "BPM", "RPGC", "RPKM"),
        active = "CPM"
    )
)

#########
# setup_bmc_experiment_flags: Control what is displayed during setup_bmc_experiment.R
# Distinct from debug configuration, this is relevant when setting up the experiment only.
#########
show_all_metadata <- TRUE
show_particular_metadata <- TRUE
category_to_show <- "antibody"
value_to_show <- "HM1108"
values_in_category <- unlist(EXPERIMENT_CONFIG$CATEGORIES[category_to_show])

stopifnot(
    "show_all_metadata has to be logical" = is.logical(show_all_metadata),
    "flag must be logical type" = is.logical(show_particular_metadata),
    "category must be in grid" = category_to_show %in% names(EXPERIMENT_CONFIG$CATEGORIES),
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
    # Execution Mode
    #show_debug_output = TRUE,      # (formerly debug_verbose)
    #require_confirmation = FALSE,   # (formerly debug_interactive)
    #validate_extensively = TRUE,    # (formerly debug_validate)

    ## Processing Scope
    #target_comparison = "comp_1108forNoneAndWT",
    #target_chromosome = 10,
    #target_batch = 10,
    #samples_per_batch = 4,

    ## Output Control
    #plot_display_duration = 2

)

################################################################################
# Configuration Validation
################################################################################
experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
stopifnot(
    "Experiment ID must be a character string" = is.character(experiment_id),
    "Invalid experiment ID format. Expected: YYMMDD'Bel'" = grepl("^\\d{6}Bel$", experiment_id)
)
source("~/lab_utils/core_scripts/functions_for_bmc_config_validation.R")

# !! Update if you want thourough messages during validation.
validation_verbose <- FALSE  # Set to TRUE for detailed validation output

# Validate configuration structure
if (validation_verbose) cat("\nValidating EXPERIMENT_CONFIG structure...\n")

required_sections <- c("METADATA", "CATEGORIES", "INVALID_COMBINATIONS",
                      "EXPERIMENTAL_CONDITIONS", "COMPARISONS",
                      "CONTROL_FACTORS", "COLUMN_ORDER", "NORMALIZATION",
                      "SAMPLE_CLASSIFICATIONS")


missing_sections <- setdiff(required_sections, names(EXPERIMENT_CONFIG))
if (length(missing_sections) > 0) {
    stop(sprintf("Missing required config sections: %s",
                paste(missing_sections, collapse = ", ")))
}

if (validation_verbose) cat("[PASS] All required sections present\n\n")

# Validate configuration structure
stopifnot(
    "Missing required config sections" =
        all(required_sections %in% names(EXPERIMENT_CONFIG))
)

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
