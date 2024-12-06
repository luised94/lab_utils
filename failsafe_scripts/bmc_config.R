################################################################################
# BMC Experiment Configuration
################################################################################
#
# PURPOSE:
#   Defines and validates experimental design for BMC ChIP-seq experiments,
#   including sample categories, valid combinations, and comparison groups.
#
# USAGE:
#   1. Use '/!!' in vim/neovim to jump to required updates
#   2. Modify METADATA section with experiment details
#   3. Update CATEGORIES if experimental design changes
#   4. Review INVALID_COMBINATIONS and EXPERIMENTAL_CONDITIONS
#
# !! ----> REQUIRED UPDATES:
# !! EXPERIMENT_CONFIG$METADATA <- list(
# !!     EXPERIMENT_ID = "241010Bel",
# !!     EXPECTED_SAMPLES = 65,
# !!     VERSION = "1.0.0"
# !! )
#
# STRUCTURE:
#   EXPERIMENT_CONFIG/
#   +-- METADATA/
#   |   +-- EXPERIMENT_ID    # Format: YYMMDD'Bel'
#   |   +-- EXPECTED_SAMPLES # Total valid combinations
#   |   +-- VERSION         # Configuration version
#   +-- CATEGORIES/         # Valid values for each factor
#   +-- INVALID_COMBINATIONS/# Excluded experimental combinations
#   +-- EXPERIMENTAL_CONDITIONS/# Valid sample definitions
#   +-- COMPARISONS/        # Analysis groupings
#   +-- CONTROL_FACTORS/    # Control sample definitions
#   +-- COLUMN_ORDER/       # Standard column arrangement
#
# VALIDATION:
#   1. Category Values: Must be character vectors, unique
#   2. Column References: All referenced columns must exist
#   3. Column Order: Must include all category columns
#   4. Sample Count: Must match EXPECTED_SAMPLES
#
# DEPENDENCIES:
#   - R base packages only
#   - functions_for_bmc_config_validation.R for validation functions
#
# COMMON ISSUES:
#   1. Mismatched EXPERIMENT_ID -> Check format YYMMDD'Bel'
#   2. Wrong sample count -> Review INVALID_COMBINATIONS
#   3. Missing categories -> Check CATEGORIES vs COLUMN_ORDER
#
# AUTHOR: Luis
# DATE: 2024-11-27
# VERSION: 2.0.0
#
################################################################################
# !! Update EXPERIMENT_CONFIG if starting a new experiment.
EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "241122Bel",
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

    SAMPLE_CLASSIFICATIONS = list(
        is_input = quote(antibody == "Input"),

        is_negative = quote(
            antibody == "ProtG" |  # Protein G negative control
            (antibody == "V5" & rescue_allele == "NONE") | # No-tag control
            (time_after_release == "0" & antibody == "UM174") | # MCM at G2
            (auxin_treatment == "YES" & antibody == "ALFA") # Degradation of Orc4-ALFA.
        ),

        is_positive = quote(
            (antibody == "HM1108") |  # Known working condition
            (antibody == "V5" & rescue_allele == "WT") | # Test by removing rescue_allele condition.
            (antibody == "UM174" & time_after_release %in% c("1", "2") & rescue_allele == "WT")
        )
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

    COLUMN_ORDER = c("antibody", "rescue_allele", "auxin_treatment", "time_after_release"),

    NORMALIZATION = list(
        methods = c("CPM", "BPM", "RPGC", "RPKM"),
        active = "CPM"  # Set via config
    )
)


################################################################################
# Configuration Validation
################################################################################
source("~/lab_utils/failsafe_scripts/functions_for_bmc_config_validation.R")

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
    conditions = EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS,
    verbose = validation_verbose
)

validate_column_order(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
    verbose = validation_verbose
)

cat("\n[VALIDATED] Experiment configuration loaded successfully\n")
