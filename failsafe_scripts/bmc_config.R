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

# Validate configuration structure
stopifnot(
    "Missing required config sections" = 
        all(c("METADATA", "CATEGORIES", "INVALID_COMBINATIONS", 
              "EXPERIMENTAL_CONDITIONS", "COMPARISONS", "CONTROL_FACTORS", 
              "COLUMN_ORDER", "NORMALIZATION") %in% names(EXPERIMENT_CONFIG))
)

# Validate each section
validate_category_values(EXPERIMENT_CONFIG$CATEGORIES)

validate_column_references(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    comparisons = EXPERIMENT_CONFIG$COMPARISONS,
    control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS,
    conditions = EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS
)

validate_column_order(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER
)

cat("[VALIDATED] Experiment configuration loaded successfully\n")
