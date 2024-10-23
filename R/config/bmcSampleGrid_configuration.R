#source("~/lab_utils/R/functions/001_logging.R")
source("~/lab_utils/R/functions/002_table_operations.R")
# Configuration settings
configuration_settings <- list(
    current_experiment = "240808Bel",
    expected_number_of_samples = 33,
    all_categories_list = list(
        strain_source = c("lemr", "oa"),
        rescue_allele = c("none", "wt"),
        mcm_tag = c("none", "2", "7"),
        auxin_treatment = c("no", "yes"),
        cell_cycle = c("G1", "M"),
        antibody = c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")
    )
)

# Generate combinations grid
configuration_settings$combinations_grid <- do.call(expand.grid, configuration_settings$all_categories_list)

# Define impossible combinations
impossible_conditions <- list(
    lemr_mcm7 = quote(strain_source == "lemr" & mcm_tag == "7"),
    lemr_mcm2 = quote(strain_source == "lemr" & mcm_tag == "2"),
    oa_mcm2 = quote(strain_source == "oa" & mcm_tag == "2"),
    oa_wt = quote(strain_source == "oa" & rescue_allele == "wt")
)

# Determine the indexes that are not part of the experiment.
configuration_settings$impossible_settings <- Reduce(`|`, lapply(impossible_conditions, eval, envir = configuration_settings$combinations_grid))

# Filter out impossible combinations
configuration_settings$combinations_grid <- subset(configuration_settings$combinations_grid, subset = !configuration_settings$impossible_settings)

# Define experiment conditions
experiment_conditions <- list(
    is_input = quote(rescue_allele == "none" & cell_cycle == "M" & antibody == "Input" & 
                     ((strain_source == "oa" & auxin_treatment == "no") | (strain_source == "lemr" & auxin_treatment == "no"))),
    is_protg = quote(rescue_allele == "wt" & mcm_tag == "none" & cell_cycle == "M" & antibody == "ProtG" & 
                     strain_source == "lemr" & auxin_treatment == "no"),
    is_alfa = quote(rescue_allele == "none" & mcm_tag == "none" & cell_cycle == "M" & antibody == "ALFA" & 
                    ((strain_source == "oa" & auxin_treatment == "no") | (strain_source == "lemr"))),
    is_1108 = quote(rescue_allele == "none" & mcm_tag == "none" & cell_cycle == "M" & antibody == "HM1108" & 
                    ((strain_source == "oa" & auxin_treatment == "no") | strain_source == "lemr")),
    is_174 = quote(antibody == "74" & auxin_treatment == "no"),
    is_cha = quote(antibody == "CHA" & auxin_treatment == "no"),
    is_11HA = quote(antibody == "11HA" & auxin_treatment == "no" & 
                    !(strain_source == "lemr" & mcm_tag == "none" & rescue_allele == "wt" & cell_cycle == "M"))
)

# Apply experiment conditions
configuration_settings$combinations_grid <- subset(configuration_settings$combinations_grid, subset = Reduce(`|`, lapply(experiment_conditions, eval, envir = configuration_settings$combinations_grid)))

verify_expected_number_of_samples(configuration_settings$combinations_grid, configuration_settings$expected_number_of_samples)

# Define the order for sorting the columns. Must include all columns.
column_sort_order <- c("antibody", "strain_source","rescue_allele","mcm_tag","auxin_treatment","cell_cycle")

configuration_settings$combinations_grid <- sort_columns(configuration_settings$combinations_grid, column_sort_order)

configuration_settings$experimental_comparisons <- list(
        #Effect of MCM: Input, G1 and M phase for 174
        comp_cellCycle = quote(strain_source == "lemr" & antibody == "74"),

        #Effect of Alfa: Input, no auxin, yes auxin
        comp_alfa = quote(antibody == "ALFA"),

        #Effect of CHA in cell cycle: Input, G1 and M phase for CHA
        comp_HaCoa = quote(strain_source == "oa" & antibody == "CHA"),

        #Effect of CHA in cell cycle: Input, G1 and M phase for CHA
        comp_11HAoa = quote(strain_source == "oa" & antibody == "11HA")
)

configuration_settings$combinations_grid <- add_comparisons(configuration_settings$combinations_grid,
    configurations_settings$experimental_comparisons)


# Define control factors
configuration_settings$control_factors <- list(
    genotype = c("strain_source", "rescue_allele", "mcm_tag")
)

configuration_settings$combinations_grid <- add_attributes(configuration_settings$combinations_grid,
    configuration_settings$control_factors)

# Print results
print(configuration_settings$combinations_grid)
