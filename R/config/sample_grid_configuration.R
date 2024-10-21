#source("~/lab_utils/R/functions/002_table_operations.R")
#@update
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
configuration_settings$combinations_grid <- do.call(base::expand.grid,
    configuration_settings$all_categories_list)

# Define impossible combinations 
configuration_settings$impossible_settings <- with(configuration_settings$combinations_grid,
    ( strain_source == "lemr" & mcm_tag == "7" ) |
    ( strain_source == "lemr" & mcm_tag == "2" ) |
    ( strain_source == "oa" & mcm_tag == "2" ) |
    ( strain_source == "oa" & rescue_allele == "wt" )
)

configuration_settings$combinations_grid <- configuration_settings$combinations_grid[!configuration_settings$impossible_settings, ]
#@update

experiment_conditions <- list(
     is_input = with(configuration_settings$combinations_grid,
         rescue_allele == "none" &
         cell_cycle == "M" &
         antibody == "Input" &
         ((strain_source == "oa" & auxin_treatment == "no") |
         (strain_source == "lemr" & auxin_treatment == "no"))
     ),
     is_protg = with(configuration_settings$combinations_grid,
             rescue_allele == "wt" &
             mcm_tag == "none" &
             cell_cycle == "M" &
             antibody == "ProtG" &
             strain_source == "lemr" &
             auxin_treatment == "no"
     ),
     is_alfa = with(configuration_settings$combinations_grid,
         rescue_allele == "none" &
         mcm_tag == "none" &
         cell_cycle == "M" &
         antibody == "ALFA" &
         (( strain_source == "oa" & auxin_treatment == "no") | ( strain_source == "lemr"))
     ),
     is_1108 =  with(configuration_settings$combinations_grid,
             rescue_allele == "none" &
             mcm_tag == "none" &
             cell_cycle == "M" &
             antibody == "HM1108" &
            (( strain_source == "oa" &  auxin_treatment == "no") | strain_source == "lemr")
     ),
     is_174 = with(configuration_settings$combinations_grid,
         antibody == "74" &
         auxin_treatment == "no"
     ),
     is_cha = with(configuration_settings$combinations_grid,
         antibody == "CHA" &
         auxin_treatment == "no"
      ),
    is_11HA = with(configuration_settings$combinations_grid,
        antibody == "11HA" &
         auxin_treatment == "no" &
         !(strain_source == "lemr" & mcm_tag == "none" & rescue_allele == "wt" & cell_cycle == "M")
    )
)
configuration_settings$combinations_grid <- configuration_settings$combinations_grid[Reduce('|', experiment_conditions), ] 

configuration_settings$control_factors <- list(
    genotype = c("strain_source", "rescue_allele", "mcm_tag")
)


print(configuration_settings$combinations_grid)
