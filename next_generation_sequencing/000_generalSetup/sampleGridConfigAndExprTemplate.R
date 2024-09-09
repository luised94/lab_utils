# Description: Configures and generates a sample table with experiments.
# Usage: Rscript sampleGridConfigAndExprTemplate.R <experiment_id>

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    current_experiment <- validate_input(args)
    print(current_experiment)
    all_categories_list <- define_categories()
    ordered_samples <- generate_filtered_samples(all_categories_list)
    named_samples <- add_sample_names_to_table(ordered_samples)
    named_samples$experiment_id <- current_experiment
    bmc_table <- create_bmc_table(named_samples)
    print_summary(named_samples, bmc_table)

}

validate_input <- function(args) {
    if(length(args) != 1){
        cat("Error: Script requires at least one argument.\n")
        cat("Usage: Rscript sampleGridConfigAndExprTemplate.R <experiment_id>\n")
        stop()
    }
    return(args[1])
}

# Define categories
define_categories <- function() {
    cat("Defining categories...\n")
    all_categories_list <- list(
        strain_source = c("lemr", "oa"),
        rescue_allele = c("none", "wt"),
        mcm_tag = c("none", "2", "7"),
        auxin_treatment = c("no", "yes"),
        cell_cycle = c("G1", "M"),
        antibody = c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")
    )
    return(all_categories_list)
}

filter_samples <- function(combinations_grid){
    is_input <- with(combinations_grid,
        rescue_allele == "none" &
        cell_cycle == "M" &
        antibody == "Input" &
        ((strain_source == "oa" & auxin_treatment == "no") | (strain_source == "lemr" & auxin_treatment == "no")) &
        !( strain_source == "lemr" & mcm_tag == "7" ) &
        !( strain_source == "lemr" & mcm_tag == "2" ) &
        !( strain_source == "oa" & mcm_tag == "2" ) &
        !( strain_source == "oa" & rescue_allele == "wt" )
    )
    is_protg <- with(combinations_grid,
            rescue_allele == "wt" &
            mcm_tag == "none" &
            cell_cycle == "M" &
            antibody == "ProtG" &
            strain_source == "oa" &
            auxin_treatment == "no"
    )
    is_alfa <- with(combinations_grid,
        rescue_allele == "none" &
        mcm_tag == "none" &
        cell_cycle == "M" &
        antibody == "ALFA" &
        (( strain_source == "oa" & auxin_treatment == "no") | ( strain_source == "lemr"))
    )
    is_1108 <-  with(combinations_grid,
            rescue_allele == "none" &
            mcm_tag == "none" &
            cell_cycle == "M" &
            antibody == "HM1108" &
           (( strain_source == "oa" &  auxin_treatment == "no") | strain_source == "lemr")
    )
    is_174 <- with(combinations_grid,
        antibody == "74" &
        auxin_treatment == "no" &
        !( strain_source == "lemr" &  rescue_allele == "none") &
        !( strain_source == "oa" &  rescue_allele == "wt") &
        !( strain_source == "lemr" &  mcm_tag == "7") &
        !( strain_source == "oa" &  mcm_tag == "2")
    )
    is_cha <- with(combinations_grid,
        antibody == "CHA" &
        auxin_treatment == "no" &
        !(strain_source == "lemr" & rescue_allele == "none") &
        !(strain_source == "oa" & rescue_allele == "wt") &
        !(strain_source == "lemr" & mcm_tag == "7") &
        !(strain_source == "oa" & mcm_tag == "2")
     )
    is_11HA <- with(combinations_grid,
        antibody == "11HA" &
        auxin_treatment == "no" &
        !(strain_source == "lemr" & rescue_allele == "none") &
        !(strain_source == "oa" & rescue_allele == "wt") &
        !(strain_source == "lemr" & mcm_tag == "7") &
        !(strain_source == "oa" & mcm_tag == "2") &
        !(strain_source == "lemr" & mcm_tag == "none" & rescue_allele == "wt" & cell_cycle == "M")
    )
    return(combinations_grid[is_input | is_protg | is_alfa | is_1108 | is_174 | is_cha | is_11HA , ])
}

generate_filtered_samples <- function(categories_list) {
    cat("Generating and filtering samples\n")
    combinations <- expand.grid(categories_list)
    filtered_samples <- filter_samples(combinations)
    reorder_index <- with(filtered_samples, order(antibody, strain_source, rescue_allele, mcm_tag, auxin_treatment, cell_cycle))
    ordered_samples <- filtered_samples[reorder_index, ]
    return(ordered_samples)
}

add_sample_names_to_table <- function(ordered_samples_table) {
    cat("Adding names to sample table\n")
    ordered_samples_table$full_name <- apply(ordered_samples_table, 1, paste, collapse = "_")
    colnames_sans_fullname <- !grepl("full_name", colnames(ordered_samples_table))
    ordered_samples_table_sans_fullname <- ordered_samples_table[, colnames_sans_fullname]
    ordered_samples_table$short_name <- apply(ordered_samples_table_sans_fullname, 1, function(row) paste0(substr(row, 1, 1), collapse = ""))
    return(ordered_samples_table)
}

add_comparisons <- function(ordered_samples_table) {
    cat("Adding columns with comparison values\n")
    df <- ordered_samples_table
    comparisons <- list(

    comp_cellCycle = with(df, 
    )


    )

    for(comp_name in names(comparisons)) {
        df[[comp_name]] <- comparisons[[comp_name]]

    }
    return(df)
}

create_bmc_table <- function(named_samples_table) {
    cat("Making bmc_table from sample table\n")
    bmc_table <- data.frame(SampleName = named_samples_table$full_name,
       Vol..uL = 10,
       Conc = NA,
       Type = "ChIP",
       Genome = "Saccharomyces cerevisiae",
       Notes = "none",
       Pool = "A"
    )
    return(bmc_table)
}
print_summary <- function(sample_table, bmc_table) {
    cat("Printing results...\n")
    cat("Dimensions of sample_table:\n")
    print(dim(sample_table))
    cat("Dimensions of bmc_table:\n")
    print(dim(bmc_table))
    cat("Breakdown by antibody:\n")
    print(table(sample_table$antibody))
    cat("Elements of sample_table:\n")
    print(sample_table)
    cat("Ensure all elements are ordered according to sample submission.\n")
    cat("First elements of BMC table:\n")
    print(head(bmc_table))
}

if(!interactive()) {
    main()
} else {

    args <- "testBel"
    current_experiment <- validate_input(args)
    print(current_experiment)
    all_categories_list <- define_categories()
    ordered_samples <- generate_filtered_samples(all_categories_list)
    named_samples <- add_sample_names_to_table(ordered_samples)
    bmc_table <- create_bmc_table(named_samples)
    named_samples$experiment_id <- current_experiment
    print_summary(named_samples, bmc_table)
    #complete_table <- add_comparisons(named_samples)
    #print_summary(complete_table, bmc_table)
    cat("Loaded all functions and testing variables.\n")
}
