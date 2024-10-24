# Description: Configures and generates a sample table with experiments.
# Usage: Rscript sampleGridConfigAndExprTemplate.R <experiment_id>
# functions moved
source("~/lab_utils/R/functions/001_logging.R")
source("~/lab_utils/R/functions/002_table_operations.R")
source("~/lab_utils/R/config/sample_grid_configuration.R")

validate_input <- function(args) {
    if(length(args) != 2){
        cat("Error: Script requires two arguments.\n")
        cat("Usage: Rscript sampleGridConfigAndExprTemplate.R <experiment_id> <expected_number_of_samples>\n")
        stop()
    } else if (is.numeric(args[2])) {
        cat("Error: Second argument must be a number larger than zero.\n")
        cat("Usage: Rscript sampleGridConfigAndExprTemplate.R <experiment_id> <expected_number_of_samples>\n")
        stop()
    }
    validated_args <- list(
        current_experiment = args[1],
        expected_number_of_samples = args[2]
    )
    return(validated_args)
}
args <- c(current_experiment, expected_number_of_samples)
main <- function(args) {
    validated_args <- validate_input(args)
    current_experiment <- validated_args$current_experiment
    all_categories_list <- define_categories()
    ordered_samples <- generate_filtered_samples(all_categories_list)
    verify_number_of_samples(ordered_samples, validated_args$expected_number_of_samples)
    named_samples <- add_sample_names_to_table(ordered_samples)
    named_samples$experiment_id <- current_experiment
    bmc_table <- create_bmc_table(named_samples)
    table_with_comparisons <- add_comparisons(named_samples)
    complete_table <- add_attributes(table_with_comparisons, control_factors)
    print_summary(complete_table, bmc_table)
    return(list(
        sample_table = complete_table,
        bmc_table = bmc_table
    ))
}
if(!interactive()) {
    sample_config_output <- main(args)
    cat("Script run through command line.\n")
    cat("Open in text editor and modify @update tag sections then run 000_setupExperimentDir.R.\n")
} else {
    validated_args <- validate_input(args)
    current_experiment <- validated_args$current_experiment
    print(current_experiment)
    all_categories_list <- define_categories()
    ordered_samples <- generate_filtered_samples(all_categories_list)
    verify_number_of_samples(ordered_samples, validated_args$expected_number_of_samples)
    named_samples <- add_sample_names_to_table(ordered_samples)
    bmc_table <- create_bmc_table(named_samples)
    named_samples$experiment_id <- current_experiment
    table_with_comparisons <- add_comparisons(named_samples)
    complete_table <- add_attributes(table_with_comparisons, control_factors)
    print_summary(complete_table, bmc_table)
    cat("Loaded all functions and testing variables.\n")
    cat("Script sourced from repl.\n")
    cat("Open in text editor and modify @update tag sections then run 000_setupExperimentDir.R.\n")
}
