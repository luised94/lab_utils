source("~/lab_utils/R/functions/001_logging.R")
source("~/lab_utils/R/init.R")
library(assertthat)

sort_columns <- function(df, column_sort_order) {
    # Input validation
    if(!setequal(column_sort_order, colnames(df))) {
        log_error("Column must be sorted using all columns to ensure proper order.")
        stop("Modify column_sort_order variable appropriately.")
    }
    # Create a list of sorting criteria
    sort_criteria <- lapply(column_sort_order, function(col) df[[col]])
    # Sort the dataframe
    df[do.call(order, sort_criteria), ]
}

add_sample_names_to_table <- function(df) {
     df$full_name <- apply(df, 1, paste, collapse = "_")
     colnames_sans_fullname <- !grepl("full_name", colnames(df))
     df_sans_fullname <- df[, colnames_sans_fullname]
     df$short_name <- apply(df_sans_fullname, 1, function(row) paste0(substr(row, 1, 1), collapse = ""))
     return(df)
 }

verify_expected_number_of_samples <- function(df, expected_number_of_samples) {
    # Input validation
    assert_that(is.numeric(expected_number_of_samples), "expected_number_of_samples must be a numeric.")
    if(!nrow(df) == expected_number_of_samples) {
        log_error("Combinations grid does not contain the expected_number_of_samples.")
        log_info("Elements of sample_table:\n")
        print(df)
        log_info("Dimensions of sample_table:\n")
        print(dim(df))
        log_info("Breakdown by antibody:\n")
        print(table(df$antibody))
        stop("Update the impossible_settings and experiment_conditions variable.\nEnsure that it matches configuration_settings$expected_number_of_samples.")
    }
}

add_comparisons <- function(df, comparison_list) {
    log_info("Adding columns with comparison values\n")

    # Input validation
    if (!is.data.frame(df)) {
        stop("Input 'df' must be a data.frame", call. = FALSE)
    }
    if (!is.list(comparison_list) || length(comparison_list) == 0) {
        stop("Input 'comparison_list' must be a non-empty list", call. = FALSE)
    }


    ## Check if all comparison expressions are valid
    #invalid_comps <- sapply(comparison_list, function(comp) !is.language(comp))
    #if (any(invalid_comps)) {
    #    stop("Invalid comparison expressions: ", 
    #         paste(names(comparison_list)[invalid_comps], collapse = ", "), 
    #         call. = FALSE)
    #}
    #                                                                                                       
    ## Get all unique column names referenced in comparisons
    #all_cols <- unique(unlist(lapply(comparison_list, all.vars)))
    #missing_cols <- setdiff(all_cols, names(df))
    #if (length(missing_cols) > 0) {
    #    stop("Columns referenced in comparisons but not in dataframe: ", 
    #         paste(missing_cols, collapse = ", "), 
    #         call. = FALSE)
    #}
    #                                                                                                       
    ## Performance optimization: pre-allocate result list
    #results <- vector("list", length(comparison_list))
    #names(results) <- names(comparison_list)
    #                                                                                                       
    ## Evaluate comparisons
    #for (comp_name in names(comparison_list)) {
    #    tryCatch({
    #        log_info(paste("Evaluating comparison:", comp_name))
    #        results[[comp_name]] <- eval(comparison_list[[comp_name]], df)
    #        if (!is.logical(results[[comp_name]])) {
    #            stop("Comparison result must be logical", call. = FALSE)
    #        }
    #        if (length(results[[comp_name]]) != nrow(df)) {
    #            stop("Comparison result length does not match number of rows in dataframe", call. = FALSE)
    #        }
    #    }, error = function(e) {
    #        stop(paste("Error in comparison", comp_name, ":", e$message), call. = FALSE)
    #    })
    #}
    #                                                                                                       
    ## Add results to dataframe
    #df[names(results)] <- results
    # For all comparisons, create column with that name and TRUE/FALSE values for rows.
    for(comp_name in names(comparison_list)) {
        df[[comp_name]] <- eval(comparison_list[[comp_name]], df)

    }
    return(df)
}

# Function to add new comparisons easily
add_new_comparison <- function(df, name, condition) {
    df[[name]] <- eval(condition, df)
    return(df)
}

add_attributes <- function(df, control_factors) {
    cat("Adding attributes to column names\n")
    control_columns <- unlist(unname(control_factors))
    at_least_one_not_in_df_column <- !all(control_columns %in% colnames(df))
    if(at_least_one_not_in_df_column){
        control_column_not_in_df <- which(!(control_columns %in% colnames(df)))
        log_error("One of the control factors is in the dataframe.")
        stop("Verify the columns in categories and control_factors list to ensure you are assigning correctly\n")
    }

    for (factor in names(control_factors)) {
        new_column_name <- paste0("X__cf_", factor)
        df[[new_column_name]] <- paste(control_factors[[factor]], collapse = ",")
    }
    return(df)

}

create_bmc_table <- function(named_samples_table) {
    cat("Making bmc_table from sample table\n")
    bmc_table <- data.frame(SampleName = named_samples_table$full_name,
       Vol..uL = 10,
       Conc = 0,
       Type = "ChIP",
       Genome = "Saccharomyces cerevisiae",
       Notes = ifelse(named_samples_table$antibody == "Input", "Run on fragment analyzer.", "Run on femto pulse."),
       Pool = "A"
    )
    return(bmc_table)
}

table_has_ID_column <- function(sample_table){
    if(!("sample_ID" %in% colnames(sample_table))){
        cat("No sample_ID column found.\n")
        cat("Must determine sample_IDs from fastq files\n")
        return(FALSE)
    } else {
        cat("Table has sample_ID column.\n")
        return(TRUE)
    }
}
modify_and_output_table <- function(sample_table, sample_ID_array, output_file_path) {
    if(nrow(sample_table) != length(sample_ID_array)) {
        cat("Number of rows is different from length of sample_ID_array.\n")
        cat("Verify fastq file names to ensure proper number is being extracted.\n")
        cat(sprintf("Number of rows: %s\n", nrow(sample_table)))
        cat(sprintf("Length of array: %s\n", length(sample_ID_array)))
        stop()
    } else if ("sample_ID" %in% colnames(sample_table)) {
        cat("sample_ID already part of the sample table.\n")
        print(colnames(sample_table))
        stop()
    } else {
        sample_table$sample_ID <- sample_ID_array
        print(head(sample_table))
        write.table(sample_table, output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        cat("Output modified sample table with sample ID column.\n")
    }
}
