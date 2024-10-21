
filter_samples <- function(combinations_grid, impossible_combinations, experimental_combinations){

    #Filter initial grid with impossible combinations based on experiment design.
    valid_combinations <- combinations_grid[!impossible_combinations, ]

    # Combine all experimental conditions and apply final filter
    experimental_conditions <- Reduce('|', experimental_combinations)
    return(valid_combinations[experimental_conditions, ])
}

generate_filtered_samples <- function(categories_list) {
    cat("Generating and filtering samples\n")
    combinations <- expand.grid(categories_list)
    filtered_samples <- filter_samples(combinations)
    #@update
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
    #@update
    cat("Adding columns with comparison values\n")
    df <- ordered_samples_table
    comparisons <- list(
        #Effect of MCM: Input, G1 and M phase for 174
        comp_cellCycle = with(df,
                strain_source == "lemr" & antibody == "74"
        ),

        #Effect of Alfa: Input, no auxin, yes auxin
        comp_alfa = with(df,
                antibody == "ALFA"
        ),

        #Effect of CHA in cell cycle: Input, G1 and M phase for CHA
        comp_HaCoa = with(df,
                strain_source == "oa" & antibody == "CHA"
        ),

        #Effect of CHA in cell cycle: Input, G1 and M phase for CHA
        comp_11HAoa = with(df,
                strain_source == "oa" & antibody == "11HA"
        )
    )

    # For all comparisons, create column with that name and TRUE/FALSE values for rows.
    for(comp_name in names(comparisons)) {
        df[[comp_name]] <- comparisons[[comp_name]]

    }
    return(df)
}
add_attributes <- function(table_with_comparisons, control_factors) {
    cat("Adding attributes to column names\n")
    df <- table_with_comparisons
    control_columns <- unlist(unname(control_factors))
    at_least_one_not_in_df_column <- !all(control_columns %in% colnames(df))
    if(at_least_one_not_in_df_column){
        control_column_not_in_df <- which(!(control_columns %in% colnames(df)))
        #cat("Columns not in the sample table\n")
        #print(control_columns[control_column_not_in_df])
        stop("Verify the columns in categories and control_factors list to ensure you are assigning correctly\n")
    }

    for (factor in names(control_factors)) {
        #cat(sprintf("Assign column to %s\n", factor))
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
print_summary <- function(sample_table, bmc_table) {
    cat("Printing results...\n")
    cat("Elements of sample_table:\n")
    print(sample_table)
    cat("First elements of BMC table:\n")
    print(head(bmc_table))
    cat("Dimensions of sample_table:\n")
    print(dim(sample_table))
    cat("Dimensions of bmc_table:\n")
    print(dim(bmc_table))
    cat("Breakdown by antibody:\n")
    print(table(sample_table$antibody))
}

process_control_factors <- function(sample_table) {
    cat("Process control factors from __cf_ columns\n")
    df <- sample_table
    cf_cols <- grep("^X__cf_", names(df), value = TRUE)
    if(length(cf_cols) == 0) {
        cat("No columns containing X__cf_ tag found in sample table")
        stop("Verify sample table was produced with updated sampleGridConfig.")
    }
    control_factors <- lapply(df[cf_cols], function(x) strsplit(x[1], ",")[[1]])
    names(control_factors) <- sub("^X__cf_", "", cf_cols)
    df[cf_cols] <- NULL
    attr(df, "control_factors") <- control_factors
    return(df)
}

get_factors_to_match <- function(sample_table) {
    cat("Grabbing attributes from sample table\n")
    df <- sample_table
    control_factors <- attr(df, "control_factors")
    if (is.null(control_factors)) {
        stop("No control factors defined in sample data.")
    }
    all_factors <- unlist(control_factors)
    return(intersect(all_factors, colnames(df)))
}
determine_matching_control <- function(sample_row, sample_table, factors_to_match) {
    cat("Determining control row for sample row.\n")
    df <- sample_table
    comparison_row <- sample_row[factors_to_match]
    rows_with_same_factors <- apply(df[, factors_to_match], 1, function(row) {
        all(row == comparison_row)
    })
    is_input <- df$antibody == "Input"
    index <- as.numeric(unname(which(is_input & rows_with_same_factors)))
    return(index)
}
select_control_index <- function(control_indices, max_controls = 1) {
    cat("Processing control index to ensure one is used.\n")
    if (length(control_indices) == 0) {
    stop("No matching control found")
    }
    if (length(control_indices) > max_controls) {
    warning(paste("Multiple matching controls found, using first", max_controls))
    control_indices[1:max_controls]
    } else {
    control_indices
    }
}
verify_number_of_samples <- function(ordered_samples, expected_number_of_samples) {
    if (nrow(ordered_samples) != expected_number_of_samples){
        cat("Breakdown by antibody\n")
        print(table(ordered_samples$antibody))
        cat("All samples\n")
        print(ordered_samples)
        cat(sprintf("Number of samples after filtering: %s\n", nrow(ordered_samples)))
        cat(sprintf("Number of samples expected: %s\n", expected_number_of_samples))
        stop("Number of rows filtered does not match expected number of samples.")
    }
}
