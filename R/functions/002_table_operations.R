source("~/lab_utils/R/functions/001_logging.R")

sort_columns <- function(df, column_sort_order) {
    if(!setequal(column_sort_order, colnames(df))) {
        log_error("Column must be sorted using all columns to ensure proper error.")
        stop("Modify column_sort_order variable appropriately.")
    }
    # Create a list of sorting criteria
    sort_criteria <- lapply(column_sort_order, function(col) df[[col]])
    # Sort the dataframe
    df[do.call(order, sort_criteria), ]
}

add_sample_names_to_table <- function(df) {
     ordered_samples_table$full_name <- apply(ordered_samples_table, 1, paste, collapse = "_")
     colnames_sans_fullname <- !grepl("full_name", colnames(ordered_samples_table))
     ordered_samples_table_sans_fullname <- ordered_samples_table[, colnames_sans_fullname]
     ordered_samples_table$short_name <- apply(ordered_samples_table_sans_fullname, 1, function(row) paste0(substr(row, 1, 1), collapse = ""))
     return(ordered_samples_table)
 }

add_comparisons <- function(df) {
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

