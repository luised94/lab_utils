output_tables_in_list <- function(experiment_directory, list_of_tables, OUTPUT_TABLE = FALSE){
        experiment_name <- basename(experiment_directory)
        if (!(typeof(list_of_tables) == "list")){
            stop("Argument must be a list.")
        }
        names_of_tables <- names(list_of_tables)
        for (name_of_table in names_of_tables){
            output_table <- sample_config_output[[name_of_table]]
            cat("============\n")
            print(head(output_table))
            output_file_path <- file.path(experiment_directory, "documentation", paste(experiment_name, "_", name_of_table, ".tsv", sep = ""))
            cat(sprintf("Outputting to %s: \n", output_file_path))
            if(OUTPUT_TABLE) {
                write.table(output_table, file = output_file_path, sep = "\t", row.names = FALSE)
            } else {
                cat("Skip writing table. MODIFY OUTPUT_TABLE value to output.\n")
            }
        }
}
