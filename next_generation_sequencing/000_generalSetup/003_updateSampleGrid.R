#STATUS:
#PURPOSE: Process the BMC sample grid data to create a csv file with the sample number, information and short name for further processing
#USAGE: Run as a script with Rscript and positional arguments or replace commented directory to scan line. 
#Example: Rscript ./001_R_node_processBMCSampleGridDataCSV.R <directory_to_process>
#DEPENDENCY: Assumes directory structure of directory_creation.sh and requires the BMC_sample_grid_data_*.csv file. 
#TODO: How to deal with multiple BMC files. Maybe just have it as separate directory for each sequencing run. Maybe will neeed if I do a bunch of samples.
#TODO: Need to integrate with script that creates the grid in the first place which is not in the repository currently. Need to provide as argument for the sample_ID and access with args[2] then use nrow to create the last number in the seq function.
validate_input <- function(args) {
    if(length(args) != 1){
        cat("No argument provided. One argument is required.\n")
        cat("Usage: Rscript 003_updateSampleGrid.R <directory>\n")
        cat("Example: Rscript 003_updateSampleGrid.R 240630Bel\n")
        stop()
    } 
    experiment_name <- args[1]
    directory_path <- file.path(Sys.getenv("HOME"), "data", experiment_name)
    if(!dir.exists(directory_path)) {
        stop("Directory does not exist:", directory_path)
    } else {
        cat(sprintf("Using directory %s:\n", directory_path))
    }
    return(directory_path)
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

determine_sample_id <- function(directory_path) {
    fastq_directory_path <- file.path(directory_path, "fastq")
    fastq_file_paths <- list.files(fastq_directory_path, pattern = "*.fastq", full.names = TRUE)
    if(length(fastq_file_paths) == 0) {
        cat(sprintf("No fastq files found in %s\n", fastq_directory_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project and download files from BMC\n")
        stop()
    } else {
        cat(sprintf("Found %s files in %s.\n", length(fastq_file_paths), fastq_directory_path))
    }
    fastq_file_names <- basename(fastq_file_paths)
    ID_regex <- "\\d{5,6}"
    fastq_split_string_list <- strsplit(fastq_file_names, "_|-")
    sample_IDs <- lapply(fastq_split_string_list, function(fastq_split_string_list) {
        for(split_string in fastq_split_string_list) {
            if(grepl(ID_regex, split_string)) {
                return(split_string)
            }
        }
    })
    if(!all(unlist(lapply(sample_IDs, length)) == 1)) {
        cat("At least one of the files did not extract exactly one sample ID.\n")
        cat("Files with problems:\n")
        print(fastq_file_names[unlist(lapply(sample_IDs, length)) != 1])
        cat("Verify sample names. Redownload from BMC if necessary.\n")
        cat(sprintf("Regex pattern used %s:\n", ID_regex))
        stop()
    } else {
        sample_IDs <- unlist(sample_IDs)
        cat(sprintf("Found %s sample_IDs.\n", length(sample_IDs)))
        cat(sprintf("First sample_ID: %s\n",sample_IDs[1]))
        cat(sprintf("Last sample_ID: %s\n", sample_IDs[length(sample_IDs)]))
        cat("Returning sample_ID array.\n")
        return(sample_IDs)
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

load_sample_table <- function(directory_path) {
    cat("Loading sample_table from", directory_path, "\n")
    documentation_dir_path <- file.path(directory_path, "documentation")
    sample_table_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
    if(length(sample_table_path) == 0){
        cat(sprintf("No files with pattern sample_table found in %s\n", documentation_dir_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
        stop()
    } else if(length(sample_table_path) > 1){
        cat(sprintf("Multiple files with pattern sample_table found in %s\n", documentation_dir_path))
        cat(sprintf("Files found in %s\n", documentation_dir_path))
        print(sample_table)
        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
        stop()
    }
    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
    cat(sprintf("Reading %s\n", sample_table_path))
    cat("Head of sample_table\n")
    print(head(sample_table))
    return(sample_table)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    directory_path <- validate_input(args)
    sample_table <- load_sample_table(directory_path)
    if(table_has_ID_column(sample_table)){
        cat("Sample table has sample ID column")
        return(sample_table)
    } else {
        cat("Sample table doesnt have sample ID column.\n")
        cat("Determine sample IDs and overwrite sample table.\n")
        documentation_dir_path <- file.path(directory_path, "documentation")
        output_file_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
        sample_IDs <- determine_sample_id(directory_path)
        modify_and_output_table(sample_table, sample_IDs, output_file_path)
        sample_table <- load_sample_table(directory_path)
        return(sample_table)
    }
}

if(!interactive()){
    main()
}
