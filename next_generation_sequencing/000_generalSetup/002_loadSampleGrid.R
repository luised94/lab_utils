#PURPOSE: Process the BMC sample grid data to create a csv file with the sample number, information and short name for further processing
#USAGE: Run as a script with Rscript and positional arguments or replace commented directory to scan line. 
#Example: Rscript ./001_R_node_processBMCSampleGridDataCSV.R <directory_to_process>
#DEPENDENCY: Assumes directory structure of directory_creation.sh and requires the BMC_sample_grid_data_*.csv file. 
#TODO: How to deal with multiple BMC files. Maybe just have it as separate directory for each sequencing run. Maybe will neeed if I do a bunch of samples.
#TODO: Need to integrate with script that creates the grid in the first place which is not in the repository currently. Need to provide as argument for the sample_ID and access with args[2] then use nrow to create the last number in the seq function.
validate_input <- function(args) {
    if(length(args) != 1){
        cat("No argument provided. One argument is required.\n")
        cat("Usage: 002_loadSampleGrid.R <directory>\n")
        cat("Example: 002_loadSampleGrid.R 240630Bel\n")
        stop()
    } 
    experiment_name <- args[1]
    directory_path <- file.path(Sys.getenv("HOME"), "data", experiment_name)
    if(!dir.exists(directory_path)) {
        stop("Directory does not exist:", directory_path)
    } else {
        cat(sprintf("Using directory %s\n:", directory_path))
    }
    return(directory_path)
}

table_has_ID_column <- function(sample_table){
    if(!("sample_ID" %in% colnames(sample_table)){
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
    fastq_file_paths <- list.files(fastq_directory_path, pattern = "\\.fastq", full.names = TRUE)
    if(length(fastq_file_paths == 0)) {
        cat(sprintf("No fastq files found in %s\n", documentation_dir_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
        q(status = 1)
    } else {
        cat(sprintf("Found %s files in %s.\n", length(fastq_file_paths), fastq_directory_path))
    }

}
load_sample_table <- function(directory_path) {
    cat("Loading sample_table from", directory_path, "\n")
    documentation_dir_path <- file.path(directory_path, "documentation")
    sample_table_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
    if(length(sample_table_path) == 0){
        cat(sprintf("No files with pattern sample_table found in %s\n", documentation_dir_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
        q(status = 1)
    } else if(length(sample_table_path) > 1){
        cat(sprintf("Multiple files with pattern sample_table found in %s\n", documentation_dir_path))
        cat(sprintf("Files found in %s\n", documentation_dir_path))
        print(sample_table)
        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
        q(status = 1)
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
    if(table_has_ID_column){
        return(sample_table)
    } else {
        determine_sample_ID(directory_path)
        sample_table <- load_sample_table(directory_path)
        return(sample_table)
    }

    #return(sample_table)
}

if(!interactive()){
    main()
}
#if (!("sample_ID" %in% colnames(sample_table))) {
#    #Determine the sample ID from fastq files
#    fastq_files <- basename(list.files(directory_to_scan, recursive = TRUE, pattern = "\\.fastq$"))
#    sample_IDs <- unique(sub(".*D24-(\\d{6}).*", "\\1", fastq_files))
#    sample_IDs <- sample_IDs[!(grepl("unmapped", sample_IDs))]
#
#    if (length(sample_IDs) == nrow(sample_table)){
#        sample_table$sample_ID <- sample_IDs
#        write.table(sample_table, file = sample_table_path, row.names = FALSE, col.names = TRUE, sep = "\t")
#    } else {
#        cat("Sample IDs determined: ", length(sample_IDs), "\n", "Rows in sample table: ", nrow(sample_table))
#        stop("Sample ID do not have the same length as the number of samples in sample_table.tsv")
#    }
#}
