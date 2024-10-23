# functions moved
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
