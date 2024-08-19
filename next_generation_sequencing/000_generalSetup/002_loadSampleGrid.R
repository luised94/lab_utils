#PURPOSE: Process the BMC sample grid data to create a csv file with the sample number, information and short name for further processing
#USAGE: Run as a script with Rscript and positional arguments or replace commented directory to scan line. 
#Example: Rscript ./001_R_node_processBMCSampleGridDataCSV.R <directory_to_process>
#DEPENDENCY: Assumes directory structure of directory_creation.sh and requires the BMC_sample_grid_data_*.csv file. 
#TODO: How to deal with multiple BMC files. Maybe just have it as separate directory for each sequencing run. Maybe will neeed if I do a bunch of samples.
#TODO: Need to integrate with script that creates the grid in the first place which is not in the repository currently. Need to provide as argument for the sample_ID and access with args[2] then use nrow to create the last number in the seq function.
args <- commandArgs(trailingOnly = TRUE)

directory_to_scan <- paste(Sys.getenv("HOME"), "data", args[1], "documentation", sep = "/")

#directory_to_scan <- paste(Sys.getenv("HOME"), "data", "240304Bel", "documentation", sep = "/")

sample_table_path <- list.files(path = directory_to_scan, pattern = "sample_table")

sample_table <- read.table(file = sample_table_path, header = TRUE, sep = "\t")

if (!("sample_ID" %in% colnames(sample_table))) {
    #Determine the sample ID from fastq files
    fastq_files <- basename(list.files(directory_to_scan, recursive = TRUE, pattern = "\\.fastq$"))
    sample_IDs <- unique(sub(".*D24-(\\d{6}).*", "\\1", fastq_files))
    sample_IDs <- sample_IDs[!(grepl("unmapped", sample_IDs))]

    if (length(sample_IDs) == nrow(sample_table)){
        sample_table$sample_ID <- sample_IDs
        write.table(sample_table, file = sample_table_path, row.names = FALSE, col.names = TRUE, sep = "\t")
    } else {
        cat("Sample IDs determined: ", length(sample_IDs), "\n", "Rows in sample table: ", nrow(sample_table))
        stop("Sample ID do not have the same length as the number of samples in sample_table.tsv")
    }
}
