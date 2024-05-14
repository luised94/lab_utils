#Run with: Rscript lab_utils/next_generation_sequencing/004_bamProcessing/003_R_node_checkQcBam.R dirname > output.log 2>&1
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 ) { print("No argument provided. Provide a directory to process") ; q() }

get_current_datetime_string <- function() {
                  return(format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
}

find_the_directory <- function(dirname){
    data_directory <- paste0(Sys.getenv("HOME"), "/data")
    contains_dirname <- grepl(dirname, list.dirs(data_directory, recursive = FALSE))
    directory_path <- list.dirs(data_directory, recursive = FALSE)[contains_dirname]
    if (length(directory_path) > 0 && dir.exists(directory_path)) { return(directory_path) }
    else { print(paste0("Directory ", dirname, " not found")) ; return(NULL) }
}

directory_to_process <- find_the_directory(args[1])
if (is.null(directory_to_process)) { print("Directory to process is null. Exiting script.") ; q() }

message <- sprintf("Directory to process is %s.", directory_to_process)
print(message)
setwd(directory_to_process)

fastqc_file_path <- list.files("./qualityControl", pattern = "bam", recursive = FALSE, full.names = TRUE)
print(fastqc_file_path)
print(paste0("Number of files to process ", length(fastqc_file_path)))

