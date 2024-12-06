#STATUS:
# DESCRIPTION: Parse FastQC files in a given directory and output summarized data
# USAGE: Rscript lab_utils/next_generation_sequencing/003_fastqProcessing/006_R_node_parseFastqc.R dirname > output.log 2>&1
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

fastqc_file_path <- list.files("./qualityControl", pattern = "fastqc_data", recursive = TRUE, full.names = TRUE)
print(fastqc_file_path)
print(paste0("Number of files to process ", length(fastqc_file_path)))

number_of_files <- length(fastqc_file_path)
for (file_idx in 1:number_of_files) {
    current_time <- get_current_datetime_string()
    lines <- readLines(fastqc_file_path[file_idx])
    modules_starts <- which(grepl(">>", lines))
    modules_ends <- which(grepl(">>END_MODULE", lines))
    modules_starts <- modules_starts[!(modules_starts %in% modules_ends)]
    output_dir <- dirname(fastqc_file_path[file_idx])
    fastqc_summary <- list()
    for (module_idx in seq_along(modules_starts)) {
        module_lines <- lines[modules_starts[module_idx]:modules_ends[module_idx]]
        module_summary <- gsub(">>", "", module_lines[1])
        fastqc_summary <- append(fastqc_summary, module_summary)
        module_filename <- gsub(" ", "", strsplit(module_summary, "\t")[[1]][1])
        potential_headers <- which(grepl("^#", module_lines[2:length(module_lines)-1]))
        last_potential_header <- potential_headers[length(potential_headers)] 
        for (header_idx in potential_headers) {
            number_of_elements_in_header <- length(strsplit(module_lines[header_idx], "\t")[[1]])
            number_of_elements_in_line <- length(strsplit(module_lines[last_potential_header], "\t")[[1]])
            if(number_of_elements_in_header == number_of_elements_in_line) {
                header <- gsub("#", "", module_lines[header_idx])
                data <- read.table(text = module_lines[(header_idx+1):length(module_lines)-1],
                                   header = FALSE,
                                   col.names = strsplit(header, "\t")[[1]],
                                   sep = "\t")
            }

        }
        output_file_name <- paste0(current_time, "_", "fastqc_", module_filename, ".tab")
        output_file_path <- file.path(output_dir, output_file_name)
        write.table(data, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    fastqc_summary <- read.table(text = unlist(fastqc_summary),
                       header = FALSE,
                       col.names = c("Stat", "Value"),
                       sep = "\t")
    print(head(fastqc_summary))
    output_file_name <- paste0(current_time, "_", "fastqc_", "summary", ".tab")
    output_file_path <- file.path(output_dir, output_file_name)
    write.table(fastqc_summary, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
#Files can be deleted using: find "dir" -type f -name "*.tab" -exec rm {} + 
#Use a dry-run or wc command to see what files and how many will be deleted.
