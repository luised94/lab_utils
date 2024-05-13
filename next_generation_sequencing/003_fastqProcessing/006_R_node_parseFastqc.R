#Run with: Rscript lab_utils/next_generation_sequencing/003_fastqProcessing/006_R_node_parseFastqc.R dirname > output.log 2>&1
args <- commandArgs(trailingOnly = TRUE)

get_current_datetime_string <- function() {
                  return(format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
}

find_the_directory <- function(dirname){
	data_directory <- paste0(Sys.getenv("HOME"), "/data")
	contains_dirname <- grepl(dirname, list.dirs(data_directory, recursive = FALSE))
	return(list.dirs(data_directory, recursive = FALSE)[contains_dirname])
}
print(args[1])
directory_to_process <- find_the_directory(args[1])
setwd(directory_to_process)

fastqc_file_path <- list.files("./qualityControl", pattern = "fastqc_data", recursive = TRUE, full.names = TRUE)	

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
		
		head(data)
		output_file_name <- paste0(current_time, "_", "fastqc_", module_filename, ".tab")
		output_file_path <- file.path(output_dir, output_file_name)
		write.table(data, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	}
	fastqc_summary <- read.table(text = unlist(fastqc_summary),
					   header = FALSE,
					   col.names = c("Stat", "Value"),
					   sep = "\t")
	output_file_name <- paste0(current_time, "_", "fastqc_", "summary", ".tab")
	output_file_path <- file.path(output_dir, output_file_name)
	write.table(fastqc_summary, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
#Files can be deleted using: find "dir" -type f -name "*.tab" -exec rm {} + 
#Use a dry-run or wc command to see what files and how many will be deleted.
