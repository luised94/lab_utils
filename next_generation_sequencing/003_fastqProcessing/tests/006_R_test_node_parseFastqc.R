#Run with: Rscript lab_utils/next_generation_sequencing/003_fastqProcessing/006_R_node_readInFastqc.R > output.log 2>&1

get_current_datetime_string <- function() {
                  return(format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
}
#Run with source("/home/luised94/lab_utils/next_generation_sequencing/003_fastqProcessing/006_R_node_readInFastqc.R")
setwd("/net/bmc-pub14/data/bell/users/luised94/240304Bel")
#df_sample_info <- read.delim("./documentation/sample_info_table.csv", header = TRUE, sep = ",")
#sample_ids <- df_sample_info$sample_ID
#base_dir <- "./qualityControl"

#qualityControl_dir <- list.dirs("./qualityControl", recursive = FALSE)
#contains_sample_id <- grepl(sample_ids[1], list.dirs("./qualityControl", recursive = FALSE))
#subset_qualityControl_dir <- qualityControl_dir[contains_sample_id][1]
fastqc_file_path <- list.files("./qualityControl", pattern = "fastqc_data", recursive = TRUE, full.names = TRUE)	
summary_file_path <- list.files("./qualityControl", pattern = "summary", recursive = TRUE, full.names = TRUE)	
#head(fastqc_file_path)
#head(summary_file_path)
#delimiters <- c(",", "\t", ";", " ")
max_index <- length(fastqc_file_path)
for (i in 1:max_index) {
current_time <- get_current_datetime_string()
	lines <- readLines(fastqc_file_path[i])
#	print(lines[1:5])
	data_blocks <- list()
	modules_starts <- which(grepl(">>", lines))
	modules_ends <- which(grepl(">>END_MODULE", lines))
	modules_starts <- modules_starts[!(modules_starts %in% modules_ends)]
#	print(lines[modules_starts])
#	print(lines[modules_ends])
#	seq_along(modules_starts)
output_dir <- dirname(fastqc_file_path[i])
fastqc_summary <- list()
for (i in seq_along(modules_starts)) {
module_lines <- lines[modules_starts[i]:modules_ends[i]]
module_summary <- gsub(">>", "", module_lines[1])
fastqc_summary <- append(fastqc_summary, module_summary)
module_filename <- gsub(" ", "", strsplit(module_summary, "\t")[[1]][1])
potential_headers <- which(grepl("^#", module_lines[2:length(module_lines)-1]))
last_potential_header <- potential_headers[length(potential_headers)] 
for (header_idx in potential_headers) {
	number_of_elements_in_header <- length(strsplit(module_lines[header_idx], "\t")[[1]])
	number_of_elements_in_line <- length(strsplit(module_lines[last_potential_header], "\t")[[1]])
	if(number_of_elements_in_header == number_of_elements_in_line) {
	#	message <- sprintf("Proper header is %s", i)
	#	print(message)
	header <- gsub("#", "", module_lines[header_idx])
data <- read.table(text = module_lines[(header_idx+1):length(module_lines)-1],
				   header = FALSE,
				   col.names = strsplit(header, "\t")[[1]],
				   sep = "\t")
}

}

print(data)
#message <- sprintf("#CURRENT_INDEX: %s", i)
#print(message)
output_file_name <- paste0(current_time, "_", "fastqc_", module_filename, ".tab")
output_file_path <- file.path(output_dir, output_file_name)
print(output_file_path)
#write.table(data, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
cat("\n")
#for (delim in delimiters) {
#    split_header <- strsplit(header, delim)[[1]]
#	split_first_line <-strsplit(module_lines[3], delim)[[1]]
#	if (length(split_header) > 1) {
#      message <- sprintf("Delimiter %s produces more than one column for the header. Total length is %s." , delim, length(split_header))
#	  print(message)
#	  print("Header is:")
#	  print(header)
#	  print(split_header)
#    }
#
#	if (length(split_first_line) > 1) {
#      message <- sprintf("Delimiter %s produces more than one column for the header. Total length is %s." , delim, length(split_first_line))
#	  print(message)
#	  print("Line is:")	
#      print(module_lines[3])
#	  print(split_first_line)
#    }
#
#}
#
	}
print(unlist(fastqc_summary))
output_file_name <- paste0(current_time, "_", "fastqc_", "summary", ".tab")
output_file_path <- file.path(output_dir, output_file_name)
print(output_file_path)
#write.table(fastqc_summary, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
#write.table(data, file = output_file_name, sep = "\t", quotes = FALSE, col.names = TRUE)
#if ((strsplit(module_summary, "\t")[[1]][1] == "Sequence Length Distribution") == FALSE) {
#} else {
#print("#SKIPPED Sequence Length Distribution")
#}
