#PURPOSE: Process the BMC sample grid data to create a csv file with the sample number, information and short name for further processing
#USAGE: Run as a script with Rscript and positional arguments or replace commented directory to scan line. 
#DEPENDENCY: Assumes directory structure of directory_creation.sh and requires the BMC_sample_grid_data_*.csv file. 
args <- commandArgs(trailingOnly = TRUE)

directory_to_scan <- paste(Sys.getenv("HOME"), "data", args[1], "documentation", sep = "/")

#directory_to_scan <- paste(Sys.getenv("HOME"), "data", "240304Bel", "documentation", sep = "/")
# TODO: How to deal with multiple BMC files. Maybe just have it as separate directory for each sequencing run. Maybe will neeed if I do a bunch of samples.
sample_grid_data_paths <- list.files(path = directory_to_scan, pattern = "BMC")

sample_info_string_names <- read.csv(sample_grid_data_paths)[,1]

sample_info_table <- as.data.frame(
				   do.call(rbind,
	                                   stringi::stri_split_regex(sample_info_string_names, pattern = "_")
                                   )
                     ) 



# Can determine in bash by using find ../240304Bel/ -type f -name "processed_*.fastq" | cut -d_ -f3 | cut -d- -f2 | sort | uniq
colnames(sample_info_table) <- c("ORC4_Rescue", "Suppressor", "Cell_Cycle", "Auxin_treatment", "antibody") 

sample_info_table$sample_ID <- seq(148001, 148049, by = 1)

#Take the first string of the first five columns to create the short name column. Easier for filtering the columns. 
sample_info_table$short_name <- apply(sample_info_table[, 1:5], 1, function(x) paste0(substr(x, 1,1), collapse = ""))

write.csv(sample_info_table, paste(directory_to_scan, "sample_info_table.csv", sep="/"), row.names = FALSE)
